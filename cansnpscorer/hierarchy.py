"""The hierarchy of a canSNP scheme.

A scheme is a rooted tree of named nodes, each the branch one canSNP marks.
Which genomes sit at which node is a separate file, the genome map, as NCBI
Taxonomy keeps nodes.dmp apart from its accession table and FlexTaxD keeps its
taxonomy apart from --genomeid2taxid: a node may hold no genome yet, and a
genome may sit at an internal node.

``Hierarchy`` is built from a child and parent table or from Newick.
"""

import re
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

SUPPORT_LABEL = re.compile(r"^\d+(\.\d+)?(/\d+(\.\d+)?)*$")
#: Marks a token that was a quoted label, so that a quoted "(" or ":" is a
#: name and never structure.
QUOTED = "\x00"


class HierarchyError(ValueError):
    """The hierarchy cannot be read, or is not a tree."""


def _newick_tokens(text: str) -> List[str]:
    """Split Newick text into brackets, commas, the closing semicolon, labels
    and branch lengths, the latter kept with their leading colon."""
    tokens: List[str] = []
    index = 0
    text = text.strip()
    while index < len(text):
        character = text[index]
        if character.isspace():
            index += 1
        elif character in "(),;":
            tokens.append(character)
            index += 1
        elif character == "'":
            end = index + 1
            label = []
            while end < len(text):
                if text[end] == "'" and end + 1 < len(text) and text[end + 1] == "'":
                    label.append("'")
                    end += 2
                    continue
                if text[end] == "'":
                    break
                label.append(text[end])
                end += 1
            else:
                raise HierarchyError("the Newick string holds an unclosed quoted label")
            tokens.append(QUOTED + "".join(label))
            index = end + 1
        elif character == "[":
            end = text.find("]", index)
            if end < 0:
                raise HierarchyError("the Newick string holds an unclosed comment")
            index = end + 1
        else:
            end = index
            if character == ":":
                end += 1
            while end < len(text) and text[end] not in "(),;:'[" and not text[end].isspace():
                end += 1
            tokens.append(text[index:end])
            index = end
    return tokens


def _support_value(label: str) -> Optional[float]:
    """The support a numeric internal label carries, or None for a name."""
    if not SUPPORT_LABEL.match(label):
        return None
    return float(label.split("/")[-1])


class Hierarchy:
    """A rooted tree of node names.

    Built through one of the ``from_`` constructors, which all check that the
    result is a single rooted tree: no node with two parents, no cycle and no
    second root.
    """

    def __init__(self, pairs: Iterable[Tuple[str, Optional[str]]]) -> None:
        self.parents: Dict[str, Optional[str]] = {}
        self.children: Dict[str, List[str]] = {}
        self.order: List[str] = []
        # Read from a tree that carries them: the support of each internal
        # node as a fraction, and the length of the branch above each node.
        self.support: Dict[str, float] = {}
        self.lengths: Dict[str, float] = {}

        for child, parent in pairs:
            child = str(child).strip()
            parent = None if parent is None else str(parent).strip()
            if not child:
                raise HierarchyError("the hierarchy holds a node without a name")
            if child == parent:
                raise HierarchyError(f"node '{child}' is its own parent")
            if child in self.parents and self.parents[child] != parent:
                raise HierarchyError(
                    f"node '{child}' is given more than one parent: "
                    f"{self.parents[child]}, {parent}")
            if child not in self.parents:
                self.parents[child] = parent
                self.order.append(child)

        for child, parent in self.parents.items():
            if parent is not None:
                self.children.setdefault(parent, []).append(child)
        for node in list(self.children):
            self.parents.setdefault(node, None)
            if node not in self.order:
                self.order.append(node)
        for node in self.order:
            self.children.setdefault(node, [])

        self._check_is_a_tree()

    def _check_is_a_tree(self) -> None:
        roots = self.roots
        if not roots:
            raise HierarchyError("the hierarchy has no root; every node has a parent")
        if len(roots) > 1:
            raise HierarchyError(
                f"the hierarchy has {len(roots)} roots where one is expected: "
                f"{', '.join(sorted(roots))}")

        reached: Set[str] = set()
        stack = list(roots)
        while stack:
            node = stack.pop()
            if node in reached:
                raise HierarchyError(f"the hierarchy is not a tree; node '{node}' is reached twice")
            reached.add(node)
            stack.extend(self.children[node])

        unreached = set(self.order) - reached
        if unreached:
            raise HierarchyError(
                f"these nodes are not connected to the root: {', '.join(sorted(unreached))}")

    # -- constructors --------------------------------------------------------

    @classmethod
    def from_child_parent(cls, pairs: Iterable[Tuple[str, Optional[str]]]) -> "Hierarchy":
        """Build from (child, parent) pairs, as a repgenr hierarchy file lists them.

        The root is the node that appears only as a parent, so it needs no row
        of its own; a row that leaves the parent empty names a root explicitly.
        """
        return cls([
            (child, parent if parent and str(parent).strip() else None)
            for child, parent in pairs
        ])

    @classmethod
    def from_newick(cls, text: str, support_labels: bool = False) -> "Hierarchy":
        """Build from a Newick string.

        Branch lengths are kept in ``lengths``. Newick does not say whether an
        internal label is a name or a support value, so the caller does: in a
        scheme's hierarchy every label is a name, as a node called 2.1 must
        stay 2.1. With ``support_labels``, as for a genome tree, an internal
        label that is a number, or numbers separated by slashes as IQ-TREE
        writes SH-aLRT and ultrafast bootstrap together, is a support value: it
        is kept in ``support``, as a fraction, from the last number of the
        label. Unnamed internal nodes are given a generated name, since a canSNP
        scheme addresses its branches by name. Single-quoted labels may hold
        any character.
        """
        tokens = _newick_tokens(text)
        if not tokens:
            raise HierarchyError("the Newick string holds no nodes")
        structure = {"(", ")", ",", ";"}

        def is_label(token: Optional[str]) -> bool:
            return token is not None and (token.startswith(QUOTED) or (
                token not in structure and not token.startswith(":")))

        # Generated names must not take a name the tree already uses.
        taken = {token.lstrip(QUOTED) for token in tokens if is_label(token)}
        names_seen: Set[str] = set()

        pairs: List[Tuple[str, Optional[str]]] = []
        support: Dict[str, float] = {}
        lengths: Dict[str, float] = {}
        generated = 0
        position = 0

        def peek() -> Optional[str]:
            return tokens[position] if position < len(tokens) else None

        def parse_node() -> str:
            nonlocal position, generated
            children: List[str] = []
            if peek() == "(":
                position += 1
                while True:
                    children.append(parse_node())
                    if peek() == ",":
                        position += 1
                        continue
                    if peek() == ")":
                        position += 1
                        break
                    if peek() is None:
                        raise HierarchyError("the Newick string ends inside a group")
                    raise HierarchyError(f"unexpected '{peek()}' in the Newick string")

            label = ""
            quoted = False
            if is_label(peek()):
                quoted = tokens[position].startswith(QUOTED)
                label = tokens[position].lstrip(QUOTED)
                position += 1
            length = None
            if peek() is not None and not peek().startswith(QUOTED) and peek().startswith(":"):
                try:
                    length = float(tokens[position][1:])
                except ValueError:
                    raise HierarchyError(
                        f"'{tokens[position][1:]}' is not a branch length") from None
                position += 1

            value = (_support_value(label) if children and support_labels and not quoted
                     else None)
            name = label if value is None else ""
            if not name:
                if not children:
                    raise HierarchyError("the Newick string holds a leaf without a name")
                generated += 1
                while f"node{generated}" in taken:
                    generated += 1
                name = f"node{generated}"
            if name in names_seen:
                raise HierarchyError(
                    f"the Newick string names '{name}' twice; every node of a tree needs its "
                    f"own name")
            names_seen.add(name)
            if value is not None:
                support[name] = value
            if length is not None:
                lengths[name] = length

            for child in children:
                pairs.append((child, name))
            return name

        root = parse_node()
        if peek() == ";":
            position += 1
        if peek() is not None:
            raise HierarchyError(
                "the Newick string holds more than one tree, or text after its ';'")
        pairs.append((root, None))

        hierarchy = cls(pairs)
        # Supports written as percentages are read as fractions, so that a
        # threshold means the same whichever the tree builder wrote.
        if any(value > 1.0 for value in support.values()):
            support = {name: value / 100.0 for name, value in support.items()}
        hierarchy.support = support
        hierarchy.lengths = lengths
        return hierarchy

    # -- questions about the tree -------------------------------------------

    @property
    def roots(self) -> List[str]:
        return [node for node, parent in self.parents.items() if parent is None]

    @property
    def root(self) -> str:
        return self.roots[0]

    @property
    def nodes(self) -> List[str]:
        """Every node, parents before their children."""
        ordered: List[str] = []
        stack = [self.root]
        while stack:
            node = stack.pop(0)
            ordered.append(node)
            stack = self.children[node] + stack
        return ordered

    @property
    def leaves(self) -> List[str]:
        return [node for node in self.nodes if not self.children[node]]

    def descendants(self, node: str) -> List[str]:
        """Every node below this one, excluding itself."""
        if node not in self.parents:
            raise HierarchyError(f"the hierarchy has no node '{node}'")
        found: List[str] = []
        stack = list(self.children[node])
        while stack:
            current = stack.pop(0)
            found.append(current)
            stack = self.children[current] + stack
        return found

    def leaves_below(self, node: str) -> List[str]:
        """The leaves at or below this node."""
        below = [node, *self.descendants(node)]
        return [name for name in below if not self.children[name]]

    def path(self, node: str) -> List[str]:
        """The nodes from the root down to this one, inclusive."""
        if node not in self.parents:
            raise HierarchyError(f"the hierarchy has no node '{node}'")
        path = [node]
        parent = self.parents[node]
        while parent is not None:
            path.append(parent)
            parent = self.parents[parent]
        return list(reversed(path))

    def isolates(self, clusters: Mapping[str, Sequence[str]]) -> Dict[str, List[str]]:
        """The isolates below each node.

        Args:
            clusters: The isolates of each leaf, as a genome map lists them.

        Returns:
            The isolates at or below every node, so that an internal node
            carries everything below it. Leaves with no isolates are reported
            as empty rather than dropped, so that the caller can say which
            nodes of the scheme cannot be typed.
        """
        unknown = sorted(set(clusters) - set(self.parents))
        if unknown:
            raise HierarchyError(
                f"the genome map names clusters that are not in the hierarchy: "
                f"{', '.join(unknown)}")

        isolates: Dict[str, List[str]] = {}
        for node in reversed(self.nodes):
            own = list(clusters.get(node, []))
            for child in self.children[node]:
                own.extend(isolates[child])
            # Keep the first mention of each isolate, so that a cluster named at
            # two depths does not inflate the count.
            seen: Set[str] = set()
            isolates[node] = [name for name in own if not (name in seen or seen.add(name))]
        return isolates

    def __len__(self) -> int:
        return len(self.parents)

    def __contains__(self, node: str) -> bool:
        return node in self.parents
