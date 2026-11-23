"""The hierarchy of a canSNP scheme.

A scheme is a rooted tree of named nodes, each the branch one canSNP marks.
Which genomes sit at which node is a separate file, the genome map, as NCBI
Taxonomy keeps nodes.dmp apart from its accession table and FlexTaxD keeps its
taxonomy apart from --genomeid2taxid: a node may hold no genome yet, and a
genome may sit at an internal node.

``Hierarchy`` is built from a child and parent table or from Newick.
"""

import re
from collections.abc import Iterable, Mapping, Sequence

#: A non-negative decimal number, as tree builders write support values:
#: 95, 95.0, 0.95, .95, 1. or 1e-1.
NUMBER = r"(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
#: One number, or several separated by slashes as IQ-TREE writes SH-aLRT and
#: ultrafast bootstrap together; a field may be empty where a value is missing.
SUPPORT_LABEL = re.compile(rf"^(?=.*\d)(?:{NUMBER})?(?:/(?:{NUMBER})?)*$")
#: Marks a token that was a quoted label, so that a quoted "(" or ":" is a
#: name and never structure.
QUOTED = "\x00"


class HierarchyError(ValueError):
    """The hierarchy cannot be read, or is not a tree."""


def _newick_tokens(text: str) -> list[str]:
    """Split Newick text into brackets, commas, the closing semicolon, labels
    and branch lengths, the latter kept with their leading colon."""
    tokens: list[str] = []
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


def _support_value(label: str) -> float | None:
    """The support a numeric internal label carries, or None for a name.

    Of several values separated by slashes the last is taken. A label whose
    last field is empty, as ``80/``, is a support label without that value,
    and gives NaN, which the caller reads as no support.
    """
    if not SUPPORT_LABEL.match(label):
        return None
    last = label.split("/")[-1]
    return float(last) if last else float("nan")


class Hierarchy:
    """A rooted tree of node names.

    Built through one of the ``from_`` constructors, which all check that the
    result is a single rooted tree: no node with two parents, no cycle and no
    second root.
    """

    def __init__(self, pairs: Iterable[tuple[str, str | None]]) -> None:
        self.parents: dict[str, str | None] = {}
        self.children: dict[str, list[str]] = {}
        self.order: list[str] = []
        # Read from a tree that carries them: the support of each internal
        # node as a fraction, and the length of the branch above each node.
        self.support: dict[str, float] = {}
        self.lengths: dict[str, float] = {}

        for child, parent in pairs:
            child = str(child).strip()
            parent = None if parent is None else str(parent).strip()
            if not child:
                raise HierarchyError("the hierarchy holds a node without a name")
            if child == parent:
                raise HierarchyError(f"node '{child}' is its own parent")
            if child in self.parents and self.parents[child] != parent:
                raise HierarchyError(
                    f"node '{child}' is given more than one parent: {self.parents[child]}, {parent}"
                )
            if child not in self.parents:
                self.parents[child] = parent
                self.order.append(child)

        for child, parent in self.parents.items():
            if parent is not None:
                self.children.setdefault(parent, []).append(child)
        for node in list(self.children):
            if node not in self.parents:
                self.parents[node] = None
                self.order.append(node)
        for node in self.order:
            self.children.setdefault(node, [])

        self._root, self._nodes, self._leaves = self._check_is_a_tree()

    def _check_is_a_tree(self) -> tuple[str, tuple[str, ...], tuple[str, ...]]:
        """Check the nodes form one rooted tree; return its root, its nodes
        with parents before their children, and its leaves in that order."""
        roots = self.roots
        if not roots:
            raise HierarchyError("the hierarchy has no root; every node has a parent")
        if len(roots) > 1:
            raise HierarchyError(
                f"the hierarchy has {len(roots)} roots where one is expected: {', '.join(sorted(roots))}"
            )

        reached: set[str] = set()
        ordered: list[str] = []
        stack = list(roots)
        while stack:
            node = stack.pop()
            if node in reached:
                raise HierarchyError(f"the hierarchy is not a tree; node '{node}' is reached twice")
            reached.add(node)
            ordered.append(node)
            stack.extend(reversed(self.children[node]))

        unreached = set(self.order) - reached
        if unreached:
            raise HierarchyError(f"these nodes are not connected to the root: {', '.join(sorted(unreached))}")
        return (roots[0], tuple(ordered), tuple(node for node in ordered if not self.children[node]))

    # -- constructors --------------------------------------------------------

    @classmethod
    def from_child_parent(cls, pairs: Iterable[tuple[str, str | None]]) -> "Hierarchy":
        """Build from (child, parent) pairs, as a repgenr hierarchy file lists them.

        The root is the node that appears only as a parent, so it needs no row
        of its own; a row that leaves the parent empty names a root explicitly.
        """
        return cls([(child, parent if parent and str(parent).strip() else None) for child, parent in pairs])

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
        label; a label whose last field is empty gives no support value.
        Unnamed internal nodes are given a generated name, since a canSNP
        scheme addresses its branches by name. Single-quoted labels may hold
        any character, and comments in square brackets are skipped. The tree
        is read without recursion, so its depth is not limited.
        """
        tokens = _newick_tokens(text)
        if not tokens:
            raise HierarchyError("the Newick string holds no nodes")
        structure = {"(", ")", ",", ";"}

        def is_label(token: str | None) -> bool:
            return token is not None and (
                token.startswith(QUOTED) or (token not in structure and not token.startswith(":"))
            )

        # Generated names must not take a name the tree already uses.
        taken = {token.lstrip(QUOTED) for token in tokens if is_label(token)}
        names_seen: set[str] = set()

        pairs: list[tuple[str, str | None]] = []
        support: dict[str, float] = {}
        lengths: dict[str, float] = {}
        generated = 0
        position = 0

        def peek() -> str | None:
            return tokens[position] if position < len(tokens) else None

        def finish_node(children: list[str]) -> str:
            """Read the label and branch length after a node, and name it."""
            nonlocal position, generated
            label = ""
            quoted = False
            if is_label(peek()):
                quoted = tokens[position].startswith(QUOTED)
                label = tokens[position].lstrip(QUOTED)
                position += 1
            length = None
            following = peek()
            if following is not None and not following.startswith(QUOTED) and following.startswith(":"):
                try:
                    length = float(tokens[position][1:])
                except ValueError:
                    raise HierarchyError(f"'{tokens[position][1:]}' is not a branch length") from None
                position += 1

            value = _support_value(label) if children and support_labels and not quoted else None
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
                    f"the Newick string names '{name}' twice; every node of a tree needs its own name"
                )
            names_seen.add(name)
            if value is not None and value == value:
                support[name] = value
            if length is not None:
                lengths[name] = length

            for child in children:
                pairs.append((child, name))
            return name

        # Read without recursion, so that a deep tree, such as a caterpillar
        # of thousands of genomes, is read as any other: each open group is a
        # list of the children read so far.
        groups: list[list[str]] = []
        while True:
            while peek() == "(":
                position += 1
                groups.append([])
            name = finish_node([])
            while groups:
                if peek() == ",":
                    position += 1
                    groups[-1].append(name)
                    break
                if peek() == ")":
                    position += 1
                    name = finish_node([*groups.pop(), name])
                    continue
                if peek() is None:
                    raise HierarchyError("the Newick string ends inside a group")
                raise HierarchyError(f"unexpected '{peek()}' in the Newick string")
            else:
                break
        root = name
        if peek() == ";":
            position += 1
        if peek() is not None:
            raise HierarchyError("the Newick string holds more than one tree, or text after its ';'")
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
    def roots(self) -> list[str]:
        return [node for node, parent in self.parents.items() if parent is None]

    # The root, the node order and the leaves are found once, when the tree
    # is checked; a Hierarchy is not changed after it is built.

    @property
    def root(self) -> str:
        return self._root

    @property
    def nodes(self) -> list[str]:
        """Every node, parents before their children."""
        return list(self._nodes)

    @property
    def leaves(self) -> list[str]:
        return list(self._leaves)

    def descendants(self, node: str) -> list[str]:
        """Every node below this one, excluding itself, parents before their
        children."""
        if node not in self.parents:
            raise HierarchyError(f"the hierarchy has no node '{node}'")
        found: list[str] = []
        stack = list(reversed(self.children[node]))
        while stack:
            current = stack.pop()
            found.append(current)
            stack.extend(reversed(self.children[current]))
        return found

    def leaves_below(self, node: str) -> list[str]:
        """The leaves at or below this node."""
        below = [node, *self.descendants(node)]
        return [name for name in below if not self.children[name]]

    def path(self, node: str) -> list[str]:
        """The nodes from the root down to this one, inclusive."""
        if node not in self.parents:
            raise HierarchyError(f"the hierarchy has no node '{node}'")
        path = [node]
        parent = self.parents[node]
        while parent is not None:
            path.append(parent)
            parent = self.parents[parent]
        return list(reversed(path))

    def isolates(self, clusters: Mapping[str, Sequence[str]]) -> dict[str, list[str]]:
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
                f"the genome map names clusters that are not in the hierarchy: {', '.join(unknown)}"
            )

        isolates: dict[str, list[str]] = {}
        for node in reversed(self.nodes):
            own = list(clusters.get(node, []))
            for child in self.children[node]:
                own.extend(isolates[child])
            # Keep the first mention of each isolate, so that a cluster named at
            # two depths does not inflate the count.
            isolates[node] = list(dict.fromkeys(own))
        return isolates

    def __len__(self) -> int:
        return len(self.parents)

    def __contains__(self, node: str) -> bool:
        return node in self.parents
