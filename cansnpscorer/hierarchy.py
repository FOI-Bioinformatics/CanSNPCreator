"""The hierarchy of a canSNP scheme.

A scheme is a rooted tree whose leaves are clusters of isolates and whose
internal nodes are the branches a canSNP marks. The hierarchy arrives in one of
several shapes, and the tools used to treat them differently: the nested form
was walked recursively so that an internal node received the isolates of
everything below it, while the child and parent table was handed straight to
that walker, which could not read it and left every internal node empty. No
canSNP could then be assigned to a branch, which is what the tools exist to do.

``Hierarchy`` is the one representation they all become.
"""

import json
import re
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

NEWICK_TOKEN = re.compile(r"[(),;]|[^(),;:\s]+(?::[^(),;\s]+)?")


class HierarchyError(ValueError):
    """The hierarchy cannot be read, or is not a tree."""


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
    def from_nested(cls, data) -> "Hierarchy":
        """Build from the nested form, as the YAML and JSON files hold it.

        A mapping is a node with its children below it; a sequence holds
        siblings; a bare value is a leaf. A node with nothing below it, written
        either as ``name`` or as ``name:``, is a leaf.
        """
        pairs: List[Tuple[str, Optional[str]]] = []

        def walk(node, parent: Optional[str]) -> None:
            if isinstance(node, Mapping):
                for name, below in node.items():
                    pairs.append((name, parent))
                    walk(below, name)
            elif isinstance(node, (list, tuple, set)):
                for item in node:
                    walk(item, parent)
            elif node is not None:
                pairs.append((node, parent))

        walk(data, None)
        if not pairs:
            raise HierarchyError("the hierarchy holds no nodes")
        return cls(pairs)

    @classmethod
    def from_json(cls, text: str) -> "Hierarchy":
        try:
            return cls.from_nested(json.loads(text))
        except json.JSONDecodeError as error:
            raise HierarchyError(f"the hierarchy is not valid JSON: {error}") from None

    @classmethod
    def from_newick(cls, text: str) -> "Hierarchy":
        """Build from a Newick string whose internal nodes are named.

        Branch lengths are accepted and discarded. Unnamed internal nodes are
        given a generated name, since a canSNP scheme addresses its branches by
        name.
        """
        tokens = NEWICK_TOKEN.findall(text.strip())
        if not tokens:
            raise HierarchyError("the Newick string holds no nodes")

        pairs: List[Tuple[str, Optional[str]]] = []
        generated = 0

        def name_of(token: str) -> str:
            return token.split(":")[0].strip()

        position = 0

        def parse_node() -> str:
            nonlocal position, generated
            children: List[str] = []
            if position < len(tokens) and tokens[position] == "(":
                position += 1
                while True:
                    children.append(parse_node())
                    if position >= len(tokens):
                        raise HierarchyError("the Newick string ends inside a group")
                    if tokens[position] == ",":
                        position += 1
                        continue
                    if tokens[position] == ")":
                        position += 1
                        break
                    raise HierarchyError(f"unexpected '{tokens[position]}' in the Newick string")

            name = ""
            if position < len(tokens) and tokens[position] not in {"(", ")", ",", ";"}:
                name = name_of(tokens[position])
                position += 1
            if not name:
                generated += 1
                name = f"node{generated}"

            for child in children:
                pairs.append((child, name))
            return name

        root = parse_node()
        pairs.append((root, None))
        return cls(pairs)

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
            clusters: The isolates of each leaf, as a group file lists them.

        Returns:
            The isolates at or below every node, so that an internal node
            carries everything below it. Leaves with no isolates are reported
            as empty rather than dropped, so that the caller can say which
            nodes of the scheme cannot be typed.
        """
        unknown = sorted(set(clusters) - set(self.parents))
        if unknown:
            raise HierarchyError(
                f"the group file names clusters that are not in the hierarchy: "
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
