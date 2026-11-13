"""What an audit reports, and how.

A finding is either a problem, which makes part of a scheme wrong or
impossible to type, or a recommendation, which leaves it usable. Each also
carries the action it asks for: modify an existing SNP or node, extend the
scheme with SNPs or nodes it does not have, or check something the tool cannot
settle on its own.
"""

from dataclasses import dataclass

PROBLEM = "problem"
RECOMMENDATION = "recommendation"

MODIFY = "modify"
EXTEND = "extend"
CHECK = "check"


@dataclass(frozen=True)
class Finding:
    """One thing the audit has to say about a scheme."""

    HEADER = ("severity", "action", "category", "node", "detail")

    severity: str
    action: str
    category: str
    node: str
    detail: str

    def row(self) -> list[str]:
        return [self.severity, self.action, self.category, self.node, self.detail]

    def __str__(self) -> str:
        where = f" [{self.node}]" if self.node else ""
        return f"{self.action}{where}: {self.detail}"


def problem(action: str, category: str, node: str, detail: str) -> Finding:
    return Finding(PROBLEM, action, category, node, detail)


def recommendation(action: str, category: str, node: str, detail: str) -> Finding:
    return Finding(RECOMMENDATION, action, category, node, detail)
