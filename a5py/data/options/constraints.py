"""Tools to enforce constraints on option parameters.

The parameters in options are subject to constraints. In order to have a single
source of truth, these constraints are defined in the docstring of each
dataclass attribute.
"""
import re

def parse(doc):
    """Parse constraints from attribute docstring.

    Assuming that the attribute's docstring has the following format:

    "<info>: <constraint>, default=<default>.",

    this function extracts the constraint and the default value along with the
    property name. These can be used to generate XML file storing the options.

    Parameters
    ----------
    doc : str
        Attribute docstring.

    Returns
    -------
    dict : {"name", "constraint", "default"}
        The extracted summary.
    """
    parts = []
    for line in doc.strip().splitlines():
        if not line.strip():
            break
        parts.append(line.strip())
    summary_in_one_line = " ".join(parts)
    match_regex = re.match(
        r"^(?P<info>[^:]+): (?P<constraint>.*), default=(?P<default>[^.]+)\.",
        summary_in_one_line
        )
    if not match_regex:
        raise ValueError(f"Docstring format invalid: {summary_in_one_line}")
    return {
        "constraint": match_regex.group("constraint").strip(),
        "default": match_regex.group("default").strip(),
    }


def convert(constraint: str):
    """Convert human-readable constraint string to machine-readable.

    Parameters
    ----------
    constraint : str
        A string representation of a constraint.

        Acceptable formats are: "{value1, value2, ...}" for sets, "(condition)"
        for scalars, and "[condition, condition, ...]" for lists.

    Returns
    -------
    kind : str
        One of "set", "scalar", or "list".
    constraints : list
        A list of constraints.
    """
    constraint = constraint.strip()
    if constraint.startswith("{") and constraint.endswith("}"):
        return ("set", [eval(x) for x in constraint[1:-1].split(",")])

    if constraint.startswith("(") and constraint.endswith(")"):
        return ("scalar", constraint[1:-1].strip())

    if constraint.startswith("[") and constraint.endswith("]"):
        return ("list", [c.strip() for c in constraint[1:-1].split(",")])

    raise ValueError(f"Unknown constraint format: {constraint}")


def enforce(value, constraint: str):
    """Check if a value satisfies the constraint.

    Parameters
    ----------
    value : Any
        The value to check.
    constraint : str
        A string representation of a constraint.

        The following types of constraints are enforced:

        - {a, b, c, ...}: value must be in the set.
        - (a), (a > 0): a scalar value must satisfy a condition.
        - [a, b], [a > 0, ...] value must either be a list of fixed length
          (or a scalar or a list if "..." is used) where each element must
          satisfy a condition.

    Returns
    -------
    accepted : bool
        True if the value satisfies the constraint, False otherwise.
    """
    kind, constraints = convert(constraint)

    if kind == "set":
        items_in_the_set = constraints
        return value in items_in_the_set

    if kind == "scalar":
        if constraints == "x":
            return True
        return eval(f"{value} {constraints}")

    if kind == "list":
        if "..." in constraints[-1]:
            constraint_to_apply_all = constraints[0].replace("a", "x")
            if not isinstance(value, (list, tuple)):
                # Scalars are accepted for non-fixed lists
                value = [value]
            return all(eval(constraint_to_apply_all, {"x": v}) for v in value)
        else:
            value_has_expected_length = ( isinstance(value, (list, tuple)) and
                                          len(value) == len(constraints) )
            if not value_has_expected_length:
                return False
            for list_item, item_constraint in zip(value, constraints):
                expr = item_constraint.replace("a", "x").replace("b", "x")
                if not eval(expr, {"x": list_item}):
                    return False
            return True
