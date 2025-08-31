"""Tools to enforce constraints on option parameters.

The parameters in options are subject to constraints. In order to have a single
source of truth, these constraints are defined in the docstrings of the
properties.
"""
import re

def summarize_property(property):
    """Get a summary of property.

    Assuming that the property's docstring has the following format:

    "<info>: <constraint>, default=<default>.",

    this function extracts the constraint and the default value along with the
    property name. These can be used to generate XML file storing the options.

    Parameters
    ----------
    property : Property
        The property to summarize.

    Returns
    -------
    dict : {"name":name, "constraint":constraint, "default":default}
        The extracted summary.
    """
    parts = []
    for line in property.__doc__.strip().splitlines():
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
        "name": property.__name__,
        "constraint": match_regex.group("constraint").strip(),
        "default": match_regex.group("default").strip(),
    }


def parse_constraint(constraint: str):
    """Parse a human-readable constraint string.

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


def enforce_constraint(value, constraint: str):
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
    kind, constraints = parse_constraint(constraint)

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
