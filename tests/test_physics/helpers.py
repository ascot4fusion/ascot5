import numpy as np

def assert_isclose(name, result, expected, tolerance, relative=False):
    """Asserts that the result is close to the expected value and prints
    informative message on failure.

    Parameters
    ----------
    name : str
        Name of the quantity (for the error message).
    result : float
        Value to check.
    expected : float
        Expected value.
    tolerance : float
        Tolerance for the quantity.
    relative : bool, optional
        If True, the tolerance is relative to the expected value.

        By default the tolerance is absolute.
    """
    __tracebackhide__ = True
    rtol = tolerance if relative else 0
    atol = tolerance if not relative else 0
    difference = np.abs(result - expected)
    if relative:
        difference /= np.abs(expected)
    assert np.isclose(
        result,
        expected,
        rtol=rtol,
        atol=atol,
    ), (
        f"{name} mismatch:\n"
        f"  numerical:  {result:.6e}\n"
        f"  expected: {expected:.6e}\n"
        f"  difference: {difference:.6e}\n"
        f"  tolerance:  {tolerance:.2e}"
    )