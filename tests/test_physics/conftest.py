import pytest

from a5py import Ascot

def pytest_addoption(parser):
    parser.addoption(
        "--inspect",
        action="store_true",
        help="Load test data from HDF5 file",
    )

    parser.addoption(
        "--dump",
        action="store_true",
        help="Dump test data even on success",
    )

    parser.addoption(
        "--plot",
        action="store_true",
        help="Plot test results",
    )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    rep = outcome.get_result()
    setattr(item, f"rep_{rep.when}", rep)


@pytest.fixture
def inspect(request):
    return request.config.getoption("--inspect")


@pytest.fixture
def dump(request):
    return request.config.getoption("--dump")


@pytest.fixture
def plot(request):
    return request.config.getoption("--plot")

@pytest.fixture
def ascot(request, inspect, dump):

    fn = request.node.nodeid.split("::")[1] + ".h5"
    fn = fn.replace("[", "-").replace("]", "")
    out = fn if inspect else None
    ascot = Ascot(out)
    yield ascot

    failed = request.node.rep_call.failed
    if not inspect and (failed or dump):
        ascot.data.save(fn)
