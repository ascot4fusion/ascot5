import os
import tempfile
import pytest

from a5py import Ascot

@pytest.fixture
def ascot(request):
    file = getattr(request, "param", False)
    if file:
        fd, filepath = tempfile.mkstemp()
        os.close(fd)
        try:
            yield Ascot(filepath)
        finally:
            os.remove(filepath)
    else:
        yield Ascot()
