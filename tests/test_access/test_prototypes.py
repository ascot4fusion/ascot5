import os
import pytest

import numpy as np
from .conftest import QIDS, INPUTS, OUTPUTS, DATES, NOTE, FNEMPTY
from a5py.exceptions import AscotIOException

from .prototypes import InputPrototype, OutputPrototype, TreePrototype


CATEGORY = "group"

@pytest.fixture()
def tree():
    """Create empty tree."""
    yield TreePrototype([CATEGORY], (FNEMPTY, False))
    if os.path.exists(FNEMPTY):
        os.unlink(FNEMPTY)


def test_input_create_write_read(tree):
    """"""
    data = np.array([0., 1., 2.])
    tree.create_inputprototype(bval=data)
    assert all(tree[CATEGORY].active.bval == data)
    assert not tree[CATEGORY].active.bval is data
    with pytest.raises(ValueError):
        tree[CATEGORY].active.bval[0] = 1

    tree = TreePrototype([CATEGORY], (FNEMPTY, True))
    assert all(tree[CATEGORY].active.bval == data)


def test_input_stage_unstage(tree):
    """"""
    data = np.array([0., 1., 2.])
    leaf = tree.create_inputprototype(bval=data, save=False)
    with pytest.raises(AscotIOException):
        leaf.unstage()
    with pytest.raises(AscotIOException):
        leaf.stage()
    leaf.save()
    leaf.unstage()
    assert all(leaf.bval == data)
    leaf.stage()
    assert all(leaf.bval == data)


def test_output_create_write_read(tree):
    """"""
    data = np.array([0., 1., 2.])
    leafin = tree.create_inputprototype(bval=data, save=False)
    leafout = OutputPrototype(QIDS[1], DATES[0], NOTE, "runprototype", {CATEGORY:leafin})
    tree._treemanager.enter_output(leafout, [leafin.qid], save=False)
    assert tree.active[CATEGORY].qid == leafin.qid

    params = {"diag1": True, "dval": data}
    leafout._setup(params)
    assert all(leafout.get_diag1() == data)

    leafin.save()
    tree.active.save()

    assert all(leafout.get_diag1() == data)
    tree = TreePrototype([CATEGORY], (FNEMPTY, True))
    assert all(tree.active.get_diag1() == data)
