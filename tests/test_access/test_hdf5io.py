# pylint: disable=protected-access, no-member
"""Tests that cover the HDF5 IO operations performed by the tree structure."""
import os
from unittest.mock import MagicMock

import unyt
import numpy as np
import pytest

from a5py.data.access.hdf5io import (
    TreeFile, TreeFileManager, DataAccess, RESULTGROUP
    )
from a5py.data.access.tree import ROOT

from .conftest import (
    DATE, NOTE, INPUTS, OUTPUTS, CATEGORIES,
    FNTEST, FNEMPTY,
    )


@pytest.fixture()
def treefile():
    """Open a file for reading and writing.

    The file is closed and removed once a test is completed.
    """
    file = TreeFile(FNTEST, "a")
    yield file
    file.close()
    os.unlink(FNTEST)


@pytest.fixture()
def testfile():
    """Create a file with one input and one output."""
    with TreeFile(FNTEST, "a") as h5:
        h5.set_node(CATEGORIES[0], f"{INPUTS[0]}_1")
        h5.set_node(RESULTGROUP, f"{OUTPUTS[0]}_1")
        h5.set_datagroup(
            CATEGORIES[0], f"{INPUTS[0]}_1", date=DATE, note=NOTE,
            )
        h5.set_datagroup(
            RESULTGROUP, f"{OUTPUTS[0]}_1", date=DATE, note=NOTE,
            **{CATEGORIES[0]:f"{INPUTS[0]}_1"},
        )
        h5[f"{CATEGORIES[0]}/{INPUTS[0]}_1"]
    yield FNTEST
    os.unlink(FNTEST)


@pytest.fixture()
def emptyfile():
    yield FNEMPTY
    try:
        os.unlink(FNEMPTY)
    except FileNotFoundError:
        pass


@pytest.fixture()
def testmanager(testfile):
    return TreeFileManager(
        ROOT, testfile, file_exists=True, input_categories=CATEGORIES,
        )


@pytest.fixture()
def emptymanager(emptyfile):
    return TreeFileManager(
        ROOT, emptyfile, file_exists=False, input_categories=CATEGORIES,
        )


@pytest.fixture()
def treefilemanager():
    """Prepare a mock of TreeFileManager that manages a (immutable) file with
    fixed datasets.
    """
    treefilemanager = MagicMock()

    def read_node(node):
        if node == ROOT:
            return (f"{OUTPUTS[0]}_1", [f"{OUTPUTS[0]}_1"], [OUTPUTS[0]])
        elif node == CATEGORIES[0]:
            return (f"{INPUTS[0]}_1", [f"{INPUTS[0]}_1"], [INPUTS[0]])
        return (None, [], [])
    treefilemanager.read_node.side_effect = read_node

    def read_input(qid, var, category):
        return (DATE, NOTE)
    treefilemanager.read_input.side_effect = read_input

    def read_output(qid, var):
        return (DATE, NOTE, [f"{INPUTS[0]}_1"])
    treefilemanager.read_output.side_effect = read_output
    return treefilemanager


def test_treefile_set_node(treefile):
    """Test TreeFile.set_node.

    1. Setting a node creates node if it doesn't exist.
    2. The active attribute is also created but it's value is empty string.
    3. Setting the active attribute updates the value.
    4. Setting existing node without specifying active attribute retains
       previous value.
    """
    treefile.set_node("node")
    assert "node" in treefile

    node = treefile["node"]
    assert node.attrs["active"] == np.bytes_("")

    treefile.set_node("node", active="value")
    assert node.attrs["active"] == np.bytes_("value")

    treefile.set_node("node")
    assert node.attrs["active"] == np.bytes_("value")


def test_treefile_get_node(treefile):
    """Test TreeFile.get_node.

    After setting the active attribute, it can be retrieved as a string.
    """
    treefile.set_node("node", active="value")
    value = treefile.get_node(node="node")
    assert value == "value"


def test_treefile_set_datagroup(treefile):
    """Test TreeFile.set_datagroup.

    1. Setting a data group with no attributes creates the group if it doesn't
       exist. The date and note attributes are created but they are empty
       strings.
    2. Setting the date and note attributes updates the values.
    3. Setting existing data group without specifying date and note retains
       the values.
    4. Note can be updated.
    5. Extra attributes can be set.
    """
    treefile.set_datagroup("/", "dataset")
    assert "dataset" in treefile
    attrs = treefile["dataset"].attrs
    assert attrs["date"] == np.bytes_("")
    assert attrs["note"] == np.bytes_("")

    treefile.set_datagroup("/", "dataset", date=DATE, note=NOTE)
    assert attrs["date"] == np.bytes_(DATE)
    assert attrs["note"] == np.bytes_(NOTE)

    treefile.set_datagroup("/", "dataset")
    assert attrs["date"] == np.bytes_(DATE)
    assert attrs["note"] == np.bytes_(NOTE)

    treefile.set_datagroup("/", "dataset", note="New note")
    assert attrs["date"] == np.bytes_(DATE)
    assert attrs["note"] == np.bytes_("New note")

    treefile.set_datagroup("/", "dataset", extra_attribute="value")
    assert "extra_attribute" in attrs
    assert attrs["extra_attribute"] == np.bytes_("value")


def test_treefile_get_datagroup(treefile):
    """Test TreeFile.get_datagroup.

    1. Test that correct date and note are returned when there are no extra
       attributes.
    2. Test that date and note are returned along with the extra attributes
       when they are present. Test that the extra attribute has a correct value.
    """
    treefile.set_node("node")
    treefile.set_datagroup("node", "dataset", date=DATE, note=NOTE)
    date, note, _ = treefile.get_datagroup("node", "dataset")
    assert date == DATE
    assert note == NOTE

    treefile.set_datagroup("node", "dataset", extra_attribute="value")
    date, note, extra = treefile.get_datagroup("node", "dataset")
    assert "extra_attribute" in extra
    assert date == DATE
    assert note == NOTE
    assert extra["extra_attribute"] == "value"


def test_treefile_write_dataset(treefile):
    """Test TreeFile.write_datasets."""
    data = np.array([1, 2], dtype=np.float64).reshape((2,1))
    qnt = np.array([1.0])*unyt.m
    treefile.write_dataset("/", "data", data, False)
    group = treefile["/"]
    assert np.all( np.isclose(group["data"][:], data) )
    assert "units" not in group["data"].attrs

    treefile.write_dataset("/", "qnt", qnt, False)
    assert np.all( np.isclose(group["qnt"][:], qnt.v) )
    assert "units" in group["qnt"].attrs
    assert group["qnt"].attrs["units"], np.bytes_("m")


def test_treefile_read_dataset(treefile):
    """Test TreeFile.read_datasets."""
    data = np.array([1, 2], dtype=np.float64).reshape((2,1))
    qnt = np.array([1.0])*unyt.m
    treefile.write_dataset("/", "data", data, False)
    treefile.write_dataset("/", "qnt", qnt, False)

    dataout = treefile.read_dataset("/", "data")
    assert np.all( np.isclose(dataout, data) )
    dataout = treefile.read_dataset("/", "qnt")
    assert np.all( np.isclose(dataout, qnt, atol=0*unyt.m) )


@pytest.mark.parametrize("isempty", [True, False])
def test_manager_init(emptyfile, testfile, isempty):
    """Test TreeFileManager initialization from a file.

    1. If file is assumed to exist when it doesn't, an error is raised or
       if file is assumed to not exist when it does, raise an error.
    2. Otherwise, the file has nodes for results and input categories present.
       For empty file each node has active attribute set to empty string. For
       non-empty file no data is overwritten (each node has active attribute
       set to the corresponding QID and output has reference to input).
    """
    file = emptyfile if isempty else testfile
    with pytest.raises((FileNotFoundError, FileExistsError)):
        TreeFileManager(
            root=ROOT, filename=file, file_exists=isempty,
            input_categories=CATEGORIES,
            )

    TreeFileManager(
        root=ROOT, filename=file, file_exists=not isempty,
        input_categories=CATEGORIES,
        )
    with TreeFile(file, "r") as h5:
        for node in (CATEGORIES + [RESULTGROUP]):
            assert node in h5
            if not isempty and node == CATEGORIES[0]:
                assert h5.get_node(node) == f"{INPUTS[0]}_1"
            elif not isempty and node == RESULTGROUP:
                assert h5.get_node(node) == f"{OUTPUTS[0]}_1"
                outputname = f"{OUTPUTS[0]}_1"
                _, _, inputqids = h5.get_datagroup(node, outputname)
                assert CATEGORIES[0] in inputqids
                assert inputqids[CATEGORIES[0]] == f"{INPUTS[0]}_1"
                assert not CATEGORIES[1] in inputqids
            else:
                assert h5.get_node(node) == ""


@pytest.mark.parametrize("category", (CATEGORIES[0], ROOT))
def test_manager_read_node(testmanager, category):
    """Test TreeFileManager.read_node.

    Both input and result nodes should return active field, and QIDs and
    variants belonging to it.
    """
    active, names = testmanager.read_node(category)
    assert active == f"{OUTPUTS[0]}_1" if category == ROOT else f"{INPUTS[0]}_1"
    assert names == [f"{OUTPUTS[0]}_1" if category == ROOT else f"{INPUTS[0]}_1"]


@pytest.mark.parametrize("category", (CATEGORIES[0], ROOT))
def test_manager_set_active(testmanager, category):
    """Test TreeFileManager.set_active.

    The TreeFileManager does not care about consistency, so we can change the
    active QID to whatever and check that it was changed.
    """
    name = f"{OUTPUTS[0]}_1" if category == "root" else f"{INPUTS[0]}_1"
    testmanager.set_active(category, name)
    active, _ = testmanager.read_node(category)
    assert active == name


def test_manager_read_input(testmanager):
    """Test TreeFileManager.read_input.

    This should return variant date and note.
    """
    date, note, _ = testmanager.read_variant(f"{INPUTS[0]}_1", CATEGORIES[0])
    assert date == DATE
    assert note == NOTE


def test_manager_read_output(testmanager):
    """Test TreeFileManager.read_output.

    This should return variant date, note, and a dictionary of inputs.
    """
    date, note, inputs = testmanager.read_variant(f"{OUTPUTS[0]}_1", ROOT)
    assert date == DATE
    assert note == NOTE
    assert inputs == {CATEGORIES[0]: f"{INPUTS[0]}_1"}


def test_manager_write_input(testmanager):
    """Test TreeFileManager.read_input.

    Use QID that does not exist yet in the file and read the data afterwards
    for verification.
    """
    name = INPUTS[0]+ "_1"
    testmanager.write_variant(name, (DATE, NOTE, None), CATEGORIES[0])
    date, note, _ = testmanager.read_variant(name, CATEGORIES[0])
    assert date == DATE
    assert note == NOTE


def test_manager_write_output(emptymanager):
    """Test TreeFileManager.write_output.

    Use QID that does not exist yet in the file and read the data afterwards
    for verification.
    """
    name = OUTPUTS[0]+ "_1"
    emptymanager.write_variant(
        name, (DATE, NOTE, {CATEGORIES[0]:f"{INPUTS[0]}_1"}), ROOT
        )
    date, note, inputs = emptymanager.read_variant(name, ROOT)
    assert date == DATE
    assert note == NOTE
    assert inputs == {CATEGORIES[0]: f"{INPUTS[0]}_1"}


@pytest.mark.parametrize("isroot", (True, False))
def test_manager_set_note(testmanager, isroot):
    """Test TreeFileManager.set_note."""
    newnote = "New note"
    if isroot:
        testmanager.set_note(f"{OUTPUTS[0]}_1", ROOT, newnote)
        _, note, _ = testmanager.read_variant(f"{OUTPUTS[0]}_1", ROOT)
    else:
        testmanager.set_note(f"{INPUTS[0]}_1", CATEGORIES[0], newnote)
        _, note, _ = testmanager.read_variant(f"{INPUTS[0]}_1", CATEGORIES[0])
    assert note == newnote


@pytest.mark.parametrize(
        "variant, category",
        [(INPUTS[0], CATEGORIES[0]), (OUTPUTS[0], RESULTGROUP)],
        ids=["input", "output"]
        )
def test_manager_remove_data(testfile, testmanager, variant, category):
    """Test TreeFileManager.remove_data."""
    name = f"{variant}_1"
    testmanager.remove_variant(name, category)
    with TreeFile(testfile, "r") as h5:
        assert name not in h5[category]


def test_manager_repack(testfile, testmanager):
    """Test TreeFileManager.repack.

    Add some data to a group, then remove the group and compare the size
    a) before removal, b) after removal but before repacking, and c) after
    repacking.
    """
    with TreeFile(testfile, "a") as h5:
        name = f"{INPUTS[0]}_1"
        h5[CATEGORIES[0]][name].create_dataset("data", data=np.ones(10000,))

    originalsize = os.path.getsize(testfile)
    testmanager.remove_variant(f"{INPUTS[0]}_1", CATEGORIES[0])
    unpackedsize = os.path.getsize(testfile)
    testmanager.repack()
    repackedsize = os.path.getsize(testfile)

    assert unpackedsize <= originalsize
    assert repackedsize < unpackedsize
