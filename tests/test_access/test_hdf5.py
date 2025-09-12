# pylint: disable=protected-access, no-member
"""Tests that cover the HDF5 IO operations performed by the tree structure."""
import os
from unittest.mock import MagicMock

import unyt
import numpy as np
import pytest

from a5py.data.access.leaf import MetaData
from a5py.data.access.hdf5 import (
    HDF5Interface, HDF5Manager, HDF5MiniManager, RESULTGROUP
    )
from a5py.data.access.tree import ROOT

from .conftest import (
    QIDS, DATE, NOTE, INPUTS, OUTPUTS, CATEGORIES,
    FNTEST, FNEMPTY,
    )


@pytest.fixture()
def interface():
    """Open an interface to a file for reading and writing.

    The file is closed and removed once a test is completed.
    """
    interface = HDF5Interface(FNTEST, "a")
    yield interface
    interface.close()
    os.unlink(FNTEST)


@pytest.fixture()
def testfile():
    """Create a file with one input and one output."""
    with HDF5Interface(FNTEST, "a") as h5:
        h5.set_node(CATEGORIES[0], QIDS[0])
        h5.set_node(RESULTGROUP, QIDS[1])
        h5.set_datagroup(
            CATEGORIES[0], f"{INPUTS[0]}_{QIDS[0]}", date=DATE, note=NOTE,
            )
        h5.set_datagroup(
            RESULTGROUP, f"{OUTPUTS[0]}_{QIDS[1]}", date=DATE, note=NOTE,
            **{CATEGORIES[0]:QIDS[0]},
        )
        h5[f"{CATEGORIES[0]}/{INPUTS[0]}_{QIDS[0]}"]
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
    return HDF5Manager(
        ROOT, testfile, file_exists=True, input_categories=CATEGORIES,
        )


@pytest.fixture()
def emptymanager(emptyfile):
    return HDF5Manager(
        ROOT, emptyfile, file_exists=False, input_categories=CATEGORIES,
        )


@pytest.fixture()
def hdf5manager():
    """Prepare a mock of HDF5Manager that manages a (immutable) file with
    fixed datasets.
    """
    hdf5manager = MagicMock()

    def read_node(node):
        if node == ROOT:
            return (QIDS[1], [QIDS[1]], [OUTPUTS[0]])
        elif node == CATEGORIES[0]:
            return (QIDS[0], [QIDS[0]], [INPUTS[0]])
        return (None, [], [])
    hdf5manager.read_node.side_effect = read_node

    def read_input(qid, var, category):
        return (DATE, NOTE)
    hdf5manager.read_input.side_effect = read_input

    def read_output(qid, var):
        return (DATE, NOTE, [QIDS[0]])
    hdf5manager.read_output.side_effect = read_output
    return hdf5manager


def test_interface_set_node(interface):
    """Test HDF5Interface.set_node.

    1. Setting a node creates node if it doesn't exist.
    2. The active attribute is also created but it's value is empty string.
    3. Setting the active attribute updates the value.
    4. Setting existing node without specifying active attribute retains
       previous value.
    """
    interface.set_node("node")
    assert "node" in interface

    node = interface["node"]
    assert node.attrs["active"] == np.bytes_("")

    interface.set_node("node", active="value")
    assert node.attrs["active"] == np.bytes_("value")

    interface.set_node("node")
    assert node.attrs["active"] == np.bytes_("value")


def test_interface_get_node(interface):
    """Test HDF5Interface.get_node.

    After setting the active attribute, it can be retrieved as a string.
    """
    interface.set_node("node", active="value")
    value = interface.get_node(node="node")
    assert value == "value"


def test_interface_set_datagroup(interface):
    """Test HDF5Interface.set_datagroup.

    1. Setting a data group with no attributes creates the group if it doesn't
       exist. The date and note attributes are created but they are empty
       strings.
    2. Setting the date and note attributes updates the values.
    3. Setting existing data group without specifying date and note retains
       the values.
    4. Note can be updated.
    5. Extra attributes can be set.
    """
    interface.set_datagroup("/", "dataset")
    assert "dataset" in interface
    attrs = interface["dataset"].attrs
    assert attrs["date"] == np.bytes_("")
    assert attrs["desc"] == np.bytes_("")

    interface.set_datagroup("/", "dataset", date=DATE, note=NOTE)
    assert attrs["date"] == np.bytes_(DATE)
    assert attrs["desc"] == np.bytes_(NOTE)

    interface.set_datagroup("/", "dataset")
    assert attrs["date"] == np.bytes_(DATE)
    assert attrs["desc"] == np.bytes_(NOTE)

    interface.set_datagroup("/", "dataset", note="New note")
    assert attrs["date"] == np.bytes_(DATE)
    assert attrs["desc"] == np.bytes_("New note")

    interface.set_datagroup("/", "dataset", extra_attribute="value")
    assert "extra_attribute" in attrs
    assert attrs["extra_attribute"] == np.bytes_("value")


def test_interface_get_datagroup(interface):
    """Test HDF5Interface.get_datagroup.

    1. Test that correct date and note are returned when there are no extra
       attributes.
    2. Test that date and note are returned along with the extra attributes
       when they are present. Test that the extra attribute has a correct value.
    """
    interface.set_node("node")
    interface.set_datagroup("node", "dataset", date=DATE, note=NOTE)
    date, note = interface.get_datagroup("node", "dataset")
    assert date == DATE
    assert note == NOTE

    interface.set_datagroup("node", "dataset", extra_attribute="value")
    date, note, extra = interface.get_datagroup("node", "dataset")
    assert "extra_attribute" in extra
    assert date == DATE
    assert note == NOTE
    assert extra["extra_attribute"] == "value"


def test_interface_write_datasets(interface):
    """Test HDF5Interface.write_datasets."""
    data = {
        "data":np.array([1, 2], dtype=np.float64).reshape((2,1)),
        "qnt":np.array([1.0])*unyt.m,
    }
    interface.write_datasets("/", data)
    group = interface["/"]
    assert np.all( np.isclose(group["data"][:], data["data"]) )
    assert "units" not in group["data"].attrs
    assert np.all( np.isclose(group["qnt"][:], data["qnt"].v) )
    assert "units" in group["qnt"].attrs
    assert group["qnt"].attrs["units"], np.bytes_("m")


def test_interface_read_datasets(interface):
    """Test HDF5Interface.read_datasets."""
    data = {
        "data":np.array([1, 2], dtype=np.float64).reshape((2,1)),
        "qnt":np.array([1.0])*unyt.m,
    }
    interface.write_datasets("/", data)
    dataout = interface.read_datasets("/")
    assert np.all( np.isclose(dataout["data"], data["data"]) )
    assert np.all( np.isclose(dataout["qnt"], data["qnt"], atol=0*unyt.m) )
    dataout = interface.read_datasets("/", "data")
    assert np.all( np.isclose(dataout, data["data"]) )


@pytest.mark.parametrize("isempty", [True, False])
def test_manager_init(emptyfile, testfile, isempty):
    """Test HDF5Manager initialization from a file.

    1. If file is assumed to exist when it doesn't, an error is raised or
       if file is assumed to not exist when it does, raise an error.
    2. Otherwise, the file has nodes for results and input categories present.
       For empty file each node has active attribute set to empty string. For
       non-empty file no data is overwritten (each node has active attribute
       set to the corresponding QID and output has reference to input).
    """
    file = emptyfile if isempty else testfile
    with pytest.raises((FileNotFoundError, FileExistsError)):
        HDF5Manager(
            root=ROOT, filename=file, file_exists=isempty,
            input_categories=CATEGORIES,
            )

    HDF5Manager(
        root=ROOT, filename=file, file_exists=not isempty,
        input_categories=CATEGORIES,
        )
    with HDF5Interface(file, "r") as h5:
        for node in (CATEGORIES + [RESULTGROUP]):
            assert node in h5
            if not isempty and node == CATEGORIES[0]:
                assert h5.get_node(node) == QIDS[0]
            elif not isempty and node == RESULTGROUP:
                assert h5.get_node(node) == QIDS[1]
                outputname = f"{OUTPUTS[0]}_{QIDS[1]}"
                _, _, inputqids = h5.get_datagroup(node, outputname)
                assert CATEGORIES[0] in inputqids
                assert inputqids[CATEGORIES[0]] == QIDS[0]
                assert not CATEGORIES[1] in inputqids
            else:
                assert h5.get_node(node) == ""


@pytest.mark.parametrize("category", (CATEGORIES[0], ROOT))
def test_manager_read_node(testmanager, category):
    """Test HDF5Manager.read_node.

    Both input and result nodes should return active field, and QIDs and
    variants belonging to it.
    """
    active, qids, variants = testmanager.read_node(category)
    assert active == QIDS[1] if category == ROOT else QIDS[0]
    assert qids == [QIDS[1] if category == ROOT else QIDS[0]]
    assert variants == [OUTPUTS[0] if category == ROOT else INPUTS[0]]


@pytest.mark.parametrize("category", (CATEGORIES[0], ROOT))
def test_manager_set_active(testmanager, category):
    """Test HDF5Manager.set_active.

    The HDF5Manager does not care about consistency, so we can change the
    active QID to whatever and check that it was changed.
    """
    testmanager.set_active(category, QIDS[2])
    active, _, _ = testmanager.read_node(category)
    assert active == QIDS[2]


def test_manager_read_input(testmanager):
    """Test HDF5Manager.read_input.

    This should return variant date and note.
    """
    date, note = testmanager.read_input(QIDS[0], INPUTS[0], CATEGORIES[0])
    assert date == DATE
    assert note == NOTE


def test_manager_read_output(testmanager):
    """Test HDF5Manager.read_output.

    This should return variant date, note, and list of input QIDs.
    """
    date, note, inputqids = testmanager.read_output(QIDS[1], OUTPUTS[0])
    assert date == DATE
    assert note == NOTE
    assert inputqids == [QIDS[0]]


def test_manager_write_input(testmanager):
    """Test HDF5Manager.read_input.

    Use QID that does not exist yet in the file and read the data afterwards
    for verification.
    """
    meta = MetaData(qid=QIDS[3], date=DATE, note=NOTE, variant=INPUTS[0])
    testmanager.write_input(meta, CATEGORIES[0])
    date, note = testmanager.read_input(meta.qid, meta.variant, CATEGORIES[0])
    assert date == DATE
    assert note == NOTE


def test_manager_write_output(emptymanager):
    """Test HDF5Manager.write_output.

    Use QID that does not exist yet in the file and read the data afterwards
    for verification.
    """
    meta = MetaData(qid=QIDS[3], date=DATE, note=NOTE, variant=OUTPUTS[0])
    emptymanager.write_output(meta, {CATEGORIES[0]:QIDS[0]})
    out = emptymanager.read_output(meta.qid, meta.variant)
    date, note, inputqids = out
    assert date == DATE
    assert note == NOTE
    assert inputqids == [QIDS[0]]


@pytest.mark.parametrize("isroot", (True, False))
def test_manager_set_note(testmanager, isroot):
    """Test HDF5Manager.set_note."""
    newnote = "New note"
    if isroot:
        testmanager.set_note(QIDS[0], INPUTS[0], CATEGORIES[0], newnote)
        _, note = testmanager.read_input(QIDS[0], INPUTS[0], CATEGORIES[0])
    else:
        testmanager.set_note(QIDS[1], OUTPUTS[0], RESULTGROUP, newnote)
        _, note, _ = testmanager.read_output(QIDS[1], OUTPUTS[0])
    assert note == newnote


@pytest.mark.parametrize(
        "qid, variant, category",
        [(QIDS[0], INPUTS[0], CATEGORIES[0]),
         (QIDS[1], OUTPUTS[0], RESULTGROUP)],
        ids=["input", "output"]
        )
def test_manager_remove_data(testfile, testmanager, qid, variant, category):
    """Test HDF5Manager.remove_data."""
    testmanager.remove_variant(qid, variant, category)
    with HDF5Interface(testfile, "r") as h5:
        name = f"{variant}_{qid}"
        assert name not in h5[category]


def test_manager_repack(testfile, testmanager):
    """Test HDF5Manager.repack.

    Add some data to a group, then remove the group and compare the size
    a) before removal, b) after removal but before repacking, and c) after
    repacking.
    """
    with HDF5Interface(testfile, "a") as h5:
        name = f"{INPUTS[0]}_{QIDS[0]}"
        h5[CATEGORIES[0]][name].create_dataset(
            "data", data=np.ones(10000,),
            )

    originalsize = os.path.getsize(testfile)
    testmanager.remove_variant(QIDS[0], INPUTS[0], CATEGORIES[0])
    unpackedsize = os.path.getsize(testfile)
    testmanager.repack()
    repackedsize = os.path.getsize(testfile)

    assert unpackedsize <= originalsize
    assert repackedsize < unpackedsize
