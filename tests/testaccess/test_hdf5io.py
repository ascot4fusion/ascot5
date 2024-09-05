# pylint: disable=protected-access, no-member
"""Tests that cover the HDF5 IO operations performed by the tree structure."""
import os
from unittest.mock import MagicMock

import unyt
import numpy as np
import pytest

import a5py.data.access.metadata as metadata
from a5py.data.access.metadata import MetaData
from a5py.data.access.hdf5 import HDF5Interface, HDF5Manager
from a5py.data.access.tree import Tree
from a5py import AscotIOException

from .conftest import (
    QID1, QID2, QID3, QID4, DATE, NOTE, INPUTVAR, OUTPUTVAR, CATEGORY,
    FNTEST, FNEMPTY,
    )


@pytest.fixture
def interface():
    "Create temporary file."
    interface = HDF5Interface(FNTEST)
    yield interface
    interface.close()
    os.unlink(FNTEST)


@pytest.fixture
def testfile():
    """Write a file with one input and one output."""
    with HDF5Interface(FNTEST) as h5:
        h5.set_node(CATEGORY, QID1)
        h5.set_node("results", QID2)
        h5.set_datagroup(
            CATEGORY, f"{INPUTVAR}_{QID1}", date=DATE, note=NOTE,
            )
        h5.set_datagroup(
            "results", f"{OUTPUTVAR}_{QID2}", date=DATE, note=NOTE,
            **{CATEGORY:QID1},
        )
        g = h5[f"{CATEGORY}/{INPUTVAR}_{QID1}"]
    yield FNTEST
    os.unlink(FNTEST)


@pytest.fixture
def emptyfile():
    yield FNEMPTY
    try:
        os.unlink(FNEMPTY)
    except FileNotFoundError:
        pass


@pytest.fixture
def hdf5manager():
    """Prepare a mock of HDF5Manager that manages a (immutable) file with
    fixed datasets."""
    hdf5manager = MagicMock()

    def read_node(node):
        if node == "root":
            return (QID2, [QID2], [OUTPUTVAR])
        elif node == CATEGORY:
            return (QID1, [QID1], [INPUTVAR])
        return (None, [], [])
    hdf5manager.read_node.side_effect = read_node

    def read_input(qid, var):
        return (DATE, NOTE)
    hdf5manager.read_input.side_effect = read_input

    def read_output(qid, var):
        return (DATE, NOTE, [QID1])
    hdf5manager.read_output.side_effect = read_output
    return hdf5manager


@pytest.fixture
def tree():
    """Create an empty tree."""
    return Tree()


def test_interface_set_node(interface):
    """Test setting a node with an active attribute."""
    interface.set_node("node")
    assert "node" in interface

    node = interface["node"]
    assert node.attrs["active"] == np.bytes_("")

    interface.set_node("node", active="value")
    assert node.attrs["active"] == np.bytes_("value")

    interface.set_node("node")
    assert node.attrs["active"] == np.bytes_("value")


def test_interface_get_node(interface):
    """Test getting the active attribute of a node."""
    interface.set_node("node", active="value")
    value = interface.get_node(node="node")
    assert value == "value"


def test_interface_set_datagroup(interface):
    """Test setting metadata to a dataset."""
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
    """Test getting metadata from a data group."""
    interface.set_node("node")
    interface.set_datagroup("node", "dataset")
    interface.set_datagroup("node", "dataset", date=DATE, note=NOTE)
    date, note = interface.get_datagroup("node", "dataset")
    assert date == DATE
    assert note == NOTE

    interface.set_datagroup(
        "node", "dataset", date=DATE, note=NOTE, extra_attribute="value",
        )
    _, _, extra = interface.get_datagroup("node", "dataset")
    assert "extra_attribute" in extra
    assert extra["extra_attribute"] == "value"


def test_interface_write_datasets(interface):
    """Test writing actual data."""
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
    """Test reading actual data."""
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


def test_manager_init(testfile, emptyfile):
    """Test initialization from a file with data."""
    with pytest.raises(FileExistsError):
        HDF5Manager(testfile, file_exists=False)

    HDF5Manager(testfile, file_exists=True)
    with HDF5Interface(testfile) as h5:
        for node in list(metadata.input_categories) + ["results"]:
            assert node in h5
            if node == CATEGORY:
                assert h5.get_node(node) == QID1
            elif node == "results":
                assert h5.get_node(node) == QID2
            else:
                assert h5.get_node(node) == ""


def test_manager_init_empty(testfile, emptyfile):
    """Test initialization from an empty file."""
    with pytest.raises(FileNotFoundError):
        HDF5Manager(emptyfile, file_exists=True)

    HDF5Manager(emptyfile, file_exists=False)
    with HDF5Interface(emptyfile) as h5:
        for node in list(metadata.input_categories) + ["results"]:
            assert node in h5
            assert h5.get_node(node) == ""


def test_manager_read_node(testfile):
    """Test reading metadata field from a node."""
    manager = HDF5Manager(testfile, file_exists=True)
    active, qids, variants = manager.read_node(CATEGORY)
    assert active == QID1
    assert qids == [QID1]
    assert variants == [INPUTVAR]

    active, qids, variants = manager.read_node("root")
    assert active == QID2
    assert qids == [QID2]
    assert variants == [OUTPUTVAR]


def test_manager_set_active(testfile):
    """Test setting active QID to a node."""
    manager = HDF5Manager(testfile, file_exists=True)

    manager.set_active("root", QID1)
    active, _, _ = manager.read_node("root")
    assert active == QID1

    manager.set_active(CATEGORY, QID2)
    active, _, _ = manager.read_node(CATEGORY)
    assert active == QID2


def test_manager_read_input(testfile):
    """Test reading input variant."""
    manager = HDF5Manager(testfile, file_exists=True)
    date, note = manager.read_input(QID1, INPUTVAR)
    assert date == DATE
    assert note == NOTE


def test_manager_read_run(testfile):
    """Test reading run variant."""
    manager = HDF5Manager(testfile, file_exists=True)
    out = manager.read_run(QID2, OUTPUTVAR)
    date, note, inputqids = out
    assert date == DATE
    assert note == NOTE
    assert inputqids == [QID1]


def test_manager_write_input(emptyfile):
    """Test writing input variant."""
    manager = HDF5Manager(emptyfile, file_exists=False)
    meta = metadata.MetaData(
        qid=QID1, date=DATE, note=NOTE,
        variant=INPUTVAR,
        )
    manager.write_input(meta)
    date, note = manager.read_input(meta.qid, meta.variant)
    assert date == DATE
    assert note == NOTE


def test_manager_write_run(emptyfile):
    """Test writing run variant."""
    manager = HDF5Manager(emptyfile, file_exists=False)
    meta = metadata.MetaData(
        qid=QID2, date=DATE, note=NOTE,
        variant=OUTPUTVAR,
        )
    manager.write_run(meta, {CATEGORY:QID1})
    out = manager.read_run(meta.qid, meta.variant)
    date, note, inputqids = out
    assert date == DATE
    assert note == NOTE
    assert inputqids == [QID1]


def test_manager_set_note(testfile):
    """Test modifying the note on the file."""
    newnote = "New note"
    manager = HDF5Manager(testfile, file_exists=True)

    manager.set_note(QID1, INPUTVAR, newnote)
    _, note = manager.read_input(QID1, INPUTVAR)
    assert note == newnote

    manager.set_note(QID2, OUTPUTVAR, newnote)
    _, note, _ = manager.read_run(QID2, OUTPUTVAR)
    assert note == newnote


def test_manager_remove_data(testfile):
    """Test removing a data variant."""
    manager = HDF5Manager(testfile, file_exists=True)

    manager.remove_variant(QID1, INPUTVAR)
    with HDF5Interface(testfile) as h5:
        name = f"{INPUTVAR}_{QID1}"
        assert name not in h5[CATEGORY]

    manager.remove_variant(QID2, OUTPUTVAR)
    with HDF5Interface(testfile) as h5:
        name = f"{OUTPUTVAR}_{QID2}"
        assert name not in h5["results"]


def test_manager_repack(testfile):
    """Test repacking the file."""
    manager = HDF5Manager(testfile, file_exists=True)
    with HDF5Interface(testfile) as h5:
        name = f"{INPUTVAR}_{QID1}"
        h5[CATEGORY][name].create_dataset(
            "data", data=np.ones(10000,),
            )

    originalsize = os.path.getsize(testfile)
    manager.remove_variant(QID1, INPUTVAR)
    unpackedsize = os.path.getsize(testfile)
    manager.repack()
    repackedsize = os.path.getsize(testfile)

    assert unpackedsize <= originalsize
    assert repackedsize < unpackedsize


def test_manager_read_input_variant(testfile):
    """Test reading input variant."""
    manager = HDF5Manager(testfile, file_exists=True)
    manager.write_input
    meta = metadata.MetaData(
        qid=QID1, date=DATE, note=NOTE,
        variant=INPUTVAR,
        )
    manager.write_input(meta)


def test_write_input_variant(emptyfile):
    """Test writing input variant."""
    manager = HDF5Manager(emptyfile, file_exists=False)
    manager.write_input
    meta = metadata.MetaData(
        qid=QID1, date=DATE, note=NOTE,
        variant=INPUTVAR,
        )
    manager.write_input(meta)


def test_initialization(tree, hdf5manager):
    """Test tree initialization from HDF5 file."""
    tree._treemanager.hdf5manager = hdf5manager
    tree._treemanager.init_from_hdf5()

    assert tree[CATEGORY].active.qid == QID1
    assert tree[CATEGORY].active.date == DATE
    assert tree[CATEGORY].active.variant == INPUTVAR
    assert tree[CATEGORY].active.note == NOTE

    with pytest.raises(AscotIOException):
        _ = tree.efield.active

    assert tree.active.qid == QID2
    assert tree.active.date == DATE
    assert tree.active.variant == OUTPUTVAR
    assert tree.active.note == NOTE
    assert tree.active[CATEGORY].qid == QID1

    with pytest.raises(AscotIOException):
        _ = tree.active.efield


def test_add_data_no_manager(tree):
    """Test adding input or output when there is no hdf5 manager."""
    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
    with pytest.raises(AscotIOException) as e:
        tree._treemanager.enter_input(meta, store_hdf5=True)
    assert "No HDF5 file was provided." in str(e.value)
    assert tree._treemanager.get_leaf(QID1).qid == QID1

    meta = MetaData(qid=QID3, date=DATE, note=NOTE, variant=OUTPUTVAR)
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_run(meta, [], [QID1], store_hdf5=True)
    assert "No HDF5 file was provided." in str(e.value)
    assert tree._treemanager.get_leaf(QID3).qid == QID3


def test_add_data(tree, hdf5manager):
    """Test adding data on tree and file."""
    tree._treemanager.hdf5manager = hdf5manager
    tree._treemanager.init_from_hdf5()

    meta = MetaData(qid=QID3, date=DATE, note=NOTE, variant=INPUTVAR)
    tree._treemanager.enter_input(meta)
    hdf5manager.write_input.assert_called_once_with(meta)

    meta = MetaData(qid=QID4, date=DATE, note=NOTE, variant=OUTPUTVAR)
    tree._treemanager.enter_run(meta, [], [QID1])
    hdf5manager.write_run.assert_called_once_with(meta, {CATEGORY:QID1})


def test_add_data_failure(tree, hdf5manager):
    """Test that no data is written if it fails to be included in the tree."""
    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
    try:
        tree._treemanager.enter_input(meta)
    except:
        pass
    hdf5manager.write_input.assert_not_called()

    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=OUTPUTVAR)
    try:
        tree._treemanager.enter_run(meta, [], [])
    except:
        pass
    hdf5manager.write_run.assert_not_called()


@pytest.mark.parametrize("repack", [False, True])
@pytest.mark.parametrize("method", ["leaf", "node", "tree"])
def test_destroy_data(tree, hdf5manager, repack, method):
    """Test removing data."""
    tree._treemanager.hdf5manager = hdf5manager
    tree._treemanager.init_from_hdf5()

    match method:
        case "leaf":
            tree[f"q{QID2}"].destroy(repack=repack)
        case "node":
            tree.destroy(repack=repack)
        case "tree":
            tree.destroy(QID2, repack=repack)
    hdf5manager.remove_variant.assert_called_once_with(QID2, OUTPUTVAR)
    if repack:
        hdf5manager.repack.assert_called_once()
    hdf5manager.reset_mock()

    match method:
        case "leaf":
            tree[CATEGORY][f"q{QID1}"].destroy(repack=repack)
        case "node":
            tree[CATEGORY].destroy(repack=repack)
        case "tree":
            tree.destroy(QID1, repack=repack)
    hdf5manager.remove_variant.assert_called_once_with(QID1, INPUTVAR)
    if repack:
        hdf5manager.repack.assert_called_once()


def test_destroy_input_inuse(tree, hdf5manager):
    """Test that input which is used cannot be removed from file."""
    tree._treemanager.hdf5manager = hdf5manager
    tree._treemanager.init_from_hdf5()
    try:
        tree.destroy(QID1)
    except AscotIOException:
        pass
    hdf5manager.remove_variant.assert_not_called()


def test_set_active(tree, hdf5manager):
    """Test that the active field is changed on disk."""
    tree._treemanager.hdf5manager = hdf5manager
    tree._treemanager.init_from_hdf5()
    hdf5manager.reset_mock()
    meta = MetaData(qid=QID3, date=DATE, note=NOTE, variant=OUTPUTVAR)
    leaf = tree._treemanager.enter_run(meta, [], [])
    leaf.activate()
    hdf5manager.set_active.assert_called_once_with("root", leaf.qid)


def test_set_active_not_on_file(tree, hdf5manager):
    """Test that we don't store active QID of data which is not on disk."""
    tree._treemanager.hdf5manager = hdf5manager
    tree._treemanager.init_from_hdf5()
    hdf5manager.reset_mock()
    meta = MetaData(qid=QID3, date=DATE, note=NOTE, variant=OUTPUTVAR)
    leaf = tree._treemanager.enter_run(meta, [], [], store_hdf5=False)
    leaf.activate()
    hdf5manager.set_active.assert_not_called()


def test_set_note(tree, hdf5manager):
    """Test that changing the note on data is reflected on disk."""
    tree._treemanager.hdf5manager = hdf5manager
    tree._treemanager.init_from_hdf5()
    tree.active.note = "Get to the choppa!"
    hdf5manager.set_note.assert_called_once_with(
        QID2, OUTPUTVAR, "Get to the choppa!"
        )
