"""Module for creating and modifying the HDF5 file and for accessing meta data.

The file has a following structure:

    /
    ├── bfield
    │   ├── B_TC_0123456789
    |   |   ├── bxyz
    │   │   └── Other *datasets*
    │   └── Other input *variants* in <variant>_<qid> format
    ├── Other input *categories*
    └── results
        ├── run_2345678901
        │   ├── inistate
        │   │   ├── weight
        │   │   └── Other *datasets*
        │   └── Other *diagnostics*
        └── Other run *variants* in <variant>_<qid> format

Input categories have fixed names and they are always present in the file even
though there would be no data within them. Each data group (input and output
variants) are named as <variant>_<qid>. They are the only wants (adn their
possible subgroups) that contain any datasets.

Attributes are used to store metadata within the groups:

- Each data group has `date` which is the date when the data was created and
  `desc` which is an user-specified description of the data.
- Groups representing results and input categories have `active` which is
  the QID of the currently active data within that group.
- Run variants have attributes named <input category> for each input category
  used in the simulation, and it contains QID of the input variant used.
"""
import os
import tempfile
from typing import Optional, TYPE_CHECKING

import h5py
import unyt
import numpy as np

from .leaf import QIDLEN, MetaData

if TYPE_CHECKING:
    from .leaf import MetaData


RESULTGROUP = "results"
"""Name of the results group as it appears in HDF5 file."""


class HDF5Interface(h5py.File):
    """Low level interface to HDF5.

    Helper class for the `HDF5Manager`. Does not traverse the tree but just
    creates and removes groups and sets the given groups attributes.
    """

    def __init__(self, filename, mode):
        """Initialize HDF5Interface.

        Parameters
        ----------
        filename : str
            HDF5 file name.
        mode : {"r", "a"}
            File open mode ("r"ead or "a"ppend).
        """
        super().__init__(filename, mode)

    def set_node(self, node: str, active: Optional[str]=None) -> None:
        """Set metadata on a given node, creating the group if none exists.

        Parameters
        ----------
        node : str, optional
            Name of the node.
        active : str, optional
            The value of the active attribute.
        """
        if node not in self:
            self.create_group(node)
        group = self[node]

        group.attrs["active"] = (
            np.bytes_(active) if not active is None
            else group.attrs.get("active", np.bytes_(""))
            )

    def get_node(self, node: str) -> str:
        """Get metadata from a node.

        Parameters
        ----------
        node : str, optional
            Name of the node or None for the root.

        Returns
        -------
        active : str
            The value of the active attribute.
        """
        return self[node].attrs["active"].decode("utf-8")

    def set_datagroup(
            self,
            node: str,
            name: str,
            date: Optional[str]=None,
            note: Optional[str]=None,
            **additional_attrs,
            ) -> None:
        """Set metadata on a data group, creating the group if none exists.

        If any attributes specified in the parameters are None, then empty
        string is used when creating the group or the value is not changed if
        the group exists already.

        Parameters
        ----------
        node : str
            Node to which the data group belongs to.
        name : str
            Name of the data group as it appears in the file.
        date : str, optional
            Date when the data was created.
        note : str, optional
            User-defined note of the data.
        **additional_attrs : dict[str, str]
            Any additional attributes to store.
        """
        if name not in self[node]:
            self[node].create_group(name)

        group = self[node][name]
        group.attrs["date"] = (
            np.bytes_(date) if date else group.attrs.get("date", np.bytes_(""))
            )
        group.attrs["desc"] = (
            np.bytes_(note) if note else group.attrs.get("desc", np.bytes_(""))
            )
        for attr, value in additional_attrs.items():
            group.attrs[attr] = (
                np.bytes_(value) if value
                else group.attrs.get(attr, np.bytes_(""))
            )

    def get_datagroup(
            self, node: str, name: str,
            ) -> tuple[str, str] | tuple[str, str, dict[str, str]]:
        """Get metadata from a data group.

        Parameters
        ----------
        node : str
            Node to which this data group belongs to.
        name : str
            Name of the data group as it appears in the file.

        Returns
        -------
        date : str
            Date when the data was created.
        note : str
            User-defined note of the data.
        additional_attrs : Dict[str, str]
            Any additional attributes.

        Raises
        ------
        IOError
            If failed to read date or note from the file.
        """
        group = self[node][name]
        additional_attrs = {}
        date, note = None, None
        for attr, value in group.attrs.items():
            if attr == "date":
                date = value.decode("utf-8")
            elif attr == "desc":
                note = value.decode("utf-8")
            else:
                additional_attrs[attr] = value.decode("utf-8")
        if date is None or note is None:
            raise IOError("Failed to read date and note from the file.")

        if not additional_attrs:
            return date, note
        return date, note, additional_attrs

    def write_datasets(
            self, path: str, data: dict[str, np.ndarray | unyt.unyt_array],
            ) -> None:
        """Write datasets to a group.

        If the data has physical quantities, those are stored in "units"
        attribute of the corresponding dataset.

        Parameters
        ----------
        path : str
            Path to the group where the datasets will be stored.
        data : dict
            Name of the data on file and the value to be written.
        """
        group = self[path]
        for field, value in data.items():
            dset = group.create_dataset(field, data=value)
            if hasattr(value, "units"):
                dset.attrs["units"] = np.bytes_(value.units)

    def read_datasets(
            self, path: str, name: Optional[str]=None,
            ) -> dict[str, np.ndarray | unyt.unyt_array]:
        """Read datasets from a group.

        Parameters
        ----------
        path : str
            Path to the group where the datasets belong to.
        name : str, optional
            Instead of reading all datasets, read only this one.

        Returns
        -------
        data : Dict[str, np.ndarray | unyt.unyt_array]
            Datasets that were read.
        """
        def read_dataset(dataset):
            """Read dataset with associated units."""
            try:
                units = unyt.unyt_quantity.from_string(
                    dataset.attrs["units"].decode("utf-8")
                )
            except KeyError:
                units = 1
            if dataset.shape == ():
                return dataset[()] * units
            return dataset[:] * units

        group = self[path]
        if name is None:
            data = {}
            for dataname, dataset in group.items():
                if isinstance(dataset, h5py.Dataset):
                    data[dataname] = read_dataset(dataset)
            return data

        dataset = group[name]
        return read_dataset(dataset)


class HDF5Manager():
    """Manages contents of the ASCOT5 HDF5 file.

    Attributes
    ----------
    filename : str
        Name of the file which is managed.
    """

    def __init__(
            self, root: str,
            filename: str,
            file_exists: bool,
            input_categories: str,
            ) -> None:
        """Initializes the managed file.

        1. Sets the active flag in the root to "" (if none exists).
        2. Creates all input category groups (for those that do not exist) and
           initializes the active flag to "" (for groups which were created).
        3. Sets input references to "" in simulation output groups (for unused
           inputs).

        Parameters
        ----------
        filename : str
            Name of the HDF5 file.
        file_exists : bool
            Assume that the file exists, otherwise create a new file.

        Raises
        ------
        Ascot5IOException
            If creating a file when one already exists or not creating a file
            when one does not exist.
        """
        self.filename = filename
        self.root = root
        if file_exists:
            if not os.path.isfile(filename):
                raise FileNotFoundError(
                    f"The file '{filename}' does not exist."
                    )
            with h5py.File(filename, "a"):
                pass
        else:
            if os.path.isfile(filename):
                raise FileExistsError(
                    f"The file '{filename}' already exists."
                    )
            with h5py.File(filename, "w"):
                pass

        with HDF5Interface(self.filename, "a") as h5:
            h5.set_node(node=RESULTGROUP)

            for category in input_categories:
                h5.set_node(node=category)

    def read_node(self, node: str) -> tuple[str, list[str], list[str]]:
        """Return list of QIDs, variants, and the active QID for a given node.

        Parameters
        ----------
        node : str
            Name of the node.

        Returns
        -------
        active : str
            QID of the active data.
        qids : List[str]
            List of QIDs of the data belonging to this node.
        variants : List[str]
            List of variants of the data belonging to this node in same order as
            `qids`.
        """
        node = RESULTGROUP if node == self.root else node
        qids, variants = [], []
        with HDF5Interface(self.filename, "r") as h5:
            active = h5.get_node(node)
            for datagroup in h5[node]:
                name = datagroup.split("/")[-1]

                split = name.split("_")
                is_variant = (
                    len(split) == 2 and
                    len(split[1]) == QIDLEN and
                    split[1].isdigit()
                )
                if is_variant:
                    variant, qid = split
                    qids.append(qid)
                    variants.append(variant)

        return active, qids, variants

    def set_active(self, node: str, active: str) -> None:
        """Store the active QID.

        Parameters
        ----------
        node : str
            Name of the node.
        active : str, optional
            QID of the active data.
        """
        node = RESULTGROUP if node == self.root else node
        with HDF5Interface(self.filename, "a") as h5:
            h5.set_node(node, active)

    def read_input(
            self, qid: str,
            variant: str,
            category: str,
            ) -> tuple[str, str]:
        """Return date and note of this input variant.

        Parameters
        ----------
        qid : str
            QID of the data.
        variant : str
            Data variant.

        Returns
        -------
        date : str
            Date when the data was created.
        note : str
            User-defined note on the data.
        """
        name = f"{variant}_{qid}"
        with HDF5Interface(self.filename, "r") as h5:
            return h5.get_datagroup(category, name)

    def write_input(self, meta: MetaData, category: str) -> None:
        """Create an input variant with the given metadata.

        If this is the first input in this node, the active flag will be set to
        the QID of the input.

        Parameters
        ----------
        meta : MetaData
            Input metadata.
        category : str
            Input category.
        """
        name = f"{meta.variant}_{meta.qid}"
        with HDF5Interface(self.filename, "a") as h5:
            h5.set_datagroup(
                node=category, name=name, date=meta.date, note=meta.note,
                )
            if h5.get_node(category) == "":
                h5.set_node(category, meta.qid)

    def read_output(self, qid: str, variant: str) -> tuple[str, str, list[str]]:
        """Return date, note, and list of input QIDs for this run variant.

        Parameters
        ----------
        qid : str
            QID of the run.
        variant : str
            Run variant.

        Returns
        -------
        date : str
            Date when the data was created.
        note : str
            User-defined note of the data.
        inputqids : List[str]
            List of QIDs for the inputs used to produce the output.
        """
        name = f"{variant}_{qid}"
        with HDF5Interface(self.filename, "r") as h5:
            date, note, usedinputs = h5.get_datagroup(RESULTGROUP, name)
            if not usedinputs:
                raise IOError(
                    "Could not read output metadata: the file does not contain "
                    "references to inputs."
                )
        return date, note, list(usedinputs.values())

    def write_output(self, meta: MetaData, usedinputs: dict[str, str]) -> None:
        """Create a run variant with the given metadata.

        Parameters
        ----------
        meta : MetaData
            Output metadata.
        usedinputs : Dict[str, str]
            Key-value pairs of input category and QID of the input used from
            that category in the simulation.
        """
        name = f"{meta.variant}_{meta.qid}"
        with HDF5Interface(self.filename, "a") as h5:
            h5.set_datagroup(
                node=RESULTGROUP, name=name, date=meta.date, note=meta.note,
                **usedinputs,
                )

    def set_note(self, qid: str, variant: str, node: str, note: str) -> None:
        """Set the note attribute in a data.

        Parameters
        ----------
        qid : str
            QID of the data.
        variant : str
            Data variant.
        note : str
            The new note.
        """
        name = f"{variant}_{qid}"
        node = RESULTGROUP if node == self.root else node
        with HDF5Interface(self.filename, "a") as h5:
            h5.set_datagroup(node, name, note=note)

    def remove_variant(self, qid: str, variant: str, node: str) -> None:
        """Remove data from the file.

        Does not repack the file.

        Parameters
        ----------
        qid : str
            QID of the data.
        variant : str
            Data variant.
        """
        name = f"{variant}_{qid}"
        node = RESULTGROUP if node == self.root else node
        with HDF5Interface(self.filename, "a") as h5:
            del h5[node][name]

    def repack(self) -> None:
        """Repack the file.

        Removing data from a HDF5 file using h5py only removes the references to
        it. To free memory, repacking the data is necessary. This method copies
        the data to a new HDF5 file and replaces the original with it.
        """
        permissions = os.stat(self.filename).st_mode
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            with (
                h5py.File(self.filename, "r") as oldfile,
                h5py.File(temp_file, "a") as newfile,
                ):
                for group in oldfile.keys():
                    oldfile.copy(group, newfile, group)
                for attr, value in oldfile.attrs.items():
                    newfile.attrs[attr] = value
                tempfilename = temp_file.name

        os.rename(tempfilename, self.filename)
        os.chmod(self.filename, permissions)

    def get_minimanager(self, node, variant, qid, subpath=None):
        node = RESULTGROUP if node == self.root else node
        path = f"/{node}/{variant}_{qid}"
        if subpath:
            path = f"{path}/{subpath}"
        return HDF5MiniManager(self.filename, path)


class HDF5MiniManager():
    """Provides access to a single HDF5 group only."""

    def __init__(self, filename, path):
        self.filename = filename
        self.path = path

    def get_minimanager(self, group):
        with HDF5Interface(self.filename, "a") as h5:
            try:
                h5[self.path].create_group(group)
            except ValueError:
                pass
        return HDF5MiniManager(self.filename, f"{self.path}/{group}")

    def get_children(self):
        """Return names of any groups belonging to this group."""
        with HDF5Interface(self.filename, "r") as h5:
            return list(h5[self.path].keys())

    def add_children(self, name):
        """Add a group to this group."""
        with HDF5Interface(self.filename, "a") as h5:
            h5.add_group(self.path, name)

    def write(
            self, name, value,
            ) -> None:
        """Write datasets to a data group.

        Parameters
        ----------
        data : dict[str, np.array or unyt.array]
            Name-value pairs of data to be written.
        """
        with HDF5Interface(self.filename, "a") as h5:
            h5.write_datasets(self.path, {name: value})

    def read(
            self,
            name: Optional[str] = None,
            ) -> dict[str, np.ndarray | unyt.unyt_array]:
        """Read the data from the HDF5 file.

        Parameters
        ----------
        qid : str
            QID of the data.
        variant : str
            Data variant.
        name : str, optional
            Instead of reading all datasets, read only this one.

        Returns
        -------
        data : dict[str, np.array or unyt.array] or np.array or unyt.array
            Name-value pairs of read data or just the value of the data if
            `name` was specified.
        """
        with HDF5Interface(self.filename, "r") as h5:
            return h5.read_datasets(self.path, name=name)
