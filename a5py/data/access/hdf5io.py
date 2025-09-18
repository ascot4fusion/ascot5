"""Module for creating and modifying the HDF5 file and for accessing meta data.

The file has the following structure: ::

    /
    ├── <input category>
    │   ├── <input variant name>
    |   |   ├── <dataset>
    │   │   └── <additional datasets>
    │   └── <additional input variants>
    ├── <additional input categories>
    └── results
        ├── <output variant name>
        │   ├── <diagnostic>
        │   │   ├── <dataset>
        │   │   └── <additional datasets>
        │   └── <additional diagnostics>
        └── <additional output variants>

Groups representing input categories and results have fixed names and they are
always present in the file even though there would be no data within them. Each
data group (input and output variants) are named as `<variant>_<id>`. Only input
variants and subgroups within output variants (diagnostics) contain actual data.

Variant groups contain following metadata:

- `date`: Date when the variant was created. This is stored as a string in
  format .
- `note`: User defined note describing the data.
- Output variants also contain attributes `<input category>` for each input
  category used in the simulation. The value of the attribute is the name of
  the variant that was used.

Every group corresponding to results or any input category has an attribute
`active` which is the name of the active variant within that group or empty
string if the group has no variants.
"""
from __future__ import annotations

import os
import tempfile
from typing import Optional

import h5py
import unyt
import numpy as np


RESULTGROUP = "results"
"""Name of the results group as it appears in HDF5 file."""


class TreeFile(h5py.File):
    """Low level interface to the HDF5 file.

    Extends :class:`h5py.File` with additional methods for reading and writing
    data assuming that the file stores the tree structure.
    """

    def __init__(self, filename: str, mode: str) -> None:
        """Open the file for reading or writing.

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
        node : str, *optional*
            Name of the node.
        active : str, *optional*
            The value of the active attribute.
        """
        if node not in self:
            self.create_group(node)
        group = self[node]

        group.attrs["active"] = (
            np.bytes_(active) if active is not None
            else group.attrs.get("active", np.bytes_(""))
            )

    def get_node(self, node: str) -> str:
        """Get metadata from a node.

        Parameters
        ----------
        node : str, *optional*
            Name of the node or None for the root.

        Returns
        -------
        active : str
            The value of the active attribute.
        """
        return self[node].attrs["active"].decode("utf-8")

    def set_datagroup(
            self, node: str,
            name: str,
            date: Optional[str]=None,
            note: Optional[str]=None,
            **additional_attrs: dict[str, str],
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
        date : str, *optional*
            Date when the data was created.
        note : str, *optional*
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
        group.attrs["note"] = (
            np.bytes_(note) if note else group.attrs.get("note", np.bytes_(""))
            )
        for attr, value in additional_attrs.items():
            group.attrs[attr] = (
                np.bytes_(value) if value
                else group.attrs.get(attr, np.bytes_(""))
            )

    def get_datagroup(
            self, node: str, name: str,
            ) -> tuple[str, str, dict[str, str]]:
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
            Names and values of any additional attributes.

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
            elif attr == "note":
                note = value.decode("utf-8")
            else:
                additional_attrs[attr] = value.decode("utf-8")
        if date is None or note is None:
            raise IOError("Failed to read date and/or note from the file.")

        return date, note, additional_attrs

    def write_dataset(
            self, path: str,
            name: str,
            dataset: np.ndarray | unyt.unyt_array,
            compress: bool,
            ) -> None:
        """Write a dataset.

        Parameters
        ----------
        path : str
            Path to the group where the datasets will be stored.

            The group should exist.
        name : str
            Name of the dataset on file.
        dataset : np.ndarray or unyt.array
            The value to be written.

            If the dataset is a unyt.array, the units will be stored
            as `units` attribute on the dataset.
        compress : bool
            Whether to compress the dataset.
        """
        group = self[path]
        if compress:
            dset = group.create_dataset(
                name, data=dataset, compression="gzip", compression_opts=9
                )
        else:
            dset = group.create_dataset(name, data=dataset)
        if hasattr(dataset, "units"):
            dset.attrs["units"] = np.bytes_(dataset.units)

    def read_dataset(
            self, path: str, name: Optional[str]=None,
            ) -> (dict[str, np.ndarray | unyt.unyt_array] |
                  np.ndarray | unyt.unyt_array):
        """Read datasets from a group.

        Parameters
        ----------
        path : str
            Path to the group where the datasets belong to.
        name : str, *optional*
            Instead of reading all datasets, read only this one.

        Returns
        -------
        data : Dict[str, np.ndarray | unyt.array]
            Datasets that were read.
        """
        dataset = self[path][name]
        try:
            units = unyt.unyt_quantity.from_string(
                dataset.attrs["units"].decode("utf-8")
            )
        except KeyError:
            units = 1
        if dataset.shape == ():
            return dataset[()] * units
        return dataset[:] * units


class DataAccess():
    """Provides access to datasets within a group.

    Attributes
    ----------
    filename : str
        Name of the file.
    path : str
        Path to the group where data is accessed.
    """

    def __init__(self, filename: str, path: str) -> None:
        self.filename: str = filename
        self.path: str = path

    @property
    def children(self) -> list[str]:
        """Names of groups belonging to this group."""
        with TreeFile(self.filename, "r") as h5:
            return list(h5[self.path].keys())

    def access_data(self, group: str) -> DataAccess:
        """Access data within a group that is current group's child.

        Parameters
        ----------
        group : str
            Name of the group within the current group.
        """
        with TreeFile(self.filename, "a") as h5:
            try:
                h5[self.path].create_group(group)
            except ValueError:
                pass
        return DataAccess(self.filename, f"{self.path}/{group}")

    def add_child(self, name: str) -> None:
        """Add a group to this group.

        Parameters
        ----------
        name : str
            Name of the group to be added.
        """
        with TreeFile(self.filename, "a") as h5:
            h5[self.path].create_group(name)

    def write(
            self, name: str,
            value: np.ndarray | unyt.unyt_array,
            compress: bool=True,
            ) -> None:
        """Write a dataset.

        Parameters
        ----------
        name : str
            Name of the dataset as it appears in the file.
        value : np.ndarray or unyt.array
            Data to be written.
        """
        with TreeFile(self.filename, "a") as h5:
            h5.write_dataset(self.path, name, value, compress)

    def read(self, name: str) -> np.ndarray | unyt.unyt_array:
        """Read a datasets

        Parameters
        ----------
        name : str, *optional*
            Name of the dataset as it appears in the file.

        Returns
        -------
        value : np.ndarray or unyt.array
            Data that was read.
        """
        with TreeFile(self.filename, "r") as h5:
            return h5.read_dataset(self.path, name)


class TreeFileManager():
    """Manages the structure of the HDF5 file that stores the data in a
    tree-like format.

    During initialization all nodes are created with empty meta data.
    """

    def __init__(
            self, root: str,
            filename: str,
            file_exists: bool,
            input_categories: list[str],
            ) -> None:
        """
        Parameters
        ----------
        root : str
            Name of the root node.
        filename : str
            Name of the HDF5 file.
        file_exists : bool
            Assume that the file exists, otherwise try to create a new file.
        input_categories : list[str]
            Names of the input categories.

        Raises
        ------
        FileNotFoundError
            If the file is assumed to exists when it does not exist.
        FileExistsError
            If the file is assumed to not exist when it does exist.
        """
        self.filename: str = filename
        """Name of the file."""

        self.root: str = root
        """Name of the root node.

        This node is stored as :attr:`RESULTGROUP` in the file.
        """

        if file_exists:
            if not os.path.isfile(filename):
                raise FileNotFoundError(
                    f"The file '{filename}' does not exist."
                    )
            with h5py.File(filename, "a"):
                # Just check that the file is HDF5 file.
                pass
        else:
            if os.path.isfile(filename):
                raise FileExistsError(
                    f"The file '{filename}' already exists."
                    )
            with h5py.File(filename, "w"):
                # Creates the file
                pass

        with TreeFile(self.filename, "a") as h5:
            h5.set_node(node=RESULTGROUP)
            for category in input_categories:
                h5.set_node(node=category)

    def access_data(
            self, name: str, node: str, subpath: Optional[str]=None
            ) -> DataAccess:
        """Access data within given variant's group.

        Parameters
        ----------
        name : str
            Name of the data.
        node : str
            Name of the node.
        subpath : str, *optional*
            Path to the child group if the data is stored there instead.
        """
        node = RESULTGROUP if node == self.root else node
        path = f"/{node}/{name}"
        if subpath:
            path = f"{path}/{subpath}"
        return DataAccess(self.filename, path)

    def set_active(self, node: str, active: str) -> None:
        """Store name of the active data.

        Parameters
        ----------
        node : str
            Name of the node.
        active : str, *optional*
            Name of the active data.
        """
        node = RESULTGROUP if node == self.root else node
        with TreeFile(self.filename, "a") as h5:
            h5.set_node(node, active)

    def set_note(self, name: str, node: str, note: str) -> None:
        """Set the note attribute for a given variant.

        Parameters
        ----------
        name : str
            Name of the data.
        node : str
            Name of the node where the data belongs.
        note : str
            The new note.
        """
        node = RESULTGROUP if node == self.root else node
        with TreeFile(self.filename, "a") as h5:
            h5.set_datagroup(node, name, note=note)

    def read_variant(
            self, name: str,
            node: str,
            ) -> tuple[str, str, dict[str, str] | None]:
        """Return variant's metadata.

        Parameters
        ----------
        name : str
            Name of the variant to be read.
        node : str
            Name of the node.

        Returns
        -------
        date : str
            Date when the data was created.
        note : str
            User-defined note on the data.
        inputs : Dict[str, str] | None
            Names of the inputs used by the variant or None for input variants.

        Raises
        ------
        IOError
            If failed to read metadata from the file.
        """
        node = RESULTGROUP if node == self.root else node
        with TreeFile(self.filename, "r") as h5:
            metadata = h5.get_datagroup(node, name)

        if node != RESULTGROUP:
            date, note, _ = metadata
            return date, note, _

        if len(metadata) < 3:
            raise IOError(
                "Could not read output metadata: the file does not contain "
                "references to inputs."
            )
        date, note, inputs = metadata
        return date, note, inputs

    def write_variant(
            self, name: str,
            metadata: tuple[str, str, dict[str, str] | None],
            node: str,
            ) -> None:
        """Create group for storing variant's metadata.

        Parameters
        ----------
        name : str
            Name of the data.
        metadata : tuple[str, str, dict[str, str] | None]
            Date when the data was created, user-defined note, and, for output
            variants, any used input category along with the name of the
            corresponding input variant.
        node : str, *optional*
            Node where the data belongs.

            Mandatory for input variants.
        """
        node = RESULTGROUP if node == self.root else node
        date, note, inputs = metadata
        with TreeFile(self.filename, "a") as h5:
            if inputs is None:
                inputs = {}
            h5.set_datagroup(
                node=node, name=name, date=date, note=note, **inputs,
                )
            if h5.get_node(node) == "":
                h5.set_node(node, name)

    def remove_variant(self, name: str, node: str) -> None:
        """Remove data from the file.

        Does not repack the file.

        Parameters
        ----------
        name : str
            Name of the data.
        node : str
            Name of the node where the data belongs.
        """
        node = RESULTGROUP if node == self.root else node
        with TreeFile(self.filename, "a") as h5:
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

    def read_node(self, node: str) -> tuple[str, list[str]]:
        """Return name of the active data and all names of data that belong to
        a given node.

        Parameters
        ----------
        node : str
            Name of the node.

        Returns
        -------
        active : str
            Name of the active data.
        names : List[str]
            Names of the data belonging to this node.
        """
        node = RESULTGROUP if node == self.root else node
        with TreeFile(self.filename, "r") as h5:
            active = h5.get_node(node)

            names = []
            for datagroup in h5[node]:
                name = datagroup.split("/")[-1]
                split = name.split("_")
                is_variant = (
                    len(split) == 2 and
                    len(split[0]) > 0 and
                    len(split[1]) > 0
                )
                if is_variant:
                    names.append(name)

        return active, names
