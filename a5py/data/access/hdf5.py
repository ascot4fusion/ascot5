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
from __future__ import annotations

import os
import tempfile
from typing import List, Dict, Tuple, Optional, Union

import h5py
import unyt
import numpy as np

from . import metadata
from . metadata import QIDLEN, MetaData

TreeContent = Tuple[ List[str], List[str], List[str], List[ List[str] ] ]

class HDF5Interface(h5py.File):
    """Low level interface to HDF5.

    Helper class for the `HDF5Manager`. Does not traverse the tree but just
    creates and removes groups and sets the given groups attributes.
    """

    def __init__(self, filename):
        """Initialize HDF5Interface.

        Parameters
        ----------
        filename : str
            HDF5 file name.
        """
        super().__init__(filename, "a")

    def set_node(self, node: str, active: Optional[str] = None) -> None:
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
            np.bytes_(active) if active
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
            date: Optional[str] = None,
            note: Optional[str] = None,
            **additional_attrs,
            ) -> None:
        """Set metadata on a data group, creating the group if none exists.

        If any attributes specified in the parameters are None, the default
        value is used when creating the group or the value is not changed if
        the group exists.

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
            ) -> Union[Tuple[str, str], Tuple[str, str, Dict[str, str]]]:
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
        additional_attrs** : Dict[str, str]
            Any additional attributes.
        """
        group = self[node][name]
        additional_attrs = {}
        for attr, value in group.attrs.items():
            if attr == "date":
                date = value.decode("utf-8")
            elif attr == "desc":
                note = value.decode("utf-8")
            else:
                additional_attrs[attr] = value.decode("utf-8")
        if not additional_attrs:
            return date, note
        return date, note, additional_attrs

    def write_datasets(
            self, path: str, data: Union[np.array, unyt.unyt_array],
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
            self, path: str, name: Optional[str] = None,
            ) -> Dict[str, Union[np.array, unyt.unyt_array]]:
        """Read datasets from a group.

         Parameters
        ----------
        path : str
            Path to the group where the datasets belong to.
        name : str, optional
            Instead of reading all datasets, read only this one.
        """
        def read_dataset(dataset):
            """Read dataset with associated units."""
            try:
                units = unyt.unyt_quantity.from_string(
                    dataset.attrs["units"].decode("utf-8")
                )
            except KeyError:
                units = 1
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

    def __init__(self, filename: str, file_exists: bool) -> None:
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
        if file_exists:
            if not os.path.isfile(filename):
                raise FileNotFoundError(f"The file '{filename}' does not exist.")
            with h5py.File(filename, "a"):
                pass
        else:
            if os.path.isfile(filename):
                raise FileExistsError(f"The file '{filename}' already exists.")
            with h5py.File(filename, "w"):
                pass

        with HDF5Interface(self.filename) as h5:
            h5.set_node(node="results")

            input_dict = {}
            for category in metadata.input_categories:
                h5.set_node(node=category)
                input_dict[category] = ""

            for run in h5:
                if not run.startswith(metadata.run_variants):
                    continue
                h5.set_datagroup(node="results", name=run, **input_dict)

    def read_node(self, node: str) -> Tuple[str, List[str], List[str]]:
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
        node = "results" if node == "root" else node
        qids, variants = [], []
        with HDF5Interface(self.filename) as h5:
            active = h5.get_node(node)
            for datagroup in h5[node]:
                name = datagroup.split("/")[-1]
                if name not in metadata.input_categories:
                    qids.append(name[-QIDLEN:])
                    variants.append(name[:-QIDLEN-1])

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
        with HDF5Interface(self.filename) as h5:
            if node == "root":
               h5.set_node("results", active)
            else:
                h5.set_node(node, active)

    def read_input(self, qid: str, variant: str) -> Tuple[str, str]:
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
        category = metadata.get_input_category(variant)
        with HDF5Interface(self.filename) as h5:
            return h5.get_datagroup(category, name)

    def write_input(self, meta: MetaData) -> None:
        """Create an input variant with the given metadata.

        Parameters
        ----------
        meta : MetaData
            Input's metadata.
        """
        name = f"{meta.variant}_{meta.qid}"
        category = metadata.get_input_category(meta.variant)
        with HDF5Interface(self.filename) as h5:
            h5.set_datagroup(
                node=category, name=name, date=meta.date, note=meta.note,
                )

    def read_run(self, qid: str, variant: str) -> Tuple[str, str, List[str]]:
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
        with HDF5Interface(self.filename) as h5:
            date, note, usedinputs = h5.get_datagroup("results", name)
        usedinputs = {
            category: qid for category, qid in usedinputs.items() if qid
            }
        return date, note, list(usedinputs.values())

    def write_run(self, meta: MetaData, usedinputs: Dict[str, str]) -> None:
        """Create a run variant with the given metadata.

        Parameters
        ----------
        meta : MetaData
            Run's metadata.
        usedinputs : Dict[str, str]
            Key-value pairs of input category and QID of the input used from
            that category in the simulation.
        """
        name = f"{meta.variant}_{meta.qid}"
        for category in metadata.input_categories:
            if category not in usedinputs:
                usedinputs[category] = ""

        with HDF5Interface(self.filename) as h5:
            h5.set_datagroup(
                node="results", name=name, date=meta.date, note=meta.note,
                **usedinputs,
                )

    def set_note(self, qid: str, variant: str, note: str) -> None:
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
        if variant in metadata.run_variants:
            node = "results"
        else:
            node = metadata.get_input_category(variant)

        with HDF5Interface(self.filename) as h5:
            h5.set_datagroup(node, name, note=note)

    def remove_variant(self, qid: str, variant: str) -> None:
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
        if variant in metadata.run_variants:
            node = "results"
        else:
            node = metadata.get_input_category(variant)

        with HDF5Interface(self.filename) as h5:
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

    def write_datasets(
            self,
            qid: str,
            variant: str,
            data: Dict[str, Union[np.array, unyt.unyt_array]],
            subpath: str = None,
            ) -> None:
        """Write datasets to a data group.

        Parameters
        ----------
        qid : str
            QID of the data.
        variant : str
            Data variant.
        data : dict[str, np.array or unyt.array]
            Name-value pairs of data to be written.
        subpath : str, optional
            Path within the data group if the data is located within a subgroup.
        """
        if variant in metadata.run_variants:
            path = f"results/{variant}_{qid}"
        else:
            category = metadata.get_input_category(variant)
            path = f"{category}/{variant}_{qid}"
        if subpath:
            path += f"/{subpath}"
        with HDF5Interface(self.filename) as h5:
            h5.require_group(path)
            h5.write_datasets(path, data)

    def read_datasets(
            self,
            qid: str,
            variant: str,
            name: Optional[str] = None,
            subpath: Optional[str] = None,
            ) -> Dict[str, Union[np.array, unyt.unyt_array]]:
        """Read the data from the HDF5 file.

        Parameters
        ----------
        qid : str
            QID of the data.
        variant : str
            Data variant.
        name : str, optional
            Instead of reading all datasets, read only this one.
        subpath : str, optional
            Path within the data group if the data is located within a subgroup.

        Returns
        -------
        data : dict[str, np.array or unyt.array] or np.array or unyt.array
            Name-value pairs of read data or just the value of the data if
            `name` was specified.
        """
        if variant in metadata.run_variants:
            path = f"results/{variant}_{qid}"
        else:
            category = metadata.get_input_category(variant)
            path = f"{category}/{variant}_{qid}"
        if subpath:
            path += f"/{subpath}"
        with HDF5Interface(self.filename) as h5:
            return h5.read_datasets(path, name=name)
