"""Module for creating and modifying the HDF5 file and for accessing meta data.

ASCOT5 HDF5 file is built as follows. At top level are parent groups, one (and
only one) for each type of input (e.g. bfield for magnetic field input) and one
parent group that contains all the results.

Parent groups contain groups that hold the actual datasets. Single parent can
have multiple groups e.g. bfield can store multiple magnetic field inputs.
The result group is a different in a sense that it contains run groups, one
for each simulation, which then contain the groups that store inistate,
endstate, and whatever diagnostics were used. However, here we refer to the
top level groups as parents, and the groups that are directly below them
as data containers.

A typical structure of a file can be like this:

> /
> /bfield
> /bfield/B_2DS_1234567890
> /bfield/B_2DS_2345678901
> /bfield/B_3DS_3456789012
> ...
> <other input groups>
> ...
> /results
> /results/run_4567890123
> /results/run_5678901234

The number in group name is unique identifier (QID) which is generated when a
group is created. It is generated from a 32 bit random integer and as such it
is almost quaranteed to be unique.

QID is used to determine which input group is active, i.e., to be used in
a simulation, and to which run group the simulation results are written.
QID of the active group is stored as an attribute (in a string format) in the
parent. For the results parent this is the QID of the most recent run (unless
manually changed). Run groups also store QIDs of all input fields used in that
particular run.

All groups also store date and time at which they where created and also a
note field that user can modify at will.
"""
from __future__ import annotations

import os
import tempfile
from typing import List, Dict, Tuple, Optional, Union

import h5py
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

    def set_dataset(
            self,
            node: str,
            name: str,
            date: Optional[str] = None,
            note: Optional[str] = None,
            **additional_attrs,
            ) -> None:
        """Set metadata on a dataset, creating the group if none exists.

        If any attributes specified in the parameters are None, the default
        value is used when creating the group or the value is not changed if
        the group exists.

        Parameters
        ----------
        node : str
            Node to which this dataset belongs to.
        name : str
            Name of the dataset as it appears in the file.
        date : str, optional
            Date when the dataset was created.d
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

    def get_dataset(
            self, node: str, name: str,
            ) -> Union[Tuple[str, str], Tuple[str, str, Dict[str, str]]]:
        """Get metadata from a dataset.

        Parameters
        ----------
        node : str
            Node to which this dataset belongs to.
        name : str
            Name of the dataset as it appears in the file.

        Returns
        -------
        date : str
            Date when the dataset was created.
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
            h5.set_node(node="/")

            input_dict = {}
            for category in metadata.input_categories:
                h5.set_node(node=category)
                input_dict[category] = ""

            for run in h5:
                if not run.startswith(metadata.run_variants):
                    continue
                h5.set_dataset(node="/", name=run, **input_dict)

    def read_node(self, node: str) -> Tuple[str, List[str], List[str]]:
        """Return list of QIDs, variants, and the active QID for datasets within
        this node.

        Parameters
        ----------
        node : str
            Name of the node.

        Returns
        -------
        active : str
            QID of the active dataset.
        qids : List[str]
            List of QIDs belonging to this node.
        variants : List[str]
            List of variants belonging to this node in same order as `qids`.
        """
        node = "/" if node == "root" else node
        qids, variants = [], []
        with HDF5Interface(self.filename) as h5:
            active = h5.get_node(node)
            for dataset in h5[node]:
                name = dataset.split("/")[-1]
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
            QID of the active dataset.
        """
        with HDF5Interface(self.filename) as h5:
            if node == "root":
               h5.set_node("/", active)
            else:
                h5.set_node(node, active)

    def read_input(self, qid: str, variant: str) -> Tuple[str, str]:
        """Return date and note of this input dataset.

        Parameters
        ----------
        qid : str
            Quasi-unique identifier for the dataset.
        variant : str
            Dataset variant.

        Returns
        -------
        date : str
            Date when the dataset was created.
        note : str
            User-defined note of the data.
        """
        name = f"{variant}_{qid}"
        category = metadata.get_input_category(variant)
        with HDF5Interface(self.filename) as h5:
            return h5.get_dataset(category, name)

    def read_output(self, qid: str, variant: str) -> Tuple[str, str, List[str]]:
        """Return date, note, and list of input QIDs for this output dataset.

        Parameters
        ----------
        qid : str
            Quasi-unique identifier for the dataset.
        variant : str
            Dataset variant.

        Returns
        -------
        date : str
            Date when the dataset was created.
        note : str
            User-defined note of the data.
        inputqids : List[str]
            List of QIDs for the inputs used to produce the output.
        """
        name = f"{variant}_{qid}"
        with HDF5Interface(self.filename) as h5:
            date, note, usedinputs = h5.get_dataset("/", name)
        usedinputs = {
            category: qid for category, qid in usedinputs.items() if qid
            }
        return date, note, list(usedinputs.values())

    def write_input(self, meta: MetaData) -> None:
        """Create an input dataset with the given metadata.

        Parameters
        ----------
        meta : MetaData
            Dataset's metadata.
        """
        name = f"{meta.variant}_{meta.qid}"
        category = metadata.get_input_category(meta.variant)
        with HDF5Interface(self.filename) as h5:
            h5.set_dataset(
                node=category, name=name, date=meta.date, note=meta.note,
                )

    def write_output(self, meta: MetaData, usedinputs: Dict[str, str]) -> None:
        """Create an output dataset with the given metadata.

        Parameters
        ----------
        meta : MetaData
            Dataset's metadata.
        usedinputs : Dict[str, str]
            Key value pairs of input category and QID of the input used from
            that category in the simulation.
        """
        name = f"{meta.variant}_{meta.qid}"
        for category in metadata.input_categories:
            if category not in usedinputs:
                usedinputs[category] = ""

        with HDF5Interface(self.filename) as h5:
            h5.set_dataset(
                node="/", name=name, date=meta.date, note=meta.note,
                **usedinputs,
                )

    def set_note(self, qid: str, variant: str, note: str):
        """Set the note attribute in a dataset."""
        name = f"{variant}_{qid}"
        if variant in metadata.run_variants:
            node = "/"
        else:
            node = metadata.get_input_category(variant)

        with HDF5Interface(self.filename) as h5:
            h5.set_dataset(node, name, note=note)

    def remove_dataset(self, qid: str, variant: str) -> None:
        """Removes dataset from the file.

        Does not repack the file.
        """
        name = f"{variant}_{qid}"
        if variant in metadata.run_variants:
            node = "/"
        else:
            node = metadata.get_input_category(variant)

        with HDF5Interface(self.filename) as h5:
            del h5[node][name]

    def contains_dataset(self, dataset):
        """Checks if a dataset exists on the file.

        Parameters
        ----------
        

        Returns
        -------
        """
        with HDF5Interface.contains_dataset(self.filename):

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
