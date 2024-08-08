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
description field that user can modify at will.

The groups should be created, removed, and their metadata accessed via this
module. The actual datasets are accessed via corresponding modules. This module
assumes the HDF5 file is open when these routines are called, and closed
afterwards.
"""
from __future__ import annotations

from importlib.metadata import version as importlib_version
from collections import namedtuple
from typing import List, Dict, Tuple

import h5py
import numpy as np
import unyt

from . treestructure import (
    input_categories, run_variants, get_input_category, get_qid, QIDLEN,
    )
from a5py.exceptions import AscotNoDataException, AscotIOException

VERSION = importlib_version("a5py")
"""Current version of the code."""

def add_group():
    pass

class HDF5Interface(h5py.File):

    def __init__(self, filename, create=False):
        """Initialize HDF5Interface.

        Groups corresponding to input categories are created if they do not
        already exist and active attributes are initialized empty.

        Parameters
        ----------
        filename : str
            HDF5 file name.
        create : bool, optional
            Create the HDF5 file if it does not exist.
        """
        super().__init__(filename, "a")
        with h5py.File(filename, "a") as h5:
            if "active" not in self.attrs:
                self.attrs["active"] = np.bytes_("")
            for category in input_categories:
                if category not in self:
                    group = self.create_group(category)
                    group.attrs["active"] = np.bytes_("")

            for name, group in self.items():
                if name.startswith(run_variants):
                    for category in input_categories:
                        group.attrs[category] = np.bytes_("")

    def set_input(self, metadata):
        """Set metadata of an input data creating the group if one doesn't
        exist.
        """
        category = get_input_category(metadata.variant)
        group = self[category].create_group(f"{metadata.variant}_{metadata.qid}")
        group.attrs["date"] = np.bytes_(metadata.date)
        group.attrs["description"] = np.bytes_(metadata.description)

    def set_node(self, category: str, active: str) -> None:
        """Set metadata of the root or an input category creating the group if
        one doesn't exist.

        This routine does not touch the simulation output or input data
        belonging to this category.

        Parameters
        ----------
        category : str
            Name of the input category or empty string for root.
        """
        if category:
            try:
                group = self.create_group(category)
            except Exception:
                group = self[category]
        else:
            group = self
        group.attrs["active"] = np.bytes_(active)

    def set_simulation_output(self, simulation: str):
        """Set metadata of a simulation output data creating the group if one
        doesn't exist.

        Diagnostics data is not stored by this routine.
        """
        group = self["results"].create_group(
            f"{simulation.variant}_{simulation.qid}",
            )
        group.attrs["date"] = np.bytes_(simulation.date)
        group.attrs["description"] = np.bytes_(simulation.description)

        for category in input_categories:
            try:
                qid = simulation[category].qid
                group.attrs[category] = np.bytes_(qid)
            except AscotIOException:
                continue

    def read_input(self, ):
        """
        """
        return None

    def read_input_category(
            self,
            category: str,
            ) -> Tuple[str, List[Dict[str, str]]]:
        """Read metadata from an input category group.

        Parameters
        ----------
        category : str
            Name of the input category.

        Returns
        -------
        active : str
            QID of the active group in the input category.
        metadata : list of dict [str, str]
            Metadata ("qid", "date", "description", "variant") of the input
            groups in the category.
        """
        metadata = []
        for name, group in self[category].items():
            metadata.append({
                "qid":get_qid(name),
                "variant":name[:-(QIDLEN + 1)],
                "description":group.attrs["desc"].decode("utf-8"),
                "date":group.attrs["date"].decode("utf-8"),
                })
        return self[category].attrs["active"].decode("utf-8"), metadata

    def read_simulation_output(self):
        """Read metadata from a simulation output.
        """
        metadata = []
        for name, group in self.items():
            if not name.startswith(run_variants):
                continue
            inputs = []
            for category in input_categories:
                if category in group.attrs:
                    inputs.append(group.attrs[category].decode("utf-8"))

            metadata.append({
                "qid":get_qid(name),
                "variant":name[:-(QIDLEN + 1)],
                "description":group.attrs["desc"].decode("utf-8"),
                "date":group.attrs["date"].decode("utf-8"),
                "inputs":inputs
                })
        return self.attrs["active"].decode("utf-8"), metadata

    def read_tree(self):
        """
        """
        return inputs, runs, active


def write_data(group, name, data, shape, dtype, unit=None, compress=False):
    """Write a dataset.

    The shape of the written dataset is same as the input array.

    Parameters
    ----------
    group : `h5py.Group`
        HDF5 group where the dataset will be written.
    name : str
        Name of the new dataset.
    data : `np.array`
        Data to be written.
    dtype : str
        Data type.
    unit : str, optional
        Unit string if the data has units.
    compress : boolean, optional
        Use gzip compression
    """
    if compress:
        g = group.create_dataset(
            name  = name,
            shape = shape,
            data  = data,
            dtype = dtype,
            compression="gzip", compression_opts=9
        )
    else:
        g = group.create_dataset(
            name  = name,
            shape = shape,
            data  = data,
            dtype = dtype
        )
    if unit is not None:
        g.attrs.create("unit", np.bytes_(unit))

def read_data(group, name):
    """Read a dataset including its units if present.

    Parameters
    ----------
    group : `h5py.Group`
        HDF5 group containing the dataset.
    name : str
        Name of the dataset.

    Returns
    -------
    data : `np.array` or `unyt.array`, `(N,)`
        Dataset with units read from the dataset attribute "unit".

        If dataset has no "unit" attribute, ordinary `np.array` is returned.
    """
    unit = 1
    if "unit" in group[name].attrs.keys():
        unit_str = group[name].attrs["unit"]
        unit     = unyt.Unit(unit_str)

    return group[name][:].ravel() * unit
