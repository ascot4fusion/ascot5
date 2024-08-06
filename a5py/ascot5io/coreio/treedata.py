from __future__ import annotations
"""Abstract classes for objects containing data.
"""
import h5py
import random
import datetime
from enum import Enum
from ... import utils

import numpy as np

from a5py.exceptions import AscotIOException
from . import fileapi
from . import tools

class Address():
    """Data class that contains information on how to access the data.
    """

    class Format(Enum):

        IMAS_IDS = 1
        HDF5_FILE = 2
        IN_MEMORY = 3

    def __init__(self,
                 hdf5_filename=None,
                 path_within_hdf5=None,
                 imas_ids=None,
                 pointer_to_c_struct=None,
                 **kwargs
                 ) -> None:
        super().__init__(**kwargs)
        self.imas_ids = imas_ids
        self.hdf5_filename = hdf5_filename
        self.path_within_hdf5 = path_within_hdf5
        self.pointer_to_c_struct = pointer_to_c_struct

        self.format = None
        inconsistent_address = False
        if imas_ids is not None:
            self.format = Address.Format.IMAS_IDS
        if hdf5_filename is not None and path_within_hdf5 is not None:
            self.format = Address.Format.HDF5_FILE
        if pointer_to_c_struct is not None:
            self.format = Address.Format.IN_MEMORY
        no_address = self.format is None

        if inconsistent_address:
            raise ValueError("Inconsistent address specified.")
        if no_address:
            raise ValueError("No address specified.")


    @classmethod
    def from_hdf5(cls, hdf5_filename, path_within_hdf5) -> Address:
        """Create an Address instance from an HDF5 file.

        Parameters
        ----------
        hdf5_filename : str
            Path to HDF5 file.

        Returns
        -------
        address : Address
            An instance of Address.
        """
        return cls(
            hdf5_filename=hdf5_filename,
            path_within_hdf5=path_within_hdf5,
            )

    @classmethod
    def from_imas(cls, imas_ids) -> Address:
        """Create an Address instance from IMAS IDS.

        Parameters
        ----------
        imas_ids : list of str
            List of IMAS IDS.

        Returns
        -------
        address : Address
            An instance of Address.
        """
        return cls(imas_ids=imas_ids)

class DataHolder():
    """Object for accessing data via its attributes but the actual data can be
    stored in memory or on disk.

    Note that if the data is stored on disk, it is written there in the moment
    the attributes are set.
    """

    def __init__(self, **kwargs):
        """
        """
        self._address = common.Address()
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __getattr__(self, key):
        """Retrieve value of the attribute whether it is stored in the memory or
        on disk.
        """
        notdata = key.startswith("_")
        if notdata or self._address.in_memory:
            return super().__getattr__(key)
        else:
            return self._read_from_hdf5(key)

    def __setattr__(self, key, value):
        """Set value of the attribute whether it is stored in the memory or on
        disk.
        """
        notdata = key.startswith("_")
        if notdata or self._address.in_memory:
            super().__setattr__(key, value)
        else:
            self._write_to_hdf5(key, value)

    @classmethod
    def _init_from_hdf5(cls):
        """Initialize the object from a HDF5 file.

        With this method the attributes and the values are read from the HDF5
        file instead of writing them there on the moment the attributes are set.
        Use this to initialize object whose data already exists in the file.
        """
        pass

class DataContainer():
    """Interface for providing access to underlying HDF5 data.

    This class provides methods to change and access the raw data. The
    raw data is accessed (via :obj:`h5py`) as

    .. code-block:: python

       with datacontainer as d:
           d["data"][:]

    where ``d`` is the :obj:`h5py.Group` corresponding to this container.

    Attributes
    ----------
    _root : :class:`RootNode`
        The rootnode this object belongs to.
    _path : str
        Path to this data within the HDF5 file.
    _opened : list [obj]
        List containing a single object, which is this data's HDF5 group if
        the data is being accessed and `None` otherwise.

        We are storing single item list, instead of storing the item directly,
        since this class is inherited by `Note` objects whose attributes cannot
        be changed directly. However, if the attribute is a list, the contents
        of that list can be changed ;).
    """

    def __init__(self, root, path, **kwargs):
        """Initialize data container.

        Parameters
        ----------
        root : :class:`RootNode`
            The root node this data container belongs to.
        path : str
            Path to this data within the HDF5 file.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        self._root   = root
        self._path   = path
        self._opened = [None]
        super().__init__(**kwargs)

    def __enter__(self):
        """Call _open when entering ``with`` clause.

        Returns
        -------
        data : :class:`h5py.Group`
            HDF5 group corresponding to this data.
        """
        return self._open()

    def __exit__(self, exception_type, exception_value, traceback):
        """Call _close when exiting ``with`` clause.
        """
        self._close()

    def _open(self, mode="r"):
        """Open and return the HDF5 group corresponding to this data.

        Returns
        -------
        data : :class:`h5py.Group`
            HDF5 group corresponding to this data.
        mode : {"r", "a"}
            Is the file opened for (r)eading or (a)ppending.

        Raises
        ------
        AscotIOException
            Raised if this instance has already opened the file.
        """
        if self._opened[0] is not None:
            raise AscotIOException(
                "File already opened by this instance")

        fn = self._root._ascot.file_getpath()
        self._opened[0] = h5py.File(fn, mode)[self._path]
        return self._opened[0]

    def _close(self):
        """Close the HDF5 group corresponding to this data.
        """
        self._opened[0].file.close()
        self._opened[0] = None

class DataGroup(DataContainer):
    """Data container that also has meta data (QID, date, description).
    """

    def set_desc(self, desc):
        """Set description for this group.

        Note that the first word in the description is this group's tag which
        can be used to refer to this group.

        Parameters
        ----------
        desc : str
            Short description for the user to document this group.
        """
        f = self._open("a")
        fileapi.set_desc(f.file, self._path, desc)
        self._close()
        self._root._build(self._root._ascot.file_getpath())

    def get_desc(self):
        """Get this group's description.

        Returns
        -------
        desc : str
            Documentation the user has used to describe this group.
        """
        f = self._open()
        val = fileapi.get_desc(f.file, self._path)
        self._close()
        return val

    def get_date(self):
        """Get the date when this group was created.

        Returns
        -------
        date : str
            The date in YYYY-MM-DD hh:mm:ss format.
        """
        f = self._open()
        val = fileapi.get_date(f.file, self._path)
        self._close()
        return val

    def get_qid(self):
        """Get QID of this group.

        Returns
        -------
        qid : str
            String with 10 characters containing numbers from 0-9 which is
            an unique identifier for this group.
        """
        val = fileapi.get_qid(self._path)
        return val

    def get_type(self):
        """Return type of this group.

        Returns
        -------
        gtype : str
            Group type.
        """
        path = self._path.split("/")[-1]
        return path[:-11]

    def get_name(self):
        """Return name of this group.

        Returns
        -------
        name : str
            The name of the group as "<group type>_<group qid>".
        """
        return self._path.split("/")[2]

    def activate(self):
        """Set this group as active.

        Active inputs are used when the simulation is run. Active groups are
        also used during post-processing.
        """
        self._root._activate_group(self.get_qid())

    def destroy(self, repack=True):
        """Remove this group from the HDF5 file.

        Parameters
        ----------
        repack : bool, optional
            If True, repack the HDF5 file.

            Removing data from the HDF5 file only removes references to it and
            repacking is required to free the disk space. Repacking makes a copy
            of the HDF5 file and destroys the original, and it also requires
            3rd party tool `h5repack` which is why it's use is optional here.
        """
        self._root._destroy_group(self.get_qid(), repack)

    def export(self, target_file, newgroup=False):
        """Copy this group with its contents to another HDF5 file.

        Parameters
        ----------
        target_file : str
            Path to the file where the data is copied to.
        newgroup : bool, optional
            If True, a new QID and date is generated for the copied group.
        """
        group = tools.copygroup(self._root.file_getpath(), target_file,
                                self._path, newgroup=newgroup)
        return group
