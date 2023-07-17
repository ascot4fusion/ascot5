"""Abstract classes for objects containing data.
"""
import h5py

from a5py.exceptions import AscotIOException
from . import fileapi

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
        val = fileapi.set_desc(f.file, self._path, desc)
        self._close()
        self._root._build(self._root._ascot.file_getpath())
        return val

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

    def export(self,target_file,newgroup=False):
        """
        """

        import a5py.ascot5io.ascot5tools as tools
        group = tools.copygroup(self._root.file_getpath(), target_file,
                                self._path, newgroup=newgroup)
        return group
