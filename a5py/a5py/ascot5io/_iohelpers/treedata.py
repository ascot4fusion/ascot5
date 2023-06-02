"""
Contains definition of DataContainer class.
"""

import h5py

from . import ascot5file

class DataGroup():
    """
    Class inherited by the input and result groups in Ascot treeview.

    This class provides methods to change and view the meta data (QID,
    name, type, description, date) and access to the raw data. The
    raw data is accessed (via h5py) as
    with datacontainer as d:
      d["data"][:]
    where d is the h5py.group corresponding to this container.

    Attributes:
      _root : RootNode
        The rootnode this object belongs to.
      _path : str
        Path to this data within the HDF5 file.
      _opened : bool
        If True, the raw data can be accessed via h5py.
    """

    def __init__(self, root, hdf5group):
        """
        Initialize container.

        Args:
          root : RootNode
            The rootnode this object belongs to.
          hdf5group : h5py.group
            The h5py group this data corresponds to.
        """
        self._root   = root
        self._path   = hdf5group.name
        self._opened = None

    def __enter__(self):
        """
        Call _open when entering with clause.

        Returns:
          h5py.group
            HDF5 group corresponding to this data.
        """
        return self._open()

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Call _close when exiting with clause.
        """
        self._close()

    def _open(self):
        """
        Open and return the HDF5 group corresponding to this data.

        Returns:
          h5py.group
            HDF5 group corresponding to this data.
        Raises:
          AscotIOException
            When this instance has already opened the file.
        """
        if self._opened is not None:
            raise AscotIOException(
                """The file has already been opened by this instance and
                   not closed by _close()""")

        self._opened = h5py.File(self._file, "a")[self._path]
        return self._opened

    def _close(self):
        """
        Close the HDF5 group corresponding to this data.
        """
        self._opened.file.close()
        self._opened = None

    def set_desc(self, desc):
        """
        Set description for this group.

        Args:
          desc : str
            Short description for the user to document this group.
        """
        f = self._open()
        val = ascot5file.set_desc(f.file, self._group,desc)
        self._close()
        self._root._build(self._root.file_getpath())
        return val

    def get_desc(self):
        """
        Get this group's description.

        Returns:
          str
            Documentation the user has used to describe this group.
        """
        f = self._open()
        val = ascot5file.get_desc(f.file, self._group)
        self._close()
        return val

    def get_date(self):
        """
        Get the date when this group was created.

        Returns:
          The date in XXX format.
        """
        f = self._open()
        val = ascot5file.get_date(f.file, self._group)
        self._close()
        return val

    def get_qid(self):
        """
        Get QID of this group.

        Returns:
          str
            String with 10 characters containing numbers from 0-9 which is
            an unique identifier for this group.
        """
        return ascot5file.get_qid(self._group)

    def get_type(self):
        """
        Return type of this group.

        Returns:
          str
            Group type.
        """
        path = self._path.split("/")[-1]
        return path[:-11]

    def get_name(self):
        """
        Return name of this group.

        Returns:
          str
            The name of the group as "<group type>_<group qid>".
        """
        return self._path.split("/")[2]

    def activate(self):
        """
        Set this group as active.
        """
        self._root._activate(self.get_qid())

    def destroy(self, repack=True):
        """
        Remove this group from the HDF5 file.
        """
        self._root._destroy(self.get_qid(), repack)

    def copy_to_hdf5file(self,target_file,newgroup=False):

        import a5py.ascot5io.ascot5tools as tools
        group = tools.copygroup(self._root.file_getpath(), target_file,
                                self._path, newgroup=newgroup)
        return group
