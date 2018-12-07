import h5py

from . ascot5file import get_desc, get_qid, get_date

class AscotInput():

    def __init__(self, hdf5):
        self._file = hdf5.file.filename
        self._path = hdf5.name
        self._isopen = False

    def desc(self):
        f = self._open()
        val = get_desc(f.file, f)
        self._close(f)
        return val

    def date(self):
        f = self._open()
        val = get_date(f.file, f)
        self._close(f)
        return val

    def qid(self):
        return get_qid(self._path)

    def type(self):
        return get_qid(self._path)

    def _open(self):
        if self._isopen == True:
            print("da")
            return None

        self._isopen = True
        return h5py.File(self._file, "r")[self._path]

    def _close(self, f):
        self._isopen = False
        f.file.close()
