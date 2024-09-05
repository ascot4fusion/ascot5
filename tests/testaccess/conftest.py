"""Set up the coreio tests."""
from unittest.mock import MagicMock

from a5py.data.access import metadata
from a5py.data.access.treeparts import Leaf

FNTEST = "test.h5"
FNEMPTY = "empty.h5"
QID1 = "2991786571"
QID2 = "9753987342"
QID3 = "4404229430"
QID4 = "0963810214"
QID5 = "5960585966"
DATE_FRI = "1997-08-29 02:14:00"
DATE_SAT = "1997-08-30 02:14:00"
DATE_SUN = "1997-08-31 02:14:00"
INPUTVAR = "input"
INPUTVAR2 = "E_TC"
OUTPUTVAR = "output"
CATEGORY = "wall"
DIAGNOSTIC = "diagnostic"
DATE = DATE_FRI
NOTE = "Let off some steam Bennett"

# Make Ascot recognize our dummy variants
metadata.data_variants[CATEGORY] = (
    metadata.data_variants[CATEGORY] + (INPUTVAR,)
)
metadata.run_variants += (OUTPUTVAR,)
# Make Ascot recognize our dummy variants
metadata.data_variants[CATEGORY] = (
    metadata.data_variants[CATEGORY] + (INPUTVAR,)
)
metadata.simulation_diagnostics += (DIAGNOSTIC,)


def create_leaf(
        qid: str, date: str = DATE, variant: str = INPUTVAR, note: str = NOTE,
        ) -> Leaf:
    """Create a `Leaf` instance with mock tree manager."""
    leaf = Leaf(qid=qid, date=date, variant=variant, note=note)
    leaf._treemanager = MagicMock()
    return leaf
