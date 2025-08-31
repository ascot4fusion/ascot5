from mpi4py import MPI
from a5py import Ascot
from a5py.templates import PremadeMagneticField


rank, size = MPI.COMM_WORLD.Get_rank(), MPI.COMM_WORLD.Get_size()
print("Rank {rank} started.")
if rank == 0:
    a5 = Ascot()
    template = PremadeMagneticField(a5, field="iter-baseline")
    template.create_input()
else:
    a5 = Ascot()

a5.init(rank)
print(a5.data.bfield.active.rmajor, rank)
print("Rank {rank} completed.")
