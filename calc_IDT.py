
import cantera as ct
import numpy as np
from mpi4py import MPI
import h5py
import subprocess


fuel = 'DME'
oxid = 'O2:1, AR:3.76'
oxid = 'O2:1, N2:3.76'

min_temp    = 800.0
max_temp    = 1000.0
del_temp    = 100.0
min_presAtm = 40.0
max_presAtm = 50.0
del_presAtm = 10.0

equivalence_ratio = 1.0

endTime_Max = 5.0

flag_save_history = True

array_temp    = np.arange(min_temp, max_temp+del_temp/10.0, del_temp)
array_presAtm = np.arange(min_presAtm, max_presAtm+del_presAtm/10.0, del_presAtm)
array_pres    = ct.one_atm*array_presAtm

# Set case number in each process
caseNumber_temp = array_temp.shape[0]
caseNumber_pres = array_pres.shape[0]
caseNumber_all  = caseNumber_temp*caseNumber_pres

# MPI settings
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

work1 = (caseNumber_all)//size
work2 = (caseNumber_all)%size
caseNumber_sta = rank*work1 + min(rank, work2)
caseNumber_end = caseNumber_sta + work1-1
if (work2 > rank): caseNumber_end += 1
print('Rank:{:>2}, Case:{:>3}-{}, Total case: {}'.format(\
       rank, caseNumber_sta, caseNumber_end, caseNumber_all))



