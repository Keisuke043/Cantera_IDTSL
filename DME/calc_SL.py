
import cantera as ct
import numpy as np
from mpi4py import MPI
import h5py
import subprocess


fuel = 'CH3OCH3'
oxid = 'O2:1, N2:3.76'

chemPath_DTU_cti = 'KineticModels/Chem/DME/DTUmechanism/chem.cti'
chemPath_Lu_cti  = 'KineticModels/Chem/DME/Lu/sk39/chem.cti'
chemPath_dic = {'DTU': chemPath_DTU_cti, 'Lu':chemPath_Lu_cti}

key_chem = 'DTU'
key_chem = 'Lu'
chemPath_cti = chemPath_dic[key_chem]

min_temp    = 300.0
max_temp    = 1000.0
del_temp    = 50.0
min_presAtm = 1.0
max_presAtm = 10.0
max_presAtm = 1.0
del_presAtm = 9.0

equivalence_ratio = 1.0

endTime_Max = 5.0

savedir = 'Results_SL'
subprocess.call(['mkdir','-p', savedir])

flag_save_profile = True
if flag_save_profile:
    savedir_profile = 'Results_1Dprofile/profile_{}_hdf'.format(key_chem)
    subprocess.call(['mkdir','-p', savedir_profile])

array_temp    = np.arange(min_temp, max_temp+del_temp/10.0, del_temp)
array_presAtm = np.arange(min_presAtm, max_presAtm+del_presAtm/10.0, del_presAtm)
array_pres    = ct.one_atm*array_presAtm

array_width = np.array([0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002])  # m
array_width = np.array([0.05])  # m

# Set case number in each process
caseNumber_width = array_width.shape[0]
caseNumber_temp  = array_temp.shape[0]
caseNumber_pres  = array_pres.shape[0]
caseNumber_all   = caseNumber_width*caseNumber_temp*caseNumber_pres

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

caseNumber = 0
casePTdict = {}
for cur_width in array_width:
    for cur_pres in array_pres:
        for cur_temp in array_temp:
            casePTdict[caseNumber] = [cur_width, cur_pres, cur_temp]
            caseNumber += 1

filename_sl_hdf = '{}/Results_SL_{}_phi{}.h5'.format(savedir, key_chem, equivalence_ratio)
f_sl = h5py.File(filename_sl_hdf, 'w', driver='mpio', comm=MPI.COMM_WORLD)

dsetwidth = f_sl.create_dataset('width',(caseNumber_width,),dtype="f")
dsetT     = f_sl.create_dataset('T',    (caseNumber_temp,), dtype="f")
dsetp     = f_sl.create_dataset('p',    (caseNumber_pres,), dtype="f")
dsetpAtm  = f_sl.create_dataset('p_atm',(caseNumber_pres,), dtype="f")
dsetwidth[:] = array_width
dsetT[:]     = array_temp
dsetp[:]     = array_pres
dsetpAtm[:]  = array_presAtm

grp = f_sl.create_group('SL')
dsetSL = grp.create_dataset('SL', (caseNumber_all,5), dtype="f")

gas = ct.Solution(chemPath_cti)

for caseNumber in range(caseNumber_sta, caseNumber_end+1):

    cur_width = float(casePTdict[caseNumber][0])
    cur_p     = float(casePTdict[caseNumber][1])
    cur_p_atm = cur_p/ct.one_atm
    cur_T     = float(casePTdict[caseNumber][2])

    gas.TP = (cur_T, cur_p)
    gas.set_equivalence_ratio(equivalence_ratio, fuel=fuel, oxidizer=oxid)

    # Set up flame object
    f = ct.FreeFlame(gas, width=cur_width)
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)

    # Solve with mixture-averaged transport model
    f.transport_model = 'Mix'
    loglevel = 0  # amount of diagnostic output (0 to 8)

    print('Case = {:>4}, phi={}, width={}, p={}, T={}'.format(caseNumber, \
                              equivalence_ratio, cur_width, cur_p, cur_T))

    try:

        f.solve(loglevel=loglevel, auto=False)

        # Save to HDF container file if h5py is installed

        dsetSL[caseNumber]  = [cur_width, cur_p, cur_p_atm, cur_T, f.velocity[0]]

        if flag_save_profile:

            try:
                f.write_hdf('{}/adiabaticFlame_width{}_phi{}_p{}_T{}.h5'\
                 .format(savedir_profile, cur_width, equivalence_ratio, cur_p, cur_T),
                 group='mix', mode='w', description='solution with mixture-averaged transport')
            except ImportError:
                f.write_csv('{}/adiabaticFlame_width{}_phi{}_p{}_T{}.csv'\
                 .format(savedir_profile, cur_width, equivalence_ratio, cur_p, cur_T), quiet=False)

        print('mixture-averaged flamespeed = {0:7f} m/s, maximum HRR = {1} W/m^3'\
                                  .format(f.velocity[0], max(f.heat_release_rate)))

    except:

        dsetSL[caseNumber]  = [cur_width, cur_p, cur_p_atm, cur_T, 0]
        f.write_hdf('{}/adiabaticFlame_width{}_phi{}_p{}_T{}.h5_failure'\
                    .format(savedir_profile, cur_width, equivalence_ratio, cur_p, cur_T), \
                    group='mix', mode='w', description='solution with mixture-averaged transport')

        print('Failure')

    del f

f_sl.close()

