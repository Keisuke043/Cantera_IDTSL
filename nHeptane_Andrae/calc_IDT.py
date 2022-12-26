
import numpy as np
from mpi4py import MPI
import h5py
import subprocess
import matplotlib.pyplot as plt
import cantera as ct



fuel = 'C7H16'
oxid = 'O2:1, N2:3.76'
key_chem = 'Andrae'
chemPath_cti = 'Andrae2013_C7H16/chem.yaml'

min_temp    = 500.0
max_temp    = 2000.0
del_temp    = 50.0
min_presAtm = 1.0
max_presAtm = 1.0
del_presAtm = 1.0

equivalence_ratio = 1.0

endTime_Max = 5.0

savedir = 'Results_IDT'
subprocess.call(['mkdir','-p', savedir])

flag_save_history = False
if flag_save_history:
    savedir_timeHistory = 'Results_0DtimeHistory/timeHistory_{}_hdf'.format(key_chem)
    subprocess.call(['mkdir','-p', savedir_timeHistory])

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

caseNumber = 0
casePTdict = {}
for cur_pres in array_pres:
    for cur_temp in array_temp:
        casePTdict[caseNumber] = [cur_pres, cur_temp]
        caseNumber += 1

filename_idt_hdf = '{}/Results_IDT_{}_phi{}.h5'.format(savedir, key_chem, equivalence_ratio)
f_idt = h5py.File(filename_idt_hdf, 'w', driver='mpio', comm=MPI.COMM_WORLD)

dsetT    = f_idt.create_dataset('T',     (caseNumber_temp,), dtype="f")
dsetp    = f_idt.create_dataset('p',     (caseNumber_pres,), dtype="f")
dsetpAtm = f_idt.create_dataset('p_atm', (caseNumber_pres,), dtype="f")
dsetT[:] = array_temp
dsetp[:] = array_pres
dsetpAtm[:] = array_presAtm

grp = f_idt.create_group('tau')

dset_T_grad  = grp.create_dataset('tau_T_grad',  (caseNumber_all, 4), dtype="f")
dset_oh_peak = grp.create_dataset('tau_oh_peak', (caseNumber_all, 5), dtype="f")

gas = ct.Solution(chemPath_cti)
gasEquil = ct.Solution(chemPath_cti)

for caseNumber in range(caseNumber_sta, caseNumber_end+1):

    T = float(casePTdict[caseNumber][1])
    p = float(casePTdict[caseNumber][0])
    p_atm = p/ct.one_atm

    gas.TP = (T, p)
    gas.set_equivalence_ratio(equivalence_ratio, \
                              fuel=fuel, oxidizer=oxid)

    gasEquil.TP = (T, p)
    gasEquil.set_equivalence_ratio(equivalence_ratio, \
                                   fuel=fuel, oxidizer=oxid)
    print(gas.mole_fraction_dict())

    gasEquil.equilibrate('UV')
    equil_temp = gasEquil.T

    r = ct.IdealGasReactor(gas)
    reactorNetwork = ct.ReactorNet([r])
    timeHistory_state = ct.SolutionArray(gas, extra=['t'])

    cur_time = 0.0
    count = 0
    while(cur_time < endTime_Max):
        cur_time = reactorNetwork.step()
        if count % 10 == 0:
            timeHistory_state.append(r.thermo.state, t=cur_time)
        count += 1

    grad_list = []
    for i in range(len(timeHistory_state.t)-1):
        grad_list.append((timeHistory_state.T[i+1] - timeHistory_state.T[i]) \
        / (timeHistory_state.t[i+1] - timeHistory_state.t[i]))
    try:
        idx_T_grad = grad_list.index(max(grad_list))
        if(timeHistory_state.T[idx_T_grad] > 0.5*equil_temp):
            tau_T_grad = timeHistory_state.t[idx_T_grad]
        else:
            tau_T_grad = endTime_Max
    except:
        tau_T_grad = endTime_Max

    count_oh_peak = 0
    for i in range(len(timeHistory_state.t)-2):
        if(   timeHistory_state('OH').Y[i]   < timeHistory_state('OH').Y[i+1] \
          and timeHistory_state('OH').Y[i+1] > timeHistory_state('OH').Y[i+2] \
          and timeHistory_state.t[i] < tau_T_grad \
          and count_oh_peak == 0):
            count_oh_peak = 1
            tau_oh_peak1 = timeHistory_state.t[i]
        elif(timeHistory_state('OH').Y[i]   < timeHistory_state('OH').Y[i+1] \
            and timeHistory_state('OH').Y[i+1] > timeHistory_state('OH').Y[i+2] \
            and count_oh_peak == 1):
            tau_oh_peak2 = timeHistory_state.t[i]
            count_oh_peak = 2
    if(count_oh_peak==0):
        tau_oh_peak1 = tau_T_grad
        tau_oh_peak2 = tau_T_grad
    elif(count_oh_peak==1):
        tau_oh_peak2 = tau_T_grad

    dset_oh_peak[caseNumber] = [p, p_atm, T, tau_oh_peak1, tau_oh_peak2]
    dset_T_grad[caseNumber]  = [p, p_atm, T, tau_T_grad]

    if(flag_save_history):

        timeHistory_state.write_hdf('{}/timeHistory_p{}_T{}.h5'.format(savedir_timeHistory, p, T), \
                     cols=('t', 'T', 'P', 'X', 'heat_release_rate'), group='p{}_T{}'.format(p, T))

    # print('Case = {:>/4}, p,T = {:>7} {:>6}, tau = {:e} '.format(caseNumber, p, T, tau_T_grad))
    print(caseNumber, p, T, tau_T_grad)


f_idt.close()


