import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import subprocess
import h5py
import cantera as ct

import settings_matplotlib as setMPL


def read_IDT_hdf(read_file_hdf):

    f_h5 = h5py.File(read_file_hdf, 'r')
    dset_name = list(f_h5.keys())
    T_array = f_h5['/T'][:]
    p_array = f_h5['/p_atm'][:]
    dsetTau1_name = 'tau_oh_peak'
    dsetTau2_name = 'tau_T_grad'
    tau_array1 = f_h5['/tau/{}'.format(dsetTau1_name)][:]
    tau_array2 = f_h5['/tau/{}'.format(dsetTau2_name)][:]
    f_h5.close()
    
    tau1 = []
    tau2 = []
    count = 0
    for y in p_array:
        tmp_tau1 = []
        tmp_tau2 = []
        for x in T_array:
            tmp_tau1.append(tau_array1[count][3])
            tmp_tau2.append(tau_array2[count][3])
            count += 1
        tau1.append(tmp_tau1)
        tau2.append(tmp_tau2)
 
    df = pd.DataFrame()
    df['T'] = T_array
    df['T_inv'] = 1000.0/T_array
    tau1_p_columns = ['tau1_p{}'.format(p) for p in p_array]
    tau2_p_columns = ['tau2_p{}'.format(p) for p in p_array]
    df_tau1_p = pd.DataFrame(tau1, index=tau1_p_columns).T
    df_tau2_p = pd.DataFrame(tau2, index=tau2_p_columns).T
    df = pd.concat([df, df_tau1_p], axis=1)
    df = pd.concat([df, df_tau2_p], axis=1)

    return df


if __name__ == '__main__':

    phi = 1.0
    p_atm = 1.0 # atm
    Tini_list = [600.0, 700.0, 800.0, 900.0] # K

    chemPath_cti = 'KineticModels/Chem/nHeptane/Lu/sk88/chem.yaml'

    read_dir = 'Results_timeHist_IDT'
    save_dir = 'Results_timeHist_IDT'
    df_IDT = read_IDT_hdf('./{}/Results_IDT_Lu_phi{}.h5'.format(read_dir, phi))

    clist = ['r', 'b', 'g', 'y']
    clist = ['violet', 'cornflowerblue', 'burlywood', 'darkseagreen']
    clist = ['#CB4335', '#2E86C1', '#229954', '#CA6F1E']

    tmax = 0.0

    cols = ['phi', 'p [atm]', 'T [K]', 'tau1 [s]', 'tau2 [s]']
    df_tau_table = pd.DataFrame(columns=cols)

    fig, ax = plt.subplots(nrows=2, ncols=1, sharex='col', figsize=(11,9))

    for i, T in enumerate(reversed(Tini_list)):

        p = p_atm*ct.one_atm 
        read_timeHist_hdf = '{}/timeHist_Lu_phi{}_hdf/timeHistory_p{}_T{}.h5'.format(read_dir, phi, p, T)
        print(read_timeHist_hdf)

        gas = ct.Solution(chemPath_cti)
        states = ct.SolutionArray(gas)
        states.read_hdf(read_timeHist_hdf, group='p{}_T{}'.format(p, T))

        c = clist[i]
        ax[0].plot(states.t, states.T, alpha=1.0, lw=3, ls='-', c=c, label='{} K'.format(T))
        ax[1].plot(states.t, states.heat_release_rate, lw=3, ls='-', c=c, label='')

        tau1_s = float(df_IDT[df_IDT['T']==T]['tau1_p{}'.format(p_atm)])
        tau2_s = float(df_IDT[df_IDT['T']==T]['tau2_p{}'.format(p_atm)])
        tau_ms = tau2_s*1e+3

        if (tmax < tau2_s): tmax = tau2_s
        tlim = [0.0, tmax*1.07]
        HRRlim_log = [5e+3, 2e+11]
        ax[0].text(tau2_s, 1700.0, '{:.0f} K'.format(T), 
                   ha='center', backgroundcolor='w', c=c, fontsize=20)
        if phi == 1.0 and T == 800.0:
            ax[1].text(tau2_s, HRRlim_log[1]*4, '{:.1f} ms'.format(tau_ms), 
                       ha='center', va='bottom', c=c, fontsize=20)
        else:
            ax[1].text(tau2_s, HRRlim_log[1], '{:.1f} ms'.format(tau_ms), 
                       ha='center', va='bottom', c=c, fontsize=20)

        df_tau_table.loc[i] = [phi, p, T, tau1_s, tau2_s]
    df_tau_table = df_tau_table.sort_values(by=['T [K]'])
    file_path = '{}/tau_phi{}_{}atm.csv'.format(save_dir, phi, p_atm)
    df_tau_table.to_csv(file_path, index=False)

    ax[0].set_title(r'$n$-heptane/air, $\phi$ = {}, {} atm'.format(phi, p_atm), fontsize=24)
    ax[0].set_ylabel(r'Temperature [K]')
    ax[1].set_ylabel('Heat release rate [J/m$^3$-s]')
    ax[1].set_xlabel(r't [s]')

    ax[1].set_xlim(tlim)
    ax[1].set_ylim(HRRlim_log)
    ax[1].set_yscale('log')

    file_name = 't_THRR_phi{}_{}atm.png'.format(phi, p_atm)
    file_path = '{0}/{1}'.format(save_dir, file_name)
    fig.tight_layout()
    fig.savefig(file_path, bbox_inches = 'tight')
    subprocess.call(['imgcat', file_path])


# # read_file_hdf = './Results_IDT.h5'
# read_file_hdf = read_path_base+'/Results_IDT.h5'
# 
# dsetTau_name = 'tau_T_grad'
# tau_array = read_tau_Cantera0D(read_file_hdf, dsetTau_name)
# df_Tgrad = pd.DataFrame(tau_array, columns=['P_Pa', 'P_atm', 'T', 'tau1'])
# print(df_Tgrad)
# 
# dsetTau_name = 'tau_oh_peak'
# tau_array = read_tau_Cantera0D(read_file_hdf, dsetTau_name)
# df_OHpeak = pd.DataFrame(tau_array, columns=['P_Pa', 'P_atm', 'T', 'tau1', 'tau2'])
# print(df_OHpeak)
# 
# # dsetTau_name = dset_name_tau_list[9]
# # tau_array_HRR = read_tau_Cantera0D(read_file_hdf, dsetTau_name)
# # df_HRRpeak = pd.DataFrame(tau_array_HRR, columns=['P_Pa', 'P_atm', 'T', 'tau1', 'tau2', 'tau3'])
# # print(df_HRRpeak)
# 
# endTime_Max  = 5.0
# for p in p_list:
#     data_nout = []
#     for T in T_list:
#         states = ct.SolutionArray(gas)
#         filename = '{}/timeHistory_p{}_T{}.h5'.format(read_dir_hdf, p, T)
#         states.read_hdf(filename, group='p{}_T{}'.format(p, T))
# 
#         df_Tgrad_i = df_Tgrad[(df_Tgrad['P_Pa']==p) & (df_Tgrad['T']==T)]
#         tau_Tgrad = float(df_Tgrad_i['tau1'])
# 
#         df_OHpeak_i = df_OHpeak[(df_OHpeak['P_Pa']==p) & (df_OHpeak['T']==T)]
#         tau_OHpeak1 = float(df_OHpeak_i['tau1'])
# 
#         print('p{}_T{}'.format(p, T))
# 
#         count_oh_peak = 0
#         for i in range(len(states.t)-2):
#             if(   states('OH').Y[i]   < states('OH').Y[i+1] \
#               and states('OH').Y[i+1] > states('OH').Y[i+2] \
#               and count_oh_peak == 0):
#                 count_oh_peak = 1
#                 idx_oh_peak1 = i
#         if(count_oh_peak==0):
#             tau_oh_peak1 = endTime_Max
#             tau_oh_peak2 = endTime_Max
# 
#         def make_patch_spines_invisible(ax):
#             ax.set_frame_on(True)
#             ax.patch.set_visible(False)
#             for sp in ax.spines.values():
#                 sp.set_visible(False)
# 
#         Tmin, Tmax = 500, 2000
#         HHRmin, HHRmax = 1e+4, 2e+13
#         Tmin, Tmax = 500, 4000
#         HHRmin, HHRmax = 1e-4, 2e+15
#         cm = plt.cm.get_cmap('tab20')
#         # sp_list = ['nC7H16','nC7H15','nC7H15OO','nC7H14OOH','nHOOC7H14OO',\
#         #            'CH2O', 'H2O2', 'OH', 'HO2', 'CO', 'CO2']
#         # sp_list = ['nC7H16', 'CH2O', 'H2O2', 'OH', 'HO2', 'CO', 'CO2']
#         # sp_list = ['NXC12H26', 'CH2O', 'H2O2', 'OH', 'HO2', 'CO', 'CO2']
#         sp_list = ['nC7H16', 'iC8H18', 'CH2O', 'H2O2', 'OH', 'HO2', 'CO', 'CO2']
#         fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
# 
#         ax2 = ax1.twinx()
#         ax3 = ax1.twinx()
#         ax3.spines["right"].set_position(("axes", 1.15))
#         make_patch_spines_invisible(ax3)
#         ax3.spines["right"].set_visible(True)
# 
#         for i, sp in enumerate(sp_list):
#             ax1.plot(states.t, states(sp).X, lw=2, c=cm.colors[i], label=sp)
#         ax2.plot(states.t, states.T, lw=2, c='k', label='T')
#         ax3.plot(states.t, states.heat_release_rate, lw=2, ls='--', c='firebrick', label='HRR')
#         # ax1.set_xscale('log')
#         ax1.set_yscale('log')
#         ax3.set_yscale('log')
# 
#         ax1.set_xlim(0, tau_Tgrad+tau_Tgrad/10.0)
#         ax1.set_ylim(1e-14, 5e-1)
#         ax2.set_ylim(Tmin, Tmax)
#         ax3.set_ylim(HHRmin, HHRmax)
# 
#         ax1.set_xlabel(r'Time [s]')
#         ax1.set_ylabel(r'Mole fraction [-]')
#         ax2.set_ylabel(r'T [K]')
#         ax3.set_ylabel(r'Heat release rate [kJ/cm3-sec]')
# 
#         l1 = ax1.legend(loc='lower right', fontsize=10)
#         l1.set_zorder(10)
# 
#         h2, l2 = ax2.get_legend_handles_labels()
#         h3, l3 = ax3.get_legend_handles_labels()
#         l23 = ax2.legend(h2+h3, l2+l3, loc='upper right', fontsize=10)
#         l23.set_zorder(10)
# 
#         ax2.annotate(r'$\mathdefault{T_{in}}$='+'{:.0f}K'.format(T), 
#              xy=(0, T), xytext=(0, T),
#              ha='right', va='bottom', fontsize=11, color='k', 
#              bbox=dict(boxstyle='square', ec= 'gray', fc='w', alpha=1.0))
#         ax2.annotate('', xy=(tau_Tgrad, Tmax), xytext=(0, Tmax),
#             arrowprops=dict(arrowstyle="<|-|>", color='r', linewidth=1.5),
#             fontsize=8)
#         ax2.annotate(r'$\mathdefault{IDT_{hot}}$='+'{:1.2e}s'.format(tau_Tgrad), xy=(0, 1750), 
#              xytext=(tau_Tgrad/2.0, 1750),
#              ha='center', va='center', fontsize=11, color='k', 
#              bbox=dict(boxstyle='round', ec= 'r', fc='w', alpha=0.6))
#         ax2.annotate('', xy=(tau_OHpeak1, T-30), xytext=(0, T-30),
#             arrowprops=dict(arrowstyle="<|-|>", color='b', linewidth=1.5),
#             va='center', fontsize=8)
#         ax2.annotate(r'$\mathdefault{IDT_{OHpeak_{1st}}}$='+'{:1.2e}s'.format(tau_OHpeak1), 
#              xy=(0, T-50), xytext=(tau_OHpeak1, T-50),
#              ha='center', va='top', fontsize=11, color='k', 
#              bbox=dict(boxstyle='round', ec= 'b', fc='w', alpha=0.6))
#         try:
#             if states.T[idx_oh_peak1] < Tmax:
#                 ax2.annotate(r'$\mathdefault{T_{OHpeak_{1st}}}$='+'{:.0f}K'.format(states.T[idx_oh_peak1]), 
#                      xy=(0, states.T[idx_oh_peak1]+50), xytext=(tau_OHpeak1, states.T[idx_oh_peak1]+50),
#                      ha='left', va='bottom', fontsize=11, color='k', 
#                      bbox=dict(boxstyle='square', ec= 'gray', fc='w', alpha=0.6))
#         except:
#             pass
# 
#         save_filename = 'p{}_T{}'.format(p, T)
# 
#         # plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
#         # fig.set_figheight(6)
#         # fig.set_figwidth(8)
# 
#         fig.tight_layout()
# 
#         plt.savefig('{0}/{1}.png'.format(save_dir, save_filename))
#         plt.clf()
#         plt.close()
#         # plt.savefig('{0}/{1}.pdf'.format(save_dir, save_filename), bbox_inches = 'tight')


