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
    

    ###################################################################################

    phi = 1.0
    p_atm = 1.0 # atm
    Tini_list = [700.0] # K
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

    # fig, ax = plt.subplots(nrows=2, ncols=1, sharex='col', figsize=(11,9))

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

        if tau1_s != tau2_s:
            for time, Temp in zip(states.t, states.T):
                if tau1_s < time:
                    temp_tau1 = Temp
                    break
        else:
            temp_tau1 = np.nan
        print(T, temp_tau1, temp_tau1-T)


    # phi_list = np.arange(0.5, 2.0, 0.10)
    # phi_list = map(lambda a:a*0.1, list(range(5, 20, 1)))
    phi_list = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 
                1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    phi_list = [0.8, 2.0]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,7))
    ax2 = ax.twinx()
    # ax2 = ax.twiny()

    for i, phi in enumerate(phi_list):
        print(phi)
        c = clist[i]

        dT_list = []
        Tini_list = [600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0] # K
        for j, T in enumerate(Tini_list):
            p = p_atm*ct.one_atm 
            read_timeHist_hdf = '{}/timeHist_Lu_phi{}_hdf/timeHistory_p{}_T{}.h5'.format(read_dir, phi, p, T)
            print(read_timeHist_hdf)

            gas = ct.Solution(chemPath_cti)
            states = ct.SolutionArray(gas)
            states.read_hdf(read_timeHist_hdf, group='p{}_T{}'.format(p, T))

            # ax[0].plot(states.t, states.T, alpha=1.0, lw=3, ls='-', c=c, label='{} K'.format(T))
            # ax[1].plot(states.t, states.heat_release_rate, lw=3, ls='-', c=c, label='')

            tau1_s = float(df_IDT[df_IDT['T']==T]['tau1_p{}'.format(p_atm)])
            tau2_s = float(df_IDT[df_IDT['T']==T]['tau2_p{}'.format(p_atm)])
            tau_ms = tau2_s*1e+3

            if tau1_s != tau2_s:
                print(phi, T, tau1_s, tau1_s)
                for time, Temp in zip(states.t, states.T):
                    if tau1_s < time:
                        if Temp > 1000:
                            temp_tau1 = np.nan
                            break
                        temp_tau1 = Temp
                        break
            else:
                temp_tau1 = np.nan
            print(T, temp_tau1-T)
            dT_list.append(temp_tau1-T)

        print(Tini_list)
        print(dT_list)
        ax2.plot(Tini_list, dT_list, lw=4, c=c, ls='-')

        df_IDT = read_IDT_hdf('Results_timeHist_IDT/Results_IDT_Lu_phi{}.h5'.format(phi))
        df_IDT1 = df_IDT[(600.0<=df_IDT['T']) & (df_IDT['T']<=850.0)]
        df_IDT2 = df_IDT
        print(df_IDT1)


        df_IDT1.plot(kind='line', ax=ax, x='T', y='tau1_p{}'.format(p_atm), lw=4, c=c, ls='--', label=r'$\phi$={}'.format(phi))
        # df_IDT1.plot(kind='line', ax=ax, x='T_inv', y='tau1_p{}'.format(p), lw=4, c=c, ls='--', label='')
        # df_IDT2.plot(kind='line', ax=ax, x='T_inv', y='tau2_p{}'.format(p), lw=4, c=c, ls='-',  label=r'$n$-heptane/air, $\phi$ = {}'.format(phi))
    ax.set_ylabel(r'Ignition Delay Time [s]')
    ax.legend(loc='upper right', fontsize=18)

    xlim = [600, 800]
    tlim = [1.0e-3, 5.0]
    # ax.legend(loc='lower right', fontsize=16)
    ax.set_xlim(xlim)
    ax.set_ylim(tlim)
    ax.set_yscale('log')
    # ax.set_xlabel(r'1000/T [K$^\mathdefault{-1}$]')
    ax.set_xlabel(r'Temperature [K]')
    ax.set_ylabel(r'Ignition Delay Time [s]')
    ax.text(0.03, 0.9, r'$n$-heptane/air, {} atm'.format(p_atm), transform=ax.transAxes)

    # ticks = ax.get_xticks()
    # ax2.set_xticks(ticks)
    # ax2.set_xticklabels((1000/ticks).round(1))
    # ax2.set_xlim(ax.get_xlim())
    # ax2.set_xlabel(r'Temperature [K]')

    # # T_list = [600, 700, 900, 1000]
    # # clist_line = ['#CB4335', '#2E86C1', '#229954', '#CA6F1E']
    # # for i, T in enumerate(reversed(T_list)):
    # #     ax.vlines(1000.0/T, tlim[0], tlim[1], color=clist_line[i], ls='-.', lw=2, alpha=0.9, zorder=-1)

    # # secax = ax.secondary_xaxis('top', functions=(lambda T:1000/T, lambda Tinv:1000/Tinv))
    # # secax.set_xlabel('Temperature [K]')

    file_name = 'idt_{}atm.png'.format(p)
    file_path = '{0}/{1}'.format(save_dir, file_name)
    fig.tight_layout()
    fig.savefig(file_path, bbox_inches = 'tight')
    subprocess.call(['imgcat', file_path])




    #     if (tmax < tau2_s): tmax = tau2_s
    #     tlim = [0.0, tmax*1.07]
    #     HRRlim_log = [5e+3, 2e+11]
    #     ax[0].text(tau2_s, 1700.0, '{:.0f} K'.format(T), 
    #                ha='center', backgroundcolor='w', c=c, fontsize=20)
    #     if phi == 1.0 and T == 800.0:
    #         ax[1].text(tau2_s, HRRlim_log[1]*4, '{:.1f} ms'.format(tau_ms), 
    #                    ha='center', va='bottom', c=c, fontsize=20)
    #     else:
    #         ax[1].text(tau2_s, HRRlim_log[1], '{:.1f} ms'.format(tau_ms), 
    #                    ha='center', va='bottom', c=c, fontsize=20)

    #     df_tau_table.loc[i] = [phi, p, T, tau1_s, tau2_s]
    # df_tau_table = df_tau_table.sort_values(by=['T [K]'])
    # file_path = '{}/tau_phi{}_{}atm.csv'.format(save_dir, phi, p_atm)
    # df_tau_table.to_csv(file_path, index=False)

    # ax[0].set_title(r'$n$-heptane/air, $\phi$ = {}, {} atm'.format(phi, p_atm), fontsize=24)
    # ax[0].set_ylabel(r'Temperature [K]')
    # ax[1].set_ylabel('Heat release rate [J/m$^3$-s]')
    # ax[1].set_xlabel(r't [s]')

    # ax[1].set_xlim(tlim)
    # ax[1].set_ylim(HRRlim_log)
    # ax[1].set_yscale('log')

    # file_name = 't_THRR_phi{}_{}atm.png'.format(phi, p_atm)
    # file_path = '{0}/{1}'.format(save_dir, file_name)
    # fig.tight_layout()
    # fig.savefig(file_path, bbox_inches = 'tight')
    # subprocess.call(['imgcat', file_path])


