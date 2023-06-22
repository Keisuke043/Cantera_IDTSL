import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import subprocess
import h5py

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

    save_dir = 'Results_IDT'

    p = 1.0 # atm
    phi_list = [0.8, 1.0, 1.4, 2.0]
    clist = ['violet', 'cornflowerblue', 'burlywood', 'darkseagreen']

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,7))
    ax2 = ax.twiny()

    for i, phi in enumerate(phi_list):
        c = clist[i]

        df_IDT = read_IDT_hdf('Results_IDT/Results_IDT_Lu_phi{}.h5'.format(phi))
        df_IDT1 = df_IDT[df_IDT['T']<1200.0]
        df_IDT2 = df_IDT
        df_IDT1.plot(kind='line', ax=ax, x='T_inv', y='tau1_p{}'.format(p), lw=4, c=c, ls='--', label='')
        df_IDT2.plot(kind='line', ax=ax, x='T_inv', y='tau2_p{}'.format(p), lw=4, c=c, ls='-',  label=r'$n$-heptane/air, $\phi$ = {}'.format(phi))

    tlim = [1.0e-5, 2.0]
    ax.legend(loc='lower right', fontsize=16)
    ax.set_ylim(tlim)
    ax.set_yscale('log')
    ax.set_xlabel(r'1000/T [K$^\mathdefault{-1}$]')
    ax.set_ylabel(r'Ignition Delay Time [s]')
    ax.text(0.03, 0.9, r'$n$-heptane/air, {} atm'.format(p), transform=ax.transAxes)

    ticks = ax.get_xticks()
    ax2.set_xticks(ticks)
    ax2.set_xticklabels((1000/ticks).round(1))
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xlabel(r'Temperature [K]')

    # T_list = [600, 700, 900, 1000]
    # clist_line = ['#CB4335', '#2E86C1', '#229954', '#CA6F1E']
    # for i, T in enumerate(reversed(T_list)):
    #     ax.vlines(1000.0/T, tlim[0], tlim[1], color=clist_line[i], ls='-.', lw=2, alpha=0.9, zorder=-1)

    # secax = ax.secondary_xaxis('top', functions=(lambda T:1000/T, lambda Tinv:1000/Tinv))
    # secax.set_xlabel('Temperature [K]')

    file_name = 'idt_{}atm.png'.format(p)
    file_path = '{0}/{1}'.format(save_dir, file_name)
    fig.tight_layout()
    fig.savefig(file_path, bbox_inches = 'tight')
    subprocess.call(['imgcat', file_path])

