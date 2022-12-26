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

    df_IDT = read_IDT_hdf('Results_IDT/Results_IDT_Andrae_phi1.0.h5')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,7))
    ax2 = ax.twiny()

    df_IDT.plot(kind='line', ax=ax, x='T_inv', y='tau1_p1.0', lw=2, c='r', ls='--', label='')
    df_IDT.plot(kind='line', ax=ax, x='T_inv', y='tau2_p1.0', lw=2, c='r', ls='-',  label='Andrae')

    ax.legend(loc='upper left', fontsize=16)
    ax.set_ylim(1.0e-5, 5.0)
    ax.set_yscale('log')
    ax.set_xlabel(r'1000/T [K$^\mathdefault{-1}$]')
    ax.set_ylabel(r'Ignition Delay Time [s]')

    ticks = ax.get_xticks()
    ax2.set_xticks(ticks)
    ax2.set_xticklabels((1000/ticks).round(1))
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xlabel(r'Temperature [K]')

    # secax = ax.secondary_xaxis('top', functions=(lambda T:1000/T, lambda Tinv:1000/Tinv))
    # secax.set_xlabel('Temperature [K]')

    file_name = 'idt.png'
    file_path = '{0}/{1}'.format(save_dir, file_name)
    fig.tight_layout()
    fig.savefig(file_path, bbox_inches = 'tight')
    subprocess.call(['imgcat', file_path])

