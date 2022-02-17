
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import subprocess
import h5py

import settings_matplotlib as setMPL


def read_SL_hdf(read_file_hdf):

    print(read_file_hdf)
    f_h5 = h5py.File(read_file_hdf, 'r')
    dset_name = list(f_h5.keys())
    T_array = f_h5['/T'][:]
    p_array = f_h5['/p_atm'][:]
    sl_array = f_h5['/SL/SL'][:]
    f_h5.close()
    
    sl = []
    count = 0
    for y in p_array:
        tmp_sl = []
        for x in T_array:
            tmp_sl.append(sl_array[count][4])
            count += 1
        sl.append(tmp_sl)
 
    df = pd.DataFrame()
    df['T'] = T_array
    sl_p_columns = ['sl_p{}'.format(p) for p in p_array]
    df_sl_p = pd.DataFrame(sl, index=sl_p_columns).T
    df = pd.concat([df, df_sl_p], axis=1)

    return df


if __name__ == '__main__':

    save_dir = 'Results_SL'
    read_hdf_dic = {'DTU':'Results_SL/Results_SL_DTU_phi1.0.h5', 'Lu':'Results_SL/Results_SL_Lu_phi1.0.h5'}

    df_DTU = read_SL_hdf(read_hdf_dic['DTU'])
    df_Lu  = read_SL_hdf(read_hdf_dic['Lu'])

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9,6))

    df_DTU.plot(kind='line', ax=ax, x='T', y='sl_p1.0', lw=2, c='r', ls='-', label='DTU')
    df_Lu.plot(kind='line', ax=ax, x='T', y='sl_p1.0', lw=2, c='b', ls='-', label='Lu')

    ax.legend(loc='upper left', fontsize=16)
    ax.set_xlim(290, 900)
    ax.set_ylim(0, 4.0)
    ax.set_xlabel(r'Temperature [K]')
    ax.set_ylabel(r'SL [m/s]')

    file_name = 'sl_dme.png'
    fig.tight_layout()
    fig.savefig("{0}/{1}".format(save_dir, file_name), bbox_inches = 'tight')


