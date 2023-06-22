
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
    width_array = f_h5['/width'][:]
    phi_array   = f_h5['/phi'][:]
    T_array     = f_h5['/T'][:]
    p_array     = f_h5['/p_atm'][:]
    sl_array    = f_h5['/SL/SL'][:]
    f_h5.close()

    df = pd.DataFrame(sl_array, \
                      columns=['width', 'phi', 'p', 'p_atm', 'T', 'SL'])
    index_failure = df[df['SL'] == 0].index
    df.loc[index_failure,'SL'] = None

    return T_array, df


def graph_phi_SL():

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9,6))

    T_array = np.array([300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])
    for T in np.flip(T_array):
        df_graph = df_Lu[df_Lu['width']==width]
        df_graph = df_graph[df_graph['T']==T]
        df_graph.plot(kind='line', ax=ax, x='phi', y='SL', lw=2, ls='-', marker='.', label=f'{T} K')

    ax.grid(c='lightgray', ls=':')
    ax.legend(loc='upper right', fontsize=16, frameon=False)
    ax.set_xlim(0.5, 2.0)
    ax.set_ylim(limSL)
    ax.set_xlabel(r'Equivalence ratio, $\phi$')
    ax.set_ylabel(r'Laminar flame speed, S$_L$ [m/s]')

    file_name = 'phi_sl_width{}_c7.png'.format(width)
    save_path = '{0}/{1}'.format(save_dir, file_name)
    fig.tight_layout()
    fig.savefig(save_path, bbox_inches = 'tight')
    # subprocess.call(['imgcat', save_path])


def graph_T_SL():

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9,6))

    phi_array = np.array([0.6, 0.8, 1.0, 1.5, 2.0])
    for phi in phi_array:
        df_graph = df_Lu[df_Lu['width']==width]
        df_graph = df_graph[df_graph['phi']==phi]
        df_graph.plot(kind='line', ax=ax, x='T', y='SL', lw=2, ls='-', marker='.', label=f'$\phi$ = {phi}')

    ax.grid(c='lightgray', ls=':')
    ax.legend(loc='upper left', fontsize=16, frameon=False)
    ax.set_xlim(300.0, None)
    ax.set_ylim(limSL)
    ax.set_xlabel(r'Initial temperature [K]')
    ax.set_ylabel(r'Laminar flame speed, S$_L$ [m/s]')

    file_name = 'T_sl_width{}_c7.png'.format(width)
    save_path = '{0}/{1}'.format(save_dir, file_name)
    fig.tight_layout()
    fig.savefig(save_path, bbox_inches = 'tight')
    # subprocess.call(['imgcat', save_path])


def graph_T_normSL():

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9,6))

    phi_array = np.array([0.6, 0.8, 1.0, 1.5, 2.0])
    for phi in phi_array:
        df_graph = df_Lu[df_Lu['width']==width]
        df_graph = df_graph[df_graph['phi']==phi]
        SL_300K = float(df_graph[df_graph['T']==300.0]['SL'])
        df_graph['normSL'] = df_graph['SL']/SL_300K
        df_graph.plot(kind='line', ax=ax, x='T', y='normSL', lw=2, ls='-', marker='.', label=f'$\phi$ = {phi}')

    ax.grid(c='lightgray', ls=':')
    ax.legend(loc='upper left', fontsize=16, frameon=False)
    ax.set_xlim(300.0, None)
    ax.set_ylim(1.0, None)
    ax.set_xlabel(r'Initial temperature [K]')
    ax.set_ylabel(r'Norm. laminar flame speed [-]')

    file_name = 'T_normsl_width{}_c7.png'.format(width)
    save_path = '{0}/{1}'.format(save_dir, file_name)
    fig.tight_layout()
    fig.savefig(save_path, bbox_inches = 'tight')
    # subprocess.call(['imgcat', save_path])


if __name__ == '__main__':

    read_hdf = 'Results_SL/Results_SL_Lu_width2mm.h5'
    read_hdf = 'Results_SL/Results_SL_Lu.h5'
    save_dir = 'Results_SL'
    limSL = [0.0, 5.0]
    T_array, df_Lu = read_SL_hdf(read_hdf)

    width_list = [0.003]
    width_list = [0.02, 0.01, 0.005, 0.002, 0.001, 0.0005]
    for width in width_list:
        graph_phi_SL()
        graph_T_SL()
        graph_T_normSL()


