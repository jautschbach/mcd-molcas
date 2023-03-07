import pandas as pd
import numpy as np
import glob
import os

def _read_file(files):
    dfs = []
    for f in files:
        head_lines = 0
        with open(f, 'r') as fn:
            lines = fn.readlines()
            while lines[head_lines][0] != '#':
                head_lines += 1
        df = pd.read_csv(f, index_col=False, skiprows=head_lines+1,
                         header=None, delim_whitespace=True)
        idx = int(f.split('-')[-1])
        df['comp'] = idx
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    df.sort_values(by=['comp', 0], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df

def test_ucl6_f_d():
    test_dir = 'ucl6-f-d'
    test_files = glob.glob(os.path.join(test_dir, 'mcd-c-spectrum-?'))
    actual_files = glob.glob(os.path.join(test_dir, 'mcdspectrum-?'))
    test_df = _read_file(test_files)
    actual_df = _read_file(actual_files)
    assert np.allclose(test_df, actual_df)

def test_ucl6_lmct():
    test_dir = 'ucl6-lmct'
    test_files = glob.glob(os.path.join(test_dir, 'mcd-c-spectrum-?'))
    actual_files = glob.glob(os.path.join(test_dir, 'mcdspectrum-?'))
    test_df = _read_file(test_files)
    actual_df = _read_file(actual_files)
    assert np.allclose(test_df, actual_df)

