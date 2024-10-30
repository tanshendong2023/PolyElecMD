# Copyright (c) 2024. PEMD developers. All rights reserved.
# Distributed under the terms of the MIT License.

# ******************************************************************************
# Module Docstring
# ******************************************************************************

import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
from tqdm.auto import tqdm
from statsmodels.tsa.stattools import acovf

def distance(x0, x1, box_length):
    delta = x1 - x0
    delta -= box_length * np.round(delta / box_length)
    return delta

def load_data_traj(work_dir, tpr_file, xtc_wrap_file, xtc_unwrap_file, select_atoms, run_start, run_end, cutoff_radii):
    tpr_filepath = os.path.join(work_dir, tpr_file)
    xtc_wrap_filepath = os.path.join(work_dir, xtc_wrap_file)
    xtc_unwrap_filepath = os.path.join(work_dir, xtc_unwrap_file)

    run_wrap = mda.Universe(tpr_filepath, xtc_wrap_filepath)
    run_unwrap = mda.Universe(tpr_filepath, xtc_unwrap_filepath)

    # Selections
    cations_wrap = run_wrap.select_atoms(select_atoms['cation'])
    cations_unwrap = run_unwrap.select_atoms(select_atoms['cation'])
    anions_wrap = run_wrap.select_atoms(select_atoms['anion'])
    PEO_wrap = run_wrap.select_atoms(select_atoms['PEO'])
    SN_wrap = run_wrap.select_atoms(select_atoms['SN'])

    n_cations = len(cations_wrap)
    total_frames = run_end - run_start
    li_positions = np.zeros((total_frames, n_cations, 3), dtype=np.float64)
    li_env = np.zeros((total_frames, n_cations), dtype=np.int32)

    # Trajectory analysis for positions
    for ts in tqdm(run_unwrap.trajectory[run_start:run_end], desc='Processing positions'):
        li_positions[ts.frame - run_start] = cations_unwrap.positions.copy()

    # Trajectory analysis for environments
    for ts in tqdm(run_wrap.trajectory[run_start:run_end], desc='Processing environments'):
        box_size = ts.dimensions[0]

        li_pos = cations_wrap.positions.copy()
        PEO_pos = PEO_wrap.positions.copy()
        anions_pos = anions_wrap.positions.copy()
        SN_pos = SN_wrap.positions.copy()

        for n in range(n_cations):
            li = li_pos[n]

            distances_oe = np.linalg.norm(distance(PEO_pos, li, box_size), axis=1)
            close_oe = np.any(distances_oe < cutoff_radii['PEO'])

            distances_tfsi = np.linalg.norm(distance(anions_pos, li, box_size), axis=1)
            close_tfsi = np.any(distances_tfsi < cutoff_radii['TFSI'])

            distances_sn = np.linalg.norm(distance(SN_pos, li, box_size), axis=1)
            close_sn = np.any(distances_sn < cutoff_radii['SN'])

            # Assign environment type
            if close_oe and not (close_tfsi or close_sn):
                li_env[ts.frame - run_start, n] = 1
            elif close_oe and (close_tfsi or close_sn):
                li_env[ts.frame - run_start, n] = 2
            elif not close_oe and (close_tfsi or close_sn):
                li_env[ts.frame - run_start, n] = 3

    return run_wrap, cations_wrap, li_positions, li_env

def load_md_trajectory(work_dir, tpr_filename='nvt_prod.tpr', xtc_filename='nvt_prod.xtc'):
    data_tpr_file = os.path.join(work_dir, tpr_filename)
    data_xtc_file = os.path.join(work_dir, xtc_filename)
    u = mda.Universe(data_tpr_file, data_xtc_file)
    return u

def times_array(run, run_start, run_end, time_step=5):
    times = []
    for step, _ts in enumerate(run.trajectory[run_start:run_end]):
        times.append(step * time_step)
    return np.array(times)

def calc_acf(a_values: dict[str, np.ndarray]) -> list[np.ndarray]:
    acfs = []
    for neighbors in a_values.values():  # for _atom_id, neighbors in a_values.items():
        acfs.append(acovf(neighbors, demean=False, adjusted=True, fft=True))
    return acfs

def calc_neigh_corr(run, center_atoms, distance_dict, select_dict, run_start, run_end):
    acf_avg = {}
    # center_atoms = run.select_atoms('resname LIP and name Li')
    for kw in distance_dict:
        acf_all = []
        for atom in tqdm(center_atoms[::]):
            distance = distance_dict.get(kw)
            assert distance is not None
            bool_values = {}
            for time_count, _ts in enumerate(run.trajectory[run_start:run_end:]):
                if kw in select_dict:
                    selection = (
                            "("
                            + select_dict[kw]
                            + ") and (around "
                            + str(distance)
                            + " index "
                            + str(atom.id - 1)
                            + ")"
                    )
                    shell = run.select_atoms(selection)
                else:
                    raise ValueError("Invalid species selection")
                # 获取这些原子所属的分子ID
                mols = set(atom.residue for atom in shell)

                for mol in mols:
                    if str(mol.resid) not in bool_values:
                        bool_values[str(mol.resid)] = np.zeros(int((run_end - run_start) / 1))
                    bool_values[str(mol.resid)][time_count] = 1

            acfs = calc_acf(bool_values)
            acf_all.extend(list(acfs))
        acf_avg[kw] = np.mean(acf_all, axis=0)
    return acf_avg


if __name__ == '__main__':
    work_dir = './'

    tpr_file = 'nvt_prod.tpr'
    xtc_wrap_file = 'nvt_prod.xtc'
    xtc_unwrap_file = 'nvt_prod_unwrap.xtc'

    select_atoms = {
        'cation': 'resname LIP and name Li',
        'anion': 'resname NSC and name OBT',
        'PEO': 'resname MOL and name O',
        'SN': 'resname SN and name N',
    }

    cutoff_radii = {
        'PEO': 3.7,
        'TFSI': 3.5,
        'SN': 3.7,
    }

    run_start = 0
    run_end = 10001  # frames
    time_step = 0.001  # ps per MD step
    dt_collection = 5000  # MD steps between frames
    time_per_frame = time_step * dt_collection  # ps per frame

    # 加载数据并获取满足条件的锂离子
    run_wrap, cations_wrap, li_positions, li_env = load_data_traj(
        work_dir,
        tpr_file,
        xtc_wrap_file,
        xtc_unwrap_file,
        select_atoms,
        run_start,
        run_end,
        cutoff_radii,
    )

    fraction_env1 = np.mean(li_env == 1, axis=0)
    cation_indices = np.where(fraction_env1 >= 0.8)[0]
    center_atoms = cations_wrap[cation_indices]

    run = load_md_trajectory(work_dir, tpr_file, xtc_wrap_file)

    times = times_array(run, run_start, run_end, time_step=5)

    # 定义距离和选择字典
    distance_dict = {"anion": 3.5, 'SN':3.7}
    select_dict = {
        "cation": "resname LIP and name Li",
        "anion": "resname NSC and name OBT",
        "SN": 'resname SN and name N'
    }

    # 计算邻居自相关函数
    acf_avg = calc_neigh_corr(run, center_atoms, distance_dict, select_dict, run_start, run_end)

    # 归一化自相关函数
    acf_avg_norm = {}
    species_list = list(acf_avg.keys())
    for kw in species_list:
        if kw in acf_avg:
            acf_avg_norm[kw] = acf_avg[kw] / acf_avg[kw][0]  # 防止除以零的错误处理

    # 准备将时间和归一化后的自相关函数保存在同一个CSV文件中
    acf_df = pd.DataFrame({'Time (ps)': times})
    for key, value in acf_avg_norm.items():
        acf_df[key + ' ACF'] = pd.Series(value)

    # 将数据保存到CSV文件
    acf_df.to_csv('residence_time.csv', index=False)


