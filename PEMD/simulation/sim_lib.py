# Copyright (c) 2024. PEMD developers. All rights reserved.
# Distributed under the terms of the MIT License.

# ******************************************************************************
# Module Docstring
# ******************************************************************************

import os
import time
import subprocess
import numpy as np
from rdkit import Chem
from openbabel import openbabel as ob
from simple_slurm import Slurm
from PEMD.model import model_lib
import PEMD.model.MD_lib as MDlib
from pysimm import system, lmps, forcefield

# OpenBabel setup
obConversion = ob.OBConversion()
ff = ob.OBForceField.FindForceField('UFF')
mol = ob.OBMol()
np.set_printoptions(precision=20)

def get_slurm_job_status(job_id):
    command = f'sacct -j {job_id} --format=State --noheader'
    process = subprocess.run(command, shell=True, capture_output=True, text=True)
    # Split the output by newlines and strip whitespace
    statuses = [line.strip() for line in process.stdout.strip().split('\n')]
    # Check if all statuses indicate the job is completed
    if all(status == 'COMPLETED' for status in statuses):
        return 'COMPLETED'
    elif any(status == 'FAILED' for status in statuses):
        return 'FAILED'
    elif any(status == 'CANCELLED' for status in statuses):
        return 'CANCELLED'
    else:
        return 'RUNNING'

# Modified order_energy_xtb function
def order_energy_xtb(file_path, numconf, output_file):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    structures = []
    current_structure = []
    is_first_line = True  # Indicates if we are at the start of the file

    for line in lines:
        line = line.strip()
        if line.isdigit() and (is_first_line or current_structure):
            if current_structure and not is_first_line:
                energy_line = current_structure[1]  # Second line contains energy information
                try:
                    energy = float(energy_line.split()[-1])
                except ValueError:
                    print(f"Could not parse energy value: {energy_line}")
                    energy = float('inf')  # Assign infinite energy if parsing fails
                structures.append((energy, current_structure))
            current_structure = [line]  # Start a new structure
            is_first_line = False
        else:
            current_structure.append(line)

    # Add the last structure
    if current_structure:
        energy_line = current_structure[1]
        try:
            energy = float(energy_line.split()[-1])
        except ValueError:
            print(f"Could not parse energy value: {energy_line}")
            energy = float('inf')
        structures.append((energy, current_structure))

    # Sort structures by energy
    structures.sort(key=lambda x: x[0])

    # Select the lowest energy structures
    selected_structures = structures[:numconf]

    # Write the selected structures to the output .xyz file
    with open(output_file, 'w') as outfile:
        for idx, (energy, structure) in enumerate(selected_structures):
            for line_num, line in enumerate(structure):
                if line_num == 1:
                    # Modify the comment line to include the energy value
                    outfile.write(f"Energy = {energy}\n")
                else:
                    outfile.write(f"{line}\n")

    print(f"The lowest {numconf} energy structures have been written to {output_file}")

def read_xyz_file(file_path):
    """
    读取 .xyz 文件，返回一个结构列表。
    每个结构是一个包含原子数、注释行和原子坐标的字典。
    """
    structures = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        num_atoms_line = lines[i].strip()
        if num_atoms_line.isdigit():
            num_atoms = int(num_atoms_line)
            comment_line = lines[i + 1].strip()
            atoms = []
            for j in range(i + 2, i + 2 + num_atoms):
                atom_line = lines[j].strip()
                atoms.append(atom_line)
            structure = {
                'num_atoms': num_atoms,
                'comment': comment_line,
                'atoms': atoms
            }
            structures.append(structure)
            i = i + 2 + num_atoms
        else:
            i += 1
    return structures

def submit_job_to_slurm(command, job_name, node, core, mem, gaussian_dir, file, soft):

    file_contents = f"#!/bin/bash"
    file_contents += f"#SBATCH -J {job_name}\n"
    file_contents += f"#SBATCH -N {node}\n"
    file_contents += f"#SBATCH -n {core}\n"
    file_contents += f"#SBATCH -p {mem}\n\n"
    file_contents += f"module load {soft}\n\n"
    file_contents += f"{command} {file}/\n"

    # 将作业脚本写入文件
    script_path = os.path.join(gaussian_dir, f"sub.sh")
    with open(script_path, 'w') as f:
        f.write(file_contents)

    # 提交作业并获取作业 ID
    submit_command = ['sbatch', script_path]
    result = subprocess.run(submit_command, capture_output=True, text=True)
    if result.returncode == 0:
        # 提取作业 ID
        output = result.stdout.strip()
        # 通常，sbatch 的输出格式为 "Submitted batch job <job_id>"
        job_id = output.split()[-1]
        return job_id
    else:
        print(f"提交作业失败：{result.stderr}")
        return None

def read_energy_from_gaussian(log_file_path):
    """
    从 Gaussian 输出文件中读取能量（自由能）
    """
    with open(log_file_path, 'r') as file:
        lines = file.readlines()
    energy = None
    for line in lines:
        if 'Sum of electronic and thermal Free Energies=' in line:
            energy = float(line.strip().split()[-1])
    return energy

def read_final_structure_from_gaussian(log_file_path):
    """
    从 Gaussian 输出文件中提取优化后的结构坐标
    """
    with open(log_file_path, 'r') as file:
        lines = file.readlines()
    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        if 'Standard orientation:' in line:
            start_idx = i + 5  # 跳过标题行
        elif start_idx and '---------------------------------------------------------------------' in line:
            end_idx = i
            break
    if start_idx and end_idx:
        atoms = []
        for line in lines[start_idx:end_idx]:
            tokens = line.strip().split()
            atom_number = int(tokens[1])
            x = float(tokens[3])
            y = float(tokens[4])
            z = float(tokens[5])
            atom_symbol = Chem.PeriodicTable.GetElementSymbol(Chem.GetPeriodicTable(), atom_number)
            atoms.append(f"{atom_symbol}   {x}   {y}   {z}")
        return atoms
    else:
        return None

def order_energy_gaussian(dir_path, output_file):
    data = []
    # Traverse all files in the specified folder
    for file in os.listdir(dir_path):
        if file.endswith(".out"):
            log_file_path = os.path.join(dir_path, file)
            energy = read_energy_from_gaussian(log_file_path)
            atoms = read_final_structure_from_gaussian(log_file_path)
            if energy is not None and atoms is not None:
                data.append({"Energy": energy, "Atoms": atoms})

    # Check if data is not empty
    if data:
        # Sort the structures by energy
        sorted_data = sorted(data, key=lambda x: x['Energy'])
        # Write the sorted structures to an .xyz file
        with open(output_file, 'w') as outfile:
            for item in sorted_data:
                num_atoms = len(item['Atoms'])
                outfile.write(f"{num_atoms}\n")
                outfile.write(f"Energy = {item['Energy']}\n")
                for atom_line in item['Atoms']:
                    outfile.write(f"{atom_line}\n")
        print(f"Sorted structures have been saved to {output_file}")
    else:
        print(f"No successful Gaussian output files found in {dir_path}")

def get_gaff2(unit_name, length, relax_polymer_lmp_dir, mol, atom_typing='pysimm'):
    print("\nGenerating GAFF2 parameter file ...\n")
    # r = MDlib.get_coord_from_pdb(outfile_name + ".pdb")
    # from pysimm import system, forcefield

    file_base = relax_polymer_lmp_dir + '/' + f'{unit_name}_N{length}'

    obConversion.SetInAndOutFormats("pdb", "cml")
    if os.path.exists(file_base + '.pdb'):
        mol = ob.OBMol()
        obConversion.ReadFile(mol, file_base + '.pdb')
    else:
        try:
            count_atoms = mol.NumAtoms()
            if count_atoms > 0:
                pass
            else:
                print("ERROR: Number of atoms = ", count_atoms)
                print("Couldn't generate GAFF2 parameter file\n")
                return
        except BaseException:
            print("ERROR: pdb file not found; OBMol not provided")
            print("Couldn't generate GAFF2 parameter file\n")
            return

    obConversion.WriteFile(mol, file_base + '.cml')
    data_fname = file_base + '_gaff2.lmp'

    try:
        print("Pysimm working on {}".format(file_base + '.mol2'))
        s = system.read_cml(file_base + '.cml')

        f = forcefield.Gaff2()
        if atom_typing == 'pysimm':
            for b in s.bonds:
                if b.a.bonds.count == 3 and b.b.bonds.count == 3:
                    b.order = 4
            s.apply_forcefield(f, charges='gasteiger')
        elif atom_typing == 'antechamber':
            obConversion.SetOutFormat("mol2")
            obConversion.WriteFile(mol, file_base + '.mol2')
            print("Antechamber working on {}".format(file_base + '.mol2'))
            MDlib.get_type_from_antechamber(s, file_base + '.mol2', 'gaff2', f)
            s.pair_style = 'lj'
            s.apply_forcefield(f, charges='gasteiger', skip_ptypes=True)
        else:
            print('Invalid atom typing option, please select pysimm or antechamber.')
        s.write_lammps(data_fname)
        print("\nGAFF2 parameter file generated.")
    except BaseException:
        print('problem reading {} for Pysimm.'.format(file_base + '.cml'))


def relax_polymer_lmp(unit_name, length, relax_polymer_lmp_dir, core):
    origin_dir = os.getcwd()
    os.chdir(relax_polymer_lmp_dir)
    # 创建LAMMPS输入文件字符串
    file_base = f'{unit_name}_N{length}'
    lmp_commands = """
    units real
    boundary s s s
    dimension 3
    atom_style full
    bond_style harmonic
    angle_style harmonic
    dihedral_style fourier
    improper_style harmonic
    pair_style lj/cut 2.0
    read_data {0}_gaff2.lmp
    thermo 100
    thermo_style custom step temp pxx pyy pzz ebond eangle edihed eimp epair ecoul evdwl pe ke etotal lx ly lz vol density
    fix xwalls all wall/reflect xlo EDGE xhi EDGE
    fix ywalls all wall/reflect ylo EDGE yhi EDGE
    fix zwalls all wall/reflect zlo EDGE zhi EDGE
    velocity all create 300 3243242
    minimize 1e-8 1e-8 10000 10000
    dump 1 all custom 500 soft.lammpstrj id type mass mol x y z
    fix 1 all nvt temp 800 800 100 drag 2
    run 5000000
    write_data {0}_gaff2.data
    write_dump all xyz {0}_lmp.xyz
    """.format(file_base)

    # 将LAMMPS输入命令写入临时文件
    with open('lammps_input.in', "w") as f:
        f.write(lmp_commands)

    slurm = Slurm(J='lammps',
                  N=1,
                  n=f'{core}',
                  output=f'slurm.{Slurm.JOB_ARRAY_MASTER_ID}.out'
                  )

    slurm.add_cmd('module load LAMMPS')
    job_id = slurm.sbatch('mpirun lmp < lammps_input.in >out.lmp 2>lmp.err')
    while True:
        status = get_slurm_job_status(job_id)
        if status in ['COMPLETED', 'FAILED', 'CANCELLED']:
            print("\n", unit_name, ": MD simulation normally terminated.\n")
            model_lib.toxyz_lammps(f'{file_base}_lmp.xyz', f'{file_base}_gmx.xyz', f'{file_base}_gaff2.lmp')
            os.chdir(origin_dir)
            break
        else:
            print("polymer relax not finish, waiting...")
            time.sleep(10)

