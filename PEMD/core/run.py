# Copyright (c) 2024. PEMD developers. All rights reserved.
# Distributed under the terms of the MIT License.

# ******************************************************************************
# Module Docstring
# ******************************************************************************

from PEMD.core.model import PEMDModel
from PEMD.model.build import (
    gen_poly_smiles,
    gen_poly_3D,
)
from PEMD.simulation.qm import (
    conformer_search_xtb,
    conformer_search_gaussian,
    QMcalc_gaussian,
    calc_resp_gaussian
)
from PEMD.simulation.sim_lib import (
    order_energy_xtb,
    read_xyz_file
)

from PEMD.simulation.md import (
    gen_poly_gmx_oplsaa,
)

class ConformerSearch:

    def __init__(self, work_dir, smiles, epsilon, ):
        self.work_dir = work_dir
        self.smiles = smiles
        self.epsilon = epsilon

    def search_by_xtb(self, max_conformers=1000, top_n_MMFF=100, slurm=True, job_name='xtb', nodes=1,
                     ntasks_per_node=32, partition='interactive', ):

        return conformer_search_xtb(
            self.work_dir,
            self.smiles,
            self.epsilon,
            max_conformers,
            top_n_MMFF,
            slurm,
            job_name,
            nodes,
            ntasks_per_node,
            partition,
        )

    def order_energy_xtb(self, numconf, ):

            return order_energy_xtb(
                self.work_dir,
                numconf,
            )

    def search_by_gaussian(self, xyz_file, core, memory, function, basis_set, chg, mult, ):

        sorted_gaussian_file, lowest_energy_structure_file = conformer_search_gaussian(
            xyz_file,
            core,
            memory,
            function,
            basis_set,
            chg,
            mult,
            self.epsilon
        )

        return sorted_gaussian_file, lowest_energy_structure_file





class PEMDSimulation:

    def __init__(self, ):
        pass

    def conformer_search(self, smiles, epsilon, core, memory, function, chg, mult, basis_set, ):

        xyz_file = conformer_search_xtb(
            smiles,
            epsilon,
            max_conformers=1000,
            top_n_MMFF=100,
            top_n_xtb=10,
        )

        sorted_gaussian_file, lowest_energy_structure_file = conformer_search_gaussian(
            xyz_file,
            core,
            memory,
            function,
            basis_set,
            chg,
            mult,
            epsilon
        )

        return sorted_gaussian_file, lowest_energy_structure_file

    def QM_calculation(self, input_xyzfile, core, memory, function, basis_set, chg, mult, epsilon, gaussian_dir, filename):

        structure = read_xyz_file(input_xyzfile)

        return QMcalc_gaussian(
            structure,
            core,
            memory,
            function,
            basis_set,
            chg,
            mult,
            epsilon,
            gaussian_dir,
            filename
        )

    # 生成提交脚本



    # def calc_resp_charge(
    #         self,
    #         epsilon,
    #         core,
    #         memory,
    #         function,
    #         basis_set,
    #         method, # resp1 or resp2
    # ):
    #     sorted_df = self.conformer_search(
    #         epsilon,
    #         core,
    #         memory,
    #         function,
    #         basis_set,
    #     )
    #
    #     return calc_resp_gaussian(
    #         sorted_df,
    #         epsilon,
    #         core,
    #         memory,
    #         method,
    #     )
    #
    # def build_polymer(self,):
    #
    #     return  gen_poly_3D(
    #         self.poly_name,
    #         self.length,
    #         self.gen_poly_smiles(),
    #     )
    #
    # def gen_polymer_force_field(self,):
    #
    #     gen_poly_gmx_oplsaa(
    #         self.poly_name,
    #         self.poly_resname,
    #         self.poly_scale,
    #         self.poly_charge,
    #         self.length,
    #     )








