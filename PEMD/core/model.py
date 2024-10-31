# Copyright (c) 2024. PEMD developers. All rights reserved.
# Distributed under the terms of the MIT License.

# ******************************************************************************
# core.model Module
# ******************************************************************************


import os
import json
from PEMD.model.build import (
    gen_poly_smiles,
    gen_poly_3D,
)


class PEMDModel:
    def __init__(self, poly_name, repeating_unit, leftcap, rightcap, length_short, length_poly, ):
        self.poly_name = poly_name
        self.repeating_unit = repeating_unit
        self.leftcap = leftcap
        self.rightcap = rightcap
        self.length_short = length_short
        self.length_poly = length_poly

    @classmethod
    def from_json(cls, work_dir, json_file):

        json_path = os.path.join(work_dir, json_file)
        with open(json_path, 'r', encoding='utf-8') as file:
            model_info = json.load(file)

        poly_name = model_info['polymer']['compound']
        repeating_unit = model_info['polymer']['repeating_unit']
        leftcap = model_info['polymer']['left_cap']
        rightcap = model_info['polymer']['right_cap']
        length_short = model_info['polymer']['length'][0]
        length_poly = model_info['polymer']['length'][1]

        return cls(poly_name, repeating_unit, leftcap, rightcap, length_short, length_poly)

    def gen_poly_smiles(self, short=False):

        if short:
            return gen_poly_smiles(
                self.poly_name,
                self.repeating_unit,
                self.length_short,
                self.leftcap,
                self.rightcap,
            )
        else:
            return gen_poly_smiles(
                self.poly_name,
                self.repeating_unit,
                self.leftcap,
                self.rightcap,
                self.length_poly,
            )

    # def build_polymer(self,):
    #
    #     return  gen_poly_3D(
    #         self.poly_name,
    #         self.length,
    #         self.gen_poly_smiles(),
    #     )
















