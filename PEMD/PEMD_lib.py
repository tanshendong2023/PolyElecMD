import osimport numpy as npimport pandas as pdfrom rdkit import Chemfrom rdkit.Chem import AllChemfrom openbabel import openbabel as obfrom LigParGenPEMD import Converterimport PEMD.MD_lib as MDlib# OpenBabel setupobConversion = ob.OBConversion()ff = ob.OBForceField.FindForceField('UFF')mol = ob.OBMol()np.set_printoptions(precision=20)def rdkitmol2xyz(unit_name, m, dir_xyz, IDNum):    try:        Chem.MolToXYZFile(m, dir_xyz + unit_name + '.xyz', confId=IDNum)    except Exception:        obConversion.SetInAndOutFormats("mol", "xyz")        Chem.MolToMolFile(m, dir_xyz + unit_name + '.mol', confId=IDNum)        mol = ob.OBMol()        obConversion.ReadFile(mol, dir_xyz + unit_name + '.mol')        obConversion.WriteFile(mol, dir_xyz + unit_name + '.xyz')# This function create XYZ files from SMILES# INPUT: ID, SMILES, directory name# OUTPUT: xyz files in 'work_dir', result = DONE/NOT DONE, mol without Hydrogen atomdef smiles_xyz(unit_name, SMILES, dir_xyz):    try:        m1 = Chem.MolFromSmiles(SMILES)    # Get mol(m1) from smiles        m2 = Chem.AddHs(m1)   # Add H        AllChem.Compute2DCoords(m2)    # Get 2D coordinates        AllChem.EmbedMolecule(m2)    # Make 3D mol        m2.SetProp("_Name", unit_name + '   ' + SMILES)    # Change title        AllChem.UFFOptimizeMolecule(m2, maxIters=200)    # Optimize 3D str        rdkitmol2xyz(unit_name, m2, dir_xyz, -1)        result = 'DONE'    except Exception:        result, m1 = 'NOT_DONE', ''    return result, m1# This function indentifies row numbers of dummy atoms# INPUT: SMILES# OUTPUT: row indices of dummy atoms and nature of bond with connecting atomdef FetchDum(smiles):    m = Chem.MolFromSmiles(smiles)    dummy_index = []    if m is not None:        for atom in m.GetAtoms():            if atom.GetSymbol() == '*':                dummy_index.append(atom.GetIdx())        for bond in m.GetBonds():            if (                bond.GetBeginAtom().GetSymbol() == '*'                or bond.GetEndAtom().GetSymbol() == '*'            ):                bond_type = bond.GetBondType()                break    return dummy_index, str(bond_type)# Connection information obtained by OpenBabel# INPUT: XYZ file# OUTPUT: Connectivity informationdef connec_info(unit_name):    obConversion = ob.OBConversion()    obConversion.SetInFormat("xyz")    mol = ob.OBMol()    obConversion.ReadFile(mol, unit_name)    neigh_atoms_info = []    for atom in ob.OBMolAtomIter(mol):        neigh_atoms = []        bond_orders = []        for allatom in ob.OBAtomAtomIter(atom):            neigh_atoms.append(allatom.GetIndex())            bond_orders.append(atom.GetBond(allatom).GetBondOrder())        neigh_atoms_info.append([neigh_atoms, bond_orders])    neigh_atoms_info = pd.DataFrame(neigh_atoms_info, columns=['NeiAtom', 'BO'])    return neigh_atoms_infodef Init_info(unit_name, smiles_each_ori,length):    # Get index of dummy atoms and bond type associated with it    try:        dum_index, bond_type = FetchDum(smiles_each_ori)        if len(dum_index) == 2:            dum1 = dum_index[0]            dum2 = dum_index[1]        else:            print(                unit_name,                ": There are more or less than two dummy atoms in the SMILES string; "                "Hint: PSP works only for one-dimensional polymers.",            )            return unit_name, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'REJECT'    except Exception:        print(            unit_name,            ": Couldn't fetch the position of dummy atoms. Hints: (1) In SMILES strings, use '*' for a dummy atom,"            "(2) Check RDKit installation.",        )        return unit_name, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'REJECT'    # Assign dummy atom according to bond type    if bond_type == 'SINGLE':        dum, unit_dis = 'Cl', -0.17        # List of oligomers        oligo_list = list(set(length) - set(['n']))    elif bond_type == 'DOUBLE':        dum, unit_dis = 'O', 0.25        # List of oligomers        oligo_list = []    else:        print(            unit_name,            ": Unusal bond type (Only single or double bonds are acceptable)."            "Hints: (1) Check bonds between the dummy and connecting atoms in SMILES string"            "       (2) Check RDKit installation.",        )        return unit_name, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'REJECT'    # Replace '*' with dummy atom    smiles_each = smiles_each_ori.replace(r'*', dum)    # Convert SMILES to XYZ coordinates    convert_smiles2xyz, m1 = smiles_xyz(unit_name, smiles_each, './')    # if fails to get XYZ coordinates; STOP    if convert_smiles2xyz == 'NOT_DONE':        print(            unit_name,            ": Couldn't get XYZ coordinates from SMILES string. Hints: (1) Check SMILES string,"            "(2) Check RDKit installation.",        )        return unit_name, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'REJECT'    # Collect valency and connecting information for each atom    neigh_atoms_info = connec_info('./' + unit_name + '.xyz')    try:        # Find connecting atoms associated with dummy atoms.        # dum1 and dum2 are connected to atom1 and atom2, respectively.        atom1 = neigh_atoms_info['NeiAtom'][dum1].copy()[0]        atom2 = neigh_atoms_info['NeiAtom'][dum2].copy()[0]    except Exception:        print(            unit_name,            ": Couldn't get the position of connecting atoms. Hints: (1) XYZ coordinates are not acceptable,"            "(2) Check Open Babel installation.",        )        return unit_name, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'REJECT'    return (        unit_name,        dum1,        dum2,        atom1,        atom2,        m1,        neigh_atoms_info,        oligo_list,        dum,        unit_dis,        '',    )def gen_oligomer_smiles(    unit_name,    dum1,    dum2,    atom1,    atom2,    input_smiles,    ln,    # loop,    smiles_LCap_,    LCap_,    smiles_RCap_,    RCap_,):    input_mol = Chem.MolFromSmiles(input_smiles)    edit_m1 = Chem.EditableMol(input_mol)    edit_m1.RemoveAtom(dum1)    if dum1 < dum2:        edit_m1.RemoveAtom(dum2 - 1)    else:        edit_m1.RemoveAtom(dum2)    monomer_mol = edit_m1.GetMol()    inti_mol = monomer_mol    if atom1 > atom2:        atom1, atom2 = atom2, atom1    if dum1 < atom1 and dum2 < atom1:        second_atom = atom1 - 2    elif (dum1 < atom1 and dum2 > atom1) or (dum1 > atom1 and dum2 < atom1):        second_atom = atom1 - 1    else:        second_atom = atom1    if dum1 < atom2 and dum2 < atom2:        first_atom = atom2 - 2    elif (dum1 < atom2 and dum2 > atom2) or (dum1 > atom2 and dum2 < atom2):        first_atom = atom2 - 1    else:        first_atom = atom2    for i in range(1, ln):        combo = Chem.CombineMols(inti_mol, monomer_mol)        edcombo = Chem.EditableMol(combo)        edcombo.AddBond(            second_atom + (i - 1) * monomer_mol.GetNumAtoms(),            first_atom + i * monomer_mol.GetNumAtoms(),            order=Chem.rdchem.BondType.SINGLE,        )        inti_mol = edcombo.GetMol()    # if loop is True and LCap_ is False and RCap_ is False:    #     edcombo.AddBond(    #         first_atom,    #         second_atom + i * monomer_mol.GetNumAtoms(),    #         order=Chem.rdchem.BondType.SINGLE,    #     )    #     inti_mol = edcombo.GetMol()    if LCap_ is True or RCap_ is True:        inti_mol = gen_smiles_with_cap(            unit_name,            0,            0,            first_atom,            second_atom + i * monomer_mol.GetNumAtoms(),            inti_mol,            smiles_LCap_,            smiles_RCap_,            LCap_,            RCap_,            WithDum=False,        )        return inti_mol    return Chem.MolToSmiles(inti_mol)def gen_smiles_with_cap(    unit_name,    dum1,    dum2,    atom1,    atom2,    smiles_each,    smiles_LCap_,    smiles_RCap_,    LCap_,    RCap_,    WithDum=True,):    # Main chain    # Check if there are dummy atoms in the chain    if WithDum is True:        main_mol = Chem.MolFromSmiles(smiles_each)        main_edit_m1 = Chem.EditableMol(main_mol)        # Remove dummy atoms        main_edit_m1.RemoveAtom(dum1)        if dum1 < dum2:            main_edit_m1.RemoveAtom(dum2 - 1)        else:            main_edit_m1.RemoveAtom(dum2)        # Mol without dummy atom        main_mol_noDum = main_edit_m1.GetMol()        # Get linking atoms        if atom1 > atom2:            atom1, atom2 = atom2, atom1        if dum1 < atom1 and dum2 < atom1:            first_atom = atom1 - 2        elif (dum1 < atom1 and dum2 > atom1) or (dum1 > atom1 and dum2 < atom1):            first_atom = atom1 - 1        else:            first_atom = atom1        if dum1 < atom2 and dum2 < atom2:            second_atom = atom2 - 2        elif (dum1 < atom2 and dum2 > atom2) or (dum1 > atom2 and dum2 < atom2):            second_atom = atom2 - 1        else:            second_atom = atom2    else:        main_mol_noDum = smiles_each        first_atom, second_atom = atom1, atom2    LCap_add = 0    # Left Cap    if LCap_ is True:        (unit_name, dum_L, atom_L, m1L, neigh_atoms_info_L, flag_L) = Init_info_Cap(            unit_name, smiles_LCap_        )        # Reject if SMILES is not correct        if flag_L == 'REJECT':            return unit_name, 'REJECT', 0        # Editable Mol for LeftCap        LCap_m1 = Chem.MolFromSmiles(smiles_LCap_)        LCap_edit_m1 = Chem.EditableMol(LCap_m1)        # Remove dummy atoms        LCap_edit_m1.RemoveAtom(dum_L)        # Mol without dummy atom        LCap_m1 = LCap_edit_m1.GetMol()        LCap_add = LCap_m1.GetNumAtoms()        # Linking atom        if dum_L < atom_L:            LCap_atom = atom_L - 1        else:            LCap_atom = atom_L        # Join main chain with Left Cap        combo = Chem.CombineMols(LCap_m1, main_mol_noDum)        edcombo = Chem.EditableMol(combo)        edcombo.AddBond(            LCap_atom, first_atom + LCap_add, order=Chem.rdchem.BondType.SINGLE        )        main_mol_noDum = edcombo.GetMol()    # Right Cap    if RCap_ is True:        (unit_name, dum_R, atom_R, m1L, neigh_atoms_info_R, flag_R) = Init_info_Cap(            unit_name, smiles_RCap_        )        # Reject if SMILES is not correct        if flag_R == 'REJECT':            return unit_name, 'REJECT', 0        # Editable Mol for RightCap        RCap_m1 = Chem.MolFromSmiles(smiles_RCap_)        RCap_edit_m1 = Chem.EditableMol(RCap_m1)        # Remove dummy atoms        RCap_edit_m1.RemoveAtom(dum_R)        # Mol without dummy atom        RCap_m1 = RCap_edit_m1.GetMol()        # Linking atom        if dum_R < atom_R:            RCap_atom = atom_R - 1        else:            RCap_atom = atom_R        # Join main chain with Left Cap        combo = Chem.CombineMols(main_mol_noDum, RCap_m1)        edcombo = Chem.EditableMol(combo)        edcombo.AddBond(            LCap_add + second_atom,            RCap_atom + main_mol_noDum.GetNumAtoms(),            order=Chem.rdchem.BondType.SINGLE,        )        main_mol_noDum = edcombo.GetMol()    return Chem.MolToSmiles(main_mol_noDum)def Init_info_Cap(unit_name, smiles_each_ori):    # Get index of dummy atoms and bond type associated with it    try:        dum_index, bond_type = FetchDum(smiles_each_ori)        if len(dum_index) == 1:            dum1 = dum_index[0]        else:            print(                unit_name,                ": There are more or less than one dummy atoms in the SMILES string; ",            )            return unit_name, 0, 0, 0, 0, 'REJECT'    except Exception:        print(            unit_name,            ": Couldn't fetch the position of dummy atoms. Hints: (1) In SMILES string, use '*' for a dummy atom,"            "(2) Check RDKit installation.",        )        return unit_name, 0, 0, 0, 0, 'REJECT'    # Replace '*' with dummy atom    smiles_each = smiles_each_ori.replace(r'*', 'Cl')    # Convert SMILES to XYZ coordinates    convert_smiles2xyz, m1 = smiles_xyz(unit_name, smiles_each, './')    # if fails to get XYZ coordinates; STOP    if convert_smiles2xyz == 'NOT_DONE':        print(            unit_name,            ": Couldn't get XYZ coordinates from SMILES string. Hints: (1) Check SMILES string,"            "(2) Check RDKit installation.",        )        return unit_name, 0, 0, 0, 0, 'REJECT'    # Collect valency and connecting information for each atom    neigh_atoms_info = connec_info('./' + unit_name + '.xyz')    try:        # Find connecting atoms associated with dummy atoms.        # dum1 and dum2 are connected to atom1 and atom2, respectively.        atom1 = neigh_atoms_info['NeiAtom'][dum1].copy()[0]    except Exception:        print(            unit_name,            ": Couldn't get the position of connecting atoms. Hints: (1) XYZ coordinates are not acceptable,"            "(2) Check Open Babel installation.",        )        return unit_name, 0, 0, 0, 0, 'REJECT'    return (        unit_name,        dum1,        atom1,        m1,        neigh_atoms_info,        '',    )# This function try to create a directorydef build_dir(path):    try:        os.mkdir(path)    except OSError:        passdef is_nan(x):    return x != xdef gen_conf_xyz_vasp(unit_name, m1, out_dir, ln, Nconf, NCores_opt, OPLS,):    m2 = Chem.AddHs(m1)    NAttempt = 10000    if NCores_opt != 1:        NAttempt = 1000000    for i in range(10):        cids = AllChem.EmbedMultipleConfs(            m2,            numConfs=Nconf + 10,            numThreads=NCores_opt,            randomSeed=i,            maxAttempts=NAttempt,        )        if len(cids) > 0:            break    n = 0    for cid in cids:        n += 1        AllChem.UFFOptimizeMolecule(m2, confId=cid)        # AllChem.MMFFOptimizeMolecule(m2, confId=cid)        # 使用 os.path.join 来正确地拼接路径和文件名        file_base = '{}_N{}_C{}'.format(unit_name, ln, n)        pdb_filename = os.path.join(out_dir, file_base + '.pdb')        xyz_filename = os.path.join(out_dir, file_base + '.xyz')        Chem.MolToPDBFile(m2, pdb_filename, confId=cid)  # Generate pdb file        Chem.MolToXYZFile(m2, xyz_filename, confId=cid)  # Generate xyz file        # outfile_name = out_dir + unit_name + '_N' + str(ln) + '_C' + str(n)        #        # # if IrrStruc is False:        # Chem.MolToPDBFile(m2, outfile_name + '.pdb', confId=cid)  # Generate pdb file        # Chem.MolToXYZFile(m2, outfile_name + '.xyz', confId=cid)  # Generate xyz file        # Generate OPLS parameter file        if n == 1 and OPLS is True:            print(unit_name, ": Generating OPLS parameter file ...")            if os.path.exists(pdb_filename):                try:                    Converter.convert(                        pdb=pdb_filename,                        resname=unit_name,                        charge=0,                        opt=0,                        outdir= out_dir,                        ln = ln,                    )                    print(unit_name, ": OPLS parameter file generated.")                except BaseException:                    print('problem running LigParGen for {}.pdb.'.format(pdb_filename))        if n == Nconf:            break    return len(cids)