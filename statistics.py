from Bio.PDB import *

pdbl = PDBList()
filename = '6gov'
pdbl.retrieve_pdb_file(filename, pdir='.', file_format='mmCif')
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure(filename, filename + '.cif')


class ResiduesAngles:

    def __init__(self, residue_name, num_of_dihedral_angels, arr_of_atoms_names, ):
        """Constructor"""
        self.residue_name = residue_name
        self.num_of_dihedral_angels = num_of_dihedral_angels
        self.arr_of_atoms_names = arr_of_atoms_names
        arr = []
        for i in range(0, num_of_dihedral_angels):
            arr.append([])
        self.arr_of_values = arr


full = []
full.append(ResiduesAngles("ARG", 5, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "NE"],
                                      ["CG", "CD", "NE", "CZ"], ["CD", "NE", "CZ", "Nh1"]]))
full.append(ResiduesAngles("ARG", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]))
full.append(ResiduesAngles("ASP", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]))
full.append(ResiduesAngles("CYS", 1, [["N", "CA", "CB", "SG"]]))
full.append(ResiduesAngles("GLN", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]))
full.append(ResiduesAngles("GLU", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]))
full.append(ResiduesAngles("HIS", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]]))
full.append(ResiduesAngles("IlE", 2, [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD"]]))
full.append(ResiduesAngles("LEU", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]))
full.append(ResiduesAngles("LYS", 4, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "CE"],
                                      ["CG", "CD", "CE", "NZ"]]))
full.append(ResiduesAngles("MET", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "SD", "CE"]]))
full.append(ResiduesAngles("PHE", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]))
full.append(ResiduesAngles("PRO", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"]]))
full.append(ResiduesAngles("SER", 1, [["N", "CA", "CB", "CG"]]))
full.append(ResiduesAngles("THR", 1, [["N", "CA", "CB", "OG1"]]))
full.append(ResiduesAngles("TRP", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]))
full.append(ResiduesAngles("TYR", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]))
full.append(ResiduesAngles("VAL", 1, [["N", "CA", "CB", "CG1"]]))

for residue in structure.get_residues():
    for i in range(0, 18):
        if residue.get_resname() == full[i].residue_name:
            for j in range(full[i].num_of_dihedral_angels):
                if not residue.is_disordered():
                # try:
                    vector1 = residue[full[i].arr_of_atoms_names[j][0]].get_vector()
                    vector2 = residue[full[i].arr_of_atoms_names[j][1]].get_vector()
                    vector3 = residue[full[i].arr_of_atoms_names[j][2]].get_vector()
                    vector4 = residue[full[i].arr_of_atoms_names[j][3]].get_vector()
                    angle = calc_dihedral(vector1, vector2, vector3, vector4)
                    full[i].arr_of_values[j].append(angle)
                # except KeyError:
                #     print(residue.get_resname())
print(full[11].arr_of_values)
