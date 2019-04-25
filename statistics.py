from Bio.PDB import *
import numpy as np
import math as math
from scipy.signal import argrelextrema

pdbl = PDBList()
filename = '6gov'
pdbl.retrieve_pdb_file(filename, pdir='.', file_format='mmCif')
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure(filename, filename + '.cif')


class ResiduesAngles:

    def __init__(self, residue_name, num_of_dihedral_angels, arr_of_atoms_names):
        """Constructor"""
        n = 90
        self.residue_name = residue_name
        self.num_of_dihedral_angels = num_of_dihedral_angels
        self.arr_of_atoms_names = arr_of_atoms_names
        arr = []
        arr2 = []
        self.arr_result_values_of_angels = arr
        for i in (0, num_of_dihedral_angels):
            arr.append([])
            arr2.append(n * [0])
        self.arr_of_values = arr
        self.arr_of_incidence = arr2


full = [ResiduesAngles("ARG", 5, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "NE"],
                                  ["CG", "CD", "NE", "CZ"], ["CD", "NE", "CZ", "Nh1"]]),
        ResiduesAngles("ARG", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]),
        ResiduesAngles("ASP", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]),
        ResiduesAngles("CYS", 1, [["N", "CA", "CB", "SG"]]),
        ResiduesAngles("GLN", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]),
        ResiduesAngles("GLU", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]),
        ResiduesAngles("HIS", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]]),
        ResiduesAngles("IlE", 2, [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD"]]),
        ResiduesAngles("LEU", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]),
        ResiduesAngles("LYS", 4, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "CE"],
                                  ["CG", "CD", "CE", "NZ"]]),
        ResiduesAngles("MET", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "SD", "CE"]]),
        ResiduesAngles("PHE", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]),
        ResiduesAngles("PRO", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"]]),
        ResiduesAngles("SER", 1, [["N", "CA", "CB", "CG"]]), ResiduesAngles("THR", 1, [["N", "CA", "CB", "OG1"]]),
        ResiduesAngles("TRP", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]),
        ResiduesAngles("TYR", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]),
        ResiduesAngles("VAL", 1, [["N", "CA", "CB", "CG1"]])]

# n - width of histogram column
n = 1
for residue in structure.get_residues():
    for i in range(0, 18):
        if residue.get_resname() == full[i].residue_name:
            for j in range(full[i].num_of_dihedral_angels):
                if not residue.is_disordered():
                    try:
                        vector1 = residue[full[i].arr_of_atoms_names[j][0]].get_vector()
                        vector2 = residue[full[i].arr_of_atoms_names[j][1]].get_vector()
                        vector3 = residue[full[i].arr_of_atoms_names[j][2]].get_vector()
                        vector4 = residue[full[i].arr_of_atoms_names[j][3]].get_vector()
                        angle = calc_dihedral(vector1, vector2, vector3, vector4)
                        full[i].arr_of_values[j].append(angle)
                        full[i].arr_of_incidence[j][math.floor(angle / n)] += 1
                    except KeyError:
                        pass


for i in (0, len(full)):
    for j in (0, full[i].num_of_dihedral_angels):
        full[i].arr_result_values_of_angels.append(
            argrelextrema(full[i].arr_of_incidence[j], np.greater))
