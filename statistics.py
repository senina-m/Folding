from Bio.PDB import *
from matplotlib.pyplot import ylabel, plot, legend, show, title, text
import math as math

pdbl = PDBList()
filename = '6gov'
pdbl.retrieve_pdb_file(filename, pdir='.', file_format='mmCif')
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure(filename, filename + '.cif')


class ResiduesAngles:

    def __init__(self, residue_name, num_of_dihedral_angels, arr_of_atoms_names):
        """Constructor"""
        n = 360
        self.residue_name = residue_name
        self.num_of_dihedral_angels = num_of_dihedral_angels
        self.arr_of_atoms_names = arr_of_atoms_names
        arr = []
        arr2 = []
        self.arr_result_values_of_angels = []
        for i in range(0, num_of_dihedral_angels):
            arr.append([])
            arr2.append(n * [0])
        # self.num_of_such_res = 0
        self.arr_of_incidence = arr2


full = [ResiduesAngles("ARG", 5, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "NE"],
                                  ["CG", "CD", "NE", "CZ"], ["CD", "NE", "CZ", "NH1"]]),
        ResiduesAngles("ASN", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]),
        ResiduesAngles("ASP", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]),
        ResiduesAngles("CYS", 1, [["N", "CA", "CB", "SG"]]),
        ResiduesAngles("GLN", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]),
        ResiduesAngles("GLU", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]),
        ResiduesAngles("HIS", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]]),
        ResiduesAngles("IlE", 2, [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD"]]),
        ResiduesAngles("LEU", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]),
        ResiduesAngles("LYS", 4, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "CE"],
                                  ["CG", "CD", "CE", "NZ"]]),
        ResiduesAngles("MET", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "SD"], ["CB", "CG", "SD", "CE"]]),
        ResiduesAngles("PHE", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]),
        ResiduesAngles("PRO", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"]]),
        ResiduesAngles("SER", 1, [["N", "CA", "CB", "OG"]]),
        ResiduesAngles("THR", 1, [["N", "CA", "CB", "OG1"]]),
        ResiduesAngles("TRP", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]),
        ResiduesAngles("TYR", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]),
        ResiduesAngles("VAL", 1, [["N", "CA", "CB", "CG1"]])]
# n - width of histogram column
n = 360


def get_index(angle):
    return int(round((angle + math.pi) * (n - 1) * 0.5 / math.pi))


# number_of_residues_in_str = 0
models = list(structure.get_models())
for residue in models[0].get_residues():
    # number_of_residues_in_str = 1 + number_of_residues_in_str
    for i in range(0, 18):
        if residue.get_resname() == full[i].residue_name:
            full[i].num_of_such_res += 1
            for j in range(0, full[i].num_of_dihedral_angels):
                if not residue.is_disordered():
                    vector1 = residue[full[i].arr_of_atoms_names[j][0]].get_vector()
                    vector2 = residue[full[i].arr_of_atoms_names[j][1]].get_vector()
                    vector3 = residue[full[i].arr_of_atoms_names[j][2]].get_vector()
                    vector4 = residue[full[i].arr_of_atoms_names[j][3]].get_vector()
                    angle = calc_dihedral(vector1, vector2, vector3, vector4)
                    full[i].arr_of_incidence[j][get_index(angle)] = 1 + full[i].arr_of_incidence[j][get_index(angle)]
# print(number_of_residues_in_str)


def get_peaks(arr):
    emission = 5
    average = sum(arr) / len(arr) + emission
    e = 3
    result = []
    for i in range(e, len(arr) - e):
        peak = True
        for j in range(i - e, i + e + 1):
            if not (average <= arr[i] and arr[i] >= arr[j]) or (arr[i] == arr[i + 1]):
                peak = False

        if peak:
            result.append(i)

    return result


for i in range(0, len(full)):
    for j in range(0, full[i].num_of_dihedral_angels):
        full[i].arr_result_values_of_angels.append(get_peaks(full[i].arr_of_incidence[j]))
        title(full[i].residue_name)
        plot(full[i].arr_of_incidence[j], label="angel {}".format(j + 1))
        ylabel(str(sum(full[i].arr_of_incidence[j])))
        legend()
    show()

for i in range(0, len(full)):
    print()
    print(full[i].residue_name)
    # print(full[i].num_of_such_res)
    for j in range(0, full[i].num_of_dihedral_angels):
        print(full[i].arr_result_values_of_angels[j])
        # print(len(full[i].arr_of_values[j]))
    print(sum(full[i].arr_of_incidence[j]))
