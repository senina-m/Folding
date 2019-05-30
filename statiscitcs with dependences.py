from Bio.PDB import *
import matplotlib.pyplot as plt
import math as math

pdbl = PDBList()
filename = '6gov'
pdbl.retrieve_pdb_file(filename, pdir='.', file_format='mmCif')
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure(filename, filename + '.cif')


class DependentAngle:

    def __init__(self, angle):
        self.angle = angle
        self.dependent_angle = []
        self.has_already_seen = 1


class ResiduesAngles:

    # incidence: Dict[int, list]

    def __init__(self, residue_name, num_of_dihedral_angels, arr_of_atoms_names):
        n = 360
        self.residue_name = residue_name
        self.num_of_dihedral_angels = num_of_dihedral_angels
        self.arr_of_atoms_names = arr_of_atoms_names
        self.incidence = {}
        self.result = {}


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


def get_angel(i, j, residue):
    vector1 = residue[full[i].arr_of_atoms_names[j][0]].get_vector()
    vector2 = residue[full[i].arr_of_atoms_names[j][1]].get_vector()
    vector3 = residue[full[i].arr_of_atoms_names[j][2]].get_vector()
    vector4 = residue[full[i].arr_of_atoms_names[j][3]].get_vector()
    angle = get_index(calc_dihedral(vector1, vector2, vector3, vector4))
    return angle


models = list(structure.get_models())
for residue in models[0].get_residues():
    for i in range(0, 18):
        if residue.get_resname() == full[i].residue_name:
            if not residue.is_disordered():
                angle1 = get_angel(i, 0, residue)
                if angle1 not in full[i].incidence:
                    full[i].incidence[angle1] = [1, {}]
                    if full[i].num_of_dihedral_angels > 1:
                        angle2 = get_angel(i, 1, residue)
                        if angle2 not in full[i].incidence[angle1][1]:
                            full[i].incidence[angle1][1][angle2] = [1, {}]
                            if full[i].num_of_dihedral_angels > 2:
                                angle3 = get_angel(i, 2, residue)
                                if angle3 not in full[i].incidence[angle1][1][angle2][1]:
                                    full[i].incidence[angle1][1][angle2][1][angle3] = [1, {}]
                                    if full[i].num_of_dihedral_angels > 3:
                                        angle4 = get_angel(i, 3, residue)
                                        if angle4 not in full[i].incidence[angle1][1][angle2][1][angle3][1]:
                                            full[i].incidence[angle1][1][angle2][1][angle3][1][angle4] = [1, {}]
                                            if full[i].num_of_dihedral_angels == 5:
                                                angle5 = get_angel(i, 4, residue)
                                                if angle5 not in \
                                                        full[i].incidence[angle1][1][angle2][1][angle3][1][angle4][1]:
                                                    full[i].incidence[angle1][1][angle2][1][angle3][1][angle4][1][
                                                        angle5] = [1]
                                                else:
                                                    full[i].incidence[angle1][1][angle2][1][angle3][1][angle4][1][
                                                        angle5][0] += 1
                                        else:
                                            full[i].incidence[angle1][1][angle2][1][angle3][1][angle4][0] += 1
                                else:
                                    full[i].incidence[angle1][1][angle2][1][angle3][0] += 1
                        else:
                            full[i].incidence[angle1][1][angle2][0] += 1
                else:
                    full[i].incidence[angle1][0] += 1


def get_peaks(dict):
    emission = 5
    average = 0
    for key in dict.keys():
        average += dict[key][0]
    average = emission + average/len(dict)
    e = 5
    result = {}
    for key in sorted(dict.keys()):
        peak = True
        for j in range(key - e, key + e + 1):
            if j in dict:
                if not (average <= dict[key] and dict[key] >= dict[j]):
                    peak = False
        if peak:
            result = {key: {}}
    return result


dictionary = {}
for i in range(0, len(full)):
    for key in full[i].incidence.keys():
        dictionary = {key: full[i].incidence[key][0]}
        full[i].result = get_peaks(dictionary)
    for key1 in full[i].result.keys():
        dictionary.clear()
        for key2 in full[i].incidence[key1][1].keys():
            dictionary = {key1: full[i].incidence[key1][1][key2][0]}
            full[i].result[key1] = get_peaks(dictionary)
            dictionary.clear()
            for key3 in full[i].incidence[key1][1][key2][1].keys():
                dictionary = {key2: full[i].incidence[key1][1][key2][1][key3][0]}
                full[i].result[key1][key2] = get_peaks(dictionary)
                dictionary.clear()
                for key4 in full[i].incidence[key1][1][key2][1][key3][1].keys():
                    dictionary = {key3: full[i].incidence[key1][1][key2][1][key3][1][key4][0]}
                    full[i].result[key1][key2][key3] = get_peaks(dictionary)
                    dictionary.clear()
                    for key5 in full[i].incidence[key1][1][key2][1][key3][1][key4][1].keys():
                        dictionary = {key4: full[i].incidence[key1][1][key2][1][key3][1][key4][1][key5][0]}
                        full[i].result[key1][key2][key3][key4] = get_peaks(dictionary)

for i in range(0, len(full)):
    pass
# plt.title(full[i].residue_name)
# plt.bar(range(len(full[i].incidence)), list(full[i].incidence.values()[0]))
# plt.xticks(range(len(full[i].incidence)), list(full[i].incidence.keys()))
# plt.show()
