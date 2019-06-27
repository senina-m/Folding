from Bio.PDB import *
import matplotlib.pyplot as plt
import math as math


class DependentAngle:

    def __init__(self, angle):
        self.angle = angle
        self.dependent_angle = []
        self.has_already_seen = 1


class ResiduesAngles:

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


def get_angle_of_next_level(dict, residue, my_res, level, a):
    a.append(get_angel(i, 0, residue))
    next_level = level + 1
    if a[level] not in dict:
        dict[a[level]] = [1, {}]
        if my_res.num_of_dihedral_angels > level + 1:
            get_angle_of_next_level(dict[tuple(a)[level]][1], residue, my_res, next_level, a)
    else:
        dict[tuple(a)[level]][0] += 1
    return dict


pdbl = PDBList()
# names_file = input()
# f = open(names_file, "r")
# for line in f:
#     filename = f.readline()
filename = '6gov'
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure(filename, filename + '.cif')
models = list(structure.get_models())
for residue in models[0].get_residues():
    for i in range(0, 18):
        if residue.get_resname() == full[i].residue_name:
            if not residue.is_disordered():
                angle = []
                full[i].incidence = get_angle_of_next_level(full[i].incidence, residue, full[i], 0, angle)


def get_peaks(dict):
    emission = 1
    average = sum(list(dict.values())) / (emission * len(dict))
    e = 5
    result = {}
    for key in sorted(dict.keys()):
        peak = True
        for j in range(key - e, key + e + 1):
            if j in dict:
                if not (average <= dict[key] and dict[key] >= dict[j]):
                    peak = False
        if peak:
            result[key] = {}
    return result


def find_result_angles_of_this_level(output_dict, input_dict):
    dict = {}
    for angle in input_dict.keys():
        dict[angle] = input_dict[angle][0]
    # lists = sorted(dict.items())
    # x, y = zip(*lists)
    # plt.plot(x, y)
    # plt.show()
    if len(dict) != 0:
        output_dict = get_peaks(dict)
    else:
        return output_dict
    dict.clear()
    for peak_angle in output_dict.keys():
        find_result_angles_of_this_level(output_dict[peak_angle], input_dict[peak_angle][1])
    return output_dict


for i in range(0, len(full)):
    full[i].result = find_result_angles_of_this_level(full[i].result, full[i].incidence)


def printing_of_angles(num_of_angles, num, result_dict):
    if len(result_dict.keys()) != 0:
        print('Values of  angle:')
        print(result_dict.keys())
        print('Enter value of ', num, ' angle to see the next angle dependeces:')
        angle_val = input()
        if num_of_angles > num and angle_val in result_dict:
            printing_of_angles(num_of_angles, num + 1, result_dict[angle_val])
        else:
            print('There no common values for', num + 1, 'angle :(')
    else:
        print('There no common values for this angle :(')


print('Enter residue name')
input_res_name = input()
for s in full:
    if input_res_name == s.residue_name:
        printing_of_angles(s.num_of_dihedral_angels, 1, s.result)
