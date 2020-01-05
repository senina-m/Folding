import shutil

from Bio.PDB import *
import matplotlib.pyplot as plt
import math as math
import os
import pickle


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
        self.frequency = {}
        self.result = {}


full = [ResiduesAngles("ARG", 5, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "NE"],
                                  ["CG", "CD", "NE", "CZ"], ["CD", "NE", "CZ", "NH1"]]),
        ResiduesAngles("ASN", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]),
        ResiduesAngles("ASP", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]),
        ResiduesAngles("CYS", 1, [["N", "CA", "CB", "SG"]]),
        ResiduesAngles("GLN", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]),
        ResiduesAngles("GLU", 3, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]),
        ResiduesAngles("HIS", 2, [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]]),
        ResiduesAngles("ILE", 2, [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD"]]),
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


# Takes: i: residue # (1..18)
#        level: level
#        residue: structure residue
def get_angel(i, level, residue):
    try:
        vector1 = residue[full[i].arr_of_atoms_names[level][0]].get_vector()
        vector2 = residue[full[i].arr_of_atoms_names[level][1]].get_vector()
        vector3 = residue[full[i].arr_of_atoms_names[level][2]].get_vector()
        vector4 = residue[full[i].arr_of_atoms_names[level][3]].get_vector()
        angle = get_index(calc_dihedral(vector1, vector2, vector3, vector4))
    except KeyError:
        angle = -1
    return angle


# Takes: empty frequency dictionary
#        my_res: my residue class
#        level: level 1..5
#        angle list: list with angles of this residue
#        res_num: 1..18
#        residue: structure residue
# Returns: frequency dictionary angle -> [#, {}]
def get_angle_of_next_level(dict, my_res, level, angle_list, res_num, residue):
    angle_list.append(get_angel(res_num, level, residue))
    if angle_list[level] not in dict:
        dict[angle_list[level]] = [1, {}]
        if my_res.num_of_dihedral_angels > level + 1:
            get_angle_of_next_level(dict[tuple(angle_list)[level]][1], my_res, level + 1, angle_list, res_num, residue)
    else:
        dict[tuple(angle_list)[level]][0] += 1
    return dict


path = "C:/Практика Biocad/PDB/resultbase"
names_file = os.listdir(path)
for filename in names_file:
    if filename[4:8] == ".cif":
        # print(filename)
        shutil.move(path + "/" + filename, "C:/Практика Biocad/Python Projects/" + filename)
        parser = MMCIFParser(QUIET=True)
        # print(filename[0:4] + " " + filename)
        # structure = parser.get_structure(filename[0:3], filename)
        structure = parser.get_structure(filename, filename)
        models = list(structure.get_models())
        for residue in models[0].get_residues():
            for i in range(0, 18):
                if residue.get_resname() == full[i].residue_name:
                    if not residue.is_disordered():
                        full[i].frequency = get_angle_of_next_level(full[i].frequency, full[i], 0, [], i, residue)
        shutil.move("C:/Практика Biocad/Python Projects/" + filename, path + "/" + filename)

# Saving python object to file
file = open("statistics_savings.dat", "w+")
path = "C:/Практика Biocad/Python Projects/statistics_savings.dat"
with open(path, "wb") as f:
    pickle.dump(full, f)
with open(path, "rb") as f:
    print(pickle.load(f))


# Takes: dict:(angle -> frequency)
# Returns: dict:(angle(peak) -> [angles of the peak])
def get_peaks(dict):
    average = sum(list(dict.values())) / len(dict)
    e = 20
    ev = 10
    result = {}
    for key in sorted(dict.keys()):
        if dict[key] > average:
            peak = True
            for j in [c for c in range(key - e, key + e + 1) if c != key]:
                if j in dict:
                    if not (dict[key] > dict[j] or (key < j and dict[key] == dict[j])):
                        peak = False
            if peak:
                result[key] = []
                i = key
                not_bound = True
                while (not_bound and i < 360) and i > 0:
                    not_bound = True
                    i -= 1
                    for k in range(i - ev, i + ev + 1):
                        if k in dict and i in dict:
                            if not dict[i] <= dict[k]:
                                not_bound = False
                for p in range(i, key):
                    if p in dict:
                        result[key].append(p)
                i = key
                not_bound = True
                while (not_bound and i < 360) and i > 0:
                    not_bound = True
                    i += 1
                    for k in range(i - ev, i + ev + 1):
                        if k in dict and i in dict:
                            if not dict[i] <= dict[k]:
                                not_bound = False
                for p in range(key, i):
                    if p in dict:
                        result[key].append(p)
    return result


# Takes:
# output dictionary -  result dict of this level
# input dictionary -  dictionary like: {angle -> [frecuency, {}]}
# Returns dictionary with freaquent angles like: {a -> {a -> {}, a -> {a -> {}}}
def find_result_angles_of_this_level(output_dict, input_dict, level, prev_angle):
    dict = {}
    for angle in input_dict.keys():
        if not angle == -1:
            dict[angle] = input_dict[angle][0]
        dict[angle] = input_dict[angle][0]
    if len(dict) != 0:
        if len(dict) > 3:
            string = 'prev angle value ' + prev_angle + ', this angle # ' + str(level)
            col = ['b', 'g', 'r', 'c', 'm']
            plt.plot(dict.keys(), dict.values(), '.', color=col[level], label=string)
            plt.xlabel('angle')
            plt.ylabel('frequency')
            plt.legend()
            plt.show()
        peaks = get_peaks(dict)
    else:
        return output_dict
    dict.clear()
    for peak_angle in peaks.keys():
        f = make_res_dict(peak_angle, input_dict, peaks[peak_angle])
        output_dict[peak_angle] = f[peak_angle]
        find_result_angles_of_this_level(output_dict[peak_angle][1], output_dict[peak_angle][1], level + 1,
                                         str(peak_angle))
    return output_dict


def sum_dicts(general_dict, dict):
    for angle in dict.keys():
        if angle in general_dict:
            general_dict[angle] = sum_arrs(general_dict[angle], dict[angle])
        else:
            general_dict[angle] = dict[angle]
    return general_dict


def sum_arrs(general_arr, arr):
    res_arr = [general_arr[0] + arr[0]]
    if len(general_arr[1]) > 0 or len(arr[1]) > 0:
        res_arr.append(sum_dicts(general_arr[1], arr[1]))
    else:
        res_arr.append({})
    return res_arr


def make_res_dict(general_angle, dict, angles_list):
    res_array = dict[general_angle]
    for angle in angles_list:
        res_array = sum_arrs(res_array, dict[angle])
    result = {general_angle: res_array}
    return result


plt.figure(figsize=(30, 30))
for i in range(0, len(full)):
    plt.subplot(3, 6, i + 1)
    plt.title(str(full[i].residue_name))
    full[i].result = find_result_angles_of_this_level(full[i].result, full[i].frequency, 1, '')
    plt.legend()
plt.subplots_adjust(top=0.95, bottom=0.001, left=0.10, right=0.95, hspace=0.25, wspace=0.35)
plt.show()


def printing_of_angles(num_of_angles, num, result_dict):
    if len(result_dict.keys()) != 0:
        print('Values of  angle:')
        print(result_dict.keys())
        print('Enter value of ', num, ' angle to see the next angle dependeces:')
        angle_val = input()
        if num_of_angles > num and int(angle_val) in result_dict:
            printing_of_angles(num_of_angles, num + 1, result_dict[int(angle_val)][1])
        else:
            print('There no common values for', num + 1, 'angle :(')
    else:
        print('There no common values for this angle :(')


while True:
    print('Enter residue name')
    input_res_name = input()
    for s in full:
        if input_res_name == s.residue_name:
            printing_of_angles(s.num_of_dihedral_angels, 1, s.result)
