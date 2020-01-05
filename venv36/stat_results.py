# import matplotlib.pyplot as plt
import pickle
import numpy


class ResiduesAngles:

    def __init__(self, residue_name, num_of_dihedral_angels, arr_of_atoms_names):
        n = 360
        self.residue_name = residue_name
        self.num_of_dihedral_angels = num_of_dihedral_angels
        self.arr_of_atoms_names = arr_of_atoms_names
        self.frequency = {}
        self.result = {}


# Saving python object to file
file = open("statistics_savings.dat", "w+")
path = "C:/Практика Biocad/Python Projects/statistics_savings.dat"
with open(path, "rb") as f:
    full = pickle.load(f)

for i in full:
    i.result.clear()


# Takes: dict:(angle -> frequency)
# Returns: dict:(angle(peak) -> [angles of the peak])
def get_peaks(dict):
    average = sum(list(dict.values())) / len(dict)
    edge = 10
    e = 20
    ev = 10
    b = ev / 2
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
                not_bound_left = True
                while (not_bound_left and i <= 360) and i > 0:
                    not_bound_left = True
                    i -= 1
                    count = 0
                    for k in range(0, ev):
                        if (i - k) in dict and dict[i] < dict[i - k]:
                            count += 1
                    if count >= b:
                        not_bound_left = False
                for p in range(i, key):
                    if p in dict:
                        result[key].append(p)
                i = key
                not_bound_right = True
                while (not_bound_right and i < 360) and i >= 0:
                    not_bound_right = True
                    i += 1
                    count = 0
                    for k in range(0, ev):
                        if (i + k) in dict and dict[i] < dict[i + k]:
                            count += 1
                    if count >= b:
                        not_bound_right = False
                for p in range(key + 1, i + 1):
                    if p in dict:
                        result[key].append(p)
    # delete 0 or 360 angle if both are stable
    list = result.keys()
    if list[0] <= edge and list[len(list) - 1] >= 360 - edge:
        if dict[list[0]] > dict[list[len(list) - 1]]:
            del result[list[len(list) - 1]]
        else:
            del result[list[0]]
    return result


# Takes:
# output dictionary -  result dict of this level
# input dictionary -  dictionary like: {angle -> [frecuency, {}]}
# Returns dictionary with freaquent angles like: {a -> {a -> {}, a -> {a -> {}}}
def find_result_angles_of_this_level(output_dict, input_dict, level, prev_angle, num_of_levels):
    dict = {}
    for angle in range(0, 361):
        if angle in input_dict.keys():
            dict[angle] = input_dict[angle][0]
        else:
            dict[angle] = 0
    if len(dict) != 0:
        if level <= num_of_levels:
            # string = "First dihedral angle"
            d = {}
            for j in dict:
                if dict[j] != 0:
                    d[j] = dict[j]
            # if len(dict) > 3:
            #     string = '' + prev_angle + ',#   ' + str(level)
            #     col = ['b', 'g', 'r', 'c', 'm']
            #     lists = sorted(d.items())
            #
            #     x, y = zip(*lists)
            #     print(prev_angle + ', # ' + str(level))
            #     plt.plot(x, y, '.', color=numpy.random.rand(3,), label=string)
            #
            #     # Printing of the whole results of the first angle
            #     # plt.plot(x, y, '.', color=col[level - 1])
            #     plt.xlabel('angle')
            #     plt.ylabel('frequency')
            #     plt.legend()
            #     # if level == num_of_levels:
            #     #     plt.show()
        peaks = get_peaks(dict)
    else:
        return output_dict
    dict.clear()
    for peak_angle in peaks.keys():
        f = make_res_dict(peak_angle, input_dict, peaks[peak_angle])
        output_dict[peak_angle] = find_result_angles_of_this_level({}, f[peak_angle][1], level + 1,
                                                                   prev_angle + "-" + str(peak_angle), num_of_levels)
    return output_dict


def sum_dicts(dict1, dict2):
    res_dict = dict(dict1)
    for angle in dict2.keys():
        if angle in res_dict:
            res_dict[angle] = sum_arrs(res_dict[angle], dict2[angle])
        else:
            res_dict[angle] = dict2[angle]
    return res_dict


def sum_arrs(arr1, arr2):
    res_arr = [arr1[0] + arr2[0]]
    if len(arr1[1]) > 0 or len(arr2[1]) > 0:
        res_arr.append(sum_dicts(arr1[1], arr2[1]))
    else:
        res_arr.append({})
    return res_arr


def make_res_dict(peak_angle, input_dict, neighbor_angles_list):
    res_array = input_dict[peak_angle]
    for angle in neighbor_angles_list:
        if angle in input_dict:
            res_array = sum_arrs(res_array, input_dict[angle])
    return {peak_angle: res_array}


for i in range(0, len(full)):
    if i != 7:
        # plt.title(str(full[i].residue_name))
        print(str(full[i].residue_name))
        full[i].result = find_result_angles_of_this_level(full[i].result, full[i].frequency, 1, 'no',
                                                          full[i].num_of_dihedral_angels)
        # plt.legend()
        # plt.show()


# Printing of the whole results of the first angle for 17 residues
# plt.figure(figsize=(30, 30))
# for i in range(0, len(full) + 1):
#     if i < 7:
#         plt.subplot(3, 6, i + 1)
#         plt.title(str(full[i].residue_name))
#         full[i].result = find_result_angles_of_this_level(full[i].result, full[i].frequency, 1, '')
#         plt.legend()
#     elif i > 8:
#         plt.subplot(3, 6, i - 1)
#         plt.title(str(full[i - 1].residue_name))
#         full[i-1].result = find_result_angles_of_this_level(full[i-1].result, full[i-1].frequency, 1, '')
#         plt.legend()
# plt.subplots_adjust(top=0.9, bottom=0.1, left=0.10, right=0.95, hspace=0.9, wspace=0.9)
# plt.show()


def count_leavs(angle_dict, num):
    if len(angle_dict) == 0:
        num = num + 1
        return num
    else:
        for angle in angle_dict.keys():
            num = count_leavs(angle_dict[angle], num)
        return num


for res in full:
    print(res.residue_name + ": " + str(count_leavs(res.result, 0)))


def printing_of_angles(num_of_angles, num, result_dict):
    if len(result_dict.keys()) != 0:
        print('Values of  angle:')
        print(result_dict.keys())
        print('Enter value of ', num, ' angle to see the next angle dependeces:')
        angle_val = input()
        if num_of_angles > num and (int(angle_val) in result_dict):
            printing_of_angles(num_of_angles, num + 1, result_dict[int(angle_val)])
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
