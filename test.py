a = {1: {2: {5: {}, 6: {}}, 3: {}}, 4: {}}


def count_leavs(angle_dict, num):
    if len(angle_dict) == 0:
        num = num + 1
        return num
    else:
        for angle in angle_dict.keys():
            num = count_leavs(angle_dict[angle], num)
        return num


print(count_leavs(a, 0))
