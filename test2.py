import matplotlib.pylab as plt


def get_peaks(dict):
    emission = 5
    average = emission + (sum(list(dict.values())) / (3 * len(dict)))
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


d = {244: 4, 356: 6, 9: 5, 3: 8, 117: 5, 359: 7, 10: 7, 112: 16, 0: 8, 115: 12, 243: 3, 114: 7, 116: 8, 113: 10, 107: 8, 5: 12, 250: 3, 105: 5, 1: 6, 110: 15, 109: 7, 354: 3, 2: 10, 118: 4, 108: 4, 357: 1, 248: 3, 99: 1, 12: 4, 119: 7, 111: 7, 4: 7, 106: 5, 358: 3, 18: 3, 11: 4, 104: 4, 15: 1, 16: 1, 13: 5, 247: 3, 6: 8, 351: 1, 253: 1, 353: 2, 7: 5, 14: 1, 249: 2, 123: 3, 17: 2, 122: 3, 121: 2, 258: 1, 254: 2, 237: 2, 240: 2, 8: 3, 246: 2, 232: 2, 355: 2, 245: 3, 252: 1, 103: 4, 88: 1, 124: 1, 120: 1, 349: 1, 102: 1, 242: 3, 101: 1, 352: 1}
print(sorted(d.keys()))
a = get_peaks(d)
print(a)
lists = sorted(d.items())
x, y = zip(*lists)
plt.plot(x, y)
plt.show()
