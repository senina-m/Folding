a = [1, 2, 3, 4, 4, 3, 2, 1]


def get_peaks(arr):
    emission = 0
    average = sum(arr) / len(arr) + emission
    e = 1
    result = []
    for i in range(e, len(arr) - e):
        peak = True
        for j in range(i - e, i + e + 1):
            if not (average <= arr[i] and arr[i] >= arr[j]) or (arr[i] == arr[i+1]):
                peak = False
        if peak:
            result.append(i)
    return result


print(get_peaks(a))
