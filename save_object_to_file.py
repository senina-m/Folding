import pickle

full = {1: 2, 2: 3}
file = open("statistics_savings.dat", "w+")
path = "statistics_savings.dat"

with open(path, "wb") as f:
    pickle.dump(full, f)
with open(path, "rb") as f:
    print(pickle.load(f))
