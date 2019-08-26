import os

# path = 'C:/Практика Biocad/PDB/all/mmCIF'
# path_2_destination_folder = "C:/Практика Biocad/PDB/mmCIF"
# list_of_names = (os.listdir(path))
#
# for name in list_of_names:
#     path_2_current_folder = path + "/" + name
#     for i in os.listdir(path_2_current_folder):
#         os.rename(path_2_current_folder + "/" + i, path_2_destination_folder + "/" + i)

key = 3
e = 1
for j in [c for c in range(key - e, key + e + 1) if c != key]:
    print(j)
