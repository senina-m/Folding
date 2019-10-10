list = []
file_path = "C:/Школа/Программирование/IdeaProjects/HW/src/ru/senina/school/RSA_simple_numbers.txt"
f = open(file_path, 'r+')
lines = f.readlines()
for line in lines:
    list += line.split()
string = ' '.join(list)
f.truncate(0)
f.write(str(len(list)))
f.write(string)
print(list)
f.close()
