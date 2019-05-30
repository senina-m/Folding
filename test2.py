import matplotlib.pyplot as plt

D = {u'Label1': 26, u'Label2': 17, u'Label3':30}

plt.bar(range(len(D)), list(D.values()))

plt.show()

# plot(full[i].arr_of_incidence[j], label="angel {}".format(j + 1))
# plt.ylabel(str(sum(full[i].arr_of_incidence[j])))
# plt.legend()
