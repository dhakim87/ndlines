from matplotlib import pyplot as plt

species_A = [.7, .3, 0]
species_B = [.3, .7, 0]
species_C = [.2, .2, .6]

Xs = []
Ys = []
for i in range(100):
    for j in range(200):
        Xs.append(species_A[0] * i + species_C[0] * j)
        Ys.append(species_A[1] * i + species_C[1] * j)
plt.scatter(Xs, Ys, c='gray')

Xs = []
Ys = []
for i in range(100):
    Xs.append(species_A[0] * i)
    Ys.append(species_A[1] * i)
    plt.scatter(Xs, Ys, c='r')

Xs = []
Ys = []
for i in range(100):
    Xs.append(species_B[0] * i)
    Ys.append(species_B[1] * i)
    plt.scatter(Xs, Ys, c='g')

Xs = []
Ys = []
for i in range(200):
    Xs.append(species_C[0] * i)
    Ys.append(species_C[1] * i)
    plt.scatter(Xs, Ys, c='b')



plt.xlabel("Species A")
plt.ylabel("Species B")
plt.title("Read Coupling")
plt.show()