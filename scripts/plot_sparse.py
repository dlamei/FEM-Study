import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("TkAgg")

path = "../build/galerkin_mat.txt"

with open(path) as sparse_mat:
    for line in sparse_mat:
        if '$' in line:
            next(sparse_mat)
            break

    mat = []
    for row in sparse_mat:
        row = np.array(row.strip().split(' '))
        if len(row) <= 1:
            break
        mat.append(row.astype(np.float32))

    matrix = np.array(mat)
    matrix[np.isclose(matrix, 0) == False] = 1

    plt.matshow(matrix)
    plt.show()

