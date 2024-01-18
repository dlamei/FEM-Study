import scipy.sparse as sparse
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("TkAgg")


class SparseMat:

    def __init__(self, data, name = ''):
        self.data = data
        self.name = name

    def spy(self, axis):
        axis.spy(self.data, markersize=0.1, color='black')
        axis.set_title(self.name)
        axis.get_xaxis().set_ticks([])
        axis.get_yaxis().set_ticks([])


def parse_sparse(path):
    with open(path) as sparse_file:
        shape = [int(n) for n in next(sparse_file).split()]

        row = []
        col = []
        val = []
        
        for line in sparse_file:
            (x, y, v) = line.split()
            col.append(int(x))
            row.append(int(y))
            val.append(float(v))

        data = sparse.coo_matrix((val, (row, col)), shape)
        return SparseMat(data)


m = parse_sparse("build/m.txt")
l = parse_sparse("build/l.txt")
u = parse_sparse("build/u.txt")

m.name = "M"
l.name = "L"
u.name = "U"

mats = [m, l, u]

f, axs = plt.subplots(2, 2)

axs = axs.flatten()

for i, m in enumerate(mats):
    m.spy(axs[i])

for i in range(len(axs) - len(mats)):
    f.delaxes(axs[i + len(mats)])


plt.show()
