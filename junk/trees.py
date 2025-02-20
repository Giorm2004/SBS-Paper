import numpy as np
import tskit 
import matplotlib.pyplot as plt

ts = tskit.load("./overlay_II.trees")
ts = ts.simplify()

height_for_pos = np.zeros(int(ts.sequence_length))

for tree in ts.trees():
    mean_height = np.mean([tree.time(root) for root in tree.roots])
    left, right = map(int, tree.interval)
    height_for_pos[left: right] = mean_height

print(height_for_pos.shape)
np.savetxt("heights.gz", height_for_pos)
plt.plot(height_for_pos)
plt.savefig('trees.png')
