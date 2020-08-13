import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import os


r"""
Plots ray traing data (midplane)
"""


plt.figure(figsize=(6,6))
ax = plt.gca()

scale = 5
plt.xlim(-scale,scale)
plt.ylim(-scale,scale)


def main(filename,color="black",label='data',linestyle="-"):
    data_type = np.dtype([
        ('block',bool),
        ('x',float),
        ('y',float),
        ('z',float)
        ])

    data = np.fromfile(filename,dtype=data_type)

    indexesToSplit = np.where(data["block"] == True)[0]
    splittedBlocks = []
    previousIndex = 0

    # Separation of different geodesics
    for indexToSplit in indexesToSplit:
        newData = data[previousIndex+1:indexToSplit]
        newData = newData[["x","y","z"]]
        splittedBlocks.append(newData)
        previousIndex = indexToSplit


    nSplittedBlocks = len(splittedBlocks)

    for datablock in splittedBlocks[::8]:
        ax.plot(datablock["x"], datablock["y"],linestyle,alpha=0.5,color=color)

    ax.plot([],[],alpha=0.5,label=label,color=color)



# ==========================================================================
# Configuration of the script
# --------------------------------------------------------------------------
# * (ray-tracing) Data ID (directory)
dataid = "raytr_test1"
# * Root directory where Data ID is located
rootdir = "data/raytr"
# ===========================================================================


for filename in os.listdir(rootdir+"/"+dataid):
    if filename.endswith('.dat'):
        main(rootdir+"/"+dataid+"/"+filename,label=filename[:-4])


ax.legend(ncol=2,loc=4)
plt.show()
