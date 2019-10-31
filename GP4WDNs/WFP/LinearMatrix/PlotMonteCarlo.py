import scipy.io
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
import random

# https://stackoverflow.com/questions/23880138/display-a-3d-bar-graph-using-transparency-and-multiple-colors-in-matplotlib
# https://stackoverflow.com/questions/35210337/can-i-plot-several-histograms-in-3d
# https://matplotlib.org/3.1.0/gallery/color/colormap_reference.html#sphx-glr-gallery-color-colormap-reference-py
# https://stackoverflow.com/questions/45239801/loading-a-cell-array-mat-file-into-python
# https://stackoverflow.com/questions/43831123/python-best-way-to-draw-3d-function-with-random-x-and-y

# https://matplotlib.org/users/usetex.html

# add best fit check here:
#https://matplotlib.org/gallery/statistics/histogram_features.html


# Prepare for the colormap
def get_cmap(N):
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    #scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

def plotDistribution(FigName,Index,SourceData,xlabelName,titleName,elev,azim):
    fig = plt.figure()
    # USE LATEX
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rcParams.update({'font.size': 14})
    ax = fig.add_subplot(111, projection='3d')
    #print(type(ax))

    colorInd = 0
    z = np.arange(Index.size);

    print('start')
    Index = Index[0,:]
    #print(HeadIndex)
    #print(ID_Index)
    for i in Index:
        print(i)
        c=colr_list[colorInd]
        #print(c)
        #ys = np.random.normal(loc=100, scale=10, size=2000)
        ys = SourceData[i-1,:]
        #print(ys)
        hist, bins = np.histogram(ys,density=True)
        xs = (bins[:-1] + bins[1:])/2
        ax.bar(xs, hist, zs=z[colorInd]-0.5, zdir='y', color=c, ec=c, alpha=0.8)
        #print(z)
        colorInd = colorInd + 1

    # Ylabel
    ylab= (ID_Index[i-1][0][0] for i in Index)
    #print(ylab)
    # convert to tuplex
    ylab = tuple(ylab)
    print(ylab)
    plt.yticks(z,ylab)

    ax.set_xlabel(xlabelName)
    ax.set_ylabel('ID')
    ax.set_zlabel('Probility')
    ax.view_init(elev, azim)
    plt.title(titleName)

    #plt.show()
    plt.subplots_adjust(left=0.1, right=0.90, top=0.90, bottom=0.1)
    plt.savefig(FigName,dpi=300)



# Load mat file provided by Matlab Monte_Carlo.m


mat = scipy.io.loadmat('MonteCarloData.mat')
mat1 = scipy.io.loadmat('MonteCarloDemand.mat')
'''
print(type(mat))
print(type(mat['deterministic']))
print(mat['deterministic'])
print(mat['deterministic'].shape)
print(type(mat['FlowIndex']))
print(mat['FlowIndex'])
print(mat['FlowIndex'].size)
print(type(mat['HeadIndex']))
print(mat['HeadIndex'])
print(mat['HeadIndex'].size)
print(type(mat['ID_Index']))
print(mat['ID_Index'])
print(mat['ID_Index'].shape)
print(type(mat['MCSolution']))
print(mat['MCSolution'])
print(mat['MCSolution'].shape)
'''

MCSolution = mat['MCSolution']
HeadIndex= mat['HeadIndex']
FlowIndex = mat['FlowIndex']
ID_Index= mat['ID_Index']
Deterministic= mat['deterministic']

DemandIndex = mat1['DemandIndex']
DemandMC = mat1['demand_MC']




colr_list = []
num_colr = HeadIndex.size + FlowIndex.size
cmap = get_cmap(num_colr*2)

# We create a list with 10 points in the range 1-20 to match our
# cmap_list. Then we put those 10 RGB tuples into colr_list
tmp_list=random.sample(range(1,num_colr*2),num_colr)
for x in tmp_list:
	colr_list.append(cmap(float(x)))


titleName = 'Head Distribution'
FigName = "/Users/shenwang/Dropbox/Probabilistic State Estimation/Figures/HeadDistribution.eps"
Index = HeadIndex
xlabelName = 'Head (ft)'
elev = 10
azim = -132
plotDistribution(FigName,Index,MCSolution,xlabelName,titleName,elev,azim)

titleName = 'Flow Distribution'
FigName = "/Users/shenwang/Dropbox/Probabilistic State Estimation/Figures/FlowDistribution.eps"
Index = FlowIndex
xlabelName = 'Flow (GPM)'
elev = 31
azim = -100
plotDistribution(FigName,Index,MCSolution,xlabelName,titleName,elev,azim)


titleName = 'Demand Distribution'
FigName = "/Users/shenwang/Dropbox/Probabilistic State Estimation/Figures/DemandDistribution.eps"
Index = DemandIndex
xlabelName = 'Demand (GPM)'
elev = 10
azim = -132
plotDistribution(FigName,Index,DemandMC,xlabelName,titleName,elev,azim)
