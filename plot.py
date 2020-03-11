#rocket-nozzle-tools
__author__ = "Ben Appleby"
__email__ = "ben.appleby@sky.com","b6040585@my.shu.ac.uk"
__copyright__ = "Copyright 2019, Ben Appleby"

from numpy import array
from matplotlib import pyplot as plt

def plotGraph(coordinates):
    data = array(coordinates)
    x,y = data.T
    plt.scatter(x,y)
    plt.axis('equal')
    plt.ylim(bottom=0)
    plt.show()