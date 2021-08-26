
"""
Module to create topo and qinit data files for this example.
"""

from __future__ import absolute_import
from clawpack.geoclaw.topotools import Topography
from numpy import *

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 301
    nypoints =301
    xlower = 0.e0
    xupper = 1.e0
    ylower = 0.e0
    yupper = 1.e0
    outfile= "Lbar.tt3"     

    topography = Topography(topo_func=topo)
    topography.x = linspace(xlower,xupper,nxpoints)
    topography.y = linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=3, Z_format="%22.15e")

def makeqinit():
    """
    Create qinit data file
    """
    nxpoints = 301
    nypoints = 301
    xlower = 0.e0
    xupper = 1.e0
    yupper = 0.e0
    ylower = 1.e0
    outfile= "qinit.xyz"     

    topography = Topography(topo_func=qinit)
    topography.x = linspace(xlower,xupper,nxpoints)
    topography.y = linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=1)

def topo(x,y):
    """
    Island
    """
    ze = -((x-40.)**2 + (y-35.)**2)/20.
    
    #z_island = where(ze>-10., 100.*exp(ze), 0.)
    z_island = where(ze>-10., 150.*exp(ze), 0.)
    z = -50 + z_island
    return z


def qinit(x,y):
    """
    Dam break
    """
    from numpy import where
    eta = where(x<10, 40., 0.)
    return eta

if __name__=='__main__':
    maketopo()
    makeqinit()
