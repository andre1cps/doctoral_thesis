"""
Grids
=====

Tools for creating and dealing with regular grids for Earth surface:

- grid: 
    Make a regular grid.
    
- grid_area: 
    Area value for each piece of the grid. 

- grid_size: 
    Number of pieces in the grid.

- coordinates:
    For a given a grid number, finds the coordinates (lat, lon) for 
    the center of the corresponging grid piece. 

- position:
    For a given grid number, returns the rectangular coordinates 
    (x, y, z) for the center of the corresponging grid piece.

- geodesic:
    Calculates the geodesic distance between the two grid points on 
    the Earth grid.    

"""

from numpy import pi, sin, cos, arccos, exp
from numpy import abs, array, arange, shape, reshape, meshgrid
from numpy import random  
from networkx import Graph

                     
def grid(Dtheta, Dphi):
    """

    """

    theta = (pi/180)*(Dtheta/2 + arange(-90, 90, Dtheta))
    phi = (pi/180)*(Dphi/2 + arange(0, 360, Dphi))
   
    LAT, LON = meshgrid(theta, phi)

    return LAT, LON 


def grid_area(LAT, LON, R_earth=6371e3):
    """

    """
    
    Dtheta = (180/pi)*(LAT[0, 1] - LAT[0, 0])
    Dphi = (180/pi)*(LON[1, 0] - LON[0, 0])
    
    Area = (R_earth**2)*cos(LAT)*Dtheta*Dphi

    return Area


def grid_size(grid):
    """

    """

    return grid.shape[0]*grid.shape[1]





def coordinates(i, LAT, LON):
    """

    """

    size = grid_size(LAT)

    x = reshape(LAT, size)
    y = reshape(LON, size)

    return x[i-1], y[i-1]


def position(i, LAT, LON, R_earth=6371e3):
    """

    """

    theta, phi = coordinates(i, LAT, LON)
    
    p = R_earth*array([cos(theta)*cos(phi), +
                       cos(theta)*sin(phi), +
                       sin(theta)])

    return p


def correct(alpha):
    '''
    This is necessary because of numerical problems
    in calculating the dot product.'''
    
    eps = 1e-15
 
    if (alpha > 1):
        alpha = alpha - eps

    if (alpha < -1):
        alpha = alpha + eps

    return alpha


def geodesic(i, j, LAT, LON, R_earth=6371e3):
    ''' 
    Calculate the geodesic distance between the i and j grid points
    on the Earth grid.'''

    pi = points(i, LAT, LON, R_earth)
    pj = points(j, LAT, LON, R_earth)

    alpha = (pi @ pj)/R_earth**2

    if (abs(alpha) > 1):
        alpha = correct(alpha)

    return R_earth*arccos(alpha)
