"""
Grids
=====

Tools for creating and dealing with regular grids for Earth surface.
Earth radius is considered to be 6371 km. These are the functions 
defined here:

- grid: 
    Make a regular grid.
    
- grid_size: 
    Number of pieces in the grid.

- grid_area: 
    Area value for each piece of the grid. 

- coordinates_latlon:
    For a given a grid number, finds the geographic coordinates 
    for the center of the corresponding grid piece. 

- coordinates_xyz:
    For a given grid number, returns the rectangular coordinates 
    for the center of the corresponding grid piece.

- geodesic:
    Calculates the geodesic distance between the two grid points 
    localizations on the Earth.    

"""

import numpy as np


def grid(Dtheta, Dphi):
    """
    Create a meshgrid. The parameters Dtheta and Dphi are the size 
    intervals in degrees for latitude and longitude, respectively.
    Latitude goes from -90 to 90 degrees and longitude does from 
    0 to 360 degrees. Latitude and longitude for each grid point 
    correspond to the center of that grid piece, both given in 
    radians.

    Example:

    >>> LAT, LON = grid(5, 5)
    >>> LAT.shape
    (72, 36)
    >>> LON.shape
    (72, 36)

    """

    theta = (np.pi/180)*(Dtheta/2 + np.arange(-90, 90, Dtheta))
    phi = (np.pi/180)*(Dphi/2 + np.arange(0, 360, Dphi))
   
    LAT, LON = np.meshgrid(theta, phi)

    return LAT, LON 


def grid_size(grid):
    """
    For a given meshgrid, calculates the number of grid pieces.

    Example:

    >>> LAT, LON = grid(5, 5)
    >>> SIZE = grid_size(LAT)
    >>> SIZE
    2592
  
    """

    return grid.shape[0]*grid.shape[1]


def grid_area(LAT, R_earth=6371e3):
    """
    For given latitude meshgrid, calculates the area in squared meters 
    for each grid piece:

    Area = (R_earth**2)*cos(LAT)*Dtheta*Dphi

    Note that area doesn't depend on the longitude degree.

    Example:

    >>> LAT, LON = grid(5, 5)
    >>> AREA = grid_area(LAT, LON)
    >>> AREA.shape
    (72, 36)

    """
    
    Dtheta = (180/np.pi)*(LAT[0, 1] - LAT[0, 0])
    Dphi = (180/np.pi)*(LON[1, 0] - LON[0, 0])
    
    return (R_earth**2)*np.cos(LAT)*Dtheta*Dphi


def coordinates_latlon(i, LAT, LON):
    """
    For a given grid number i, calculates geographic coordinates
    (lat_i, lon_i) in radians for the center of that grid piece.
    Grid number i goes from 1 to the grid size.

    Example:

    >>> LAT, LON = grid(5, 5)
    >>> lat_i, lon_i = coordinates_latlon(1, LAT, LON)
    >>> lat_i
    -1.5271630954950384
    >>> lon_i
    0.043633231299858237

    """

    size = grid_size(LAT)

    lat_i = np.reshape(LAT, size)[i-1]
    lon_i = np.reshape(LON, size)[i-1]

    return np.array([lat_i, lon_i])


def coordinates_xyz(i, LAT, LON, R_earth=6371e3):
    """
    For a given grid number i, calculates rectangular coordinates
    (x, y, z) in meters for the center of that grid piece:

    x = R_earth*cos(lat_i)*cos(lon_i)
    y = R_earth*cos(lat_i)*sin(lon_i)
    z = R_earth*sin(lon_i)

    Grid number i goes from 1 to the grid size.

    Example:

    >>> LAT, LON = grid(5, 5)
    >>> x, y, z = coordinates_xyz(1, LAT, LON)
    >>> x
    277634.61852266517
    >>> y
    12121.789228744608
    >>> z
    -6364936.2196980156

    """

    lat_i, lon_i = coordinates_latlon(i, LAT, LON)
    
    x = R_earth*np.cos(lat_i)*np.cos(lon_i)
    y = R_earth*np.cos(lat_i)*np.sin(lon_i)
    z = R_earth*np.sin(lat_i)
    
    return np.array([x, y, z])


def geodesic(i, j, LAT, LON, R_earth=6371e3):
    """ 
    Calculate the geodesic distance in meters between the i and j grid
    points on the Earth grid.
    
    Example:

    >>> LAT, LON = grid(5, 5)
    >>> d12 = geodesic(1, 2, LAT, LON)
    >>> d12
    555974.63322280894
    
    """

    vi = coordinates_xyz(i, LAT, LON)
    vj = coordinates_xyz(j, LAT, LON)

    alpha = (vi@vj)/R_earth**2

    # This is necessary because of numerical problems in calculating
    # the dot product.
    if (alpha > 1):
        alpha = alpha - 1e-15

    if (alpha < -1):
        alpha = alpha + 1e-15

    return R_earth*np.arccos(alpha)
