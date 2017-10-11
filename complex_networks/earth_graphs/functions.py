
from numpy import pi, sin, cos, arccos, exp
from numpy import abs, array, arange, shape, reshape, meshgrid
from numpy import random  
from networkx import Graph

                     
def grid(Dtheta, Dphi):
    '''
    Make a regular grid for Earth surface. Dtheta and Dphi are 
    the interval sizes (in degrees) for latitude and longitude, 
    respectively. Latitude varies between -90 and 90 degrees, and 
    longitude goes from 0 to 360 degrees. For each grid piece this 
    function returns values of latitude and longitude (in radians) 
    for the middle point of that grid piece.'''

    theta = (pi/180)*(Dtheta/2 + arange(-90, 90, Dtheta))
    phi = (pi/180)*(Dphi/2 + arange(0, 360, Dphi))
   
    LAT, LON = meshgrid(theta, phi)

    return LAT, LON 


def grid_area(LAT, LON, R_earth=6371e3):
    '''
    Given a regular grid for Earth in radians (LAT, LON), calculates
    the area value in squared meters for each piece of the grid.''' 
    
    Dtheta = (180/pi)*(LAT[0, 1] - LAT[0, 0])
    Dphi = (180/pi)*(LON[1, 0] - LON[0, 0])
    
    Area = (R_earth**2)*cos(LAT)*Dtheta*Dphi

    return Area


def grid_size(grid):
    '''
    '''

    return grid.shape[0]*grid.shape[1]


def coordinates(i, LAT, LON):
    '''
    Given a regular grid for Earth in radians (LAT, LON), this 
    function returns the corresponding coordinates for the given 
    grid point i, where i goes from 1 to the grid size.'''

    size = grid_size(LAT)

    x = reshape(LAT, size)
    y = reshape(LON, size)

    return x[i-1], y[i-1]


def points(i, LAT, LON, R_earth=6371e3):
    '''
    Given a regular grid for Earth in radians (LAT, LON) and a grid
    point number (i), this function returns the coordinate vector
    p = (x, y, z) for the center of that grid piece.'''

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


def probability(i, j, LAT, LON, Lambda=1110e3, R_earth=6371e3):
    '''
    '''
 
    gij = geodesic(i, j, LAT, LON)
    pij = exp(-gij/Lambda)

    return pij


def make_nodes(G, LAT, LON):
    '''
    '''

    size = grid_size(LAT)

    for i in range(size):
        node_number = i + 1
        G.add_node(node_number)

    return G


def make_edges(G, LAT, LON, verbose=True):
    '''
    '''

    # Number of nodes.
    N = grid_size(LAT)

    # Number of edges created.
    m = 0

    # Message to be printed in verbose mode.
    message = 'm = %u\t(i,j)=(%u,%u)\tgij=%.0fkm\tpij=%.6f'

    # Loop over the nodes. 
    for i in range(1, N+1):
        for j in range(i+1, N+1):
        
            gij = geodesic(i, j, LAT, LON)
            pij = probability(i, j, LAT, LON)    
            
            # Create or not a edge.
            if (random.rand() < pij):
                m = m + 1
                G.add_edge(i, j, weight=pij)

                # Print details about the edge in verbose mode.
                if (verbose == True):
                    print(message %(m, i, j, gij/1000, pij))

    return G
        

        























