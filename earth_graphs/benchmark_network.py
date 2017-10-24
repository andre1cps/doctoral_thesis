"""
Benchmark network
=================

Here is an attempt to reproduce the benchmark complex network 
described in the reference below. Earth is represented as a regular
grid, where the middle point of each piece in the grid is a node
of the complex network, and the edges are constructed 
probabilistically as a function of the geodesic distance among the
nodes. Earth radius is considered to be 6371 km. These are the 
functions defined here:

- probability:
    The probability that two nodes are linked is a function
    of the geodesic distance between them.
 
- make_network:
    Create the benchmark network and returns the corresponding graph.

Reference:
Zemp, D., Wiedermann, M., Kurths, J., Rammig, A., and Donges, J. F. 
(2014). Node-weighted measures for complex networks with directed and
weighted edges for studying continental moisture recycling. EPL 
(Europhysics Letters), 107(5):58005.

"""

import numpy as np
import networkx as nx
from .grids import geodesic, grid_size


def probability(i, j, LAT, LON, Lambda=1110e3, R_earth=6371e3):
    """
    The probability pij that node i and node j are linked is a function
    of the geodesic distance gij between them:

    pij = exp(-gij/Lambda)

    with Lambda representing the typical length scale. In this case 
    it was chosen Lambda = 1110 km to ensure an edge density of about 
    0.02.

    """
   
    gij = geodesic(i, j, LAT, LON)
    
    return np.exp(-gij/Lambda)


def make_network(LAT, LON, verbose=False):
    """
    Create the benchmark network according to the Reference and 
    returns the corresponging graph. In the verbose mode it prints
    out information about the edges that were created.

    Example:

    >>> LAT, LON = grid(10, 10)
    >>> G = make_network(LAT, LON)

    Reference:
    Zemp, D., Wiedermann, M., Kurths, J., Rammig, A., and Donges, J. F. 
    (2014). Node-weighted measures for complex networks with directed 
    and weighted edges for studying continental moisture recycling.
    EPL (Europhysics Letters), 107(5):58005.

    """

    # Create the graph.
    G = nx.Graph()

    # Number of nodes.
    N = grid_size(LAT)

    # Create all the nodes.
    for i in range(N):
        node_number = i + 1
        G.add_node(node_number)

    # Number of created edges.
    m = 0

    # Format the message to be printed in verbose mode.
    message = 'edge %u;\t(i,j)=(%u,%u);\tgij=%.0fkm;\tpij=%.6f'

    # Loop over the nodes. 
    for i in range(1, N+1):
        for j in range(i+1, N+1):
        
            gij = geodesic(i, j, LAT, LON)
            pij = probability(i, j, LAT, LON)    
            
            # Create or not a edge.
            if (np.random.rand() < pij):
                m = m + 1
                G.add_edge(i, j, weight=pij)

                # Print details about the edge in verbose mode.
                if (verbose == True):
                    print(message %(m, i, j, gij/1000, pij))

    return G

