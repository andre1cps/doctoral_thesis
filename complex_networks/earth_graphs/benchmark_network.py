"""
Benchmark networks 
==================

"""


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


