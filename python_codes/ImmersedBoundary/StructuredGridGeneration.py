from __future__ import division
import numpy as np
from numbers import Number
import copy
import matplotlib.pyplot as plt
from Grids import Point, Boundary, TwoDGrid, BdType

# Univariate method (or something else, dont know)
def get_univariate_grid(bd1, bd2, bd3, bd4, stretching_f=lambda xi: xi):
    if bd1.bdtype == BdType.xi_min or bd1.bdtype == BdType.xi_max:
        jmx = bd1.n + 1
    elif bd1.bdtype == BdType.eta_min or bd1.bdtype == BdType.eta_max:
        imx = bd1.n + 1

    if bd2.bdtype == BdType.xi_min or bd2.bdtype == BdType.xi_max:
        jmx = bd2.n + 1
    elif bd2.bdtype == BdType.eta_min or bd2.bdtype == BdType.eta_max:
        imx = bd2.n + 1

    if bd3.bdtype == BdType.xi_min or bd3.bdtype == BdType.xi_max:
        jmx = bd3.n + 1
    elif bd3.bdtype == BdType.eta_min or bd3.bdtype == BdType.eta_max:
        imx = bd3.n + 1

    if bd4.bdtype == BdType.xi_min or bd4.bdtype == BdType.xi_max:
        jmx = bd4.n + 1
    elif bd4.bdtype == BdType.eta_min or bd4.bdtype == BdType.eta_max:
        imx = bd4.n + 1

    AirfoilGrid = TwoDGrid(imx, jmx)
    Grid = AirfoilGrid.Grid
    AirfoilGrid.set_boundary(bd1)
    AirfoilGrid.set_boundary(bd2)
    AirfoilGrid.set_boundary(bd3)
    AirfoilGrid.set_boundary(bd4)
    
    
    # In the given case, assuming constant xi lines are straight lines
    for xi in range(0, imx):
        # start and end are instances of class Point
        start = Grid[xi, 0]
        end = Grid[xi, jmx-1]
        for eta in range(1, jmx-1):
            s = stretching_f(eta / (jmx - 1))
            Grid[xi, eta] = ((1 - s) * start) + (s * end)

    return AirfoilGrid
            
        
# Transfinite interpolation
def get_transfinite_interpolated_grid(bd1, bd2, bd3, bd4):
    if bd1.bdtype == BdType.xi_min or bd1.bdtype == BdType.xi_max:
        jmx = bd1.n + 1
    elif bd1.bdtype == BdType.eta_min or bd1.bdtype == BdType.eta_max:
        imx = bd1.n + 1

    if bd2.bdtype == BdType.xi_min or bd2.bdtype == BdType.xi_max:
        jmx = bd2.n + 1
    elif bd2.bdtype == BdType.eta_min or bd2.bdtype == BdType.eta_max:
        imx = bd2.n + 1

    AirfoilGrid = TwoDGrid(imx, jmx)
    Grid = AirfoilGrid.Grid
    AirfoilGrid.set_boundary(bd1)
    AirfoilGrid.set_boundary(bd2)
    AirfoilGrid.set_boundary(bd3)
    AirfoilGrid.set_boundary(bd4)

    # Notations as per class notes
    ra = Grid[0,0]
    rb = Grid[imx-1, 0]
    rc = Grid[0, jmx-1]
    rd = Grid[imx-1, jmx-1]
    for xi in range(imx):
        for eta in range(jmx):
            r_i = ((1 - eta/(jmx-1)) * Grid[xi, 0]) + \
                  ((eta/(jmx-1)) * Grid[xi, jmx-1])
            r_j = ((1 - xi/(imx-1)) * Grid[0, eta]) + \
                  ((xi/(imx-1)) * Grid[imx-1, eta])
            Grid[xi, eta] = r_i + r_j - \
                                    ((1 - xi/(imx-1))*(1 - eta/(jmx-1))*ra) - \
                                        ((xi/(imx-1))*(1 - eta/(jmx-1))*rb) - \
                                        ((1 - xi/(imx-1))*(eta/(jmx-1))*rc) - \
                                        ((xi/(imx-1))*(eta/(jmx-1))*rd)
            
    return AirfoilGrid    

# Implement Elliptic grid generation
def elliptic_grid_generation(Grid_obj, eps=0.001, iter_max=80000, restartVTK = \
                                                    None, saveInterval = -1):
    # Note that the boundary points wont be iterated over. Hence, central
    # differencing schemes can be used
    # 
    # restartVTK = filename of VTK file to start from
    #
    # saveInterval: Save grid files at how many intervals. Enter negative values
    # if only final grid is needed.

    if restartVTK == None:
        Grid = Grid_obj.Grid
    else:
        Grid = Grid_obj.populate_from_vtk(filename)
        
    imx = Grid_obj.imx
    jmx = Grid_obj.jmx
    iter_no = 0
    
    Grid_obj.save_to_vtk('Elliptic-0')

    # Using grid objects will be slow
    # Hence converting to numpy arrays
    X = np.zeros((imx, jmx))
    Y = np.zeros((imx, jmx))
    Xold = np.zeros((imx, jmx))
    Yold = np.zeros((imx, jmx))
    alpha = np.zeros((imx-2, jmx-2))
    beta = np.zeros((imx-2, jmx-2))
    gamma = np.zeros((imx-2, jmx-2))

    for xi in range(imx):
        for eta in range(jmx):
            X[xi, eta] = Grid[xi, eta].x
            Y[xi, eta] = Grid[xi, eta].y

    resnorm = 1
    resnorm0 = 1
    while resnorm/resnorm0 >= eps:
        iter_no += 1
        Xold = copy.deepcopy(X)
        Yold = copy.deepcopy(Y)
        
        alpha = 0.25 * ((X[1:imx-1, 2:jmx] - X[1:imx-1, 0:jmx-2])**2 + \
                        (Y[1:imx-1, 2:jmx] - Y[1:imx-1, 0:jmx-2])**2)
        gamma = 0.25 * ((X[2:imx, 1:jmx-1] - X[0:imx-2, 1:jmx-1])**2 + \
                        (Y[2:imx, 1:jmx-1] - Y[0:imx-2, 1:jmx-1])**2)
        beta = 0.25 * (((X[2:imx, 1:jmx-1] - X[0:imx-2, 1:jmx-1]) * \
                        (X[1:imx-1, 2:jmx] - X[1:imx-1, 0:jmx-2])) + \
                       ((Y[2:imx, 1:jmx-1] - Y[0:imx-2, 1:jmx-1]) * \
                        (Y[1:imx-1, 2:jmx] - Y[1:imx-1, 0:jmx-2])))

        X[1:imx-1, 1:jmx-1] =  \
                    ((alpha[:, :]*(X[2:imx, 1:jmx-1] + X[0:imx-2, 1:jmx-1])) + \
                    (gamma[:, :]*(X[1:imx-1, 2:jmx] + X[1:imx-1, 0:jmx-2])) - \
                    (0.5 * beta[:, :] * \
                     (X[2:imx, 2:jmx] + X[0:imx-2, 0:jmx-2] - \
                      X[2:imx, 0:jmx-2] - X[0:imx-2, 2:jmx]) \
                    )) / \
                    (2*(alpha[:, :] + gamma[:, :]))
        Y[1:imx-1, 1:jmx-1] =  \
                    ((alpha[:, :]*(Y[2:imx, 1:jmx-1] + Y[0:imx-2, 1:jmx-1])) + \
                    (gamma[:, :]*(Y[1:imx-1, 2:jmx] + Y[1:imx-1, 0:jmx-2])) - \
                    (0.5 * beta[:, :] * \
                     (Y[2:imx, 2:jmx] + Y[0:imx-2, 0:jmx-2] - \
                      Y[2:imx, 0:jmx-2] - Y[0:imx-2, 2:jmx]) \
                    )) / \
                    (2*(alpha[:, :] + gamma[:, :]))

        resnorm = np.linalg.norm(np.sqrt((X - Xold)**2 + (Y - Yold)**2))
        if iter_no == iter_max:
            break
        if iter_no == 1:
            resnorm0 = resnorm
        if iter_no % 200 == 0:
            print 'Iter No.: ', iter_no, '. Resnorm_ratio, Resnorm = ', \
                                                resnorm/resnorm0, ' ', resnorm
        if saveInterval > 0:
            if iter_no % int(saveInterval) == 0:
                Grid_obj.add_to_grid(X, Y)
                Grid_obj.save_to_vtk('Elliptic-' + iter_no.__str__())

    Grid_obj.add_to_grid(X, Y)
    Grid_obj.save_to_vtk('Elliptic-' + iter_no.__str__())
    return Grid_obj
   
    
if __name__ == '__main__':
    Grid_obj = get_univariate_grid()
    #Grid_obj = get_transfinite_interpolated_grid()
    #Grid_obj.save_to_vtk('Elliptic-Initial')
    Grid_obj = elliptic_grid_generation(Grid_obj, eps=0.000001, \
                        iter_max=200000, saveInterval=-1)


