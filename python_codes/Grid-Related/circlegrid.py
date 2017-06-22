from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from Grids import Point, Boundary, TwoDGrid, BdType
import BoundaryGrids
import StructuredGridGeneration

def get_circle_boundary(radius, centre, bd_eta_max):
    # centre is an object of point data type
    I = bd_eta_max.n
    r = np.empty(I+1, dtype=object)

    # Find point on circles such that line from circle
    # to eta_max boundary is normal to circle
    # theta = arctan((y - yc) / (x - xc))
    points = bd_eta_max.Points
    for xi in range(I+1):
        point = points[xi]
        xp, yp = point.x, point.y
        dist = np.sqrt((xp - centre.x)**2 + (yp - centre.y)**2)
        costheta = (xp - centre.x) / dist
        sintheta = (yp - centre.y) / dist
        x = centre.x + radius*costheta
        y = centre.y + radius*sintheta
        r[xi] = Point(x, y)

 #  x = np.empty(I+1)
 #  y = np.empty(I+1)
 #  for i in range(I+1):
 #      x[i] = r[i].x
 #      y[i] = r[i].y

 #  plt.plot(x, y, 'k.')
 #  plt.show()
        
    bd = Boundary(I, BdType.eta_min)
    bd.populate_boundary_ndar_points(r)
    return bd


def get_xi_boundaries(start, end, I):
    # For the flow past a cylinder, the xi boundaries are the same line
    # The points are stretched closer to the start
    delta = 4.5 # trial and error?
    stretching_f = lambda xi: 1 + (np.tanh(delta*(xi - 1)) / np.tanh(delta))
    r = BoundaryGrids.fit_line(start, end, I, stretching_f)
    bd2 = Boundary(I, BdType.xi_min)
    bd2.populate_boundary_ndar_points(r)
    bd3 = Boundary(I, BdType.xi_max)
    bd3.populate_boundary_ndar_points(r)
    return (bd2, bd3)


def get_eta_max_boundary(p1, p2, p3, p4, I):
    # The eta_max boundary is a rectangle with 4 points.
    # Specify from top right point anti-clockwise
    # We also need to return the indices of the array where 
    # p = pi, i = 1 to 4

    # The actual starting point is midpoint of p1 and p4
    start_point = 0.5 * (p1 + p4)

    # Totally 5 segments
    # L1: start to p1 with I / 4 points
    r1 = BoundaryGrids.fit_line(start_point, p1, int(I/4))
    # L2: p1 to p2 with 3*I/16 points
    r2 = BoundaryGrids.fit_line(p1, p2, int(3*I/16))
    # L3: p2 to p3 with I/8 points
    r3 = BoundaryGrids.fit_line(p2, p3, int(I/8))
    # L4: p3 to p4 with 3*I/16 points
    r4 = BoundaryGrids.fit_line(p3, p4, int(3*I/16))
    # L5: p4 to start with I/4 points
    r5 = BoundaryGrids.fit_line(p4, start_point, int(I/4))

    # Combine them (numpy.append)
    # Skip p1 common
    r = np.append(r1[:-1], r2)
    # Skip p2 common
    r = np.append(r[:-1], r3)
    # Skip p3 common
    r = np.append(r[:-1], r4)
    # Skip p4 common
    r = np.append(r[:-1], r5)
    
    # I is number of number of grid points - 1
    I = len(r) - 1
    
    bd = Boundary(I, BdType.eta_max)
    bd.populate_boundary_ndar_points(r)

    for i in range(len(r)):
        if r[i] == p1:
            i1 = i
        elif r[i] == p2:
            i2 = i
        elif r[i] == p3:
            i3 = i
        elif r[i] == p4:
            i4 = i
    i_list = [i1, i2, i3, i4]

    return (bd, i_list)


def get_boundary_grids():
    # Eta_max boundary
    I = 300
    radius = 0.001
    p1 = Point(21*radius, 5*radius)
    p2 = Point(-6*radius, 5*radius)
    p3 = Point(-6*radius, -5*radius)
    p4 = Point(21*radius, -5*radius)
    bd_eta_max, i_list = get_eta_max_boundary(p1, p2, p3, p4, I)
    I = bd_eta_max.n

    # Eta_min boundary
    centre = Point(0, 0)
    bd_eta_min = get_circle_boundary(radius, centre, bd_eta_max)

    # Xi min and max boundaries
    I = 200
    start = Point(centre.x + radius, centre.y)
    end = 0.5 * (p1 + p4)
    bd_xi_min, bd_xi_max = get_xi_boundaries(start, end, I)

    print 'Index list: ', i_list
    return (bd_xi_min, bd_xi_max, bd_eta_min, bd_eta_max, i_list)


def get_circle_grid():
    bd_xi_min, bd_xi_max, bd_eta_min, bd_eta_max, i_list = \
        get_boundary_grids()

    print bd_xi_min.n, bd_xi_max.n, bd_eta_min.n, bd_eta_max.n
    print len(bd_xi_min.Points), len(bd_xi_max.Points), \
          len(bd_eta_min.Points), len(bd_eta_max.Points)

    delta = 3.0 # From Previous section. Clustering close to start
    stretching_f = lambda xi: 1 + (np.tanh(delta*(xi - 1)) / np.tanh(delta))
    Grid_obj = StructuredGridGeneration.get_univariate_grid(bd_xi_min, \
            bd_xi_max, bd_eta_min, bd_eta_max, stretching_f)
    Grid_obj.save_to_vtk('Univariate-circle.vtk', \
                          'Point indices: ' +  i_list.__str__())
    Grid_obj.save_to_txt('Univariate-circle.txt')

if __name__ == '__main__':
    get_circle_grid()
