from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from Grids import Point, Boundary, TwoDGrid, BdType
import BoundaryGrids
import StructuredGridGeneration

def get_plane_IB(delta, I, J):
    global plt
    x_min = -1.0
    x_max = 5.0
    y_min = -5.0
    y_max = 5.0

    airfoil_start = x_min + 43.5/100*(x_max - x_min)
    airfoil_end = x_min + 57.5/100*(x_max - x_min)

    airfoil_y_s = y_min + 48/100*(y_max - y_min)
    airfoil_y_e = y_min + 52/100*(y_max - y_min)
    
    dx = 0.01
    dy = 0.01

 #  x_min = -6.0 * 0.01
 #  x_max = 21.0 * 0.01
 #  y_min = -5.0 * 0.01
 #  y_max = 5.0 * 0.01

 #  airfoil_start = x_min + 15.38/100*(x_max - x_min)
 #  airfoil_end = x_min + 26.92/100*(x_max - x_min)
#   airfoil_start = -1 * 0.01
#   airfoil_end = 1 * 0.01

 #  airfoil_y_s = y_min + 40./100*(y_max - y_min)
 #  airfoil_y_e = y_min + 60./100*(y_max - y_min)
#   airfoil_y_s = -2.5 * 0.01
#   airfoil_y_e = 2.5 * 0.01
    
#   dx = 0.02 * 0.01
#   dy = 0.02 * 0.01
    length1 = abs(airfoil_start - y_min)
    length2 = abs(x_max - airfoil_end)
    delta1 = BoundaryGrids.find_hyp_tan_stretching_factor(dx, int(I/4), length1)
    delta2 = BoundaryGrids.find_hyp_tan_stretching_factor(dx, int(3*I/4), length2)
    stretching_f1 = lambda xi: 1 + (np.tanh(delta1*(xi - 1)) / np.tanh(delta1))
    stretching_f2 = lambda xi: 1 + (np.tanh(delta2*(xi - 1)) / np.tanh(delta2))

    end_point = Point(x_min, y_min)
    start_point = Point(airfoil_start, y_min)
    r1_t = BoundaryGrids.fit_line(start_point, end_point, int(I/4), stretching_f1)
    r1 = np.empty(len(r1_t), dtype=object)
    for i in range(len(r1)):
        r1[i] = r1_t[len(r1) - 1 - i]

    num = int((airfoil_end - airfoil_start) / dx) + 1
    start_point = Point(airfoil_start, y_min)
    end_point = Point(airfoil_end, y_min)
    r2 = BoundaryGrids.fit_line(start_point, end_point, num)

    start_point = Point(airfoil_end, y_min)
    end_point = Point(x_max, y_min)
    r3 = BoundaryGrids.fit_line(start_point, end_point, int(3*I/4), stretching_f2)

    r = np.append(r1[:-1], r2)
    r = np.append(r[:-1], r3)
    bd1 = Boundary(len(r)-1, BdType.eta_min)
    bd1.populate_boundary_ndar_points(r)

    f = plt.figure()
    plt = bd1.plot_boundary(plt, 'eta_min')

    end_point = Point(x_min, y_max)
    start_point = Point(airfoil_start, y_max)
    r1_t = BoundaryGrids.fit_line(start_point, end_point, int(I/4), stretching_f1)
    r1 = np.empty(len(r1_t), dtype=object)
    for i in range(len(r1)):
        r1[i] = r1_t[len(r1) - 1 - i]

    num = int((airfoil_end - airfoil_start) / dx) + 1
    start_point = Point(airfoil_start, y_max)
    end_point = Point(airfoil_end, y_max)
    r2 = BoundaryGrids.fit_line(start_point, end_point, num)

    start_point = Point(airfoil_end, y_max)
    end_point = Point(x_max, y_max)
    r3 = BoundaryGrids.fit_line(start_point, end_point, int(3*I/4), stretching_f2)

    r = np.append(r1[:-1], r2)
    r = np.append(r[:-1], r3)
    bd2 = Boundary(len(r)-1, BdType.eta_max)
    bd2.populate_boundary_ndar_points(r)
    plt = bd2.plot_boundary(plt, 'eta_max')

    I = len(r) - 1

#   start_point = Point(x_min, y_min)
#   end_point = Point(x_min, y_max)
#   A = 3
#   L = 1
#   yc = 0.5
#   stretching_f = lambda xi: (L*xi) + A*(yc - L*xi)*(1 - xi)*xi
#   r = BoundaryGrids.fit_line(start_point, end_point, int(J), stretching_f)
#   bd3 = Boundary(len(r)-1, BdType.xi_min)
#   bd3.populate_boundary_ndar_points(r)
#   plt = bd3.plot_boundary(plt, 'xi_min')

#   start_point = Point(x_max, y_min)
#   end_point = Point(x_max, y_max)
#   A = 3
#   L = 1
#   yc = 0.5
#   stretching_f = lambda xi: (L*xi) + A*(yc - L*xi)*(1 - xi)*xi
#   r = BoundaryGrids.fit_line(start_point, end_point, int(J), stretching_f)
#   bd4 = Boundary(len(r)-1, BdType.xi_max)
#   bd4.populate_boundary_ndar_points(r)
#   plt = bd4.plot_boundary(plt, 'xi_max')

 #  plt.show()

    num = int((airfoil_y_e - airfoil_y_s)/ dy) + 1
    K = int(J/2) + int(J/2) + num
    
    Grid_obj = TwoDGrid(I+1, K+1)
    Grid_obj.set_boundary(bd1)
    Grid_obj.set_boundary(bd2)
    Grid = Grid_obj.Grid

    length1 = abs(airfoil_y_s - y_min)
    length2 = abs(y_max - airfoil_y_e)
    delta1 = BoundaryGrids.find_hyp_tan_stretching_factor(dy, int(J/2), length1)
    delta2 = BoundaryGrids.find_hyp_tan_stretching_factor(dy, int(J/2), length2)
    stretching_f1 = lambda xi: 1 + (np.tanh(delta1*(xi - 1)) / np.tanh(delta1))
    stretching_f2 = lambda xi: 1 + (np.tanh(delta2*(xi - 1)) / np.tanh(delta2))

    for xi in range(0, I+1):
        P_eta_min = Grid[xi, 0]
        P_eta_max = Grid[xi, K]
        
        end_point = P_eta_min
        start_point = Point(P_eta_min.x, airfoil_y_s)
        r1_t = BoundaryGrids.fit_line(start_point, end_point, int(J/2), stretching_f1)
        r1 = np.empty(len(r1_t), dtype=object)
        for i in range(len(r1)):
            r1[i] = r1_t[len(r1) - 1 - i]
        
    #   num = int((airfoil_y_e - airfoil_y_s)/ 0.005) + 1
        start_point = Point(P_eta_min.x, airfoil_y_s)
        end_point = Point(P_eta_max.x, airfoil_y_e)
        r2 = BoundaryGrids.fit_line(start_point, end_point, num)

        start_point = Point(P_eta_max.x, airfoil_y_e)
        end_point = P_eta_max
        r3 = BoundaryGrids.fit_line(start_point, end_point, int(J/2), stretching_f2)

        r = np.append(r1[:-1], r2)
        r = np.append(r[:-1], r3)
     #  K = int(J/2) + int(J/2) + num
    #   print len(r1), len(r2), len(r3), len(r), J, K
        
        Grid[xi, :] = r

  # Grid_obj = StructuredGridGeneration.get_univariate_grid(bd1, bd2, \
  #         bd3, bd4, stretching_f)
    Grid_obj.save_to_vtk('Circle-IB.vtk', 'Circle-IB')
 
    Grid_obj.save_to_txt('Circle-IB.txt')
    

def get_exp_plane_grid(I, J):
    delta = -10 * np.pi / 180
    x_min = -1
    x_corner = 0
    y_max = (x_corner - x_min) * 1.5
    x_end_dist = (x_corner - x_min) * 2

    x_max = x_end_dist * np.cos(delta)
    y_min = x_end_dist * np.sin(delta)

    x_min = 0
    x_max = 1
    y_min = 0
    y_max = 0.2

    start_point = Point(x_min, y_min)
    end_point = Point(x_max, y_min)
    r = BoundaryGrids.fit_line(start_point, end_point, int(I))
    bd1 = Boundary(len(r)-1, BdType.eta_min)
    bd1.populate_boundary_ndar_points(r)
    
    start_point = Point(x_min, y_max)
    end_point = Point(x_max, y_max)
    r = BoundaryGrids.fit_line(start_point, end_point, int(I))
    bd2 = Boundary(len(r)-1, BdType.eta_max)
    bd2.populate_boundary_ndar_points(r)
    
    start_point = Point(x_min, y_min)
    end_point = Point(x_min, y_max)
    r = BoundaryGrids.fit_line(start_point, end_point, int(J))
    bd3 = Boundary(len(r)-1, BdType.xi_min)
    bd3.populate_boundary_ndar_points(r)
    
    start_point = Point(x_max, y_min)
    end_point = Point(x_max, y_max)
    r = BoundaryGrids.fit_line(start_point, end_point, int(J))
    bd4 = Boundary(len(r)-1, BdType.xi_max)
    bd4.populate_boundary_ndar_points(r)
    
    Grid_obj = StructuredGridGeneration.get_univariate_grid(bd1, bd2, \
            bd3, bd4)
 #  Grid_obj.save_to_vtk('ExpFan-IB.vtk', 'ExpFan-IB')
 
#   Grid_obj.save_to_txt('ExpFan-IB.txt')
 #  Grid_obj.save_to_vtk('RampChannel-IB.vtk', 'RampChannel-IB')
 
    Grid_obj.save_to_txt('Rchannel-IB-2.txt')
