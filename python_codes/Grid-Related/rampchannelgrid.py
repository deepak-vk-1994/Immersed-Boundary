# Ramp channel grid
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from Grids import Point, Boundary, TwoDGrid, BdType
import BoundaryGrids
import StructuredGridGeneration

def get_ramp_grid(I, J):
# Grid as per ASM2010 AIAA code verification paper
    # Ramp at I/3
    x_bl = 0.
    y_bl = 0.
    x_rs = 0.25
    y_rs = 0.0
    x_re = 0.5
    y_re = 0.067
    x_br = 1.0
    y_br = 0.067
    x_tl = 0.0
    y_tl = 0.2
    x_tr = 1.0
    y_tr = 0.2
    
    # eta_min
    start_point = Point(x_bl, y_bl)
    end_point = Point(x_rs, y_rs)
    r1 = BoundaryGrids.fit_line(start_point, end_point, int(I/4))

    start_point = Point(x_rs, y_rs)
    end_point = Point(x_re, y_re)
    r2 = BoundaryGrids.fit_line(start_point, end_point, int(I/4))

    start_point = Point(x_re, y_re)
    end_point = Point(x_br, y_br)
    r3 = BoundaryGrids.fit_line(start_point, end_point, int(I/2))
    
    r = np.append(r1[:-1], r2)
    r = np.append(r[:-1], r3)
    bd1 = Boundary(I, BdType.eta_min)
    bd1.populate_boundary_ndar_points(r)
    I = len(bd1.Points) - 1
    
    
    # eta_max
    start_point = Point(x_tl, y_tl)
    end_point = Point(x_tr, y_tr)
    r = BoundaryGrids.fit_line(start_point, end_point, I)
    bd2 = Boundary(I, BdType.eta_max)
    bd2.populate_boundary_ndar_points(r)

    # xi_min
    start_point = Point(x_bl, y_bl)
    end_point = Point(x_tl, y_tl)
    r = BoundaryGrids.fit_line(start_point, end_point, J)
    bd3 = Boundary(J, BdType.xi_min)
    bd3.populate_boundary_ndar_points(r)

    
    # xi_max
    start_point = Point(x_br, y_br)
    end_point = Point(x_tr, y_tr)
    r = BoundaryGrids.fit_line(start_point, end_point, J)
    bd4 = Boundary(J, BdType.xi_max)
    bd4.populate_boundary_ndar_points(r)
    
    print len(bd1.Points), len(bd2.Points), len(bd3.Points), len(bd4.Points)

    Grid_obj = StructuredGridGeneration.get_univariate_grid(bd1, bd2, bd3, bd4)
    Grid_obj.save_to_vtk('ramp_channel_uniform.vtk')
    Grid_obj.save_to_txt('ramp_channel_uniform.txt')