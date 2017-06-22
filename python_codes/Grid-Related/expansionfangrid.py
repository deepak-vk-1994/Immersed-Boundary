from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from Grids import Point, Boundary, TwoDGrid, BdType
import BoundaryGrids
import StructuredGridGeneration

def get_expfan_grid(delta, I):
    # Grid as per ASM2010 AIAA code verification paper
    # Ramp at I/3
    J = int(I / 2)
    x_min = -1
    x_corner = 0
    y = (x_corner - x_min) * 1.5
    x_end_dist = (x_corner - x_min) * 2

    # eta_min
    start_point = Point(x_min, 0)
    corner_point = Point(x_corner, 0)
    end_point = Point(x_end_dist*np.cos(delta), x_end_dist*np.sin(delta))

    r1 = BoundaryGrids.fit_line(start_point, corner_point, int(I/3))
    r2 = BoundaryGrids.fit_line(corner_point, end_point, int(2*I/3))
    r = np.append(r1[:-1], r2)
    bd1 = Boundary(I, BdType.eta_min)
    bd1.populate_boundary_ndar_points(r)

    # eta_max
    start_point = Point(x_min, y)
    corner_point = Point(x_corner, y)
    end_point = Point(x_end_dist*np.cos(delta), x_end_dist*np.sin(delta) + y)

    r1 = BoundaryGrids.fit_line(start_point, corner_point, int(I/3))
    r2 = BoundaryGrids.fit_line(corner_point, end_point, int(2*I/3))
    r = np.append(r1[:-1], r2)
    bd2 = Boundary(I, BdType.eta_max)
    bd2.populate_boundary_ndar_points(r)

    # xi_min
    start_point = Point(x_min, 0)
    end_point = Point(x_min, y)
    r = BoundaryGrids.fit_line(start_point, end_point, J)
    bd3 = Boundary(J, BdType.xi_min)
    bd3.populate_boundary_ndar_points(r)

    
    # xi_max
    start_point = Point(x_end_dist*np.cos(delta), x_end_dist*np.sin(delta))
    end_point = Point(x_end_dist*np.cos(delta), x_end_dist*np.sin(delta) + y)
    r = BoundaryGrids.fit_line(start_point, end_point, J)
    bd4 = Boundary(J, BdType.xi_max)
    bd4.populate_boundary_ndar_points(r)

    Grid_obj = StructuredGridGeneration.get_univariate_grid(bd1, bd2, bd3, bd4)
    Grid_obj.save_to_vtk('expansion_fan.vtk')
    Grid_obj.save_to_txt('expansion_fan.txt')


def get_oblique_shock_grid(delta, I):
    # Grid as per ASM2010 AIAA code verification paper
    # Ramp at I/2
    J = int(I / 2)
    x_min = -1
    x_corner = 0
    y = (x_corner - x_min)
    x_end_dist = (x_corner - x_min)

    # eta_min
    start_point = Point(x_min, 0)
    corner_point = Point(x_corner, 0)
    end_point = Point(x_end_dist*np.cos(delta), x_end_dist*np.sin(delta))

    r1 = BoundaryGrids.fit_line(start_point, corner_point, int(I/2))
    r2 = BoundaryGrids.fit_line(corner_point, end_point, int(I/2))
    r = np.append(r1[:-1], r2)
    bd1 = Boundary(I, BdType.eta_min)
    bd1.populate_boundary_ndar_points(r)
    I = len(bd1.Points) - 1

    # eta_max
    start_point = Point(x_min, y)
    corner_point = Point(x_corner, y)
  # end_point = Point(x_end_dist*np.cos(delta), x_end_dist*np.sin(delta) + y)
    end_point = Point(x_end_dist*np.cos(delta), y)

  # r1 = BoundaryGrids.fit_line(start_point, corner_point, int(I/2))
  # r2 = BoundaryGrids.fit_line(corner_point, end_point, int(I/2))
  # r = np.append(r1[:-1], r2)
    r = BoundaryGrids.fit_line(start_point, end_point, I)
    bd2 = Boundary(I, BdType.eta_max)
    bd2.populate_boundary_ndar_points(r)

    # xi_min
    start_point = Point(x_min, 0)
    end_point = Point(x_min, y)
    r = BoundaryGrids.fit_line(start_point, end_point, J)
    bd3 = Boundary(J, BdType.xi_min)
    bd3.populate_boundary_ndar_points(r)

    
    # xi_max
    start_point = Point(x_end_dist*np.cos(delta), x_end_dist*np.sin(delta))
    end_point = Point(x_end_dist*np.cos(delta), y)
    r = BoundaryGrids.fit_line(start_point, end_point, J)
    bd4 = Boundary(J, BdType.xi_max)
    bd4.populate_boundary_ndar_points(r)

    Grid_obj = StructuredGridGeneration.get_univariate_grid(bd1, bd2, bd3, bd4)
    Grid_obj.save_to_vtk('oblique_shock.vtk')
    Grid_obj.save_to_txt('oblique_shock.txt')


def get_shock_BL_grid(delta, I):
    # Grid as per ASM2010 AIAA code verification paper
    # Ramp at I/2
    J = int(I / 2)
    x_min = 0
    x_corner = 0.1
    x_mid = 0.2 * x_corner
    y = 0.0368
    x_end_dist = x_corner - x_min

    # eta_max
    start_point = Point(x_min, y)
    turning_point = Point(x_mid, (x_mid-x_min)*np.tan(delta) + y)
    end_point = Point(x_corner, -(x_corner-x_mid)*np.tan(delta) + \
                                 turning_point.y)
    r1 = BoundaryGrids.fit_line(start_point, turning_point, int(I*0.2))
    r2 = BoundaryGrids.fit_line(turning_point, end_point, int(I*0.8))
    r = np.append(r1[:-1], r2)
    bd2 = Boundary(I, BdType.eta_max)
    bd2.populate_boundary_ndar_points(r)
    
    I = len(bd2.Points) - 1
    # eta_min
    start_point = Point(x_min, 0.)
    end_point = Point(x_corner, .0)
    r = BoundaryGrids.fit_line(start_point, end_point, I)
    bd1 = Boundary(I, BdType.eta_min)
    bd1.populate_boundary_ndar_points(r)


    dt = 3.0
    # xi_min
    start_point = Point(x_min, 0)
    end_point = Point(x_min, y)
    stretching_f = lambda xi: 1 + (np.tanh(dt*(xi - 1)) / np.tanh(dt))
    r = BoundaryGrids.fit_line(start_point, end_point, J, stretching_f)
    bd3 = Boundary(J, BdType.xi_min)
    bd3.populate_boundary_ndar_points(r)

    
    # xi_max
    start_point = Point(x_corner, 0.)
    end_point = Point(x_corner, -(x_corner-x_mid)*np.tan(delta) + \
                                 turning_point.y)
    r = BoundaryGrids.fit_line(start_point, end_point, J, stretching_f)
    bd4 = Boundary(J, BdType.xi_max)
    bd4.populate_boundary_ndar_points(r)

    Grid_obj = StructuredGridGeneration.get_univariate_grid(bd1, bd2, bd3, bd4, stretching_f)
    Grid_obj.save_to_vtk('shock_BL.vtk')
    Grid_obj.save_to_txt('shock_BL.txt')

def get_trap_grid(I, J):
    # Grid as per ASM2010 AIAA code verification paper
    # eta_min
    start_point = Point(0., 0.)
    end_point = Point(0.3, 0.0)
    r = BoundaryGrids.fit_line(start_point, end_point, I)
    bd1 = Boundary(I, BdType.eta_min)
    bd1.populate_boundary_ndar_points(r)

    # eta_max
    start_point = Point(0., 0.1)
    end_point = Point(0.3, 0.005)
    r = BoundaryGrids.fit_line(start_point, end_point, I)
    bd2 = Boundary(I, BdType.eta_max)
    bd2.populate_boundary_ndar_points(r)

    # xi_min
    start_point = Point(0., 0.)
    end_point = Point(0., 0.1)
    r = BoundaryGrids.fit_line(start_point, end_point, J)
    bd3 = Boundary(J, BdType.xi_min)
    bd3.populate_boundary_ndar_points(r)

    # xi_max
    start_point = Point(0.3, 0.)
    end_point = Point(0.3, 0.005)
    r = BoundaryGrids.fit_line(start_point, end_point, J)
    bd4 = Boundary(J, BdType.xi_max)
    bd4.populate_boundary_ndar_points(r)

    Grid_obj = StructuredGridGeneration.get_univariate_grid(bd1, bd2, bd3, bd4)
    Grid_obj.save_to_vtk('trap.vtk')
    Grid_obj.save_to_txt('trap.txt')
