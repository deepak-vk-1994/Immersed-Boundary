# This module implements classes for grid generation
from __future__ import division
from numbers import Number
from enum import Enum
import numpy as np
import os
import matplotlib.pyplot as plt

BdType = Enum('Bdtype', 'xi_min xi_max eta_min eta_max')

class Point(object):
    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y

    def magnitude(self):
        return np.sqrt(self.x**2 + self.y**2)

    def dot(self, other):
        if isinstance(other, Point):
            return (self.x*other.x + self.y*other.y)
        else:
            return NotImplemented
    
    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        if isinstance(other, Number):
            return Point(self.x*other, self.y*other)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Number):
            return Point(self.x / other, self.y / other)
        else:
            return NotImplemented

    def __pow__(self, other):
        if  isinstance(other, Number):
            return Point(self.x ** other, self.y ** other)
        else:
            return NotImplemented

    def __eq__(self, other):
        return ((self.x == other.x) and (self.y == other.y))


class Boundary(object):
    def __init__(self, n = None, bdtype=None):
        # Number of points = n + 1
        self.n = n
        #vPoint = np.vectorize(Point)
        #self.Points = vPoint(np.empty((n+1, 1), dtype=object))
        self.Points = np.empty(n+1, dtype=object)
        self.bdtype = bdtype

    def populate_boundary_from_ind(self, r_x, r_y):
        for i in range(self.n + 1):
            #self.Points[i].x = r_x[i]
            #self.Points[i].y = r_y[i]
            self.Points[i] = Point(r_x[i], r_y[i])
            
    def populate_boundary_ndar_points(self, r):
        if type(r) == np.ndarray:        
            self.Points = r
        
    def plot_boundary(self, plt, label):
        r_x = []
        r_y = []
        for i in range(self.n + 1):
            r_x.append(self.Points[i].x)
            r_y.append(self.Points[i].y)
        plt.plot(r_x, r_y, marker='*', label = label)
        return plt
        


class TwoDGrid(object):
    def __init__(self, imx = None, jmx = None):
        self.imx = imx
        self.jmx = jmx
        #vPoint = np.vectorize(Point)
        #self.Grid = vPoint(np.empty((imx, jmx), dtype=object))
        if self.imx is not None and self.jmx is not None:
            self.Grid = np.empty((imx, jmx), dtype=object)
    
    def add_to_grid(self, X, Y):
        for xi in range(self.imx):
            for eta in range(self.jmx):
                self.Grid[xi, eta] = Point(X[xi,eta], Y[xi,eta])

    def set_boundary(self, boundary):
        if type(boundary) is not Boundary:
            print 'Invalid boundary passed. Exiting...'
            return

        # In python a[0:n] returns an array from a[0] to a[n-1]
        if boundary.bdtype == BdType.xi_min:
            self.Grid[0, 0:self.jmx] = boundary.Points
        elif boundary.bdtype == BdType.xi_max:
            self.Grid[self.imx-1, 0:self.jmx] = boundary.Points
        elif boundary.bdtype == BdType.eta_min:
            self.Grid[0:self.imx, 0] = boundary.Points
        elif boundary.bdtype == BdType.eta_max:
            self.Grid[0:self.imx, self.jmx-1] = boundary.Points

    def populate_from_vtk(self, filename):
        filename = os.path.join(os.curdir, filename)
        f = open(filename, 'r')
        while f.readline().strip() != '':
            pass
        # Next line contains Dimensions
        line = f.readline()
        line = line.split()
        imx = int(line[1])
        jmx = int(line[2])
        
        self.imx = imx
        self.jmx = jmx
        self.Grid = np.empty((imx, jmx), dtype=object)
        
        # Next line is about number of points - skip items
        f.readline()
        for j in range(jmx):
            for i in range(imx):
                line = f.readline()
                line = line.split()
                x = float(line[0])
                y = float(line[1])
                self.Grid[i, j] = Point(x, y)
        
        f.close()
    
    
    def save_to_vtk(self, filename, comment='Grid Generation'):
        if filename[-4:] != '.vtk':
            filename = os.path.join(os.curdir, filename + '.vtk')
        else:
            filename = os.path.join(os.curdir, filename)
        f = open(filename, 'w')

        f.write(r'# vtk DataFile Version 3.1')
        f.write('\n')
        f.write(comment + '\n')
        f.write('ASCII\n')
        f.write('DATASET STRUCTURED_GRID\n\n')

        f.write('DIMENSIONS ' + self.imx.__str__() + ' ' + \
                                self.jmx.__str__() + ' 1\n')
        f.write('POINTS ' + (self.imx * self.jmx).__str__() + ' FLOAT\n') 

        for j in range(self.jmx):
            for i in range(self.imx):
                f.write(self.Grid[i][j].x.__str__() + ' ' + \
                        self.Grid[i][j].y.__str__() + ' 0.0 \n')
                
        f.write('\n')        
        f.close()

    def save_to_txt(self, filename):
        f = open(filename, 'w')
        f.write(self.imx.__str__() + ' ' + self.jmx.__str__() + '\n')
        for j in range(self.jmx):
            for i in range(self.imx):
                f.write(self.Grid[i][j].x.__str__() + ' ' + \
                        self.Grid[i][j].y.__str__()  + '\n')
        f.close()


    def plot_grid(self, plt, label):
        r_x = []
        r_y = []
        for i in range(self.imx):
            for j in range(self.jmx):
                r_x.append(self.Grid[i][j].x)
                r_y.append(self.grid[i][j].y)
        plt.plot(r_x, r_y, '*', label = label)
        return plt        
