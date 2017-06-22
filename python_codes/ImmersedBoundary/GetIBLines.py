from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from Grids import Point, Boundary, TwoDGrid, BdType
from BoundaryGrids import fit_cubic_spline, extract_control_points


# This is a script to generate the required IB lines in the following format
# Each line should contain the start point S, end point E and the outward unit normal vector N
# Hence, the ith line should contain
# Si_x Si_y Ei_x Ei_y Ni_x Ni_y

# Note that E_i = S_i+1

def put_to_file(start_points, end_points, normals, filename):
    # start_points, end_points and normals are vectors -> can be lists
    # or numpy arrays. Should be of the same length.
    # If 2D, then each entry of start_points should be a vector of length 2
    f = open(filename, 'w')
    n = len(start_points)
    f.write(n.__str__() + '\n')
    for i in range(n):
        sp = start_points[i]
        ep = end_points[i]
        nv = normals[i]
        line = list(sp) + list(ep) + list(nv)
        for l in line:
            f.write(l.__str__() + ' ')
        f.write('\n')
    f.close()


def get_IB_naca_4_dig(airfoil_str, chord_len= 1.0, num_points=51):
    # Input: airfoil_str: String of 4 numbers
    #        chord_len: length of chord (float or int)
    # Only implemented for symmetric 4 digit
    #TODO: Cambered 4 digit
    #TODO: Validation of input

    t = int(airfoil_str[-2:]) * 0.01
    c = chord_len
    x = np.linspace(0, c, num_points)
    x = 0.5 * c * (1 - np.cos(x/c * np.pi))
    y = 5 * t * c * \
        (0.2969*np.sqrt(x/c) - 0.126*(x/c) - 0.3516*(x/c)**2 + \
         0.2843*(x/c)**3 - 0.1036*(x/c)**4)
    return (x, y)


def get_lines_from_upper_airfoil(x, y):    
    # Now to create the start, end and normals list
    # Given n points, the number of line segments is n-1
    # But we want a closed airfoil with the trailing edge matching
    # Hence number of line segments = 2*(n-1)
    # But total number of points = 2*n - 2
    # We need outward normal
    start = []
    end = []
    normals = []

    # First storing the top surface line segments
    for i in range(len(x) - 1):
        start.append([x[i], y[i]])
        end.append([x[i+1], y[i+1]])
        # Given a line segment, find dx = x2 - x1 and dy = y2 - y1
        # The normal vector to the line segment is given by (dy, -dx)
        # We also need to normalise it
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        ds = np.sqrt(dx**2 + dy**2)
        dx = dx / ds
        dy = dy / ds
        normals.append([-dy, dx])

    # Now we have the top surface line segments from LE to TE
    # Now the bottom surface line segments from TE to LE
    
    for i in range(len(x) - 1, 0, -1):
        start.append([x[i], -y[i]])
        end.append([x[i-1], -y[i-1]])
        # Given a line segment, find dx = x2 - x1 and dy = y2 - y1
        # The normal vector to the line segment is given by (dy, -dx)
        # We also need to normalise it
        dx = x[i-1] - x[i]
        dy = -y[i-1] + y[i]
        ds = np.sqrt(dx**2 + dy**2)
        dx = dx / ds
        dy = dy / ds
        normals.append([-dy, dx])

    return (start, end, normals)


def get_concave_polygon():
    x = [0, 0, 1, 2, 2, 0]
    y = [0, 2, 1, 2, 0, 0]
    start = []
    end = []
    normals = []

    # First storing the top surface line segments
    for i in range(len(x) - 1):
        start.append([x[i], y[i]])
        end.append([x[i+1], y[i+1]])
        # Given a line segment, find dx = x2 - x1 and dy = y2 - y1
        # The normal vector to the line segment is given by (dy, -dx)
        # We also need to normalise it
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        ds = np.sqrt(dx**2 + dy**2)
        dx = dx / ds
        dy = dy / ds
        normals.append([-dy, dx])
    return (start, end, normals)


def get_ramp():
  # x = [-5, 5, 10.5]
  # y = [-5, -3, -3]
    x = [0.25, 0.5, 1.0]
    y = [0, 0.0667, 0.05]
    start = []
    end = []
    normals = []

    # First storing the top surface line segments
    for i in range(len(x) - 1):
        start.append([x[i], y[i]])
        end.append([x[i+1], y[i+1]])
        # Given a line segment, find dx = x2 - x1 and dy = y2 - y1
        # The normal vector to the line segment is given by (dy, -dx)
        # We also need to normalise it
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        ds = np.sqrt(dx**2 + dy**2)
        dx = dx / ds
        dy = dy / ds
        normals.append([-dy, dx])
    return (start, end, normals)


def get_circle(radius, centre, I):
    r = np.empty(I+1, dtype=object)
    for xi in range(I+1):
        t =  2 * np.pi * (1 - xi/ I)
        x = centre.x + radius*np.cos(t)
        y = centre.y + radius*np.sin(t)
        r[xi] = Point(x, y)

    start = []
    end = []
    normals = []
    for i in range(len(r) - 1):
        start.append([r[i].x, r[i].y])
        end.append([r[i+1].x, r[i+1].y])
        # Given a line segment, find dx = x2 - x1 and dy = y2 - y1
        # The normal vector to the line segment is given by (dy, -dx)
        # We also need to normalise it
        dx = r[i+1].x - r[i].x
        dy = r[i+1].y - r[i].y
        ds = np.sqrt(dx**2 + dy**2)
        dx = dx / ds
        dy = dy / ds
        normals.append([-dy, dx])
    return (start, end, normals)

def get_naca_grid_gridgen():
    # Use points from naca equation as control points so that the trailing
    # edge is sharp
    x, y = get_IB_naca_4_dig('0012', chord_len= 1.0, num_points=66)
    CP = []
    for i in range(len(x)):
        CP.append(Point(x[i], y[i]))

    ################ Generating around airfoil
    I = 100
#   xi_spacing_f = lambda xi, I: 0.5 * I * (1 - np.cos(xi * np.pi / I))
    xi_spacing_f = lambda xi, I: xi
    s_CP_f = lambda j,n,I: j * I / n
    s_CP_inverse_f = lambda s,n,I: s * n / I 

    r = fit_cubic_spline(CP, I, xi_spacing_f, s_CP_f, s_CP_inverse_f)
    bd1 = Boundary(I, BdType.eta_min)
    bd1.populate_boundary_ndar_points(r)
    x = []
    y = []
    for i in range(len(r)):
        x.append(r[i].x)
        y.append(r[i].y)
    return (x, y)        


def get_airfoil_dat_gridgen():
    xt, yt = zip(*np.loadtxt('top'))
    xb, yb = zip(*np.loadtxt('bottom'))
    CPt = []
    CPb = []
    for i in range(len(xt)):
        CPt.append(Point(xt[i], yt[i]))
    for i in range(len(xb)):
        CPb.append(Point(xb[i], yb[i]))

    I = 200
#   xi_spacing_f = lambda xi, I: 0.5 * I * (1 - np.cos(xi * np.pi / I))
    xi_spacing_f = lambda xi, I: xi
    s_CP_f = lambda j,n,I: j * I / n
    s_CP_inverse_f = lambda s,n,I: s * n / I 

    rt = fit_cubic_spline(CPt, I, xi_spacing_f, s_CP_f, s_CP_inverse_f)
    rb = fit_cubic_spline(CPb, I, xi_spacing_f, s_CP_f, s_CP_inverse_f)
    bd1 = Boundary(I, BdType.eta_min)
    bd1.populate_boundary_ndar_points(rt)
    bd2 = Boundary(I, BdType.eta_min)
    bd2.populate_boundary_ndar_points(rb)
    xtnew = []
    ytnew = []
    xbnew = []
    ybnew = []
    for i in range(len(rt)):
        xtnew.append(rt[i].x)
        ytnew.append(rt[i].y)
    for i in range(len(rb)):
        xbnew.append(rb[i].x)
        ybnew.append(rb[i].y)

    np.savetxt('top-new', zip(xtnew, ytnew))
    np.savetxt('bottom-new', zip(xbnew, ybnew))

    plt.plot(xtnew, ytnew, '*-')
    plt.plot(xbnew, ybnew, '.-')
    plt.axis('equal')
    plt.show()

def get_exp_fan_IB():
    delta = -10*np.pi/180
    x_st = -1
    x_corner = 0
    x_end_dist = (x_corner - x_st) * 2
    x_end = x_end_dist*np.cos(delta)
    y_st = 0.0
    y_corner = 0.0
    y_end = x_end_dist*np.sin(delta)

    x = [x_st, x_corner, x_end]
    y = [y_st, y_corner, y_end]
    start = []
    end = []
    normals = []

    # First storing the top surface line segments
    for i in range(len(x) - 1):
        start.append([x[i], y[i]])
        end.append([x[i+1], y[i+1]])
        # Given a line segment, find dx = x2 - x1 and dy = y2 - y1
        # The normal vector to the line segment is given by (dy, -dx)
        # We also need to normalise it
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        ds = np.sqrt(dx**2 + dy**2)
        dx = dx / ds
        dy = dy / ds
        normals.append([-dy, dx])
    return (start, end, normals)


def plot_start_end_normals(start, end, normals):
    # Works only for 2D!!
    #TODO: 3D extension??
    ax = plt.axes()
    scale_to_div = 20
    for i in range(len(start)):
        st = start[i]
        en = end[i]
        no = normals[i]
        no = [(no[0]/scale_to_div), (no[1]/scale_to_div)]

        plt.plot(st[0], st[1], 'ko')
        plt.plot(en[0], en[1], 'k*')
        plt.plot([st[0], en[0]], [st[1], en[1]], 'b')
        cen = [0.5*(st[0] + en[0]), 0.5*(st[1] + en[1])]
        ax.arrow(cen[0], cen[1], no[0], no[1], \
          head_width=0.002, head_length=0.003)
    
    # Not setting axes equal will make normals look like not perpendicular to line segment
    plt.axis('equal')
    plt.show()


def read_from_file(filename):
    start = []
    end = []
    normals = []
    f = open(filename, 'r')
    # First line contains number of iblines
    num_iblines = int(f.readline())
    for i in range(num_iblines):
        line = f.readline()
        data = line.split()
        data = [float(d) for d in data]
        start.append([data[0], data[1]])
        end.append([data[2], data[3]])
        normals.append([data[4], data[5]])

    return (start, end, normals)


if __name__ == '__main__':
 #  x, y = get_naca_grid_gridgen()
 #  start, end, normals = get_lines_from_upper_airfoil(x, y)
    radius = 0.005
    centre = Point(0., 0.)
    I = 50
 #  start, end, normals = get_circle(radius, centre, I)
 #  start, end, normals = get_exp_fan_IB()
    start, end, normals = get_concave_polygon()
    plot_start_end_normals(start, end, normals)
    filename = 'concave.dat'
    put_to_file(start, end, normals, filename)
