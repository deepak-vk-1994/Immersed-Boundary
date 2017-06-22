from __future__ import division
import numpy as np
import random
import matplotlib.pyplot as plt
from Grids import Point, Boundary, TwoDGrid, BdType


#  TDMA solver using Thomas Algorithm
def Thomas_algorithm(M, b):
    # M is a square matrix which is tri-diagonal
    # b is the RHS of size (n, 1) or (1, n)
    # No checks are performed to check if given matrix is tri-diagonal or not
    # User needs to pass Tri diagonal matrix, else calculated solution will not
    # be correct
    
    # Note: Has been verified to be working
    
    n = len(M)
    a = np.zeros((n))
    d = np.zeros((n))
    c = np.zeros((n))
    
    # Solution
    X = np.zeros((n))
    
    for i in range(n):
        a[i] = M[i, i+1] if i < n-1 else 0
        d[i] = M[i, i]
        c[i] = M[i, i-1] if i > 0 else 0
        
    # Forward elimination
    for i in range(1, n):
        m = c[i] / d[i-1]
        d[i] = d[i] - (m * a[i-1])
        b[i] = b[i] - (m * b[i-1])
    
    # Reverse substitution
    X[n-1] = b[n-1] / d[n-1]
    
    for i in range(n-2, -1, -1):
        X[i] = (b[i] - (a[i] * X[i+1])) / d[i]
    
    return X
    
    
def extract_control_points(filename):
    f = open(filename, 'r')
    f.readline() # First line - name of airfoil
    
    n = int(float(f.readline().split()[0]))
    # n is the number of control points or airfoil coordinates given
    
    f.readline() # Blank line
    
    # Start reading coordinates
    p_list = []
    for i in range(n):
        x, y = f.readline().split()
        x, y = float(x), float(y)
        p_list.append(Point(x,y))

    return p_list
    
    
def fit_cubic_spline(CP, I, xi_spacing_f=lambda xi: xi, s_CP_f=lambda j,n,I: \
                     j*I/n, s_CP_inverse_f=lambda s,n,I: s*n/I):
    # CP is the list of control points
    # It should be a list of dictionaries with keys 'x' and 'y'
    
    # I + 1 is the number of grid points
    # Note: xi ranges from 0 to I and is an integer only
    
    # xi_spacing_f is a function handle that transforms the spacing between 
    # consecutive xi. This should return values in the same range as that of xi.
    # The default argument is equal spacing.
    # Note that now, s = xi_spacing(xi). Hence the grid generation parameter
    # is 's' and not xi. However since xi is an integer only, it can be used
    # as array index.
    
    # s_CP_f is a function handle which gives the value of s given the index of 
    # control point j, i.e., it return s_j. The values of s at the intermediate
    # control points s_j can be got from equal spacing between s_0 and s_I.

    # s_CP_inverse_f is the function that returns the value of 'j' given a value
    # s. The s_CP_f function can be though of a map between the 'j' space and
    # the 's' space (0<=j<=n and 0<=s<=I). This function is the reverse map.
    # Non integer values of 'j' dont have any physical meaning. But the floor
    # and ceil of j correspond to s_j and s_j+1 which give the location of the
    # j and j+1 th control point respectively. This can be used to give the 
    # equation of the jth cubic spline that was generated
    
    n = len(CP)
    
    # Solving for second derivative of r, implicitly
    r_dd_CP_x = np.zeros((n))
    r_dd_CP_y = np.zeros((n))
    r_dd_CP_x[0] = 0
    r_dd_CP_x[n-1] = 0
    r_dd_CP_y[0] = 0
    r_dd_CP_y[n-1] = 0

    # M is the same for x or y coordinate
    M = np.zeros((n-2, n-2)) # No need to solve for first and last r_dd
    
    b_x = np.zeros((n-2)) # RHS of the implicit equation
    b_y = np.zeros((n-2)) # RHS of the implicit equation
    
    delta = lambda j: s_CP_f(j+1, n-1, I) - s_CP_f(j, n-1, I)
    
    for j in range(n-2):
        # The delta function is written taking into account that the value of j
        # ranges from 0 to n-1. The inner points hence take value of j such that
        # 1 <= j <= n-2. But the computational domain is from 0 <= j <= n-3 
        # (array index).
        # Hence all delta function calls should be made with j -> j+1
        # Also, all calls to CP should also be made with j -> j+1
        if j > 0:
            M[j, j-1] = delta(j-1 + 1) / 6.0                 # c_i
        M[j, j] = (delta(j-1 + 1) + delta(j + 1)) / 3.0      # d_i
        if j < n-3:
            M[j, j+1] = delta(j + 1) / 6.0                   # a_i
        
        b_x[j] = (-1 * (CP[j + 1].x - CP[j-1 + 1].x) / delta(j-1 + 1)) + \
                 (+1 * (CP[j+1 + 1].x - CP[j + 1].x) / delta(j + 1))
        b_y[j] = (-1 * (CP[j + 1].y - CP[j-1 + 1].y) / delta(j-1 + 1)) + \
                 (+1 * (CP[j+1 + 1].y - CP[j + 1].y) / delta(j + 1))
    
    # Now use TDMA solver
    r_dd_CP_x[1:n-1] = Thomas_algorithm(M, b_x)
    r_dd_CP_y[1:n-1] = Thomas_algorithm(M, b_y)


    r_dd_CP = np.empty(n, dtype=object)
    for i in range(n):
        r_dd_CP[i] = Point(r_dd_CP_x[i], r_dd_CP_y[i])

    
    # Now we can construct the cubic spline
    # Note, we will construct the cubic spline at integer values of xi.
    # The equivalent 's' is found as xi_spacing_f(xi).
    # The values of s at the control points are found using the function
    # s_CP_f. 
    
    # Cubic polynomials have been constructed between the control points
    # But, given a value of s, it is necessary to know which interval that s
    # belongs to, i.e., which 'j'
    # That shall be the inverse of the function s_CP_f which is s_CP_inverse_f. 
    # But the inverse will return a fractional value. Hence, the floor and ceil 
    # of that number will give the required interval    
    r = np.empty(I+1, dtype=object)    
    r[0] = Point(CP[0].x, CP[0].y)
    r[I] = Point(CP[n-1].x, CP[n-1].y)        

    for xi in range(1, I):
        s = xi_spacing_f(xi, I)
        j = s_CP_inverse_f(s, n-1, I)
        if int(round(j*10)) % 10 == 0:
            j_l = int(j)
            if j_l == n:
                j_l = n-1
            j_u = j_l + 1
        else:
            j_l = int(np.floor(j))
            j_u = int(np.ceil(j))
        s_l = s_CP_f(j_l, n-1, I)
        s_u = s_CP_f(j_u, n-1, I)
        
        # xi is used as the array index since it is always an integer
        r[xi] = (s_u - s)**3 * r_dd_CP[j_l] / (6.0 * (s_u - s_l)) + \
                 (s - s_l)**3 * r_dd_CP[j_u] / (6.0 * (s_u - s_l)) + \
                 ( (s_u - s) * \
                   ((CP[j_l] / (s_u - s_l)) - ((s_u - s_l) * \
                                                      r_dd_CP[j_l]/6.0)) ) + \
                 ( (s - s_l) * \
                   ((CP[j_u] / (s_u - s_l)) - ((s_u - s_l) * \
                                                      r_dd_CP[j_u]/6.0)) )

    return r


def fit_line_old(start, end, I, stretching_f=lambda xi: xi):
    # I + 1 is the number of grid points
    # start and stop are 2 dictionaries have x and y as keys
    r = np.empty(I+1, dtype=object)

    for xi in range(I+1):
        s = stretching_f(xi/I)
        print start['x']
        print xi
        x = ((1 - s) * start['x']) + (s * end['x'])
        y = ((1 - s) * start['y']) + (s * end['y'])
        r[xi] = Point(x, y)

    return r

def fit_line(start, end, I, stretching_f=lambda xi: xi):
    # I + 1 is the number of grid points
    # start and stop are objects of point type
    r = np.empty(I+1, dtype=object)

    for xi in range(I+1):
        s = stretching_f(xi/I)
        r[xi] = ((1 - s) * start) + (s * end)

    return r

def find_hyp_tan_stretching_factor(ds, I, length):
    # ds = mag(r(xi = 1) - r(xi = 0))
    # length = mag(r2 - r1)
    # Brute force method - Can improve
    # deltaS = dr/dxi at xi = 0
    delta = 0.01
    dd = 0.01 # Change in delta
    
    # As delta increases, ds calculated from the stretching function
    # decreases
    ds_func = 1

    xi = 1
    while ds_func >= ds:
        delta += dd
        ds_func = length * delta / (np.tanh(delta) * I) * (1 - \
                                    (np.tanh(delta*(-1)))**2)

    return delta


def fit_semicircle(Radius, t1, t2, centre, I, stretching_f=lambda xi: xi):
    # Assuming 2D. centre is a dictionary with 'x' and 'y' as keys
    # t1 and t2 - theta ranges - in radians
    r = np.empty(I+1, dtype=object)

    for xi in range(I+1):
        s = stretching_f(xi)
        t = (s * (t2 - t1) / I) + t1
        x = centre['x'] + Radius*np.cos(t)
        y = centre['y'] + Radius*np.sin(t)
        r[xi] = Point(x, y)
    
    return r
    

def get_Boundaries_Problem1():
    # Grid points go from xi = 0 to I. Hence I+1 grid points
    # n + 1 = number of control points

    global plt
    f = plt.figure()
    
    ################ Extract grid
    filename = 'n0012.dat'
    CP = extract_control_points(filename)
    CP_x = []
    CP_y = []
    
    for p in CP:
        CP_x.append(p.x)
        CP_y.append(p.y)
    

    ################ Generating around airfoil
    I = 100
    xi_spacing_f = lambda xi, I: 0.5 * I * (1 - np.cos(xi * np.pi / I))
    s_CP_f = lambda j,n,I: j * I / n
    s_CP_inverse_f = lambda s,n,I: s * n / I 

    r = fit_cubic_spline(CP, I, xi_spacing_f, s_CP_f, s_CP_inverse_f)
    bd1 = Boundary(I, BdType.eta_min)
    bd1.populate_boundary_ndar_points(r)
    
    #plt = bd1.plot_boundary(plt, '1 - Control points')    
    #plt.plot(CP_x, CP_y, marker='.', markevery=3, label = '1 - Control points')

    chord = np.sqrt((CP_x[-1] - CP_x[0])**2 + (CP_y[-1] - CP_y[0])**2)


    ################ Generating the lines
    # Clustering near the start
    I = 100
    start = {'x':1, 'y': bd1.Points[-1].y}
    end = {'x': 21, 'y': 0}
    length = np.sqrt((end['x'] - start['x'])**2 + (end['y'] - start['y'])**2)
    ds1 = 0.001 * chord
    delta = find_hyp_tan_stretching_factor(ds1, I, length)
    stretching_f = lambda xi: 1 + (np.tanh(delta*(xi - 1)) / np.tanh(delta))
    r = fit_line(start, end, I, stretching_f)
    bd2 = Boundary(I, BdType.xi_max)
    bd2.populate_boundary_ndar_points(r)
    #plt = bd2.plot_boundary(plt, 'Boundary 2')    
    
    # Clustering near the end
    I = 100
    start = {'x': 0, 'y': 0}
    end = {'x':-20, 'y': 0}
    length = np.sqrt((end['x'] - start['x'])**2 + (end['y'] - start['y'])**2)
    ds1 = 0.001 * chord
    delta = find_hyp_tan_stretching_factor(ds1, I, length)
    stretching_f = lambda xi: 1 + (np.tanh(delta*(xi - 1)) / np.tanh(delta))
    r = fit_line(start, end, I, stretching_f)
    bd3 = Boundary(I, BdType.xi_min)
    bd3.populate_boundary_ndar_points(r)
    #plt = bd3.plot_boundary(plt, 'Boundary 3')        

    
    ################ Plotting the semicircle
    I = 100
    Radius = 20.5
    t1 = np.pi
    t2 = 0
    centre = {'x': 0.5, 'y':0}
    r = fit_semicircle(Radius, t1, t2, centre, I)
    bd4 = Boundary(I, BdType.eta_max)
    bd4.populate_boundary_ndar_points(r)
    #plt = bd4.plot_boundary(plt, 'Boundary 4')        
    
    #plt.legend(loc='best')
    #plt.show()

    return (bd1, bd2, bd3, bd4)

