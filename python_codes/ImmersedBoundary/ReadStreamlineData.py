from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def extract_data_along_streamline(filename):
    f = open(filename, 'r')

    # Ignore random lines
    for i in range(4):
        f.readline()

    # next line contains number of points
    l = f.readline()
    l = l.split()
    num_points = int(l[1])

    # From the next line that many points are to be read
    points_read = 0
    points = np.empty((num_points, 3))

    while(points_read < num_points):
        l = f.readline()
        l = l.split()
        n = int(len(l) / 3)
        for i in range(n):
            points[points_read, 0] = float(l[i*3 + 0])
            points[points_read, 1] = float(l[i*3 + 1])
            points[points_read, 2] = float(l[i*3 + 2])
            points_read += 1

    # Read line connectivity?? Right now skip
    while not (f.readline().startswith('LINES')):
        pass
    f.readline()

    # blank line
    f.readline()

    # Point data. We know it is same as num_points
    pressure = np.empty(num_points)
    while not (f.readline().startswith('POINT_DATA')):
        pass
    for i in range(2):
        f.readline()

    pressure_read = 0
    while (pressure_read < num_points):
        l = f.readline()
        l = l.split()
        n = len(l)
        for i in range(n):
            pressure[pressure_read] = float(l[i])
            pressure_read += 1

    print points_read, pressure_read, num_points

    return (pressure, points)


def plot_cp_vs_x(pressure, points, p_inf, v_inf, rho_inf):
    x = points[:, 0]
    Cp = (pressure - p_inf) / (0.5 * rho_inf * v_inf * v_inf)

    x_new = []
    Cp_new = []
    for i in range(len(x)):
        if ((x[i] >= 0) and (x[i] <= 1)):
            x_new.append(x[i])
            Cp_new.append(Cp[i])
    
    plt.plot(x_new, Cp_new, 'b.-', label='Cp')
    plt.title('NACA 0012 Cp vs x')
    plt.xlabel('x (m)')
    plt.ylabel('Cp')
    plt.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    filename = 'bottom.vtk'
    p_inf = 90748.127
    rho_inf = 1.1322
    u_inf = 133.99

    pressure, points = extract_data_along_streamline(filename)
    plot_cp_vs_x(pressure, points, p_inf, u_inf, rho_inf)


