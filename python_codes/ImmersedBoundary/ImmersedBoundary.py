from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import GetIBLines as gb

# This module serves as the platform to test any and all modules
# related to immersed boundary method. Mostly, just geometric
# modules can be tested since the solver is not integrated via
# python (yet ;-) )


def read_grid_2D(gridfile):
    f1 = open(gridfile, 'r')
    [imx, jmx] = [int(i) for i in f1.readline().split()]
    gridx = np.zeros((imx,jmx))
    gridy = np.zeros((imx,jmx))
    for j in range(jmx):
        for i in range(imx):
            [gridx[i][j], gridy[i][j]] = [float(k) for k in f1.readline().split()]
    f1.close()
    return [imx, jmx, gridx, gridy]


def convert_to_dict(start, end, normals):
    lines = []
    for i in range(len(start)):
        line = {'st': [start[i][0], start[i][1]], \
                'en': [end[i][0], end[i][1]], \
                'no': [normals[i][0], normals[i][1]]}
        lines.append(line)
    return lines


def first_cut_elimination():
    # To make life easier, make iblines as global. This is just a testing
    # module anyway
    imx, jmx = gridx.shape
    indices = []
    num_iblines = len(iblines)

    # IB need not be a closed loop. Hence include start and end
    centroid = [0, 0]
    for line in iblines:
        centroid[0] += line['st'][0] + line['en'][0]
        centroid[1] += line['st'][1] + line['en'][1]
    centroid[0] = centroid[0] / (num_iblines * 2)
    centroid[1] = centroid[1] / (num_iblines * 2)

    # Find max radius
    max_radius = 0.0
    for line in iblines:
        radius = np.sqrt((centroid[0] - line['st'][0])**2 + \
                         (centroid[1] - line['st'][1])**2)
        if radius > max_radius:
            max_radius = radius

    # Check if cell centre of the grid cells lie within the circle
    # of radius = max_radius

    for j in range(jmx-1):
        for i in range(imx-1):
            cell_centre = [0.25 * (\
                          gridx[i][j] + gridx[i+1][j] + \
                          gridx[i+1][j+1] + gridx[i][j+1]), \
                          0.25 * (\
                          gridy[i][j] + gridy[i+1][j] + \
                          gridy[i+1][j+1] + gridy[i][j+1])]
            dist = np.sqrt((centroid[0] - cell_centre[0])**2 + \
                           (centroid[1] - cell_centre[1])**2)
            
            #TODO: 20% allowance... 
            if dist < max_radius*1.2:
                indices.append([i, j])

#   # Plotting to check
#   # Plot IB lines as blue lines
#   for line in iblines:
#       x = [line['st'][0], line['en'][0]]
#       y = [line['st'][1], line['en'][1]]
#       plt.plot(x, y, 'b')
#   # Plot a circle with the centroid as green circle
#   t = np.linspace(0, 2*np.pi, 101)
#   x = centroid[0] + max_radius*np.cos(t)
#   y = centroid[1] + max_radius*np.sin(t)
#   plt.plot(x, y, 'g', label='Outer circle of IB')
#   plt.plot(gridx.flatten(), gridy.flatten(), 'k.', markersize=2)
#   # Plot the cell centres within the circle as red stars
#   # Plot grid points as black dots
#   for ind in indices:
#       i, j = ind[0], ind[1]
#       plt.plot(gridx[i][j], gridy[i][j], 'k.', markersize=4)
#       plt.plot(gridx[i+1][j], gridy[i+1][j], 'k.', markersize=4)
#       plt.plot(gridx[i+1][j+1], gridy[i+1][j+1], 'k.', markersize=4)
#       plt.plot(gridx[i][j+1], gridy[i][j+1], 'k.', markersize=4)
#       cell_centre = [0.25 * (\
#                     gridx[i][j] + gridx[i+1][j] + \
#                     gridx[i+1][j+1] + gridx[i][j+1]), \
#                     0.25 * (\
#                     gridy[i][j] + gridy[i+1][j] + \
#                     gridy[i+1][j+1] + gridy[i][j+1])]

#       plt.plot(cell_centre[0], cell_centre[1], 'r*')
#   plt.legend(loc='best')
#   plt.axis('equal')
#   plt.show()

    return indices


def find_intersection(x0, y0, ibline):
    # ibline is a dictionary of a line with st, en and no as keys
    # Essentially solving a set of 2 linear equations
    s_x, s_y = ibline['st'][0], ibline['st'][1]
    e_x, e_y = ibline['en'][0], ibline['en'][1]
    n_x, n_y = ibline['no'][0], ibline['no'][1]
    a11 = s_x - e_x
    a12 = n_x
    a21 = s_y - e_y
    a22 = n_y
    b1 = s_x - x0
    b2 = s_y - y0

    det = a11*a22 - a12*a21
    # Inverse of [[a b]; [c d]] = 1/(ad - bc) * [[d -b]; [-c a]]
    alpha = (a22*b1 - a12*b2) / det
    beta = (-a21*b1 + a11*b2) / det
    return (alpha, beta)


def find_nearest_ibline(x, y):  
    # Distance between point and line
    # Using equation 14 from mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    pt_line = {}
    pt_line['min_d'] = 1e10

    for i in range(len(iblines)):
        temp_dict = {}
        ibline = iblines[i]
        s_x, s_y = ibline['st'][0], ibline['st'][1]
        e_x, e_y = ibline['en'][0], ibline['en'][1]
        alpha, beta = find_intersection(x, y, ibline)
        xnc = x + beta*ibline['no'][0]
        ync = y + beta*ibline['no'][1]
        if 0 < alpha < 1:
            d = np.sqrt((x - xnc)**2 + (y - ync)**2)
            temp_dict['int_pt'] = [xnc, ync]
        else:
            ds = np.sqrt((x - s_x)**2 + (y - s_y)**2)
            de = np.sqrt((x - e_x)**2 + (y - e_y)**2)
            if ds < de:
                d = ds
                temp_dict['int_pt'] = [s_x, s_y]
            else:
                d = de
                temp_dict['int_pt'] = [e_x, e_y]
        if d <= pt_line['min_d']:
            pt_line['min_d'] = d
            pt_line['min_d_ind'] = i
            pt_line['alpha'] = alpha
            pt_line['beta'] = beta
            pt_line['int_pt'] = temp_dict['int_pt']
    
    return pt_line


#ef classify_inside_outside(first_cut_indices):
#   # Indices are the indices of the 2D array which are got from first cut elimination
#   # Idea is that for each cell centre in the first elimination, the closest immersed
#   # boundary/line is found out. Then it is found if the cell centre lies internal or
#   # external to the line

#   # Using cell centres to classify inside or outside
#   inside_cell_indices = []
#   outside_cell_indices = []

#   for ind in first_cut_indices:
#       i, j = ind[0], ind[1]
#       # Using cell centres. first_cut_indices are indices of cells
#       x, y = [0.25 * (\
#              gridx[i][j] + gridx[i+1][j] + \
#              gridx[i+1][j+1] + gridx[i][j+1]), \
#              0.25 * (\
#              gridy[i][j] + gridy[i+1][j] + \
#              gridy[i+1][j+1] + gridy[i][j+1])]
#       cell_dict = find_nearest_ibline(x, y)
#       cell_dict['cell_ind'] = [i, j]
#       
#   #   # Construct line to debug
#   #   xc, yc = cell_dict['int_pt'][0], cell_dict['int_pt'][1]
#   #   xp, yp = x, y
#   #   beta = cell_dict['beta']
#   #   alpha = cell_dict['alpha']
#   #   if 0 < alpha < 1:
#   #       pass
#   #   else:
#   #       if beta < 0:
#   #           continue

#   #       ibline = iblines[cell_dict['min_d_ind']]
#   #     # xc = xp + beta*ibline['no'][0]
#   #     # yc = yp + beta*ibline['no'][1]
#   #       plt.plot([x, xc], [y, yc], 'g')
#   #       for line in iblines:
#   #           x = [line['st'][0], line['en'][0]]
#   #           y = [line['st'][1], line['en'][1]]
#   #           plt.plot(x, y, 'b')

#   #       s_x, s_y = ibline['st'][0], ibline['st'][1]
#   #       e_x, e_y = ibline['en'][0], ibline['en'][1]
#   #       plt.plot([s_x, e_x, xp], [s_y, e_y, yp], 'k*')
#   #       plt.plot([s_x, e_x], [s_y, e_y], 'k-', linewidth=4)
#   #       plt.axis('equal')
#   #       if beta > 0:
#   #           loc = 'internal'
#   #       else:
#   #           loc = 'external'
#   #       titlestr = 'alpha = ' + alpha.__str__() + ', beta = ' + beta.__str__() + '. Marked as: ' + loc
#   #       plt.title(titlestr)
#   #       plt.show()

#       if cell_dict['beta'] > 0:
#           # Cell centre is inside
#           # If beta = 0, then the cell centre is on the ibline, which we
#           # term to be external / outside
#           if not 0 < cell_dict['alpha'] < 1:
#               # It means that the intersection point has been chosen as one of the
#               # vertices of the line segment
#               # Possible that it is actually interior, checking next line segment
#               # which shares the same vertex
#               ib_ind = cell_dict['min_d_ind']
#               xnc, ync = cell_dict['int_pt']
#               s_x, s_y = iblines[ib_ind]['st']
#               e_x, e_y = iblines[ib_ind]['en']
#               if [xnc, ync] == iblines[ib_ind]['st']:
#                   # The end of the previous line is same
#                   # Border case: ib_ind = 0
#                   if ib_ind == 0:
#                       common_v_ind = len(iblines) - 1
#                   else:
#                       common_v_ind = ib_ind - 1
#               else:
#                   # The start of the next line is the same
#                   # Border case: ib_ind = len(iblines) - 1
#                   if ib_ind == len(iblines) - 1:
#                       commmon_v_ind = 0
#                   else:
#                       common_v_ind = ib_ind + 1
#               # Check status of neighbour
#               alpha_new, beta_new = find_intersection(x, y, iblines[common_v_ind])
#               if beta_new > 0:
#                   # Neighbour also internal, so confirmed
#                   inside_cell_indices.append(cell_dict)
#               else:
#                   outside_cell_indices.append(cell_dict)
#           else:
#               # 0 < alpha < 1. Not a border condition. No issue
#               inside_cell_indices.append(cell_dict)
#       else:
#           # Outside only...
#           outside_cell_indices.append(cell_dict)

#   return(inside_cell_indices, outside_cell_indices)


def classify_inside_outside(first_cut_indices):
    # Indices are the indices of the 2D array which are got from first cut elimination
    # Idea is that for each cell centre in the first elimination, the closest immersed
    # boundary/line is found out. Then it is found if the cell centre lies internal or
    # external to the line

    # Using cell centres to classify inside or outside
    inside_cell_indices = []
    outside_cell_indices = []
    num_common = 0
    print 'Length of first cut indices ', len(first_cut_indices)

    for ind in first_cut_indices:
        cell_dict = {}
        i, j = ind[0], ind[1]
        # Using cell centres. first_cut_indices are indices of cells
        x, y = [0.25 * (\
               gridx[i][j] + gridx[i+1][j] + \
               gridx[i+1][j+1] + gridx[i][j+1]), \
               0.25 * (\
               gridy[i][j] + gridy[i+1][j] + \
               gridy[i+1][j+1] + gridy[i][j+1])]
        cell_dict = find_nearest_ibline(x, y)
        cell_dict['cell_ind'] = [i, j]
        
        if ((not (0 < cell_dict['alpha'] < 1)) and (cell_dict['beta'] >= 0)):
            # It means that the intersection point has been chosen as one of the
            # vertices of the line segment
            # Possible that it is actually interior, checking next line segment
            # which shares the same vertex
            ib_ind = cell_dict['min_d_ind']
            xnc, ync = cell_dict['int_pt']
            s_x, s_y = iblines[ib_ind]['st']
            e_x, e_y = iblines[ib_ind]['en']

     #      if [xnc, ync] == iblines[ib_ind]['st']:
     #          # The end of the previous line is same
     #          # Border case: ib_ind = 0
     #          if ib_ind == 0:
     #              common_v_ind = len(iblines) - 1
     #          else:
     #              common_v_ind = ib_ind - 1
     #      elif [xnc, ync] == iblines[ib_ind]['en']:
     #          # The start of the next line is the same
     #          # Border case: ib_ind = len(iblines) - 1
     #          if ib_ind == len(iblines) - 1:
     #              common_v_ind = 0
     #          else:
     #              common_v_ind = ib_ind + 1
     #      else:
     #          print 'Idiot! something is wrong'

            common_v_ind = -1
            for h in range(len(iblines)):
                if [xnc, ync] == iblines[h]['st']:
                    if not h == ib_ind:
                        common_v_ind = h
                elif [xnc, ync] == iblines[h]['en']:
                    if not h == ib_ind:
                        common_v_ind = h
            
            if not common_v_ind == -1:
                num_common =  num_common + 1
                # Check status of neighbour
                alpha_new, beta_new = find_intersection(x, y, iblines[common_v_ind])
                cell_dict['alpha'] = alpha_new
                cell_dict['beta'] = beta_new
                cell_dict['min_d_ind'] = common_v_ind

        if cell_dict['beta'] >= 0:
            inside_cell_indices.append(cell_dict)
        else:
            outside_cell_indices.append(cell_dict)
    
    print 'Num common ', num_common
    return(inside_cell_indices, outside_cell_indices)


def classify_int_ext_band():
    # Now that we have classified each cell centre as inside or outside the domain,
    # we use 'neighbour' logic to determine if the cell is purely internal or external
    # to the immersed boundary

    # Logic: Among the exterior cells, if at least one neighbour is internal, then band 
    #        cell. Else, exterior. Interior cells are designated as interior
    # This is the part of the logic that needs to be tested

    # Create list of interior cell indices for comparison
    global band_cell_indices
    band_cell_indices = []
    inside_cell_ind_list = [ins['cell_ind'] for ins in inside_cell_indices]
    for out_cell in outside_cell_indices:
        i,j = out_cell['cell_ind']
        
        if (([i+1, j] in inside_cell_ind_list) or ([i-1, j] in inside_cell_ind_list) or \
            ([i, j+1] in inside_cell_ind_list) or ([i, j-1] in inside_cell_ind_list)):
            band_cell_indices.append(out_cell)

    # Remove band cells from outside_cells
    for band_cell in band_cell_indices:
        if band_cell in outside_cell_indices:
            outside_cell_indices.remove(band_cell)


def classify_int_ext_band_new():
    # Now that we have classified each cell centre as inside or outside the domain,
    # we use 'neighbour' logic to determine if the cell is purely internal or external
    # to the immersed boundary

    # Logic: Among the exterior cells, if at least one neighbour is internal, then band 
    #        cell. Else, exterior. Interior cells are designated as interior
    # This is the part of the logic that needs to be tested

    # Create list of interior cell indices for comparison
    global band_cell_indices
    band_cell_indices = []
    outside_cell_ind_list = [ins['cell_ind'] for ins in outside_cell_indices]
    for in_cell in inside_cell_indices:
        i,j = in_cell['cell_ind']
        
        if (([i+1, j] in outside_cell_ind_list) or ([i-1, j] in outside_cell_ind_list) or \
            ([i, j+1] in outside_cell_ind_list) or ([i, j-1] in outside_cell_ind_list)):
            band_cell_indices.append(in_cell)

    # Remove band cells from outside_cells
    for band_cell in band_cell_indices:
        if band_cell in inside_cell_indices:
            inside_cell_indices.remove(band_cell)


def plot_and_check():
    # Plot the entire grid
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(gridx.flatten(), gridy.flatten(), 'k.', markersize=2)

    # Plot the immersed line
    # Plot IB lines as blue lines
    for line in iblines:
        x = [line['st'][0], line['en'][0]]
        y = [line['st'][1], line['en'][1]]
        plt.plot(x, y, 'r')

    # Plotting rectanges: This seems close to impossible as control
    # over colour
    # Plot interior cells as red plus
    # Plot interior cells as red traingles
    for int_cell in inside_cell_indices:
        i, j = int_cell['cell_ind']
        rect = np.array([\
                [gridx[i][j], gridy[i][j]], \
                [gridx[i+1][j], gridy[i+1][j]], \
                [gridx[i+1][j+1], gridy[i+1][j+1]], \
                [gridx[i][j+1], gridy[i][j+1]]])
        ax.add_patch(Polygon(rect, True, color='red', alpha=0.4))
        x, y = [0.25 * (\
                gridx[i][j] + gridx[i+1][j] + \
                gridx[i+1][j+1] + gridx[i][j+1]), \
                0.25 * (\
                gridy[i][j] + gridy[i+1][j] + \
                gridy[i+1][j+1] + gridy[i][j+1])]
        plt.plot(x, y, 'k*', markersize=3)

    # Plot exterior cells as yellow dots
#   for ext_cell in outside_cell_indices:
#       i, j = ext_cell['cell_ind']
#       rect = np.array([\
#               [gridx[i][j], gridy[i][j]], \
#               [gridx[i+1][j], gridy[i+1][j]], \
#               [gridx[i+1][j+1], gridy[i+1][j+1]], \
#               [gridx[i][j+1], gridy[i][j+1]]])
#       ax.add_patch(Polygon(rect, True, color='yellow', alpha=0.4))
       #x, y = [0.25 * (\
       #        gridx[i][j] + gridx[i+1][j] + \
       #        gridx[i+1][j+1] + gridx[i][j+1]), \
       #        0.25 * (\
       #        gridy[i][j] + gridy[i+1][j] + \
       #        gridy[i+1][j+1] + gridy[i][j+1])]
       #plt.plot(x, y, 'y.', markersize=10)

    # Plot band cells as blue stars
    for band_cell in band_cell_indices:
        i, j = band_cell['cell_ind']
        rect = np.array([\
                [gridx[i][j], gridy[i][j]], \
                [gridx[i+1][j], gridy[i+1][j]], \
                [gridx[i+1][j+1], gridy[i+1][j+1]], \
                [gridx[i][j+1], gridy[i][j+1]]])
        ax.add_patch(Polygon(rect, True, color='blue', alpha=0.4))
      
        x, y = [0.25 * (\
                gridx[i][j] + gridx[i+1][j] + \
                gridx[i+1][j+1] + gridx[i][j+1]), \
                0.25 * (\
                gridy[i][j] + gridy[i+1][j] + \
                gridy[i+1][j+1] + gridy[i][j+1])]
        plt.plot(x, y, 'k^', markersize=3)

    plt.axis('equal')
    plt.show()


def plot_from_file():
    # Plot the entire grid
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(gridx.flatten(), gridy.flatten(), 'k.', markersize=2)

    # Plot the immersed line
    # Plot IB lines as blue lines
    for line in iblines:
        x = [line['st'][0], line['en'][0]]
        y = [line['st'][1], line['en'][1]]
        plt.plot(x, y, 'r', linewidth=2)

    band_indices = []
    with open('band_faces', 'r') as f:
        for line in f:
            s = line.split()
            s[1] = int(s[1])
            s[2] = int(s[2])
            s[3] = int(s[3])
            s[4] = int(s[4])
            s[5] = int(s[5])
            s[6] = int(s[6])
            band_indices.append(s)
    for ind in band_indices:
        fdir, i, j, i_f, j_f, i_b, j_b = ind
        i = i-1
        j = j-1
        i_f = i_f - 1
        j_f = j_f - 1
        i_b = i_b - 1
        j_b = j_b - 1
        i_face = i + i_f - i_b
        j_face = j + j_f - j_b
        if fdir == 'x':
            x = [gridx[i, j], gridx[i, j+1]]
            y = [gridy[i, j], gridy[i, j+1]]
            x2 = [gridx[i_face, j_face], gridx[i_face, j_face+1]]
            y2 = [gridy[i_face, j_face], gridy[i_face, j_face+1]]
        else:
            x = [gridx[i, j], gridx[i+1, j]]
            y = [gridy[i, j], gridy[i+1, j]]
            x2 = [gridx[i_face, j_face], gridx[i_face+1, j_face]]
            y2 = [gridy[i_face, j_face], gridy[i_face+1, j_face]]
        plt.plot(x, y, 'k-', linewidth=2)
#       plt.plot(x2, y2, 'k-', linewidth=4)
            
    bandcells_ind = np.loadtxt('band_cells')
    for ind in bandcells_ind:
        i, j = ind
        i = i-1
        j = j-1
        rect = np.array([\
                [gridx[i][j], gridy[i][j]], \
                [gridx[i+1][j], gridy[i+1][j]], \
                [gridx[i+1][j+1], gridy[i+1][j+1]], \
                [gridx[i][j+1], gridy[i][j+1]]])
        ax.add_patch(Polygon(rect, True, color='blue', alpha=0.4))
        x, y = [0.25 * (\
                gridx[i][j] + gridx[i+1][j] + \
                gridx[i+1][j+1] + gridx[i][j+1]), \
                0.25 * (\
                gridy[i][j] + gridy[i+1][j] + \
                gridy[i+1][j+1] + gridy[i][j+1])]
        plt.plot(x, y, 'k^', markersize=3)
    
    interiorcells_ind = np.loadtxt('interior_cells')
    for ind in interiorcells_ind:
        i, j = ind
        i = i-1
        j = j-1
        rect = np.array([\
                [gridx[i][j], gridy[i][j]], \
                [gridx[i+1][j], gridy[i+1][j]], \
                [gridx[i+1][j+1], gridy[i+1][j+1]], \
                [gridx[i][j+1], gridy[i][j+1]]])
        ax.add_patch(Polygon(rect, True, color='red', alpha=0.4))
        x, y = [0.25 * (\
                gridx[i][j] + gridx[i+1][j] + \
                gridx[i+1][j+1] + gridx[i][j+1]), \
                0.25 * (\
                gridy[i][j] + gridy[i+1][j] + \
                gridy[i+1][j+1] + gridy[i][j+1])]
        plt.plot(x, y, 'k*', markersize=3)

    plt.xlim([-0.05, 1.05])
    plt.axis('equal')
    plt.show()



def write_out():
    f = open('band_cells_py', 'w')
    for cell in band_cell_indices:
        i, j = cell['cell_ind']
        f.write(i.__str__() + '    ' + j.__str__() + '\n')
    f.close()

    f = open('interior_cells_py', 'w')
    for cell in inside_cell_indices:
        i, j = cell['cell_ind']
        f.write(i.__str__() + '    ' + j.__str__() + '\n')
    f.close()


if __name__ == '__main__':
#   gridfile = 'Circle-IB-3.txt'
#   gridfile = 'Plane-IB-5.txt'
    gridfile = 'ExpFan-IB.txt'
   #gridfile = 'Rchannel-IB.txt'
#   gridfile = 'Circle-IB.txt'
#   IBfile = 'naca0012-3.dat'
#   IBfile = 'circle.dat'
#   IBfile = 'ramp2.dat'
 #  IBfile = 'naca0012-2.dat'
#   IBfile = 'concave.dat'
    IBfile = 'exp_fan.dat'
  # IBfile = 'ramp.dat'
    global gridx
    global gridy
    global iblines
    global inside_cell_indices
    global outside_cell_indices
    global band_cell_indices
    
    imx, jmx, gridx, gridy = read_grid_2D(gridfile)
    start, end, normals = gb.read_from_file(IBfile)
    iblines = convert_to_dict(start, end, normals)
    
#   first_cut_indices = first_cut_elimination()
#   inside_cell_indices, outside_cell_indices = classify_inside_outside(first_cut_indices)
#   classify_int_ext_band_new()
#   plot_and_check()

    plot_from_file()
#   write_out()
