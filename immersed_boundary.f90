module immersed_boundary

    use global, only: FILE_NAME_LENGTH, IB_FILE_UNIT, &
                DESCRIPTION_STRING_LENGTH
    use grid, only: grid_x, grid_y, imx, jmx
    use geometry, only: xnx, xny, ynx, yny, xA, yA, volume
    use utils, only: alloc, dealloc, dmsg
    use state, only: density, x_speed, y_speed, pressure, mu_ref, gm, &
                     density_inf, pressure_inf, x_speed_inf, y_speed_inf, R_gas
    use face_interpolant
    use scheme, only: residue, F_p, G_p
    use viscous, only: T_ref, Sutherland_Temp, Pr

    implicit none
    private

    type line
        real, dimension(2) :: st, en, no
        real :: u, v
    end type

    type domain
        integer, dimension(2):: cell_ind
        real, dimension(2) :: int_point
        real :: alpha, beta, min_d
        integer :: nearest_ibline_ind
    end type

    type face
        integer, dimension(2) :: face_ind, field_ind, band_ind
        character :: face_type
        integer :: nearest_ibline_ind
        real :: band_min_d, field_min_d, face_min_d
        real, dimension(2) :: int_point
    end type

    integer, dimension(:, :), allocatable, public :: cell_type
    type(line), dimension(:), allocatable, public :: iblines
    type(domain), dimension(:), allocatable, public :: bandcells, interiorcells, &
                                            tempcells
    type(face), dimension(:), allocatable, public :: bandfaces
    integer, dimension(:, :), allocatable :: first_cut_indices
    character(len=FILE_NAME_LENGTH), public :: IBfilename
    integer, public :: ind_count, num_iblines, num_bandcells, &
                       num_interiorcells, num_tempcells
    integer, public :: num_bandfaces
    real, public, dimension(:, :), allocatable :: signed_distance
    character(len=DESCRIPTION_STRING_LENGTH), public :: signed_dist_varname

    public :: setup_IB, IB_step, reset_states_at_interface_faces, &
              destroy_IB, update_band_interior_cells, reset_gradients_at_interfaces

    contains

    subroutine read_iblines()
        ! This module reads the immersed boundary coordinates from a 
        ! given file

        implicit none

        integer :: i, mem_stat
        
        call dmsg(1, 'immersed_boundary', 'read_iblines')

        open(IB_FILE_UNIT, file=IBfilename)
        read(IB_FILE_UNIT, *) num_iblines
        
        allocate(iblines(1:num_iblines), stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not allocate memory for imlines'
            stop
        end if

        do i = 1, num_iblines
            read(IB_FILE_UNIT, *) iblines(i)%st(1), iblines(i)%st(2), &
                        iblines(i)%en(1), iblines(i)%en(2), &
                        iblines(i)%no(1), iblines(i)%no(2) 

            !TODO: Stationary IB. Make generic
            iblines(i)%u = 0.0
            iblines(i)%v = 0.0
        end do
        close(IB_FILE_UNIT)

    end subroutine read_iblines

    subroutine setup_IB()

        implicit none
        
        call dmsg(1, 'immersed_boundary', 'setup_IB')
        
        call read_iblines()
        ! Change alloc size if shifting to node based selection
        call alloc(cell_type, 1, imx-1, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for first_cut_blanket')
        call alloc(first_cut_indices, 1, (imx-1)*(jmx-1), 1, 2, &
                    errmsg='Error: Unable to allocate memory for first cut indices')
        call alloc(signed_distance, 1, imx-1, 1, jmx - 1, &
                    errmsg='Error: Unable to allocate memory for signed_distance')
        
        signed_distance = -1.0e8
        signed_dist_varname = 'Signed_distance'

        num_bandcells = 0
        num_interiorcells = 0
        num_tempcells = 0

    end subroutine setup_IB

    subroutine first_cut_elimination()
        ! -----------------------------------------------------------------------
        ! Eliminate a first round of cell centres to be classified by constructing
        ! a rectangular bounding box
        ! Have a blanket array cell_type
        ! cell_type(i, j) = 0, if i,j is an external cell
        ! cell_type(i, j) = -1, if i,j is an internal cell
        ! All first cut elimination cells are first marked as -1. 
        ! The next subroutine, which classifies them properly, resets
        ! the actual external cells as 0
        ! -----------------------------------------------------------------------
    
        implicit none

        integer :: i, j
        type(line) :: ibline
        real :: xmin, xmax, ymin, ymax, x, y, l, b
        real, dimension(2) :: cell_centre
        
        call dmsg(1, 'immersed_boundary', 'first_cut_elimination')

        xmin = 1e10
        xmax = -1e10
        ymin = 1e10
        ymax = -1e10
        ind_count = 0

        do i = 1, num_iblines
            ibline = iblines(i)
            x = ibline%st(1)
            y = ibline%st(2)
            if (x < xmin) then
                xmin = x
            end if
            if (x > xmax) then
                xmax = x
            end if
            if (y < ymin) then
                ymin = y
            end if
            if (y > ymax) then
                ymax = y
            end if

            x = ibline%en(1)
            y = ibline%en(2)
            if (x < xmin) then
                xmin = x
            end if
            if (x > xmax) then
                xmax = x
            end if
            if (y < ymin) then
                ymin = y
            end if
            if (y > ymax) then
                ymax = y
            end if
        end do

        ! Add cushioning of 10%
        l = xmax - xmin
        b = ymax - ymin
        xmin = xmin - 0.2*l
        xmax = xmax + 0.2*l
        ymin = ymin - 0.2*b
        ymax = ymax + 0.2*b

        ! Check if cell centre of the grid cells lie within the rectangle
        do j = 1, jmx - 1
         do i = 1, imx - 1
            cell_centre(1) = 0.25 * (grid_x(i, j) + grid_x(i+1, j) + &
                                     grid_x(i+1, j+1) + grid_x(i, j+1))
            cell_centre(2) = 0.25 * (grid_y(i, j) + grid_y(i+1, j) + &
                                     grid_y(i+1, j+1) + grid_y(i, j+1))
            
            if ((cell_centre(1) >= xmin) .and. (cell_centre(1) <= xmax) .and. &
                (cell_centre(2) >= ymin) .and. (cell_centre(2) <= ymax)) then
                ! Could be internal cell. Mark as -1
                cell_type(i, j) = -1
                ind_count = ind_count + 1
            else
                ! Sure as hell is an external cell. Mark as 0
                cell_type(i, j) = 0
            end if
         end do
        end do

    end subroutine first_cut_elimination

    function get_intersection(ibline_ind, xc, yc) result(res_matrix)
        ! -----------------------------------------------------------------------
        ! Given a ibline and a cell centre, find the intersection point, signed distance,
        ! alpha, beta
        ! -----------------------------------------------------------------------

        implicit none

        integer, intent(in) :: ibline_ind
        real, intent(in) :: xc, yc
        real, dimension(5) :: res_matrix
        
        type(line) :: ibline
        real :: a11, a12, a21, a22, b1, b2, det
        real :: xnc, ync, d, de, ds
        real :: alpha, beta
        real, dimension(2) :: int_point

        ibline = iblines(ibline_ind)
        a11 = ibline%st(1) - ibline%en(1)
        a12 = ibline%no(1)
        a21 = ibline%st(2) - ibline%en(2)
        a22 = ibline%no(2)
        b1 = ibline%st(1) - xc
        b2 = ibline%st(2) - yc

        det = a11*a22 - a12*a21
        ! Inverse of [[a b]; [c d]] = 1 / (ad - bc) * [[d -b]; [-c a]]
        alpha = (a22*b1 - a12*b2) / det
        beta = (-a21*b1 + a11*b2) / det

        xnc = xc + beta*ibline%no(1)
        ync = yc + beta*ibline%no(2)
        if ((alpha > 0) .and. (alpha < 1)) then
            d = sqrt((xc - xnc)**2 + (yc - ync)**2)
            int_point = (/ xnc, ync /)
        else
            ds = sqrt((xc - ibline%st(1))**2 + (yc - ibline%st(2))**2)
            de = sqrt((xc - ibline%en(1))**2 + (yc - ibline%en(2))**2)
            if (ds < de) then
                d = ds
                int_point = (/ ibline%st(1), ibline%st(2) /)
            else
                d = de
                int_point = (/ ibline%en(1), ibline%en(2) /)
            end if
        end if

        res_matrix(1) = alpha
        res_matrix(2) = beta
        res_matrix(3) = int_point(1)
        res_matrix(4) = int_point(2)
        res_matrix(5) = d
        
    end function get_intersection

    function in_or_out(i, j) result(cell)
        ! -----------------------------------------------------------------------
        ! Given a cell centre, classify it as internal or external
        ! The nearest IB line is first found, then the point is classified as internal
        ! or external to the IB. Band cell classification is not done here
        ! A list of internal cells are stored, but the external cells are not stored
        ! separately
        ! -----------------------------------------------------------------------

        implicit none

        type(domain) :: cell
        integer, intent(in) :: i, j
        type(line) :: ibline
        integer :: n, common_ind, l
        real :: xc, yc, d, alpha, beta
        real, dimension(5) :: res_matrix
        real, dimension(2) :: int_point
        real :: dot_prod

        ! Choosing centroid
        xc = 0.25 * (grid_x(i, j) + grid_x(i+1, j) + &
                     grid_x(i+1, j+1) + grid_x(i, j+1))
        yc = 0.25 * (grid_y(i, j) + grid_y(i+1, j) + &
                     grid_y(i+1, j+1) + grid_y(i, j+1))

        cell%min_d = 1.0e10
        cell%cell_ind = (/ i, j /)
        ! Finding nearest ib line to this centroid
        do l = 1, num_iblines
            ibline = iblines(l)
            res_matrix = get_intersection(l, xc, yc)
            alpha = res_matrix(1) 
            beta = res_matrix(2) 
            int_point(1) = res_matrix(3) 
            int_point(2) = res_matrix(4) 
            d = res_matrix(5) 
            if (abs(d) <= abs(cell%min_d)) then
                ! beta >= 0 -> interior or band
                ! beta < 0 -> exterior
                cell%min_d = d * beta / abs(beta)
                cell%nearest_ibline_ind = l
                cell%alpha = alpha
                cell%beta = beta
                cell%int_point = int_point
            end if
        end do

        ! Nearest ib line found
        int_point = cell%int_point
        beta = cell%beta
        if ( (cell%alpha < 0.) .or. (cell%alpha > 1.) ) then
            ! Could be a concave issue. Check neighbour 
            common_ind = -1
            do n = 1, num_iblines
                if ((int_point(1) .eq. iblines(n)%st(1)) .and. &
                    (int_point(2) .eq. iblines(n)%st(2))) then
                    if (n /= cell%nearest_ibline_ind) then
                        common_ind = n
                        exit
                    end if
                else if ((int_point(1) .eq. iblines(n)%en(1)) .and. &
                         (int_point(2) .eq. iblines(n)%en(2))) then
                    if (n /= cell%nearest_ibline_ind) then
                        common_ind = n
                        exit
                    end if
                end if
            end do
            if (common_ind /= -1) then
                dot_prod = ((xc - int_point(1))*(iblines(cell%nearest_ibline_ind)%no(1) + &
                                                iblines(common_ind)%no(1))) + &
                           ((yc - int_point(2))*(iblines(cell%nearest_ibline_ind)%no(2) + &
                                                iblines(common_ind)%no(2)))
                if (dot_prod >= 0.0) then
                    cell%beta = -1.0
                else
                    cell%beta = 1.0
                end if
                cell%min_d = abs(cell%min_d) * cell%beta / abs(cell%beta)
            end if 
        end if
        
    end function in_or_out

    subroutine classify_cells()
        ! -----------------------------------------------------------------------
        ! Classify all cells as internal or external cell
        ! Store the set of internal cells (includes band)
        ! Blanket array cell_type is reset for cells found to be external
        ! ind_count is the number of cell centres inside bounding box
        ! -----------------------------------------------------------------------
        
        implicit none

        integer :: mem_stat, i, j
        type(domain) :: cell
        
        call dmsg(1, 'immersed_boundary', 'classify_cells')
        
        allocate(tempcells(1:ind_count), stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not allocate memory for interiorcells'
            stop
        end if

        do j = 1, jmx-1
         do i = 1, imx-1
            if (cell_type(i, j) .eq. -1) then
                cell = in_or_out(i, j)
                if (cell%beta >= 0) then
                    num_tempcells = num_tempcells + 1
                    tempcells(num_tempcells) = cell
                else
                    ! Ah ha! Found another external cell. Mark as 0
                    cell_type(i, j) = 0
                end if
                signed_distance(i, j) = cell%min_d
            end if
         end do
        end do
        
    end subroutine classify_cells

    subroutine band_cells_as_int()
        ! -----------------------------------------------------------------------
        ! The subset of interior cells which have at least one exterior cell
        ! as its neighbour are termed as band cells.
        ! The interface between a band cell and a field cell is called a band face.
        ! A single band cell can have multiple band faces, hence all the neighbours
        ! are found one by one.
        ! tempcells included both internal and band cells. This subroutine extracts
        ! the actual interior cells separately.
        ! -----------------------------------------------------------------------

        implicit none

        integer :: i, j, k, mem_stat
        type(domain) :: cell, cell2
        type(face) :: bandface
        real :: xb, yb, xe, ye, xf, yf, beta
        real, dimension(5) :: res_matrix
        logical :: n1, n2, n3, n4, n_avail

        call dmsg(1, 'immersed_boundary', 'band_cells_as_int')
        
        allocate(bandcells(1:num_tempcells), stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not allocate memory for bandcells'
            stop
        end if
        allocate(interiorcells(1:num_tempcells), stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not allocate memory for interior cells'
            stop
        end if
        allocate(bandfaces(1:2*num_tempcells), stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not allocate memory for bandfaces'
            stop
        end if

        num_bandcells = 0
        num_interiorcells = 0
        num_bandfaces = 0

        ! Using cell_type to check if any neighbour of temp cell is an external or
        ! field cell. But cell_type has an index range. If cell to be classified at
        ! boundary, then correspondng common face should not be considered as band face
        ! and should be ignored. Using the logical flag n_avail to check if neighbouring
        ! cell is available inside the domain and is not a ghost cell
        do k = 1, num_tempcells
            cell = tempcells(k)
            i = cell%cell_ind(1)
            j = cell%cell_ind(2)
            n1 = .false.
            n2 = .false.
            n3 = .false.
            n4 = .false.
            
            ! Checking the 4 neighbours for field cells
            ! Right neighbour
            n_avail = .true.
            if (i .eq. imx-1) n_avail = .false.
            if (n_avail) then
            if (cell_type(i+1, j) .eq. 0) then
                n1 = .true.
                num_bandfaces = num_bandfaces + 1
                bandface%face_ind = (/ i+1, j /)
                bandface%field_ind = (/ i+1, j /)
                bandface%band_ind = (/ i, j /)
                bandface%face_type = 'x'
                bandface%band_min_d = cell%min_d
                cell2 = in_or_out(i+1, j)
                bandface%field_min_d = cell2%min_d
                ! Nearest ibline: it will be the nearest ibline of the band cell 
                ! or the field cell 
                xf = 0.5 * (grid_x(i+1, j) + grid_x(i+1, j+1))
                yf = 0.5 * (grid_y(i+1, j) + grid_y(i+1, j+1))
                ! Checking the band cell
                xb = cell%int_point(1)
                yb = cell%int_point(2)
                ! Checking the exterior cell. Above loop breaks and hence value of l can be
                ! used to access the appropriate exterior cell
                xe = cell2%int_point(1)
                ye = cell2%int_point(2)
                ! Find distance from the midpoint of face to the above intersection points.
                ! Choose the ibline which is closer
                if (sqrt((xf - xb)**2 + (yf - yb)**2) < sqrt((xf - xe)**2 + (yf - ye)**2)) then
                    ! Choose that of band cell
                    bandface%nearest_ibline_ind = cell%nearest_ibline_ind
                else
                    bandface%nearest_ibline_ind = cell2%nearest_ibline_ind
                end if
                res_matrix = get_intersection(bandface%nearest_ibline_ind, xf, yf)
                beta = res_matrix(2) 
                bandface%face_min_d = res_matrix(5) * beta / abs(beta)
                bandface%int_point = (/ res_matrix(3), res_matrix(4) /)
                ! Store to list of bandfaces
                bandfaces(num_bandfaces) = bandface
            end if
            end if
            
            ! Left neighbour    
            n_avail = .true.
            if (i .eq. 1) n_avail = .false.
            if (n_avail) then
            if (cell_type(i-1, j) .eq. 0) then
                n2 = .true.
                ! Adding the interface between band and field cell
                num_bandfaces = num_bandfaces + 1
                bandface%face_ind = (/ i, j /)
                bandface%field_ind = (/ i-1, j /)
                bandface%band_ind = (/ i, j /)
                bandface%face_type = 'x'
                bandface%band_min_d = cell%min_d
                cell2 = in_or_out(i-1, j)
                bandface%field_min_d = cell2%min_d
                ! Nearest ibline: it will be the nearest ibline of the band cell 
                ! or the field cell
                xf = 0.5 * (grid_x(i, j) + grid_x(i, j+1))
                yf = 0.5 * (grid_y(i, j) + grid_y(i, j+1))
                ! Checking the band cell
                xb = cell%int_point(1)
                yb = cell%int_point(2)
                ! Checking the exterior cell. Above loop breaks and hence value of l can be
                ! used to access the appropriate exterior cell
                xe = cell2%int_point(1)
                ye = cell2%int_point(2)
                ! Find distance from the midpoint of face to the above intersection points.
                ! Choose the ibline which is closer
                if (sqrt((xf - xb)**2 + (yf - yb)**2) < sqrt((xf - xe)**2 + (yf - ye)**2)) then
                    ! Choose that of band cell
                    bandface%nearest_ibline_ind = cell%nearest_ibline_ind
                else
                    bandface%nearest_ibline_ind = cell2%nearest_ibline_ind
                end if
                res_matrix = get_intersection(bandface%nearest_ibline_ind, xf, yf)
                beta = res_matrix(2) 
                bandface%face_min_d = res_matrix(5) * beta / abs(beta)
                bandface%int_point = (/ res_matrix(3), res_matrix(4) /)
                ! Store to list of bandfaces
                bandfaces(num_bandfaces) = bandface
            end if
            end if

            ! Top neighbour
            n_avail = .true.
            if (j .eq. jmx-1) n_avail = .false.
            if (n_avail) then
            if (cell_type(i, j+1) .eq. 0) then
                n3 = .true.
                ! Adding the interface between band and field cell
                num_bandfaces = num_bandfaces + 1
                bandface%face_ind = (/ i, j+1 /)
                bandface%field_ind = (/ i, j+1 /)
                bandface%band_ind = (/ i, j /)
                bandface%face_type = 'y'
                bandface%band_min_d = cell%min_d
                cell2 = in_or_out(i, j+1)
                bandface%field_min_d = cell2%min_d
                ! Nearest ibline: it will be the nearest ibline of the band cell 
                ! or the field cell 
                xf = 0.5 * (grid_x(i, j+1) + grid_x(i+1, j+1))
                yf = 0.5 * (grid_y(i, j+1) + grid_y(i+1, j+1))
                ! Checking the band cell
                xb = cell%int_point(1)
                yb = cell%int_point(2)
                ! Checking the exterior cell. Above loop breaks and hence value of l can be
                ! used to access the appropriate exterior cell
                xe = cell2%int_point(1)
                ye = cell2%int_point(2)
                ! Find distance from the midpoint of face to the above intersection points.
                ! Choose the ibline which is closer
                if (sqrt((xf - xb)**2 + (yf - yb)**2) < sqrt((xf - xe)**2 + (yf - ye)**2)) then
                    ! Choose that of band cell
                    bandface%nearest_ibline_ind = cell%nearest_ibline_ind
                else
                    bandface%nearest_ibline_ind = cell2%nearest_ibline_ind
                end if
                res_matrix = get_intersection(bandface%nearest_ibline_ind, xf, yf)
                beta = res_matrix(2) 
                bandface%face_min_d = res_matrix(5) * beta / abs(beta)
                bandface%int_point = (/ res_matrix(3), res_matrix(4) /)
                ! Store to list of bandfaces
                bandfaces(num_bandfaces) = bandface
            end if
            end if

            ! Bottom neighbour
            n_avail = .true.
            if (j .eq. 1) n_avail = .false.
            if (n_avail) then
            if (cell_type(i, j-1) .eq. 0) then
                n4 = .true.
                ! Adding the interface between band and field cell
                num_bandfaces = num_bandfaces + 1
                bandface%face_ind = (/ i, j /)
                bandface%field_ind = (/ i, j-1 /)
                bandface%band_ind = (/ i, j /)
                bandface%face_type = 'y'
                bandface%band_min_d = cell%min_d
                cell2 = in_or_out(i, j-1)
                bandface%field_min_d = cell2%min_d
                ! Nearest ibline: it will be the nearest ibline of the band cell 
                ! or the field cell 
                xf = 0.5 * (grid_x(i, j) + grid_x(i+1, j))
                yf = 0.5 * (grid_y(i, j) + grid_y(i+1, j))
                ! Checking the band cell
                xb = cell%int_point(1)
                yb = cell%int_point(2)
                ! Checking the exterior cell. Above loop breaks and hence value of l can be
                ! used to access the appropriate exterior cell
                xe = cell2%int_point(1)
                ye = cell2%int_point(2)
                ! Find distance from the midpoint of face to the above intersection points.
                ! Choose the ibline which is closer
                if (sqrt((xf - xb)**2 + (yf - yb)**2) < sqrt((xf - xe)**2 + (yf - ye)**2)) then
                    ! Choose that of band cell
                    bandface%nearest_ibline_ind = cell%nearest_ibline_ind
                else
                    bandface%nearest_ibline_ind = cell2%nearest_ibline_ind
                end if
                res_matrix = get_intersection(bandface%nearest_ibline_ind, xf, yf)
                beta = res_matrix(2) 
                bandface%face_min_d = res_matrix(5) * beta / abs(beta)
                bandface%int_point = (/ res_matrix(3), res_matrix(4) /)
                ! Store to list of bandfaces
                bandfaces(num_bandfaces) = bandface
            end if  
            end if

            ! Multiple neighbours of the temp cell can be an exterior cell, but the
            ! cell should be counted as band cell only once. Hence if at least one of the
            ! neighbour is an exterior cell, then count it as band cell
            if (n1 .or. n2 .or. n3 .or. n4) then
                num_bandcells = num_bandcells + 1
                bandcells(num_bandcells) = cell
            end if

            ! If no neighbour of temp cells is an exterior cell, then it is
            ! interior cell
            if (.not. (n1 .or. n2 .or. n3 .or. n4)) then
                num_interiorcells = num_interiorcells + 1
                interiorcells(num_interiorcells) = cell
            end if
        end do

    end subroutine band_cells_as_int

    subroutine IB_step()

        implicit none

        call dmsg(1, 'immersed_boundary', 'ib_step')
        
        call first_cut_elimination()
        call classify_cells()
        call band_cells_as_int()

!       call dmsg(5, 'immersed_boundary', 'ib_step', 'IB classification complete')
!       call test_bandfaces()
!       call test_classification()

    end subroutine IB_step

    subroutine set_slip_velocity_face()
        ! -----------------------------------------------------------------------
        ! Let F1, F2 and B be 3 consecutive cells, where F stands for
        ! field cells and B for band cells. In addition to setting the
        ! state at the interface between F2 and B, the interface between
        ! F1 and F2 also needs to be reset since higher order methods
        ! use the cell centered value of B at this interface

        ! In structured grid if i is the index of a cell, then i is the 
        ! index of the left face of that cell

        ! From a simple schematic, it is seen that the index of the
        ! interface between F1 and F2 is (i + (i_f - i_b))
        ! Index of F1 is (i_f + (i_f - i_b))
        ! -----------------------------------------------------------------------
        
        implicit none

        integer :: i, j, i_f, j_f, i_b, j_b, n
        character :: f_type
        real :: u_temp, v_temp, nx, ny, df, dface
        type(face) :: bandface
        
        call dmsg(1, 'immersed_boundary', 'set_slip_velocity')

        do n = 1, num_bandfaces
            bandface = bandfaces(n)
            i = bandface%face_ind(1)
            j = bandface%face_ind(2)
            i_f = bandface%field_ind(1)
            j_f = bandface%field_ind(2)
            i_b = bandface%band_ind(1)
            j_b = bandface%band_ind(2)
            f_type = bandface%face_type
            df = bandface%field_min_d
            df = df + sign(1.0, df)*1e-15
            dface = bandface%face_min_d
            dface = min(0.0, dface)
            nx = iblines(bandface%nearest_ibline_ind)%no(1)
            ny = iblines(bandface%nearest_ibline_ind)%no(2)
            u_temp = x_speed(i_f, j_f)
            v_temp = y_speed(i_f, j_f)

            if (f_type .eq. 'x') then
                x_x_speed_left(i, j) = u_temp - ((1 - dface/df) * &
                                      (u_temp * nx + v_temp * ny) * nx)
                x_y_speed_left(i, j) = v_temp - ((1 - dface/df) * &
                                      (u_temp * nx + v_temp * ny) * ny)
                x_x_speed_right(i, j) = x_x_speed_left(i, j)                                      
                x_y_speed_right(i, j) = x_y_speed_left(i, j)                                      
            else
                y_x_speed_left(i, j) = u_temp - ((1 - dface/df) * &
                                      (u_temp * nx + v_temp * ny) * nx)
                y_y_speed_left(i, j) = v_temp - ((1 - dface/df) * &
                                      (u_temp * nx + v_temp * ny) * ny)
                y_x_speed_right(i, j) = y_x_speed_left(i, j)                                      
                y_y_speed_right(i, j) = y_y_speed_left(i, j)                                      
            end if
        end do

    end subroutine set_slip_velocity_face

    subroutine set_no_slip_velocity_face()
        ! -----------------------------------------------------------------------
        ! Let F1, F2 and B be 3 consecutive cells, where F stands for
        ! field cells and B for band cells. In addition to setting the
        ! state at the interface between F2 and B, the interface between
        ! F1 and F2 also needs to be reset since higher order methods
        ! use the cell centered value of B at this interface

        ! In structured grid if i is the index of a cell, then i is the 
        ! index of the left face of that cell

        ! From a simple schematic, it is seen that the index of the
        ! interface between F1 and F2 is (i + (i_f - i_b))
        ! Index of F1 is (i_f + (i_f - i_b))
        ! -----------------------------------------------------------------------
        
        implicit none

        character :: f_type
        integer :: n, i, j, i_f, j_f, i_b, j_b
        type(face) :: bandface
        type(line) :: ibline
        real :: d_f, d_ib, df, dface
        real :: x0, y0, x1, y1
        real :: u_f, v_f, u_ib, v_ib, nx, ny, u_temp, v_temp
        
        call dmsg(1, 'immersed_boundary', 'set_no_slip_velocity')

        do n = 1, num_bandfaces
            bandface = bandfaces(n)
            ibline = iblines(bandface%nearest_ibline_ind)
            i = bandface%face_ind(1)
            j = bandface%face_ind(2)
            i_f = bandface%field_ind(1)
            j_f = bandface%field_ind(2)
            i_b = bandface%band_ind(1)
            j_b = bandface%band_ind(2)
            f_type = bandface%face_type
            d_ib = abs(bandface%face_min_d)
            d_ib = max(1e-6, d_ib)

            if (f_type .eq. 'x') then
                x0 = 0.5 * (grid_x(i, j) + grid_x(i, j+1))
                y0 = 0.5 * (grid_y(i, j) + grid_y(i, j+1))
            else
                x0 = 0.5 * (grid_x(i, j) + grid_x(i+1, j))
                y0 = 0.5 * (grid_y(i, j) + grid_y(i+1, j))
            end if
            
            x1 = 0.25 * (grid_x(i_f, j_f) + grid_x(i_f + 1, j_f) + &
                         grid_x(i_f, j_f + 1) + grid_x(i_f + 1, j_f + 1))
            y1 = 0.25 * (grid_y(i_f, j_f) + grid_y(i_f + 1, j_f) + &
                         grid_y(i_f, j_f + 1) + grid_y(i_f + 1, j_f + 1))
            d_f = abs(sqrt((x1 - x0)**2 + (y1 - y0)**2))
            d_f = max(1e-6, d_f)

            nx = iblines(bandface%nearest_ibline_ind)%no(1)
            ny = iblines(bandface%nearest_ibline_ind)%no(2)
            u_temp = x_speed(i_f, j_f)
            v_temp = y_speed(i_f, j_f)
            u_f = u_temp - (1 * (u_temp * nx + v_temp * ny) * nx)
            v_f = v_temp - (1 * (u_temp * nx + v_temp * ny) * ny)
            
            u_temp = iblines(bandface%nearest_ibline_ind)%u
            v_temp = iblines(bandface%nearest_ibline_ind)%v
            u_ib = u_temp - (1 * (u_temp * nx + v_temp * ny) * nx)
            v_ib = v_temp - (1 * (u_temp * nx + v_temp * ny) * ny)

            if (f_type .eq. 'x') then
                x_x_speed_left(i, j) = ((u_f/d_f) + (u_ib/d_ib)) / (1/d_f + 1/d_ib)
                x_y_speed_left(i, j) = ((v_f/d_f) + (v_ib/d_ib)) / (1/d_f + 1/d_ib)
                x_x_speed_right(i, j) = x_x_speed_left(i, j)
                x_y_speed_right(i, j) = x_y_speed_left(i, j)
            else
                y_x_speed_left(i, j) = ((u_f/d_f) + (u_ib/d_ib)) / (1/d_f + 1/d_ib)
                y_y_speed_left(i, j) = ((v_f/d_f) + (v_ib/d_ib)) / (1/d_f + 1/d_ib)
                y_x_speed_right(i, j) = y_x_speed_left(i, j)
                y_y_speed_right(i, j) = y_y_speed_left(i, j)
            end if
            
            ! Adding perpendicular velocity booundary condition 
            df = bandface%field_min_d
            df = df + sign(1.0, df)*1e-15
            dface = bandface%face_min_d
            dface = min(0.0, dface)
            
            u_temp = x_speed(i_f, j_f)
            v_temp = y_speed(i_f, j_f)
            u_f = (u_temp * nx + v_temp * ny) * nx
            v_f = (u_temp * nx + v_temp * ny) * ny
            
            u_temp = iblines(bandface%nearest_ibline_ind)%u
            v_temp = iblines(bandface%nearest_ibline_ind)%v
            u_ib = (u_temp * nx + v_temp * ny) * nx
            v_ib = (u_temp * nx + v_temp * ny) * ny

            if (f_type .eq. 'x') then
                x_x_speed_left(i, j) = x_x_speed_left(i, j)+u_ib+dface/df*(u_f-u_ib)
                x_y_speed_left(i, j) = x_y_speed_left(i, j)+v_ib+dface/df*(v_f-v_ib)
                x_x_speed_right(i, j) = x_x_speed_left(i, j)
                x_y_speed_right(i, j) = x_y_speed_left(i, j)
            else
                y_x_speed_left(i, j) = y_x_speed_left(i, j)+u_ib+dface/df*(u_f-u_ib)
                y_y_speed_left(i, j) = y_y_speed_left(i, j)+v_ib+dface/df*(v_f-v_ib)
                y_x_speed_right(i, j) = y_x_speed_left(i, j)
                y_y_speed_right(i, j) = y_y_speed_left(i, j)
            end if
        end do

    end subroutine set_no_slip_velocity_face

    subroutine set_pressure_and_density_face()

        implicit none

        integer :: i, j, i_f, j_f, i_b, j_b, n
        character :: f_type
        type(face) :: bandface
        
        call dmsg(1, 'immersed_boundary', 'set_pressure_and_density_at_interfaces')
        
        do n = 1, num_bandfaces
            bandface = bandfaces(n)
            i = bandface%face_ind(1)
            j = bandface%face_ind(2)
            i_f = bandface%field_ind(1)
            j_f = bandface%field_ind(2)
            i_b = bandface%band_ind(1)
            j_b = bandface%band_ind(2)
            f_type = bandface%face_type
            ! Set velocity
            if (f_type .eq. 'x') then
                ! xi face
                x_density_left(i,j) = density(i_f, j_f)
                x_density_right(i,j) = density(i_f, j_f)
                x_pressure_left(i,j) = pressure(i_f, j_f)
                x_pressure_right(i,j) = pressure(i_f, j_f)
            else
                ! eta face
                y_density_left(i,j) = density(i_f, j_f)
                y_density_right(i,j) = density(i_f, j_f)
                y_pressure_left(i,j) = pressure(i_f, j_f)
                y_pressure_right(i,j) = pressure(i_f, j_f)
            end if
        end do

    end subroutine set_pressure_and_density_face

    subroutine reset_states_at_interface_faces()

        implicit none

        call dmsg(1, 'immersed_boundary', 'reset_states_at_interface_faces')
        
        if (mu_ref .eq. 0.0) then
            call set_slip_velocity_face()
        else
            call set_no_slip_velocity_face()
        end if

        call set_pressure_and_density_face()

    end subroutine reset_states_at_interface_faces

    subroutine update_band_interior_cells()

        implicit none

        type(face) :: bandface
        integer :: n, i, j, i_f, j_f, i_b, j_b
        real :: u_temp, v_temp
        character :: f_type

        call dmsg(1, 'immersed_boundary', 'update_band_interior_cells')
        
        do n = 1, num_interiorcells
            i = interiorcells(n)%cell_ind(1)
            j = interiorcells(n)%cell_ind(2)
            density(i, j) = density_inf
            pressure(i, j) = pressure_inf
            x_speed(i, j) = 0.
            y_speed(i, j) = 0.
            residue(i, j, :) = (/ 0., 0., 0., 0. /)
        end do

        do n = 1, num_bandfaces
            bandface = bandfaces(n)
            i = bandface%face_ind(1)
            j = bandface%face_ind(2)
            i_f = bandface%field_ind(1)
            j_f = bandface%field_ind(2)
            i_b = bandface%band_ind(1)
            j_b = bandface%band_ind(2)
            f_type = bandface%face_type
            
            density(i_b, j_b) = density(i_f, j_f)
            pressure(i_b, j_b) = pressure(i_f, j_f)
            if (f_type .eq. 'x') then
                u_temp = x_x_speed_left(i, j)
                v_temp = x_y_speed_left(i, j)
            else
                u_temp = y_x_speed_left(i, j)
                v_temp = y_y_speed_left(i, j)
            end if
            x_speed(i_b, j_b) = (2 * u_temp) - x_speed(i_f, j_f)
            y_speed(i_b, j_b) = (2 * v_temp) - y_speed(i_f, j_f)
            residue(i_b, j_b, :) = (/ 0., 0., 0., 0. /)
        end do

    end subroutine update_band_interior_cells

    subroutine reset_gradients_at_interfaces()
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the 
    ! interfaces between band and field cells. Using the same
    ! scheme as that in viscous.f90.
    !-----------------------------------------------------------------

        implicit none
        
        real :: Tau_xx, Tau_xy, Tau_yy
        real :: dudx, dudy, dvdx, dvdy, dTdx, dTdy 
        real :: K_heat, mu, Qx, Qy
        character :: f_type
        integer :: i, j, i_f, j_f, i_b, j_b, n, dir_flag
        real :: nx, ny, A, Tfield, uface, vface, Tface
        real, dimension(:), pointer :: Flux_p
        type(face) :: bandface

        call dmsg(1, 'immersed_boundary', 'reset_gradients')

        ! Calculating the fluxes at the faces
        ! NOTE: The normal vector should be in the direction from band to field
        ! dir_flag ensures the above
        do n = 1, num_bandfaces
            bandface = bandfaces(n)
            i = bandface%face_ind(1)
            j = bandface%face_ind(2)
            i_f = bandface%field_ind(1)
            j_f = bandface%field_ind(2)
            i_b = bandface%band_ind(1)
            j_b = bandface%band_ind(2)
            f_type = bandface%face_type
            Tfield = pressure(i_f, j_f) / (density(i_f, j_f) * R_gas)
            dir_flag = (i_b - i_f) + (j_b - j_f)

            if (f_type .eq. 'x') then
                nx = xnx(i, j)
                ny = xny(i, j)
                A = xA(i, j)
                Flux_p => F_p(i, j, :)
                uface = 0.5 * (x_x_speed_left(i, j) + x_x_speed_right(i, j))
                vface = 0.5 * (x_y_speed_left(i, j) + x_y_speed_right(i, j))
                Tface = 0.5 * ((x_pressure_left(i, j) / (x_density_left(i, j) * R_gas)) + &
                               (x_pressure_right(i, j) / (x_density_right(i, j) * R_gas)))
            else
                nx = ynx(i, j)
                ny = yny(i, j)
                A = yA(i, j)
                Flux_p => G_p(i, j, :)
                uface = 0.5 * (y_x_speed_left(i, j) + y_x_speed_right(i, j))
                vface = 0.5 * (y_y_speed_left(i, j) + y_y_speed_right(i, j))
                Tface = 0.5 * ((y_pressure_left(i, j) / (y_density_left(i, j) * R_gas)) + &
                               (y_pressure_right(i, j) / (y_density_right(i, j) * R_gas)))
            end if
            
            dudx = 4 * dir_flag * (uface - x_speed(i_f, j_f)) * A * nx / (volume(i_f, j_f) + volume(i_b, j_b))
            dudy = 4 * dir_flag * (uface - x_speed(i_f, j_f)) * A * ny / (volume(i_f, j_f) + volume(i_b, j_b))
            dvdx = 4 * dir_flag * (vface - y_speed(i_f, j_f)) * A * nx / (volume(i_f, j_f) + volume(i_b, j_b))
            dvdy = 4 * dir_flag * (vface - y_speed(i_f, j_f)) * A * ny / (volume(i_f, j_f) + volume(i_b, j_b))
            dTdx = 4 * dir_flag * (Tface - Tfield) * A * nx / (volume(i_f, j_f) + volume(i_b, j_b))
            dTdy = 4 * dir_flag * (Tface - Tfield) * A * ny / (volume(i_f, j_f) + volume(i_b, j_b))
            
            ! mu requires T at the face. Hence:
            ! T_L and T_R are the left and right states of the face i,j,k
            ! The values at face used instead of element values
            mu = mu_ref * (Tface / T_ref)**1.5 * ((T_ref + &
                          Sutherland_temp) / (Tface + Sutherland_temp))

            ! Using lambda = -2 * mu / 3
            ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
            ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz
            Tau_xx = 2. * mu * (dudx - ((dudx + dvdy) / 3.)) 
            Tau_yy = 2. * mu * (dvdy - ((dudx + dvdy) / 3.)) 
            Tau_xy = mu * (dvdx + dudy)

            ! Pr: Prandtl Number
            ! Qx, Qy, Qz: Conduction fluxes
            K_heat = mu * gm * R_gas / ((gm - 1) * Pr)
            Qx = K_heat*dTdx
            Qy = K_heat*dTdy

            ! Note that the xi-direction faces only need the following quantities:
            ! Tau_xx, Tau_xy, Tau_xz -> dudx, dudy, dudz, dvdx, dvdy, dwdx, dwdz
            ! Qx -> dTdx
            ! The mass flux has no viscous component
            ! momentum for xi-face:
            Flux_p(2) = - ( ((Tau_xx * nx) + &
                         (Tau_xy * ny)) * A)
            Flux_p(3) = - ( ((Tau_xy * nx) + &
                         (Tau_yy * ny)) * A)
           
            ! Energy flux
            ! (TijVj + Qi)ni
            Flux_p(4) = - (A * ( &
                            ((Tau_xx*uface + Tau_xy*vface + Qx) * &
                             nx) + &
                            ((Tau_xy*uface + Tau_yy*vface + Qy) * &
                             ny) ) )
        end do

    end subroutine reset_gradients_at_interfaces

! ---------------------------------------------------------------------------------
! NOTE: Testing functions which validate cell classification
! Use the output files band_cells, band_faces, interior cells and put them in
! python codes/Immersed Boundary folder
! Open ImmersedBoundary.py and give file names of grid file and IB file.
! ensure that plot_from_file() function is un-commented
! Running ImmersedBoundary.py will plot the band cells, interior cells, band faces
! ---------------------------------------------------------------------------------
!   subroutine test_classification()

!       implicit none

!       integer :: i
!       open(75, file='band_cells')
!       do i = 1, num_bandcells
!           write(75, *) bandcells(i)%cell_ind(1), bandcells(i)%cell_ind(2)
!       end do
!       close(75)
!       open(75, file='interior_cells')
!       do i = 1, num_interiorcells
!           write(75, *) interiorcells(i)%cell_ind(1), interiorcells(i)%cell_ind(2)
!       end do
!       close(75)

!   end subroutine test_classification

!   subroutine test_bandfaces()

!       implicit none

!       integer :: i, j, i_f, j_f, i_b, j_b, n
!       character :: f_type
!       type(face) :: bandface
!       
!       call dmsg(1, 'immersed_boundary', 'set_pressure_and_density_at_interfaces')
!       
!       open(75, file='band_faces')
!       do n = 1, num_bandfaces
!           bandface = bandfaces(n)
!           i = bandface%face_ind(1)
!           j = bandface%face_ind(2)
!           i_f = bandface%field_ind(1)
!           j_f = bandface%field_ind(2)
!           i_b = bandface%band_ind(1)
!           j_b = bandface%band_ind(2)
!           f_type = bandface%face_type
!           write(75, *) f_type, i, j, i_f, j_f, i_b, j_b
!       end do
!       close(75)

!   end subroutine test_bandfaces
    
    subroutine destroy_IB()

        implicit none

        integer :: mem_stat
        
        call dmsg(1, 'immersed_boundary', 'destroy_IB')
        
        deallocate(iblines, stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not deallocate iblines'
        end if

        deallocate(bandcells, stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not deallocate bandcells'
        end if

        deallocate(interiorcells, stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not deallocate interiorcells'
        end if

        deallocate(tempcells, stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not deallocate tempcells'
        end if

        deallocate(bandfaces, stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not deallocate bandfaces'
        end if

        deallocate(first_cut_indices, stat=mem_stat)
        if (mem_stat /= 0) then
            print *, 'Error: Could not deallocate first_cut_indices'
        end if

    end subroutine destroy_IB

end module immersed_boundary
