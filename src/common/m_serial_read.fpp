module m_serial_read

    use m_derived_types

    use m_global_parameters

    use m_variables_conversion

    use m_mpi_common

    use m_compile_specific

    use m_boundary_common

    private; public :: s_read_serial_data_files, &
        s_read_serial_grid_binary

    logical :: file_exist !< checks if file exists
    character(LEN=15) :: FMT !< Format string for writing data files
    character(LEN=3) :: status
    integer :: i, j, k, l, r, c !< Generic loop iterator
    logical, parameter :: serial_txt_io = .false.

contains

    !> The goal of this subroutine is to read in any preexisting
        !!      initial condition data files so that they may be used by
        !!      the pre-process as a starting point in the creation of an
        !!      all new initial condition.
        !! @param q_cons_vf Conservative variables
        !! @param ib_markers track if a cell is within the immersed boundary
    subroutine s_read_serial_data_files(t_step_dir, q_cons_vf, ib_markers, pb, mv, bc_type, ib_levelset, &
                                        ib_levelset_norm, airfoil_grid_u, airfoil_grid_l)

        character(len=*), intent(in) :: t_step_dir
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), intent(inout) :: ib_markers
        type(pres_field), intent(inout), optional :: pb, mv
        type(integer_field), dimension(1:num_dims, -1:1), intent(inout), optional :: bc_type
        type(levelset_field), intent(IN), optional :: ib_levelset
        type(levelset_norm_field), intent(IN), optional :: ib_levelset_norm
        type(vec3_dt), dimension(:), allocatable, intent(inout), optional :: airfoil_grid_u, airfoil_grid_l

        call s_read_serial_conservative_Variables_binary(t_step_dir, q_cons_vf)
        if (ib) call s_read_serial_ib_binary(t_step_dir, ib_markers, ib_levelset, ib_levelset_norm, &
                                             airfoil_grid_u, airfoil_grid_l)

#ifndef MFC_POST_PROCESS
        if (qbmm .and. .not. polytropic) call s_read_serial_nonpolytropic_qbmm_binary(t_step_dir, pb, mv)
#endif

        if (bc_io .and. present(bc_type)) then
            call s_read_serial_boundary_condition_files(t_step_dir, bc_type)
        else
            call s_assign_default_bc_type(bc_type)
        end if

    end subroutine s_read_serial_data_files

    !> The goal of this subroutine is to read in any preexisting
        !!      grid data as well as based on the imported grid, complete
        !!      the necessary global computational domain parameters.
    impure subroutine s_read_serial_grid_binary(step_dirpath)

        ! Generic string used to store the address of a particular file
        character(len=*), intent(in) :: step_dirpath
        character(len=len_trim(step_dirpath) + 2*name_len) :: file_loc

        ! Inquiring as to the existence of the time-step directory
        file_loc = trim(step_dirpath)//'/.'
        call my_inquire(file_loc, file_exist)

        ! If the time-step directory is missing, the pre-process exits
        if (.not. file_exist) then
            call s_mpi_abort('Time-step folder '//trim(step_dirpath)// &
                             ' is missing. Exiting.')
        end if

        ! Reading the Grid Data File for the x-direction
        #:for VAR, XYZ in [('m', 'x'), ('n', 'y'), ('p', 'z')]
            if (${VAR}$ > 0) then
                ! Checking whether x[y,z]_cb.dat exists
                file_loc = trim(step_dirpath)//'/${XYZ}$_cb.dat'
                call my_inquire (trim(file_loc), file_exist)

                ! If it exists, x[y,z]_cb.dat is read
                if (file_exist) then
                    open (1, FILE=trim(file_loc), FORM='unformatted', &
                          STATUS='old', ACTION='read')
                    read (1) ${XYZ}$_cb(-1:${VAR}$)
                    close (1)
                else
                    call s_mpi_abort('File ${XYZ}$_cb.dat is missing in '// &
                                     trim(step_dirpath)//'. Exiting.')
                end if

#ifndef MFC_PRE_PROCESS
                d${XYZ}$(0:${VAR}$) = ${XYZ}$_cb(0:${VAR}$) - ${XYZ}$_cb(-1:${VAR}$ - 1)
                ${XYZ}$_cc(0:${VAR}$) = ${XYZ}$_cb(-1:${VAR}$ - 1) + d${XYZ}$(0:${VAR}$)/2._wp
#else
                ! Computing cell-center locations
                ${XYZ}$_cc(0:${VAR}$) = (${XYZ}$_cb(0:${VAR}$) + ${XYZ}$_cb(-1:(${VAR}$ - 1)))/2._wp

                ! Computing minimum cell-width
                d${XYZ}$ = minval(${XYZ}$_cb(0:${VAR}$) - ${XYZ}$_cb(-1:${VAR}$ - 1))
                if (num_procs > 1) call s_mpi_reduce_min(d${XYZ}$)

                ! Setting locations of domain bounds
                ${XYZ}$_domain%beg = ${XYZ}$_cb(-1)
                ${XYZ}$_domain%end = ${XYZ}$_cb(${VAR}$)
#endif
            end if
        #:endfor

    end subroutine s_read_serial_grid_binary

    subroutine s_read_serial_ib_binary(step_dirpath, ib_markers, ib_levelset, ib_levelset_norm, &
                                       airfoil_grid_u, airfoil_grid_l)

        character(len=*), intent(in) :: step_dirpath
        type(integer_field), intent(inout) :: ib_markers
        type(levelset_field), intent(IN), optional :: ib_levelset
        type(levelset_norm_field), intent(IN), optional :: ib_levelset_norm
        type(vec3_dt), dimension(:), allocatable, intent(inout), optional :: airfoil_grid_u, airfoil_grid_l
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc !<

        ! Read Markers
        file_loc = trim(step_dirpath)//'/ib.dat'
        call my_inquire(trim(file_loc), file_exist)
        if (file_exist) then
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS='old', ACTION='read')
            read (1) ib_markers%sf(0:m, 0:n, 0:p)
            close (1)
        else
            call s_mpi_abort('File ib.dat is missing in ' &
                             //trim(step_dirpath)// &
                             '. Exiting.')
        end if

        ! Read Levelset
        if (present(ib_levelset)) then
            write (file_loc, '(A)') &
                trim(step_dirpath)//'/levelset.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist) then
                open (2, FILE=trim(file_loc), FORM='unformatted', ACTION='read')
                read (2) ib_levelset%sf(0:m, 0:n, 0:p, 1:num_ibs);
                close (2)
            else
                call s_mpi_abort(trim(file_loc)//' is missing. Exiting.')
            end if
        end if

        ! Read Levelset Norm
        if (present(ib_levelset_norm)) then
            write (file_loc, '(A)') trim(step_dirpath)//'/levelset_norm.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist) then
                open (2, FILE=trim(file_loc), FORM='unformatted', ACTION='read')
                read (2) ib_levelset_norm%sf(0:m, 0:n, 0:p, 1:num_ibs, 1:3); close (2)
            else
                call s_mpi_abort(trim(file_loc)//' is missing. Exiting.')
            end if
        end if

#ifndef MFC_POST_PROCESS
        ! Read Airfoil Variables
        if (present(airfoil_grid_u) .and. present(airfoil_grid_L)) then
            do i = 1, num_ibs
                if (patch_ib(i)%geometry == 4) then
#ifdef MFC_SIMULATION
                    ! In pre_process Np is calculated in m_ib_patches
                    Np = int((patch_ib(i)%p*patch_ib(i)%c/dx(0))*20) + int(((patch_ib(i)%c - patch_ib(i)%p*patch_ib(i)%c)/dx(0))*20) + 1
#endif
                    #:for VAR in ['u', 'l']
                        allocate (airfoil_grid_${VAR}$(1:Np))
                        write (file_loc, '(A)') trim(step_dirpath)//'/airfoil_${VAR}$.dat'
                        inquire (FILE=trim(file_loc), EXIST=file_exist)
                        if (file_exist) then
                            open (2, FILE=trim(file_loc), FORM='unformatted', ACTION='read')
                            read (2) airfoil_grid_${VAR}$; close (2)
                        else
                            call s_mpi_abort(trim(file_loc)//' is missing. Exiting.')
                        end if
                    #:endfor
                end if
            end do
        end if
#endif

    end subroutine s_read_serial_ib_binary

    subroutine s_read_serial_conservative_variables_binary(step_dirpath, q_cons_vf)

        character(len=*), intent(in) :: step_dirpath
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc
        character(LEN=int(floor(log10(real(sys_size, wp)))) + 1) :: file_num

        do i = 1, sys_size
            ! Checking whether data file associated with variable position
            ! of the currently manipulated conservative variable exists
            write (file_num, '(I0)') i
            file_loc = trim(step_dirpath)//'/q_cons_vf'// &
                       trim(file_num)//'.dat'
            call my_inquire(trim(file_loc), file_exist)

            ! If it exists, the data file is read
            if (file_exist) then
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS='old', ACTION='read')
                read (1) q_cons_vf(i)%sf(0:m, 0:n, 0:p)
                close (1)
            else
                call s_mpi_abort('File q_cons_vf'//trim(file_num)// &
                                 '.dat is missing in '//trim(step_dirpath)// &
                                 '. Exiting.')
            end if
        end do

    end subroutine s_read_serial_conservative_variables_binary

    subroutine s_read_serial_nonpolytropic_qbmm_binary(step_dirpath, pb, mv)

        character(len=*), intent(in) :: step_dirpath
        type(pres_field), intent(in) :: pb, mv
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc
        character(LEN=int(floor(log10(real(sys_size, wp)))) + 1) :: file_num

        do i = 1, nb
            do r = 1, nnode
                ! Checking whether data file associated with variable position
                ! of the currently manipulated bubble variable exists
                write (file_num, '(I0)') sys_size + r + (i - 1)*nnode
                file_loc = trim(step_dirpath)//'/pb'// &
                           trim(file_num)//'.dat'
                call my_inquire(trim(file_loc), file_exist)

                ! If it exists, the data file is read
                if (file_exist) then
                    open (1, FILE=trim(file_loc), FORM='unformatted', &
                          STATUS='old', ACTION='read')
                    read (1) pb%sf(0:m, 0:n, 0:p, r, i)
                    close (1)
                else
                    call s_mpi_abort('File pb'//trim(file_num)// &
                                     '.dat is missing in '//trim(step_dirpath)// &
                                     '. Exiting.')
                end if
            end do

        end do

        do i = 1, nb
            do r = 1, 4
                ! Checking whether data file associated with variable position
                ! of the currently manipulated bubble variable exists
                write (file_num, '(I0)') sys_size + r + (i - 1)*4
                file_loc = trim(step_dirpath)//'/mv'// &
                           trim(file_num)//'.dat'
                call my_inquire (trim(file_loc), file_exist)

                ! If it exists, the data file is read
                if (file_exist) then
                    open (1, FILE=trim(file_loc), FORM='unformatted', &
                          STATUS='old', ACTION='read')
                    read (1) mv%sf(0:m, 0:n, 0:p, r, i)
                    close (1)
                else
                    call s_mpi_abort('File mv'//trim(file_num)// &
                                     '.dat is missing in '//trim(step_dirpath)// &
                                     '. Exiting.')
                end if
            end do

        end do

    end subroutine s_read_serial_nonpolytropic_qbmm_binary

end module m_serial_read
