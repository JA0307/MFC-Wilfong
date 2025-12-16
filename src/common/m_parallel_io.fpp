module m_parallel_io

    use m_derived_types

    use m_global_parameters

    use m_mpi_common

    private; public :: s_read_parallel_grid_data_files

contains

    !> Cell-boundary data are checked for consistency by looking
        !!      at the (non-)uniform cell-width distributions for all the
        !!      active coordinate directions and making sure that all of
        !!      the cell-widths are positively valued
    impure subroutine s_read_parallel_grid_data_files()

#ifdef MFC_MPI

        real(wp), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        #:for VAR, IDX, XYZ in [('m', 1, 'x'), ('n', 2, 'y'), ('p', 3, 'z')]
            if (${VAR}$ > 0) then
                ! Read in cell boundary locations in x-direction
                file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'${XYZ}$_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = ${VAR}$_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ_ALL(ifile, ${XYZ}$_cb_glb, data_size, mpi_p, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting. ')
                end if

#ifdef MFC_PRE_PROCESS
                ! Assigning local cell boundary locations
                ${XYZ}$_cb(-1:${VAR}$) = ${XYZ}$_cb_glb((start_idx(${IDX}$) - 1):(start_idx(${IDX}$) + ${VAR}$))
                ! Computing cell center locations
                ${XYZ}$_cc(0:${VAR}$) = (${XYZ}$_cb(0:${VAR}$) + ${XYZ}$_cb(-1:(${VAR}$ - 1)))/2._wp
                ! Computing minimum cell width
                d${XYZ}$ = minval(${XYZ}$_cb(0:${VAR}$) - ${XYZ}$_cb(-1:(${VAR}$ - 1)))
                if (num_procs > 1) call s_mpi_reduce_min(d${XYZ}$)
                ! Setting locations of domain bounds
                ${XYZ}$_domain%beg = ${XYZ}$_cb(-1)
                ${XYZ}$_domain%end = ${XYZ}$_cb(m)
#else
                ! Assigning local cell boundary locations
                ${XYZ}$_cb(-1:${VAR}$) = ${XYZ}$_cb_glb((start_idx(${IDX}$) - 1):(start_idx(${IDX}$) + ${VAR}$))
                ! Computing the cell width distribution
                d${XYZ}$(0:${VAR}$) = ${XYZ}$_cb(0:${VAR}$) - ${XYZ}$_cb(-1:${VAR}$ - 1)
                ! Computing the cell center locations
                ${XYZ}$_cc(0:${VAR}$) = ${XYZ}$_cb(-1:${VAR}$ - 1) + d${XYZ}$(0:${VAR}$)/2._wp
#endif

            end if
        #:endfor

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

#endif

    end subroutine s_read_parallel_grid_data_files

end module m_parallel_io
