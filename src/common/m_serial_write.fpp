
#:include 'case.fpp'
#:include 'macros.fpp'

module m_serial_write

    use m_derived_types

    use m_global_parameters

    use m_variables_conversion

    use m_mpi_common

    use m_compile_specific

    use m_boundary_common

    logical :: file_exist !< checks if file exists
    character(LEN=15) :: FMT !< Format string for writing data files
    character(LEN=3) :: status
    integer :: i, j, k, l, r, c !< Generic loop iterator
    logical, parameter :: serial_txt_io = .true.

contains

    !>  Writes grid and initial condition data files to the "0"
        !!  time-step directory in the local processor rank folder
        !! @param q_cons_vf Conservative variables
        !! @param ib_markers track if a cell is within the immersed boundary
        !! @param levelset closest distance from every cell to the IB
        !! @param levelset_norm normalized vector from every cell to the closest point to the IB
    impure subroutine s_write_serial_data_files(t_step_dir, save_count, q_cons_vf, q_prim_vf, bc_type, pb, mv, ib_markers, q_T_sf)

        character(len=*), intent(in) :: t_step_dir
        integer, intent(in) :: save_count
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf, q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type
        type(pres_field), intent(in) :: pb, mv
        type(integer_field), intent(in) :: ib_markers
        type(scalar_field), intent(inout), optional :: q_T_sf

        character(LEN=len_trim(t_step_dir) + name_len) :: file_loc, temp_dir
        integer :: t_step !< Time step number
#ifndef MFC_POST_PROCESS

        call my_inquire(t_step_dir, file_exist)
        if (file_exist) call s_delete_directory(trim(t_step_dir))
        call s_create_directory(trim(t_step_dir))

#ifdef MFC_PRE_PROCESS
        status = merge('old', 'new', old_grid)
        if (bc_io) then
            if (igr) then
                call s_write_serial_boundary_condition_files(q_cons_vf, bc_type, t_step_dir, status)
            else
                call s_write_serial_boundary_condition_files(q_prim_vf, bc_type, t_step_dir, status)
            end if
        end if
#else
        status = 'new'
#endif

        call s_write_serial_grid_binary(t_step_dir)

        call s_write_serial_conservative_variables_binary(t_step_dir, q_cons_vf)
        if (qbmm .and. .not. polytropic) call s_write_serial_nonpolytropic_qbmm_binary(t_step_dir, pb, mv)
        if (ib) call s_write_serial_ib_binary(t_step_dir, ib_markers)

        if (serial_txt_io) then
            ! Query if time-step directory exists, if not create it
            write (file_loc, '(A,I0,A,I0)') trim(case_dir)//'/D'
            if (proc_rank == 0) then
                inquire (FILE=trim(file_loc), EXIST=file_exist)
                if (.not. file_exist) call s_create_directory(trim(file_loc))
            end if

            if (num_dims == 1) then
                if (precision == 1) then
                    FMT = "(2F30.7)"
                else
                    FMT = "(2F40.14)"
                end if
            elseif (num_dims == 2) then
                if (precision == 1) then
                    FMT = "(3F30.7)"
                else
                    FMT = "(3F40.14)"
                end if
            elseif (num_dims == 3) then
                if (precision == 1) then
                    FMT = "(4F30.7)"
                else
                    FMT = "(4F40.14)"
                end if
            end if

#ifdef MFC_PRE_PROCESS
            if (num_dims == 1) call s_write_prim_variables_txt_pre_process(file_loc, save_count, q_cons_vf)
#else
            call s_write_prim_variables_txt(file_loc, save_count, q_cons_vf, q_prim_vf, q_T_sf)
#endif

            call s_write_serial_conservative_variables_txt(file_loc, save_count, q_cons_vf)
            if (qbmm .and. .not. polytropic) call s_write_serial_nonpolytropic_qbmm_txt(file_loc, save_count, pb, mv)
            if (ib) call s_write_serial_ib_txt(file_loc, save_count, ib_markers)
        end if
#endif
    end subroutine s_write_serial_data_files

    subroutine s_write_serial_grid_binary(step_dirpath)

        character(LEN=*), intent(in) :: step_dirpath
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc
#ifndef MFC_POST_PROCESS
        #:for VAR, XYZ in [('m', 'x'), ('n', 'y'), ('p', 'z')]
            if (${VAR}$ > 0) then
                file_loc = trim(step_dirpath)// '/${XYZ}$_cb.dat'
                open  (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
                write (1) ${XYZ}$_cb(-1:${VAR}$)
                close (1)
            endif
        #:endfor
#endif
    end subroutine s_write_serial_grid_binary

    subroutine s_write_serial_ib_binary(step_dirpath, ib_markers)

        character(len=*), intent(in) :: step_dirpath
        type(integer_field), intent(in) :: ib_markers
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc !<
#ifndef MFC_POST_PROCESS
        ! Outputting IB Markers
        file_loc = trim(step_dirpath)// '/ib.dat'
        open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
        write (1) ib_markers%sf
        close (1)

        ! Write airfoil specific data
        do i = 1, num_ibs
            if (patch_ib(i)%geometry == 4) then
                #:for VAR in ['u', 'l']
                    file_loc = trim(step_dirpath)//'/airfoil_${VAR}$.dat'
                    open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
                    write (1) airfoil_grid_${VAR}$(1:Np)
                    close (1)
                #:endfor
            end if
        end do
#endif

    end subroutine s_write_serial_ib_binary

    subroutine s_write_serial_ib_txt(step_dirpath, save_count, ib_markers)

        character(len=*), intent(in) :: step_dirpath
        integer, intent(in) :: save_count
        type(integer_field), intent(in) :: ib_markers
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc !<
        !! Generic string used to store the address of a particular file
#ifndef MFC_POST_PROCESS
        ! Write IB Markers
        write (file_loc, '(A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/ib_markers.', proc_rank, '.', save_count, '.dat'
        open (2, FILE=trim(file_loc))
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    if (p > 0) then
                        write (2, FMT) x_cc(j), y_cc(k), z_cc(l), real(ib_markers%sf(j, k, l))
                    else
                        write (2, FMT) x_cc(j), y_cc(k), real(ib_markers%sf(j, k, l))
                    end if
                end do
            end do
        end do
        close (2)

        ! Write airfoil specific data
        do i = 1, num_ibs
            if (patch_ib(i)%geometry == 4) then
                write (file_loc, '(A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/airfoil_u.', proc_rank, '.', save_count, '.dat'
                open (2, FILE=trim(file_loc))
                do j = 1, Np
                    write (2, FMT) airfoil_grid_u(j)%x, airfoil_grid_u(j)%y
                end do
                close (2)

                write (file_loc, '(A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/airfoil_l.', proc_rank,  '.', save_count, '.dat'
                open (2, FILE=trim(file_loc))
                do j = 1, Np
                    write (2, FMT) airfoil_grid_l(j)%x, airfoil_grid_l(j)%y
                end do
                close (2)
            end if
        end do
#endif
    end subroutine s_write_serial_ib_txt

    subroutine s_write_serial_conservative_variables_binary(step_dirpath, q_cons_vf)

        character(len=*), intent(in) :: step_dirpath
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc !<
        character(LEN=int(floor(log10(real(sys_size, wp)))) + 1) :: file_num
#ifndef MFC_POST_PROCESS
        do i = 1, sys_size
            write (file_num, '(I0)') i
            file_loc = trim(step_dirpath)//'/q_cons_vf'//trim(file_num) &
                       //'.dat'
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS=status)
            write (1) q_cons_vf(i)%sf(0:m, 0:n, 0:p)
            close (1)
        end do
#endif
    end subroutine s_write_serial_conservative_variables_binary

    subroutine s_write_serial_conservative_variables_txt(step_dirpath, save_count, q_cons_vf)

        character(len=*), intent(in) :: step_dirpath
        integer :: save_count
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc !<
        !! Generic string used to store the address of a particular file
#ifndef MFC_POST_PROCESS
        do i = 1, sys_size
            write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/cons.', i, '.', proc_rank, '.', save_count , '.dat'
            open (2, FILE=trim(file_loc))
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        if (num_dims == 1) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        elseif (num_dims == 2) then
                            write (2, FMT) x_cb(j), y_cb(k), q_cons_vf(i)%sf(j, k, 0)
                        elseif (num_dims == 3) then
                            write (2, FMT) x_cb(j), y_cb(k), z_cb(l), q_cons_vf(i)%sf(j, k, l)
                        end if
                    end do
                    if (p > 0) write (2, *)
                end do
                if (n > 0) write (2, *)
            end do
            close (2)
        end do
#endif
    end subroutine s_write_serial_conservative_variables_txt

    subroutine s_write_serial_nonpolytropic_qbmm_binary(step_dirpath, pb, mv)

        character(len=*), intent(in) :: step_dirpath
        type(pres_field), intent(in) :: pb, mv
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc
        character(LEN=int(floor(log10(real(sys_size, wp)))) + 1) :: file_num
#ifndef MFC_POST_PROCESS
        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do r = 1, nnode
                    write (file_num, '(I0)') r + (i - 1)*nnode + sys_size
                    file_loc = trim(step_dirpath) // '/pb' // trim(file_num) // '.dat'
                    open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
                    write (1) pb%sf(0:m, 0:n, 0:p, r, i)
                    close (1)
                end do
            end do

            do i = 1, nb
                do r = 1, nnode
                    write (file_num, '(I0)') r + (i - 1)*nnode + sys_size
                    file_loc = trim(step_dirpath) // '/mv' // trim(file_num) // '.dat'
                    open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
                    write (1) mv%sf(0:m, 0:n, 0:p, r, i)
                    close (1)
                end do
            end do
        end if
#endif
    end subroutine s_write_serial_nonpolytropic_qbmm_binary

    subroutine s_write_serial_nonpolytropic_qbmm_txt(step_dirpath, save_count, pb, mv)

        character(len=*), intent(in) :: step_dirpath
        integer :: save_count
        type(pres_field), intent(in) :: pb, mv
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc
#ifndef MFC_POST_PROCESS
        do i = 1, nb
            do r = 1, nnode
                write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/pres.', i, '.', r, '.', proc_rank, '.', save_count, '.dat'
                open (2, FILE=trim(file_loc))
                do j = 0, m
                    do k = 0, n
                        do l = 0, p
                            if (num_dims == 1) then
                                write (2, FMT) x_cb(j), pb%sf(j, 0, 0, r, i)
                            elseif  (num_dims == 2) then
                                write (2, FMT) x_cb(j), y_cb(k), pb%sf(j, k, 0, r, i)
                            elseif (num_dims == 3) then
                                write (2, FMT) x_cb(j), y_cb(k), z_cb(l), pb%sf(j, k, l, r, i)
                            end if
                        end do
                        if (p > 0) write (2, *)
                    end do
                    if (n > 0) write (2, *)
                end do
                close (2)
            end do
        end do

        do i = 1, nb
            do r = 1, nnode
                write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/mv.', i, '.', r, '.', proc_rank, '.', save_count, '.dat'
                open (2, FILE=trim(file_loc))
                do j = 0, m
                    do k = 0, n
                        do l = 0, p
                            if (num_dims == 1) then
                                write (2, FMT) x_cb(j), mv%sf(j, 0, 0, r, i)
                            elseif (num_dims == 2) then
                                write (2, FMT) x_cb(j), y_cb(k), mv%sf(j, k, 0, r, i)
                            elseif (num_dims == 3) then
                                write (2, FMT) x_cb(j), y_cb(k), z_cb(l), mv%sf(j, k, l, r, i)
                            end if
                        end do
                        if (p > 0) write (2, *)
                    end do
                    if (n > 0) write (2, *)
                end do
                close (2)
            end do
        end do
#endif
    end subroutine s_write_serial_nonpolytropic_qbmm_txt

#ifdef MFC_PRE_PROCESS
    subroutine s_write_prim_variables_txt_pre_process(step_dirpath, save_count, q_cons_vf)

        character(len=*), intent(in) :: step_dirpath
        integer :: save_count
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc

        real(wp), dimension(nb) :: nRtmp         !< Temporary bubble concentration
        real(wp) :: nbub                         !< Temporary bubble number density
        real(wp) :: gamma, lit_gamma, pi_inf, qv !< Temporary EOS params
        real(wp) :: rho                          !< Temporary density
        real(wp) :: pres, T                         !< Temporary pressure
        real(wp) :: rhoYks(1:num_species) !< Temporary species mass fractions
        real(wp) :: pres_mag

        pres_mag = 0._wp

        T = dflt_T_guess

        gamma = gammas(1)
        lit_gamma = gs_min(1)
        pi_inf = pi_infs(1)
        qv = qvs(1)

        if (model_eqns == 2) then
            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/prim.', i, '.', proc_rank, '.', save_count, '.dat'
                print*, proc_rank, file_loc
                open (2, FILE=trim(file_loc))
                do j = 0, m

                    if (chemistry) then
                        do c = 1, num_species
                            rhoYks(c) = q_cons_vf(chemxb + c - 1)%sf(j, 0, 0)
                        end do
                    end if

                    call s_convert_to_mixture_variables(q_cons_vf, j, 0, 0, rho, gamma, pi_inf, qv)

                    lit_gamma = 1._wp/gamma + 1._wp

                    if ((i >= chemxb) .and. (i <= chemxe)) then
                        write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)/rho
                    else if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                             .or. &
                             ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
                             .or. &
                             ((i >= chemxb) .and. (i <= chemxe)) &
                             ) then
                        write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                    else if (i == mom_idx%beg) then !u
                        write (2, FMT) x_cb(j), q_cons_vf(mom_idx%beg)%sf(j, 0, 0)/rho
                    else if (i == stress_idx%beg) then !tau_e
                        write (2, FMT) x_cb(j), q_cons_vf(stress_idx%beg)%sf(j, 0, 0)/rho
                    else if (i == E_idx) then !p
                        if (mhd) then
                            pres_mag = 0.5_wp*(Bx0**2 + q_cons_vf(B_idx%beg)%sf(j, 0, 0)**2 + q_cons_vf(B_idx%beg + 1)%sf(j, 0, 0)**2)
                        end if

                        call s_compute_pressure( &
                            q_cons_vf(E_idx)%sf(j, 0, 0), &
                            q_cons_vf(alf_idx)%sf(j, 0, 0), &
                            0.5_wp*(q_cons_vf(mom_idx%beg)%sf(j, 0, 0)**2._wp)/rho, &
                            pi_inf, gamma, rho, qv, rhoYks, pres, T, pres_mag=pres_mag)
                        write (2, FMT) x_cb(j), pres
                    else if (mhd) then
                        if (i == mom_idx%beg + 1) then ! v
                            write (2, FMT) x_cb(j), q_cons_vf(mom_idx%beg + 1)%sf(j, 0, 0)/rho
                        else if (i == mom_idx%beg + 2) then ! w
                            write (2, FMT) x_cb(j), q_cons_vf(mom_idx%beg + 2)%sf(j, 0, 0)/rho
                        else if (i == B_idx%beg) then ! By
                            write (2, FMT) x_cb(j), q_cons_vf(B_idx%beg)%sf(j, 0, 0)/rho
                        else if (i == B_idx%beg + 1) then ! Bz
                            write (2, FMT) x_cb(j), q_cons_vf(B_idx%beg + 1)%sf(j, 0, 0)/rho
                        end if
                    else if ((i >= bub_idx%beg) .and. (i <= bub_idx%end) .and. bubbles_euler) then

                        if (qbmm) then
                            nbub = q_cons_vf(bubxb)%sf(j, 0, 0)
                        else
                            if (adv_n) then
                                nbub = q_cons_vf(n_idx)%sf(j, 0, 0)
                            else
                                do k = 1, nb
                                    nRtmp(k) = q_cons_vf(bub_idx%rs(k))%sf(j, 0, 0)
                                end do

                                call s_comp_n_from_cons(real(q_cons_vf(alf_idx)%sf(j, 0, 0), kind=wp), nRtmp, nbub, weight)
                            end if
                        end if
                        write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)/nbub
                    else if (i == n_idx .and. adv_n .and. bubbles_euler) then
                        write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                    else if (i == damage_idx) then
                        write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                    end if
                end do
                close (2)
            end do
        end if

    end subroutine s_write_prim_variables_txt_pre_process
#elif defined(MFC_SIMULATION)
    subroutine s_write_prim_variables_txt(step_dirpath, save_count, q_cons_vf, q_prim_vf, q_T_sf)

        character(len=*), intent(in) :: step_dirpath
        integer :: save_count
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf, q_prim_vf
        type(scalar_field), intent(inout), optional :: q_T_sf

        character(LEN=len_trim(step_dirpath) + name_len) :: file_loc

        if ((prim_vars_wrt .or. (n == 0 .and. p == 0)) .and. (.not. igr)) then
            call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_vf, idwint)
            do i = 1, sys_size
                $:GPU_UPDATE(host='[q_prim_vf(i)%sf(:,:,:)]')
            end do
            ! q_prim_vf(bubxb) stores the value of nb needed in riemann solvers, so replace with true primitive value (=1._wp)
            if (qbmm) then
                q_prim_vf(bubxb)%sf = 1._wp
            end if
        end if

        if (num_dims == 1) then
            if (model_eqns == 2 .and. (.not. igr)) then
                do i = 1, sys_size
                    write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/prim.', i, '.', proc_rank, '.', save_count, '.dat'

                    open (2, FILE=trim(file_loc))
                    do j = 0, m
                        ! todo: revisit change here
                        if (((i >= adv_idx%beg) .and. (i <= adv_idx%end))) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        else
                            write (2, FMT) x_cb(j), q_prim_vf(i)%sf(j, 0, 0)
                        end if
                    end do
                    close (2)
                end do
            end if
        elseif (num_dims == 2) then
            if (prim_vars_wrt .and. (.not. igr)) then
                do i = 1, sys_size
                    write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/prim.', i, '.', proc_rank, '.', save_count, '.dat'

                    open (2, FILE=trim(file_loc))

                    do j = 0, m
                        do k = 0, n
                            if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                                .or. &
                                ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
                                ) then
                                write (2, FMT) x_cb(j), y_cb(k), q_cons_vf(i)%sf(j, k, 0)
                            else
                                write (2, FMT) x_cb(j), y_cb(k), q_prim_vf(i)%sf(j, k, 0)
                            end if
                        end do
                        write (2, *)
                    end do
                    close (2)
                end do
            end if
        elseif (num_dims == 3) then
            if (prim_vars_wrt .and. (.not. igr)) then
                do i = 1, sys_size
                    write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(step_dirpath)//'/prim.', i, '.', proc_rank, '.', save_count, '.dat'

                    open (2, FILE=trim(file_loc))

                    do j = 0, m
                        do k = 0, n
                            do l = 0, p
                                if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                                    .or. &
                                    ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
                                    .or. &
                                    ((i >= chemxb) .and. (i <= chemxe)) &
                                    ) then
                                    write (2, FMT) x_cb(j), y_cb(k), z_cb(l), q_cons_vf(i)%sf(j, k, l)
                                else
                                    write (2, FMT) x_cb(j), y_cb(k), z_cb(l), q_prim_vf(i)%sf(j, k, l)
                                end if
                            end do
                            write (2, *)
                        end do
                        write (2, *)
                    end do
                    close (2)
                end do
            end if
        end if

    end subroutine s_write_prim_variables_txt
#endif

end module m_serial_write
