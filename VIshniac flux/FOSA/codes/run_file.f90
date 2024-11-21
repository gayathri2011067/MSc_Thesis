program run_all

    use alpha_profile
    use velocity_profile
    use initial_field
    use seed
    use parameters
    use time_grid
    use physical_grid
    use make_a_grid
    use eta_profile
    use omega_profile
    use equations
    use spatial_derivatives
    use timestepping
  
    implicit none

    integer :: kk, j
    character(len=30) :: data_path, filename, xfile, omegafile, alphafile, Br_ini_file,&
    B_r_final_file, B_phi_ini_file, B_phi_final_file, time_file, alpham_final_file,turb_vel_file,&
    vishniac_term_file,alp_tot_file, alp_k_file

    ! Call subroutine to construct the physical grid
    call construct_grid

    ! Call subroutine to construct eta profile
    call construct_eta_profile
    call construct_omega_profile
    call construct_alpha_profile
    call field_initialization


    ! Define the output file name
    data_path = '../run_files/'
    filename =  'eta_fz_values.txt'
    xfile=  'z_values.txt'
    omegafile=  'omega_values.txt'
    alphafile=  'alpha_values.txt'
    Br_ini_file=  'Br_ini.txt'
    B_phi_ini_file=  'B_phi_ini.txt'
    B_r_final_file=  'Br_final.txt'
    B_phi_final_file=  'B_phi_final.txt'
    time_file=  'time.txt'
    alpham_final_file=  'alpham_final.txt'
    turb_vel_file=  'turb_vel.txt'
    vishniac_term_file=  'vishniac_term.txt'
    alp_tot_file = 'alp_tot.txt'
    alp_k_file = 'alp_k.txt'
    ! print *, "Computational time in Gyr = ", total_t

    ! Open the file for writing
    open(unit=10, file=trim(data_path) // filename)
    open(unit=17, file=trim(data_path) // xfile)
    open(unit=19, file=trim(data_path) // alphafile)
    open(unit=20, file=trim(data_path) // Br_ini_file)
    open(unit=21, file=trim(data_path) // B_phi_ini_file)
    open(unit=22, file=trim(data_path) // B_r_final_file)
    open(unit=23, file=trim(data_path) // B_phi_final_file)
    open(unit=24, file=trim(data_path) // time_file)
    open(unit=25, file=trim(data_path) // alpham_final_file)
    open(unit=26, file=trim(data_path) // turb_vel_file)
    open(unit=27, file=trim(data_path) // vishniac_term_file)
    open(unit=28, file=trim(data_path) // alp_tot_file)
    open(unit=29, file=trim(data_path) // alp_k_file)

    ! open(unit=10, file=filename)
    ! open(unit=17, file=xfile)
    ! open(unit=19, file=alphafile)
    ! open(unit=20, file=Br_ini_file)
    ! open(unit=21, file=B_phi_ini_file)
    ! open(unit=22, file=B_r_final_file)
    ! open(unit=23, file=B_phi_final_file)
    ! open(unit=24, file=time_file)
    ! Write the values to the file
    do i = 1, nx
        write(10, '(F12.8)') eta_fz(i)
        write(17, '(F12.8)') x(i) 
        write(19, '(F12.8)') alpha_k(i)
        write(20, '(F12.8)') B_r(i)
        write(21, '(F12.8)') B_phi(i)
        write(26, '(F12.8)') small_u(i)

    end do

    ! Close the file
    close(10)
    close(17)
    close(19)
    close(20)
    close(21)
    close(26)

    ! print *, 'z values have been saved to ', xfile
    ! print *, 'omega values have been saved to ', omegafile
    ! print *, 'alpha values have been saved to ', alphafile
    ! print *, 'Br_ini values have been saved to ', Br_ini_file
    print *, 't_d ', t_d_dim
    

  call field_initialization

  ! ************************************************************************
  !RK4 without chain rule for derivatives
  ! do kk = 1, n1 ! for n1 iterations
  !   do j = 1, n2 ! for n2 time steps
  !     call RK4_new
  !     ! print*, 't=', t
  !     ! print*, 'B_r=', B_r

  !   end do
  !   ! print*, 'B_r=', B_r
  !   write (22, *) B_r
  !   write (23, *) B_phi
  !   write (24, *) t
  !   write (25, *) alpha_m
  ! end do
  ! ************************************************************************
  ! Forward differencing for time-stepping
  ! do kk = 1, n1 ! for n1 iterations
  !   do j = 1, n2 ! for n2 time steps
  !     call forward_difference
  !     ! print*, 't=', t
  !     ! print*, 'B_r=', B_r
  !     ! print*, 'B_phi=', B_phi
  !   end do
  !   ! print*, 'B_r=', B_r
  !   write (22, *) B_r
  !   write (23, *) B_phi
  !   write (24, *) t
  ! end do
  ! ************************************************************************
  ! Backward differencing for time-stepping
  ! do kk = 1, n1 ! for n1 iterations
  !   do j = 1, n2 ! for n2 time steps
  !     call backward_difference
  !     ! print*, 't=', t
  !     ! print*, 'B_r=', B_r
  !     ! print*, 'B_phi=', B_phi
  !   end do
  !   ! print*, 'B_r=', B_r
  !   write (22, *) B_r
  !   write (23, *) B_phi
  !   write (24, *) t
  ! end do
  ! ************************************************************************
  ! Central differencing for time-stepping
  ! old_Br = B_r
  ! old_Bphi = B_phi
  ! call RK4_new    !NOTE: using RK4 to get the first time step
  ! first = 0.  !resetting
  ! t=0.
  ! do kk = 1, n1 ! for n1 iterations
  !   do j = 1, n2 ! for n2 time steps
  !     call central_difference
  !     ! print*, 't=', t
  !     ! print*, 'B_r=', B_r
  !     ! print*, 'B_phi=', B_phi
  !   end do
  !   ! print*, 'B_r=', B_r
  !   write (22, *) B_r
  !   write (23, *) B_phi
  !   write (24, *) t
  ! end do
  ! ************************************************************************
  !RK3 implicit
  do kk = 1, n1 ! for n1 iterations
    ! print*, 'iteration ', kk -1, 'completed'
    do j = 1, n2 ! for n2 time steps
      call RK3_implicit
      ! print*, 't=', t
      ! print*, 'B_r=', B_r

    end do
    ! print*, 'B_r=', B_r
    write (22, *) B_r/B_0
    write (23, *) B_phi/B_0
    write (24, *) t*t_d_dim
    write (25, *) alpha_m
    write (28, *) alpha
    write (29, *) alpha_k
    ! write (27, *) vishniac_term


  end do
  ! ************************************************************************



    close(22)
    close(23)
    close(24)
    close(25)
    close(28)
    close(29)
    ! close(27)
  ! print*, 'n2=', n2
  ! print*, 'dt=', dt
  ! print*, 'Nt=', Nt
  ! print*, 'n1=', n1
  ! print*, 'total_t=', total_t
  print*, 't=', t
  ! print*, 'first=', first
    print *, 'File sucessfully run'

end program run_all
