module eta_profile

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
!
  implicit none
!
  double precision, dimension(nx) :: eta_dyn_dim, eta_dyn

  

! 
contains
    subroutine construct_eta_profile

    ! eta_dyn_dim=(1.0/3.0 * tau_dim )*u_zero_dim**2 *exp((x*h_dim)**2)
    ! eta_dyn=eta_dyn_dim/(h_dim**2/t_d_dim)
      eta_dyn=(1.0/3.0 * tau )*u_zero**2 *exp(x**2)
      


    end subroutine construct_eta_profile
!
  end module eta_profile