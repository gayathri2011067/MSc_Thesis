module alpha_profile

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  ! use omega_profile
  ! use spatial_derivatives
!
  implicit none
!
 
  double precision, dimension(nx) :: alpha_k,alpha_m,alpha
  double precision, dimension(nx) :: alpha_cap,d_alpha_cap,d2alpha_cap, alpha_cap2

 

! 
contains
    subroutine construct_alpha_profile
        character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
        alpha_cap = sin(pi*x)
        d_alpha_cap = cos(pi*x)
        alpha_k = R_alpha*alpha_cap  !NOTE: kinetic helicity set to 0 now
        alpha_m = 0.
        alpha = alpha_m + alpha_k
  

    end subroutine construct_alpha_profile
!
  end module alpha_profile