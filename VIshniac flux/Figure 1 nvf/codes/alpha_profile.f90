module alpha_profile

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid

  implicit none

  double precision, dimension(nx) :: alpha_k,alpha_m,alpha
  double precision, dimension(nx) :: alpha_cap,d_alpha_cap,d2alpha_cap, alpha_cap2

  contains
    subroutine construct_alpha_profile
        alpha_cap = sin(pi*x)
        alpha_k = 0d0!alpha_0*alpha_cap  !NOTE: kinetic helicity set to 0 now
        alpha_m = 0d0
        d_alpha_cap = cos(pi*x)
        alpha = alpha_m !+alpha_k
      end subroutine construct_alpha_profile

end module alpha_profile