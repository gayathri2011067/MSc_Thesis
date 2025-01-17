module initial_field


  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use eta_profile
  use alpha_profile
  use velocity_profile
  use seed
implicit none


double precision, dimension(nx) :: B_r, B_phi, dBr, d2Br, dBphi, d2Bphi, B_eq, rho, rho_0
double precision, dimension(nx) :: alpha_Br, alpha_Bphi, Uz_Br, Uz_Bphi,d_alpha_Br
double precision, dimension(nx) :: d2_alpha_Br,d_alpha_Bphi,d2_alpha_Bphi,d_Uz_Br
double precision, dimension(nx) :: d2_Uz_Br,d_Uz_Bphi,d2_Uz_Bphi, old_Br, old_Bphi
double precision,parameter ::  B_0=0.81873075307798182!=exp(-radius/R)
!


!



! 
contains
!NOTE: I set rho_0 to 10, after running many trials, this seemes to be correct. Do not change it for now.

    subroutine field_initialization
        call construct_alpha_profile
        call construct_velocity_profile
        !print xseed
        ! print*, "xseed=",xseed
        ! B_0=exp(-radius/R)
        ! print*, "B_0=",B_0
        B_r = 0.001*(1.0-x**2)*exp(-x**2)*B_0
        B_phi = 0.0
        B_eq = exp(-radius/R - x**2./2.)
 
        
        
        ! B_r=10**(-3)*B_0*exp(-x**2)*(1.-x**2)
        ! B_phi=0.0
        

    end subroutine field_initialization
!
  end module initial_field