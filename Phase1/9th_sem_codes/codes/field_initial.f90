module initial_field


  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use eta_profile
  use alpha_profile
  use velocity_profile
  use seed
  use spatial_derivatives
implicit none


double precision, dimension(nx) :: B_r, B_phi, dBr, d2Br, dBphi, d2Bphi, B_eq
double precision, dimension(nx) :: Fr, Fphi, Er, Ephi, Tee, Phii, dphii, d2phii,small_u!,small_u_0
double precision, dimension(nx) :: alpha_Br, alpha_Bphi, Uz_Br, Uz_Bphi,d_alpha_Br, B_0,eta
double precision, dimension(nx) :: d2_alpha_Br,d_alpha_Bphi,d2_alpha_Bphi,d_Uz_Br
double precision, dimension(nx) :: d2_Uz_Br,d_Uz_Bphi,d2_Uz_Bphi, old_Br, old_Bphi
! double precision, dimension(nx) :: rho,rho_0
double precision :: Bseed
! double precision, dimension(nx) :: tau
double precision, dimension(nx) :: l


contains

    subroutine field_initialization
        call construct_alpha_profile
        call construct_velocity_profile
        call init_random_seed
        call spatial_derivative(phii, 6, dphii, d2phii)

        !NOTE: to make rho stratified, uncomment the following line, and set B_eq and B_o constant.
        !also, uncomment rho in parameters_constants.f90, and uncomment line 21 here.
        ! rho_0     = 1.
        ! rho       = rho_0*exp(-x**2/2.)
        ! B_eq      = 4.*pi*rho_0*small_u_0**2
        ! B_0       = 4.*pi*rho_0*small_u_0**2

  
        Bseed     = 100

        small_u   = small_u_0*exp(x**2/2.)
        ! small_u   = small_u_0
        ! small_u  = small_u_0*(x)
        l = tau*small_u
        ! tau = l/small_u


        eta       = (1.*tau/(3.))*small_u**2
        !new expressionfrom ss21
        ! small_u   = small_u_0*(abs(x**(1./2.)))
        ! eta       = (1.*tau/(3.))*small_u**2


        
        B_eq      = 4.*pi*rho*small_u**2
        B_0       = 4.*pi*rho*small_u_0**2

        B_r       = xseed*exp(-x**2)*(1.-x**2)*Bseed
        B_phi     = xseed*exp(-x**2)*(1.-x**2)*Bseed




    end subroutine field_initialization
!
  end module initial_field