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
double precision :: Bseed
double precision, dimension(nx) :: l


double precision, dimension(nx) :: rho,rho_0


contains

    subroutine field_initialization
        call construct_alpha_profile
        call construct_velocity_profile
        call init_random_seed
        call spatial_derivative(phii, 6, dphii, d2phii)
        Bseed     = 100.  !CONCERN: Bseed is 100 because we are getting 10e-3 data for this number

        !NOTE: to make ONLY rho stratified, uncomment the following lines
        ! l         = l_dim/h_dim
        ! small_u   = l/tau
        ! rho_0     = 1.
        ! rho       = rho_0*exp(-x**2/2.)
        ! B_eq      = 4.*pi*rho*small_u**2
        ! B_0       = 4.*pi*rho*small_u_0**2
        ! eta       = (1.*tau/(3.))*small_u**2





        !NOTE: to make  ONLY u stratified, uncomment the following lines----!CONCERN: Old gaussian u profile

        rho_0     = 1.
        rho       = 1.
        small_u   = small_u_0*exp(x**2/2.)
        l = tau*small_u
        eta       = (1.*tau/(3.))*small_u**2
        B_eq      = 4.*pi*rho*small_u**2
        B_0       = 4.*pi*rho_0*small_u_0**2





        !NOTE: to make rho and u stratified, uncomment the following lines, here, B_eq is a constant----!CONCERN: Old gaussian u profile
        ! rho_0     = 1.
        ! rho       = rho_0*exp(-x**2/2.)
        ! small_u   = small_u_0*exp(x**2/2.)
        ! l = tau*small_u
        ! eta       = (1.*tau/(3.))*small_u**2
        ! B_eq      = 4.*pi*rho_0*small_u_0**2
        ! B_0       = 4.*pi*rho_0*small_u_0**2


        !NOTE: to make  ONLY u stratified, uncomment the following lines----!CONCERN: New power law u profile

        ! rho_0     = 1.
        ! rho       = 1.
        ! small_u   = 45.*km_kpc/s_Gyr * abs(x)**(0.5)  !CONCERN: small_U_0 is 45 not 10
        ! l = tau*small_u
        ! eta       = (1.*tau/(3.))*small_u**2
        ! B_eq      = 4.*pi*rho*small_u**2
        ! B_0       = 4.*pi*rho*small_u_0**2


        !NOTE: to make rho and u stratified, uncomment the following lines, here, B_eq is a constant----!CONCERN: New power law u profile

        ! rho_0     = 1.
        ! rho       = rho_0*exp(-x**2/2.)
        ! small_u   = 45.*km_kpc/s_Gyr * abs(x)**(0.5)  !CONCERN: small_U_0 is 45 not 10
        ! l = tau*small_u
        ! eta       = (1.*tau/(3.))*small_u**2
        ! B_eq      = 4.*pi*rho_0*small_u_0**2
        ! B_0       = 4.*pi*rho_0*small_u_0**2






  

        B_r       = xseed*exp(-x**2)*(1.-x**2)*Bseed
        B_phi     = xseed*exp(-x**2)*(1.-x**2)*Bseed




    end subroutine field_initialization
!
  end module initial_field