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


double precision, dimension(nx) :: B_r, B_phi, dBr, d2Br, dBphi, d2Bphi, B_eq,small_u
double precision, dimension(nx) :: Fr, Fphi, Er, Ephi, Tee, Phii, dphii, d2phii
double precision, dimension(nx) :: alpha_Br, alpha_Bphi, Uz_Br, Uz_Bphi,d_alpha_Br, B_0
double precision, dimension(nx) :: d2_alpha_Br,d_alpha_Bphi,d2_alpha_Bphi,d_Uz_Br
double precision, dimension(nx) :: d2_Uz_Br,d_Uz_Bphi,d2_Uz_Bphi, old_Br, old_Bphi
double precision :: Bseed
!


!



! 
contains
!NOTE: I set rho_0 to 10, after running many trials, this seemes to be correct. Do not change it for now.

    subroutine field_initialization

        call construct_alpha_profile
        call construct_velocity_profile
        call init_random_seed
        call spatial_derivative(phii, 6, dphii, d2phii)

        Bseed = 0.00001
        small_u = u_zero*exp(x**2/2.)
        B_eq = 4.*pi*rho*small_u**2
        B_0  = 4.*pi*rho*u_zero**2

        B_r = xseed*exp(-x**2)*(1.-x**2)*Bseed
        B_phi = xseed*exp(-x**2)*(1.-x**2)*Bseed

        ! phii = xseed*exp((-x**2)/2)*(1./sqrt(2*pi))*Bseed
        ! Tee = xseed*exp((-x**2)/2)*(1./sqrt(2*pi))*Bseed

        ! B_phi = Tee/(x+0.01)
        ! B_r = -dphii/(x+0.01)

        Fr=0.0
        Fphi=0.0
        Er=0.0
        Ephi=0.0




    end subroutine field_initialization
!
end module initial_field