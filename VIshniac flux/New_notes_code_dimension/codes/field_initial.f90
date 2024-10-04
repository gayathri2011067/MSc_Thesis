module initial_field


  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use eta_profile
  use alpha_profile
  use velocity_profile
  use spatial_derivatives
!
  implicit none
!
  double precision, dimension(nx) :: B_r, B_phi, dBr, d2Br, dBphi, d2Bphi, B_eq
  double precision, dimension(nx) :: alpha_Br, alpha_Bphi, Uz_Br, Uz_Bphi,d_alpha_Br
  double precision, dimension(nx) :: d2_alpha_Br,d_alpha_Bphi,d2_alpha_Bphi,d_Uz_Br
  double precision, dimension(nx) :: d2_Uz_Br,d_Uz_Bphi,d2_Uz_Bphi, old_Br, old_Bphi
  double precision, dimension(nx) :: T_torr,phi, Fr,F_phi,Fz, B0,rho,rho_0,B_0

!CONCERN:

! 
contains


    subroutine field_initialization
        call construct_alpha_profile
        call construct_velocity_profile

       
        rho_0 = 1.0
        rho = rho_0*exp(-x**2)

        B_eq = sqrt(4.0*pi*rho)*U_0_dim
        B_0 = sqrt(4.0*pi*rho_0)*U_0_dim
      
        T_torr=10.**(-2)*(1d0-radius)*exp(-radius - x**2.)
        phi=10.**(-2)*(1d0-radius)*exp(-radius - x**2.)
        Fr=0.
        F_phi=0.
        Fz=0.


        




    end subroutine field_initialization
!
  end module initial_field