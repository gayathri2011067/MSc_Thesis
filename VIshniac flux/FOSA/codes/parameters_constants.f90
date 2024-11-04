module parameters
  !this module contains the all the parameter values used for simulation

      implicit none
  
  !CONSTANT VALUES ♥ for conversions, all the dimensional values are in kpc and gyr terms
  
      double precision, parameter :: pi= 3.14156295358
      double precision, parameter :: g_mp=   1.67262158d-24
      double precision, parameter :: cm_kpc= 3.08567758d21
      double precision, parameter :: km_kpc= 3.08567758d16
      double precision, parameter :: cm_km=  1.d5
      double precision, parameter :: s_Gyr=  1.d9*365.25d0*24*3600
      double precision, parameter :: G_muG= 1.d6
  
  !***************************************************************************************************************
  
  !WITH DIMENSION ♥

      double precision, parameter :: small_u_0_dim = 10.*km_kpc/s_Gyr ! km/s.kpc --> 1/Gyr
      ! double precision, parameter :: small_u_0_dim = 45.*km_kpc/s_Gyr ! km/s.kpc --> 1/Gyr
!above expression for new u from ss21
      double precision, parameter :: radius_dim = 4.d0 ! kpc
      double precision, parameter :: r_d_dim = 10.d0 ! kpc
      double precision, parameter :: h_d_dim = 0.5d0 ! kpc
      double precision, parameter :: h_dim = 0.5!h_d_dim*(sqrt(1+(radius_dim/r_d_dim)**2)) !kpc !disc flaring
      double precision, parameter :: omega_0_dim = 127.*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr
      double precision, parameter :: r__omega_dim = 2. ! kpc
      double precision, parameter :: l_dim = 0.1 ! kpc
      double precision, parameter :: omega_dim = omega_0_dim*(1.+(radius_dim/r__omega_dim)**2.)**(-0.5) ! 1/Gyr
      double precision, parameter :: G_dim = -45.6*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr
      double precision, parameter :: alpha_0_dim = 1.50 ! kpc/Gyr
      double precision, parameter :: k_dim = 0.3*s_Gyr/km_kpc !km.kpc/s --> kpc**2/Gyr
      double precision, parameter :: R_dim = 20.!kpc
      double precision, parameter :: z_i_dim = -h_dim!kpc
      double precision, parameter :: z_f_dim = +h_dim !kpc
      double precision, parameter :: tau_dim = 10.*0.001 !Gyr ! TODO_LATER: Find correct value
      double precision, parameter :: B_0_dim = 0.82/G_muG !muG to G 
      double precision, parameter :: t_d_dim = h_dim**2/((tau_dim*small_u_0_dim**2)/3.)! Gyr
  !***************************************************************************************************************
  
  !DIMENSIONLESS ♥

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !NOTE: For models A,B and D, rho is constant and small U is stratified, which makes eta z dependent.
      !Thus, the F_NV term changes from NVF_pdf, as given
  
      double precision, parameter :: radius = radius_dim/h_dim
      double precision, parameter :: r_d = r_d_dim/h_dim
      double precision, parameter :: h_d = h_d_dim/h_dim
      double precision, parameter :: h = h_dim/h_dim
      double precision, parameter :: t_d = t_d_dim/(t_d_dim)
      double precision, parameter :: omega_0 = omega_0_dim*h_dim/(h_dim**2/t_d_dim) 
      double precision, parameter :: r__omega = r__omega_dim/h_dim
      ! double precision, parameter :: l = l_dim/h_dim
      double precision, parameter :: omega = omega_dim*h_dim/(h_dim**2/t_d_dim)
      double precision, parameter :: G = G_dim*(h_dim**2/(h_dim**2/t_d_dim))
      double precision, parameter :: alpha_0 = alpha_0_dim*h_dim/(h_dim**2/t_d_dim)
      double precision, parameter :: rho = 1. !NOTE:!comment for rho stratified
      double precision, parameter :: k = k_dim/(h_dim**2/t_d_dim) 
      double precision, parameter :: R = R_dim/h_dim
      double precision, parameter :: z_i = z_i_dim/h_dim
      double precision, parameter :: z_f = z_f_dim/h_dim
      double precision, parameter :: R_m_inv = 0. 
      double precision, parameter :: tau = tau_dim/t_d_dim
      double precision, parameter :: small_u_0 = small_u_0_dim*h_dim/(h_dim**2/t_d_dim)
 


!For the rest of the parameters, different trials were run, details as given below:    
!********************************************************************************************************************************

!********************************************************************************************************************************
!FOSA!figures
!Trials to recreate the figures from the 2017 new notes
! TRIAL:1____A
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.
      ! double precision, parameter :: f_para = 1.
!********************************************************************************************************************************
!TRIAL:2____1B
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 1.
!********************************************************************************************************************************
!TRIAL:3____1C !CONCERN:!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 1.
!********************************************************************************************************************************
!TRIAL:4____1D
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 0.1
!********************************************************************************************************************************
!TRIAL:10_____New expression for U with l constant-A
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 1.
      ! double precision, parameter :: R_U = 0.
      ! double precision, parameter :: f_para = 1.
!********************************************************************************************************************************
!TRIAL:11_____New expression for U-B
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 1.
!********************************************************************************************************************************
!TRIAL:12_____New expression for U-D
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 0.1
!********************************************************************************************************************************
! 13,14,15------>New expressions for U with tau constant, A,B,D
!********************************************************************************************************************************


!NOTE: The following trials are for recreating 2014 toolbox paper results
! ********************************************************************************************************************************
!TRIAL:21 
double precision, parameter :: R_alpha = 1.71  !alpha_0
double precision, parameter :: R_omega = -19.8     !-G
double precision, parameter :: R_k = 0.0
double precision, parameter :: R_U = 1.0
double precision, parameter :: f_para = 0.
! ********************************************************************************************************************************
!TRIAL:22
! double precision, parameter :: R_alpha = 1.71  !alpha_0
! double precision, parameter :: R_omega = -19.8     !-G
! double precision, parameter :: R_k = 0.3
! double precision, parameter :: R_U = 0.0
! double precision, parameter :: f_para = 0.
! ********************************************************************************************************************************
!TRIAL:23
! double precision, parameter :: R_alpha = 1.71  !alpha_0
! double precision, parameter :: R_omega = -19.8     !-G
! double precision, parameter :: R_k = 0.3
! double precision, parameter :: R_U = 1.0
! double precision, parameter :: f_para = 0.
! ********************************************************************************************************************************
!TRIAL:24
! double precision, parameter :: R_alpha = 1.71  !alpha_0
! double precision, parameter :: R_omega = -19.8     !-G
! double precision, parameter :: R_k = 0.0
! double precision, parameter :: R_U = 0.0
! double precision, parameter :: f_para = 0.


!python3 codes/compare_plots.py 1 "model A" 2 "model B" 3 "model C" 4 "model D" 













!********************************************************************************************************************************
      double precision, parameter :: Dynamo_number = R_alpha * R_omega

    end module parameters