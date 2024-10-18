module parameters
  !this module contains the all the parameter values used for simulation
  
    implicit none
  
  !CONSTANT VALUES
  
      double precision, parameter :: pi= 3.14156295358
      double precision, parameter :: g_mp=   1.67262158d-24
      double precision, parameter :: cm_kpc= 3.08567758d21
      double precision, parameter :: km_kpc= 3.08567758d16
      double precision, parameter :: cm_km=  1.d5
      double precision, parameter :: s_Gyr=  1.d9*365.25d0*24*3600
      double precision, parameter :: G_muG= 1.d6
  
  !***************************************************************************************************************
  
  !WITH DIMENSION
  
      double precision, parameter :: radius_dim = 4.d0 ! kpc
      double precision, parameter :: r_d_dim = 10.d0 ! kpc
      double precision, parameter :: h_d_dim = 0.5d0 ! kpc
      double precision, parameter :: eta_dim = 0.1d27*s_Gyr/(cm_kpc**2.) ! cm2/s --> kpc2/Gyr
      double precision, parameter :: h_dim = 0.5!h_d_dim*(sqrt(1+(radius_dim/r_d_dim)**2)) !kpc !disc flaring
      double precision, parameter :: t_d_dim = 0.73 ! Gyr
      double precision, parameter :: omega_0_dim = 127.*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr
      double precision, parameter :: r__omega_dim = 2. ! kpc
      double precision, parameter :: l_dim = 0.1 ! kpc
      double precision, parameter :: omega_dim = omega_0_dim*(1.+(radius_dim/r__omega_dim)**2.)**(-0.5) ! 1/Gyr
      ! double precision, parameter :: G_dim= -omega_dim ! 1/Gyr !NOTE: from paper!
      double precision, parameter :: G_dim = -45.6*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr!REVIEW: why not use the above line
      double precision, parameter :: alpha_0_dim = 1.50 ! kpc/Gyr
      double precision, parameter :: small_u_dim = 10*s_Gyr/km_kpc ! km/s
      double precision, parameter :: k_dim = 0.1*s_Gyr/km_kpc !km.kpc/s --> kpc**2/Gyr
      double precision, parameter :: R_dim = 20.!kpc
      double precision, parameter :: z_i_dim = -h_dim!kpc
      double precision, parameter :: z_f_dim = +h_dim !kpc
      double precision, parameter :: tau_c_dim = 10.*0.001 !Gyr ! TODO_LATER: Find correct value
      double precision, parameter :: B_0_dim = 8.2/G_muG
      !B_0_dim=8.2e-6*u.G,B_0=1
  !***************************************************************************************************************
  
  !DIMENSIONLESS
  
      double precision, parameter :: radius = radius_dim/h_dim
      double precision, parameter :: r_d = r_d_dim/h_dim
      double precision, parameter :: h_d = h_d_dim/h_dim
      double precision, parameter :: eta = eta_dim/(eta_dim) !NOTE: use profile for thick disc
      double precision, parameter :: h = h_dim/h_dim
      double precision, parameter :: t_d = t_d_dim/(h_dim**2./eta_dim)
      double precision, parameter :: omega_0 = omega_0_dim*h_dim/(eta_dim) 
      double precision, parameter :: r__omega = r__omega_dim/h_dim
      double precision, parameter :: l = l_dim/h_dim
      double precision, parameter :: omega = omega_dim*h_dim/(eta_dim)
      double precision, parameter :: G = G_dim*(h_dim**2/(eta_dim))
      double precision, parameter :: alpha_0 = alpha_0_dim*h_dim/eta_dim
      double precision, parameter :: small_u = small_u_dim*h_dim/eta_dim !NOT_SURE!CHECK: old code might be wrong
      double precision, parameter :: k = k_dim/(eta_dim) !NOT_SURE
      double precision, parameter :: R = R_dim/h_dim
      double precision, parameter :: z_i = z_i_dim/h_dim
      double precision, parameter :: z_f = z_f_dim/h_dim
      double precision, parameter :: R_m_inv = 0. ! TODO_LATER: Find correct value
      double precision, parameter :: tau_c = tau_c_dim/(h_dim**2./eta_dim)
      double precision, parameter :: B_0 = B_0_dim/(B_0_dim)
 

! double precision, parameter :: B_0 = 1.0


!For the rest of the parameters, different trials were run, details as given below:    
!********************************************************************************************************************************
!MTA! figures>old folder
      
      !TRIAL:1--->High resolution,1A
      !TRIAL:2--->Low resolution,1A
      !TRIAL:3--->f=0,1G
      !TRIAL:4--->2
      !TRIAL:5--->1D
      !TRIAL:6--->1A
!********************************************************************************************************************************
!MTA!figures

!TRIAL:1____1A
      double precision, parameter :: R_alpha = 1.5  
      double precision, parameter :: R_omega = -20.
      double precision, parameter :: R_k = 0.
      double precision, parameter :: R_U = 0.45
      double precision, parameter :: f_para = 1.
      double precision, parameter :: c_tau = 1.
!********************************************************************************************************************************
!TRIAL:2____1B
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 1.16
      ! double precision, parameter :: f_para = 1.
      ! double precision, parameter :: c_tau = 1.
!********************************************************************************************************************************
!TRIAL:3____1C
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -30.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 1.
      ! double precision, parameter :: c_tau = 1.
!********************************************************************************************************************************
!TRIAL:4____1D
      ! double precision, parameter :: R_alpha = 4.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 1.
      ! double precision, parameter :: c_tau = 1.
!********************************************************************************************************************************
!TRIAL:5____1E
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 1.
      ! double precision, parameter :: c_tau = 2.
!********************************************************************************************************************************
!TRIAL:6____1F
      ! double precision, parameter :: R_alpha = 0.  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.3
      ! double precision, parameter :: R_U = 0.
      ! double precision, parameter :: f_para = 1.
      ! double precision, parameter :: c_tau = 1.
!********************************************************************************************************************************
!TRIAL:7____1G
      ! double precision, parameter :: R_alpha = 1.5  
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 0.
      ! double precision, parameter :: c_tau = 1.
!********************************************************************************************************************************
!TRIAL:8____1H
      ! double precision, parameter :: R_alpha = 0.
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0.
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para = 0.34
      ! double precision, parameter :: c_tau = 1.
!********************************************************************************************************************************

!python3 codes/compare_plots.py 1 "1A" 2 "1B" 3 "1C" 4 "1D" 5 "1E" 6 "1F" 7 "1G" 8 "1H" 













!********************************************************************************************************************************
      double precision, parameter :: tau = c_tau*tau_c
      double precision, parameter :: Dynamo_number = R_alpha * R_omega

    end module parameters