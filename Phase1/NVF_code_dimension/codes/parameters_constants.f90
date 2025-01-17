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
      double precision, parameter :: s_Myr=  1.d6*365.25d0*24*3600
      double precision, parameter :: G_muG= 1.d6
  
  !***************************************************************************************************************
  
  !WITH DIMENSION
  
      double precision, parameter :: radius_dim = 4.d0 ! kpc
      double precision, parameter :: r_d_dim = 10.d0 ! kpc
      double precision, parameter :: h_d_dim = 0.35d0 ! kpc
      double precision, parameter :: eta_dim = 0.1d27*s_Gyr/(cm_kpc**2.) ! cm2/s --> kpc2/Gyr
      double precision, parameter :: h_dim = 0.5!h_d_dim*(sqrt(1+(radius_dim/r_d_dim)**2)) !kpc !disc flaring
      double precision, parameter :: t_d_dim = 0.73!h_dim**2./eta_dim ! Gyr
      double precision, parameter :: omega_0_dim = 127.*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr
      double precision, parameter :: r__omega_dim = 2. ! kpc
      double precision, parameter :: l_dim = 0.1 ! kpc
      double precision, parameter :: omega_dim = omega_0_dim*(1.+(radius_dim/r__omega_dim)**2.)**(-0.5) ! 1/Gyr
      ! double precision, parameter :: G_dim= -omega_dim ! 1/Gyr !NOTE: from paper!
      double precision, parameter :: G_dim = -45.6*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr!REVIEW: why not use the above line
      double precision, parameter :: alpha_0_dim = (l_dim**2.)*omega_dim/(h_dim) ! kpc/Gyr
      double precision, parameter :: U_0_dim = 10.*s_Gyr/km_kpc ! km/s --> kpc/Gyr
      double precision, parameter :: k_dim = 0.1*s_Gyr/km_kpc !km.kpc/s --> kpc**2/Gyr
      double precision, parameter :: R_dim = 20.!kpc
      double precision, parameter :: z_i_dim = -h_dim!kpc
      double precision, parameter :: z_f_dim = +h_dim !kpc
      double precision, parameter :: tau_c_dim =10.*s_Myr ! in gyr
   
  
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
      double precision, parameter :: U_0 = U_0_dim*h_dim/eta_dim !NOT_SURE!CHECK: old code might be wrong
      double precision, parameter :: k = k_dim/(eta_dim) !NOT_SURE
      double precision, parameter :: R = R_dim/h_dim
      double precision, parameter :: z_i = z_i_dim/h_dim
      double precision, parameter :: z_f = z_f_dim/h_dim
      double precision, parameter :: R_m_inv = 0. ! TODO_LATER: Find correct value
      double precision, parameter :: tau_c = t_d/75.


      
      
      
      !For the rest of the parameters, different trials were run, details as given below:    
      !********************************************************************************************************************************
      ! !TRIAL: 1
      ! double precision, parameter :: R_alpha = 0
      ! double precision, parameter :: R_omega = -20
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para=1.
      !********************************************************************************************************************************
      ! !TRIAL: 2
      ! double precision, parameter :: R_alpha = alpha_0
      ! double precision, parameter :: R_omega = -G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para=0.
      !********************************************************************************************************************************
      ! !TRIAL: 3
      ! double precision, parameter :: R_alpha = alpha_0
      ! double precision, parameter :: R_omega = -20
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = U_0
      ! double precision, parameter :: f_para=0.
      !********************************************************************************************************************************
      ! !TRIAL: 4
      ! double precision, parameter :: R_alpha = alpha_0
      ! double precision, parameter :: R_omega = -G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = U_0
      ! double precision, parameter :: f_para=1.
      !********************************************************************************************************************************
      ! !TRIAL: 5
      ! double precision, parameter :: R_alpha = 4
      ! double precision, parameter :: R_omega = -20
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para=1.
      !********************************************************************************************************************************
      ! !TRIAL: 6
      ! double precision, parameter :: R_alpha = 4
      ! double precision, parameter :: R_omega = -20
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para=-1.
      !********************************************************************************************************************************
      ! !TRIAL: 7
      ! double precision, parameter :: R_alpha = alpha_0
      ! double precision, parameter :: R_omega = -G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = U_0
      ! double precision, parameter :: f_para=1.
      !********************************************************************************************************************************
      ! !TRIAL: 8,9,10: different equations
      !********************************************************************************************************************************
      ! !TRIAL: 11 1A
      ! double precision, parameter :: R_alpha = 0
      ! double precision, parameter :: R_omega = -20
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para=1.
      !********************************************************************************************************************************
      ! !TRIAL: 12 1B
      ! double precision, parameter :: R_alpha = 0
      ! double precision, parameter :: R_omega = -20
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 1.16
      ! double precision, parameter :: f_para=1.
      !********************************************************************************************************************************
      ! !TRIAL: 13 1D
      ! double precision, parameter :: R_alpha = 4.
      ! double precision, parameter :: R_omega = -20.
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0.45
      ! double precision, parameter :: f_para=1.
      !********************************************************************************************************************************
      ! !TRIAL: 14: using sharanya's equations, no vishniac term
      ! double precision, parameter :: R_alpha = alpha_0
      ! double precision, parameter :: R_omega = -G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para=0
      !********************************************************************************************************************************
      !TRIAL: 15: using sharanya's equations, vishniac term, added term in dBrdt
      ! double precision, parameter :: R_alpha = alpha_0
      ! double precision, parameter :: R_omega = -G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = 1.
      ! ********************************************************************************************************************************
      !TRIAL: 16: using sharanya's equations, vishniac term, added term in dBrdt
      ! double precision, parameter :: R_alpha = alpha_0
      ! double precision, parameter :: R_omega = -G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = -1.
      ! ********************************************************************************************************************************
      !TRIAL: 17: using sharanya's equations, vishniac term, added term in dBrdt, added 3rd eqn to RK3 implicit, b0 1e-5
      ! double precision, parameter :: R_alpha = alpha_0
      ! double precision, parameter :: R_omega = -G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = -1.
      ! ********************************************************************************************************************************
      !TRIAL: 18: using sharanya's equations, vishniac term, paper values, corrected alpha_k = 0
      ! double precision, parameter :: R_alpha = alpha_0  !alpha_0
      ! double precision, parameter :: R_omega = -20.     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = -0.1
      ! ********************************************************************************************************************************
      !TRIAL: 19: using sharanya's equations, vishniac term, paper values, corrected Beq = 0.00001*exp(x**2), R_alpha = 1
      ! double precision, parameter :: R_alpha = 1.  !alpha_0
      ! double precision, parameter :: R_omega = -20.     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = -0.1
      ! ********************************************************************************************************************************
      !TRIAL: 20: using sharanya's equations, vishniac term, paper values
      ! double precision, parameter :: R_alpha = 1.  !alpha_0
      ! double precision, parameter :: R_omega = -50.     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = -0.1
      ! ********************************************************************************************************************************
      !TRIAL: 21: using sharanya's equations, vishniac term, paper values
      ! double precision, parameter :: R_alpha = 1.  !alpha_0
      ! double precision, parameter :: R_omega = -30.     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = -0.1
      ! ********************************************************************************************************************************
      !TRIAL: 22: using sharanya's equations, vishniac term, paper values
      ! double precision, parameter :: R_alpha = 1.  !alpha_0
      ! double precision, parameter :: R_omega = -40.     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = -0.1
      ! ********************************************************************************************************************************
      ! !TRIAL: 23: using sharanya's equations, vishniac term, paper values
      ! double precision, parameter :: R_alpha = 1.  !alpha_0
      ! double precision, parameter :: R_omega = -30.     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = 0.1
      ! ********************************************************************************************************************************
      !TRIAL: 24: using sharanya's equations, vishniac term, paper values
      ! double precision, parameter :: R_alpha = 1.  !alpha_0
      ! double precision, parameter :: R_omega = -40.     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = 0.1
      ! ********************************************************************************************************************************
      !TRIAL: 25: using sharanya's equations, vishniac term, paper values
      ! double precision, parameter :: R_alpha = 1.  !alpha_0
      ! double precision, parameter :: R_omega = -50.     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = 0.1
      ! ********************************************************************************************************************************
      !TRIAL: 26: using sharanya's equations, vishniac term, paper values
      ! double precision, parameter :: R_alpha = 1.  !alpha_0
      ! double precision, parameter :: R_omega = -20.     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = -0.1
      ! ********************************************************************************************************************************  
      !unist,dongsu ryu
      !takuya akahori,naoj
      ! ********************************************************************************************************************************
      !TRIAL: 27: using sharanya's equations, vishniac term, paper values
      ! double precision, parameter :: R_alpha = 1d0  !alpha_0
      ! double precision, parameter :: R_omega = -30d0     !-G
      ! double precision, parameter :: R_k = 0
      ! double precision, parameter :: R_U = 0
      ! double precision, parameter :: f_para = 0.1
      ! ********************************************************************************************************************************
      !TRIAL: 28-46
      !same as above, trials to find out the small discrepancy
      !removed the 0.01 factor added to avoid zero errors.
      ! ********************************************************************************************************************************
      !TRIAL: 28-46
      !same as above, trials to find out the small discrepancy
      !removed the 0.01 factor added to avoid zero errors.
      ! ********************************************************************************************************************************
      !CONCERN: due to some unfortunate reason, I did a git reset, thus lost all coded trials, starting again from trial 1.
      ! ********************************************************************************************************************************
      !TRIAL: 1- repeat trial 23 with Nt= 50000d0,nxphys= 101
      ! ********************************************************************************************************************************
      !TRIAL: 2- repeat trial 24
      ! ********************************************************************************************************************************
      !TRIAL: 3- repeat trial 25
      ! ********************************************************************************************************************************
      !TRIAL: 4- repeat trial 19 
      ! ********************************************************************************************************************************
      !TRIAL: trial 20
      double precision, parameter :: R_alpha = 0.  !alpha_0
      double precision, parameter :: R_omega = -20.     !-G
      double precision, parameter :: R_k = 0.
      double precision, parameter :: R_U = 0.45
      double precision, parameter :: f_para = 1.
      ! ********************************************************************************************************************************
      !TRIAL: 6- repeat trial 21
      ! ********************************************************************************************************************************



      !TRIAL: 7- catastrophic quenching test
      ! ********************************************************************************************************************************
      !TRIAL: 8- repeat trial 23,alpha_k=0
      ! ********************************************************************************************************************************
      !NOTE: with the above trials, we can now generate all plots in sharanya's notes.
      !REFER: https://github.com/gayathri2011067/MSc_Thesis/blob/main/resources/sharanya_notes.pdf
      ! ********************************************************************************************************************************
      
      
      double precision, parameter :: R_alpha_dim = R_alpha * eta_dim/h_dim
      double precision, parameter :: R_omega_dim = R_omega * eta_dim/(h_dim**2.)
      double precision, parameter :: R_k_dim = R_k * eta_dim
      double precision, parameter :: R_U_dim = R_U * eta_dim/h_dim
      double precision, parameter :: Dynamo_number = R_alpha * R_omega
    end module parameters