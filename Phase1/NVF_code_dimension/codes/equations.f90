module equations
    use parameters
    use eta_profile
    use velocity_profile
    use alpha_profile
    use omega_profile
    use initial_field
    use time_grid
    use physical_grid
    use make_a_grid

    implicit none
    double precision, dimension(nx) :: dBrdt, dBphidt, dalpdt,dTdt,dphi_dt,d_alpha_m_dt
    double precision, dimension(nx) :: dFr_dt, dFphi_dt, dFz_dt
     !,d_alpha_Bphi, d_alpha_Br, d_Uz_Bphi, d2_Uz_Br
    !   NOTE: Assumed to be predefined: B_r, B_phi,dBr, d2Br, dBphi, d2Bphi, d_alpha
    contains
    ! subroutine diff_equations_split(B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy)
    !   double precision, intent(in), dimension(nx) :: B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy
    !   dBrdt = R_alpha*B_phi_dummy*d_alpha_cap + R_alpha*alpha_cap2*dBphi_dummy &
    !   + d2Br_dummy - R_U*U_z_cap*dBr_dummy - R_U*B_r_dummy*d_U_z_cap
    !   dBphidt = R_omega*B_r_dummy + R_alpha*B_r_dummy*d_alpha_cap + R_alpha*alpha_cap2*dBr_dummy &
    !   + d2Bphi_dummy - R_U*U_z_cap*dBphi_dummy - R_U*B_phi_dummy*d_U_z_cap

    ! end subroutine diff_equations_split
  
!     subroutine diff_equations_no_split(B_r_dummy, B_phi_dummy)
!       double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy
!       double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
!       double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
!       double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
!       double precision,  dimension(nx) :: Uz_Bphi, Uz_Br
!       double precision,  dimension(nx) :: alpha_Bphi, alpha_Br

!       character(len=30) :: ghost_zone_type = 'anti-symmetric'
!       call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
!       call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)

!       alpha_cap2 = alpha_cap / (1.0 + (B_r_dummy**2 + B_phi_dummy**2) / B_eq**2)
!       alpha_Bphi = alpha_cap2 * B_phi_dummy
!       alpha_Br = alpha_cap2 * B_r_dummy
!       Uz_Br = B_r_dummy * U_z_cap
!       Uz_Bphi = B_phi_dummy * U_z_cap

!       call spatial_derivative(B_r_dummy, 2, dBr, d2Br)
!       call spatial_derivative(B_phi_dummy, 2, dBphi, d2Bphi)
!       call spatial_derivative(alpha_Br, 2, d_alpha_Br, d2_alpha_Br)
!       call spatial_derivative(alpha_Bphi, 2, d_alpha_Bphi, d2_alpha_Bphi)
!       call spatial_derivative(Uz_Br, 2, d_Uz_Br, d2_Uz_Br)
!       call spatial_derivative(Uz_Bphi, 2, d_Uz_Bphi, d2_Uz_Bphi)


!       dBrdt = - R_alpha*d_alpha_Bphi &
!       + d2Br - R_U*d_Uz_Br
!       dBphidt = R_omega*B_r_dummy + R_alpha*d_alpha_Br &
!       + d2Bphi - R_U*d_Uz_Bphi

!     end subroutine diff_equations_no_split
  
!     subroutine dyn_quenching_nvf(B_r_dummy, B_phi_dummy, alpha_m_dummy)
!       double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
!       double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
!       double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi, d2_alpha_Br, d2_alpha_Bphi
!       double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
!       double precision,  dimension(nx) :: Uz_Bphi, Uz_Br,der_u, d_sq_u, d2_Uz_Br,d2_Uz_Bphi
!       double precision,  dimension(nx) :: alpha_Bphi, alpha_Br
!       double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
!       double precision, dimension(nx) :: alpha_total
!       double precision, dimension(nx) :: vishniac_term


!       character(len=30) :: ghost_zone_type = 'anti-symmetric'
!       call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
!       call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
!       call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type)
!       alpha_total = R_alpha*alpha_k + alpha_m_dummy
!       ! alpha_cap2 = alpha_cap / (1.0 + (B_r_dummy**2 + B_phi_dummy**2) / B_eq**2)
!       alpha_Bphi = alpha_total * B_phi_dummy
!       alpha_Br = alpha_total * B_r_dummy
!       Uz_Br = B_r_dummy * U_z_cap
!       Uz_Bphi = B_phi_dummy * U_z_cap


!       call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
!       call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
!       call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
!       call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
!       call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
!       call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)
!       call spatial_derivative(alpha_m_dummy, 6, d_alpha_m, d2_alpha_m)
!       call spatial_derivative(U_z_cap, 6, der_u, d_sq_u)


!       ! vishniac_term = x* f_para * R_omega*(2*l/3*h)**2 
!       ! vishniac_term = -f_para * R_omega *x

      
!       dBrdt = - d_alpha_Bphi &
!       + d2Br - R_U*d_Uz_Br

!       dBphidt = R_omega*B_r_dummy + d_alpha_Br &
!       + d2Bphi - R_U*d_Uz_Bphi

!       ! dalpdt = -2./3.*t_d(alpha_total*(B_r_dummy**2. + B_phi_dummy**2.)/B_eq**2. &
!       ! (B_phi_dummy*dBphi_dummy-B_r_dummy*dBr_dummy)/B_eq**2. + vishniac_term) &
!       ! +R_k*d2_alpha_m - R_U*U_z_cap*d_alpha_m-R_U*alpha_m_dummy*d_U_z_cap

!       ! dalpdt = -2.0 / 3.0 * t_d * (alpha_total * (B_r_dummy**2.0 + B_phi_dummy**2.0) / B_eq**2.0 +&
!       !      (B_phi_dummy * dBphi - B_r_dummy * dBr) / B_eq**2.0 +vishniac_term)&
!       !      + R_k * d2_alpha_m - R_U * U_z_cap * d_alpha_m - R_U * alpha_m_dummy * d_U_z_cap

!       dalpdt =    (-2/((3*t_d)*(B_eq**2)))*((alpha_total)*(B_r_dummy**2+B_phi_dummy**2)&
!       - (B_phi_dummy*dbr-B_r_dummy*dBphi))&
!       - R_u*alpha_m_dummy*der_u-R_u*U_z*d_alpha_m + R_k*d2_alpha_m

!     end subroutine dyn_quenching_nvf
  
  
!     subroutine sharanyas_notes_eqxn(B_r_dummy, B_phi_dummy, alpha_m_dummy)
!       double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
!       double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
!       double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
!       double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
!       double precision,  dimension(nx) :: Uz_Bphi, Uz_Br
!       double precision,  dimension(nx) :: alpha_Bphi, alpha_Br
!       double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
!       double precision, dimension(nx) :: alpha_total
!       double precision, dimension(nx) :: vishniac_term
!       double precision, dimension(nx) :: der_u, d_sq_u
!       double precision, dimension(nx) :: d2_alpha_Br, d2_alpha_Bphi
!       double precision, dimension(nx) :: d2_Uz_Br, d2_Uz_Bphi

!       character(len=30) :: ghost_zone_type = 'anti-symmetric'
!       character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'

!       call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
!       call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
!       call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type2)
!       alpha_total = alpha_m_dummy!R_alpha*alpha_k + 
!       ! alpha_cap2 = alpha_cap / (1.0 + (B_r_dummy**2 + B_phi_dummy**2) / B_eq**2)
!       alpha_Bphi = alpha_total * B_phi_dummy
!       alpha_Br = alpha_total * B_r_dummy
!       Uz_Br = B_r_dummy * U_z_cap
!       Uz_Bphi = B_phi_dummy * U_z_cap


!       call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
!       call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
!       call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
!       call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
!       call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
!       call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)
!       call spatial_derivative(alpha_m_dummy, 6, d_alpha_m, d2_alpha_m)
!       call spatial_derivative(U_z_cap, 6, der_u, d_sq_u)


!       ! vishniac_term = -f_para * R_omega*(2.*l/3)**2 *x
!       vishniac_term = f_para * x* exp(-x**2)

      
!       dBrdt = (- R_alpha*d_alpha_Bphi) + d2Br !- R_U*d_Uz_Br ! CHANGE

!       dBphidt = (R_omega*B_r_dummy) + d2Bphi !+ d_alpha_Br &
!       !- R_U*d_Uz_Bphi

!       ! dalpdt = -2./3.*t_d(alpha_total*(B_r_dummy**2. + B_phi_dummy**2.)/B_eq**2. &
!       ! (B_phi_dummy*dBphi_dummy-B_r_dummy*dBr_dummy)/B_eq**2. + vishniac_term) &
!       ! +R_k*d2_alpha_m - R_U*U_z_cap*d_alpha_m-R_U*alpha_m_dummy*d_U_z_cap

!       ! dalpdt = -2.0 / 3.0 * t_d * (alpha_total * (B_r_dummy**2.0 + B_phi_dummy**2.0) / B_eq**2.0 +&
!       !      (B_phi_dummy * dBphi - B_r_dummy * dBr) / B_eq**2.0 +vishniac_term)&
!       !      + R_k * d2_alpha_m - R_U * U_z_cap * d_alpha_m - R_U * alpha_m_dummy * d_U_z_cap

!       dalpdt =    (-2d0*(h**2/l**2))*((alpha_total)*((B_r_dummy**2)+(B_phi_dummy**2))&
!       - (((B_phi_dummy*dbr)-(B_r_dummy*dBphi))/R_alpha) + (alpha_m_dummy*R_m_inv))+ &
!       (R_omega/R_alpha) * vishniac_term

!     end subroutine sharanyas_notes_eqxn


! !NOTE:phi dummy is dphi/dz and dphi_term_dz is second derivative of phi
    
!     subroutine potential_equations(Phi_dummy, T_dummy, alpha_m_dummy)
!       double precision,intent(inout), dimension(nx) :: Phi_dummy,T_dummy, alpha_m_dummy
!       double precision, dimension(nx) :: dphi_dz,dT_dz,d_alpha_m_dz
!       double precision, dimension(nx) :: d2phi_dz,d2T_dz,d2_alpha_m_dz,alpha_total
!       double precision, dimension(nx) :: U_z_T, U_z_Phi,U_z_alm,vishniac_term
!       double precision, dimension(nx) :: d_U_z_T, d_U_z_Phi,d_U_z_alm 
!       double precision, dimension(nx) :: d2_U_z_T, d2_U_z_Phi,d2_U_z_alm 
!       double precision, dimension(nx) :: dalp_dz, d2alp_dz  
      
!       character(len=30) :: ghost_zone_type = 'anti-symmetric'
!       character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
!       character(len=30) :: ghost_zone_type3 = 'symmetric'

!       call impose_boundary_conditions(T_dummy, ghost_zone_type)
!       call impose_boundary_conditions(Phi_dummy, ghost_zone_type3)
!       call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type2)

!       alpha_total = alpha_m_dummy !+ alpha_k
!       U_z_T = T_dummy * U_z_cap
!       U_z_alm = alpha_m_dummy * U_z_cap
!       vishniac_term = f_para * x* exp(-x**2)




!       call spatial_derivative(Phi_dummy, 6, dphi_dz, d2phi_dz)
!       call spatial_derivative(T_dummy, 6, dT_dz, d2T_dz)
!       call spatial_derivative(alpha_m_dummy, 6, d_alpha_m_dz, d2_alpha_m_dz)
!       call spatial_derivative(U_z_T, 6, d_U_z_T, d2_U_z_T)
!       call spatial_derivative(U_z_Phi, 6, d_U_z_Phi, d2_U_z_Phi)
!       call spatial_derivative(U_z_alm, 6,dalp_dz, d2alp_dz)

!       dTdt = - R_U*d_U_z_T - R_omega*dphi_dz - alpha_total*d2phi_dz -dalp_dz*dphi_dz + d2T_dz 

!       dphi_dt = - R_U*U_z_cap*dphi_dz + alpha_total*T_dummy + d2phi_dz

!       d_alpha_m_dt = (-2d0*(h**2/l**2)) * (alpha_total*((T_dummy**2)+(dphi_dz**2)) &
!       - (-T_dummy*d2phi_dz + dphi_dz*dT_dz)+(R_omega/R_alpha) * vishniac_term) &
!       - R_U*d_U_z_alm + R_k*d2_alpha_m_dz 


!     end subroutine potential_equations

    subroutine nvf(Phi_dummy, T_dummy, Fr_dummy, F_phi_dummy,Fz_dummy, alpha_m_dummy)

      double precision,intent(inout), dimension(nx) :: Phi_dummy,T_dummy, alpha_m_dummy
      double precision,intent(inout),  dimension(nx) :: Fr_dummy, F_phi_dummy,Fz_dummy

      double precision, dimension(nx) :: dphi_dz,dT_dz,d_alpha_m_dz
      double precision, dimension(nx) :: d2phi_dz,d2T_dz,d2_alpha_m_dz,alpha_total
      double precision, dimension(nx) :: U_z_T, U_z_Phi,U_z_alm,vishniac_term
      double precision, dimension(nx) :: d_U_z_T, d_U_z_Phi,d_U_z_alm 
      double precision, dimension(nx) :: d2_U_z_T, d2_U_z_Phi,d2_U_z_alm 
      double precision, dimension(nx) :: dFr_dz, dF_phi_dz,dFz_dz
      double precision, dimension(nx) :: d2Fr_dz, d2F_phi_dz,d2Fz_dz
      
      character(len=30) :: ghost_zone_type = 'anti-symmetric'
      character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
      character(len=30) :: ghost_zone_type3 = 'symmetric'

      call impose_boundary_conditions(T_dummy, ghost_zone_type)
      call impose_boundary_conditions(Phi_dummy, ghost_zone_type3)
      call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type2)
      alpha_total = alpha_m_dummy + R_alpha_dim*alpha_cap 
      U_z_T = T_dummy * U_z ! NOTE: U_z HAS dimensions
      U_z_alm = alpha_m_dummy * U_z
      vishniac_term = -((2.*l/3.)**2)*f_para *R_omega*x *h_dim/t_d_dim

    !   vishniac_term = -((2.*l/3.)**2)*f_para *R_omega*x


      call spatial_derivative(Phi_dummy, 6, dphi_dz, d2phi_dz)
      call spatial_derivative(T_dummy, 6, dT_dz, d2T_dz)
      call spatial_derivative(alpha_m_dummy, 6, d_alpha_m_dz, d2_alpha_m_dz)
      call spatial_derivative(U_z_T, 6, d_U_z_T, d2_U_z_T)
      call spatial_derivative(U_z_Phi, 6, d_U_z_Phi, d2_U_z_Phi)
      call spatial_derivative(U_z_alm, 6,d_U_z_alm, d2_U_z_alm)
      call spatial_derivative(Fr_dummy, 6, dFr_dz, d2Fr_dz)
      call spatial_derivative(F_phi_dummy, 6, dF_phi_dz, d2F_phi_dz)
      call spatial_derivative(Fz_dummy, 6, dFz_dz, d2Fz_dz)

    !   dphi_dt = -U_z*dphi_dz + F_phi_dummy

    !   dTdt = - R_omega*dphi_dz - d_U_z_T  + dFr_dz + (Fz_dummy/radius)
    !   !NOTE: flat rotn curve, q = 1

    !   dFr_dt = (1./tau_c)*(-(alpha_total*dphi_dz) + dT_dz -Fr_dummy)
    !   dFphi_dt = (1./tau_c)*((alpha_total*T_dummy)+ d2phi_dz - F_phi_dummy )
    !   dFz_dt = (1./tau_c)*(-Fz_dummy)

      
    !   d_alpha_m_dt = (-2./(3.*tau_c)) * ( (-Fr_dummy*dphi_dz + F_phi_dummy*T_dummy)*(exp(x**2/h**2)/((radius**2)*(B0**2)))+&
    !   vishniac_term) - d_U_z_alm + R_k*d2_alpha_m_dz


      !CONCERN: DIMENSIONS!? especially alpha

!       dphi_dt = -U_z*dphi_dz + t_d_dim*F_phi_dummy

!       dTdt = - R_omega*dphi_dz - d_U_z_T  + (t_d_dim/h_dim)*dFr_dz + (t_d_dim/h_dim)*(Fz_dummy/radius)
!       !NOTE: flat rotn curve, q = 1

!       dFr_dt = (1./tau_c)*(-((alpha_m_dummy/h_dim + R_alpha*alpha_cap/t_d_dim)*dphi_dz) + (h_dim/t_d_dim)*dT_dz -Fr_dummy)

!       dFphi_dt = (1./tau_c)*(((alpha_m_dummy+ (R_alpha*alpha_cap*h_dim/t_d_dim))*T_dummy)+ (1./t_d_dim)*d2phi_dz - F_phi_dummy )
!       dFz_dt = (1./tau_c)*(-Fz_dummy)

!       d_alpha_m_dt = (-2./(3.*tau_c)) * ( (-Fr_dummy*dphi_dz/h_dim +&
!        F_phi_dummy*T_dummy)*(exp(x**2/h**2)/((h_dim**2)*(radius**2)*(B0**2)))+&
! vishniac_term) - d_U_z_alm + R_k*d2_alpha_m_dz


    dphi_dt = -U_z*dphi_dz/h_dim + F_phi_dummy

      dTdt = - R_omega*dphi_dz/h_dim - d_U_z_T/h_dim  + dFr_dz/h_dim + (Fz_dummy/radius_dim)
      !NOTE: flat rotn curve, q = 1

      dFr_dt = (1./tau_c_dim)*(-(alpha_total*dphi_dz/h_dim) + dT_dz/h_dim -Fr_dummy)
      dFphi_dt = (1./tau_c_dim)*((alpha_total*T_dummy)+ d2phi_dz/(h_dim**2.) - F_phi_dummy )
      dFz_dt = (1./tau_c_dim)*(-Fz_dummy)

      
      d_alpha_m_dt = (-2./(3.*tau_c_dim)) * ( (-Fr_dummy*dphi_dz/h_dim + &
      F_phi_dummy*T_dummy)*(exp(x**2)/((radius_dim**2)*(B0**2)))+&
      vishniac_term) - d_U_z_alm/h_dim + R_k_dim*d2_alpha_m_dz/(h_dim**2.)






    end subroutine nvf







end module equations
!NOTE: look at ss21 and see profiles of tau,Uz,rho,u,eta.
!REFER: 2020 chamandy,shukurov 
!tue10