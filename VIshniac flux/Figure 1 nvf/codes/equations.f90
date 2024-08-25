module equations
  
  
  use parameters
  use velocity_profile
  use alpha_profile
  use initial_field
  use time_grid
  use physical_grid
  use make_a_grid
  
  implicit none
  double precision, dimension(nx) :: dBrdt, dBphidt, dalpdt 

  contains
  !contains all the equations used for simulation.
  !Currently included:
    !1! diff_equations_no_split - simple alp-omg dynamo, no quenching
    !2! dyn_quenching_nvf - dynamic quenching, no vishniac flux
    !3! sharanyas_notes_eqxn - the equations with nvf, 0 alp_k, as per notes of sharanya
  !TODO: Add the equations fron nvf paper draft.

    subroutine diff_equations_no_split(B_r_dummy, B_phi_dummy)
      double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy
      double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
      double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
      double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
      double precision,  dimension(nx) :: Uz_Bphi, Uz_Br
      double precision,  dimension(nx) :: alpha_Bphi, alpha_Br

      character(len=30) :: ghost_zone_type = 'anti-symmetric'
      call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
      call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)

      alpha_cap2 = alpha_cap / (1.0 + (B_r_dummy**2 + B_phi_dummy**2) / B_eq**2)
      alpha_Bphi = alpha_cap2 * B_phi_dummy
      alpha_Br = alpha_cap2 * B_r_dummy
      Uz_Br = B_r_dummy * U_z_cap
      Uz_Bphi = B_phi_dummy * U_z_cap

      call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
      call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
      call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
      call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
      call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
      call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)


      dBrdt = - R_alpha*d_alpha_Bphi &
      + d2Br - R_U*d_Uz_Br
      dBphidt = R_omega*B_r_dummy + R_alpha*d_alpha_Br &
      + d2Bphi - R_U*d_Uz_Bphi

    end subroutine diff_equations_no_split


    subroutine dyn_quenching_nvf(B_r_dummy, B_phi_dummy, alpha_m_dummy)
      double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
      double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
      double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
      double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
      double precision,  dimension(nx) :: Uz_Bphi, Uz_Br,der_u, d_sq_u
      double precision,  dimension(nx) :: alpha_Bphi, alpha_Br
      double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
      double precision, dimension(nx) :: alpha_total
      double precision, dimension(nx) :: vishniac_term

      character(len=30) :: ghost_zone_type = 'anti-symmetric'
      call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
      call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
      call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type)
      alpha_total = R_alpha*alpha_k + alpha_m_dummy
      ! alpha_cap2 = alpha_cap / (1.0 + (B_r_dummy**2 + B_phi_dummy**2) / B_eq**2)
      alpha_Bphi = alpha_total * B_phi_dummy
      alpha_Br = alpha_total * B_r_dummy
      Uz_Br = B_r_dummy * U_z_cap
      Uz_Bphi = B_phi_dummy * U_z_cap


      call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
      call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
      call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
      call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
      call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
      call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)
      call spatial_derivative(alpha_m_dummy, 6, d_alpha_m, d2_alpha_m)
      call spatial_derivative(U_z_cap, 6, der_u, d_sq_u)


      ! vishniac_term = x* f_para * R_omega*(2*l/3*h)**2 
      ! vishniac_term = -f_para * R_omega *x

      
      dBrdt = - d_alpha_Bphi &
      + d2Br - R_U*d_Uz_Br

      dBphidt = R_omega*B_r_dummy + d_alpha_Br &
      + d2Bphi - R_U*d_Uz_Bphi

      ! dalpdt = -2./3.*t_d(alpha_total*(B_r_dummy**2. + B_phi_dummy**2.)/B_eq**2. &
      ! (B_phi_dummy*dBphi_dummy-B_r_dummy*dBr_dummy)/B_eq**2. + vishniac_term) &
      ! +R_k*d2_alpha_m - R_U*U_z_cap*d_alpha_m-R_U*alpha_m_dummy*d_U_z_cap

      ! dalpdt = -2.0 / 3.0 * t_d * (alpha_total * (B_r_dummy**2.0 + B_phi_dummy**2.0) / B_eq**2.0 +&
      !      (B_phi_dummy * dBphi - B_r_dummy * dBr) / B_eq**2.0 +vishniac_term)&
      !      + R_k * d2_alpha_m - R_U * U_z_cap * d_alpha_m - R_U * alpha_m_dummy * d_U_z_cap

      dalpdt =    (-2/((3*t_d)*(B_eq**2)))*((alpha_total)*(B_r_dummy**2+B_phi_dummy**2)&
      - (B_phi_dummy*dbr-B_r_dummy*dBphi))&
      - R_u*alpha_m_dummy*der_u-R_u*U_z*d_alpha_m + R_k*d2_alpha_m

    end subroutine dyn_quenching_nvf
  
  
    subroutine sharanyas_notes_eqxn(B_r_dummy, B_phi_dummy, alpha_m_dummy)
      double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
      double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
      double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
      double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
      double precision,  dimension(nx) :: Uz_Bphi, Uz_Br,der_u, d_sq_u
      double precision,  dimension(nx) :: alpha_Bphi, alpha_Br
      double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
      double precision, dimension(nx) :: alpha_total
      double precision, dimension(nx) :: vishniac_term

      character(len=30) :: ghost_zone_type = 'anti-symmetric'
      character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'

      call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
      call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
      call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type2)
      alpha_total = alpha_m_dummy!R_alpha*alpha_k + 
      ! alpha_cap2 = alpha_cap / (1.0 + (B_r_dummy**2 + B_phi_dummy**2) / B_eq**2)
      alpha_Bphi = alpha_total * B_phi_dummy
      alpha_Br = alpha_total * B_r_dummy
      Uz_Br = B_r_dummy * U_z_cap
      Uz_Bphi = B_phi_dummy * U_z_cap


      call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
      call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
      call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
      call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
      call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
      call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)
      call spatial_derivative(alpha_m_dummy, 6, d_alpha_m, d2_alpha_m)
      call spatial_derivative(U_z_cap, 6, der_u, d_sq_u)


      ! vishniac_term = -f_para * R_omega*(2.*l/3)**2 *x
      vishniac_term = f_para * x* exp(-x**2)

      
      dBrdt = (- R_alpha*d_alpha_Bphi) + d2Br !- R_U*d_Uz_Br ! CHANGE

      dBphidt = (R_omega*B_r_dummy) + d2Bphi !+ d_alpha_Br &
      !- R_U*d_Uz_Bphi

      ! dalpdt = -2./3.*t_d(alpha_total*(B_r_dummy**2. + B_phi_dummy**2.)/B_eq**2. &
      ! (B_phi_dummy*dBphi_dummy-B_r_dummy*dBr_dummy)/B_eq**2. + vishniac_term) &
      ! +R_k*d2_alpha_m - R_U*U_z_cap*d_alpha_m-R_U*alpha_m_dummy*d_U_z_cap

      ! dalpdt = -2.0 / 3.0 * t_d * (alpha_total * (B_r_dummy**2.0 + B_phi_dummy**2.0) / B_eq**2.0 +&
      !      (B_phi_dummy * dBphi - B_r_dummy * dBr) / B_eq**2.0 +vishniac_term)&
      !      + R_k * d2_alpha_m - R_U * U_z_cap * d_alpha_m - R_U * alpha_m_dummy * d_U_z_cap

      dalpdt =    (-2d0*(h**2/l**2))*((alpha_total)*((B_r_dummy**2)+(B_phi_dummy**2))&
      - (((B_phi_dummy*dbr)-(B_r_dummy*dBphi))/R_alpha) + (alpha_m_dummy*R_m_inv))+ &
      (R_omega/R_alpha) * vishniac_term

    end subroutine sharanyas_notes_eqxn




end module equations
