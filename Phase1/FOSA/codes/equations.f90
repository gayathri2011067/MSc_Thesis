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
    double precision, dimension(nx) :: dBrdt, dBphidt, dalpdt, dFrdt, dFphidt, dErdt, dEphidt 
    contains
    


!     subroutine FOSA(B_r_dummy, B_phi_dummy, alpha_m_dummy)
!       double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
!       double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
!       double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
!       double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
!       double precision,  dimension(nx) :: Uz_Bphi, Uz_Br
!       double precision,  dimension(nx) :: alpha_Bphi, alpha_Br
!       double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
!       double precision, dimension(nx) :: alpha_total
!       double precision, dimension(nx) :: vishniac_term,new_diffusive_flux
!       double precision, dimension(nx) :: der_u, d_sq_u
!       double precision, dimension(nx) :: d2_alpha_Br, d2_alpha_Bphi
!       double precision, dimension(nx) :: d2_Uz_Br, d2_Uz_Bphi
!       double precision :: R_alpha_inv = 0.
!       character(len=30) :: ghost_zone_type = 'anti-symmetric'
!       character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'

!       ! call construct_eta_profile
!       ! call construct_alpha_profile
!       ! call construct_velocity_profile
!       call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
!       call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
!       call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type2)

!       alpha_total = alpha_m_dummy + R_alpha*alpha_cap

!       alpha_Bphi = alpha_total * B_phi_dummy
!       alpha_Br = alpha_total * B_r_dummy
!       Uz_Br = B_r_dummy *U_z
!       Uz_Bphi = B_phi_dummy *U_z


!       call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
!       call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
!       call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
!       call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
!       call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
!       call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)
!       call spatial_derivative(alpha_m_dummy, 6, d_alpha_m, d2_alpha_m)
!       call spatial_derivative(U_z, 6, der_u, d_sq_u)


! !-----------------------------DIFFERENT EXPRESSIONS FOR VISHNIAC TERM---------------------------------
! !NOTE: vishniac term with old u stratification and no rho stratification
!         vishniac_term = (2./9.)**2 * f_para * R_omega * tau**2 * small_u_0**4 * 4*x * exp(2*x**2)

! !NOTE: vishniac term with new u stratification and no rho stratification
!       ! vishniac_term = (2./9.)**2 * f_para * R_omega * tau**2 * small_u_0**4 * 2*x 


! !NOTE: vishniac term with old u stratification and rho stratification
!         ! vishniac_term = (2./9.)**2 * f_para * R_omega * tau**2 * small_u_0**4 * 2*x * exp(2*x**2)


! !NOTE: vishniac term with new u stratification and rho stratification
!         ! vishniac_term = (2./9.)**2 * f_para * R_omega * tau**2 * small_u_0**4 * (1-x**2) * 2*x         


! !NOTE: vishniac term with rho stratification and no u stratification
!         ! vishniac_term = -(2./3.)**2 *l**2  * small_u_0**2 * f_para * R_omega * x
              
! !----------------------------------------------------------------------------------------------------    
      
!       !CONCERN: there is adisparity between the exponential term 


! !-----------------------------OLD EQUATION OF DIFFUSIVE FLUX---------------------------------      

!         dBrdt = -d_alpha_Bphi +eta*d2Br - d_Uz_Br

!         dBphidt = R_omega*B_r_dummy+eta*d2Bphi -d_Uz_Bphi +d_alpha_Br


!         dalpdt =    (-2./(3.*tau))*((alpha_total)*(B_r_dummy**2+B_phi_dummy**2)/B_eq**2&
!         - eta*(B_phi_dummy*dbr-B_r_dummy*dBphi)/B_eq**2+ R_m_inv*alpha_m + vishniac_term)&
!          + eta*R_k*d2_alpha_m -U_z*d_alpha_m- alpha_m_dummy*der_u

! !----------------------------------------------------------------------------------------------------


! !-----------------------------NEW EQUATION OF DIFFUSIVE FLUX--------------------------------- 
        
        
!         ! new_diffusive_flux = 7 *eta**2*rho*d_alpha_m*(xi-1) !------------------->>NEW DIFFUSIVE FLUX !CONCERN: ONLY VALID WHEN RHO IS ONSTANT

!         ! dBrdt = -d_alpha_Bphi +eta*d2Br - d_Uz_Br

!         ! dBphidt = R_omega*B_r_dummy+eta*d2Bphi-d_Uz_Bphi +d_alpha_Br


!         ! dalpdt =    (-2./(3.*tau))*((alpha_total)*(B_r_dummy**2+B_phi_dummy**2)/B_eq**2&
!         ! - eta*(B_phi_dummy*dbr-B_r_dummy*dBphi)/B_eq**2+ R_m_inv*alpha_m + vishniac_term)&
!         ! - alpha_m_dummy*der_u-U_z*d_alpha_m + new_diffusive_flux


! !----------------------------------------------------------------------------------------------------



!     end subroutine FOSA




    subroutine MTA_nvf(B_r_dummy, B_phi_dummy, Fr_dummy, Fphi_dummy, Er_dummy, Ephi_dummy, alpha_m_dummy)
        double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
        double precision, intent(inout), dimension(nx) :: Fr_dummy, Fphi_dummy, Er_dummy, Ephi_dummy
        double precision, dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
        double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
        double precision, dimension(nx) :: Uz_Bphi, Uz_Br
        double precision, dimension(nx) :: d_Uz_Bphi, d_Uz_Br
        double precision, dimension(nx) :: d2_Uz_Br, d2_Uz_Bphi
        double precision, dimension(nx) :: alpha_Bphi, alpha_Br
        double precision, dimension(nx) :: d_alpha_Bphi, d_alpha_Br
        double precision, dimension(nx) :: d2_alpha_Br, d2_alpha_Bphi
        double precision, dimension(nx) :: alpha_total
        double precision, dimension(nx) :: vishniac_term
        double precision, dimension(nx) :: der_u, d_sq_u
  
        character(len=30) :: ghost_zone_type = 'anti-symmetric'
        character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
        character(len=30) :: ghost_zone_type3 = 'symmetric'
        character(len=30) :: ghost_zone_type4 = 'none'
  
        call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
        call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
        call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type2)
        call impose_boundary_conditions(Fr_dummy, ghost_zone_type4)
        call impose_boundary_conditions(Fphi_dummy, ghost_zone_type4)
        call impose_boundary_conditions(Er_dummy, ghost_zone_type4)
        call impose_boundary_conditions(Ephi_dummy, ghost_zone_type4)
  
        alpha_total = alpha_m_dummy + R_alpha*alpha_cap
  
        alpha_Bphi = alpha_total * B_phi_dummy
        alpha_Br = alpha_total * B_r_dummy
        Uz_Br = B_r_dummy *U_z
        Uz_Bphi = B_phi_dummy *U_z
  
  
        call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
        call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
        call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
        call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
        call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
        call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)
        call spatial_derivative(alpha_m_dummy, 6, d_alpha_m, d2_alpha_m)
        call spatial_derivative(U_z, 6, der_u, d_sq_u)
  
        vishniac_term = (2./9.)**2 * f_para * R_omega * tau**2 * small_u_0**4 * 4*x * exp(2*x**2)
        
  
        dBrdt = - d_Uz_Br + Fr_dummy + R_m_inv*d2Br
  
        dBphidt = - d_Uz_Bphi + Fphi_dummy + R_m_inv*d2Bphi + R_omega*B_r_dummy !+ d_alpha_Br
  
        dFrdt = (1./tau) * ((-d_alpha_Bphi) + (d2Br) - Fr_dummy)
  
        dFphidt = (1./tau) * ((d_alpha_Br) + (d2Bphi) - Fphi_dummy)
  
        dErdt = (1./tau) * ((alpha_total*B_r_dummy) + (dBphi) - Er_dummy)
  
        dEphidt = (1./tau) * ((alpha_total*B_phi_dummy) - (dBr) - Ephi_dummy)
  
        dalpdt = (-2./(tau_c*3.)) * ((Er_dummy*B_r_dummy + Ephi_dummy*B_phi_dummy)/(B_eq**2) +&
        alpha_m_dummy*R_m_inv + vishniac_term) - alpha_m_dummy*der_u - U_z*d_alpha_m + R_k*d2_alpha_m
        
  
  
  
  
      end subroutine MTA_nvf
  
  

end module equations

! run 119---D
! rest in numerical order