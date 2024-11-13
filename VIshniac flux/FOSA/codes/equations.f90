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
!       alpha_total = alpha_m_dummy!+R_alpha*alpha_cap
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


!       vishniac_term = f_para * x* exp(-x**2)

      
!       dBrdt = (- R_alpha*d_alpha_Bphi) + d2Br !- R_U*d_Uz_Br ! CHANGE

!       dBphidt = (R_omega*B_r_dummy) + d2Bphi !+ d_alpha_Br &
!       !- R_U*d_Uz_Bphi

     
!       dalpdt =    (-2d0*(h**2/l**2))*((alpha_total)*((B_r_dummy**2)+(B_phi_dummy**2))&
!       - (((B_phi_dummy*dbr)-(B_r_dummy*dBphi))/R_alpha) + (alpha_m_dummy*R_m_inv))!+ &
!       !(!R_omega/R_alpha) * vishniac_term 

! !CONCERN: this last uncommented part of alpha avolution is needed for sharanyas equation plots. I uncommented becauase of an unrelated div by 0 error.

!     end subroutine sharanyas_notes_eqxn



     
    ! subroutine nvf_pdf_equations(B_r_dummy, B_phi_dummy, alpha_m_dummy)
    !   double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
    !   double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
    !   double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
    !   double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
    !   double precision,  dimension(nx) :: Uz_Bphi, Uz_Br
    !   double precision,  dimension(nx) :: alpha_Bphi, alpha_Br
    !   double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
    !   double precision, dimension(nx) :: alpha_total
    !   double precision, dimension(nx) :: vishniac_term
    !   double precision, dimension(nx) :: der_u, d_sq_u
    !   double precision, dimension(nx) :: d2_alpha_Br, d2_alpha_Bphi
    !   double precision, dimension(nx) :: d2_Uz_Br, d2_Uz_Bphi

    !   character(len=30) :: ghost_zone_type = 'anti-symmetric'
    !   character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'

    !   call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
    !   call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
    !   call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type2)
    !   alpha_total = alpha_m_dummy+R_alpha*alpha_cap
    !   alpha_Bphi = alpha_total * B_phi_dummy
    !   alpha_Br = alpha_total * B_r_dummy
    !   Uz_Br = B_r_dummy * U_z_cap
    !   Uz_Bphi = B_phi_dummy * U_z_cap


    !   call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
    !   call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
    !   call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
    !   call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
    !   call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
    !   call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)
    !   call spatial_derivative(alpha_m_dummy, 6, d_alpha_m, d2_alpha_m)
    !   call spatial_derivative(U_z_cap, 6, der_u, d_sq_u)

    !  vishniac_term = (((2.*l)/(3.))**2) * f_para * R_omega * x!(-2d0/3d0*tau_c)*((2d0*l)/(3d0*h))**2*
      

    !   dBrdt = - d_alpha_Bphi &
    !   + d2Br - R_U*d_Uz_Br 

    !   dBphidt = R_omega*B_r_dummy + d_alpha_Br &
    !   + d2Bphi - R_U*d_Uz_Bphi


    !   dalpdt =    (-2./((l**2)*(B_eq**2)))*((alpha_total)*(B_r_dummy**2+B_phi_dummy**2)&
    !   - (B_phi_dummy*dbr-B_r_dummy*dBphi)+vishniac_term)&
    !   - R_u*alpha_m_dummy*der_u-R_u*U_z*d_alpha_m + R_k*d2_alpha_m





    ! end subroutine nvf_pdf_equations

 

    ! subroutine MTA_nvf(B_r_dummy, B_phi_dummy, Fr_dummy, Fphi_dummy, Er_dummy, Ephi_dummy, alpha_m_dummy)
    !   double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
    !   double precision, intent(inout), dimension(nx) :: Fr_dummy, Fphi_dummy, Er_dummy, Ephi_dummy
    !   double precision, dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
    !   double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
    !   double precision, dimension(nx) :: Uz_Bphi, Uz_Br
    !   double precision, dimension(nx) :: d_Uz_Bphi, d_Uz_Br
    !   double precision, dimension(nx) :: d2_Uz_Br, d2_Uz_Bphi
    !   double precision, dimension(nx) :: alpha_Bphi, alpha_Br
    !   double precision, dimension(nx) :: d_alpha_Bphi, d_alpha_Br
    !   double precision, dimension(nx) :: d2_alpha_Br, d2_alpha_Bphi
    !   double precision, dimension(nx) :: alpha_total
    !   double precision, dimension(nx) :: vishniac_term
    !   double precision, dimension(nx) :: der_u, d_sq_u

    !   character(len=30) :: ghost_zone_type = 'anti-symmetric'
    !   character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
    !   character(len=30) :: ghost_zone_type3 = 'symmetric'
    !   character(len=30) :: ghost_zone_type4 = 'none'

    !   call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
    !   call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
    !   call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type2)
    !   call impose_boundary_conditions(Fr_dummy, ghost_zone_type4)
    !   call impose_boundary_conditions(Fphi_dummy, ghost_zone_type4)
    !   call impose_boundary_conditions(Er_dummy, ghost_zone_type4)
    !   call impose_boundary_conditions(Ephi_dummy, ghost_zone_type4)

    !   alpha_total = alpha_m_dummy + R_alpha*alpha_cap

    !   alpha_Bphi = alpha_total * B_phi_dummy
    !   alpha_Br = alpha_total * B_r_dummy
    !   Uz_Br = B_r_dummy *U_z
    !   Uz_Bphi = B_phi_dummy *U_z


    !   call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
    !   call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
    !   call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
    !   call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
    !   call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
    !   call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)
    !   call spatial_derivative(alpha_m_dummy, 6, d_alpha_m, d2_alpha_m)
    !   call spatial_derivative(U_z, 6, der_u, d_sq_u)

    !  vishniac_term = -(((2.*l)/(3.))**2) * f_para * R_omega * x!(-2d0/3d0*tau_c)*((2d0*l)/(3d0*h))**2*
      

    !   dBrdt = - d_Uz_Br + Fr_dummy + R_m_inv*d2Br

    !   dBphidt = - d_Uz_Bphi + Fphi_dummy + R_m_inv*d2Bphi + R_omega*B_r_dummy

    !   dFrdt = (1./tau) * ((-d_alpha_Bphi) + (d2Br) - Fr_dummy)

    !   dFphidt = (1./tau) * ((d_alpha_Br) + (d2Bphi) - Fphi_dummy)

    !   dErdt = (1./tau) * ((alpha_total*B_r_dummy) + (dBphi) - Er_dummy)

    !   dEphidt = (1./tau) * ((alpha_total*B_phi_dummy) - (dBr) - Ephi_dummy)

    !   dalpdt = (-2./(tau_c*3.)) * ((Er_dummy*B_r_dummy + Ephi_dummy*B_phi_dummy)/(B_eq**2) +&
    !   alpha_m_dummy*R_m_inv + vishniac_term) - alpha_m_dummy*der_u - U_z*d_alpha_m + R_k*d2_alpha_m
      




    ! end subroutine MTA_nvf

    subroutine FOSA(B_r_dummy, B_phi_dummy, alpha_m_dummy)
      double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
      double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
      double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
      double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
      double precision,  dimension(nx) :: Uz_Bphi, Uz_Br
      double precision,  dimension(nx) :: alpha_Bphi, alpha_Br
      double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
      double precision, dimension(nx) :: alpha_total
      double precision, dimension(nx) :: vishniac_term
      double precision, dimension(nx) :: der_u, d_sq_u
      double precision, dimension(nx) :: d2_alpha_Br, d2_alpha_Bphi
      double precision, dimension(nx) :: d2_Uz_Br, d2_Uz_Bphi
      double precision :: R_alpha_inv = 0.
      character(len=30) :: ghost_zone_type = 'anti-symmetric'
      character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'

      ! call construct_eta_profile
      ! call construct_alpha_profile
      ! call construct_velocity_profile
      call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
      call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
      call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type2)

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

!NOTE: vishniac term for 2017 new notes with old u stratification
        vishniac_term = (8./9.)*l**2 * f_para * R_omega * x*exp(x**2)!(((2.*l)/(3.))**2) * f_para * R_omega * x!(-2d0/3d0*tau_c)*((2d0*l)/(3d0*h))**2*
!NOTE: vishniac term for 2017 new notes with new u stratification
      !a)tau = not constant
        ! vishniac_term = (4./27.)*l**2 * f_para * R_omega * small_u_0**2
      !b)tau = constant
        ! vishniac_term = (8./81.)*tau**2 * f_para*x * R_omega * small_u_0**4!l=tau*u
!NOTE: vishniac term derived from GS+23
      ! vishniac_term = ((7.*f_para*tau**2 * R_omega * x * small_u_0**4)/5400.)-((f_para*tau* R_omega *  small_u_0**2)/(324.*pi*rho))

        dBrdt = -d_alpha_Bphi +eta*d2Br - d_Uz_Br

        dBphidt = R_omega*B_r_dummy+eta*d2Bphi-d_Uz_Bphi +d_alpha_Br


        dalpdt =    (-2./(3.*tau))*((alpha_total)*(B_r_dummy**2+B_phi_dummy**2)/B_eq**2&
        - eta*(B_phi_dummy*dbr-B_r_dummy*dBphi)/B_eq**2+ R_m_inv*alpha_m + vishniac_term)&
        - alpha_m_dummy*der_u-U_z*d_alpha_m + eta*R_k*d2_alpha_m





    end subroutine FOSA


end module equations

! run 119---D
! rest in numerical order