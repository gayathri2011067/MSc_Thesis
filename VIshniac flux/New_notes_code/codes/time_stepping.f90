module timestepping
    use parameters
    use eta_profile
    use alpha_profile
    use omega_profile
    use initial_field
    use time_grid
    use physical_grid
    use make_a_grid
    use equations

    implicit none
    

    double precision, dimension(nx) :: k1r,k1phi,k2r,k2phi,k3r,k3phi,k4r,k4phi
    double precision, dimension(nx) :: k1Fr,k2FR,k3Fr,k4Fr
    double precision, dimension(nx) :: k1Fphi,k2Fphi,k3Fphi,k4Fphi
    double precision, dimension(nx) :: k1Er,k2Er,k3Er,k4Er
    double precision, dimension(nx) :: k1Ephi,k2Ephi,k3Ephi,k4Ephi
    double precision, dimension(nx) :: k1alpha,k2alpha,k3alpha,k4alpha
    double precision,dimension(nx) :: g_1, g_2, g_3, z_1, z_2
    contains
        

        ! subroutine forward_difference
        !     double precision, dimension(nx) :: temp_Br, temp_Bphi


        !     do i = nxghost+1,nxphys+nxghost
        !         temp_Br(i) = B_r(i) - (beta/2)*(alpha_Bphi(i+1)-alpha_Bphi(i-1)) &
        !         +alp*(B_r(i+1)-2*B_r(i)+B_r(i-1)) - (beta/2)*(Uz_Br(i+1) - Uz_Br(i-1))
        !         temp_Bphi(i) = B_phi(i) - (beta/2)*(alpha_Br(i+1)-alpha_Br(i-1)) &
        !         +alp*(B_phi(i+1)-2*B_phi(i)+B_phi(i-1)) - (beta/2)*(Uz_Bphi(i+1) - Uz_Bphi(i-1))&
        !         +G*B_r(i)

                
        !     end do
        !     B_r = temp_Br
        !     B_phi = temp_Bphi
        !     t = t + dt
        ! end subroutine forward_difference

        ! subroutine backward_difference
        !     ! character(len=30) :: ghost_zone_type = 'anti-symmetric'

        !     double precision, dimension(nx) :: temp_Br, temp_Bphi
        !     ! call impose_boundary_conditions(B_r, ghost_zone_type)
        !     ! call impose_boundary_conditions(B_phi, ghost_zone_type)
        !     ! call construct_grid
        !     ! alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
        !     ! alpha_Bphi = alpha_cap2*B_phi
        !     ! alpha_Br = alpha_cap2*B_r
        !     ! Uz_Br = B_r*U_z_cap
        !     ! Uz_Bphi = B_phi*U_z_cap

        !     ! call spatial_derivative(B_r,2,dBr,d2Br)
        !     ! call spatial_derivative(B_phi,2,dBphi,d2Bphi)
        !     ! call spatial_derivative(alpha_Br,2,d_alpha_Br,d2_alpha_Br)
        !     ! call spatial_derivative(alpha_Bphi,2,d_alpha_Bphi,d2_alpha_Bphi)
        !     ! call spatial_derivative(Uz_Br,2,d_Uz_Br,d2_Uz_Br)
        !     ! call spatial_derivative(Uz_Bphi,2,d_Uz_Bphi,d2_Uz_Bphi)

        !     do i = nxghost+1,nxphys+nxghost
        !         temp_Br(i) = B_r(i) + (beta/2)*(alpha_Bphi(i)-alpha_Bphi(i-2)) &
        !         -alp*(B_r(i+1)-2*B_r(i)+B_r(i-1)) + (beta/2)*(Uz_Br(i) - Uz_Br(i-2))
        !         temp_Bphi(i) = B_phi(i) - (beta/2)*(alpha_Br(i)-alpha_Br(i-2)) &
        !         -alp*(B_phi(i+1)-2*B_phi(i)+B_phi(i-1)) + (beta/2)*(Uz_Bphi(i) - Uz_Bphi(i-2))&
        !         -G*B_r(i)

        !     end do
        !     B_r = temp_Br
        !     B_phi = temp_Bphi
        !     t = t + dt

        ! end subroutine backward_difference

        ! subroutine central_difference
        !     double precision, dimension(nx) :: temp_Br, temp_Bphi
        !     ! character(len=30) :: ghost_zone_type = 'anti-symmetric'

        !     ! call impose_boundary_conditions(B_r, ghost_zone_type)
        !     ! call impose_boundary_conditions(B_phi, ghost_zone_type)
        !     ! call impose_boundary_conditions(old_Br, ghost_zone_type)
        !     ! call impose_boundary_conditions(old_Bphi, ghost_zone_type)
        !     ! call construct_grid
        !     ! alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
        !     ! alpha_Bphi = alpha_cap2*B_phi
        !     ! alpha_Br = alpha_cap2*B_r
        !     ! Uz_Br = B_r*U_z_cap
        !     ! Uz_Bphi = B_phi*U_z_cap


        !     do i = nxghost+1,nxphys+nxghost
        !         temp_Br(i) = (1+2*alp)*((beta)*(alpha_Bphi(i+1)-alpha_Bphi(i-1)) &
        !         +2*alp*(B_r(i+1)+B_r(i-1)) - (beta)*(Uz_Br(i+1) - Uz_Br(i-1)))+old_Br(i)

        !         temp_Bphi(i) = (1+2*alp+G)*((beta)*(alpha_Br(i+1)-alpha_Br(i-1)) &
        !         +2*alp*(B_phi(i+1)+B_phi(i-1)) - (beta)*(Uz_Bphi(i+1) - Uz_Bphi(i-1)))+old_Bphi(i) !NOT_SURE
                

        !     end do
        !     B_r = temp_Br
        !     B_phi = temp_Bphi
        !     t = t + dt

        ! end subroutine central_difference

        ! subroutine RK4_new
        !     character(len=30) :: ghost_zone_type = 'anti-symmetric'
        !     ! call impose_boundary_conditions(B_r, ghost_zone_type)
        !     ! call impose_boundary_conditions(B_phi, ghost_zone_type)

        !     ! alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
        !     ! alpha_Bphi = alpha_cap2*B_phi
        !     ! alpha_Br = alpha_cap2*B_r
        !     ! Uz_Br = B_r*U_z_cap
        !     ! Uz_Bphi = B_phi*U_z_cap

        !     ! call spatial_derivative(B_r,2,dBr,d2Br)
        !     ! call spatial_derivative(B_phi,2,dBphi,d2Bphi)
        !     ! call spatial_derivative(alpha_Br,2,d_alpha_Br,d2_alpha_Br)
        !     ! call spatial_derivative(alpha_Bphi,2,d_alpha_Bphi,d2_alpha_Bphi)
        !     ! call spatial_derivative(Uz_Br,2,d_Uz_Br,d2_Uz_Br)
        !     ! call spatial_derivative(Uz_Bphi,2,d_Uz_Bphi,d2_Uz_Bphi)




        !     call MTA_nvf(B_r_dummy, B_phi_dummy, Fr_dummy, Fphi_dummy,alpha_m_dummy, Er_dummy, Ephi_dummy)
        !     k1r = dt*dBrdt
        !     k1phi = dt*dBphidt
        !     k1alpha = dt*dalpdt
        !     B_r = B_r+0.5*k1r
        !     B_phi = B_phi+0.5*k1phi
        !     alpha_m = alpha_m+0.5*k1alpha

        !     call MTA_nvf(B_r_dummy, B_phi_dummy, Fr_dummy, Fphi_dummy,alpha_m_dummy, Er_dummy, Ephi_dummy)
        !     k2r = dt*dBrdt
        !     k2phi = dt*dBphidt
        !     k2alpha = dt*dalpdt
        !     B_r = B_r+0.5*k2r
        !     B_phi = B_phi+0.5*k2phi
        !     alpha_m = alpha_m+0.5*k2alpha

        !     call MTA_nvf(B_r_dummy, B_phi_dummy, Fr_dummy, Fphi_dummy,alpha_m_dummy, Er_dummy, Ephi_dummy)
        !     k3r = dt*dBrdt
        !     k3phi = dt*dBphidt
        !     k3alpha = dt*dalpdt
        !     B_r = B_r+k3r
        !     B_phi = B_phi+k3phi
        !     alpha_m = alpha_m+k3alpha

        !     call MTA_nvf(B_r_dummy, B_phi_dummy, Fr_dummy, Fphi_dummy,alpha_m_dummy, Er_dummy, Ephi_dummy)
        !     k4r = dt*dBrdt
        !     k4phi = dt*dBphidt
        !     k4alpha = dt*dalpdt

        !     B_r = B_r + (k1r + 2.*k2r + 2.*k3r + k4r)/6.
        !     B_phi = B_phi + (k1phi + 2.*k2phi + 2.*k3phi + k4phi)/6.

        !     t = t + dt



        ! end subroutine RK4_new

        subroutine RK3_implicit

            implicit none
            
            ! character(len=30) :: ghost_zone_type = 'anti-symmetric'
            double precision, dimension(nx) :: Br_g, Bphi_g, alpha_m_g
            double precision, dimension(nx) :: Fr_g, Fphi_g, Er_g, Ephi_g

            g_1 = 8.0 / 15.0        !
            g_2 = 5.0 / 12.0        !
            g_3 = 3.0 / 4.0         !----> REFER: Brandenburg, A. (2001).ArXiv. 
            z_1 = -17.0 / 60.0      !             Computational aspects of astrophysical MHD and turbulence.  
            z_2 = -5.0 / 12.0       !             https://doi.org/10.48550/arXiv.astro-ph/0109497
            

            ! call impose_boundary_conditions(B_r, ghost_zone_type)
            ! call impose_boundary_conditions(B_phi, ghost_zone_type)

            ! ! alpha_cap2 = alpha_cap / (1.0 + (B_r**2 + B_phi**2) / B_eq**2)
            ! ! alpha_Bphi = alpha_cap2 * B_phi
            ! ! alpha_Br = alpha_cap2 * B_r
            ! ! Uz_Br = B_r * U_z_cap
            ! ! Uz_Bphi = B_phi * U_z_cap

            ! call spatial_derivative(B_r, 2, dBr, d2Br)
            ! call spatial_derivative(B_phi, 2, dBphi, d2Bphi)
            ! call spatial_derivative(alpha_Br, 2, d_alpha_Br, d2_alpha_Br)
            ! call spatial_derivative(alpha_Bphi, 2, d_alpha_Bphi, d2_alpha_Bphi)
            ! call spatial_derivative(Uz_Br, 2, d_Uz_Br, d2_Uz_Br)
            ! call spatial_derivative(Uz_Bphi, 2, d_Uz_Bphi, d2_Uz_Bphi)




            
            ! call diff_equations_no_split(B_r, B_phi)
            ! !STEP 1 !NOTE: have to call spatial derivatives with this new f and g
            ! k1r = dt * dBrdt
            ! k1phi = dt * dBphidt
            ! Br_f = Br_f + g_1 * k1r
            ! Bphi_f = Bphi_f + g_1 * k1phi
            ! Br_g = Br_f + z_1 * k1r
            ! Bphi_g = Bphi_f + z_1 * k1phi
            ! B_r=Br_f       !----> !NOTE: added this step because i wanted f to be carried to next step
            ! B_phi=Bphi_f  ! 
            
            ! call diff_equations_no_split(B_r, B_phi)
            ! k2r = dt * dBrdt
            ! k2phi = dt * dBphidt
            ! Br_f = Br_g + g_2 * k2r
            ! Bphi_f = Bphi_g + g_2 * k2phi
            ! Br_g = Br_f + z_2 * k2r
            ! Bphi_g = Bphi_f + z_2 * k2phi
            ! B_r=Br_f       !----> !NOTE: added this step because i wanted f to be carried to next step
            ! B_phi=Bphi_f  ! 

            ! ! 3rd step
            ! call diff_equations_no_split(B_r, B_phi)
            ! k3r = dt * dBrdt
            ! k3phi = dt * dBphidt
            ! Br_f = Br_g + g_3 * k1r
            ! Bphi_f = Bphi_g + g_3 * k3phi

            ! B_r=Br_f       !----> !NOTE: added this step because i wanted f to be carried to next step
            ! B_phi=Bphi_f  ! 
            ! t = t +  dt

            ! call MTA_nvf(B_r, B_phi, Fr, Fphi, Er, Ephi, alpha_m)
            ! !STEP 1 !NOTE: have to call spatial derivatives with this new f and g
            ! k1r = dt * dBrdt
            ! k1phi = dt * dBphidt
            ! k1Fr = dt * dFrdt
            ! k1Fphi = dt * dFphidt
            ! k1Er = dt * dErdt
            ! k1Ephi = dt * dEphidt
            ! k1alpha = dt * dalpdt

            ! B_r = B_r + g_1 * k1r
            ! B_phi = B_phi + g_1 * k1phi
            ! Fr = Fr + g_1 * k1Fr
            ! Fphi = Fphi + g_1 * k1Fphi
            ! Er = Er + g_1 * k1Er
            ! Ephi = Ephi + g_1 * k1Ephi
            ! alpha_m = alpha_m + g_1 * k1alpha

            ! Br_g = B_r + z_1 * k1r
            ! Bphi_g = B_phi + z_1 * k1phi
            ! Fr_g = Fr + z_1 * k1Fr
            ! Fphi_g = Fphi + z_1 * k1Fphi
            ! Er_g = Er + z_1 * k1Er
            ! Ephi_g = Ephi + z_1 * k1Ephi
            ! alpha_m_g = alpha_m + z_1 * k1alpha

            ! ! B_r=Br_f       !----> !NOTE: added this step because I wanted f to be carried to next step
            ! ! B_phi=Bphi_f  ! 
            
            ! call MTA_nvf(B_r, B_phi, Fr, Fphi,Er, Ephi, alpha_m)

            ! k2r = dt * dBrdt
            ! k2phi = dt * dBphidt
            ! k2Fr = dt * dFrdt
            ! k2Fphi = dt * dFphidt
            ! k2Er = dt * dErdt
            ! k2Ephi = dt * dEphidt
            ! k2alpha = dt * dalpdt

            ! B_r = Br_g + g_2 * k2r
            ! B_phi = Bphi_g + g_2 * k2phi
            ! Fr = Fr_g + g_2 * k2Fr
            ! Fphi = Fphi_g + g_2 * k2Fphi
            ! Er = Er_g + g_2 * k2Er
            ! Ephi = Ephi_g + g_2 * k2Ephi
            ! alpha_m = alpha_m_g + g_2 * k2alpha

            ! Br_g = B_r + z_2 * k2r
            ! Bphi_g = B_phi + z_2 * k2phi
            ! Fr_g = Fr + z_2 * k2Fr
            ! Fphi_g = Fphi + z_2 * k2Fphi
            ! Er_g = Er + z_2 * k2Er
            ! Ephi_g = Ephi + z_2 * k2Ephi
            ! alpha_m_g = alpha_m + z_2 * k2alpha
            ! ! B_r=Br_f       !----> !NOTE: added this step because I wanted f to be carried to next step
            ! ! B_phi=Bphi_f  !

            ! ! 3rd step
            ! call MTA_nvf(B_r, B_phi, Fr, Fphi,Er, Ephi, alpha_m)
            ! k3r = dt * dBrdt
            ! k3phi = dt * dBphidt
            ! k3Fr = dt * dFrdt
            ! k3Fphi = dt * dFphidt
            ! k3Er = dt * dErdt
            ! k3Ephi = dt * dEphidt
            ! k3alpha = dt * dalpdt

            ! B_r = Br_g + g_3 * k3r
            ! B_phi = Bphi_g + g_3 * k3phi
            ! Fr = Fr_g + g_3 * k3Fr
            ! Fphi = Fphi_g + g_3 * k3Fphi
            ! Er = Er_g + g_3 * k3Er
            ! Ephi = Ephi_g + g_3 * k3Ephi
            ! alpha_m = alpha_m_g + g_3 * k3alpha
            ! ! B_r=Br_f       !----> !NOTE: added this step because I wanted f to be carried to next step
            ! ! B_phi=Bphi_f  !
            ! t = t +  dt
!FOSA:
            call FOSA(B_r, B_phi, alpha_m)
            !STEP 1 !NOTE: have to call spatial derivatives with this new f and g
            k1r = dt * dBrdt
            k1phi = dt * dBphidt
            k1alpha = dt * dalpdt
            B_r = B_r + g_1 * k1r
            B_phi = B_phi + g_1 * k1phi
            alpha_m = alpha_m + g_1 * k1alpha
            Br_g = B_r + z_1 * k1r
            Bphi_g = B_phi + z_1 * k1phi
            alpha_m_g = alpha_m + z_1 * k1alpha
            ! B_r=Br_f       !----> !NOTE: added this step because I wanted f to be carried to next step
            ! B_phi=Bphi_f  ! 
            
            call FOSA(B_r, B_phi, alpha_m)
            k2r = dt * dBrdt
            k2phi = dt * dBphidt
            k2alpha = dt * dalpdt
            B_r = Br_g + g_2 * k2r
            B_phi = Bphi_g + g_2 * k2phi
            alpha_m = alpha_m_g + g_2 * k2alpha
            Br_g = B_r + z_2 * k2r
            Bphi_g = B_phi + z_2 * k2phi
            alpha_m_g = alpha_m + z_2 * k2alpha
            ! B_r=Br_f       !----> !NOTE: added this step because I wanted f to be carried to next step
            ! B_phi=Bphi_f  ! 

            ! 3rd step
            call FOSA(B_r, B_phi, alpha_m)
            k3r = dt * dBrdt
            k3phi = dt * dBphidt
            k3alpha = dt * dalpdt
            B_r = Br_g + g_3 * k3r
            B_phi = Bphi_g + g_3 * k3phi
            alpha_m = alpha_m_g + g_3 * k3alpha
        
            ! B_r=Br_f       !----> !NOTE: added this step because I wanted f to be carried to next step
            ! B_phi=Bphi_f  ! 
            t = t +  dt

            


        end subroutine RK3_implicit
         









    


end module timestepping