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
    

    double precision, dimension(nx) :: k1phi,k1T_torr,k2T_torr,k2phi,k3T_torr,k3phi,k4T_torr,k4phi
    double precision, dimension(nx) :: k1fr, k1fphi, k2fr, k2fphi, k3fr, k3fphi, k4fr, k4fphi
    double precision, dimension(nx) :: k1fz, k2fz, k3fz, k4fz
    double precision, dimension(nx) :: k1alpha,k2alpha,k3alpha,k4alpha
    double precision,dimension(nx) :: g_1, g_2, g_3, z_1, z_2
    contains
        

        subroutine forward_difference
            double precision, dimension(nx) :: temp_Br, temp_Bphi


            do i = nxghost+1,nxphys+nxghost
                temp_Br(i) = B_r(i) - (beta/2)*(alpha_Bphi(i+1)-alpha_Bphi(i-1)) &
                +alp*(B_r(i+1)-2*B_r(i)+B_r(i-1)) - (beta/2)*(Uz_Br(i+1) - Uz_Br(i-1))
                temp_Bphi(i) = B_phi(i) - (beta/2)*(alpha_Br(i+1)-alpha_Br(i-1)) &
                +alp*(B_phi(i+1)-2*B_phi(i)+B_phi(i-1)) - (beta/2)*(Uz_Bphi(i+1) - Uz_Bphi(i-1))&
                +G*B_r(i)

                
            end do
            B_r = temp_Br
            B_phi = temp_Bphi
            t = t + dt
        end subroutine forward_difference

        subroutine backward_difference
            ! character(len=30) :: ghost_zone_type = 'anti-symmetric'

            double precision, dimension(nx) :: temp_Br, temp_Bphi
            ! call impose_boundary_conditions(B_r, ghost_zone_type)
            ! call impose_boundary_conditions(B_phi, ghost_zone_type)
            ! call construct_grid
            ! alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
            ! alpha_Bphi = alpha_cap2*B_phi
            ! alpha_Br = alpha_cap2*B_r
            ! Uz_Br = B_r*U_z_cap
            ! Uz_Bphi = B_phi*U_z_cap

            ! call spatial_derivative(B_r,2,dBr,d2Br)
            ! call spatial_derivative(B_phi,2,dBphi,d2Bphi)
            ! call spatial_derivative(alpha_Br,2,d_alpha_Br,d2_alpha_Br)
            ! call spatial_derivative(alpha_Bphi,2,d_alpha_Bphi,d2_alpha_Bphi)
            ! call spatial_derivative(Uz_Br,2,d_Uz_Br,d2_Uz_Br)
            ! call spatial_derivative(Uz_Bphi,2,d_Uz_Bphi,d2_Uz_Bphi)

            do i = nxghost+1,nxphys+nxghost
                temp_Br(i) = B_r(i) + (beta/2)*(alpha_Bphi(i)-alpha_Bphi(i-2)) &
                -alp*(B_r(i+1)-2*B_r(i)+B_r(i-1)) + (beta/2)*(Uz_Br(i) - Uz_Br(i-2))
                temp_Bphi(i) = B_phi(i) - (beta/2)*(alpha_Br(i)-alpha_Br(i-2)) &
                -alp*(B_phi(i+1)-2*B_phi(i)+B_phi(i-1)) + (beta/2)*(Uz_Bphi(i) - Uz_Bphi(i-2))&
                -G*B_r(i)

            end do
            B_r = temp_Br
            B_phi = temp_Bphi
            t = t + dt

        end subroutine backward_difference

        subroutine central_difference
            double precision, dimension(nx) :: temp_Br, temp_Bphi
            ! character(len=30) :: ghost_zone_type = 'anti-symmetric'

            ! call impose_boundary_conditions(B_r, ghost_zone_type)
            ! call impose_boundary_conditions(B_phi, ghost_zone_type)
            ! call impose_boundary_conditions(old_Br, ghost_zone_type)
            ! call impose_boundary_conditions(old_Bphi, ghost_zone_type)
            ! call construct_grid
            ! alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
            ! alpha_Bphi = alpha_cap2*B_phi
            ! alpha_Br = alpha_cap2*B_r
            ! Uz_Br = B_r*U_z_cap
            ! Uz_Bphi = B_phi*U_z_cap


            do i = nxghost+1,nxphys+nxghost
                temp_Br(i) = (1+2*alp)*((beta)*(alpha_Bphi(i+1)-alpha_Bphi(i-1)) &
                +2*alp*(B_r(i+1)+B_r(i-1)) - (beta)*(Uz_Br(i+1) - Uz_Br(i-1)))+old_Br(i)

                temp_Bphi(i) = (1+2*alp+G)*((beta)*(alpha_Br(i+1)-alpha_Br(i-1)) &
                +2*alp*(B_phi(i+1)+B_phi(i-1)) - (beta)*(Uz_Bphi(i+1) - Uz_Bphi(i-1)))+old_Bphi(i) !NOT_SURE
                

            end do
            B_r = temp_Br
            B_phi = temp_Bphi
            t = t + dt

        end subroutine central_difference

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




        !     call nvf(B_r, B_phi, alpha_m)
        !     k1r = dt*dBrdt
        !     k1phi = dt*dBphidt
        !     k1alpha = dt*dalpdt
        !     B_r = B_r+0.5*k1r
        !     B_phi = B_phi+0.5*k1phi
        !     alpha_m = alpha_m+0.5*k1alpha

        !     call nvf(B_r, B_phi, alpha_m)
        !     k2r = dt*dBrdt
        !     k2phi = dt*dBphidt
        !     k2alpha = dt*dalpdt
        !     B_r = B_r+0.5*k2r
        !     B_phi = B_phi+0.5*k2phi
        !     alpha_m = alpha_m+0.5*k2alpha

        !     call nvf(B_r, B_phi, alpha_m)
        !     k3r = dt*dBrdt
        !     k3phi = dt*dBphidt
        !     k3alpha = dt*dalpdt
        !     B_r = B_r+k3r
        !     B_phi = B_phi+k3phi
        !     alpha_m = alpha_m+k3alpha

        !     call nvf(B_r, B_phi, alpha_m)
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
            double precision, dimension(nx) :: Br_g, Bphi_g, alpha_m_g, phi_g, T_torr_g
            double precision, dimension(nx) :: Fr_g, F_phi_g, Fz_g
            ! double precision, dimension(nx) :: Br_f, Bphi_f, alpha_m_f, phi_f, T_torr_f

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

            call nvf(Phi, T_torr, Fr, F_phi, Fz, alpha_m)
            !STEP 1 !NOTE: have to call spatial derivatives with this new f and g
            k1phi = dt * dphi_dt
            k1T_torr = dt * dTdt
            k1alpha = dt * d_alpha_m_dt
            k1fr = dt * dFr_dt
            k1fphi = dt * dFphi_dt
            k1fz = dt * dFz_dt
            Phi = Phi + g_1 * k1phi
            T_torr = T_torr + g_1 * k1T_torr
            alpha_m = alpha_m + g_1 * k1alpha
            Fr = Fr + g_1 * k1fr
            F_phi = F_phi + g_1 * k1fphi
            Fz = Fz + g_1 * k1fz
            phi_g = Phi + z_1 * k1phi
            T_torr_g = T_torr + z_1 * k1T_torr
            alpha_m_g = alpha_m + z_1 * k1alpha
            Fr_g = Fr + z_1 * k1fr
            F_phi_g = F_phi + z_1 * k1fphi
            Fz_g = Fz + z_1 * k1fz
            ! Phi=phi_f       !----> !NOTE: added this step because I wanted f to be carried to next step
            ! T_torr=T_torr_f  ! 
            
            call nvf(Phi, T_torr, Fr, F_phi, Fz, alpha_m)
            k2phi   = dt * dphi_dt
            k2T_torr = dt * dTdt
            k2alpha = dt * d_alpha_m_dt
            k2fr = dt * dFr_dt
            k2fphi = dt * dFphi_dt
            k2fz = dt * dFz_dt
            Phi = phi_g + g_2 *k2phi
            T_torr = T_torr_g + g_2 *  k2T_torr
            alpha_m = alpha_m_g + g_2 * k2alpha
            Fr = Fr_g + g_2 * k2fr
            F_phi = F_phi_g + g_2 * k2fphi
            Fz = Fz_g + g_2 * k2fz
            phi_g = Phi + z_2 * k2phi
            T_torr_g = T_torr + z_2 * k2T_torr
            alpha_m_g = alpha_m + z_2 * k2alpha
            Fr_g = Fr + z_2 * k2fr
            F_phi_g = F_phi + z_2 * k2fphi
            Fz_g = Fz + z_2 * k2fz
            ! Phi=phi_f       !----> !NOTE: added this step because I wanted f to be carried to next step
            ! T_torr=T_torr_f  ! 

            ! 3rd step
            call nvf(Phi, T_torr, Fr, F_phi,Fz, alpha_m)
            k3T_torr = dt * dphi_dt
            k3phi = dt * dTdt
            k3alpha = dt * d_alpha_m_dt
            k3fr = dt * dFr_dt
            k3fphi = dt * dFphi_dt
            k3fz = dt * dFz_dt
            Phi = phi_g + g_3 * k3phi
            T_torr = T_torr_g + g_3 * k3T_torr
            alpha_m = alpha_m_g + g_3 * k3alpha
            Fr = Fr_g + g_3 * k3fr
            F_phi = F_phi_g + g_3 * k3fphi
            Fz = Fz_g + g_3 * k3fz
        
            ! Phi=phi_f       !----> !NOTE: added this step because I wanted f to be carried to next step
            ! T_torr=T_torr_f  ! 
            t = t +  dt


        end subroutine RK3_implicit
         









    


end module timestepping