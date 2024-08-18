module constants
  implicit none
  !  CONSTANTS
  double precision, parameter :: pi= 3.14156295358
  double precision, parameter :: g_mp=   1.67262158d-24
  double precision, parameter :: cm_kpc= 3.08567758d21
  double precision, parameter :: km_kpc= 3.08567758d16
  double precision, parameter :: cm_km=  1.d5
  double precision, parameter :: s_Gyr=  1.d9*365.25d0*24*3600
  double precision, parameter :: G_muG= 1.d6
end module constants
!*****************************************************
module modules
  implicit none
  !!!!!!!!!!!!!
  !! MODULES !!
  !!!!!!!!!!!!!
  logical :: Damp=       .false.  !Set to 0 for FOSA, 1 for minimal tau approximation
  logical :: Alg_quench= .false.  !Works with dyn_quench=F; Set to 1 for algebraic quenching (alpha= alpha_k/(1+Emag/Beq^2))
  logical :: Dyn_quench= .true.  !Works with alg_quench=F; Set to 1 for dynamical quenching (d_alpha_m/dt eqn incl in sim)
  logical :: Alp_sin=    .true.  !Set to 1 to get a sinusoidal alpha profile; Set to 0 for a linear alpha profile.
  logical :: Alp_squared=.true.  !Set to 1 to include alpha^2 effect; set to 0 to use alpha-omega approximation equations
  logical :: Shear=      .true.  !Set to 1 to include Omega effect in dynamo, 0 for alpha^2 dynamo
end module modules
!*****************************************************
module physical_params
  use constants
  implicit none
!  DIMENSIONLESS UNITS
  double precision, parameter :: etat=1.
  double precision, parameter :: h0=   1.
  double precision, parameter :: B0=  1.
  double precision, parameter :: Bseed=1.d-3  !Amplitude of seed magnetic field, as a fraction of Beq
!  FIDUCIAL DISC PARAMETERS
  double precision, parameter :: r0_kpc=10.  !Fiducial radius
  double precision, parameter :: V0_kms=250.  !Circular rotation speed at r=r0_kpc in km/s
  double precision, parameter :: r_kpc=4.  !Radius at which simulation is performed
!  TURBULENCE
  double precision, parameter :: l_kpc= 0.1  !Used in Krause's law; Size of largest turbulent eddies, in parsecs
  double precision, parameter :: h0_kpc= 0.5  !The designated unit of length, corresp to a half-disc thickness at fiducial radius r=r0
  double precision, parameter :: h0_km=h0_kpc*km_kpc  !The designated unit of length, corresp to a typical half-disk thickness, in km
  double precision, parameter :: v_turb_kms= 10.  !Typical turbulent velocity; together with h0_kpc it determines the unit of time for the simulation
  double precision, parameter :: etat_cm2s= l_kpc*cm_kpc*v_turb_kms*cm_km/3  !Turbulent diffusivity in units of cm^2/s
  double precision, parameter :: td0_Gyr= h0_kpc**2/etat_cm2s/s_Gyr*cm_kpc*cm_kpc  !Typical vertical turbulent diffusion timescale in units of Gyr
  double precision, parameter :: td0_s= td0_Gyr*s_Gyr  !Typical vertical turbulent diffusion timescale in units of seconds
  double precision, parameter :: nH0_cm3=0.1  !Number density in cm^{-3} at r=0 (only used to get physical value for B0)
  double precision, parameter :: rho0_gcm3=nH0_cm3*g_mp  !Density in g/cm^3 (only used to get physical value for B0)
  double precision, parameter :: B0_muG= sqrt(4*pi*rho0_gcm3)*v_turb_kms*cm_km*G_muG
  double precision, parameter :: req_kpc= 20.  !Scale length of equiparition field in kpc
!  VERTICAL VELOCITY
  double precision, parameter :: U0_kms=0.  !Vertical mean velocity in km/s
  !  SSSB06 uses U0=0.3-3 corresponding to 0.2-2 km/s; CSS12 uses R_Uz=0.45
!   SCALE HEIGHT
  double precision, parameter :: r_D_kpc= 10.  !Relevant only if Flaring=1; characteristic radius of the hyperbolically varying scale height in kpc
  double precision, parameter :: h_D_kpc= h0_kpc/(1.+(r0_kpc/r_D_kpc)**2)**(1./2)  !Relevant only if Flaring=1; amplitude of the scale height in kpc
  double precision, parameter :: h_kpc= h_D_kpc*(1.+(r_kpc/r_D_kpc)**2)**(1./2)  !Disk scale height at r_kpc where simulation is performed
  double precision, parameter :: td_Gyr= h_kpc**2/etat_cm2s/s_Gyr*cm_kpc*cm_kpc  !Typical vertical turbulent diffusion timescale in units of Gyr
  double precision, parameter :: td_s= td_Gyr*s_Gyr  !Typical vertical turbulent diffusion timescale in units of seconds
!  ROTATION CURVE
  double precision, parameter :: r_om_kpc= 2.  !Relevant only if Om_Brandt =1; r_om is the characteristic radius of the Brandt profile in kpc
  double precision, parameter :: om0_kmskpc= (1.+(r0_kpc/r_om_kpc)**2)**(1./2)*V0_kms/r0_kpc  !om0 of Brandt profile in physical units
  double precision, parameter :: om_kmskpc= om0_kmskpc/(1. +(r_kpc/r_om_kpc)**2)**(1./2)
  double precision, parameter :: G_kmskpc=-om0_kmskpc*(r_kpc/r_om_kpc)**2/(1. +(r_kpc/r_om_kpc)**2)**(3./2)
!  ALPHA EFFECT
  double precision, parameter :: C_alp= 1.
  double precision, parameter :: alp_k_mean_kms= C_alp*l_kpc**2/h_kpc*om_kmskpc  !maximum alpha_k in km/s
!
!  DIMENSIONLESS PARAMETERS
  double precision, parameter :: l= l_kpc/h0_kpc*h0
  double precision, parameter :: h= h_kpc/h0_kpc*h0
  double precision, parameter :: tau=(l/h0)**2/3  !tau in units of td=h^2/etat
  double precision, parameter :: tautilde=1.  !Ratio of tau (eddy turnover time) to correlation time of the turbulence
  double precision, parameter :: U0=U0_kms/h0_km*td0_s  !Vertical mean velocity in normalized units h0/td0; 
  double precision, parameter :: G=G_kmskpc/km_kpc*td0_s
!  double precision, parameter :: R_omega= -18.75*25./26  !R_omega=Gh^2/etat
  double precision, parameter :: R_omega= G*h**2/etat  !R_omega=Gh^2/etat
  double precision :: tau_mta=0.
!  -18.75 corresponds to 25 km/s/kpc used in SSSB06. Value in Chamandy et al 2012 (paper I) is -18.75*25./26. 
!  double precision, parameter :: R_alpha= 1.5  
  double precision, parameter :: alp_k_mean= alp_k_mean_kms/h0_km*td0_s
  double precision, parameter :: R_alpha= alp_k_mean*h/etat
!  R_alpha=alp0*h/etat; Ampl of alpha effect, in normalized units h/td; 0.75 corresponds to 0.5 km/s used in SSSB06
  double precision, parameter :: Dyn=R_omega*R_alpha
  integer :: n_alp=  1  !Relevant only if Module Alp_sin eq 1. Sets the rate of variation of alpha in z
  double precision :: factau= 1.0 !1.  !Multiplication factor: tau_mta=factau*tau
  double precision :: Rm_inv= 0.  !1.e-5	;Inverse magnetic Reynolds number
  !  DYNAMICAL QUENCHING
  !  New Vishniac flux
  double precision, parameter :: fac_NV=0.  !fudge factor (f in Sharanya's notes) in the new Vishniac flux <b^2>=|fac|*Beq^2; 
  !  Set to 0. for no New Vishniac flux. May be +/-.
  !  Fickian diffusive flux
  double precision, parameter :: kappa=0.3  !Fickian diffusion coefficient
!
contains
  subroutine mta
    use modules
!  CLOSURE APPROXIMATION
    if (.not.Damp) then
      tau_mta=0.
    else
      tau_mta=factau*tau  !Damping time of the minimal tau approximation
    endif
  end subroutine mta
end module physical_params
!*****************************************************
module ts_params
!
  use physical_params
!
  implicit none
  !  NUMERICS
  integer, parameter :: n1= 800
  integer, parameter :: n2= 1000  !long(1./dt)+1;	;Replot n1 times (after n2 timesteps)
  double precision, parameter :: tsnap= 0.025/td0_Gyr !0.1d0  !Time between successive snapshots
  double precision, parameter :: dt= tsnap/n2 !0.0339/30 !0.0339/3500 !1.d-3!1.d-4!5.d-4  !Timestep in units of td=h^2/C_etat
  double precision :: t=0.
  double precision :: first=0.  !for Runge-Kutta routine
end module ts_params
!*****************************************************
module grid_params
!
  use physical_params
!
  implicit none
!
  integer, parameter :: nxphys= 201  !Resolution in z
  integer, parameter :: nxghost= 3  !Number of ghost cells at each end in z
  integer, parameter :: nx= nxphys +2*nxghost  !Resolution in z
  double precision, dimension(nx) :: x
  double precision :: dx
end module grid_params
!*****************************************************
module grid
!
  implicit none
!
contains
  subroutine construct_grid
    use grid_params
!
    integer :: i
    double precision, parameter :: len= 2.*h0  !Simulation domain (in units of h0)
    double precision, dimension(nx) :: spac
!
    dx=len/(nxphys-1)  !x corresponds to z coordinate
    do i=1,nx
      x(i)= -(h0 +nxghost*dx) +(i-1)*dx
    enddo
  endsubroutine construct_grid
!
end module grid
!*****************************************************
module var
! Variables
!
  use modules
!
  implicit none
  integer :: nvar
!
  contains
    subroutine init_var
    !  DEFINE ARRAYS FOR PHYSICAL VARIABLES
      if (.not.Dyn_quench) then
        if (.not.Damp) then
          nvar=2  !vars are Br, Bp
        else
          nvar=4  !vars are Br, Bp, Fr, Fp
        endif
      else
        if (.not.Damp) then
          nvar=3  !vars are Br, Bp, alp_m
        else
          nvar=7  !vars are Br, Bp, Fr, Fp, Er, Ep, alp_m
        endif
      endif
    end subroutine init_var
end module var
!*****************************************************
module ts_arrays
!
  use modules
  use grid_params
  use ts_params
  use var
!
  implicit none
!
  double precision, dimension(n1+1) :: ts_t  !time array
  double precision, dimension(n1+1,nx) :: ts_Br, ts_Bp, ts_alp_m, ts_alp_k  !Br array, Bp array, alp_m, alp_k array
!
  contains
    subroutine make_ts_arrays(it,t,f,alp_k)
!
      integer, intent(in) :: it
      double precision, intent(in) :: t
      double precision, dimension(nx,nvar), intent(in) :: f
      double precision, dimension(nx), intent(in) :: alp_k
!
      ts_t(it+1)=       t
      ts_Br(it+1,:)=    f(:,1)
      ts_Bp(it+1,:)=    f(:,2)
      if (Dyn_quench) then   
        if (.not.Damp) then
          ts_alp_m(it+1,:)=   f(:,3)
        else
          ts_alp_m(it+1,:)=   f(:,7)
        endif
      endif
      ts_alp_k(it+1,:)=    alp_k(:)
    end subroutine make_ts_arrays
end module ts_arrays
!*****************************************************
module Uz_profile
!  VERTICAL VELOCITY PROFILE
  use constants
  use modules
  use physical_params
  use grid_params
  use grid
!
  implicit none
!
  double precision, dimension(nx) :: Uz, dUzdz 
! 
  contains
    subroutine construct_Uz_profile
      Uz=U0*x/h  !Vertical mean velocity in normalized units h/td
      dUzdz=U0/h  !Vertical mean velocity in normalized units h/td
    end subroutine construct_Uz_profile
!
  end module Uz_profile
!*****************************************************
module alp_profile
!  ALPHA PROFILE      =Beq_muG,
  use constants
  use modules
  use physical_params
  use grid_params
  use grid
!
  implicit none
!
  double precision, dimension(nx) :: alp_k, dalp_kdz 
! 
contains
  subroutine construct_alp_profile
    if (Alp_sin) then
      alp_k=pi/2*alp_k_mean*sin(n_alp*pi*x/h)
      dalp_kdz=n_alp*pi/h*pi/2*alp_k_mean*cos(n_alp*pi*x/h)
    else
      alp_k=pi/2*alp_k_mean*x/h
      dalp_kdz=pi/2*alp_k_mean/h
    endif
  end subroutine construct_alp_profile
!
end module alp_profile
!*****************************************************
module Beq_profile
!  Beq PROFILE
!
  use physical_params
  use grid_params
  use grid
!
  implicit none
!
  double precision, dimension(nx) :: Beq, Beq_muG
!
  contains
    subroutine construct_Beq_profile
      Beq= B0*exp(-r_kpc/req_kpc)*exp(-x**2/2/h**2)
      Beq_muG= Beq*B0_muG/B0
    end subroutine construct_Beq_profile
    !
end module Beq_profile
!*****************************************************
module boundary_conditions
!  BOUNDARY CONDITIONS
!
  use modules
  use grid_params
  use var
!  
  implicit none
!
  contains
    subroutine impose_bc(f)
!
      double precision, dimension(nx,nvar), intent(inout) :: f
      integer :: ix
!
      if (nxghost/=0) then
        do ix=1,nxghost
!          f(ix     ,1:2)=  f(2*(nxghost+1)-ix     ,1:2)  !Symmetric     about z=-h    Neumann   BC on Br, Bp: dBrdz=dBpdz=0 at z=-h
!          f(nx+1-ix,1:2)=  f(nx+1-2*(nxghost+1)+ix,1:2)  !Symmetric     about z= h    Neumann   BC on Br, Bp: dBrdz=dBpdz=0 at z= h
            f(ix     ,1:2)= -f(2*(nxghost+1)-ix     ,1:2)  !Antisymmetric about z=-h    Dirichlet BC on Br, Bp
            f(nx+1-ix,1:2)= -f(nx+1-2*(nxghost+1)+ix,1:2)  !Antisymmetric about z= h    Dirichlet BC on Br, Bp
          if (Dyn_quench) then
            f(ix     ,nvar)=  2*f(nxghost+1       ,nvar) -f(2*(nxghost+1)-ix     ,nvar)  !Relative antisymmetric   (Specify alpha_m in ghost zones)
            f(nx+1-ix,nvar)=  2*f(nx+1-(nxghost+1),nvar) -f(nx+1-2*(nxghost+1)+ix,nvar)  !Relative antisymmetric   (Specify alpha_m in ghost zones)
          endif
        enddo
      endif
    end subroutine impose_bc
end module boundary_conditions  
  
!*****************************************************
module seed
! Seed field
!
  use physical_params
  use grid_params
  use var
!  
  implicit none
!
  double precision, dimension(nx) :: xseed
contains
  subroutine init_random_seed
!
    call random_number(xseed)
!
  end subroutine init_random_seed
end module seed
!*****************************************************
module deriv
!
  use grid_params
  use grid
!
contains
  function xder(f)
!
    implicit none
!
    integer :: ix
    double precision, dimension(nx), intent(in) :: f
    double precision :: fac
    double precision, dimension(nx) :: xder
!
!  assume uniform mesh
!
!
    fac=1./(60.*(x(2)-x(1)))
!
    do ix=4,nx-3
      xder(ix)=fac*(+45.*(f(ix+1)-f(ix-1)) & 
                    - 9.*(f(ix+2)-f(ix-2)) &
                    +    (f(ix+3)-f(ix-3)) &
                   )
    enddo
!
!  The following does the same thing but is slightly slower
!
!    xder=fac*(+45.*(cshift(f,1)-cshift(f,-1)) & 
!              - 9.*(cshift(f,2)-cshift(f,-2)) & 
!              +    (cshift(f,3)-cshift(f,-3)) & 
!             )
!
!      do boundary points (also 6th order)
!
    xder(1)=fac*(-147.*f(1)+360.*f(2)-450.*f(3)+400.*f(4)-225.*f(5)+72.*f(6)-10.*f(7))
    xder(2)=fac*( -10.*f(1)- 77.*f(2)+150.*f(3)-100.*f(4)+ 50.*f(5)-15.*f(6)+ 2.*f(7))
    xder(3)=fac*(   2.*f(1)- 24.*f(2)- 35.*f(3)+ 80.*f(4)- 30.*f(5)+ 8.*f(6)-    f(7))
!
!      outer points
!
    xder(nx  )=fac*(147.*f(nx)-360.*f(nx-1)+450.*f(nx-2)-400.*f(nx-3)+225.*f(nx-4)-72.*f(nx-5)+10.*f(nx-6))
    xder(nx-1)=fac*( 10.*f(nx)+ 77.*f(nx-1)-150.*f(nx-2)+100.*f(nx-3)- 50.*f(nx-4)+15.*f(nx-5)- 2.*f(nx-6))
    xder(nx-2)=fac*( -2.*f(nx)+ 24.*f(nx-1)+ 35.*f(nx-2)- 80.*f(nx-3)+ 30.*f(nx-4)- 8.*f(nx-5)+    f(nx-6))
!
  end function xder
!
  function xder2(f)
!
    implicit none
!
    integer :: ix
    double precision, dimension(nx), intent(in) :: f
    double precision :: fac
    double precision, dimension(nx) :: xder2
!
!  assume uniform mesh
!
    fac=1./(180.*(x(2)-x(1))**2)
!
    do ix=4,nx-3
      xder2(ix)=fac*(-490.*f(ix)         &
                 +270.*(f(ix+1)+f(ix-1)) &
                 - 27.*(f(ix+2)+f(ix-2)) &
                 +  2.*(f(ix+3)+f(ix-3)) &
                )
    enddo
!
!  The following does the same thing but is slightly slower
!
!    xder2=fac*(-490.*f                          &
!               +270.*(cshift(f,1)+cshift(f,-1)) &
!               - 27.*(cshift(f,2)+cshift(f,-2)) &
!               +  2.*(cshift(f,3)+cshift(f,-3)) &
!              )
!
!  do boundary points (also 6th order)
!
    xder2(1)=fac*(812.*f(1)-3132.*f(2)+5265.*f(3)-5080.*f(4)+2970.*f(5)-972.*f(6)+137.*f(7))
    xder2(2)=fac*(137.*f(1)- 147.*f(2)- 255.*f(3)+ 470.*f(4)- 285.*f(5)+93. *f(6)- 13.*f(7))
    xder2(3)=fac*(-13.*f(1)+ 228.*f(2)- 420.*f(3)+ 200.*f(4)+  15.*f(5)-12. *f(6)+  2.*f(7))
!
!  outer points
!
    xder2(nx  )=fac*(812.*f(nx)-3132.*f(nx-1)+5265.*f(nx-2)-5080.*f(nx-3)+2970.*f(nx-4)-972.*f(nx-5)+137.*f(nx-6))
    xder2(nx-1)=fac*(137.*f(nx)- 147.*f(nx-1)- 255.*f(nx-2)+ 470.*f(nx-3)- 285.*f(nx-4)+ 93.*f(nx-5)- 13.*f(nx-6))
    xder2(nx-2)=fac*(-13.*f(nx)+ 228.*f(nx-1)- 420.*f(nx-2)+ 200.*f(nx-3)+  15.*f(nx-4)- 12.*f(nx-5)+  2.*f(nx-6))
!
  end function xder2
!
end module deriv
!*****************************************************
module equ
!
  use constants
  use physical_params
  use ts_params
  use grid_params
  use grid
  use modules
  use Uz_profile
  use alp_profile
  use Beq_profile
  use boundary_conditions
  use seed
  use var
  use deriv
!  
  implicit none
!
  double precision, dimension(nx) :: Br, Bp, Fr, Fp, Er, Ep, dBrdz, d2Brdz2, dBpdz, d2Bpdz2
  double precision, dimension(nx) :: alp_m, alp, dalp_mdz, d2alp_mdz2, dalpdz
  double precision, dimension(nx) :: Emag, DivVishniac
!
contains
  subroutine pde(f,dfdt)
!
  double precision, dimension(nx,nvar) :: f, dfdt
!  intent(in) :: f
  intent(inout) :: f
  intent(out) :: dfdt
!
Br= f(:,1)
Bp= f(:,2)
!
if (Damp) then
  Fr=f(:,3)
  Fp=f(:,4)
  if (Dyn_quench) then
    Er=f(:,5)
    Ep=f(:,6)
    alp_m=f(:,7)
  endif
else
  Fr=0.
  Fp=0.
  if (Dyn_quench) then
    alp_m=f(:,3)
  endif
endif
!
dBrdz=  xder(Br)
d2Brdz2=xder2(Br)
dBpdz=  xder(Bp)
d2Bpdz2=xder2(Bp)
!
! KINETIC ALPHA PROFILE
!
DivVishniac = -fac_NV*R_omega*(2.*l/3)**2*x
!
! CALCULATE Emag
!
Emag= Br**2 +Bp**2
!
      call construct_Uz_profile
!
      call construct_Beq_profile
!
      call construct_alp_profile
!
!      call construct_etat_profile
!
      if (Alg_quench) then
        alp= alp_k/(1. +Emag/Beq**2)  !Formula for simple alpha quenching
        dalpdz=   xder(alp)
      elseif (Dyn_quench) then
        dalp_mdz=   xder(alp_m)
        d2alp_mdz2= xder2(alp_m)
        alp= alp_k +alp_m  !Total alpha is equal to the sum of kinetic and magnetic parts
        dalpdz=     dalp_kdz +dalp_mdz
      else
        alp= alp_k
        dalpdz=   dalp_kdz
      endif
!
! LIST OF VARIABLE NAMES FOR f ARRAY
!
! UNDER FOSA				UNDER TAU APPROXIMATION
! f(:,:,1)=Br				f(:,:,1)=Br
! f(:,:,2)=Bp				f(:,:,2)=Bp
! f(:,:,3)=alp_m			f(:,:,3)=Fr
! 					f(:,:,4)=Fp
! 					f(:,:,5)=Er
!      	  			        f(:,:,6)=Ep
!      				        f(:,:,7)=alp_m
!
  call impose_bc(f)
!
! METHOD: ALL ARRAYS ARE 2-DIMENSIONAL
!
!mark
if (.not.Damp) then
!
!  CASE 1: FOSA (tau_mta-->0 LIMIT)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
!
  dfdt(:,1) =             -dUzdz*Br -Uz*dBrdz -dalpdz*Bp -alp*dBpdz +etat*d2Brdz2 +Rm_inv*etat*d2Brdz2
  dfdt(:,2) = +R_omega*Br -dUzdz*Bp -Uz*dBpdz +dalpdz*Br +alp*dBrdz +etat*d2Bpdz2 +Rm_inv*etat*d2Bpdz2
  if (.not.Alp_squared) then
    dfdt(:,2) = +R_omega*Br -dUzdz*Bp -Uz*dBpdz                     +etat*d2Bpdz2 -Rm_inv*d2Bpdz2
  endif
  if (Dyn_quench) then
    dfdt(:,3) = -2./3/tau*( alp*(Br**2 +Bp**2)/Beq**2 -etat*(Bp*dBrdz-Br*dBpdz)/Beq**2 +Rm_inv*alp_m +DivVishniac) &
                -dalp_mdz*Uz -alp_m*dUzdz +kappa*etat*d2alp_mdz2
  endif
!
else
!
!  CASE 2: MTA (FINITE tau_mta)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
!
  dfdt(:,1) =             -dUzdz*Br -Uz*dBrdz +Fr +Rm_inv*etat*d2Brdz2
  dfdt(:,2) = +R_omega*Br -dUzdz*Bp -Uz*dBpdz +Fp +Rm_inv*etat*d2Bpdz2
  dfdt(:,3) = (-dalpdz*Bp -alp*dBpdz +etat*d2Brdz2 -Fr)/tau_mta
  dfdt(:,4) = (+dalpdz*Br +alp*dBrdz +etat*d2Bpdz2 -Fp)/tau_mta
  if (.not.Alp_squared) then
    dfdt(:,4) = (                    +etat*d2Bpdz2 -Fp)/tau_mta
  endif
  if (Dyn_quench) then
    dfdt(:,5) = (+alp*Br+etat*dBpdz -Er)/tau_mta
    dfdt(:,6) = (+alp*Bp-etat*dBrdz -Ep)/tau_mta
    dfdt(:,7) = -2./3/tau*( (Er*Br +Ep*Bp)/Beq**2 +Rm_inv*alp_m +DivVishniac) &
                -dalp_mdz*Uz -alp_m*dUzdz +kappa*etat*d2alp_mdz2
  endif
!
endif
!
! IMPOSE BOUNDARY CONDITIONS (IF INITIAL VALUE IS ZERO AND dfdt=0 THIS IS EQUIVALENT TO SETTING f=0 AT ALL TIMES)
!
!dfdt(1 ,1:2)=0.d0 !Ensure that Br, Bp remain zero at z=-h
!dfdt(nx,1:2)=0.d0 !Ensure that Br, Bp remain zero at z= h
  end subroutine pde
end module equ
!*****************************************************
module timestep 
!
  use var
  use grid_params
  use ts_params
  use equ
!
contains
  subroutine rk(f)
!
  implicit none
!
  double precision :: gam1, gam2, gam3, zet1, zet2
  double precision, dimension(nx,nvar) :: f, dfdt, pdef, ftmp
  double precision :: ttmp
!
  intent(inout) :: f
!
!  Runge Kutta 3rd order time advance
!  f = f(exact) - 0.0046*dt^3*d3f/dt3
!
    gam1=8./15. !gam1, gam2 AND gam3 ARE THE COEFFICIENTS OF THE TIMESTEPS AT WHICH dfdt IS CALCULATED
    gam2=5./12.
    gam3=3./4.
    zet1=-17./60.
    zet2=-5./12.
!
    call pde(f,dfdt)
    pdef=dfdt !pde FUNCTION CALCULATES VALUE OF TIME DERIVATIVE ACCORDING TO P.D.E.s
    f=f+dt*gam1*pdef !FIRST GET INTERMEDIATE VALUE OF f (AT t=t+dt*gam1) USING dfdt AT t_0
    t=t+dt*gam1 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet1*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1) USING dfdt AT t_0
    ttmp=t+dt*zet1 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1)
!
    call pde(f,dfdt)
    pdef=dfdt !NOW CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1
    f=ftmp+dt*gam2*pdef !USE THIS TO GET ANOTHER INTERMEDIATE VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2)
    t=ttmp+dt*gam2 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet2*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2) USING dfdt AT t=t_0+dt*gam1
    ttmp=t+dt*zet2 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2)
!
    call pde(f,dfdt)
    pdef=dfdt !CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1+dt*zet1+dt*gam2
    f=ftmp+dt*gam3*pdef !USE THIS TO GET THE FINAL VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2+dt*gam3)
    t=ttmp+dt*gam3 !THEN GO TO THAT TIMESTEP
!
    first=first+1. !COUNTS THE NUMBER OF TIMES RUNGA-KUTTA ROUTINE IS EXCUTED
  end subroutine rk
end module timestep
!*****************************************************
module start
!
  use constants
  use physical_params
  use ts_params
  use grid_params
  use grid
  use modules
  use Uz_profile
  use alp_profile
  use Beq_profile
  use seed
  use var
!
  implicit none
  double precision, allocatable, dimension(:,:) :: f, dfdt
!
contains
  subroutine init_start
!
! INITIALIZE VARIABLE ARRAY
!  
    call init_var
!
    allocate(f(nx,nvar))
    allocate(dfdt(nx,nvar))
!
! SET UP NUMERICS AND PLOTTING
!
    call construct_grid
!
!  VERTICAL MEAN VELOCITY PROFILE
!
    call construct_Uz_profile
!
!  ALPHA PROFILE
!
    call construct_alp_profile
!
!  BEQ PROFILE
!
    call construct_Beq_profile
!
!  INITIALIZE VARIABLES
!
    call init_var
!
!  DETERMINE WHETHER TO USE FOSA OR MINIMAL TAU APPROXIMATION
!
    call mta
!
!  SEED FIELD
!
    call init_random_seed
!
    f(:,:)=0.
    f(:,1)=Bseed*exp(-x**2/h**2)*(1.-x**2/h**2)
!
!default,seed,seed_id
!f(*,0)=0.;Bseed*gaunoise(nx,seed=seed)
!f(*,1)=Bseed*gaunoise(nx,seed=seed)
!Br_init=f[*,0]
!  f(:,:)=Bseed  !seed field
    f(1, 1:2)=0.  !boundary condition 1 Br=Bp=0 at z=-1
    f(nx,1:2)=0.  !boundary condition 2 Br=Bp=0 at z=+1
!  
!  PRINT INFO TO SCREEN
!
!mark
    print*,'UNITS'
    print*,'td0_Gyr=',td0_Gyr
    print*,'h0_kpc=',h0_kpc
    print*,'etat_cm2s=',etat_cm2s
    print*,''
    print*,'NUMERICS:'
    print*,'dt=',dt
    print*,'nvar=',nvar
    print*,'dx=',dx
    print*,'nx=',nx
    print*,'n1=',n1
    print*,'n2=',n2
    print*,'xmin, xmax',minval(x),maxval(x)
    print*,''
    print*,'MODULES:'
    print*,'Damp= ',Damp
    print*,'Alp_squared=',Alp_squared
    print*,''
    print*,'PARAMETERS:'
    print*,'Bseed=   ',Bseed
    print*,'tau_mta= ',tau_mta
    print*,'h=   ',h
    print*,'etat=',etat
    print*,'R_alpha=alpha*h0/etat=',R_alpha
    print*,'R_omega=G*h0^2/etat =',R_omega
    print*,'Dynamo number',Dyn
    print*,'G_kmskpc =',G_kmskpc
    print*,'G =',G
    print*,'U0_kms =',U0_kms
    print*,'U0 =',U0
    print*,'B0 =',B0
    print*,'B0_muG =',B0_muG
    print*,'Beq((nx+1)/2)=',Beq((nx+1)/2)
    print*,'Beq_muG((nx+1)/2)=',Beq_muG((nx+1)/2)
  !print*,'Directory ending   =',s0
!
!  WRITE INFO TO FILE
!
    print*,''
    print*,'Writing parameters to file param.out'
    open(10,file='/home/test/fortran_pde/1D/telegraph/z/data/AQ_012/param.out',status="replace")
    write(10,*)n1, n2, dt, tsnap, nx, nxphys, nxghost, nvar, tau, B0, td0_Gyr, h0, h, l, U0, &
               B0_muG, td_Gyr, h0_kpc, h_kpc, l_kpc, U0_kms, v_turb_kms
    close(10)
    print*,'Finished writing parameters to file param.out'
    print*,'Writing initial data to file init.out'
    open(11,file='/home/test/fortran_pde/1D/telegraph/z/data/AQ_012/init.out',status="replace")
    write(11,*)x
    write(11,*)Uz
    write(11,*)alp_k
    write(11,*)Beq
    write(11,*)f(:,1)
    write(11,*)f(:,2)
    close(11)
  end subroutine init_start
  !
end module start
!*****************************************************
program run
  use constants
  use physical_params
  use ts_params
  use ts_arrays
  use grid_params
  use grid
  use modules
  use Uz_profile
  use alp_profile
  use Beq_profile
  use seed
  use var
  use start
  use equ
  use timestep
!
  implicit none 
!
integer :: it,jt
double precision :: cpu_time_start, cpu_time_finish
!
  call cpu_time(cpu_time_start)
!
  call init_var
!
  call init_start
!
  call make_ts_arrays(0,t,f,alp_k)
!
  do it=1,n1
    do jt=1,n2
      call rk(f) !(f,t,dt,first)
      if (isnan(f((nx+1)/2,1))) stop 'NANs detected, change time step'
    enddo
  call make_ts_arrays(it,t,f,alp_k)
  print*,'it=',it,'t=',t,'Br=',f((nx+1)/2,1),'Bp=',f((nx+1)/2,2),'p_B(deg)=',180./pi*atan(f((nx+1)/2,1)/f((nx+1)/2,2))
  print*,'t(Gyr)=',t/td0_Gyr,'Br/B0=',f((nx+1)/2,1),'Bp/B0=',f((nx+1)/2,2),'B/B0=',sqrt(Emag((nx+1)/2))/B0
  enddo    
!
    call impose_bc(f)
!
!  WRITE INFO TO FILE
!
  print*,'Writing data for final timestep to file run.out'
  open(12,file='/home/test/fortran_pde/1D/telegraph/z/data/AQ_012/run.out',status="replace")
  write(12,*)t
  write(12,*)f(:,1)
  write(12,*)f(:,2)
  if (Dyn_quench) then
    if (.not.Damp) then
      write(12,*)f(:,3)
    else
      write(12,*)f(:,7)
    endif
  endif
  close(12)
  print*,'Writing time series data to file ts.out'
  open(13,file='/home/test/fortran_pde/1D/telegraph/z/data/AQ_012/ts.out',status="replace")
  write(13,*)ts_t
  write(13,*)ts_Br
  write(13,*)ts_Bp
  if (Dyn_quench) then
    write(13,*)ts_alp_m
  endif
  write(13,*)ts_alp_k
  close(13)
  call cpu_time(cpu_time_finish)
  print *, "simulation time in minutes: ", (cpu_time_finish -cpu_time_start)/60
end program run
