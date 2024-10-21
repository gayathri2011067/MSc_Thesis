COMMON cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0
COMMON ctimest,first
@common
;
;  COMPILE ESSENTIAL ROUTINES
;
@xder_6th_bc
@xder2_6th_bc
@pde
;
;  COLOUR SCHEME
loadct,39
;
;  CONSTANTS
mp_g=1.67d-24
kpc_cm=3.1d21
km_cm=1.d5
Gyr_s=1.d9*365.25d0*24*3600
;
;  DIMENSIONLESS UNITS
etat0=1.
h=1.
B0=1.
;
;  RUN NAME
s0='NewVishniacFlux_115'
;
;  DIRECTORY NAMES
program_direc=strjoin(['/Users/luke/pde/1D/telegraph/z/program/',s0])
data_direc=strjoin(['/Users/luke/pde/1D/telegraph/z/data/',s0])
plots_direc=strjoin(['/Users/luke/pde/1D/telegraph/z/plots/',s0])
;
;  NUMERICS
len=2.		;Simulation domain (in units of h)
dt=1.d-5	;Timestep in units of td=h^2/etat0
n1=1010 & n2=fix(0.01d0/dt);long(1./dt)+1;	;Replot n1 times (after n2 timesteps)
print,'n1,n2',n1,n2
nx=201		;Resolution in z
ti_alp_mod=1000.;15.;5.	;Time after which Ralpha=Ralpha_mod
tf_alp_mod=1000.;30.	;Time after which Ralpha returns to its original value
ncurve=1		;Number of curves to plot (=max value of iDyn+1)
Rom=fltarr(ncurve)
;
;  SEED
seed_id=123	;Arbitrary input code for random number generator
Bseed=1.e-2	;Amplitude of seed magnetic field, as a fraction of Beq
;
;  DIMENSIONAL PARAMETERS
l0_kpc=0.1	;Used in Krause's law; Size of largest turbulent eddies, in parsecs
h_kpc=0.5	;The designated unit of length, corresponding to a typical half-disc thickness, in parsecs; determines time unit (vertical diffusion time)
v_turb0_kms=10.	;Typical turbulent velocity; together with h_kpc it determines the unit of time for the simulation (typical vertical diffusion time)
n_ion_cm3=0.5	;number density of ions in cm^-3
n_neu_cm3=0.5	;number density of neutrals in cm^-3
mass_g=mp_g	;mass of ion/neutral particle in grams
B0_Gauss=sqrt(4*!pi*(n_ion_cm3+n_neu_cm3)*mass_g*(v_turb0_kms*km_cm)^2)		;Equipartition field in units of Gauss
etat0_cm2s=l0_kpc*kpc_cm*v_turb0_kms*km_cm/3				;Turbulent diffusivity in units of cm^2/s
td0_Gyr=h_kpc^2*kpc_cm^2/etat0_cm2s/Gyr_s					;Vertical turbulent diffusion timescale in units of years
;
;  DIMENSIONLESS PARAMETERS
l0=l0_kpc/h_kpc*h
v_turb0= v_turb0_kms/etat0_cm2s*etat0*h_kpc/h*kpc_cm*km_cm
tau=(l0/h)^2/3			;tau in units of td=h^2/etat
;R_omega=-18.75	;R_omega=Gh^2/etat; Shear in vertical turbulent diffusion times; -18.75 corresponds to 25 km/s/kpc used in SSSB06
;Rom[0]=-50. ;& Rom[1]=-30. & Rom[2]=-40. & Rom[3]=-50.	;fac = -0.1, 0.001
;Rom[0]=-30. & Rom[1]=-40. & Rom[2]=-50. 			;fac =  0.1
Rom[0]=-20. ;& Rom[1]=-30. 			;fac =  0.1, alpha_m, Br and Bp plots
R_omega=Rom[0];-18.75*25./26	;Value in Chamandy et al 2012 (paper I) is -18.75*25./26
R_alpha=0.;1.0;1.5;0.75		;R_alpha=alp0*h/etat; Amplitude of alpha effect, in normalized units h/td; 0.75 corresponds to 0.5 km/s used in SSSB06
n_alp=1;2			;Relevant only if Module Alp_sin eq 1. Sets the rate of variation of alpha in z
factau=1.		;Multiplication factor: tau_mta=factau*tau
Rm_inv=0.		;1.e-5	;Inverse magnetic Reynolds number
ion_frac=n_ion_cm3/(n_ion_cm3 +n_neu_cm3)
;
;;;;;;;;;;;;;
;; MODULES ;;
;;;;;;;;;;;;;
Damp=0 			;Set to 0 for FOSA, 1 for minimal tau approximation
Dyn_quench=1		;Works with Alp_quench=0; Set to 1 for dynamical quenching (d_alpha_m/dt equation included in simulation)
Alp_sin=1		;Set to 1 to get a sinusoidal alpha profile; Set to 0 for a linear alpha profile.
Alp_squared=1		;Set to 1 to include alpha^2 effect; set to 0 to use alpha-omega approximation equations
Shear=1			;Set to 1 to include Omega effect in dynamo, 0 for alpha^2 dynamo
;
;  DYNAMICAL QUENCHING
;
;  Advective flux
U0=0.45;0.1		;Vertical mean velocity in normalized units h/td; SSSB06 uses U0=0.3-3 corresponding to 0.2-2 km/s; CSS12 uses R_Uz=0.45
;  New Vishniac flux
fac=1.;0.		;fudge factor (f in Sharanya's notes) in the new Vishniac flux <b^2>=|fac|*Beq^2; Set to 0. for no New Vishniac flux. May be +/-.
;  Fickian diffusive flux
kappa=0.;0.3	;Fickian diffusion coefficient
;  Ambipolar drift flux
C_cm2sGauss2=0; 1.5e10*(n_ion_cm3+n_neu_cm3)^(-2)*(ion_frac*(1.-ion_frac)/0.25)^(-1)*kpc_cm	;Ambipolar drift const V=C*JxB in cgs units cm^2/s/Gauss^2
C_ad=0;C_cm2sGauss2/h_kpc^2/kpc_cm^2*td_Gyr*Gyr_s*Beq_Gauss^2			;Ambipolar drift constant C_ad, where V=C_ad*JxB, in dimensionless units
fac2=0; 10.
;
;mark2
;;;;;;;;;;;;;;;;;;;;;;;;;
;; PHYSICAL PARAMETERS ;;
;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  TURBULENCE	;Aside from determining dimensionfull values, the following are relevant only if Dec_alp=1. Then they affect the dynamo number.
v_turb_M31=8.
v_turb_M51=10.
v_turb_NGC6946=15.
;
;  CLOSURE APPROXIMATION
if (Damp eq 0) then begin
  tau_mta=0.
endif else begin
  tau_mta=factau*tau		;Damping time of the minimal tau approximation
  tautilde=1.		;Ratio of tau (eddy turnover time) to correlation time of the turbulence (which may be determined by frequency of SN explosions)
endelse
;
;  FIDUCIAL DISC PARAMETERS
r_0_kpc=10.	;Fiducial radius at which values of h_0, R_alpha_0, R_omega_0, etc. are calculated in kpc
V_0_kms=250.	;Circular rotation speed at r=r_0 in km/s
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; SET UP NUMERICS AND PLOTTING ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
dx=len/(nx-1)			;x corresponds to z coordinate
x=-1.+findgen(nx)*dx
;
if (Dyn_quench eq 0) then begin
  if (Damp eq 0) then begin
    nvar=2		;vars are Br, Bp
  endif else begin
    nvar=4		;vars are Br, Bp, Fr, Fp
  endelse
endif else begin
  if (Damp eq 0) then begin
    nvar=3		;vars are Br, Bp, alp_m
  endif else begin
    nvar=7		;vars are Br, Bp, Fr, Fp, Er, Ep, alp_m
  endelse
endelse
;
;  DEFINE ARRAYS FOR PHYSICAL VARIABLES
f=fltarr(nx,nvar)
;
;  INITIALIZE VARIABLES AND ARRAYS
first=0.	;for Runge-Kutta routine
t=0.		;time
tt=0.		;time array
Beq=fltarr(nx)
BBr=fltarr(nx)
BBp=fltarr(nx)
AM=fltarr(nx)
BBr[*]=0.		;Br array
BBp[*]=0.		;Bphi array
AM[*]=0.		;alpha_m array
EE=0.		;magnetic energy array
alp_m=0.	;Initial condition on alp_m
;
;  RADIUS ARRAY IN PHYSICAL UNITS
;r_kpc=y*r_kpc_disc
;
;  Beq PROFILE
;
;Beq=B0*exp(-x^2/2)
Beq=B0*exp((x/h)^2/2)
;
;  ALPHA PROFILE
;
if (Alp_sin eq 1) then begin
  alp_k=R_alpha*sin(n_alp*!pi*x/h)
  dalp_kdz=n_alp*!pi*R_alpha*cos(n_alp*!pi*x/h)
endif else begin
  alp_k=R_alpha*x/h
  dalp_kdz=R_alpha/h+x*0
endelse
;
;  ETAT PROFILE
;
etat= etat0*exp((x/h)^2)	;basically make v_turb a Guassian and rho, tau constant
;
;  V_TURB PROFILE
;
v_turb= v_turb0*exp((x/h)^2/2)
;
;  L PROFILE
;
l= l0*exp((x/h)^2/2)
;
etat= etat0*exp((x/h)^2)
;
;  VERTICAL MEAN VELOCITY PROFILE
;
Uz=U0*x/h		;Vertical mean velocity in normalized units h/td
dUzdz=U0/h+x*0	;Vertical mean velocity in normalized units h/td
;
;  SEED FIELD
;
default,seed,seed_id
f(*,0)=Bseed*gaunoise(nx,seed=seed)
f(*,1)=Bseed*gaunoise(nx,seed=seed)
Br_init=f[*,0]
f(0,0:1)=0.		;boundary condition 1 Br=Bp=0 at z=-1
f(nx-1,0:1)=0.		;boundary condition 2 Br=Bp=0 at z=+1
;
;  CALCULATED QUANTITIES
;
Dyn=R_alpha*R_omega
;
print,'UNITS'
print,'td0_Gyr=',td0_Gyr
print,'h_kpc=',h_kpc
print,'B0_Gauss=',B0_Gauss
print,'In addition, etat0_cm2s=',etat0_cm2s
print,''
print,'NUMERICS:'
print,'nvar',nvar
print,''
print,'MODULES:'
print,'Damp ',Damp
print,'alp^2',Alp_squared
print,''
print,'PARAMETERS:'
print,'tau= ',tau
print,'h=   ',h
print,'etat0=',etat0
print,'R_alpha=alp0*h/etat=',R_alpha
print,'R_omega=G*h^2/etat =',R_omega
print,'Dynamo number      =',Dyn
print,'C_cm2sGauss2       =',C_cm2sGauss2
print,'C_ad               =',C_ad
print,'ion_frac           =',ion_frac
print,'Directory ending   =',s0
!p.multi=[0,2,2]
!p.charsize=1.5
plot,x,alp_k,xtit='z',ytit='alp_k'
plot,x,dalp_kdz,xtit='z',ytit='dalp_kdz'
plot,x,Uz,xtit='z',ytit='Uz'
plot,x,dUzdz,xtit='z',ytit='dUzdz'
;
;  VIEWING WINDOWS
;
window,0,xsize=500,ysize=700
;window,3,xsize=1200,ysize=500
;
end
