; $Id: pde.pro,v 1.3 2011-02-12 07:46:52 luke Exp $
;***********************************************************************
FUNCTION pde, f
COMMON cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0
@common
;
      dfdt = make_array(size=size(f), /nozero)   ; make a similar array
;
;  compute df/dt
;
;-----------------------------------------------------------------------
;
;  TIME-DEPENDENT ALPHA
;
;  USE EXPLICIT NAMES
;
Br    = f(*,0)
Bp    = f(*,1)
;
if (Damp eq 1) then begin
  Fr=f(*,2)
  Fp=f(*,3)
  if (Dyn_quench eq 1) then begin
    Er=f(*,4)
    Ep=f(*,5)
    alp_m=f(*,6)
  endif
endif else begin
  Fr=0*Br
  Fp=0*Bp
  if (Dyn_quench eq 1) then begin
    alp_m=f(*,2)
  endif
endelse
;
dBrdz=xder(Br)
d2Brdz2=xder2(Br)
dBpdz=xder(Bp)
d2Bpdz2=xder2(Bp)
;
;dBrdp=yder(Br)
;d2Brdp2=yder2(Br)
;dBpdp=yder(Bp)
;d2Bpdp2=yder2(Bp)
;
if (Dyn_quench eq 1) then begin
  dalp_mdz = xder(alp_m)
  d2alp_mdz2 = xder2(alp_m)
  alp      = alp_k +alp_m
  dalpdz   = dalp_kdz +dalp_mdz
endif else begin
  dalp_mdz = 0.
  alp      = alp_k
  dalpdz   = dalp_kdz
endelse
  
;  KINETIC ALPHA PROFILE
;
if (Alp_sin eq 1) then begin
  alp_k=R_alpha*sin(n_alp*!pi*x/h)
  dalp_kdz=n_alp*!pi*R_alpha*cos(n_alp*!pi*x/h)
endif else begin
  alp_k=R_alpha*x/h
  dalp_kdz=R_alpha/h+x*0
endelse
;
;DivVishniac = -fac*R_omega*(2.*l/3)^2*x
DivVishniac = +8./9*fac*R_omega*l^2*x
;
;  CALCULATE Emag
;
Emag=Br^2+Bp^2
;
;  LIST OF VARIABLE NAMES FOR f ARRAY
;
;  UNDER FOSA				UNDER TAU APPROXIMATION
;  f(*,*,0)=Br				f(*,*,0)=Br
;  f(*,*,1)=Bp				f(*,*,1)=Bp
;  f(*,*,2)=alp_m			f(*,*,2)=Fr
;  					f(*,*,3)=Fp
;  					f(*,*,4)=Er
;		  			f(*,*,5)=Ep
;					f(*,*,6)=alp_m
;
;  METHOD: ALL ARRAYS ARE 2-DIMENSIONAL
;
;r(0,*)=0.01
;mark
if (Damp eq 0) then begin
;
;  CASE 1: FOSA (tau_mta-->0 LIMIT)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
;
  dfdt(*,0) =             -dUzdz*Br -Uz*dBrdz -dalpdz*Bp -alp*dBpdz +etat*d2Brdz2 +Rm_inv*etat*d2Brdz2
  dfdt(*,1) = +R_omega*Br -dUzdz*Bp -Uz*dBpdz +dalpdz*Br +alp*dBrdz +etat*d2Bpdz2 +Rm_inv*etat*d2Bpdz2
  if (Alp_squared eq 0) then begin
    dfdt(*,1) = +R_omega*Br -dUzdz*Bp -Uz*dBpdz                     +etat*d2Bpdz2 -Rm_inv*d2Bpdz2
  endif
  if (Dyn_quench eq 1) then begin
    dfdt(*,2) = -2./3/tau*( alp*(Br^2 +Bp^2)/Beq^2 -etat*(Bp*dBrdz-Br*dBpdz)/Beq^2 +Rm_inv*alp_m +DivVishniac) $
                -dalp_mdz*Uz -alp_m*dUzdz +kappa*etat*d2alp_mdz2 +2*C_ad*xder(alp_m*xder(Br^2+Bp^2))
  endif
  ;
endif else begin
  ;
  ;  CASE 2: MTA (FINITE tau_mta)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
  ;
  dfdt(*,0) =             -dUzdz*Br -Uz*dBrdz +Fr +Rm_inv*etat*d2Brdz2
  dfdt(*,1) = +R_omega*Br -dUzdz*Bp -Uz*dBpdz +Fp +Rm_inv*etat*d2Bpdz2
  dfdt(*,2) = (-dalpdz*Bp -alp*dBpdz +etat*d2Brdz2 -Fr)/tau_mta
  dfdt(*,3) = (+dalpdz*Br +alp*dBrdz +etat*d2Bpdz2 -Fp)/tau_mta
  if (Alp_squared eq 0) then begin
    dfdt(*,3) = (                    +etat*d2Bpdz2 -Fp)/tau_mta
  endif
  if (Dyn_quench eq 1) then begin
    dfdt(*,4) = (+alp*Br+etat*dBpdz -Er)/tau_mta
    dfdt(*,5) = (+alp*Bp-etat*dBrdz -Ep)/tau_mta
    dfdt(*,6) = -2./3/tau*( (Er*Br +Ep*Bp)/Beq^2 +Rm_inv*alp_m +DivVishniac) $
                -dalp_mdz*Uz -alp_m*dUzdz +kappa*etat*d2alp_mdz2 +2*C_ad*xder(alp_m*xder(Br^2+Bp^2))
  endif
  ;
endelse
;
;  IMPOSE BOUNDARY CONDITIONS (IF INITIAL VALUE IS ZERO AND dfdt=0 THIS IS EQUIVALENT TO SETTING f=0 AT ALL TIMES)
;
;if (Neumann_inner eq 0) then begin
    ;dfdt(0,*,0:1)   =0.d	;Ensure that Br, Bp remain zero at r=0
    ;dfdt(nx-1,*,0:1)=0.d	;Ensure that Br, Bp remain zero at r=R
  ;if (Dyn_quench eq 1) then begin
    ;if (Damp eq 0) then begin
      ;dfdt(0,*,2)=0.d		;Put BC alp_m=0 at r=0 to avoid having it blow up at the origin (compression term)
    ;endif else begin
      ;dfdt(0,*,6)=0.d		;Put BC alp_m=0 at r=0 to avoid having it blow up at the origin (compression term)
    ;endelse
  ;endif
;endif else begin
  ;dfdt(0,*,0:1)   =0.d	;Ensure that Br, Bp remain zero at r=0



  ;dfdt(nx-1,*,0:1)=0.d	;Ensure that Br, Bp remain zero at r=R
  ;f(1,*,0:1)=f(0,*,0:1)		;Make derivative zero at inner boundary--wait!! Must modify eqn for Axel's 6th order scheme



  ;print,'ERROR: Neumann boundary conditions must be specified'
  ;stop
;endelse
;
return, dfdt
END
