;1LC: This version has been modified by commenting out openW and printF commands, for use with trial version of IDL.
;	For original version see r_old.pro.	2012/05/22
;
col1=250
col2=150
col3=50
col4=20
;
for iDyn=0,ncurve-1 do begin
;
R_omega=Rom(iDyn)			;Set R_omega to next value
;
;  RECORD REAL TIME WHEN SIMULATION STARTS
;
start_time=systime(1)
;
Bz=fltarr(nx)
Fz=fltarr(nx)
;
;default,w,.1
;n1=10 & n2=800	;REPLOT N1 TIMES (AFTER N2 TIMESTEPS)
;
for i=0,n1-1 do begin
  for j=0L,n2-1 do begin
    f(0,0)=0.		;IMPOSE BOUNDARY CONDITIONS FOR RK ROUTINE
    f(0,1)=0.		;IMPOSE BOUNDARY CONDITIONS FOR RK ROUTINE
    f(nx-1,0)=0.
    f(nx-1,1)=0.
    rk,f,t,dt
    f(0,0)=0.		;IMPOSE BOUNDARY CONDITIONS AGAIN FOR PLOTTING PURPOSES
    f(0,1)=0.		;IMPOSE BOUNDARY CONDITIONS AGAIN FOR PLOTTING PURPOSES
    f(nx-1,0)=0.
    f(nx-1,1)=0.
  endfor
  ;
  ;  PERTURBATION THEORY
  ;
  C0=sqrt(2.)		;*sqrt(2./(1.+4./!pi))
  Br_pert=R_alpha*C0*(cos(!pi*x/2)+3./4/!pi^(3./2)*(-Dyn)^(1./2)*cos(3.*!pi*x/2))
  Bp_pert=-2*C0*(-Dyn/!pi)^(1./2)*cos(!pi*x/2)
  ;
  ;  PERFORM INTEGRAL TO GET Bz, Ez
  ;
  ;for j=0L,nx-1 do begin
    ;sumB=0.
    ;sumF=0.
    ;for k=j,nx-1 do begin
      ;sumB=sumB+Br(k)
      ;sumF=sumF+Fr(k)
    ;endfor
    ;if (j eq 0) then begin
      ;sumB0=sumB
      ;sumF0=sumF
    ;end
    ;Bz(j)=len/r*(sumB-sumB0/2)		;-sumB0/2 IS TO ENSURE THAT Bz(0)=0 (BOUNDARY CONDITION)
    ;Fz(j)=len/r*(sumF-sumF0/2)		;-sumB0/2 IS TO ENSURE THAT Fz(0)=0 (BOUNDARY CONDITION)
  ;endfor
;  Bz(0)=0. & Bz(nx-1)=0.
;  Fz(0)=0. & Fz(nx-1)=0.
  ;
  ;
  maxB=max([abs(max(Br))>abs(max(Bp)),abs(min(Br))>abs(min(Bp))])
  maxB_pert=max([abs(max(Br_pert))>abs(max(Bp_pert)),abs(min(Br_pert))>abs(min(Bp_pert))])
  ;
  Emag=(-Br*maxB_pert/maxB)^2+(-Bp*maxB_pert/maxB)^2
  Emag_pert=Br_pert^2+Bp_pert^2
  E=mean(Br^2+Bp^2)
  tt=[tt,t]
  BBr=[[BBr],[Br]]
  BBp=[[BBp],[Bp]]
  AM=[[AM],[alp_m]]
  EE=[EE,E]
  print,s0,'      time=',t,'       <B>/B0=',sqrt(E)/B0
  simulation_time=systime(1)-start_time			;DETERMINE REAL TIME TAKEN TO RUN SIMULATION IN SECONDS
  print,'simulation_time in minutes=',simulation_time/60
  if (E ne E) then begin
    print,'INFINITY!! CHANGE TIMESTEP.'
    stop
  endif
  ;
  nvert=3
  !p.multi=[0,1,nvert]
  !p.charthick=0.5*nvert & !p.thick=1 & !x.thick=1 & !y.thick=1
  !p.charsize=0.75*nvert
  ;
  xr_B=[-1.,1.]
  yr_B=max(abs(Bp/B0)>abs(Br/B0))*[-1,1]
  ;yr=max(abs(Bp*maxB_pert/maxB))*[-1.,0.5]
  ;mark1
  ;1)
   plot,x,Bp/B0,xr=xr_B,yr=yr_B,tit='B_phi/B_eq(____), B_r/B_eq(_ _ _)'
  oplot,x,Bp/B0,li=0,col=col2,th=2
  oplot,x,Br/B0,li=2,col=col1,th=2
  oplot,x,0*x,li=1
  ;oplot,x,Bp_pert/B0,li=3,col=col4,th=2
  ;oplot,x,Br_pert/B0,li=1,col=col3,th=2
  XYoutS,0.4,3.,'D='
  XYoutS,0.4,3.,Dyn
;  plot,x,Bz,xr=[0., 1.],yr=1.3*[min(Bz),max(Bz)],xtit='z',ytit='B_z'
;  oplot,x,zero_arr,li=1				;PLOT HORIZONTAL LINE FOR REFERENCE (zero_arr DEFINED IN start.pro)
;  oplot,[0,0],1.3*[min(Bz),max(Bz)],li=1		;PLOT VERTICAL LINE FOR REFERENCE
  ;
  ;2)
  xr_alp=[-1,1]
  yr_alp=max(abs(alp_m)>abs(alp_k)>abs(alp_k+alp_m))*[-1,1]
  plot,x,alp_k+alp_m,xr=xr_alp,yr=yr_alp,li=0,th=2,xtit='z/h',tit='alpha(___), alpha_k(_ _), alpha_m(_._)'
  oplot,x,alp_k,col=col2,li=2,th=2
  oplot,x,alp_m,col=col1,li=3,th=2
  
  ;3)
  ;xr_ev=[0,20]
  ;yr_ev=[-15.,1.]
  ;plot,tt,alog10(sqrt(EE)/B0),xr=xr_ev,yr=yr_ev,xtit='time/t!Dd!N',tit='Log(<B>/B!Deq!N)'
   plot,tt,alog10(sqrt(EE)/B0),xtit='time/t!Dd!N',tit='Log(<B>/B!Deq!N)'
  oplot,tt,alog10(sqrt(EE)/B0),col=col1
  wait,0.1
;
endfor
;
if (iDyn eq 0) then begin
  sizett=size(tt)
  neltstt=sizett[3]
  BBrstore=fltarr(ncurve,nx,neltstt)
  BBpstore=fltarr(ncurve,nx,neltstt)
  AMstore= fltarr(ncurve,nx,neltstt)
  EEstore= fltarr(ncurve,neltstt)
  ttstore= fltarr(ncurve,neltstt)
endif
BBrstore[iDyn,*,*]=BBr
BBpstore[iDyn,*,*]=BBp
AMstore[iDyn,*,*]=AM
EEstore[iDyn,*]=EE
ttstore[iDyn,*]=tt
;
;  RESET
;
t=0.
f(*,*)=0.				;Reset simulation
f(*,0)=Br_init				;Use same seed field
tt=0.
BBr=fltarr(nx)
BBp=fltarr(nx)
AM=fltarr(nx)
BBr[*]=0.		;Br array
BBp[*]=0.		;Bphi array
AM[*]=0.		;alpha_m array
EE=0.
first=0
;
endfor
;
cd,data_direc
;
openw,1,'numerical.dat'
openw,2,'physical.dat'
printf,1,nx,ncurve,neltstt
printf,2,B0,Beq
for ifile=1,2 do begin
  close,ifile
endfor
;
openw,1,'ttstore.dat'
openw,2,'BBrstore.dat'
openw,3,'BBpstore.dat'
openw,4,'AMstore.dat'
openw,5,'EEstore.dat'
openw,6,'z.dat'
printf,1,ttstore
printf,2,BBrstore
printf,3,BBpstore
printf,4,AMstore
printf,5,EEstore
printf,6,x
for ifile=1,6 do begin
  close,ifile
endfor
;
cd,program_direc
;
set_plot,'x'
;
;
END
