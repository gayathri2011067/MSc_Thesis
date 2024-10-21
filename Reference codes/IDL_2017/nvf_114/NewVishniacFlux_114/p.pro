loadct,39
col1=250
col2=150
col3=90
col4=45
!p.multi=0
!p.charthick=1. & !p.thick=2. & !x.thick=2. & !y.thick=2.
!p.charsize=1.0
;
;  RUN NUMBER
s0='NewVishniacFlux_114'
;
;  DIRECTORY NAMES
program_direc=strjoin(['/Users/luke/pde/1D/telegraph/z/program/',s0])
data_direc=strjoin(['/Users/luke/pde/1D/telegraph/z/data/',s0])
plots_direc=strjoin(['/Users/luke/pde/1D/telegraph/z/plots/',s0])
;
cd,data_direc
;
openr,1,'numerical.dat'
openr,2,'physical.dat'
;
readf,1,nx,ncurve,neltstt
;
nx=fix(nx) & ncurve=fix(ncurve) & neltstt=fix(neltstt)
;
readf,2,B0,Beq
print,'nx,ncurve,neltstt,Beq',nx,ncurve,neltstt,Beq
;
for ifile=1,2 do begin
  close,ifile
endfor
;
BBrstore=fltarr(ncurve,nx,neltstt)
BBpstore=fltarr(ncurve,nx,neltstt)
AMstore= fltarr(ncurve,nx,neltstt)
EEstore= fltarr(ncurve,neltstt)
ttstore= fltarr(ncurve,neltstt)
z=fltarr(nx)
;
openr,1,'ttstore.dat'
openr,2,'BBrstore.dat'
openr,3,'BBpstore.dat'
openr,4,'AMstore.dat'
openr,5,'EEstore.dat'
openr,6,'z.dat'
readf,1,ttstore
readf,2,BBrstore
readf,3,BBpstore
readf,4,AMstore
readf,5,EEstore
readf,6,z
print,'size(AMstore)',size(AMstore)
print,'size(ttstore)',size(ttstore)
for ifile=1,6 do begin
  close,ifile
endfor
;
cd,plots_direc
;
set_plot,'ps'
;
device,filename='alpha_m.ps',/color,xsize=16,ysize=10
!p.charsize=1.
xr_alp=[-1,1]
yr_alp=100*[-1,1]
 plot,z,AMstore[0,*,1000],xr=xr_alp,yr=yr_alp,xtit='z/h',ytit='alpha_m';,tit='R_alpha=0, R_omega=-20, U0=0, with alpha^2 effect'
oplot,z,AMstore[0,*,1000],col=-1
oplot,z,AMstore[0,*,500] ,col=col3,li=0
oplot,z,AMstore[0,*,250] ,col=col4,li=0
h=1.
v_turb= 15*exp(z^2/h^2/2)
oplot,z,v_turb,li=2
alp_k=0*1.5*sin(!pi*z);1.5*sin(!pi*z)
oplot,z,z*0,li=1
al_legend,['alpha_m','rms turb veloc'],linestyle=[0,2],colors=[-1,-1],position=[0.2*xr_alp[1],0.6*yr_alp[0]],box=0,charsize=0.9
device,/close
;
device,filename='B.ps',/color,xsize=16,ysize=10
!p.charsize=0.8
xr_B=[-1,1]
yr_B=[-0.8,1.6]
print,'ttstore[0,1000]', ttstore[0,1000]
  plot,z,BBrstore[0,*,1000],xr=xr_B,yr=yr_B,xtit='z/h',ytit='Bi';,tit='R_alpha=0, R_omega=-20, U0=0, with alpha^2 effect, t/t_d=20.0-20.5'
 oplot,z,BBrstore[0,*,1000],col=-1  ,li=0
 oplot,z,BBpstore[0,*,1000],col=-1  ,li=2
 oplot,z,BBrstore[0,*,500] ,col=col3,li=0
 oplot,z,BBpstore[0,*,500] ,col=col3,li=2
 oplot,z,BBrstore[0,*,250] ,col=col4,li=0
 oplot,z,BBpstore[0,*,250] ,col=col4,li=2
 oplot,z,0*z,li=1
al_legend,['Br','Bphi'],linestyle=[0,2],colors=[-1,-1],position=[0.4*xr_B[1],0.9*yr_B[1]],box=0,charsize=0.9
device,/close
;
device,filename='p.ps',/color,xsize=16,ysize=10
!p.charsize=0.8
xr_p=[-1,1]
yr_p=[-90,90]
pdeg= 180./!pi*atan(BBrstore/BBpstore)
  plot,z,-pdeg[0,*,1000],xr=xr_p,yr=yr_p,xtit='z/h',ytit='-p (deg)';,tit='R_alpha=0, R_omega=-20, U0=0, with alpha^2 effect, t/t_d=20.0-20.5'
 oplot,z,-pdeg[0,*,1000],col=-1,li=0
 oplot,z,-pdeg[0,*,500] ,col=col3,li=0
 oplot,z,-pdeg[0,*,250] ,col=col4,li=0
device,/close
;
;
device,filename='time_evolution.ps',/color,xsize=16,ysize=10
!p.charsize=1.
xr=[0,10]
yr=[2d-3,2]
plot,ttstore[0,*],sqrt(EEstore[0,*])/B0,xr=xr,yr=yr,xtit='time/t!Dd!N',ytit='Log(<B>/B!D0!N)',/ylog;,tit='R_alpha=0, R_omega=-20, U0=0, with alpha^2 effect'
oplot,ttstore[0,*],sqrt(EEstore[0,*])/B0,col=col1
device,/close
cd,program_direc
;
end
