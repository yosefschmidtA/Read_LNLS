pro holoxwin, option,channeltron

;This program plots holograms resulting from unpolarized light experiments
;by adding the calculated intensities for the two components. This program
;is only useful if the hologram has 4 and 3 fold symetry and the hologram of each
;comnponent was calculated with theta and phi between 0 and 90 degrees for
;4-fold and between 0 and 120 degrees for p3. 
;The program automatically generates the whole hologram from symmetry 
;considerations.
;
;arqin1 = first hologram file 
;arqin2 = second hologram file
;option = define which kind of projection to be use in the plot
;       = 'stereo' -> x=2*tan(theta/2.)*cos(phi)
;       = 'theta'  -> x=(theta)*cos(phi)
;       = 'k'      -> x=k*sin(theta)*cos(phi)

;Device, pseudo_color=16
;Device, Get_Visual_Depth=thisDepth
;IF thisDepth GT 8 THEN Device, Decomposed=0


Device, pseudo_color=1
Device, Get_Visual_Depth=thisDepth
IF thisDepth GT 8 THEN Device, Decomposed=0

;rfactortoplot_lnls,'pesaida1_Ba4d_BaOterm_in_plane_down_\[0-10\].out',0.0,0.0,thetatotal,phitotal,intensity,rafactor
read_LNLS_v3,'AB18_-',18,69,3,3,120,3,thetatotal,phitotal,intensity,channeltron
;read_lnls_v3,'MA16_-',18,69,3,0,120,3,thetatotal,phitotal,intensity,channeltron
;read_LNLS_v3,'MA11_',33,69,3,0,120,3,thetatotal,phitotal,intensity,channeltron


max_theta=max(thetatotal,min=min_theta)
print,'theta minimo=',min_theta
print,'theta maximo=',max_theta

pi2=!pi/2. & density=600

if (option eq 'stereo') then begin
 x=2*tan(thetatotal/2.)*cos(phitotal) & y=2*tan(thetatotal/2.)*sin(phitotal)
endif

if (option eq 'theta') then begin
 ;correction for 'const fild'
 ;phitotal=phitotal+thetatotal*sin(17*!pi/180)
 ;
 x=(thetatotal)*cos(phitotal) & y=(thetatotal)*sin(phitotal)
endif

if (option eq 'k') then begin
 k=10.80
 x=k*sin(thetatotal)*cos(phitotal) & y=k*sin(thetatotal)*sin(phitotal)
endif

triangulate,x,y,triang
result=trigrid(x,y,intensity,triang,nx=density,ny=density)
final=smooth(result,10)


a=60./(max(final)-min(final)) & b=(max(final)-min(final))/59. - min(final)

ll=findgen(60)/a - b

xtt=findgen(density)/(density/(max(x)-min(x))) + min(x)
ytt=findgen(density)/(density/(max(y)-min(y))) + min(y)
window, 0, xsize=800,ysize=760,retain=2
contour, final, xtt,ytt, $
         levels=ll,$
         xrange=[-1.8,1.9], xst=1,$
         yrange=[-1.9,1.8], yst=1,$
         /fill
phitt=findgen(101)/100.*2.*!pi
xfac=(4.*!pi/180.0)+max_theta


;plotando a parte externa do holograma
for i=0,1200 do begin
 k=float(i)
 plots,(k/100.+max_theta+0.5*!pi/180)*cos(phitt),(k/100.+max_theta+0.5*!pi/180)*sin(phitt),thick=4
endfor
print,min_theta
for i=0,100 do begin
 k=float(i)
 plots,(min_theta-2*!pi/180-k/100.*min_theta)*cos(phitt),(min_theta-2*!pi/180-k/100.*min_theta)*sin(phitt),thick=4
endfor
plots,max_theta*cos(phitt),max_theta*sin(phitt)
loadct,0
plots,xfac*cos(phitt),xfac*sin(phitt),color=1,thick=3
yfac=1.37
;yfac=1.5
plots,[-xfac,-xfac],[-yfac,0],color=1,thick=3
plots,[xfac,xfac],[-yfac,0],color=1,thick=3
plots,[-xfac,xfac],[-yfac,-yfac],color=1,thick=3

print,max_theta
step=max_theta/3.
for i=0,3 do begin
 x=i*step
 print,x
 plots,[x,x],[-yfac,-yfac-0.05],color=1,thick=3
 plots,[-x,-x],[-yfac,-yfac-0.05],color=1,thick=3
 strx=strcompress(fix(x*180.0/!pi),/remove_all)
 if i eq 0 then begin
  xyouts,x-0.04,-yfac-0.2,strx,color=1,chars=3,CHARTHICK=2
 endif else begin
  xyouts,x-0.09,-yfac-0.2,strx,color=1,chars=3,CHARTHICK=2
  xyouts,-x-0.09,-yfac-0.2,strx,color=1,chars=3,CHARTHICK=2
 endelse
endfor
xyouts,-0.06,-yfac-0.47,'!4h',color=1,chars=5, CHARTHICK=3
plots,[xfac,xfac-0.1],[0,0],color=1,thick=3
plots,[0,0],[xfac,xfac-0.1],color=1,thick=3
xyouts,xfac+0.08,-0.05,'0',color=1,chars=3, CHARTHICK=2
xyouts,-0.09,xfac+0.06,'90',color=1,chars=3, CHARTHICK=2
xyouts,xfac*cos(45.*!pi/180.)+0.1,xfac*sin(45.*!pi/180),'!4u',color=1,chars=5, CHARTHICK=3


;printing anysotropy
Print,'max anisotropy=',max(intensity)
Print,'min anisotropy=',min(intensity)
anisotropy=' '
anisotropy=strcompress(string(round(100*(max(intensity)+abs(min(intensity))))),/remove_all)+'%'
xyouts,-yfac+0.42,-xfac+0.035,anisotropy,color=1,chars=3,CHARTHICK=2
RA='Ra='+strmid(strcompress(string(rafactor),/remove_all),0,5)
xyouts, yfac-0.8,-xfac+0,RA,color=1,chars=3,CHARTHICK=2
end

