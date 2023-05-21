subroutine scatt_func
 use sturmians_base
 use sturmians_data
 use size_module
 use parablock
 integer i1,i2,un1,un2
 integer npt
 complex*16 psi,psi_sc,der
 character*3 itxt

un1=110+iam
un2=120+iam
write(itxt,'(i3.3)') iam

open(unit=un1,file='Re_PsiSc_TP'//itxt//'.dat')
open(unit=un2,file='Im_PsiSc_TP'//itxt//'.dat')

 npt=int(nrmax(1)/nproc)

do i1=iam*npt,(iam+1)*npt-1
do i2=0,(8+1)*npt-1
 psi_sc=(0d0,0d0)
 psi=(0d0,0d0)
! call fscatt_nsym(bs(1)%r(i1),bs(1)%r(i2),psi_sc)
 call fscatt(bs(1)%r(i1),bs(2)%r(i2),psi_sc,psi,der)
 write(un1,'(f13.6,1x,f13.6,1x,f20.10)') bs(1)%r(i1),bs(1)%r(i2),real(psi_sc)
 write(un2,'(f13.6,1x,f13.6,1x,f20.10)') bs(1)%r(i1),bs(1)%r(i2),aimag(psi_sc)
enddo
write(un1,*)
write(un2,*)
enddo

 close(un1)
 close(un2)

return
end subroutine scatt_func
!------------------------------------------------------------------------------
subroutine fscatt_nsym(r1,r2,psi_sc)
use sturmians_base
use sturmians_data
use size_module
use problemdata_module
use pi_module
use array_module
implicit none
integer in,n1,n2,i1,i2
real*8 r1,r2
complex*16 psi_sc,s1,s2

 call findi(bs(1)%r,r1,nrmax(1),i1)
 call findi(bs(1)%r,r2,nrmax(1),i2)
 i1=i1-1
 i2=i2-1

 psi_sc=(0d0,0d0)
 do in=1,imax
 call indexado(in,0,0,n1,n2)
 s1=sa1(i1,n1,0)+r1*(sb1(i1,n1,0)+r1*(sc1(i1,n1,0)+r1*sd1(i1,n1,0)))
 s2=sa2(i2,n2,0)+r2*(sb2(i2,n2,0)+r2*(sc2(i2,n2,0)+r2*sd2(i2,n2,0)))
 psi_sc=psi_sc+bn(in)*(s1*s2)
 enddo

return
end subroutine fscatt_nsym
!------------------------------------------------------------------------------
subroutine fscatt(r1,r2,psi_sc,psi,der)
use sturmians_base
use sturmians_data
use size_module
use problemdata_module
use pi_module
use array_module
implicit none
integer in,n1,n2,i1,i2,i1p,i2p
real*8 r1,r2,dh,rho,alpha
complex*16 psi,psi_sc,s1,s2,s1_exc,s2_exc,der
complex*16 psi_sc_dh

 dh=1d-3

 call findi(bs(1)%r,r1,nrmax(1),i1)
 call findi(bs(1)%r,r2,nrmax(1),i2)
 i1=i1-1
 i2=i2-1

 rho=sqrt(r1*r1+r2*r2)
 if(r1.ne.0d0)then
 alpha=atan(r2/r1)
 else
 alpha=pi/2d0
 endif

 call findi(bs(1)%r,(rho+dh)*cos(alpha),nrmax(1),i1p)
 call findi(bs(2)%r,(rho+dh)*sin(alpha),nrmax(1),i2p)
 i1p=i1p-1
 i2p=i2p-1

 psi_sc=(0d0,0d0)
 psi_sc_dh=(0d0,0d0)
 do in=1,imax
 call indexado(in,0,0,n1,n2)
 s1=sa1(i1,n1,0)+r1*(sb1(i1,n1,0)+r1*(sc1(i1,n1,0)+r1*sd1(i1,n1,0)))
 s2=sa2(i2,n2,0)+r2*(sb2(i2,n2,0)+r2*(sc2(i2,n2,0)+r2*sd2(i2,n2,0)))
 s1_exc=sa1(i2,n1,0)+r2*(sb1(i2,n1,0)+r2*(sc1(i2,n1,0)+r2*sd1(i2,n1,0)))
 s2_exc=sa2(i1,n2,0)+r1*(sb2(i1,n2,0)+r1*(sc2(i1,n2,0)+r1*sd2(i1,n2,0)))
 psi_sc=psi_sc+bn(in)*(s1*s2+(-1d0)**spin*s1_exc*s2_exc)
 s1=sa1(i1p,n1,0)+(rho+dh)*cos(alpha)*(sb1(i1p,n1,0)+(rho+dh)*cos(alpha)*(sc1(i1p,n1,0) &
    +(rho+dh)*cos(alpha)*sd1(i1p,n1,0)))
 s2=sa2(i2p,n2,0)+(rho+dh)*sin(alpha)*(sb2(i2p,n2,0)+(rho+dh)*sin(alpha)*(sc2(i2p,n2,0) &
    +(rho+dh)*sin(alpha)*sd2(i2p,n2,0)))
 s1_exc=sa1(i2p,n1,0)+(rho+dh)*sin(alpha)*(sb1(i2p,n1,0)+(rho+dh)*sin(alpha)*(sc1(i2p,n1,0) &
    +(rho+dh)*sin(alpha)*sd1(i2p,n1,0)))
 s2_exc=sa2(i1p,n2,0)+(rho+dh)*cos(alpha)*(sb2(i1p,n2,0)+(rho+dh)*cos(alpha)*(sc2(i1p,n2,0) &
    +(rho+dh)*cos(alpha)*sd2(i1p,n2,0)))
 psi_sc_dh=psi_sc_dh+bn(in)*(s1*s2+(-1d0)**spin*s1_exc*s2_exc)
 enddo

 !Derivada de la parte de scattering respecto de rho
 der=(psi_sc_dh-psi_sc)/dh

 k0=sqrt(2d0*E0)
 psi=psi_sc+sqrt(4d0*pi)*(sin(k0*r1)*2d0*r2*exp(-r2)+(-1d0)**spin*sin(k0*r2)*2d0*r1*exp(-r1))/k0

return
end subroutine fscatt
!------------------------------------------------------------------------------------------------
subroutine ffinal(final,r1,r2,k1,k2,psi,der)
!der es la derivada del estado final respecto de rho
use sturmians_data
implicit none
 real*8 r1,r2,k1,k2
 real*8 f,fp,g,gp,err,dh,rho,alpha
 complex*16 psi,der,phik1,phik2,phik1_dh,phik2_dh
 character*2 final

 dh=1d-3
 rho=sqrt(r1*r1+r2*r2)
 alpha=atan(r2/r1)

 if(final.eq.'pw')then
 phik1=sin(k1*r1)
 phik2=sin(k2*r2)
 psi=phik1*phik2
 der=k2*sin(alpha)*sin(k1*r1)*cos(k2*r2)+k1*cos(alpha)*cos(k1*r1)*sin(k2*r2)
 endif

 if(final.eq.'dw')then
 call scoul(Z(1),k1*k1/2d0,0,r1,f,fp,g,gp,err)  
 phik1=cmplx(f )
 call scoul(Z(1),k2*k2/2d0,0,r2,f,fp,g,gp,err)  
 phik2=cmplx(f)
 psi=phik1*phik2
 call scoul(Z(1),k1*k1/2d0,0,(rho+dh)*cos(alpha),f,fp,g,gp,err)  
 phik1_dh=cmplx(f) 
 call scoul(Z(1),k2*k2/2d0,0,(rho+dh)*sin(alpha),f,fp,g,gp,err)  
 phik2_dh=cmplx(f)  
 der=(phik1_dh*phik2_dh-phik1*phik2)/dh
 endif

return
end subroutine ffinal
