subroutine scatt_flux
use sturmians_base
use sturmians_data
use size_module
use problemdata_module
use array_module
 use parablock
implicit none
integer i,j,npt
complex*16 work(2)
 character*3 itxt

 write(*,*) 'iam:',iam

 write(itxt,'(i3.3)') iam
 open(unit=1000,file='Flux_TP'//itxt//'.dat')

npt=int(40d0/nproc)

 do i=iam*npt,(iam+1)*npt-1
 do j=0,(8+1)*npt-1
 	call sflux((i+1)*2d0,(j+1)*2d0,work)
 write(1000,'(f13.6,1x,f13.6,1x,es15.8,1x,es15.8)') (i+1)*2d0,(j+1)*2d0,2d0*real(work(1)),2d0*real(work(2))
 enddo
 enddo
 close(1000)

return
end subroutine scatt_flux
!-----------------------------------------------------------------------------------------------
subroutine sflux(r1,r2,fj)
use sturmians_base
use sturmians_data
use size_module
use problemdata_module
use array_module
use complex_module
implicit none
integer in,n1,n2,i1,i2,i1p,i2p
real*8 r1,r2,dh
complex*16 psi_sc,s1,s2,s1_exc,s2_exc,der
complex*16 psi_sc_dr1,psi_sc_dr2,fj(2)

 dh=1d-3

 call findi(bs(1)%r,r1,nrmax(1),i1)
 call findi(bs(2)%r,r2,nrmax(2),i2)
 i1=i1-1
 i2=i2-1

 call findi(bs(1)%r,r1+dh,nrmax(1),i1p)
 call findi(bs(2)%r,r2+dh,nrmax(2),i2p)
 i1p=i1p-1
 i2p=i2p-1

 psi_sc=(0d0,0d0)
 psi_sc_dr1=(0d0,0d0)
 psi_sc_dr2=(0d0,0d0)
 do in=1,imax
 call indexado(in,0,0,n1,n2)
 s1=sa1(i1,n1,0)+r1*(sb1(i1,n1,0)+r1*(sc1(i1,n1,0)+r1*sd1(i1,n1,0)))
 s2=sa2(i2,n2,0)+r2*(sb2(i2,n2,0)+r2*(sc2(i2,n2,0)+r2*sd2(i2,n2,0)))
 s1_exc=sa1(i2,n1,0)+r2*(sb1(i2,n1,0)+r2*(sc1(i2,n1,0)+r2*sd1(i2,n1,0)))
 s2_exc=sa2(i1,n2,0)+r1*(sb2(i1,n2,0)+r1*(sc2(i1,n2,0)+r1*sd2(i1,n2,0)))
 psi_sc=psi_sc+bn(in)*(s1*s2+(-1d0)**spin*s1_exc*s2_exc)
 s1=sa1(i1p,n1,0)+(r1+dh)*(sb1(i1p,n1,0)+(r1+dh)*(sc1(i1p,n1,0)+(r1+dh)*sd1(i1p,n1,0)))
 s2=sa2(i2,n2,0)+r2*(sb2(i2,n2,0)+r2*(sc2(i2,n2,0)+r2*sd2(i2,n2,0)))
 s1_exc=sa1(i2,n1,0)+r2*(sb1(i2,n1,0)+r2*(sc1(i2,n1,0)+r2*sd1(i2,n1,0)))
 s2_exc=sa2(i1p,n2,0)+(r1+dh)*(sb2(i1p,n2,0)+(r1+dh)*(sc2(i1p,n2,0)+(r1+dh)*sd2(i1p,n2,0)))
 psi_sc_dr1=psi_sc_dr1+bn(in)*(s1*s2+(-1d0)**spin*s1_exc*s2_exc)
 s1=sa1(i1,n1,0)+r1*(sb1(i1,n1,0)+r1*(sc1(i1,n1,0)+r1*sd1(i1,n1,0)))
 s2=sa2(i2p,n2,0)+(r2+dh)*(sb2(i2p,n2,0)+(r2+dh)*(sc2(i2p,n2,0)+(r2+dh)*sd2(i2p,n2,0)))
 s1_exc=sa1(i2p,n1,0)+(r2+dh)*(sb1(i2p,n1,0)+(r2+dh)*(sc1(i2p,n1,0)+(r2+dh)*sd1(i2p,n1,0)))
 s2_exc=sa2(i1,n2,0)+r1*(sb2(i1,n2,0)+r1*(sc2(i1,n2,0)+r1*sd2(i1,n2,0)))
 psi_sc_dr2=psi_sc_dr2+bn(in)*(s1*s2+(-1d0)**spin*s1_exc*s2_exc)
 enddo

 psi_sc_dr1=(psi_sc_dr1-psi_sc)/dh
 psi_sc_dr2=(psi_sc_dr2-psi_sc)/dh

 fj(1)=(conjg(psi_sc)*psi_sc_dr1-psi_sc*conjg(psi_sc_dr1))/2d0/ci
 fj(2)=(conjg(psi_sc)*psi_sc_dr2-psi_sc*conjg(psi_sc_dr2))/2d0/ci

return
end subroutine sflux
