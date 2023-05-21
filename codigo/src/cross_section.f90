subroutine cross_section(method,final,rf,sigma)
    !method: #1: Integral de volumen
    !        #2: Integral de superficie
    !final: 'pw': onda plana
    !       'dw': onda distorsionada
    !Si la integral es de volumen integra en el cuadrado de (r1,r2) de 0 a rf
    !Si la integral es de superficie integra en coordenada hiperesferica alfa=arctan(r2/r1) 
    !para rho=rf fijo.
    use problemdata_module
    use pi_module
    use parablock
    implicit none
    integer method,i,half,npt
    real*8 rf,interv,k1,k2,sigma(0:ninterv)
    character*2 final
    complex*16 integ

    sigma=0d0
    interv=Etot/(ninterv*1d0)
    half=int(1d0*ninterv/2d0)

    do i=iam,half,nproc
        if(i.eq.0)then
            k1=0.05d0
        else
            k1=sqrt(2d0*(1d0*i*interv))
        endif
        k2=sqrt(2d0*Etot-k1*k1)
        call amplitud(method,final,rf,k1,k2,integ)
        sigma(i)=pi*abs(integ)**2d0/k1/k2/k0/k0/8d0
        sigma(ninterv-i)=sigma(i)
    enddo

 	if(spin.eq.1)then
 	sigma=sigma*pi
 	endif

return
end subroutine cross_section
!------------------------------------------------------------------------------------------
subroutine amplitud(method,final,rf,k1,k2,result)
use array_module
use sturmians_data
use size_module
use problemdata_module
use gauleg_module
implicit none
!method: #1: Integral de volumen
!        #2: Integral de superficie
!final: 'pw': onda plana
!       'dw': onda distorsionada
!Si la integral es de volumen integra en el cuadrado de (r1,r2) de 0 a rf
!Si la integral es de superficie integra en coordenada hiperesferica alfa=arctan(r2/r1) 
!para rho=rf fijo.
 integer method,npts,i,j
 real*8 k1,k2,rf
 real*8,allocatable :: x(:),w(:)
 complex*16 result,psi_sc,psi,der,fstate,fder,v
 character*2 final

result=(0d0,0d0)

if(method.eq.1)then
 npts=int(rf/0.1d0)
 allocate(x(npts),w(npts))
 call gauleg(0d0,rf,x,w)
 !Loop para la integral en r1
 do i=1,npts
	!Loop para la integral en r2 
	do j=1,npts
	if(final.eq.'pw')then
	v=cmplx(Z(1)/x(i)+Z(1)/x(j)+1d0/max(x(i),x(j))) 
	endif
	if(final.eq.'dw')then	
	v=cmplx(1d0/max(x(i),x(j)))
	endif
	call fscatt(x(i),x(j),psi_sc,psi,der) 
	call ffinal(final,x(i),x(j),k1,k2,fstate,der)
	result=result+w(i)*w(j)*fstate*v*psi
	enddo
 enddo
 deallocate(x,w) 
endif

if(method.eq.2)then
 npts=int(rf/0.2d0)
 allocate(x(npts),w(npts))
 call gauleg(0d0,pi/2d0,x,w)
 do i=1,npts
 call fscatt(rf*cos(x(i)),rf*sin(x(i)),psi_sc,psi,der)
 call ffinal(final,rf*cos(x(i)),rf*sin(x(i)),k1,k2,fstate,fder)
 result=result+w(i)*(fstate*der-psi_sc*fder)
 enddo
 result=result*rf/2d0
 deallocate(x,w)
endif

return
end subroutine amplitud

