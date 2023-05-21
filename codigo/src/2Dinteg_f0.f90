Subroutine integ2D_f0(ijp,result)
    use sturmians_base
    use sturmians_data
    use size_module
    use problemdata_module
    implicit none
    integer j1p,j2p,ijp,i,kout
    real*8 xi,h0
    complex*16 result,integ_out_r1,integ_out_r2
    complex*16,allocatable :: integ_int_r1(:),integ_int_r2(:)

    result=(0d0,0d0)
    integ_out_r1=(0d0,0d0)
    integ_out_r2=(0d0,0d0)

    allocate(integ_int_r1(npoints),integ_int_r2(npoints))
    integ_int_r1=(0d0,0d0)
    integ_int_r2=(0d0,0d0)

    !Recupera los índices
    call indexado(ijp,0,0,j1p,j2p)

    !Integral afuera
    do kout=npoints-1,1,-1

	!Integral adentro. Resuelve con Simpson 3/8 
	h0=(bs(1)%x(kout+1)-bs(1)%x(kout))/4d0
	i=int(4d0*(kout-1d0))
	xi=bs(1)%x(kout)
	integ_int_r1(kout)=2d0*h0*(7d0*bs(1)%s(i,j1p,0)*sin(k0*xi)/xi &
             +32d0*bs(1)%s(i+1,j1p,0)*sin(k0*(xi+h0))/(xi+h0) &
             +12d0*bs(1)%s(i+2,j1p,0)*sin(k0*(xi+2d0*h0))/(xi+2d0*h0) &
             +32d0*bs(1)%s(i+3,j1p,0)*sin(k0*(xi+3d0*h0))/(xi+3d0*h0) &
             +7d0*bs(1)%s(i+4,j1p,0)*sin(k0*(xi+4d0*h0))/(xi+4d0*h0))/45d0
	integ_int_r2(kout)=2d0*h0*(7d0*bs(1)%s(i,j2p,0)*2d0*xi*exp(-xi)/xi &
             +32d0*bs(1)%s(i+1,j2p,0)*2d0*(xi+h0)*exp(-(xi+h0))/(xi+h0) &
             +12d0*bs(1)%s(i+2,j2p,0)*2d0*(xi+2d0*h0)*exp(-(xi+2d0*h0))/(xi+2d0*h0) &
             +32d0*bs(1)%s(i+3,j2p,0)*2d0*(xi+3d0*h0)*exp(-(xi+3d0*h0))/(xi+3d0*h0) &
             +7d0*bs(1)%s(i+4,j2p,0)*2d0*(xi+4d0*h0)*exp(-(xi+4d0*h0))/(xi+4d0*h0))/45d0

	integ_int_r1(kout)=integ_int_r1(kout)+integ_int_r1(kout+1)
	integ_int_r2(kout)=integ_int_r2(kout)+integ_int_r2(kout+1)

 !Calcula la integra afuera
        integ_out_r1=integ_out_r1+bs(1)%w(kout)*bs(1)%s((kout-1)*4,j1p,0)*sin(k0*xi)*integ_int_r2(kout)
        integ_out_r2=integ_out_r2+bs(1)%w(kout)*bs(2)%s((kout-1)*4,j2p,0)*2d0*xi*exp(-xi)*integ_int_r1(kout)

    enddo !kout

    deallocate(integ_int_r1,integ_int_r2)

    result=integ_out_r1+integ_out_r2

    return
end subroutine integ2D_f0
