Subroutine integ2D(ijp,ij,result)
    !Está armada para el caso de equal energy!!!!!
    !Da el elemento ijp(j1p,j2p), ij(j1,j2) de la integral 1/r>
    use sturmians_base
    use sturmians_data
    use size_module
    implicit none
    integer j1,j2,j1p,j2p,ij,ijp,i,kout
    real*8,parameter :: tol=1d-16
    real*8 time1,time2,h0
    complex*16,allocatable :: integ_int_r1(:),integ_int_r2(:)
    complex*16 result,integ_out_r1,integ_out_r2,temp

    integ_out_r1=(0d0,0d0)
    integ_out_r2=(0d0,0d0)
    result=(0d0,0d0)
    temp=(0d0,0d0)

    allocate(integ_int_r1(npoints),integ_int_r2(npoints))
    integ_int_r1=(0d0,0d0)
    integ_int_r2=(0d0,0d0)

    !Recupera los índices
    call indexado(ij,0,0,j1,j2)
    call indexado(ijp,0,0,j1p,j2p)

    !Integral afuera
    do kout=npoints-1,1,-1

	!Integral adentro. Resuelve con Simpson 3/8 
	h0=(bs(1)%x(kout+1)-bs(1)%x(kout))/4d0
	i=int(4d0*(kout-1d0))
	integ_int_r1(kout)=2d0*h0*(7d0*bs(1)%s(i,j1,0)*bs(1)%s(i,j1p,0)/bs(1)%x(kout) &
             +32d0*bs(1)%s(i+1,j1,0)*bs(1)%s(i+1,j1p,0)/(bs(1)%x(kout)+h0) &
             +12d0*bs(1)%s(i+2,j1,0)*bs(1)%s(i+2,j1p,0)/(bs(1)%x(kout)+2d0*h0) &
             +32d0*bs(1)%s(i+3,j1,0)*bs(1)%s(i+3,j1p,0)/(bs(1)%x(kout)+3d0*h0) &
             +7d0*bs(1)%s(i+4,j1,0)*bs(1)%s(i+4,j1p,0)/(bs(1)%x(kout)+4d0*h0))/45d0
	integ_int_r2(kout)=2d0*h0*(7d0*bs(1)%s(i,j2,0)*bs(1)%s(i,j2p,0)/bs(1)%x(kout) &
             +32d0*bs(1)%s(i+1,j2,0)*bs(1)%s(i+1,j2p,0)/(bs(1)%x(kout)+h0) &
             +12d0*bs(1)%s(i+2,j2,0)*bs(1)%s(i+2,j2p,0)/(bs(1)%x(kout)+2d0*h0) &
             +32d0*bs(1)%s(i+3,j2,0)*bs(1)%s(i+3,j2p,0)/(bs(1)%x(kout)+3d0*h0) &
             +7d0*bs(1)%s(i+4,j2,0)*bs(1)%s(i+4,j2p,0)/(bs(1)%x(kout)+4d0*h0))/45d0

	integ_int_r1(kout)=integ_int_r1(kout)+integ_int_r1(kout+1)
 	integ_int_r2(kout)=integ_int_r2(kout)+integ_int_r2(kout+1)


  !Calcula la integra afuera
        integ_out_r1= integ_out_r1+bs(1)%w(kout)*bs(1)%s((kout-1)*4,j1,0)*bs(1)%s((kout-1)*4,j1p,0)*integ_int_r2(kout)
        integ_out_r2= integ_out_r2+bs(1)%w(kout)*bs(2)%s((kout-1)*4,j2,0)*bs(2)%s((kout-1)*4,j2p,0)*integ_int_r1(kout)

    enddo !kout

    deallocate(integ_int_r1,integ_int_r2)

    result=integ_out_r1+integ_out_r2

    return
end subroutine integ2D
