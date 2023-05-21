Subroutine integ1D_pot(ele,n_el,i,j,result)
    !Integral para el elemento i,j <Snp|V|Sn>
    use sturmians_base
    use sturmians_data
    integer i,j,n_el,ele
    real*8 xi
    real*8,parameter :: tol=1d-16
    complex*16 result,temp,v

    result=(0d0,0d0)
    temp=(0d0,0d0)
    do k=1,npoints_in
	if(pot(n_el).eq.'yuk')then
            v=(-1d0)*exp(-alfa(n_el,ele)*bs(n_el)%x_in(k))/bs(n_el)%x_in(k)
	endif
	if(pot(n_el).eq.'exp')then
            v=(-1d0)*exp(-alfa(n_el,ele)*bs(n_el)%x_in(k))
	endif
        result=result+bs(n_el)%w_in(k)*bs(n_el)%s_in(k,i,0)*bs(n_el)%s_in(k,j,0)*v
    enddo

    return
end subroutine integ1D_pot


!------------------------------------------------------------------------------
Subroutine integ1D_ov(ele,n_el,i,j,result)
    !Integral para el elemento i,j <Snp|Sn>
    use sturmians_base
    use sturmians_data
    integer i,j,n_el,ele,k
    real*8,parameter :: tol=1d-16
    complex*16 result,temp

    result=(0d0,0d0)
    temp=(0d0,0d0)
    do k=1,npoints
        result=result+bs(n_el)%w(k)*bs(n_el)%s((k-1)*4,i,0)*bs(n_el)%s((k-1)*4,j,0)
    enddo

    return
end subroutine integ1D_ov


