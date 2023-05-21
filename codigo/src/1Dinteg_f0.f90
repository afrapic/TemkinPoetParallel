Subroutine integ1D_f0_coulomb(ele,n_el,j,result)
    !Interal para el elemento j
    use sturmians_base
    use sturmians_data
    use problemdata_module
    implicit none
    integer n_el,ele,j,k
    complex*16 result,f
    !<Snp|Z/r|sin(k0*r)>

    result=(0d0,0d0)
    do k=1,npoints
            f=bs(1)%s((k-1)*4,j,0)*Z(1)*sin(k0*bs(1)%x(k))/bs(1)%x(k)
        result=result+bs(1)%w(k)*f
    enddo

    return
end subroutine integ1D_f0_coulomb

Subroutine integ1D_f0_bound(ele,n_el,j,result)
    !Interal para el elemento j
    use sturmians_base
    use sturmians_data
    use problemdata_module
    implicit none
    integer n_el,ele,j,k
    complex*16 result,f
    !<Snp|phi0>

    result=(0d0,0d0)
    do k=1,npoints
            f=2d0*bs(1)%x(k)*exp(-bs(1)%x(k))*bs(2)%s((k-1)*4,j,0)
        result=result+bs(1)%w(k)*f
    enddo

    return
end subroutine integ1D_f0_bound

Subroutine integ1D_f0_exp(ele,n_el,j,result)
    !Interal para el elemento j
    use sturmians_base
    use sturmians_data
    use problemdata_module
    implicit none
    integer n_el,ele,j,k
    complex*16 result,f
    !<Snp|exp(-r)>

    result=(0d0,0d0)
    do k=1,npoints
            f=bs(1)%s((k-1)*4,j,0)*exp(-bs(1)%x(k))
        result=result+bs(1)%w(k)*f
    enddo

    return
end subroutine integ1D_f0_exp
