subroutine matrix_element_tp(ijp_in,ij_in,matrix,overlap,iam)
    use size_module
    use sturmians_base
    use sturmians_data
    use problemdata_module
    implicit none
    integer, optional :: iam
    integer ijp_in,ij_in
    integer ijp,ij,j1,j2,j1p,j2p
    real*8 time1,time2
    complex*16 matrix,v1,v2,o1,o2,r12
    complex*16 overlap

    ij=ij_in
    ijp=ijp_in

    call indexado(ijp,0,0,j1p,j2p)
    call indexado(ij,0,0,j1,j2)

    call integ1D_pot(0,1,j1,j1p,v1) !<S(j1p)|V|S(j1)>
    call integ1D_ov(0,1,j1,j1p,o1)  !<S(j1p)|S(j1)>

    call integ1D_pot(0,2,j2,j2p,v2) !<S(j2p)|V|S(j2)>
    call integ1D_ov(0,2,j2,j2p,o2)  !<S(j2p)|S(j2)>

    call integ2D(ijp,ij,r12)      !<S(j1p)|<S(j2p)|1/r>|S(j1)>|S(j2)>

    matrix  = bs(1)%beta(j1,0)*v1*o2+bs(2)%beta(j2,0)*v2*o1-r12-(E(1)+E(2)-Etot)*o1*o2
    overlap =o1*o2
 
    return
end subroutine matrix_element_tp
!---------------------------------------------------------------------------------------------
subroutine vector_element_tp(ijp_in,vector,iam)
    use size_module
    use sturmians_base
    use sturmians_data
    use pi_module
    use problemdata_module

    implicit none
    integer, optional :: iam
    integer ijp_in,ijp,j1p,j2p
    complex*16 vector,c1,b2,r12

    ijp=ijp_in

    call indexado(ijp,0,0,j1p,j2p)

    call integ1D_f0_coulomb(0,1,j1p,c1)  !<S(j1p)|Z/r|sin(k0*r)>
    call integ1D_f0_bound(0,2,j2p,b2)  !<S(j2p)|phi0>
    call integ2D_f0(ijp,r12)       !<S(j1p)|<S(j2p)|1/r>|sin(k0*r1)>|phi0(r2)>

    vector=sqrt(4d0*pi)*(c1*b2+r12)/k0

    return
end subroutine vector_element_tp
