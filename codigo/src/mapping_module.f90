module mapping_module
contains
!
!***********************************************************************
!
subroutine mapping(I,J,irow,icol,ipr,ipc)
    use matrixdat
    use procdata
    implicit none
    integer I,J
    integer ipr,ipc
    integer mblock
    integer iypos
    integer irow,icol
    integer lblock
    integer ixpos

    !       mapping the GLOBAL matrix A(I,J) into LOCAL A(irow,icol)

    !       ...ipr and ipc: row and column processes
    ipr = mod(irsrc + (I-1)/MB , nprow)
    ipc = mod(icsrc + (J-1)/NB , npcol)

    !       column identification (J --> icol)
    mblock = (J-1)/(npcol*NB)
    iypos  = mod(J-1,NB) + 1
    icol   = mblock*NB + iypos
    !       row identification    (I --> irow)
    lblock = (I-1)/(nprow*MB)
    ixpos  = mod(I-1,MB) + 1
    irow   = lblock*MB + ixpos

    return
end subroutine mapping


!
!***********************************************************************
!
subroutine demapping(irow,icol,I,J)
    use matrixdat
    use procdata
    implicit none
    integer irow,icol
    integer I,J
    integer ix,iy
    integer l,m

    !       demapping the local (irow,icol) in the global (I,J) position

    ix = mod(irow,nb)
    if (ix.eq.0) ix=nb
    iy = mod(icol,mb)
    if (iy.eq.0) iy=mb
    !        l and m: local block coordinate
    l = (irow-ix)/nb
    m = (icol-iy)/mb
    I = (l*nprow + myrow) * nb + ix
    J = (m*npcol + mycol) * mb + iy

    return
end subroutine demapping
end module mapping_module
