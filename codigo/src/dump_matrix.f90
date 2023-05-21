subroutine dump_matrix(cfile,A,irowsize,icolsize,N,M)
!> \brief Dumps the local Matrix matrix
!!  Note that use of pzlawrite is not recommended, since it
!!  is unefficient. This piece of code simply writes
!!  nproc files, each with the local matrix.
!!  Postprocessing can be achieved with any script.
!!
    use mapping_module
    use debug_module
    use procdata
    use parablock
    implicit none
    complex*16 :: A(irowsize,icolsize) !< Local A Matrix 
    integer irowsize                   !< Local Row Size of A
    integer icolsize                   !< Local Column Size of A
    integer N                         !< Global size of matrix A, rows
    integer M                         !< Global size of matrix A, columns
    integer un
    integer i,j,ipr,ipc,irow,icol
    character*1 cfile
    character*3 itxt
  
    un=110+iam
    write(itxt,'(i3.3)') iam
    open(unit=un,file =cfile//itxt//'.dat')
    do i=1,N
        do j=1,M
!            write(110+iam,*) i,j
            call mapping(i,j,irow,icol,ipr,ipc)
            if (ipr.eq.myrow.and.ipc.eq.mycol) then
                write(un,'(2i6.5,2e18.8)') i,j,A(irow,icol)
            endif
        enddo
    enddo


    return

end subroutine dump_matrix

