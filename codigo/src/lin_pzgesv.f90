subroutine printlin(a,b,x,irowsize,icolsize,n,nrhs)	
    use debug_module
    use parablock
    use procdata
    use matrixdat
    implicit none
    integer irowsize,icolsize         ! local  matrix size
    integer n,nrhs                    ! global matrix size
    complex*16 a(irowsize,icolsize)
    complex*16 af(irowsize,icolsize)
    complex*16 b(irowsize,nrhs)
    complex*16 x(irowsize,nrhs)
    complex*16 matrix

    logical, parameter :: check=.false.
    integer ierr
    character equed
    integer lda,info
    integer lwork, lrwork
    !     Local Arrays Descriptors
    integer, parameter :: DLEN_ = 9
    integer  DESCA( DLEN_ ), DESCB( DLEN_ ), DESCX( DLEN_ ), DESCAF( DLEN_)
    integer  DESCC(DLEN_),descb2(dlen_)
    real*8 ANORM,Xnorm,Bnorm,eps,resid
    integer locr,locc
    integer nq
    !      External Functions 
    integer  NUMROC
    real*8   PZLANGE,pdlamch

    include 'mpif.h'

    character*(MPI_MAX_PROCESSOR_NAME) name
    integer irlen

!    write(6,*) iam,irowsize,icolsize,nrhs,N

    call MPI_GET_PROCESSOR_NAME(name,irlen,ierr)
    nq = numroc(n,nb,mycol,icsrc,npcol)
    write(nout,'(a,i3,a,a)') ' In lin Compute Basis proc=:',iam,'in ',name
    write(nout,'(a,8i6)') 'In lin : ',iam,irowsize,icolsize,n,nrhs,ictxt,myrow,mycol
    call flush(nout)
    return
 end subroutine printlin





    !
    !   Wrapper for PZGESV scalapack
    !   Solves x for ax=b
    !
 subroutine lin_pzgesv(a,b,x,irowsize,icolsize,n,nrhs)	
     use debug_module
     use parablock
     use matrixdat
     use procdata
     implicit none
     integer irowsize,icolsize         ! local  matrix size
     integer n,nrhs                    ! global matrix size
     complex*16 a(irowsize,icolsize)
     complex*16 b(irowsize,nrhs)
     complex*16 x(irowsize,nrhs)
     complex*16 af(irowsize,icolsize)

     complex*16 matrix

     logical, parameter :: check=.false.
     integer ierr
     character equed
     integer, allocatable :: ipiv(:)
     real*8,  allocatable :: row(:),col(:)
     real*8,  allocatable :: rwork(:)
     real*8,  allocatable :: berr(:),ferr(:)
     complex*16,allocatable :: work(:)

     real*8 rcond
     real*8 swork(1)
     complex*16 cwork(1)
     integer lda,info
     integer lwork, lrwork
     !     Local Arrays Descriptors
     integer, parameter :: DLEN_ = 9
     integer  DESCA( DLEN_ ), DESCB( DLEN_ ), DESCX( DLEN_ ), DESCAF( DLEN_)
     integer  DESCC(DLEN_),descb2(dlen_)
     real*8 ANORM,Xnorm,Bnorm,eps,resid
     integer locr,locc
     integer nq
     !      External Functions 
     integer  NUMROC
     real*8   PZLANGE,pdlamch

     include 'mpif.h'

     character*(MPI_MAX_PROCESSOR_NAME) name
     integer irlen

     !    write(6,*) iam,irowsize,icolsize,nrhs,N

     call MPI_GET_PROCESSOR_NAME(name,irlen,ierr)
     write(nout,'(a,i3,a,a)') ' Compute Basis proc=:',iam,'in ',name
     write(nout,'(a,8i6)') 'In lin : ',iam,irowsize,icolsize,n,nrhs,ictxt,myrow,mycol
     lda = max(1,irowsize)
     !    if(iam.eq.iprint) write(nout,*) "lda = ",lda         
     call flush(nout)
     LOCr = NUMROC( N, MB, MYROW, iRSRC, NPROW )
     LOCc = NUMROC( N, NB, MYCOL, iCSRC, NPCOL )
     !    if(iam.eq.iprint) write(nout,*) "locr = ",locr,irowsize
     !    if(iam.eq.iprint) write(nout,*) "locc = ",locc,icolsize

     !
     !     Initialize the array descriptor for the matrix A and B
     !
     call DESCINIT( DESCA,  N, N  , MB, NB, irsrc, icsrc, ICTXT, lda, INFO )
     call DESCINIT( DESCAF, N, N  , MB, NB, irsrc, icsrc, ICTXT, lda, INFO )
     call DESCINIT( DESCB, N, NRHS, MB, 1, irsrc, icsrc, ICTXT, lda, INFO )
     call DESCINIT( DESCX, N, NRHS, MB, 1, irsrc, icsrc, ICTXT, lda, INFO )
     !    if(iam.eq.iprint) then
     !        write(nout,'(10i5)') iam,desca(1:9)
     !        write(nout,'(10i5)') iam,descaf(1:9)
     !        write(nout,'(10i5)') iam,descx(1:9)
     !        write(nout,'(10i5)') iam,descb(1:9)
     call flush(nout)
     !    end if

     !
     !     Compute working array
     !
     allocate(ipiv(irowsize+mb),stat=ierr)
     if(ierr.ne.0) stop "Allocation of ipiv failed"
     allocate(row(irowsize),stat=ierr)
     if(ierr.ne.0) stop "Allocation of row failed"
     allocate(col(icolsize),stat=ierr)
     if(ierr.ne.0) stop "Allocation of col failed"
     nq = numroc(n,nb,mycol,icsrc,npcol)
     allocate(ferr(nq),berr(nq),stat=ierr)
     if(ierr.ne.0) stop "Allocation of ferr failed"    
     call flush(nout)

     !write(nout,'(a,5i6)') 'In lin : ',iam,irowsize,icolsize,nq,nrhs
     call flush(nout)
     call blacs_barrier(ictxt,'A')
     lwork = -1
     lrwork = -1
     !   call mpi_barrier(mpi_comm_world,ierr)
     call PZGESVX( 'E', 'N', n, nrhs, A,  1,  1,  DESCA,     &
          AF, 1,  1,  DESCAF,    &
          IPIV, EQUED, row, col,   &
          B,  1,  1, DESCB,        &
          X,  1,  1, DESCX,        &
          RCOND, FERR, BERR, cWORK, LWORK, sWORK, LRWORK, INFO )

     lwork = int(real(cwork(1)))
     if(iam.eq.iprint) write(nout,*) 'lwork =',lwork 
     call flush(nout)

     allocate(work(2*lwork))

     lrwork = int(swork(1))
     if(iam.eq.iprint) write(nout,*) 'lrwork =',lrwork
     call flush(nout)

     allocate(rwork(lrwork))

     if(iam.eq.iprint) then
         write(nout,*)  lwork,lrwork
         write(nout,'(9i5)') desca(1:9)
         write(nout,'(9i5)') descaf(1:9)
         write(nout,'(9i5)') descx(1:9)
         write(nout,'(9i5)') descb(1:9)
         call flush(nout)
     end if

     !
     !*********************************************************************
     !     Call expert ScaLAPACK PZGESVX routine
     !*********************************************************************
     !
     write(nout,*) "Starting PZGESVX..."

     call PZGESVX( 'E', 'N', n, nrhs, A,  1,  1,  DESCA,     &
          AF, 1,  1,  DESCAF,    &
          IPIV, EQUED, row, col,       &
          B,  1,  1, DESCB,        &
          X,  1,  1, DESCX,        &
          RCOND, FERR, BERR, WORK, LWORK, RWORK, LRWORK, INFO )
     write(nout,*) "... PZGESVX ended",INFO

     if(iam.eq.iprint) write(nout,*)  "INFO =",info," in proc",iam
     write(nout,*)  "Requested lwork =",lwork,        " lrwork=",lrwork
     write(nout,*)  "Optimal   lwork =",int(real(work(1)))," lrwork=",int(rwork(1)) 

     call flush(nout)
     ! Write the solution

!     call blacs_barrier(ictxt,'A')
!     call PZLAWRITE('sol.dat', N, NRHS, x, 1,1, DESCX,0, 0, WORK )

     call blacs_barrier(ictxt,'A')
     write(nout,*)  "deallocating in lin"
     if(allocated(ipiv)) deallocate(ipiv,stat=info)
     if(info.gt.0) write(nout,*) "deallocation of ipiv failed",info

     if(allocated(row)) deallocate(row,stat=info)
     if(info.gt.0) write(nout,*) "deallocation of row failed",info

     if(allocated(col)) deallocate(col,stat=info)
     if(info.gt.0) write(nout,*) "deallocation of col failed",info

     if(allocated(berr)) deallocate(berr,stat=info)
     if(info.gt.0) write(nout,*) "deallocation of berr failed",info
     if(allocated(ferr)) deallocate(ferr,stat=info)
     if(info.gt.0) write(nout,*) "deallocation of ferr failed",info

     if(allocated(work)) deallocate(work,stat=info)
     if(info.gt.0) write(nout,*) "deallocation of work failed",info
     if(allocated(rwork)) deallocate(rwork,stat=info)
     if(info.gt.0) write(nout,*) "deallocation of rwork failed",info
     write(nout,*)  "Finished deallocating in lin"
     call flush(nout)

     return
 end subroutine lin_pzgesv
