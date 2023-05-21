!
!
!      Calculation of Temkin Poet wave function with Sturmians 
!      expansions
!
!
program sturmianTP	
    !........ Parallel Version
    !
    !	Inspired on D. Mitnik code, Jun-2009
    !   Basis calculation by A. L. Frapiccini 
    !   Scalapack, MPI debugging by F. D. Colavecchia
    !
    use debug_module
    use problemdata_module
    use parablock
    use crashblock
    use procdata
    use matrixdat
    implicit none
    real timei
    complex*16, allocatable :: H(:,:)
    complex*16, allocatable :: tind(:,:)
    complex*16, allocatable :: x(:,:)

    integer nn
    integer ierr
    integer icolsize,irowsize
    integer nrhs,nqrhs
    include 'mpif.h'

    character*(MPI_MAX_PROCESSOR_NAME) name
    integer irlen


    !  Start MPI
    !
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,iam,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)
    call MPI_GET_PROCESSOR_NAME(name,irlen,ierr)

    call set_debuglevel('a',iam)

    write(nout,'(a,i3,a,a)') ' Compute Basis proc=:',iam,'in ',name

    !      Initialization - Sturmian basis set and size
    call read_problemdata(E0,spin,ninterv,res,dump,iam)

    call BasisSet(NN,iam,nproc)
    call mpi_barrier(mpi_comm_world,ierr)

    call flush(nout)
    call timing(timei,iam,nproc,"Basis time: ")

    !      Initialization - Hamiltonian size - etc.
    if(iam.eq.iprint) write(nout,*) ' Initialization proc=:',iam

    call initialization(irowsize,icolsize,nqrhs,nrhs,NN)	

    !if(iam.eq.iprint) write(nout,*) irowsize,icolsize,nqrhs,NN 
    call flush(nout)

    !       Allocate the local Hamiltonian, x and RHS

    if(iam.eq.iprint) write(nout,*) "Allocating  Hamiltonian",irowsize,icolsize
    allocate(H(irowsize,icolsize),stat=ierr)
    if(ierr.ne.0)stop 'allocation fails - H'
    allocate(tind(irowsize,1),stat=ierr)
    if(ierr.ne.0)stop 'allocation fails - Tind'
    allocate(x(irowsize,1),stat=ierr)
    if(ierr.ne.0)stop 'allocation fails - x'
    call flush(nout)

    !       Create Hamiltonian
    if(iam.eq.iprint) write(nout,*) ' Computing Hamiltonian'

    call hamiltonian(H,irowsize,icolsize,NN)

    call timing(timei,iam,nproc,"Computing H time: ")

    !      Write Hamiltonian to disc as is
    if(dump) call dump_matrix('h',H,irowsize,icolsize,NN,NN)

    call blacs_barrier(ictxt,'A')

    call timing(timei,iam,nproc,"Writing H time: ")

    !       Create RHS
    if(iam.eq.iprint) write(nout,*),' Computing RHS'
    call flush(nout)

    call rhs(tind,irowsize,NN,nqrhs)

!      Write rhs to disc as is
    if(dump) call dump_matrix('r',tind,irowsize,nrhs,NN,nrhs)

    call blacs_barrier(ictxt,'A')

    call timing(timei,iam,nproc,"Writing R time: ")

    !      Free Sturmian basis
    call free_basis

    !       Solve the system   
    !    if(iam.eq.iprint) write(nout,*) ' Solve Linear System'

    write(nout,'(a,8i6)') 'In main : ',iam,irowsize,icolsize,NN,nrhs,ictxt,myrow,mycol
    call flush(nout)
    call blacs_barrier(ictxt,'A')

    call lin_pzgesv(H,tind,x,irowsize,icolsize,NN,nrhs)

    call timing(timei,iam,nproc,"Solve Linear System time: ")

    call dump_matrix('sol',x,irowsize,nrhs,NN,nrhs)

    if(allocated(H)) deallocate(H,stat=ierr)
    if(ierr.gt.0) write(nout,*) "deallocation of H failed",ierr

    if(allocated(x)) deallocate(x,stat=ierr)
    if(ierr.gt.0) write(nout,*) "deallocation of x failed",ierr

    if(allocated(tind)) deallocate(tind,stat=ierr)
    if(ierr.gt.0) write(nout,*) "deallocation of tind failed",ierr

    call blacs_gridexit(ictxt)
    call mpi_finalize(ierr)
    if(ierr>0) then 
        write(*,*) "MPI finishes with ierr",ierr
        stop
    end if
end program sturmianTP
!-------------------------------------------------------------------------------
!
!    Initialization of scalapack
!
!-------------------------------------------------------------------------------
subroutine initialization(irowsize,icolsize,nqrhs,nrhs,NN)	
    use debug_module
    use parablock
    use crashblock
    use procdata
    use matrixdat
    implicit real*8(a-h,o-z)
    include 'mpif.h'
    character*80 usrinfo

    !       Reading parameters for paralelization

    open(unit=10,file='scaldata.in')
    read(10,'(a)') usrinfo
    read(10,*) nrhs
    read(10,*) NB
    read(10,*) nprow
    read(10,*) npcol
    read(10,*) itype
    read(10,*) irsrc
    read(10,*) icsrc
    close(10)

    if(iam.eq.iprint) then
        write(nout,*)  'NRHS =',nrhs
        write(nout,*)  'NB = ',NB
        write(nout,*)  'NPROW = ',nprow
        write(nout,*)  'NPCOL = ',npcol
        write(nout,*)  'istype = ',itype
        write(nout,*)  'irsrc = ',irsrc
        write(nout,*)  'icsrc = ',icsrc
    end if

    t = nproc

    !write(nout,*)  'I am ',iam, 'nprow =',nprow,nn

    IF(NPROW.LE.0.OR.NPCOL.LE.0)THEN
        T=NPROC
        NPROW=NINT(SQRT(T))
        NPCOL=NPROW
        if(iam.eq.iprint) then
            write(nout,*)  'NPROW = ',nprow
            write(nout,*)  'NPCOL = ',npcol
        end if
    ENDIF

    if (nprow*npcol.ne.nproc) then
        if (iam.ne.0) stop 'USAGE: nprow*npcol must be equal to # of processors'
        istop=1
        return
    endif
    MB = 0
    if (MB.ne.NB) MB=NB                !FOR NOW

    !       initialize the BLACS context
    if(iam.eq.iprint) write(nout,*)  "Initializing BLACS context"

    call INITCTXT

    if(iam.eq.iprint) write(nout,'(a,i3,a,i3,a,i3)') 'myrow = ',myrow, " mycol = ", mycol, " on proc ",iam

    !       bails out if this process is not a part of this context
    if (myrow.eq.-1) return

    !       calculate the size of the local Hamiltonian
    if(iam.eq.iprint) write(nout,*)  "Computing  Hamiltonian local size"

    call HSIZE(irowsize,icolsize,NN)
    call flush(nout)

    !       calculate the size of the local RHS
    NQRHS =  NRHS

    return
end subroutine initialization

!
!     ******************************************************************
!
subroutine rhs(tind,irowsize,NN,nqrhs)	
    use debug_module
    use mapping_module
    use parablock
    use crashblock
    use procdata
    use matrixdat
    implicit none
    complex*16 tind(irowsize,nqrhs)
    complex*16 v

    integer nn,i
    integer ierr
    integer icolsize,irowsize
    integer nrhs,nqrhs
    integer irow,icol
    integer ipr,ipc
    include 'mpif.h'

    !       Calculates the RHS of the equation

    !        create RHS
    v = (0d0,0d0)
    do i=1,NN
        call mapping(i,1,irow,icol,ipr,ipc)
        if (ipr.eq.myrow.and.ipc.eq.mycol) then
            call vector_element_tp(i,v,iam)
            tind(irow,1) = v
        endif
    enddo
    call flush(6)
    return	
end subroutine rhs
!
!     ******************************************************************
!
subroutine hamiltonian(H,irowsize,icolsize,NN)	
    !> \brief Computes the local Hamiltonian matrix
    !
    !
    use debug_module
    use mapping_module
    use parablock
    use crashblock
    use procdata
    use matrixdat
    implicit real*8(a-h,o-z)
    complex*16 :: H(irowsize,icolsize) !< Local H Matrix 
    integer irowsize                   !< Local Row Size of H
    integer icolsize                   !< Local Column Size of H
    integer NN                         !< Global size of matrix H
    complex*16 matrix,overlap

    !       Calculates the Hamiltonian 
    include 'mpif.h'

    !        create Hamiltonian
    call random_seed
    do i=1,NN
        do j=1,NN
            call mapping(i,j,irow,icol,ipr,ipc)
            if (ipr.eq.myrow.and.ipc.eq.mycol) then
                call matrix_element_tp(i,j,matrix,overlap,iam)
                H(irow,icol) = matrix
            endif
        enddo
    enddo


    return	
end subroutine hamiltonian


!
! **********************************************************************
!
!			**** parallel ****
!--------------------------------------------------------------------------
!			subroutines for parallel code
!--------------------------------------------------------------------------

!
!     ******************************************************************
!
subroutine initctxt
    use parablock
    use crashblock
    use procdata
    use matrixdat
    implicit real*8(a-h,o-z)

    !       Initialize a single BLACS context

    include 'mpif.h'

    !       Initialize a single BLACS context
    CALL BLACS_GET( -1, 0, ictxt )
    CALL BLACS_GRIDINIT( ictxt, 'R', NPROW, NPCOL )
    CALL BLACS_GRIDINFO( ictxt, NPROW, NPCOL, MYROW, MYCOL )
    CALL BLACS_SET( ictxt, 15, 1)
    return
end subroutine initctxt

!
!***********************************************************************
!
subroutine Hsize(irowsize,icolsize,NN)
    use parablock
    use crashblock
    use procdata
    use matrixdat
    implicit real*8(a-h,o-z)

    !       calculate the size of the local Hamiltonian

    include 'mpif.h'

    !       dimension of the local matrix H
    irowsize = numroc(NN,mb,myrow,irsrc,nprow)
    icolsize = numroc(NN,nb,mycol,icsrc,npcol)
    return
end subroutine Hsize


!
!     ******************************************************************
!
subroutine printvalues(irowsize,icolsize,Z,W,Ndim)
    use parablock
    use crashblock
    use procdata
    use matrixdat
    implicit real*8(a-h,o-z)

    !       .print the eigenvalues and eigenvectors on file 

    parameter (iprocmin=1)
    parameter (ip1=100,ip2=ip1/iprocmin+1,ipch=3,np=ip1*ip2*ipch)
    data autoev/27.2113961/

    include 'mpif.h'

    dimension H(irowsize,icolsize)
    dimension Z(irowsize,icolsize)
    dimension W(Ndim)
    character*1 filec,filed,fileu

    !        open output files
    ich0 = ichar('0')
    ic = iam/100
    id = (iam - 100*ic)/10
    iu = (iam - 100*ic - 10*id)
    ic = ic + ich0
    id = id + ich0
    iu = iu + ich0
    filec = char(ic)
    filed = char(id)
    fileu = char(iu)
    open(unit=25,file='eigenvectors'//filec//filed//fileu//'.dat', &
         form='unformatted',status='unknown')
    open(unit=26,file='eigenvalues.dat',status='unknown')

    do irow=1,irowsize
        write(25) (Z(irow,icol), icol=1,icolsize)
    end do

    do i=1,Ndim
        if (iam.eq.iprint) write(26,*) i,'  ',W(i),'  ',W(i)*autoev
    end do

    return	
end subroutine printvalues


!
!***********************************************************************
!	
!			**** parallel ****



