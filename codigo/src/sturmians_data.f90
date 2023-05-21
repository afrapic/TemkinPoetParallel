module sturmians_data
    implicit none
    integer npoints,nasint,npoints_in
    integer lmax(2),nrmax(2)
    integer, parameter :: NEL=2
    integer, parameter :: NLMAX = 10
    integer,allocatable :: norder(:,:)
    real*8 rfin,rfin_in
    real*8 E(2),rmax(2),Z(2)
    real*8,allocatable :: lambda(:,:),alfa(:,:),Zi(:,:),r0(:,:)
    character*3 pot(2),wave(2)
    character*1 out,nrm

contains

    subroutine read_input(i,E_i,Z_i,lmax_i,pot_i,nrmax_i,rmax_i,wave_i,nrm_i,out_i,nasint_i)    
	use file_module
	implicit none
	character ich
	character*3 pot_i,wave_i
	real*8 E_i,Z_i,rmax_i
	integer i,lmax_i,nrmax_i,l,nasint_i
	character*1 nrm_i,out_i
	integer, parameter :: unit=1000
	integer, external :: read_int
	real*8, external :: read_double

	write(ich,'(i1.1)') i

	open(unit=unit,file="electron"//ich//".in")
 	E_i      = read_double(unit,var,dummy)
	Z_i      = read_double(unit,var,dummy)
	lmax_i   = read_int(unit,var,dummy)
	call read_string(unit,var,dummy,pot_i)
	rmax_i   = read_double(unit,var,dummy)
	nrmax_i  = read_int(unit,var,dummy)
	if(E_i.gt.0d0)then
            call read_string(unit,var,dummy,wave_i)
	else
            wave_i='xxx'
	endif
	call read_string(unit,var,dummy,nrm_i)
	call read_string(unit,var,dummy,out_i)
	nasint_i   = read_int(unit,var,dummy)
 	close(unit)

	return
    end subroutine read_input

    subroutine read_ele(i) !Reads the parameters lambda,norder and ntot for each l
	implicit none
	integer i,l
	character ich
        integer mxlmax

	if(i.eq.1)then
            mxlmax = maxval(lmax)
            allocate(lambda(NEL,0:mxlmax),norder(NEL,0:mxlmax))
            if(pot(i).eq.'yuk')then
                allocate(alfa(NEL,0:mxlmax))
                alfa=0d0
            endif
            if(pot(i).eq.'cou')then
                allocate(Zi(NEL,0:mxlmax))
                Zi=0d0
            endif
            if(pot(i).eq.'exp')then
                allocate(alfa(NEL,0:mxlmax))
                alfa=0d0
            endif
            if(pot(i).eq.'hul')then
                allocate(alfa(NEL,0:mxlmax))
                alfa=0d0
            endif
            if(pot(i).eq.'cbx')then
                allocate(r0(NEL,0:mxlmax))
                r0=0d0
            endif
            if(pot(i).eq.'vbx')then
                allocate(r0(NEL,0:mxlmax))
                r0=0d0
            endif
            lambda=0d0
            norder=0
	endif
	write(ich,'(i1.1)') i
	open(unit=100,file="lambda"//ich//".in")
	open(unit=200,file="norder"//ich//".in")
	open(unit=400,file="potdata"//ich//".in")
	do l=0,lmax(i)
            read(100,*) lambda(i,l)
            read(200,*) norder(i,l)
            if(pot(i).eq.'yuk')then
                read(400,*) alfa(i,l)
            endif
            if(pot(i).eq.'cou')then
                read(400,*) Zi(i,l)
            endif
            if(pot(i).eq.'exp')then
                read(400,*) alfa(i,l)
            endif
            if(pot(i).eq.'hul')then
                read(400,*) alfa(i,l)
            endif
            if(pot(i).eq.'cbx')then
                read(400,*) r0(i,l)
            endif
            if(pot(i).eq.'vbx')then
                read(400,*) r0(i,l)
            endif
	enddo
	close(100)
	close(200)
	close(400)

	return
    end subroutine read_ele

    subroutine print_electron_data(i,E_i,Z_i,lmax_i,pot_i,norder_i,lambda_i,wave_i)
	implicit none
	real*8 E_i,Z_i,lambda_i(0:NLMAX)
	integer i,lmax_i,norder_i(0:NLMAX),j
	character*3 pot_i,wave_i
	logical,external :: LSAMEN

	write(*,*)
	write(*,'(a)') "Finite Laguerre basis method" 
	write(*,'(a)')        "--------------"
	write(*,'(a,i1)') " Electron ",i
	write(*,'(a)')        "--------------"
	write(*,'(a,i1,a,f9.5)') "     E(",i,") = ",E_i
	write(*,'(a,i1,a,f9.5)') "     Z(",i,") = ",Z_i
	write(*,'(a,i1,a,i3)')    " lmax(",i,") = ",lmax_i
	if(pot_i.eq.'yuk')then 
            write(*,'(a)')    "Short range potential is Yukawa"
	endif
	if(pot_i.eq.'exp')then 
            write(*,'(a)')    "Short range potential is exponential"
	endif
	if(pot_i.eq.'cou')then 
            write(*,'(a)')    "Coulomb Sturmians"
 	endif
	if(pot_i.eq.'hul')then 
            write(*,'(a)')    "Short range potential is Hulten"
 	endif
	if(pot_i.eq.'cbx')then 
            write(*,'(a)')    "Short range potential is Coulomb box"
 	endif
	if(pot_i.eq.'vbx')then 
            write(*,'(a)')    "Short range potential is spherical box"
 	endif
	if(E_i.gt.0d0)then
            if(wave_i.eq.'out')then
                write(*,'(a)')    "Outgoing wave asymptotic behavior"
            endif
            if(wave_i.eq.'sta')then
                write(*,'(a)')    "Stationary asymptotic behavior"
            endif
	endif

	do j=0,lmax_i
            write(*,'(a,i1,a,i1,a,i3)')    "norder(",i,", l=",j,") = ",norder_i(j)
            write(*,'(a,i1,a,i1,a,f9.5)')  "lambda(",i,", l=",j,") = ",lambda_i(j)
            if(pot_i.eq.'yuk')then
                write(*,'(a,i1,a,i1,a,f9.5)')     "alfa(",i,",l=",j,") = ",alfa(i,j)
            endif
            if(pot_i.eq.'exp')then
                write(*,'(a,i1,a,i1,a,f9.5)')    "alfa(",i,",l=",j,") = ",alfa(i,j)
            endif
            if(pot_i.eq.'cou')then
                write(*,'(a,i1,a,i1,a,f9.5)')      "Zi(",i,",l=",j,") = ",Zi(i,j)
            endif
            if(pot_i.eq.'hul')then
                write(*,'(a,i1,a,i1,a,f9.5)')    "alfa(",i,",l=",j,") = ",alfa(i,j)
            endif
            if(pot_i.eq.'cbx')then
                write(*,'(a,i1,a,i1,a,f9.5)')    "r0(",i,",l=",j,") = ",r0(i,j)
            endif
            if(pot_i.eq.'vbx')then
                write(*,'(a,i1,a,i1,a,f9.5)')    "r0(",i,",l=",j,") = ",r0(i,j)
            endif
	enddo

   	return
    end subroutine print_electron_data

    subroutine free_sturmians_data
        implicit none
        if(allocated(lambda)) deallocate(lambda)
        if(allocated(r0)) deallocate(r0)
        if(allocated(alfa)) deallocate(alfa)
        if(allocated(Zi)) deallocate(Zi)
        if(allocated(norder)) deallocate(norder)
    end subroutine free_sturmians_data

end module sturmians_data
