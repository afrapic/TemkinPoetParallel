module problemdata_module
    implicit none

    integer ninterv,spin
    integer,parameter :: Lg=0
    real*8,parameter :: ev=27.2113962d0,Eb=-0.5d0
    real*8 E0,Etot,ratio,k0
    character*2 res
    logical dump,dumpcs,dumpf,flux

contains

    subroutine read_problemdata(En,s,interv,rs,dumpi,iam)
        use debug_module
	use file_module
	use sturmians_data
	integer s,interv
	real*8 En
	character*2 rs
        character*3 dum
        logical dumpi
	integer, parameter :: unit=1000
	integer, external :: read_int
	real*8, external :: read_double
        integer, optional :: iam

	dumpi=.false.
	open(unit=unit,file="problemdataparallel.in")
	En      = read_double(unit,var,dummy)      !Energia de incidencia del electrón
	s       = read_int(unit,var,dummy)         !spin total del sistema spin=0 (singlet) o spin=1 (triplet)
	interv  = read_int(unit,var,dummy)         !Número de puntos a calcular en la seccion eficaz
	call read_string(unit,var,dummy,rs)
	ratio      = read_double(unit,var,dummy) 
        call read_string(unit,var,dummy,dum)
        if(dum.eq.'yes') dumpi = .true.

        k0 = sqrt(2d0*En)

        if(present(iam)) then
            if(iam.eq.iprint) then
                write(nout,*) '--------------------------------'
                write(nout,'(a)') 'Temkin-Poet e-H'
                write(nout,'(a3,1x,f5.2)') 'E0:',En
                if(s.eq.0)then
                    write(nout,'(a)') 'Singlet'
                else
                    write(nout,'(a)') 'Triplet'
                endif
                if(rs.eq.'eq')then
                    write(nout,'(a)') 'Equal energy basis'
                endif
                write(nout,*) '--------------------------------'
            end if
        end if
	return
    end subroutine read_problemdata

    subroutine read_cs(dumpi,dumpf,flux,meth,fin,rr,iam)
        use debug_module
	use file_module
        integer meth
	character*2 fin
	real*8 rr
        character*3 dumc,dumf,dumg
        logical dumpi,dumpf,flux
	integer, parameter :: unit=1000
        integer, optional :: iam
	integer, external :: read_int
	real*8, external :: read_double

        dumpi=.false.
 	dumpf=.false.
        flux=.false.
	open(unit=unit,file="csparallel.in")
	call read_string(unit,var,dummy,dumf)
	call read_string(unit,var,dummy,dumg)
        call read_string(unit,var,dummy,dumc)
	meth    = read_int(unit,var,dummy)
	call read_string(unit,var,dummy,fin)
	rr      = read_double(unit,var,dummy)
        if(dumc.eq.'yes') dumpi = .true.
	if(dumf.eq.'yes') dumpf = .true.
	if(dumg.eq.'yes') flux = .true.

	return
    end subroutine read_cs


    subroutine free_basis
        use sturmians_base
        implicit none

        call free_base
        call free_sturmians_data
!        write(*,*) "Sturmians basis is free"
        return
    end subroutine free_basis


end module problemdata_module
