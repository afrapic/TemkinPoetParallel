subroutine BasisSet(nn,iam,nproc)
    use debug_module
    use sturmians_data
    use sturmians_base
    use size_module
    use problemdata_module
    implicit none
    integer i,nn,j,np
    integer, optional :: iam 
    integer, optional :: nproc
    real time1,time2

    Etot=E0+Eb

    !Calcula la base (igual energía)
    do i=1,NEL
        call read_input(i,E(i),Z(i),lmax(i),pot(i),nrmax(i),rmax(i),wave(i),nrm,out,nasint)
	if(res.eq.'eq')then
            E(1)=Etot
            E(2)=E(1)
	else
            E(1)=Etot*ratio
            E(2)=Etot-E(1)
	endif
	call read_ele(i)   
        call create_base(i,bs(i))
	if(i.eq.1)then
            if(present(iam)) then 
                if(iam.eq.iprint) then
                    call print_electron_data(i,E(i),Z(i),lmax(i),pot(i),norder(i,:),lambda(i,:),wave(i)) 
                end if
            end if
            call calcula_base(i,bs(i),0)
	else
            if(res.eq.'ne')then
                if(present(iam)) then 
                    if(iam.eq.iprint) then
                        call print_electron_data(i,E(i),Z(i),lmax(i),pot(i),norder(i,:),lambda(i,:),wave(i)) 
                    end if
                end if
                call calcula_base(i,bs(i),0)
            else
                bs(2)=bs(1)
            endif
	endif
    enddo


    if(present(iam)) then 
        if(iam.eq.iprint) write(*,*) 'Basis solved'
    end if

    call sizebase(Lg)
    if(present(iam)) then 
        if(iam.eq.iprint) write(*,'(a11,1x,i6.6)') 'Total size',imax
    end if

    nn=imax

    call cpu_time(time1)
    if(out.eq.'n')then
    np=norder(1,0)
    else
    np=norder(1,0)+nasint
    endif
    npoints=int(np*3d0+10d0)
    rfin=1d0*np
    call quad_rule
    npoints_in=int(norder(1,0)*3d0+10d0)
    rfin_in=1d0*norder(1,0)
    call quad_rule_in
    if(present(iam)) then 
        if(iam.eq.iprint) call timing(time1,iam,0,"Quadrature time: ")
    end if

    return
end subroutine BasisSet

