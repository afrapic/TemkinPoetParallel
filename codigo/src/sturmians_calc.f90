module sturmians_calc
    !Calculos sobre la base
    use sturmians_data
    use sturmians_base 
    use complex_module
contains

    !----------------------------------------------------------
    complex*16 function sturm_expansion_1D(n,cf,nf,r)
        !
        ! Calcula la expansion de una funcion radial en la base de sturmianas radiales 
        ! del electron n en un dado r.
        ! La función viene dada por el set de nf coeficientes cf, obtenidos previamente
        ! OJO: Se asume l=0 : ele=0 al llamar radial_sturmian
        !
        use complex_module
        use sturmians_base, only : NEL
        implicit none
        complex*16 cf(n)
        integer n,nf
        real*8 r
        integer i
        complex*16 st,tot
        tot = c0
        if(n.gt.NEL) stop "sturm_expansion: n > NEL"
        do i=1,n
            call radial_sturm(n,i,0,r,st)
            tot = tot+cf(i)*st
        end do
        sturm_expansion_1D = tot
        return

    end function sturm_expansion_1D

    !--------------------------------------------------------------------------------------------
    subroutine clausura(ele,filename)
	use sturmians_data
	use sturmians_base
	Implicit none
	integer ele,i,ip,n
	complex*16 cls,v,sr,srp
	character*80 filename

	open(unit=200,file='clausura_'//TRIM(filename)//'.dat')
	ip=160

	do i=1,nrmax(1)-1
            if(pot(1).eq.'hul')then
                v=exp(-bs(1)%r(ip)/alfa(1,ele))/(1d0-exp(-bs(1)%r(ip)/alfa(1,ele)))
            endif
            if(pot(1).eq.'yuk')then
                v=(-1d0)*exp(-alfa(1,ele)*bs(1)%r(ip))/(bs(1)%r(ip))
            endif
            if(pot(1).eq.'exp')then
                v=(-1d0)*exp(-alfa(1,ele)*bs(1)%r(ip))
            endif
            if(pot(1).eq.'cou')then
                v=-Zi(1,ele)/bs(1)%r(ip)
            endif
            if(pot(1).eq.'cbx')then
                if(bs(1)%r(ip).le.r0(1,ele))then
                    v=(-1d0)/(bs(1)%r(ip))
                else 
                    v=(0d0,0d0)
                endif
            endif
            if(pot(1).eq.'vbx')then
                if(bs(1)%r(ip).le.r0(1,ele))then
                    v=(-1d0)
                else 
                    v=(0d0,0d0)
                endif
            endif
            cls=(0d0,0d0)
            do n=1,norder(1,ele)
		call radial_sturm(1,n,ele,bs(1)%r(i),sr)
		call radial_sturm(1,n,ele,bs(1)%r(ip),srp)
		cls=cls+srp*sr*v
            enddo
            write(200,'(f13.6,1x,e13.6,1x,e13.6,1x,f13.6)') bs(1)%r(i),cls,bs(1)%r(ip)
	enddo

	close(200)

	return
    end subroutine clausura
    !-----------------------------------------------------------------------------------------------------
    Subroutine print_eigenvalues(num,l,filename)
	use sturmians_data
	use sturmians_base
	use inf_nan_detection
	Implicit none
	integer num,j,l,i
	integer,allocatable :: indx(:)
	complex*16 sn,norm,im
	character*2 nl
	character*80 filename
	im=(0d0,1d0)

	allocate(indx(1))
	indx=minloc(abs(bs(1)%beta))

	open(unit=200,file='beta_'//TRIM(filename)//'.dat')
	do j=1,num !-1,20,1
            write(200,'(f20.10,1x,f20.10)') bs(1)%beta(j,l) 
	enddo
	close(200)

	deallocate(indx)

	return
    end subroutine print_eigenvalues

end module sturmians_calc
