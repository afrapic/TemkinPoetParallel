!-------------------------------------------------------------------------------------------------------------
   subroutine radial_sturm(n,asnt,nu,ele,r,sturm)
	! 6 Nov. 2007
	! Calcula la funcion sturmiana de orden nu para momento angular ele para un r dado
	use sturmians_base
	use sturmians_data 
	use baseLag_module
	implicit none
	integer n,nu,ele,m
	real*8 r
	complex*16 sturm
	real*8,allocatable :: fi(:)
	character*1 asnt

	 allocate(fi(0:norder(n,ele)+nasint))
	 call basis(norder(n,ele)+nasint+1,ele,lambda(n,ele),r,fi)
	 sturm=(0d0,0d0)
	 do m=0,norder(n,ele)-1
	 sturm=sturm+bs(n)%a(m,nu,ele)*fi(m)
	 enddo

	if(asnt.eq.'y')then
	 do m=norder(n,ele),nasint+norder(n,ele)
	 sturm=sturm+bs(n)%coef_asint(m,ele)*fi(m) 
	 enddo
	endif
	 deallocate(fi)

	return
	end subroutine radial_sturm
