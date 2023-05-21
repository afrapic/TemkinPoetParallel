    !-----------------------------------------------------------
	subroutine sturmians_spline
	use sturmians_base
	use sturmians_data
	use array_module
	use problemdata_module
	implicit none
	integer i,ord,l
	real*8 interv,dh,drs0,dis0,drsf,disf
	real*8,allocatable :: ar(:),br(:),cr(:),dr(:),ai(:),bi(:),ci(:),di(:)
	real*8,allocatable :: s1_real(:),s1_imag(:)
	complex*16 sn,s0,s1,s2,im
	im=(0d0,1d0)

	nrmax(2)=nrmax(1)

	allocate(sa1(0:nrmax(1)-1,1:norder(1,0),0:lmax(1)),sb1(0:nrmax(1)-1,1:norder(1,0),0:lmax(1)))
	allocate(sc1(0:nrmax(1)-1,1:norder(1,0),0:lmax(1)),sd1(0:nrmax(1)-1,1:norder(1,0),0:lmax(1)))
	allocate(sa2(0:nrmax(2)-1,1:norder(2,0),0:lmax(2)),sb2(0:nrmax(2)-1,1:norder(2,0),0:lmax(2)))
	allocate(sc2(0:nrmax(2)-1,1:norder(2,0),0:lmax(2)),sd2(0:nrmax(2)-1,1:norder(2,0),0:lmax(2)))

	dh=1d-4
	do l=0,lmax(1)
	do ord=1,norder(1,l)
	!Calcula la derivada segunda en bs(1)%r(0) y bs(1)%r(nrmax-1)
	call radial_sturm(1,out,ord,l,bs(1)%r(0),s0)
	call radial_sturm(1,out,ord,l,bs(1)%r(0)+dh,s1)
	call radial_sturm(1,out,ord,l,bs(1)%r(0)+2d0*dh,s2)
	drs0=real(s2-2d0*s1+s0)/dh/dh
	dis0=aimag(s2-2d0*s1+s0)/dh/dh
	call radial_sturm(1,out,ord,l,bs(1)%r(nrmax(1)-1),s0)
	call radial_sturm(1,out,ord,l,bs(1)%r(nrmax(1)-1)+dh,s1)
	call radial_sturm(1,out,ord,l,bs(1)%r(nrmax(1)-1)+2d0*dh,s2)
	drsf=real(s2-2d0*s1+s0)/dh/dh
	disf=aimag(s2-2d0*s1+s0)/dh/dh
		!Guarda las sturmians para cada ord y l en funcion de r
		allocate(s1_real(0:nrmax(1)-1),s1_imag(0:nrmax(1)-1))
		do i=0,nrmax(1)-1
		call radial_sturm(1,out,ord,l,bs(1)%r(i),sn)
		s1_real(i)=real(sn)
		s1_imag(i)=aimag(sn)
		enddo
	!Spline de la parte real e imaginara para ord y l 
	allocate(ar(0:nrmax(1)-1),br(0:nrmax(1)-1),cr(0:nrmax(1)-1),dr(0:nrmax(1)-1))
	allocate(ai(0:nrmax(1)-1),bi(0:nrmax(1)-1),ci(0:nrmax(1)-1),di(0:nrmax(1)-1))
	call spline(bs(1)%r,s1_real,ar,br,cr,dr,drs0,drsf,nrmax(1))
	call spline(bs(1)%r,s1_imag,ai,bi,ci,di,dis0,disf,nrmax(1))
	sa1(:,ord,l)=ar(:)+im*ai(:)
	sb1(:,ord,l)=br(:)+im*bi(:)
	sc1(:,ord,l)=cr(:)+im*ci(:)
	sd1(:,ord,l)=dr(:)+im*di(:)
	deallocate(s1_real,s1_imag,ar,br,cr,dr,ai,bi,ci,di)
	enddo
	enddo

	if(res.eq.'eq')then
	sa2=sa1
	sb2=sb1
	sc2=sc1
	sd2=sd1
	else
		interv=rmax(2)/(nrmax(2)*1d0-1d0)
		do i=1,nrmax(2)
		bs(2)%r(i)=i*interv
		enddo
		do l=0,lmax(2)
		do ord=1,norder(2,l)
		!Calcula la derivada segunda en bs(1)%r(0) y bs(1)%r(nrmax-1)
		call radial_sturm(2,out,ord,l,bs(2)%r(0),s0)
		call radial_sturm(2,out,ord,l,bs(2)%r(0)+dh,s1)
		call radial_sturm(2,out,ord,l,bs(2)%r(0)+2d0*dh,s2)
		drs0=real(s2-2d0*s1+s0)/dh/dh
		dis0=aimag(s2-2d0*s1+s0)/dh/dh
		call radial_sturm(2,out,ord,l,bs(2)%r(nrmax(2)-1),s0)
		call radial_sturm(2,out,ord,l,bs(2)%r(nrmax(2)-1)+dh,s1)
		call radial_sturm(2,out,ord,l,bs(2)%r(nrmax(2)-1)+2d0*dh,s2)
		drsf=real(s2-2d0*s1+s0)/dh/dh
		disf=aimag(s2-2d0*s1+s0)/dh/dh
			!Guarda las sturmians para cada ord y l en funcion de r
			allocate(s1_real(0:nrmax(2)-1),s1_imag(0:nrmax(2)-1))
			do i=0,nrmax(2)-1
			call radial_sturm(2,out,ord,l,bs(2)%r(i),sn)
			s1_real(i)=real(sn)
			s1_imag(i)=aimag(sn)
			enddo
		!Spline de la parte real e imaginara para ord y l 
		allocate(ar(0:nrmax(2)-1),br(0:nrmax(2)-1),cr(0:nrmax(2)-1),dr(0:nrmax(2)-1))
		allocate(ai(0:nrmax(2)-1),bi(0:nrmax(2)-1),ci(0:nrmax(2)-1),di(0:nrmax(2)-1))
		call spline(bs(2)%r,s1_real,ar,br,cr,dr,drs0,drsf,nrmax(1))
		call spline(bs(2)%r,s1_imag,ai,bi,ci,di,dis0,disf,nrmax(1))
		sa2(:,ord,l)=ar(:)+im*ai(:)
		sb2(:,ord,l)=br(:)+im*bi(:)
		sc2(:,ord,l)=cr(:)+im*ci(:)
		sd2(:,ord,l)=dr(:)+im*di(:)
		deallocate(s1_real,s1_imag,ar,br,cr,dr,ai,bi,ci,di)
		enddo
		enddo
	endif

	return
	end subroutine sturmians_spline
