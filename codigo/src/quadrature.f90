    subroutine quad_rule
	use sturmians_base
	use sturmians_data
	use gauleg_module
        use problemdata_module
	implicit none
	integer j1,k,np,i
	real*8 xi
	real*8,allocatable :: xx(:),ww(:)
	complex*16 sj1,sj2

	allocate(xx(npoints),ww(npoints))
	call gauleg(0d0,rfin,xx,ww)
	allocate(bs(1)%x(npoints),bs(2)%x(npoints),bs(1)%w(npoints),bs(2)%w(npoints))
	bs(1)%x=xx
	bs(2)%x=xx
	bs(1)%w=ww
	bs(2)%w=ww

	np=int(npoints+(npoints-1d0)*3d0)

	allocate(bs(1)%s(0:np-1,1:norder(1,0),0:lmax(1)))
        allocate(bs(2)%s(0:np-1,1:norder(2,0),0:lmax(2)))
	do k=0,np-1
	    if(k.ne.np-1)then
            i=int(k/4d0)+1d0
            xi=xx(i)+mod(k,4)*(xx(i+1)-xx(i))/4d0
	    else
	    xi=xx(npoints)
	    endif
            do j1=1,norder(1,0)
		call radial_sturm(1,out,j1,0,xi,sj1)
		bs(1)%s(k,j1,0)=sj1
                if(res.eq.'eq')then
                    bs(2)%s(k,j1,0)=sj1
                else
                    call radial_sturm(2,out,j1,0,xi,sj2)
                    bs(2)%s(k,j1,0)=sj2
                endif
            enddo
	enddo

	deallocate(xx,ww)
	return
    end subroutine quad_rule
!------------------------------------------------------------------------------------------------------
    subroutine quad_rule_in
	use sturmians_base
	use sturmians_data
	use gauleg_module
        use problemdata_module
	implicit none
	integer j1,k,np,i
	real*8 xi
	real*8,allocatable :: xx(:),ww(:)
	complex*16 sj1,sj2

	allocate(xx(npoints_in),ww(npoints_in))
	call gauleg(0d0,rfin_in,xx,ww)
	allocate(bs(1)%x_in(npoints_in),bs(2)%x_in(npoints_in),bs(1)%w_in(npoints),bs(2)%w_in(npoints))
	bs(1)%x_in=xx
	bs(2)%x_in=xx
	bs(1)%w_in=ww
	bs(2)%w_in=ww

	allocate(bs(1)%s_in(npoints_in,1:norder(1,0),0:lmax(1)))
        allocate(bs(2)%s_in(npoints_in,1:norder(2,0),0:lmax(2)))
	do k=1,npoints_in
            do j1=1,norder(1,0)
		call radial_sturm(1,'n',j1,0,xx(k),sj1)
		bs(1)%s_in(k,j1,0)=sj1
                if(res.eq.'eq')then
                    bs(2)%s_in(k,j1,0)=sj1
                else
                    call radial_sturm(2,'n',j1,0,xx(k),sj2)
                    bs(2)%s_in(k,j1,0)=sj2
                endif
            enddo
	enddo

	deallocate(xx,ww)
	return
    end subroutine quad_rule_in
