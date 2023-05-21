module size_module
    implicit none

    integer l_max,imax
    integer l1_min,l2_min,l1_max,l2_max
    integer,allocatable :: ind_l(:,:),nmax(:)


contains
    !----------------------------------------------------------------------------------------	
    Subroutine sizebase(Lg)
 	use sturmians_base
 	use sturmians_data
        integer l1,l2,i,Lg,j
	integer pair
        integer,allocatable :: temp_l(:,:)


        allocate(temp_l(0:1,0:lmax(1)*lmax(2)))
        i=-1
  	do l1=0,lmax(1)
            do l2=0,lmax(2)
                if(abs(l1-l2).le.Lg.and.Lg.le.l1+l2.and.mod(l1+l2+Lg,2).eq.0.and.l1.le.l2)then
                    i=i+1
                    temp_l(0,i)=l1
                    temp_l(1,i)=l2
                endif
            enddo
 	enddo
        allocate(ind_l(0:1,0:i))
	l_max=size(ind_l,2)

	do j=0,i
            ind_l(:,j)=temp_l(:,j)
	enddo
	deallocate(temp_l)

	allocate(nmax(0:i))
	call sizebase_n(i,ind_l,nmax)
	imax=sum(nmax)

        l1_min=minval(ind_l(0,:),1)
        l1_max=maxval(ind_l(0,:),1)
	l2_min=minval(ind_l(1,:),1)
	l2_max=maxval(ind_l(1,:),1)

	return
    End Subroutine sizebase
    !----------------------------------------------------------------------------------------		
    Subroutine sizebase_n(lmx,ind,nmx)
	use sturmians_base
 	use sturmians_data
 	integer i,lmx,nmx(0:lmx),ind(0:1,0:lmx),l1,l2

	do i=0,lmx
            l1=ind(0,i)
            l2=ind(1,i)
            nmx(i)=norder(1,l1)*norder(2,l2)
 	enddo

 	return
    End Subroutine sizebase_n
    !----------------------------------------------------------------------------------------	
    subroutine indexado(in,l1,l2,n1,n2)
	use sturmians_data 
	use sturmians_base
	implicit none
	integer,intent(in) :: in
	integer i,n1,n2,j,l1,l2

	n1=int((1d0*in-1d0)/(norder(2,l2)*1d0)+1d0)
	n2=int(in-norder(2,l2)*(n1-1d0))

	return
    end subroutine indexado
    !----------------------------------------------------------------------------------------	
    subroutine index_j(ij,l1,l2,j1,j2)
	use sturmians_data 
	use sturmians_base
	implicit none
	integer,intent(in) :: ij
	integer l1,l2,j1,j2

	j1=int((1d0*ij-1d0)/(norder(2,l2)*1d0)+1d0)
	j2=int(ij-norder(2,l2)*(j1-1d0))

	return
    end subroutine index_j
    !----------------------------------------------------------------------------------------
end module size_module

