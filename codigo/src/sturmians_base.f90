module sturmians_base

    use sturmians_data
    !
    !   Derived type for one electron basis
    !
    type base_lag
       ! radial grid
       real*8, allocatable :: r(:),w(:),x(:),x_in(:),w_in(:)
       ! eigenvalues
       complex*16, allocatable :: beta(:,:)
       ! Laguerre coefficients
       complex*16, allocatable :: a(:,:,:),s(:,:,:),coef_asint(:,:),s_in(:,:,:)
    end type base_lag

    type(base_lag) :: bs(NEL)

contains

    subroutine create_base(i,b_i)
        use sturmians_data
        implicit none
        type (base_lag) :: b_i
        integer, intent(in) :: i !#electron
	integer nord

	nord=maxval(norder(i,:))

	allocate(b_i%r(0:nrmax(i)-1))
        allocate(b_i%beta(1:nord,0:lmax(i)))
        allocate(b_i%a(0:nord-1,1:nord,0:lmax(i)))
	allocate(b_i%coef_asint(0:nord+nasint-1,0:lmax(i)))

        return
    end subroutine create_base

    subroutine free_base
        use sturmians_data
        implicit none
        type (base_lag) :: b_i
        integer i !#electron
	integer nord
        do i=1,NEL
            b_i = bs(i)
            if(allocated(b_i%r)) deallocate(b_i%r)
            if(allocated(b_i%beta)) deallocate(b_i%beta)
            if(allocated(b_i%a)) deallocate(b_i%a)
            if(allocated(b_i%s)) deallocate(b_i%s)
	    if(allocated(b_i%s_in)) deallocate(b_i%s_in)
            if(allocated(b_i%r)) deallocate(b_i%r)
            if(allocated(b_i%w)) deallocate(b_i%w)
            if(allocated(b_i%x)) deallocate(b_i%x)
	    if(allocated(b_i%w_in)) deallocate(b_i%w_in)
            if(allocated(b_i%x_in)) deallocate(b_i%x_in)
	    if(allocated(b_i%coef_asint)) deallocate(b_i%coef_asint)
        end do

        call free_sturmians_data

        return
    end subroutine free_base



end module sturmians_base
