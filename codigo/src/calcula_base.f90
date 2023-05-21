    Subroutine calcula_base(i,bs_i,ele)
	!Calcula la base bs(i) de acuerdo al método
	use sturmians_data
	use sturmians_base
 	use baseLag_module
 	implicit none

	type(base_lag) :: bs_i
	integer i,ele

	call base_Laguerre(i,bs_i,ele)

	return
    end subroutine calcula_base
