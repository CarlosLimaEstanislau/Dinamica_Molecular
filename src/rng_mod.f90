module rng_mod
	use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    
    private

    public :: random_normal, shuffle_array

    contains

    	function random_normal(mean, std) result(r)
	        implicit none
	        real(dp), intent(in) :: mean, std
	        real(dp) :: r
	        real(dp) :: u1, u2, fac
	        real(dp), save :: spare
	        logical,  save :: has_spare = .false.

	        if (has_spare) then
	            r = mean + std * spare
	            has_spare = .false.
	        else
	            call random_number(u1)
	            call random_number(u2)

	            if (u1 < tiny(1.0_dp)) u1 = tiny(1.0_dp)

	            fac = sqrt(-2.0_dp * log(u1))
	            r     = mean + std * fac * cos(2.0_dp * pi * u2)
	            spare =               fac * sin(2.0_dp * pi * u2)

	            has_spare = .true.
	        end if
	    end function random_normal

	    subroutine shuffle_array(arr1, arr2)
	        implicit none
	        real(dp), dimension(:), intent(inout) :: arr1, arr2
	        integer :: i, j
	        real(dp) :: rtmp, temp1, temp2

	        do i = size(arr1), 2, -1
	            call random_number(rtmp)
	            j = int(rtmp * i) + 1

	            temp1   = arr1(i)
	            arr1(i) = arr1(j)
	            arr1(j) = temp1

	           	temp2   = arr2(i)
	            arr2(i) = arr2(j)
	            arr2(j) = temp2
	        end do
	    end subroutine shuffle_array
end module rng_mod