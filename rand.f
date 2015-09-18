	subroutine rand(n, values)
		implicit none
		
		integer, intent(in) :: n
		double precision, dimension(n) :: values

		integer, parameter :: idist = 1
		integer, dimension(4) :: iseed = (/0, 0, 0, 0/)

		call SYSTEM_CLOCK(iseed(4))
		iseed(4) = iseed(4) * 2 + 1

		call dlarnv(idist, iseed, n, values)

	end subroutine rand