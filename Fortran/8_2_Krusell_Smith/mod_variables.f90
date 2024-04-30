module mod_variables
!  Purpose
!	Collect all global variables.
!
!  Record of revisions:
!     Date     Programmer  Description of change
!  ==========  ==========  =====================
!  09/13/2019  T. Yamada   Original

    use mod_types
	use mod_calibration, only:ne, nz

    implicit none

	!****** grid ******
	integer,  parameter :: na = 21 ! 101
	integer,  parameter :: nd = 501 ! 1001
	integer,  parameter :: nk = 6
	real(dp), parameter :: amax  = 300.0
	real(dp), parameter :: amin  = 0.0
	real(dp), parameter :: kmax = 40.0
	real(dp), parameter :: kmin = 30.0
	real(dp), dimension(na) :: aprime
	real(dp), dimension(nd) :: grid
	real(dp), dimension(nk) :: kgrid

	!****** transition probability ******
	real(dp), dimension(ne, nz, ne, nz) :: prob ! pr(e,z,e',z')
	real(dp), dimension(ne, nz) :: pi_dist
	real(dp), dimension(nz, nz) :: tran_zz ! pr(z,z')

	!****** macro variables ******
	real(dp), dimension(nz) :: L_agg

	!****** simulation part ******
	integer, parameter :: nums = 6000 ! 11000
    integer, parameter :: numi = 1000 ! 5000

	!****** computation time ******
	real(dp) time_begin, time_end

end module mod_variables
