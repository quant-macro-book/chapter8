module mod_calibration
!  Purpose:
!	Collects all calibrated parameters.
!   All parameters are taken from Krusell and Smith (1998,JPE).
!
!  Record of revisions:
!     Date     Programmer  Description of change
!  ==========  ==========  =====================
!  09/13/2019  T. Yamada   Original

    use mod_types

	implicit none

	!****** preference ******
	real(dp), parameter :: beta  = 0.99 ! discount factor
	real(dp), parameter :: gamma = 1.0  ! relative risk aversion

	!****** production ******
	real(dp), parameter :: alpha = 0.36  ! capital share
	real(dp), parameter :: delta = 0.025 ! depreciation rate

	!****** idiosyncratic unemployment shock ******
	integer,  parameter :: ne = 2 ! employed, unemployed
	real(dp), parameter, dimension(ne) :: endow = (/1.0, 0.05/) ! ui = 5%

	!****** aggregate productivity shock and unemployment dynamics ******
	integer,  parameter :: nz = 2 ! #aggregate state = 2
	real(dp), parameter :: unempg = 0.04 ! unemployment rate when good
	real(dp), parameter :: unempb = 0.1  ! unemployment rate when bad
	real(dp), parameter :: durug  = 1.5  ! duration of unemployment when good
	real(dp), parameter :: durub  = 2.5  ! duration of unemployment when bad
	real(dp), parameter :: durgd  = 8.0  ! duration of boom (8 quarters)
	real(dp), parameter :: durbd  = 8.0  ! duration of recession (8 quarters)
	real(dp), parameter, dimension(nz) :: tfp = (/1.01, 0.99/) ! TFP shock

end module mod_calibration
