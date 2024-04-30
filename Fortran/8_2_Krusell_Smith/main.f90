program main
!  Purpose:
!	Solve Krusell and Smith (1998) model
!
!  Record of revisions:
!     Date     Programmer  Description of change
!  ==========  ==========  =====================
!  09/13/2019  T. Yamada   Original

    ! use IMSL library for regression
	include 'link_fnl_shared.h'

	use mod_types
	use mod_calibration
	use mod_variables
	use mod_generate_grid
    use mod_trans_probability
	use mod_EGM
    use mod_law_of_motion

	implicit none

	!****** variables for main loop ******
	integer :: it, z
	real(dp) :: metric1, metric2
	real(dp) :: metric = 1.0
	real(dp) :: adj = 0.5
	integer, parameter :: maxit = 1000
	real(dp), parameter :: toler = 1.0e-004
    
    !****** policy function ******
	real(dp), dimension(nd, ne, nk, nz) :: policy

    !****** simulated paths ******
	real(dp), dimension(nums) :: k_path
	integer,  dimension(nums) :: z_path

    !****** regression coefficients ******
	real(dp), dimension(nz) :: icept, icept1, icept2, R2
    real(dp), dimension(nz) :: slope, slope1, slope2 
	!-------------------------------------

	call cpu_time(time_begin)

    write (*, *)
    write (*, *) "-+-+-+- Solving Krusell and Smith model -+-+-+-"
    write (*, *) ""
    
    write (*, *) "--- calibrated parameters ---"
    write (*,"(' discount factor:            ', f10.4)") beta
    write (*,"(' relative risk aversion:     ', f10.4)") gamma
    write (*,"(' capital share:              ', f10.4)") alpha
    write (*,"(' depreciation rate:          ', f10.4)") delta
    write (*,"(' TFP level (g):              ', f10.4)") tfp(1)
    write (*,"(' TFP level (b):              ', f10.4)") tfp(2)
    write (*,"(' unemployment rate (g):      ', f10.4)") unempg
    write (*,"(' unemployment rate (b):      ', f10.4)") unempb
    write (*,"(' duration of unemp period (g):', f10.4)") durug
    write (*,"(' duration of unemp period (b):', f10.4)") durub
    write (*, *) ""

    !pause

	! uniform grid
	!call grid_uniform(amin, amax, na, aprime)
	!call grid_uniform(amin, amax, nd, grid)
    call grid_uniform(kmin, kmax, nk, kgrid)

	! exponential grid
	call grid_triple_exp(amin, amax, na, aprime)
	call grid_triple_exp(amin, amax, nd, grid)

	! transition probability matrix: prob(e,z,e',z')
	call trans_prob(prob, tran_zz)

    ! debug: simulation of transition probability matrix
    !call test_prob_zzee

	! aggragate labor supply
	pi_dist(:, 1) = (/1.0-unempg, unempg/) ! boom
	pi_dist(:, 2) = (/1.0-unempb, unempb/) ! recession
	do z = 1,nz
		L_agg(z) = dot_product(pi_dist(:, z), endow)
	end do
    !L_agg(1) = 1.0 - unempg
    !L_agg(2) = 1.0 - unempb

	! initial guess
	icept = 0.0
	slope = 1.0
	!icept(1) = 0.1346
	!icept(2) = 0.1227
	!slope(1) = 0.9631
	!slope(2) = 0.9649

	!****** main loop ******
	do it = 1, maxit

		if (metric < toler) exit

		write(*,*) "iteration counter:", it

		! get policy function by the Endogenous Gridpoints Method
		call end_grid_method(icept, slope, policy)

		! law of motion: Young (2010)
		!call law_of_motion(policy, k_path, z_path)

		! law of motion by simulation: originally in Krusell and Smith (1998)
		call law_of_motion_sim(policy, k_path, z_path)

		! new coefficients from regression
		call regress(k_path, z_path, icept1, slope1, R2)

		! metric of error
		metric1 = maxval(abs((icept-icept1)/icept))
		metric2 = maxval(abs((slope-slope1)/slope))
		metric  = max(metric1, metric2)

		! update coefficients
		icept2 = adj*icept + (1.0-adj)*icept1
		slope2 = adj*slope + (1.0-adj)*slope1
		icept  = icept2
		slope  = slope2

        write(*,"(' error:                   ', f10.4,'%')") metric*100.0
		write(*,"(' intercept:               ', f10.4)") icept
		write(*,"(' slope:                   ', f10.4)") slope
		write(*,"(' R^2:                     ', f10.4)") R2
		write(*,*) ""

	end do

	call cpu_time(time_end)

	open(10, file="result.txt")
		write(10, *) icept
		write(10, *) slope
		write(10, *) R2
		write(10, *) "program finished sucessfully in", time_end-time_begin, "seconds"
	close(10)

	!****** display results ******
	write(*,*) ""
	write(*,*) "program finished sucessfully in", time_end-time_begin, "seconds"

end program main
