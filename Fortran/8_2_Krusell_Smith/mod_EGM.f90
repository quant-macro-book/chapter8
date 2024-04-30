module mod_EGM
!  Purpose:
!	Compute a policy function of Krusell & Smith model
!	by the Endogenous Gridpoints Method.
!
!  Record of revisions:
!	  Date     Programmer  Description of change
!  ==========  ==========  =====================
!  09/13/2019  T. Yamada   Original

    use mod_types
	use mod_calibration
	use mod_variables

	implicit none

	!****** list of local variable ******
	real(dp), dimension(na, ne, nk, nz) :: coef
	!------------------------------------

	contains

	subroutine end_grid_method(icept, slope, policy)
	!  Purpose:
	!	Using the endogenous grid point method,
	!	solve households' optimization problem.
	!
	!  Record of revisions:
	!     Date     Programmer  Description of change
	!  ==========  ==========  =====================
	!  09/13/2019  T. Yamada   Original

		use mod_interpolation
		use mod_useful_functions

		!****** input *******
	    real(dp), intent(in), dimension(nz) :: icept, slope

		!****** output ******
		real(dp), intent(out), dimension(nd, ne, nk, nz) :: policy

        !****** citeria for convergence ******
        integer :: it
        real(dp) :: metric
		integer,  parameter :: maxit = 1000
		real(dp), parameter :: toler = 1.0e-004

		!****** local variables *******
		integer :: i, e, k, z
		real(dp) :: rent, wage, rhs
		real(dp), dimension(nk, nz) :: agg_cap_next
		real(dp), dimension(0:na, ne, nk, nz) :: con0, con1, coh0, coh1
		real(dp), dimension(na, ne, nk, nz) :: asset

        real(dp), dimension(nk, nz) :: tmp1, tmp2
		!------------------------------

		! initialize metric of error
		metric = 1.0

		! initial guess (hand-to-mouth)
		con0(0, :, :, :) = 0.0
		coh0(0, :, :, :) = 0.0

		do z = 1, nz
			do k = 1, nk
				call factor_price(kgrid(k), L_agg(z), tfp(z), alpha, delta, rent, wage)
                ! debug: wage is slightly higher in the recession because labor is relatively scarce in the recession by assumption
                !tmp1(k, z) = rent
                !tmp2(k, z) = wage
                do e = 1, ne
					coh0(1:na, e, k, z) = wage*endow(e) + (1.0+rent)*aprime
					con0(1:na, e, k, z) = coh0(1:na, e, k, z)
				end do
			end do
        end do

        ! given current capital and approximate aggregation, compute the next period's capital
		do z = 1, nz
			do k = 1, nk
				agg_cap_next(k, z) = approx_agg(kgrid(k), z, nz, icept, slope)
			end do
		end do

		! consumption function over cash on hand
		do it = 1, maxit

			if (metric < toler) exit

			! update cubic spline coefficient from numerical recipes
			do z = 1, nz
				do k = 1, nk
					do e = 1, ne
						call spline(coh0(1:na, e, k, z), con0(1:na, e, k, z), 1.0e30_dp, 1.0e30_dp, coef(:, e, k, z))
					end do
				end do
			end do

			! main part of EGM
			do z = 1, nz ! current TFP level
				do k = 1, nk ! current aggregate capital
					do e = 1, ne ! current employment status
						do i = 1, na ! current saving
							rhs              = RHS_Euler(aprime(i), e, agg_cap_next(k, z), z, con0, coh0)
							con1(i, e, k, z) = rhs**(-1.0/gamma)
							coh1(i, e, k, z) = con1(i, e, k, z) + aprime(i)
						end do
					end do
				end do
			end do
			con1(0, :, :, :) = 0.0
			coh1(0, :, :, :) = 0.0

			! check converegence
			metric = maxval(abs((con0(1:na, :, :, :)-con1(1:na, :, :, :))/con0(1:na, :, :, :)))

			! update policy funciton
			con0 = con1
			coh0 = coh1

        end do

        ! debug: consumption function
        !open(30, file="cash_on_hand_g.txt")
        !open(31, file="cash_on_hand_b.txt")
        !open(32, file="consumption_g.txt")
        !open(33, file="consumption_b.txt")
        !    do e = 1, ne
        !        do i = 0, na
        !            write (30, *) coh0(i, e, nk/2+1, 1)
        !            write (31, *) coh0(i, e, nk/2+1, 2)
        !            write (32, *) con0(i, e, nk/2+1, 1)
        !            write (33, *) con0(i, e, nk/2+1, 2)
        !        end do
        !    end do
        !close (30)
        !close (31)
        !close (32)
        !close (33)

		! retrieve current financial asset
		do z = 1, nz
			do k = 1, nk
				call factor_price(kgrid(k), L_agg(z), tfp(z), alpha, delta, rent, wage)
				do e = 1, ne
					do i = 1, na
						asset(i, e, k, z) = (coh0(i, e, k, z) - wage*endow(e)) / (1.0+rent)
					end do
				end do
			end do
        end do

		do z = 1, nz
			do k = 1, nk
				do e = 1, ne
					call spline(asset(:, e, k, z), aprime, 1.0e30_dp, 1.0e30_dp, coef(:, e, k, z))
					do i = 1, nd
						! out of range exception
						!if (asset(1, e, k, z) > amin .and. grid(i) < asset(i, e, k, z)) then
						!	policy(i, e, k, z) = interp1(asset(:, e, k, z), aprime, grid(i))
						if (asset(1, e, k, z) > grid(i)) then
							policy(i, e, k, z) = amin
                        else if (grid(i) > asset(na, e, k, z)) then
							policy(i, e, k, z) = interp1(asset(:, e, k, z), aprime, grid(i))
						else
							policy(i, e, k, z) = splint(asset(:, e, k, z), aprime, coef(:, e, k, z), grid(i))
							!policy(i, e, k, z) = interp1(asset(:, e, k, z), aprime, grid(i))
						end if
						if (policy(i, e, k, z) < amin) then
							policy(i, e, k, z) = amin
						end if
					end do
				end do
			end do
		end do

        ! debug: policy function
        !open(34, file="grid.txt")
        !open(35, file="policy_g.txt")
        !open(36, file="policy_b.txt")
        !    do i = 1, nd
        !        write (34, *) grid(i)
        !    end do
        !    do e = 1, ne
        !        do i = 1, nd
        !            write (35, *) policy(i, e, nk/2+1, 1)
        !            write (36, *) policy(i, e, nk/2+1, 2)
        !        end do
        !    end do
        !close (34)
        !close (35)
        !close (36)

	end subroutine end_grid_method


	function RHS_Euler(asset, e_state, K_agg, z_state, consf, xgrid)
	!  Purpose:
	!	Compute the right hand side of the Euler equation.
	!
	!  Record of revisions:
	!	  Date     Programmer  Description of change
	!  ==========  ==========  =====================
	!  09/13/2019  T. Yamada   Original

		use mod_variables
		use mod_interpolation
		use mod_useful_functions

        !****** input ******
		real(dp), intent(in) :: asset, K_agg
		integer, intent(in) :: e_state, z_state
		real(dp), intent(in), dimension(0:na, ne, nk,nz) :: consf, xgrid

        !****** output ******
		real(dp) :: RHS_Euler

		!****** local variables ******
		integer :: e, z, kloc
		real(dp) :: rent, wage, coh, exp_value, weight, cons0, cons1
		real(dp), dimension(ne, nz) :: cons, mu
		!-----------------------------

		cons = 0.0
		mu   = 0.0

        do z = 1, nz ! next period's TFP level
			do e = 1, ne ! next period's employment status

                ! factor prices in the next period
				call factor_price(K_agg, L_agg(z), tfp(z), alpha, delta, rent, wage)
				coh = wage*endow(e) + (1.0+rent)*asset

				! use linear interpolation over aggregate capital grid
				kloc = locate(kgrid, K_agg)
				if (kloc >= nk) then
					kloc = nk - 1
				else if (kloc < 1) then
					kloc = 1
				end if
				weight = (K_agg - kgrid(kloc))/(kgrid(kloc+1) - kgrid(kloc))
				if (weight > 1.0) then
					weight = 1.0
				else if (weight < 0.0) then
					weight = 0.0
				end if

				! lower capital grid
				if (coh < xgrid(1, e, kloc, z)) then
					cons0 = coh
				else if (coh > xgrid(na, e, kloc, z)) then
					cons0 = interp1(xgrid(1:na, e, kloc, z), consf(1:na, e, kloc, z), coh)
				else
					cons0 = splint(xgrid(1:na, e, kloc, z), consf(1:na, e, kloc, z), coef(:, e, kloc, z), coh)
				end if

				! upper capital grid
				if (coh < xgrid(1, e, kloc+1, z)) then
					cons1 = coh
				else if (coh > xgrid(na, e, kloc+1, z)) then
					cons1 = interp1(xgrid(1:na, e, kloc+1, z), consf(1:na, e, kloc+1, z), coh)
				else
					cons1 = splint(xgrid(1:na, e, kloc+1, z), consf(1:na, e, kloc+1, z), coef(:, e, kloc+1, z), coh)
				end if

				cons(e, z) = (1.0-weight)*cons0 + weight*cons1
    
                ! marginal utility x gross interest rate
				mu(e, z)   = marg_util(cons(e, z), gamma)*(1.0+rent)

			end do
        end do

        ! expected value
		exp_value = 0.0
		do z = 1, nz
			do e = 1, ne
				exp_value = prob(e_state, z_state, e, z)*mu(e,z) + exp_value
			end do
		end do
		RHS_Euler = beta*exp_value

	end function RHS_Euler

end module mod_EGM
