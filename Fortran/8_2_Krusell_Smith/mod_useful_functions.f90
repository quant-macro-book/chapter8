module mod_useful_functions
!  Purpose:
!	Collects several functions frequently used in economics.
!
!  Record of revisions:
!     Date     Programmer  Description of change
!  ==========  ==========  =====================
!  09/13/2019  T. Yamada   Original

    use mod_types

    implicit none	

    contains

	function CRRA(cons, risk_aversion)
	!  Purpose:
	!	CRRA type utility function
	!
	!  Record of revisions:
	!     Date     Programmer  Description of change
	!  ==========  ==========  =====================
	!  05/12/2005  T. Yamada   Original
		!****** list of local variables ******
		real(dp), intent(in) :: cons, risk_aversion
		real(dp) :: CRRA, first, second, value
		real(dp) :: epsilon = 1.0e-008
		real(dp) :: ra_lag  = 1.0e-006
		!-------------------------------------
		if (cons > epsilon) then
			if (risk_aversion-1.0 < ra_lag) then
				CRRA = log(cons)
			else
				CRRA = (1.0/(1.0-risk_aversion)) * (cons**(1.0-risk_aversion))
			endif
		else
			if (risk_aversion == 1.0) then
				value  = log(epsilon)
				first  =  1.0 / epsilon
				second = -1.0 * (epsilon**-2)
				CRRA   = value + first*(cons-epsilon) + (second*(cons-epsilon)**2) / 2.0
			else
				value  = (1.0/(1.0-risk_aversion)) * (epsilon**(1.0-risk_aversion))
				first  = epsilon**(-risk_aversion)
				second = -risk_aversion * epsilon**(-risk_aversion-1.0)
				CRRA   = value + first*(cons-epsilon) + (second*(epsilon-epsilon)**2) / 2.0
			end if
		end if
	end function CRRA


	function marg_util(cons, risk_aversion)
	!  Purpose:
	!	marginal utility of CRRA type utility function.
	!	稲田条件があるため、極端に0に近い場合や負値では
	!	うまく定義できない場合がある.
	!	それを避けるために、消費量がepsilon以下の場合には、
	!	2次関数で近似する.
	!
	!  Record of revisions:
	!    Date      Programmer  Description of change
	!  ==========  ==========  =====================
	!  05/12/2005  T. Yamada   Original
		!****** list of local variables ******
		real(dp),intent(in) :: cons, risk_aversion
		real(dp) :: marg_util
		real(dp) :: value, first, second
		real(dp) :: epsilon = 1.0e-008
		!-----------------------------
		if (cons > epsilon) then
			if (risk_aversion == 1.0) then
				marg_util = 1.0/cons
			else
				marg_util = cons**(-risk_aversion)
			endif
		else
			if (risk_aversion == 1.0) then
				value  = 1.0/epsilon
				first  = -1.0 * (epsilon**-2)
				second =  2.0 * (epsilon**-3)
				marg_util = value + first*(cons-epsilon) + (second*(cons-epsilon)**2) / 2.0
			else
				value  = epsilon**(-risk_aversion)
				first  = -risk_aversion * epsilon**(-risk_aversion-1.0)
				second = (-risk_aversion*(-risk_aversion-1.0)) * epsilon**(-risk_aversion-2.0)
				marg_util = value + first*(cons-epsilon) + (second*(cons-epsilon)**2) / 2.0
			end if
		end if
	end function marg_util


	subroutine factor_price(agg_cap, agg_lab, tfp, alpha, delta, rent, wage)
	!  Purpose
	!	Assuming the production function is of Cobb=Douglas type,
	!	compute factor prices from the first order condition.
	!
	!  Input:
	!	agg_cap :aggregate capital
	!	agg_lab :aggregate labor supply
    !   tfp     :TFP level
	!	alpha   :capital share
	!	delta   :capital depreciation rate
	!
	!  Output:
	!	rent    :interest rate
	!	wage    :wage
	!
	!  Record of revisions:
	!     Date     Programmer  Description of change
	!  ==========  ==========  =====================
	!  06/08/2007  T. Yamada   Original
		!****** input ******
		real(dp), intent(in)  :: agg_cap, agg_lab, tfp, alpha, delta
		!****** output ******
		real(dp), intent(out) :: rent, wage
		!--------------------
		rent = tfp*     alpha *agg_cap**(alpha-1.0)*agg_lab**(1.0-alpha) - delta
		wage = tfp*(1.0-alpha)*agg_cap** alpha     *agg_lab**(   -alpha)
	end subroutine factor_price


	function approx_agg(agg_cap, state, nz, icept, slope)
	!  Purpose
	!	From approximate aggregation,
    !   compute the next period's aggregate capital.
	!
	!  Input:
	!	agg_cap    :aggregate capital
	!	state      :aggregate state
	!	nz         :# of aggregate state
	!	icept      :intercept
	!	slope      :slope of prediction function
	!
	!  Output:
	!	approx_agg :next period's aggregate capital
	!
	!  Record of revisions:
	!     Date     Programmer  Description of change
	!  ==========  ==========  =====================
	!  09/13/2019  T. Yamada   Original
		!****** input ******
		integer,  intent(in) :: nz, state
		real(dp), intent(in) :: agg_cap
		real(dp), intent(in), dimension(nz) :: icept, slope
		!****** output ******
		real(dp) :: approx_agg
		!--------------------
		approx_agg = exp(icept(state) + slope(state)*log(agg_cap))
	end function approx_agg

end module mod_useful_functions
