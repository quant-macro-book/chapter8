 module mod_law_of_motion
!  Purpose:
!	Compute invariant distribution using Young(2010)'s method.
!
!  Record of revisions:
!	  Date     Programmer  Description of change
!  ==========  ==========  =====================
!  09/13/2019  T. Yamada   Original

    use mod_types
	use mod_calibration

	implicit none

    contains

    subroutine law_of_motion(policy, k_path, z_path)
    !  Purpose:
    !	given policy function, compute wealth distribution.
    !  
    !  Record of revisions:
    !	  Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  09/26/2019  T. Yamada   Original

	    use mod_variables
	    use mod_interpolation

	    !****** input *******
	    real(dp), intent(in), dimension(nd, ne, nk, nz) :: policy

	    !****** output *******
	    real(dp), intent(out), dimension(nums) :: k_path
	    integer, intent(out), dimension(nums) :: z_path

	    !****** local variables ******
	    integer :: it, ie, ie2, ik, i, e, k, z, status
	    integer,  dimension(nd, ne, nk, nz) :: loca
	    real(dp), dimension(nd, ne, nk, nz) :: wght
	    real(dp), dimension(ne, ne) :: tran_ee
	    integer :: aloc, kloc
	    real(dp) :: bunbo, bunsi, wta, wtk
	    real(dp), dimension(nd, ne) :: psi0, psi1
	    real(dp), dimension(nd) :: dist
	    !---------------------------------------------

	    k_path = 0.0

	    ! generate a sequence of aggregate shock
	    call markov(tran_zz, nz, nums, z_path)

	    ! setup for distribution iteration: for details, see Young (2010)
	    do z = 1, nz
		    do k = 1, nk
			    do e = 1, ne
				    do i = 1, nd
					    aloc = locate(grid, policy(i, e, k, z))
					    if (aloc >= nd) then
						    aloc = nd - 1
					    else if (aloc <= 1) then
						    aloc = 1
					    end if
					    loca(i, e, k, z) = aloc
					    bunsi            = policy(i, e, k, z) - grid(aloc)
					    bunbo            = grid(aloc+1) - grid(aloc)
					    wght(i, e, k, z) = bunsi / bunbo
					    if (wght(i, e, k, z) > 1.0) then
						    wght(i, e, k, z) = 1.0
					    else if (bunsi < 0.0) then
						    wght(i, e, k, z) = 0.0
					    end if
				    end do
			    end do
		    end do
	    end do

	    ! initial distribution(uniform)
	    open (10, file="psi.ini")
		    do e = 1, ne
			    do i = 1, nd
				    read(10, *) psi0(i, e)
			    end do
		    end do
	    close (10)

	    dist = 0.0
	    do e = 1, ne
		    dist = psi0(:, e) + dist
        end do
        k_path(1) = dot_product(grid, dist)

	    ! iterate distribution forwardly 
	    do it = 1, nums-1

		    ! use linear interpolation over aggregate capital grid
		    kloc = locate(kgrid, k_path(it))
		    if (kloc >= nk) then
			    kloc = nk - 1
		    else if (kloc < 1) then
			    kloc = 1
		    end if
		    wtk = (k_path(it) - kgrid(kloc)) / (kgrid(kloc+1) - kgrid(kloc))
		    if (wtk > 1.0) then
			    wtk = 1.0
		    else if (wtk < 0.0) then
			    wtk = 0.0
		    end if

		    do ie = 1, ne
			    do ik = 1, nd
				    do ie2 = 1, ne
					    ! lower bound
					    aloc    = loca(ik, ie, kloc, z_path(it))
					    wta     = wght(ik, ie, kloc, z_path(it))
					    tran_ee = prob(:, z_path(it), :, z_path(it+1)) / tran_zz(z_path(it), z_path(it+1))
					    psi1(aloc,   ie2) = (1.0-wtk)*(1.0-wta)*tran_ee(ie, ie2)*psi0(ik, ie) + psi1(aloc,   ie2)
					    psi1(aloc+1, ie2) = (1.0-wtk)*     wta *tran_ee(ie, ie2)*psi0(ik, ie) + psi1(aloc+1, ie2)

					    ! upper bound
					    aloc    = loca(ik, ie, kloc+1, z_path(it))
					    wta     = wght(ik, ie, kloc+1, z_path(it))
					    tran_ee = prob(:, z_path(it), :, z_path(it+1)) / tran_zz(z_path(it), z_path(it+1))
					    psi1(aloc,   ie2) = wtk*(1.0-wta)*tran_ee(ie, ie2)*psi0(ik, ie) + psi1(aloc,   ie2)
					    psi1(aloc+1, ie2) = wtk*     wta *tran_ee(ie, ie2)*psi0(ik, ie) + psi1(aloc+1, ie2)
				    end do
			    end do
		    end do

		    ! integerate distribution
		    dist = 0.0
		    do e = 1, ne
			    dist = psi1(:, e) + dist
            end do
            k_path(it+1) = dot_product(grid, dist)

		    ! update
		    psi0 = psi1
		    psi1 = 0.0

	    end do

        ! output simulation results
        !open(20, file="z_path.txt")
        !open(21, file="k_path.txt")
        !    do it = 1, nums
        !        write(20, *) z_path(it)
        !        write(21, *) k_path(it)
        !    end do
        !close (20)
        !close (21)

    end subroutine law_of_motion


    subroutine law_of_motion_sim(policy, k_path, z_path)
    !  Purpose:
    !	given policy function, compute wealth distribution.
    !  
    !  Record of revisions:
    !	  Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  09/26/2019  T. Yamada   Original

	    use mod_variables
	    use mod_interpolation

	    !****** input *******
	    real(dp), intent(in), dimension(nd, ne, nk, nz) :: policy

	    !****** output *******
	    real(dp), intent(out), dimension(nums) :: k_path
	    integer, intent(out), dimension(nums) :: z_path

	    !****** local variables ******
	    integer :: it, i, e, kloc, count
        real(dp) :: pr_next, savings0, savings1, wght
        integer, dimension(numi) :: emp0, emp1 ! employment status
	    real(dp), dimension(numi) :: rndw
	    real(dp), dimension(numi) :: cap0, cap1 ! savings
        real(dp), dimension(nums) :: u_path ! unemployment rate
	    !---------------------------------------------

        ! set random seed
	    call random_seed(put=(/1024/))

	    k_path = 0.0

	    ! generate a sequence of aggregate shock
	    call markov(tran_zz, nz, nums, z_path)

	    ! initial distribution: all individuals start as workers
        cap0(:) = 35.0 ! savings
        emp0(1:int(0.9*numi)) = 1 ! employed
        emp0(int(0.9*numi)+1:numi) = 2 ! unemployed
        k_path(1) = sum(cap0)/dble(numi)
        u_path(1) = 0.1

	    ! iterate distribution forwardly 
	    do it = 1, nums-1

            ! generate randum variable for each period
	        call random_number(rndw)
 
            do i = 1, numi

                ! worker         
                if (emp0(i) == 1) then

                    ! savings in the next period
		            kloc = locate(kgrid, k_path(it))
				    if (kloc >= nk) then
					    kloc = nk - 1
				    else if (kloc <= 1) then
					    kloc = 1
				    end if
                    savings0 = interp1(grid, policy(:, 1, kloc, z_path(it)), cap0(i))
                    savings1 = interp1(grid, policy(:, 1, kloc+1, z_path(it)), cap0(i))
				    wght = (k_path(it) - kgrid(kloc)) / (kgrid(kloc+1) - kgrid(kloc))
					if (wght > 1.0) then
						wght = 1.0
					else if (kgrid(kloc+1) - kgrid(kloc) < 0.0) then
						wght = 0.0
					end if
                    cap1(i) = (1.0-wght)*savings0 + wght*savings1

                    ! transition of employment status
                    pr_next = prob(1, z_path(it), 1, z_path(it+1)) / tran_zz(z_path(it), z_path(it+1))
                    do e = 1, ne
			            if (rndw(i) <= pr_next) then
				            emp1(i) = e
				            if (emp1(i) /= 0) exit
			            else
				            pr_next = pr_next + prob(1, z_path(it), e+1, z_path(it+1)) / tran_zz(z_path(it), z_path(it+1))
			            end if
                    end do

                end if

                ! unemployed
                if (emp0(i) == 2) then

                    ! savings in the next period
		            kloc = locate(kgrid, k_path(it))
				    if (kloc >= nk) then
					    kloc = nk - 1
				    else if (kloc <= 1) then
					    kloc = 1
				    end if
                    savings0 = interp1(grid, policy(:, 2, kloc, z_path(it)), cap0(i))
                    savings1 = interp1(grid, policy(:, 2, kloc+1, z_path(it)), cap0(i))
				    wght = (k_path(it) - kgrid(kloc)) / (kgrid(kloc+1) - kgrid(kloc))
					if (wght > 1.0) then
						wght = 1.0
					else if (kgrid(kloc+1) - kgrid(kloc) < 0.0) then
						wght = 0.0
					end if
                    cap1(i) = (1.0-wght)*savings0 + wght*savings1
                    
                    ! transition of employment status
                    pr_next = prob(2, z_path(it), 1, z_path(it+1)) / tran_zz(z_path(it), z_path(it+1))
                    do e = 1, ne
			            if (rndw(i) <= pr_next) then
				            emp1(i) = e
				            if (emp1(i) /= 0) exit
			            else
				            pr_next = pr_next + prob(2, z_path(it), e+1, z_path(it+1)) / tran_zz(z_path(it), z_path(it+1))
			            end if
                    end do

                end if

            end do

		    ! aggregate capital and unemployment rate
            count = 0
            do i = 1, numi
                if (emp1(i) == 2) then
                    count = count + 1 
                end if
            end do
            u_path(it+1) = dble(count)/dble(numi)
            k_path(it+1) = dble(sum(cap1))/dble(numi)

		    ! update
            emp0 = emp1
		    cap0 = cap1
            emp1 = 0
            cap1 = 0.0

	    end do

        ! output simulation results
        open(20, file="z_path.txt")
        open(21, file="k_path.txt")
        open(22, file="u_path.txt")
            do it = 1, nums
                write(20, *) z_path(it)
                write(21, *) k_path(it)
                write(22, *) u_path(it)
            end do
        close (20)
        close (21)
        close (22)

    end subroutine law_of_motion_sim


    subroutine markov(tran, ng, ns, seq_state)
    !  Purpose:
    !	Given ng*ng transition matrix and state of the grid,
    !	generate simulated sequence of ng period.
    !
    !  Input:
    !	tran     : transition matrix of (ng*ng)
    !	ng       : dimension of the matrix
    !	ns       : simulation period
    !	node     : markov state's value
    !
    !  Output:
    !	seq_state: sequence of state(node index)
    !
    !	Required external code: nrtype
    !
    !  Record of revisions:
    !	   Date     Programmer  Description of change
    !	==========  ==========  =====================
    !	04/05/2007  T. Yamada   Original

	    !****** input ******
	    integer, intent(in) :: ng
	    integer, intent(in) :: ns
	    real(dp), intent(in), dimension(ng, ng) :: tran

	    !****** output ******
	    integer, intent(out), dimension(ns) :: seq_state

	    !****** list of local variables ******
	    integer :: i, j, current
	    real(dp) :: next_pr
	    real(dp), dimension(ns) :: X
	    !-------------------------------------

    	!call random_seed( )
	    call random_seed(put=(/225/))
	    call random_number(X)

	    ! set initial state
	    seq_state(1) = 1
	    current      = seq_state(1)

	    ! main simulation part
	    do j = 2,ns
		    next_pr = tran(current, 1)
		    do i = 1,ng
			    if (X(j) <= next_pr) then
				    seq_state(j) = i
				    if (seq_state(j) /= 0) exit
			    else
				    next_pr = next_pr + tran(current,i+1)
			    end if
		    end do
		    current = seq_state(j)
        end do

    end subroutine markov


    subroutine regress(k_path, z_path, icept, slope, R2)
    ! Purpose:
    !  Given k and z path, compute regression coefficients
    !
    !  log(K(t+1)) = a(z) + b(z)*log(K(t))
    !
    ! Record of revisions:
    !    Date     Programmer  Description of change
    ! ==========  ==========  =====================
    ! 09/26/2019  T. Yamada   Original

	    use mod_calibration
	    use mod_variables, only:nums
	    use RONE_INT

	    !****** input ******
	    real(dp), intent(in), dimension(nums) :: k_path
	    integer, intent(in), dimension(nums) :: z_path

	    !****** output ******
	    real(dp), intent(out), dimension(nz) :: icept
	    real(dp), intent(out), dimension(nz) :: slope
	    real(dp), intent(out), dimension(nz) :: R2

	    !****** list of local variables ******
	    integer :: it, z, status
	    integer, dimension(nz) :: nobs, counter
	    real(dp), allocatable, dimension(:,:) :: kk_path, casevar
	    real(dp), dimension(15) :: AOV
	    real(dp), dimension(2, 5) :: COEF
	    real(dp), dimension(2, 2) :: COVB
	    real(dp), dimension(10) :: TESTLF
	    !--------------------------------------------------

	    ! count realized states(discard first 1000 periods)
	    nobs = 0
	    do it = 1001, nums-1
		    do z = 1, nz
			    if (z_path(it) == z) then
				    nobs(z) = nobs(z) + 1
			    end if
		    end do
	    end do

	    ! construct data
	    allocate(kk_path(1:nobs(1), 2), STAT=status)
	    allocate(casevar(1:nobs(1), 12), STAT=status)
	    kk_path = 0.0
	    casevar = 0.0

	    ! find regression parameters
	    counter = 1
	    do z = 1, nz

		    if (z /= 1) then
			    deallocate(kk_path, STAT=status)
			    deallocate(casevar, STAT=status)
			    allocate(kk_path(1:nobs(z), 2), STAT=status)
			    allocate(casevar(1:nobs(z), 12), STAT=status)
			    kk_path = 0.0
			    casevar = 0.0
		    end if

		    do it = 1001, nums-1
			    if (z_path(it) == z) then
				    kk_path(counter(z), 1) = log(k_path(it))
				    kk_path(counter(z), 2) = log(k_path(it+1))
				    counter(z)             = counter(z) + 1
			    end if
		    end do

		    ! regression: use IMSL
		    call RONE(kk_path, 2, 1, AOV, COEF, COVB, TESTLF, casevar)
		    icept(z) = COEF(1, 1)
		    slope(z) = COEF(2, 1)
		    R2(z)    = AOV(11) / 100.0

	    end do

	    deallocate(kk_path, STAT=status)
	    deallocate(casevar, STAT=status)

    end subroutine regress

end module mod_law_of_motion
