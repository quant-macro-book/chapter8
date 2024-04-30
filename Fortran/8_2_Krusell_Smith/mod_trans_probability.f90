module mod_trans_probability
!  Purpose:
!   Collect subroutines for computing the transition probability matrix.
!
!  Record of revisions:
!     Date     Programmer  Description of change
!  ==========  ==========  ======================
!  09/13/2019  T. Yamada   Original

    use mod_types

    implicit none

    contains

    subroutine trans_prob(prob, tran_zz)
    !  Purpose:
    !	Transition probability of idiosyncratic and aggregate shocks.
    !
    !  Input:
    !	nothing
    !
    !  Output:
    !   prob:     transition probability of (e,z)
    !   tranz_ss: transition probability of z
    !
    !	pi_dist  :stationary distribution
    !  Record of revisions:
    !	  Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  09/13/2019  T. Yamada   Original

	    use mod_calibration

	    !****** output ******
	    real(dp), dimension(ne, nz, ne, nz) :: prob
	    real(dp), dimension(nz, nz) :: tran_zz

	    !****** list of local variables ******
        integer :: i, j, k, l
	    real(dp), dimension(ne) :: pi_e
	    real(dp), dimension(nz) :: pi_z
	    real(dp), dimension(ne, ne) :: tran_gg, tran_bb, tran_gb, tran_bg
	    real(dp) :: a12, a22, b12, b22
        real(dp) :: pgg00, pbb00, pbg00, pgb00, pgg01, pbb01, pbg01, pgb01, pgg, pgb
        real(dp) :: pgg10, pbb10, pbg10, pgb10, pgg11, pbb11, pbg11, pgb11, pbg, pbb
        real(dp) :: test11, test12, test21, test22
	    real(dp), dimension(ne, nz, ne, nz) :: tmp ! pr(e,z,e',z')
	    !-------------------------------------

	    prob = 0.0

	    ! transition probability of (z,z'): business cycle
	    call gen_prob_matrix(durgd, 0.5d0, nz, pi_z, tran_zz)

	    ! transition probability of (e,e') when boom
	    !call gen_prob_matrix(durug, unempg, ne, pi_e, tran_gg)
	    !tmp(:, 1, :, 1) = tran_zz(1, 1)*tran_gg

	    ! from bad to good (with additional assumption: 0.75)
	    !a22 = 0.75*tran_gg(2,2) ! pgb00
	    !a12 = (unempg - unempb*a22) / (1.0-unempb) ! pgb01
	    !tran_gb(1, 1) = 1.0 - a12
	    !tran_gb(1, 2) = a12
	    !tran_gb(2, 1) = 1.0 - a22
	    !tran_gb(2, 2) = a22
	    !tmp(:, 2, :, 1) = tran_zz(1, 2)*tran_gb

	    ! transition probability of (e,e') when recession
	    !call gen_prob_matrix(durub, unempb, ne, pi_e, tran_bb)
	    !tmp(:,2,:,2) = tran_zz(2, 2)*tran_bb

	    ! from good to bad (with additional assumption: 1.25)
	    !b22 = 1.25*tran_bb(2,2) ! pbg00
	    !b12 = (unempb - unempg*b22) / (1.0-unempg) ! pbg01
	    !tran_bg(1, 1) = 1.0 - b12
	    !tran_bg(1, 2) = b12
	    !tran_bg(2, 1) = 1.0 - b22
	    !tran_bg(2, 2) = b22
	    !tmp(:, 1, :, 2) = tran_zz(2, 1)*tran_bg

        ! copy-n-paste from A.Smith's code
        pgg00 = (durug-1.0)/durug ! tran_ee(2,2)|good
        pbb00 = (durub-1.0)/durub ! tran_ee(2,2)|bad
        pbg00 = 1.25*pbb00
        pgb00 = 0.75*pgg00
        pgg01 = (unempg - unempg*pgg00)/(1.0-unempg)
        pbb01 = (unempb - unempb*pbb00)/(1.0-unempb)
        pbg01 = (unempb - unempg*pbg00)/(1.0-unempg)
        pgb01 = (unempg - unempb*pgb00)/(1.0-unempb)
        pgg = (durgd-1.0)/durgd       ! tran_zz(1,1)
        pgb = 1.0 - (durbd-1.0)/durbd ! tran_zz(1,2)

        pgg10 = 1.0 - (durug-1.0)/durug ! tran_ee(1,2)|good
        pbb10 = 1.0 - (durub-1.0)/durub ! tran_ee(1,2)|bad
        pbg10 = 1.0 - 1.25*pbb00
        pgb10 = 1.0 - 0.75*pgg00
        pgg11 = 1.0 - (unempg - unempg*pgg00)/(1.0-unempg)
        pbb11 = 1.0 - (unempb - unempb*pbb00)/(1.0-unempb)
        pbg11 = 1.0 - (unempb - unempg*pbg00)/(1.0-unempg)
        pgb11 = 1.0 - (unempg - unempb*pgb00)/(1.0-unempb)
        pbg = 1.0 - (durgd-1.0)/durgd ! tran_zz(2,1)
        pbb = (durbd-1.0)/durbd       ! tran_zz(2,2)

        prob(1, 1, 1, 1) = pgg*pgg11
        prob(1, 1, 1, 2) = pbg*pbg11
        prob(1, 1, 2, 1) = pgg*pgg01
        prob(1, 1, 2, 2) = pbg*pbg01
        prob(1, 2, 1, 1) = pgb*pgb11
        prob(1, 2, 1, 2) = pbb*pbb11
        prob(1, 2, 2, 1) = pgb*pgb01
        prob(1, 2, 2, 2) = pbb*pbb01
        prob(2, 1, 1, 1) = pgg*pgg10
        prob(2, 1, 1, 2) = pbg*pbg10
        prob(2, 1, 2, 1) = pgg*pgg00
        prob(2, 1, 2, 2) = pbg*pbg00
        prob(2, 2, 1, 1) = pgb*pgb10
        prob(2, 2, 1, 2) = pbb*pbb10
        prob(2, 2, 2, 1) = pgb*pgb00
        prob(2, 2, 2, 2) = pbb*pbb00

        ! debug: check consistency (==1)
        test11 = prob(1, 1, 1, 1)/tran_zz(1, 1) + prob(1, 1, 2, 1)/tran_zz(1, 1)
        test12 = prob(1, 1, 1, 2)/tran_zz(1, 2) + prob(1, 1, 2, 2)/tran_zz(1, 2)
        test21 = prob(1, 2, 1, 1)/tran_zz(2, 1) + prob(1, 2, 2, 1)/tran_zz(2, 1)
        test22 = prob(1, 2, 1, 2)/tran_zz(2, 2) + prob(1, 2, 2, 2)/tran_zz(2, 2)

        ! output transition matrix
        !open(20, file="prob.txt")
        !open(21, file="tran_zz.txt")
        !    do i = 1, nz
        !        do j = 1, ne
        !            do k = 1, nz
        !                do l = 1,ne
        !                    write(20, *) prob(l, k, j, i)
        !                end do
        !            end do
        !        end do
        !    end do
        !    do i = 1, nz
        !        do j = 1,nz
        !            write(21, *) tran_zz(i, j)
        !        end do
        !    end do
        !close (20)
        !close (21)

    end subroutine trans_prob


    subroutine gen_prob_matrix(duration, unemp, ne, pi_dist, prob)
    !  Purpose:
    !	2-state markov chain with fixed duration.
    !
    !  Input:
    !	duration :duration of unemployment/recession
    !	unemploy :unemployment rate/% of recession period
    !   ne       :#state (must be 2)
    !
    !  Output:
    !	pi_dist  :stationary distribution
    !	prob     :transition probability matrix
    !
    !  Record of revisions:
    !	  Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  09/13/2019  T. Yamada   Original

        !****** input ******
	    real(dp), intent(in) :: duration
	    real(dp), intent(in) :: unemp
        integer,  intent(in) :: ne

        !****** output ******
	    real(dp), intent(out), dimension(ne) :: pi_dist
	    real(dp), intent(out), dimension(ne, ne) :: prob
	    !--------------------

	    pi_dist(1) = 1.0 - unemp
	    pi_dist(2) = unemp

	    prob(2,2) = (duration-1.0) / duration
	    prob(2,1) = 1.0 - prob(2, 2)
	    prob(1,1) = 1.0 - (1.0-prob(2, 2))*(pi_dist(2)/pi_dist(1))
	    prob(1,2) = 1.0 - prob(1, 1)

    end subroutine gen_prob_matrix


    subroutine test_prob_zzee
    !  Purpose:
    !	Test the transition probability matrix.
    !   "Not used" in the main code.
    !  
    !  Input:
    !	nothing
    !
    !  Output:
    !   "z_path.txt" & "unemp_rate.txt"
    !
    !  Record of revisions:
    !	  Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  10/22/2019  T. Yamada   Original

        use mod_calibration
        use mod_variables

        !****** local variables ******
        integer :: t, i, zz, ee, zt, count
        real(dp) :: pr_next
        integer, parameter :: np = 11000 ! simulation period
        integer, parameter :: nw = 10000  ! #workers
        integer, dimension(np) :: z_path ! path of aggregate state
        integer, dimension(nw) :: emp0 ! emp/unemp distribution in the current period
        integer, dimension(nw) :: emp1 ! emp/unemp distribution in the next period
        real(dp), dimension(np) :: unemp_rate ! unemployment rate
	    real(dp), dimension(np) :: rndz
	    real(dp), dimension(nw) :: rndw
        !-----------------------------

        ! all individuals start with "workers"
        emp0 = 1
        unemp_rate(1) = 0.0

        ! economy start with "boom"
        z_path(1) = 1
        zt = z_path(1)

    	! set seed of randum number: use Fortran's default random number generator
	    call random_seed(put=(/1022/))
	    call random_number(rndz)

        ! business cycle dynamics
        do t = 1, np
		    pr_next = tran_zz(zt, 1)
		    do zz = 1, nz
			    if (rndz(t) <= pr_next) then
				    z_path(t) = zz
				    if (z_path(t) /= 0) exit
			    else
				    pr_next = pr_next + tran_zz(zt, zz+1)
			    end if
		    end do
		    zt = z_path(t)
        end do

        ! business cycle frequency
        count = 0
        do t = 1, np
            if (z_path(t) == 1) then
                count = count + 1
            end if
        end do

		write (*,"(' boom period:             ', i7)") count
		write (*,"(' recession period:        ', i7)") np - count

        ! unemployment rate over business cycle
        do t = 1, np-1

            ! generate randum variable for each period
	        call random_number(rndw)

            do i = 1, nw
 
                ! worker
                if (emp0(i) == 1) then
                    pr_next = prob(1, z_path(t), 1, z_path(t+1)) / tran_zz(z_path(t), z_path(t+1))
                    do ee = 1, ne
			            if (rndw(i) <= pr_next) then
				            emp1(i) = ee
				            if (emp1(i) /= 0) exit
			            else
				            pr_next = pr_next + prob(1, z_path(t), ee+1, z_path(t+1)) / tran_zz(z_path(t), z_path(t+1))
			            end if
                    end do
                end if

                ! unemployed
                if (emp0(i) == 2) then
                    pr_next = prob(2, z_path(t), 1, z_path(t+1)) / tran_zz(z_path(t), z_path(t+1))
                    do ee = 1, ne
			            if (rndw(i) <= pr_next) then
				            emp1(i) = ee
				            if (emp1(i) /= 0) exit
			            else
				            pr_next = pr_next + prob(2, z_path(t), ee+1, z_path(t+1)) / tran_zz(z_path(t), z_path(t+1))
			            end if
                    end do
                end if

            end do

            ! count unemployed workers
            count = 0
            do i = 1, nw
                if (emp1(i) == 2) then
                    count = count + 1 
                end if
            end do

            ! unemployemt rate
            unemp_rate(t+1) = dble(count)/dble(nw)

            ! update employment distribution
            emp0 = emp1
            emp1 = 0

        end do

        ! output simulation results
        open(20, file="z_path.txt")
        open(21, file="unemp_rate.txt")
            do t = 1, np
                write(20, *) z_path(t)
                write(21, *) unemp_rate(t)
            end do
        close (20)
        close (21)

    end subroutine test_prob_zzee

end module mod_trans_probability