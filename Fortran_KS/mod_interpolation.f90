MODULE mod_interpolation
!
!  Purpose:
!    Linear interpolation: interp1
!    Cubic spline interpolation: splint, spline
!    Polynomial interpolation: polint
!    Use subroutines in Numerical Recipes
!
!  Record of revisions:
!     Date     Programmer  Description of change
!  ==========  ==========  ======================
!  02/13/2005  T. Yamada   From Numerical Recipes
!  11/15/2012  T. Yamada   Collect subroutines
!  02/12/2016  T. Yamada   Add polint (from nr)

    IMPLICIT NONE

    !****** set double precision ******
    INTEGER, PARAMETER :: I4B  = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: prec = SELECTED_REAL_KIND(p = 15, r = 307)
    !----------------------------------

    !****** interface from "nrutil.f90" ******
    INTERFACE assert_eq
        MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
    END INTERFACE
    !-----------------------------------------

CONTAINS

    FUNCTION interp1(x_axis, y_axis, x) RESULT(y)
        !****** input ******
        REAL(prec), INTENT(IN), DIMENSION(:) :: x_axis
        REAL(prec), INTENT(IN), DIMENSION(:) :: y_axis
        REAL(prec), INTENT(IN) :: x
        !****** output ******
        REAL(prec) :: y
        !****** local variable ******
        INTEGER :: j
        REAL(prec) :: A, B
        !----------------------------
        j = locate(x_axis, x)
        IF (j == 0) j = 1
        IF (j == SIZE(x_axis)) j = SIZE(x_axis) - 1
        IF ((x_axis(j+1)-x_axis(j)) == 0.0) THEN
            A = 0.5
        ELSE
            A = (x_axis(j+1)-x) / (x_axis(j+1)-x_axis(j))
        END IF
        IF ((x_axis(j+1)-x_axis(j)) == 0.0) THEN
            B = 0.5
        ELSE
            B = (x-x_axis(j)) / (x_axis(j+1)-x_axis(j))
        END IF
        y = A*y_axis(j) + B*y_axis(j+1)
    END FUNCTION interp1

    ! from Numerical Recipes
    FUNCTION locate(xx,x)
        !****** input ******
        REAL(prec), DIMENSION(:), INTENT(IN) :: xx
        REAL(prec), INTENT(IN) :: x
        INTEGER(I4B) :: locate
        INTEGER(I4B) :: n,jl,jm,ju
        LOGICAL :: ascnd
        !-------------------
        n=SIZE(xx)
        ascnd = (xx(n) >= xx(1))
        jl=0
        ju=n+1
        DO
            IF (ju-jl <= 1) EXIT
            jm=(ju+jl)/2
            IF (ascnd .eqv. (x >= xx(jm))) THEN
                jl=jm
            ELSE
                ju=jm
            END IF
        END DO
        IF (x == xx(1)) THEN
            locate=1
        ELSE IF (x == xx(n)) THEN
            locate=n-1
        ELSE
            locate=jl
        END IF
    END FUNCTION locate

    ! from Numerical Recipes
    SUBROUTINE polint(xa,ya,x,y,dy)
        !****** input ******
        REAL(prec), DIMENSION(:), INTENT(IN) :: xa,ya
        REAL(prec), INTENT(IN) :: x
        !****** output ******
        REAL(prec), INTENT(OUT) :: y,dy
        !****** local ******
        INTEGER(I4B) :: m,n,ns
        REAL(prec), DIMENSION(size(xa)) :: c,d,den,ho
        !-------------------
        n=assert_eq(size(xa),size(ya),'polint')
        c=ya
        d=ya
        ho=xa-x
        ns=iminloc(abs(x-xa))
        y=ya(ns)
        ns=ns-1
        do m=1,n-1
            den(1:n-m)=ho(1:n-m)-ho(1+m:n)
            if (any(den(1:n-m) == 0.0)) &
                call nrerror('polint: calculation failure')
            den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
            d(1:n-m)=ho(1+m:n)*den(1:n-m)
            c(1:n-m)=ho(1:n-m)*den(1:n-m)
            if (2*ns < n-m) then
                dy=c(ns+1)
            else
                dy=d(ns)
                ns=ns-1
            end if
            y=y+dy
        end do
    END SUBROUTINE polint

    ! from Numerical Recipes
    FUNCTION splint(xa,ya,y2a,x)
    IMPLICIT NONE
    REAL(prec), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(prec), INTENT(IN) :: x
    REAL(prec) :: splint
    INTEGER(I4B) :: khi,klo,n
    REAL(prec) :: a,b,h
    n=assert_eq(SIZE(xa),SIZE(ya),SIZE(y2a),'splint')
    klo=MAX(MIN(locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    IF (h == 0.0) CALL nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_prec
    END FUNCTION splint

    ! from Numerical Recipes
    SUBROUTINE spline(x,y,yp1,ypn,y2)
    IMPLICIT NONE
    REAL(prec), DIMENSION(:), INTENT(IN) :: x,y
    REAL(prec), INTENT(IN) :: yp1,ypn
    REAL(prec), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER(I4B) :: n
    REAL(prec), DIMENSION(SIZE(x)) :: a,b,c,r
    n=assert_eq(SIZE(x),SIZE(y),SIZE(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0_prec*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_prec*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    IF (yp1 > 0.99e30_prec) THEN
        r(1)=0.0
        c(1)=0.0
    ELSE
        r(1)=(3.0_prec/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        c(1)=0.5
    END IF
    IF (ypn > 0.99e30_prec) THEN
        r(n)=0.0
        a(n)=0.0
    ELSE
        r(n)=(-3.0_prec/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
        a(n)=0.5
    END IF
    CALL tridag_ser(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
    END subroutine spline

    ! from Numerical Recipes
    SUBROUTINE tridag_ser(a,b,c,r,u)
    IMPLICIT NONE
    REAL(prec), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(prec), DIMENSION(:), INTENT(OUT) :: u
    REAL(prec), DIMENSION(SIZE(b)) :: gam
    INTEGER(I4B) :: n,j
    REAL(prec) :: bet
    n=assert_eq((/SIZE(a)+1,SIZE(b),SIZE(c)+1,SIZE(r),SIZE(u)/),'tridag_ser')
    bet=b(1)
    IF (bet == 0.0) CALL nrerror('tridag_ser: Error at code stage 1')
    u(1)=r(1)/bet
    DO j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j-1)*gam(j)
        IF (bet == 0.0) &
            CALL nrerror('tridag_ser: Error at code stage 2')
        u(j)=(r(j)-a(j-1)*u(j-1))/bet
    END DO
    DO j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
    END DO
    END SUBROUTINE tridag_ser

    ! From "nrutil.f90"
    FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    IF (n1 == n2) THEN
        assert_eq2=n1
    ELSE
        WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eq2'
    END IF
    END FUNCTION assert_eq2

    ! From "nrutil.f90"
    FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    IF (n1 == n2 .and. n2 == n3) THEN
        assert_eq3=n1
    ELSE
        WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eq3'
    END IF
    END FUNCTION assert_eq3

    ! From "nrutil.f90"
    FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    IF (n1 == n2 .and. n2 == n3 .and. n3 == n4) THEN
        assert_eq4=n1
    ELSE
        WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eq4'
    END IF
    END FUNCTION assert_eq4

    ! From "nrutil.f90"
    FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    IF (ALL(nn(2:) == nn(1))) THEN
        assert_eqn=nn(1)
    ELSE
        WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eqn'
    END IF
    END FUNCTION assert_eqn

    ! From "nrutil.f90"
    FUNCTION iminloc(arr)
    REAL(prec), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(1) :: imin
    INTEGER(I4B) :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
    END FUNCTION iminloc

    ! From "nrutil.f90"
    SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    WRITE (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
    END SUBROUTINE nrerror

END MODULE mod_interpolation
