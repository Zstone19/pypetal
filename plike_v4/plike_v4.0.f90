!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                  Calculating the Peak Likelihood of the ZDCF
!
!        A special version for processing the ZDCF output file *.dcf
!
!          Version 4.0 (Fortran 95 double precision, 25/1/2013)
!
!    (Bayesian prior assumption of uniform distribution in z-space)
!
!
!                               Tal Alexander
!
!   Reference to the PLIKE algorithm:
!
!   T. Alexander, 2013, 'Improved AGN light curve analysis with the
!   z-transformed discrete correlation function', arXiv:1302.1508,
!   astro-ph.IM
!
!   Reference to the ZDCF algorithm:
!
!   T. Alexander, 1997, 'Is AGN variability correlated with other AGN
!   properties? - ZDCF analysis of small samples of sparse light curves',
!   in "Astronomical Time Series", eds. D. Maoz, A. Sternberg and
!   E.M. Leibowitz, Kluwer, Dordrecht, p. 163
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SAMPLE RUN (user input marked by "<-USER". See explanations below)
! ----------
!
!PLIKE V4.0 begins.
!Enter dcf file name:
!test.dcf   <-USER
!Enter lower bound on peak location:
!-200       <-USER
!Enter upper bound on peak location:
!+200       <-USER
!
!Calculating ML in range t_lag = -2.000E+02 to  2.000E+02
!
!ZDCF peak at  +1.333E-01 r_max =  +8.811E-01 ZDCF C.O.M. at  -4.311E+00
!
!Progress meter: .......
!
!Results:
!
!# PLIKE V4.0 ZDCF peak likelihood analysis
!# for file test.dcf
!#
!# num      lag         r         -dr        +dr    likelihood
!# ---- ---------- ---------- ---------- ---------- ----------
!     1 -7.677E+01  3.552E-01  3.141E-01  2.753E-01  2.701E-03
!     2 -1.979E+01  7.796E-01  1.210E-01  9.680E-02  1.106E-01
!     3  1.333E-01  8.811E-01  8.032E-02  6.053E-02  4.577E-01
!     4  5.618E+00  8.211E-01  1.237E-01  9.322E-02  2.054E-01
!     5  2.410E+01  3.915E-01  3.061E-01  2.648E-01  3.669E-03
!     6  8.748E+01 -1.617E-01  3.189E-01  3.389E-01  2.365E-05
!     7  1.999E+02 -4.352E-01  2.513E-01  2.950E-01  7.652E-07
!#
!# ML Peak at t_lag =  +1.333E-01, Likelihood at peak =  +4.577E-01
!# 1 sigma ML t_lag interval =  +1.333E-01  +1.004E+01  -2.372E+01
!#                           = ( -2.359E+01 ,  +1.017E+01)
!
!Output written on plike.out.
!Program ended.
!
! EXPLANATION
! -----------
!
! o The program uses the *.dcf output file of the ZDCF.
!
! o The user has to give lower and upper bounds on the peak location.
!
! o The output is printed to the terminal and (over-)written on the
!   file plike.out. It consists of the likelihood function and the
!   fiducial interval, given both as t_peak + dt_up - dt_low, and as
!   the interval (t_min, t_max).
!
! FORTRAN ISSUES
! --------------
!
! o To compile:
!   gfortran -O3 -o plike plike_v4.0.f90
!
!   (no-convergence issues with the ifort 13.0 compiler have been
!    reported).
!
! o The program can be VERY slow - don't panic!
!
! o This program incorporates the functions QROMO, POLINT and MIDPNT
!   from Numerical Recipes / Press, Flannery, Teukolsky & Vetterling.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
MODULE cons
  IMPLICIT NONE
  ! Defining double precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 60)
  REAL(dp), PARAMETER :: Pi = 4*ATAN(1.0_dp)
END MODULE cons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE quad77
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NUMERICAL RECIPES "SIMPLE" QUADRATURE ROUTINES CONVERTED FROM FORTRAN77
  ! to REAL(dp) FORTRAN90
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE polint(xa,ya,n,x,y,dy)
    USE cons
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: xa(n),ya(n),x
    INTEGER, INTENT(in) :: n
    REAL(dp), INTENT(out) :: y,dy
    !
    INTEGER, PARAMETER :: NMAX=10
    INTEGER :: i,m,ns
    REAL(dp) :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    !
    ns=1
    dif=ABS(x-xa(1))
    DO i=1,n
       dift=ABS(x-xa(i))
       IF (dift < dif) THEN
          ns=i
          dif=dift
       END IF
       c(i)=ya(i)
       d(i)=ya(i)
    END DO
    y=ya(ns)
    ns=ns-1
    DO m=1,n-1
       DO i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          IF (den == 0.) THEN
             WRITE (*,*) 'failure in polint'
             STOP
          END IF
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       END DO
       IF (2*ns < n-m)THEN
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       END IF
       y=y+dy
    END DO
  END SUBROUTINE polint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE trapzd(func,a,b,s,n)
    USE cons
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    REAL(dp), INTENT(in) :: a, b
    REAL(dp), INTENT(inout) :: s
    INTERFACE
       FUNCTION func(x)
         USE cons
         IMPLICIT NONE
         REAL(dp), INTENT(in) :: x
         REAL(dp) :: func
       END FUNCTION func
    END INTERFACE
    !
    INTEGER :: it,j
    REAL(dp) :: del,sum,tnm,x
    !
    IF (n == 1) THEN
       s=0.5_dp*(b-a)*(func(a)+func(b))
    ELSE
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5_dp*del
       sum=0.0_dp
       DO j=1,it
          sum=sum+func(x)
          x=x+del
       END DO
       s=0.5_dp*(s+(b-a)*sum/tnm)
    END IF
  END SUBROUTINE trapzd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE qromb(func,a,b,ss)
    !
    ! USES polint,trapzd
    !
    USE cons
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: a,b
    REAL(dp), INTENT(out) :: ss

    INTERFACE
       FUNCTION func(x)
         USE cons
         IMPLICIT NONE
         REAL(dp), INTENT(in) :: x
         REAL(dp) :: func
       END FUNCTION func
    END INTERFACE
    REAL(dp), PARAMETER :: EPS=1.e-6
    INTEGER, PARAMETER :: JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1
    INTEGER :: j
    REAL(dp) :: dss,h(JMAXP),s(JMAXP)
    !
    h(1)=1.0_dp
    DO j=1,JMAX
       CALL trapzd(func,a,b,s(j),j)
       IF (j >= K) THEN
          CALL polint(h(j-KM),s(j-KM),K,0.0_dp,ss,dss)
          IF (ABS(dss) <= EPS*ABS(ss)) RETURN
       END IF
       s(j+1)=s(j)
       h(j+1)=0.25_dp*h(j)
    END DO
    WRITE (*,*) 'QROMB: too many steps'
    STOP
  END SUBROUTINE qromb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE qromo(func,a,b,ss,choose)
    USE cons
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: a,b
    REAL(dp), INTENT(out) :: ss
    INTERFACE
       FUNCTION func(x)
         USE cons
         IMPLICIT NONE
         REAL(dp), INTENT(in) :: x
         REAL(dp) :: func
       END FUNCTION func
       SUBROUTINE choose(funk,aa,bb,s,n)
         USE cons
         IMPLICIT NONE
         REAL(dp), INTENT(IN) :: aa,bb
         REAL(dp), INTENT(INOUT) :: s
         INTEGER, INTENT(IN) :: n
         INTERFACE
            FUNCTION funk(x)
              USE cons
              IMPLICIT NONE
              REAL(dp), INTENT(IN) :: x
              REAL(dp) :: funk
            END FUNCTION funk
         END INTERFACE
       END SUBROUTINE choose
    END INTERFACE
    !
    REAL(dp), PARAMETER :: EPS=1.0d-6
    INTEGER, PARAMETER :: JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1
    !dbg: INTEGER, PARAMETER :: JMAX=50, JMAXP=JMAX+1, K=5, KM=K-1
    INTEGER :: j
    REAL(dp) ::  dss,h(JMAXP),s(JMAXP)
    !
    h(1)=1.0_dp
    DO j=1,JMAX
       CALL choose(func,a,b,s(j),j)
       IF (j >= K) THEN
          CALL polint(h(j-KM),s(j-KM),K,0.0_dp,ss,dss)
          IF (ABS(dss) <= EPS*ABS(ss)) RETURN
       END IF
       s(j+1)=s(j)
       h(j+1)=h(j)/9.0_dp
    END DO
    WRITE (*,*) 'QROMO: too many steps'
    STOP
  END SUBROUTINE qromo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE midpnt(func,a,b,s,n)
    USE cons
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: a,b
    REAL(dp), INTENT(INOUT) :: s
    INTEGER, INTENT(IN) :: n
    INTERFACE
       FUNCTION func(x)
         USE cons, ONLY: dp
         IMPLICIT NONE
         REAL(dp), INTENT(IN) :: x
         REAL(dp) :: func
       END FUNCTION func
    END INTERFACE
    !
    INTEGER :: it,j
    REAL(dp) :: ddel,del,sum,tnm,x
    !
    IF (n == 1) THEN
       s=(b-a)*func(0.5_dp*(a+b))
    ELSE
       it=3**(n-2)
       tnm=it
       del=(b-a)/(3.0_dp*tnm)
       ddel=del+del
       x=a+0.5_dp*del
       sum=0.0_dp
       DO j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
       END DO
       s=(s+(b-a)*sum/tnm)/3.0_dp
    END IF
  END SUBROUTINE midpnt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE quad77
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variables for calculating the likelihood
MODULE like_vars
  USE cons
  IMPLICIT NONE
  INTEGER :: np ! Number of points
  ! time and correlation coefficient data
  REAL(dp), ALLOCATABLE :: t(:),dtm(:),dtp(:)
  REAL(dp), ALLOCATABLE :: r(:),drm(:),drp(:)
  ! z-transformed data and number in bin
  REAL(dp), ALLOCATABLE :: z(:)
  INTEGER, ALLOCATABLE :: n(:)
  ! the likelihood vector
  REAL(dp), ALLOCATABLE :: liklhd(:)
  ! Auxiliary integration variables
  INTEGER :: i0  ! The index of the ZDCF point for which the
                 ! likelihood is to be calculated (see LIKINT)
  INTEGER :: n0  ! Bin number (see PHIZ1)
  REAL(dp) :: z0 ! z value (see PHIZ1)
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE allocate_vars(np)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: np
    INTEGER :: ierr
    !
    ALLOCATE (t(np),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate t(',np,') Error ',ierr
       STOP
    END IF
    ALLOCATE (dtm(np),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dtm(',np,') Error ',ierr
       STOP
    END IF
    ALLOCATE (dtp(np),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dtp(',np,') Error ',ierr
       STOP
    END IF
    ALLOCATE (r(np),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate r(',np,') Error ',ierr
       STOP
    END IF
    ALLOCATE (drm(np),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate drm(',np,') Error ',ierr
       STOP
    END IF
    ALLOCATE (drp(np),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate drp(',np,') Error ',ierr
       STOP
    END IF
    ALLOCATE (z(np),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate z(',np,') Error ',ierr
       STOP
    END IF
    ALLOCATE (n(np),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate n(',np,') Error ',ierr
       STOP
    END IF
    ALLOCATE (liklhd(np),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate liklhd(',np,') Error ',ierr
       STOP
    END IF
  END SUBROUTINE allocate_vars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE deallocate_vars()
    IMPLICIT NONE
    DEALLOCATE (t)
    DEALLOCATE (dtm)
    DEALLOCATE (dtp)
    DEALLOCATE (r)
    DEALLOCATE (drm)
    DEALLOCATE (drp)
    DEALLOCATE (z)
    DEALLOCATE (n)
    DEALLOCATE (liklhd)
  END SUBROUTINE deallocate_vars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE like_vars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Functions  for calculating the likelihood
MODULE like_funcs
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate quantities associated with the Fisher transform
  SUBROUTINE calc_Fisher(r,n,zbar,sz,dzdr)
    USE cons
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: r
    INTEGER, INTENT(in) :: n
    REAL(dp), INTENT(out) :: zbar,sz,dzdr
    !
    zbar = LOG((1+r)/(1-r))/2. + r/2/(n-1) &
         & *(1+(5+r**2)/4/(n-1)+(11+2*r**2+3*r**4) &
         & /8/(n-1)**2)
    sz = SQRT(1.0_dp/(n-1) &
         & *(1+(4-r**2)/2/(n-1)+(22-6*r**2-3*r**4) &
         & /6/(n-1)**2))
    dzdr = 1.0_dp/(1+r)/(1-r)
  END SUBROUTINE calc_Fisher
!///////////////////////////////////////////////////////////////////////
! The integrand of the cumulative likelihood function
! exp(-(z(r)-zbar(rho))/sz(rho))**2/2)/sqrt(2Pi)/sz(rho)
  FUNCTION PHIZ1(Rho)
    USE cons
    USE like_vars
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: Rho
    REAL(dp) :: PHIZ1
    !
    REAL(dp) :: dzdr,zbar,sz
    ! sqrpi2 = 1/sqrt(2*Pi)
    REAL(dp), PARAMETER :: SQRPI2 = 1.0_dp/sqrt(2*Pi)
    !
    IF (ABS(Rho) >= 1.0_dp) THEN
       WRITE (*,*)  'PHIZ1: |rho| >= 1'
       STOP
    END IF
    !
    call calc_Fisher(Rho,n0,zbar,sz,dzdr)
    PHIZ1 = dzdr*SQRPI2/sz*EXP(-((z0-zbar)/sz)**2/2)
  END FUNCTION PHIZ1
!///////////////////////////////////////////////////////////////////////
! Integrating over the likelihood function from rho=-1 to rho=rho0 for a
! given value of the estimator r evaluated on n pairs.
  FUNCTION CLIKE(Rho0,Z1,N1)
    USE cons
    USE quad77
    USE like_vars
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: Rho0
    REAL(dp), INTENT(in) :: Z1
    INTEGER, INTENT(in) :: N1
    REAL(dp) :: CLIKE
    !
    IF (ABS(Rho0) > 1.0) THEN
       WRITE (*,*)  'CLIKE: |rho| > 1'
       STOP
    END IF
    ! Initializing PHIZ1
    n0 = N1
    z0 = Z1
    !
    CALL QROMO(PHIZ1,-1D0,Rho0,CLIKE,MIDPNT)
  END FUNCTION CLIKE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The integrand of the likelihood function
  ! i0 is the point whose peak likelihood is calculated
  FUNCTION LIKINT(Rho)
    USE cons
    USE like_vars
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: Rho
    REAL(dp) :: LIKINT
    !
    REAL(dp) :: cn,zbar,sz,dzdr
    INTEGER :: j
    ! cpi = 1/sqrt(2*Pi)
    REAL(dp), PARAMETER ::  CPI=1.0_dp/SQRT(2*Pi), TINY=1D-10
    !
    IF (ABS(Rho) >= 1.0) THEN
       WRITE (*,*) 'LIKINT: |rho| >= 1'
       STOP
    END IF
    !
    call calc_Fisher(Rho,n(i0),zbar,sz,dzdr)
    LIKINT = LOG(dzdr) + LOG(CPI) - LOG(sz) - ((z(i0)-zbar)/sz)**2/2
    DO j = 1 , i0-1
       cn = CLIKE(Rho,z(j),n(j))
       IF (cn <= TINY) THEN
          LIKINT = 0.0_dp
          RETURN
       END IF
       LIKINT = LIKINT + LOG(cn)
    END DO
    DO j = i0+1 , np
       cn = CLIKE(Rho,z(j),n(j))
       IF (cn <= TINY) THEN
          LIKINT = 0.0_dp
          RETURN
       END IF
       LIKINT = LIKINT + LOG(cn)
    END DO
    LIKINT = EXP(LIKINT)
  END FUNCTION LIKINT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE like_funcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM plike
  USE cons
  USE like_vars
  IMPLICIT NONE
  REAL(dp) :: mxlik,tmxlik
  REAL(dp) :: tlow,thigh,dtlikl,dtliku,tcom
  CHARACTER(len=100) :: dcf
  INTEGER :: ioerr,i,j
  REAL(dp) :: tdummy
  !
  WRITE (*,'("PLIKE V4.0 begins.")')
  WRITE (*,'("Enter dcf file name:")')
  READ (*,'(a100)') dcf
  OPEN (UNIT=11,FILE=dcf,STATUS='OLD',IOSTAT=ioerr)
  IF (ioerr /= 0) THEN
     WRITE (*,*) 'Error ',ioerr,' opening file ',TRIM(dcf)
     STOP
  END IF
  !
  WRITE (*,'("Enter lower bound on peak location:")')
  READ (*,*) tlow
  WRITE (*,'("Enter upper bound on peak location:")')
  READ (*,*) thigh
  WRITE (*,'(/,"Calculating ML in range t_lag = ",es10.3," to ",es10.3)') &
       tlow,thigh
  !
  ! Going through data first to find how many data points there are
  ! in *.dcf file in requested time-lag range.
  !
  np = 0
  DO
     READ (11,*,END=100) tdummy
     IF (tdummy >= tlow .AND. tdummy <= thigh) np = np + 1
  END DO
100 CONTINUE
  !
  WRITE (*,*)
  WRITE (*,'(i4," data points found in time-lag range")') np
  CALL allocate_vars(np)
  !
  ! Actually reading data
  !
  REWIND(11)
  i = 0
  DO WHILE (i < np)
     READ (11,*,END=200) t(i+1),dtm(i+1),dtp(i+1), &
          r(i+1),drm(i+1),drp(i+1),n(i+1)
     IF (t(i+1) >= tlow .AND. t(i+1) <= thigh) i = i+1
  END DO
200 CONTINUE
  CLOSE (UNIT=11,STATUS='KEEP',IOSTAT=ioerr)
  IF (ioerr /= 0) THEN
     WRITE (*,*) 'Error ',ioerr,' closing file ',TRIM(dcf)
     STOP
  END IF
  !
  WRITE (*,*)
  tcom = COM(i)
  WRITE (*,'("ZDCF peak at ",sp,es11.3," r_max = ",es11.3,&
       &" ZDCF C.O.M. at ",es11.3)') t(i),r(i),tcom
  WRITE (*,*)
  CALL PLIKE4(tmxlik,mxlik,dtlikl,dtliku,2)
  !
  WRITE (*,'("Results:",/)')
  OPEN (UNIT=11,FILE='plike.out',STATUS='UNKNOWN',IOSTAT=ioerr)
  IF (ioerr /= 0) THEN
     WRITE (*,*) 'Error ',ioerr,' opening file plike.out'
     STOP
  END IF
  DO j = 6,11,5
     WRITE (j,'("# PLIKE V4.0 ZDCF peak likelihood analysis",/,&
          &"# for file ",a)') TRIM(dcf)
     WRITE (j,'("#")')
     WRITE (j, &
          '("# num      lag         r         -dr        +dr    likelihood ")')
     WRITE (j, &
          '("# ---- ---------- ---------- ---------- ---------- ---------- ")')
     WRITE (j,'(2x,i4,1x,es10.3,1x,es10.3,1x,es10.3,1x,es10.3,1x,es10.3)') &
          (i,t(i),r(i),drm(i),drp(i),liklhd(i),i=1,np)
     !
     WRITE (j,'("#")')
     WRITE (j,'("# ML Peak at t_lag = ",sp,es11.3,", Likelihood at peak = ",es11.3)') &
          tmxlik,mxlik
     WRITE (j,'("# 1 sigma ML t_lag interval = ",sp,es11.3,1x,es11.3,1x,&
          &es11.3,/,"#",26x," = (",es11.3," , ",es11.3,")")') &
          tmxlik,dtliku,-dtlikl,tmxlik-dtlikl,tmxlik+dtliku
  END DO
  CLOSE (UNIT=11,STATUS='KEEP',IOSTAT=ioerr)
  IF (ioerr /= 0) THEN
     WRITE (*,*) 'Error ',ioerr,' closing file plike.out'
     STOP
  END IF
  !
  WRITE (*,*)
  WRITE (*,'("Output written on plike.out.")')
  CALL deallocate_vars()
  WRITE (*,'("Program ended.")')
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! com calculates the center-of-mass (centroid) around the ZDCF peak.
  ! The sum is taken over all the points lying between the two nearest
  ! the maximum with r < r_max/2.
  !
  FUNCTION COM(Imax)
    USE cons
    USE like_vars
    IMPLICIT NONE
    INTEGER, INTENT(out) :: Imax
    REAL(dp) :: com
    !
    REAL :: sumr,sumtr,rhalf
    INTEGER :: i,i1,i2
    ! Locating the maximum
    Imax = SUM(MAXLOC(R))
    !
    rhalf = R(Imax)/2.
    ! Locating the peak between half max to the left (i1) and right (i2)
    i1 = Imax - 1
    DO WHILE (i1 >= 1 .AND. R(i1) >= rhalf)
       i1 = i1 - 1
    END DO
    i1 = i1 + 1
    ! Warning
    IF (i1 == 1 .AND. R(1) >= rhalf) THEN
       WRITE (*,*) 'COM: Warning - Peak extends beyond range to the left.'
    END IF
    !
    i2 = Imax + 1
    DO WHILE (i2 <= np .AND. R(i2) >= rhalf)
       i2 = i2 + 1
    END DO
    i2 = i2 - 1
    ! Warning
    IF (i2 == np .AND. R(np) >= rhalf) THEN
       WRITE (*,*) 'COM: Warning - Peak extends beyond range to the right.'
    END IF
    ! Calculating the center of mass
    sumtr = SUM(T(i1:i2)*R(i1:i2))
    sumr = SUM(R(i1:i2))
    IF (sumr <= 0.0) THEN
       WRITE (*,*)  'COM: Problems. np = ',np ,' imax = ', Imax ,      &
            &' rmax = ',R(Imax),' i1 = ',i1,' i2 = ' , i2 ,   &
            &' sumr = ' , sumr
       WRITE (*,'(1x,f10.3,1x,f10.3)') (T(i),R(i),i=1,np)
       STOP
    ELSE
       COM = sumtr/sumr
    END IF
  END FUNCTION COM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE plike4(tmxlik,mxlik,dtlikl,dtliku,zrflag)
    USE cons
    USE quad77
    USE like_vars
    USE like_funcs
    IMPLICIT NONE
    INTEGER, INTENT(in) :: zrflag
    REAL(dp), INTENT(out) :: tmxlik,mxlik,dtlikl,dtliku
    !
    REAL(dp) :: rpeak
    INTEGER :: imxlik,ipeak
    REAL(dp), PARAMETER :: EPSLIK = 1e-10_dp
    !
    ipeak = SUM(MAXLOC(R))
    rpeak = R(ipeak)
    liklhd = 0.0_dp
    ! Input given in z-space (zrflag = 1): vector r holds z
    IF (zrflag == 1) THEN ! NOT IN USE!
       z(1:np) = R(1:np)
    ! Input given in r-space (rzflag <> 1)
    ELSE
       z(1:np) = LOG((1.+R(1:np))/(1.-R(1:np)))/2
    END IF
    !
    ! Calculating the likelihood function around the peak
    ! Saving time by stopping when negligibly low values reached.
    !
    WRITE (*,'("Progress meter: ",$)')
    DO i = ipeak-1,1,-1
       i0 = i ! Initializing LIKINT
       WRITE (*,'(".",$)')
       CALL QROMO(LIKINT,-1D0,1D0,liklhd(i),MIDPNT)
       IF (liklhd(i) < EPSLIK) EXIT
    END DO
    DO i = ipeak,np
       i0 = i ! Initializing LIKINT
       WRITE (*,'(".",$)')
       CALL QROMO(LIKINT,-1D0,1D0,liklhd(i),MIDPNT)
       IF (liklhd(i) < EPSLIK) EXIT
    END DO
    WRITE (*,'(/)')
    !
    ! Finding the maximal likelihood value and error interval
    !
    CALL CONF1(imxlik,tmxlik,mxlik,dtlikl,dtliku)
    !
  END SUBROUTINE plike4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finding the +/- 1 sigma intervals on the interpolated likelihood
! function which (after normalization) is the "Fiducial Probability
! Function" (Kendall & Stuart, Vol 2, 3rd ed. p. 141)
  SUBROUTINE CONF1(imxlik,tmxlik,mxlik,dtlikl,dtliku)
    !
    USE cons
    USE like_vars
    IMPLICIT NONE
    INTEGER, INTENT(out) :: imxlik
    REAL(dp), INTENT(out) :: tmxlik,mxlik,dtlikl,dtliku
    !
    INTEGER :: i
    REAL(dp) :: cn,s,p,t1sigl,t1sigu
    REAL(dp) :: Pcum(np)
    LOGICAL :: edge
    !
    ! Find the maximal likelihood
    !
    imxlik = SUM(MAXLOC(liklhd(1:np)))
    tmxlik = T(imxlik)
    mxlik = liklhd(imxlik)
    !
    ! Find the +/- 1 sigma limits
    !
    cn = 0.0
    Pcum(1) = 0.0
    DO i = 2 , Np
       cn = cn + (liklhd(i)+liklhd(i-1))/2.0*(T(i)-T(i-1))
       Pcum(i) = cn
    END DO
    Pcum(1:np) = Pcum(1:np)/cn
    !
    edge = .TRUE.
    p = Pcum(imxlik) + 0.3413*(Pcum(Np)-Pcum(imxlik))*2
    t1sigu = T(Np)
    DO i = imxlik,Np
       IF (Pcum(i) >= p) THEN
          s = (liklhd(i)-liklhd(i-1))/(T(i)-T(i-1))
          t1sigu = T(i-1) &
               & + (-liklhd(i-1)+SQRT(liklhd(i-1)**2+2*(p-Pcum(i-1))&
               & *s*cn))/s
          edge = .FALSE.
          EXIT
       END IF
    END DO
    IF (edge) WRITE (*,*) &
         'CONF1: Warning - 1 sigma interval may extend beyond range to the right.'
    !
    edge = .TRUE.
    p = Pcum(imxlik) - 0.3413*(Pcum(imxlik)-Pcum(1))*2
    t1sigl = T(1)
    DO i = imxlik ,1,-1
       IF (Pcum(i) <= p) THEN
          s = (liklhd(i+1)-liklhd(i))/(T(i+1)-T(i))
          t1sigl = T(i) &
               & + (-liklhd(i)+SQRT(liklhd(i)**2+2*(p-Pcum(i))*s*cn)&
               & )/s
          edge = .FALSE.
          EXIT
       END IF
    END DO
    IF (edge) WRITE (*,*) &
         'CONF1: Warning - 1 sigma interval may extend beyond range to the left.'
    !
    dtliku = MAX(t1sigu-tmxlik,dtp(imxlik))
    dtlikl = max(tmxlik-t1sigl,dtm(imxlik))
  END SUBROUTINE CONF1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM plike
