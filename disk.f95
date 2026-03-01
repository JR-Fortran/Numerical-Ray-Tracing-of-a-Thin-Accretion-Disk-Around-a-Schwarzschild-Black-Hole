! -----------------------------------------------------------------------------
! disk.f95 / MODULE disk
!
! Disk intersection utilities for a thin equatorial accretion disk.
!
! The disk surface is the equatorial plane:
!   theta = pi/2
!
! Provided routines:
!   disk_crossing   : detect and linearly interpolate a plane crossing
!   disk_refine_hit : refine the crossing by bisection + RK4 substeps
! -----------------------------------------------------------------------------

MODULE disk
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none
  private
  public :: disk_crossing, disk_refine_hit

contains

! -----------------------------------------------------------------------------
! disk_crossing(y_prev, y_new, hit, r_hit, phi_hit, s_hit)
!
! Detect whether the segment between two accepted integration states crosses
! the equatorial plane theta = pi/2.
!
! If a sign change in (theta - pi/2) occurs:
!   - s_hit is the interpolation fraction in [0,1]
!   - r_hit and phi_hit are linearly interpolated at the crossing
!   - phi interpolation unwraps across +/-pi to avoid discontinuities.
!
! Inputs:
!   y_prev, y_new : two consecutive accepted states (size 8)
!
! Outputs:
!   hit    : true if a crossing is detected
!   r_hit  : interpolated r at crossing
!   phi_hit: interpolated phi at crossing (unwrapped)
!   s_hit  : interpolation fraction along the segment
! -----------------------------------------------------------------------------

  subroutine disk_crossing(y_prev, y_new, hit, r_hit, phi_hit, s_hit)
    ! Detect crossing of equatorial plane theta = pi/2 between two accepted states.
    ! Also interpolate r and phi at the crossing.

    real(wp), intent(in)  :: y_prev(8), y_new(8)
    logical,  intent(out) :: hit
    real(wp), intent(out) :: r_hit, phi_hit, s_hit

    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp) :: th0, th1, a0, a1, denom
    real(wp) :: phi0, phi1, dphi

    th0 = y_prev(3)
    th1 = y_new(3)

    a0 = th0 - 0.5_wp*pi
    a1 = th1 - 0.5_wp*pi

    hit = .false.
    s_hit = 0.0_wp
    r_hit = 0.0_wp
    phi_hit = 0.0_wp

    ! Exact hit at endpoints (rare)
    if (a0 == 0.0_wp) then
      hit = .true.
      s_hit = 0.0_wp
      r_hit = y_prev(2)
      phi_hit = y_prev(4)
      return
    end if
    if (a1 == 0.0_wp) then
      hit = .true.
      s_hit = 1.0_wp
      r_hit = y_new(2)
      phi_hit = y_new(4)
      return
    end if

    ! Sign change indicates crossing
    if (a0*a1 < 0.0_wp) then
      denom = (a0 - a1)
      s_hit = a0 / denom

      r_hit = y_prev(2) + s_hit*(y_new(2) - y_prev(2))

      ! phi interpolation with unwrap across +/-pi
      phi0 = y_prev(4)
      phi1 = y_new(4)
      dphi = phi1 - phi0
      if (dphi >  pi) dphi = dphi - 2.0_wp*pi
      if (dphi < -pi) dphi = dphi + 2.0_wp*pi
      phi_hit = phi0 + s_hit*dphi

      hit = .true.
    end if

  end subroutine disk_crossing

  ! -----------------------------------------------------------------------------
  ! disk_refine_hit(rhs, lam0, y0, lam1, y1, niter, lam_hit, y_hit)
  !
  ! Refine an equatorial plane crossing inside [lam0, lam1] assuming that
  ! (theta - pi/2) changes sign across the interval.
  !
  ! Method:
  !   - Perform niter bisection iterations on lambda.
  !   - Use a small RK4 stepper (rk4_step) to propagate from the left endpoint
  !     to the midpoint to evaluate theta there.
  !   - Return lam_hit ~ crossing lambda and y_hit ~ state near the crossing.
  !
  ! This is used to make the disk edge sharper than simple linear interpolation.
  ! -----------------------------------------------------------------------------

  subroutine disk_refine_hit(rhs, lam0, y0, lam1, y1, niter, lam_hit, y_hit)
    implicit none

    interface
      subroutine rhs(lam, y, dydlam)
        import :: wp
        real(wp), intent(in)  :: lam
        real(wp), intent(in)  :: y(:)
        real(wp), intent(out) :: dydlam(size(y))
      end subroutine rhs
    end interface

    real(wp), intent(in)  :: lam0, lam1
    real(wp), intent(in)  :: y0(:), y1(:)
    integer,  intent(in)  :: niter
    real(wp), intent(out) :: lam_hit
    real(wp), intent(out) :: y_hit(size(y0))

    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp) :: aL, aR, aM, lamL, lamR, lamM, h
    integer  :: k, n

    real(wp) :: yL(size(y0)), yR(size(y0)), yM(size(y0))

    n = size(y0)
    if (size(y1) /= n) error stop "disk_refine_hit: size mismatch y0/y1"

    lamL = lam0; yL = y0
    lamR = lam1; yR = y1

    aL = yL(3) - 0.5_wp*pi
    aR = yR(3) - 0.5_wp*pi

    ! assumes aL*aR < 0 on entry
    do k = 1, niter
      lamM = 0.5_wp*(lamL + lamR)
      h    = lamM - lamL

      call rk4_step(lamL, yL, h, yM)
      aM = yM(3) - 0.5_wp*pi

      if (aL*aM <= 0.0_wp) then
        lamR = lamM; yR = yM; aR = aM
      else
        lamL = lamM; yL = yM; aL = aM
      end if
    end do

    lam_hit = 0.5_wp*(lamL + lamR)
    y_hit   = yR

  contains

    subroutine rk4_step(lam, y, h, yout)
        real(wp), intent(in)  :: lam, y(:), h
        real(wp), intent(out) :: yout(size(y))
        real(wp) :: k1(size(y)), k2(size(y)), k3(size(y)), k4(size(y))
        real(wp) :: yt(size(y))

        call rhs(lam, y, k1)
        yt = y + 0.5_wp*h*k1; call rhs(lam + 0.5_wp*h, yt, k2)
        yt = y + 0.5_wp*h*k2; call rhs(lam + 0.5_wp*h, yt, k3)
        yt = y + h*k3;        call rhs(lam + h, yt, k4)
        yout = y + (h/6.0_wp)*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
    end subroutine rk4_step

  end subroutine disk_refine_hit

END MODULE disk
