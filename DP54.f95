! -----------------------------------------------------------------------------
! DP54.f95 / MODULE DP54
!
! Adaptive Dormand–Prince Runge–Kutta integrator of order 5 with embedded order 4
! (often called "RK45" / "DP5(4)").
!
! Public routine:
!   DP54_try_step : attempt one step, return accept/reject + new step suggestion.
!
! Features:
!   - FSAL ("First Same As Last") reuse of the last stage slope on accepted steps
!   - Scaled RMS error norm using atol + rtol*max(|y|,|y5|)
!   - Basic finite-check guards against NaNs/Infs
! -----------------------------------------------------------------------------

MODULE DP54

    use, intrinsic :: iso_fortran_env, only: wp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    IMPLICIT NONE

    private
    public :: DP54_try_step

    ! User must provide RHS: dydt = f(t, y)
    abstract interface
        subroutine rhs_fun(t, y, dydt)
            import :: wp
            real(wp), intent(in) :: t
            real(wp), intent(in) :: y(:)
            real(wp), intent(out) :: dydt(size(y))
        end subroutine rhs_fun
    end interface

    CONTAINS

        ! ------------------------------------------------------------
        ! Helper: returns true if every element of v is finite
        ! (not NaN, not +Inf, not -Inf)
        ! ------------------------------------------------------------
        pure logical function all_finite(v) result(ok)
            real(wp), intent(in) :: v(:)
            ok = all(ieee_is_finite(v))
        end function all_finite

        subroutine DP54_try_step(f, t, y, h, rtol, atol, &
                                have_fsal, k_fsal, yt, &
                                k1, k2, k3, k4, k5, k6, &
                                k7, y5, y4, sc, e, err_norm, &
                                h_new, accepted)

            ! Attempt a single adaptive Dorman-Prince 5(4) step of size h.
            !
            ! Inputs
            !   f   : RHS function, dydt = f(t, y)
            !   t   : Current independent variable (time / affine parameter / etc.)
            !   y   : Current state vector at t
            !   h   : Proposed step size
            !   rtol,atol : Error tolerances (scalar version)
            !
            ! FSAL inputs/outputs:
            !   have_fsal   : if true, k_fsal already equals f(t,y) from previous accepted step
            !   k_fsal      : cache for FSAL slope (stores previous k7)
            !
            ! Workspace arrays (preallocated by caller; all same size as y):
            !   yt          : temporary "trail state" used to compute each stage
            !   k1 ... k7   : stage slope vectors
            !   y5          : 5th order candidate solution at t+h
            !   y4          : 4th order embedded solution at t+h (for error estimation)
            !   sc          : pre-component scaling for error (atol + rtol * something)
            !   e           : pre-component error estimate (y5 - y4)
            !
            ! Outputs:
            !   err_norm    : scalar error measure (scaled RMS)
            !   h_new       : suggested next step size
            !   accepted    : true if err_norm <= 1, and y is updated to y5

            procedure(rhs_fun) :: f
            real(wp), intent(in) :: t, h, rtol, atol
            real(wp), intent(inout) :: y(:)

            logical, intent(inout) :: have_fsal
            real(wp), intent(inout) :: k_fsal(:)

            real(wp), intent(inout) :: yt(:)
            real(wp), intent(inout) :: k1(:), k2(:), k3(:), k4(:), k5(:), k6(:), k7(:)
            real(wp), intent(inout) :: y5(:), y4(:), sc(:), e(:)

            real(wp), intent(out) :: err_norm, h_new
            logical, intent(out) :: accepted

            integer :: n
            real(wp) :: fac

            ! -----------------------------
            ! Dormand-Prince 5(4) constants
            ! -----------------------------
            ! c_i are the fractional time locations inside the step: t + c_i * h
            real(wp), parameter :: c2 = 1._wp/5._wp, c3 = 3._wp/10._wp, c4 = 4._wp/5._wp
            real(wp), parameter :: c5 = 8._wp/9._wp, c6 = 1._wp, c7 = 1._wp

            ! a_ij coefficients define the trial state for stage i:
            ! yt = y + h * (sum_j a_ij ) k_j)   for j < i
            real(wp), parameter :: a21 = 1._wp/5._wp
            real(wp), parameter :: a31 = 3._wp/40._wp, a32 = 9._wp/40._wp
            real(wp), parameter :: a41 = 44._wp/45._wp, a42 = -56._wp/15._wp, a43 = 32._wp/9._wp
            real(wp), parameter :: a51 = 19372._wp/6561._wp, a52 = -25360._wp/2187._wp, &
                                    a53 = 64448._wp/6561._wp, a54 = -212._wp/729._wp
            real(wp), parameter :: a61 = 9017._wp/3168._wp, a62 = -355._wp/33._wp, &
                                    a63 = 46732._wp/5247._wp, a64 = 49._wp/176._wp, a65 = -5103._wp/18656._wp
            real(wp), parameter :: a71 = 35._wp/384._wp, a72 = 0._wp, &
                                    a73 = 500._wp/1113._wp, a74 = 125._wp/192._wp, &
                                    a75 = -2187._wp/6784._wp, a76 = 11._wp/84._wp

            ! b coefficients for the 5th-order solution (main output).
            ! NOTE: DP has the nice property that this combination is ALSO the stage-7 trial state
            real(wp), parameter :: b1 = 35._wp/384._wp, b3 = 500._wp/1113._wp, b4 = 125._wp/192._wp
            real(wp), parameter :: b5 = -2187._wp/6784._wp, b6 = 11._wp/84._wp

            ! bh coefficients for the embedded 4th-order solution (used only for error estimation)
            real(wp), parameter :: bh1 = 5179._wp/57600._wp, bh3 = 7571._wp/16695._wp, bh4 = 393._wp/640._wp
            real(wp), parameter :: bh5 = -92097._wp/339200._wp, bh6 = 187._wp/2100._wp, bh7 = 1._wp/40._wp

            ! -----------------------------
            ! Step-size controller settings
            ! -----------------------------
            ! safety   : be conservative to avoid repeated rejections
            ! fac_min  : don't shrink step by more than this factor at once
            ! fac_max  : don't grow step by more than this factor at once
            ! p        : order used in h update (use 5 for DP5)
            real(wp), parameter :: safety = 0.9_wp, fac_min = 0.2_wp, fac_max = 5.0_wp, p = 5.0_wp

            n = size(y)

            ! ------------------------------------------------------------
            ! Stage 1: k1 = f(t, y)
            ! If have_fsal is true, reuse k_fsal from previous accepted step
            ! (FSAL = "First Same As Last").
            ! ------------------------------------------------------------
            if (have_fsal) then
                k1 = k_fsal
            else
                call f(t, y, k1)
                if (.not. all_finite(k1)) then
                    accepted = .false.
                    have_fsal = .false.
                    err_norm = huge(1.0_wp)
                    h_new = h * fac_min
                    return
                end if
            end if

            ! ------------------------------------------------------------
            ! Stage 2: evaluate f at t + c2 * h, with trail state y + h * a21 * k1
            ! ------------------------------------------------------------
            yt = y + h * (a21 * k1)
            call f(t + c2 * h, yt, k2)
            if (.not. all_finite(k2)) then
                accepted = .false.
                have_fsal = .false.
                err_norm = huge(1.0_wp)
                h_new = h * fac_min
                return
            end if

            ! Stage 3
            yt = y + h * (a31 * k1 + a32 * k2)
            call f(t + c3 * h, yt, k3)
            if (.not. all_finite(k3)) then
                accepted = .false.
                have_fsal = .false.
                err_norm = huge(1.0_wp)
                h_new = h * fac_min
                return
            end if

            ! Stage 4
            yt = y + h * (a41 * k1 + a42 * k2 + a43 * k3)
            call f(t + c4 * h, yt, k4)
            if (.not. all_finite(k4)) then
                accepted = .false.
                have_fsal = .false.
                err_norm = huge(1.0_wp)
                h_new = h * fac_min
                return
            end if

            ! Stage 5
            yt = y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
            call f(t + c5 * h, yt, k5)
            if (.not. all_finite(k5)) then
                accepted = .false.
                have_fsal = .false.
                err_norm = huge(1.0_wp)
                h_new = h * fac_min
                return
            end if

            ! Stage 6
            yt = y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
            call f(t + c6 * h, yt, k6)
            if (.not. all_finite(k6)) then
                accepted = .false.
                have_fsal = .false.
                err_norm = huge(1.0_wp)
                h_new = h * fac_min
                return
            end if

            ! ------------------------------------------------------------
            ! 5th-order candiate y5 at t + h:
            ! This is the main solution you'd keep if the step is accepted.
            ! Also equals the "stage 7 trail state" in Dormand-Prince.
            ! ------------------------------------------------------------
            y5 = y + h * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)
            if (.not. all_finite(y5)) then
                accepted = .false.
                have_fsal = .false.
                err_norm = huge(1.0_wp)
                h_new = h * fac_min
                return
            end if

            ! Stage 7 slope, evaluated at the accepted-time location t + h,
            ! but at the trail state y5:
            call f(t + c7 * h, y5, k7)
            if (.not. all_finite(k7)) then
                accepted = .false.
                have_fsal = .false.
                err_norm = huge(1.0_wp)
                h_new = h * fac_min
                return
            end if

            ! ------------------------------------------------------------
            ! 4th-order embedded solution y4 at t + h:
            ! It's not used as the final result; it exists to estimate error.
            ! ------------------------------------------------------------
            y4 = y + h * (bh1 * k1 + bh3 * k3 + bh4 * k4 + bh5 * k5 + bh6 * k6 + bh7 * k7)
            if (.not. all_finite(y4)) then
                accepted = .false.
                have_fsal = .false.
                err_norm = huge(1.0_wp)
                h_new = h * fac_min
                return
            end if

            ! Error estimate vector (per component)
            e = y5 - y4

            ! ------------------------------------------------------------
            ! Convert vector error e into a single number err_norm
            ! sc(i) sets how big an error is "allowed" in component i.
            !
            ! This default scaling works well if your variables are reasonably scaled.
            ! Later you may replace atol with a vector or do nondimensionalization.
            ! ------------------------------------------------------------
            sc = atol + rtol * max(abs(y), abs(y5))

            ! Scaled RMS norm (dimensionless):
            ! err_norm ~ 1 means "right at tolerance"
            err_norm = sqrt(sum((e/sc)**2) / real(n, wp))

            ! Accept/rejections
            if (err_norm <= 1.0_wp) then
                accepted = .true.

                !if accepted, comit the step
                y = y5

                ! Cache FSAL slope for next step:
                ! Next step starts at (t + h, y5), and k7 = f(t + h, y5),
                ! so k7 can be reused as k1 on the next accepted step.
                k_fsal = k7
                have_fsal = .true.
            else
                accepted = .false.

                ! If we reject the step, we will retry with a new (smaller) h.
                ! It's safest to disable FSAL reuse because (t, y) is unchanged,
                ! but k_fsal corresponds to f(t + h, y5) which we did not accept.
                have_fsal = .false.
            end if

            ! ------------------------------------------------------------
            ! Choose new step size h_new for the next attempt.
            ! Standard controller: h_new = h * safety * err^(1/p)
            ! ------------------------------------------------------------
            if (err_norm == 0.0_wp) then
                fac = fac_max
            else
                fac = safety * err_norm**(-1.0_wp/p)
                fac = min(fac_max, max(fac_min, fac))
            end if

            h_new = h * fac

        end subroutine DP54_try_step

END MODULE DP54
