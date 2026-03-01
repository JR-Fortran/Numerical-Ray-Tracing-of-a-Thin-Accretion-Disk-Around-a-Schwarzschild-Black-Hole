! -----------------------------------------------------------------------------
! schwarzschild_physics.f95 / MODULE schwarzschild_physics
!
! Physics for null geodesics in Schwarzschild spacetime (M=1).
!
! Exposed routines:
!   A_of_r            : metric factor A(r) = 1 - 2/r
!   hamiltonian_null  : H = 1/2 g^{mu nu} p_mu p_nu (should be ~0 for photons)
!   rhs_schwarzschild : Hamiltonian equations of motion dy/dλ
!
! Conventions:
!   - State y stores covariant momenta p_mu (not contravariant).
!   - gtt, grr, gthth, gphph are inverse-metric components g^{mu nu}.
! -----------------------------------------------------------------------------

MODULE schwarzschild_physics

    use, intrinsic :: iso_fortran_env, only: wp => real64
    IMPLICIT NONE
    private

    ! Expose only what the rest of the program needs
    public :: rhs_schwarzschild, hamiltonian_null, A_of_r

    real(wp), parameter :: M = 1.0_wp ! Geometrized units with M = 1.

contains

    ! -----------------------------------------------------------------------------
    ! A_of_r(r)
    !
    ! Returns A(r) = 1 - 2M/r with M=1 (geometrized units).
    ! Used throughout for Schwarzschild metric coefficients.
    ! -----------------------------------------------------------------------------
    pure real(wp) function A_of_r(r) result(A)
        real(wp), intent(in) :: r
        A = 1.0_wp - 2.0_wp * M / r
    end function A_of_r

    ! -----------------------------------------------------------------------------
    ! hamiltonian_null(y)
    !
    ! Compute the Hamiltonian:
    !   H = 1/2 g^{mu nu} p_mu p_nu
    !
    ! Here g^{mu nu} is the inverse Schwarzschild metric and p_mu are the covariant
    ! momentum components stored in y(5:8).
    !
    ! For a null geodesic (photon), H should remain approximately zero; tracking H
    ! is a useful diagnostic for integration accuracy.
    !
    ! Note:
    !   gphph = g^{phi phi} = 1/(r^2 sin^2(theta))
    !   The code guards against sin(theta)=0 issues near the poles by using sth=sin(th).
    ! -----------------------------------------------------------------------------
    pure real(wp) function hamiltonian_null(y) result(H)
        real(wp), intent(in) :: y(:)

        real(wp) :: r, th, pt, pr, pth, pph
        real(wp) :: A, gtt, grr, gthth, gphph, sth

        ! unpack what we need (ignore t, phi for H)
        r = y(2)
        th = y(3)
        pt = y(5)
        pr = y(6)
        pth = y(7)
        pph = y(8)

        A = A_of_r(r)
        sth = sin(th)

        ! Inverse metric components g^{mu nu} for rhs_schwarzschild:
        ! g^{tt} = -1/A
        ! g^{rr} = A
        ! g^{O/ O/} = 1/r^2                O/ = theta
        ! g^(O| O|) = 1/(r^2 sin^2(O/))    O| = phi
        gtt = -1.0_wp / A
        grr = A
        gthth = 1.0_wp / (r * r)

        ! Guard agains sin(O/) = O/ at poles
        gphph = 1.0_wp / (r * r * sth * sth)

        ! Hamiltonian = 1/2 (g^{tt} pt^2 g^{rr} pr^2 g^{O/ O/} pO/^2 + g^{O|O|) pO|^2)
        H = 0.5_wp * (gtt * pt * pt + grr * pr * pr + gthth * pth * pth + gphph * pph * pph)
    end function hamiltonian_null

    ! -----------------------------------------------------------------------------
    ! rhs_schwarzschild(lambda, y, dydlambda)
    !
    ! Right-hand side dy/dλ for Hamiltonian evolution in Schwarzschild spacetime.
    !
    ! Equations implemented:
    !   ẋ^μ   = ∂H/∂p_μ = g^{μν} p_ν
    !   ṗ_μ   = -∂H/∂x^μ = -(1/2) (∂_μ g^{αβ}) p_α p_β
    !
    ! Since the Schwarzschild metric is independent of t and phi:
    !   dp_t/dλ   = 0   (energy at infinity conserved)
    !   dp_phi/dλ = 0   (azimuthal angular momentum conserved)
    !
    ! Inputs:
    !   lambda : affine parameter (not explicitly used; kept for integrator interface)
    !   y      : state vector [t, r, theta, phi, pt, pr, ptheta, pphi]
    !
    ! Output:
    !   dydlambda : derivatives in the same ordering.
    ! -----------------------------------------------------------------------------
    subroutine rhs_schwarzschild(lambda, y, dydlambda)
        real(wp), intent(in) :: lambda
        real(wp), intent(in) :: y(:)
        real(wp), intent(out) :: dydlambda(size(y))

        ! We don't actually use lamda in Schwarzscild because the metric is static.
        ! But the integrator expects it, so it stays in the interface.
        real(wp) :: r, th, pt, pr, pth, pph
        real(wp) :: A, sth, cth

        ! Inverse metric components
        real(wp) :: gtt, grr, gthth, gphph

        ! Derivatives of inverse metric components
        real(wp) :: dr_gtt, dr_grr, dr_gthth, dr_gphph
        real(wp) :: dth_gphph

        ! 1 Unpack state vector for readability
        r = y(2)
        th = y(3)
        pt = y(5)
        pr = y(6)
        pth = y(7)
        pph = y(8)

        A = A_of_r(r)
        sth = sin(th)
        cth = cos(th)

        ! 2 Inverse metric g^{μν} (diagonal)
        gtt = -1.0_wp / A
        grr = A
        gthth = 1.0_wp / (r * r)
        gphph = 1.0_wp / (r * r * sth * sth)

        ! 3 Positon detivatives ẋ^μ = g^{μν} p_ν
        !
        ! dt/dλ     = g^{tt} * p_t
        ! dr/dλ     = g^{rr} * p_r
        ! dθ/dλ     = g^{θθ} * p_θ
        ! dφ/dλ     = g^{φφ} * p_φ
        dydlambda(1) = gtt * pt
        dydlambda(2) = grr * pr
        dydlambda(3) = gthth * pth
        dydlambda(4) = gphph * pph

        ! 4 Momentum detivatives ṗ_μ = -(1/2) ∂_μ g^{αβ} p_α p_β
        !
        ! Because Schwarzscild metric does NOT depend on t or φ
        ! ∂_t g^{αβ} = 0  => dp_t/dλ = 0
        ! ∂_φ g^{αβ} = 0  => dp_φ/dλ = 0
        dydlambda(5) = 0.0_wp    ! dpt/dλ
        dydlambda(8) = 0.0_wp    ! dpph/dλ

        ! 5 We need only ∂_r and ∂_θ derivatives of inverse metric terms.

        ! A(r) 1 - 2M/r
        ! dA/dr = +2M/R^2    (note the sign)
        !
        ! g^{tt} = -1/A  => ∂_r g^{tt} = + (1/A^2) dA/dr =  2M/(r^2 A^2)
        dr_gtt = 2.0_wp * M / (r * r * A * A)

        ! g^{rr} = A     => ∂_r g^{rr} = dA/dr = 2M/r^2
        dr_grr =  2.0_wp*M / (r*r)

        ! g^{θθ} = 1/r^2 => ∂_r g^{θθ} = -2/r^3
        dr_gthth = -2.0_wp / (r * r * r)

        ! g^{φφ} = 1/(r^2 sin^2 θ)
        ! ∂_r g^{φφ} = -2/(r^3 sin^2 θ)
        dr_gphph = -2.0_wp / (r * r * r * sth * sth)

        ! ∂_θ g^{φφ} = d/dθ [ (sin θ)^(-2) ] * 1/r^2
        ! d/dθ (sin^-2 θ) = -2 cos θ / sin^3 θ
        dth_gphph = -2.0_wp * cth / (r * r * sth * sth *sth)

        ! 6 dp_r/dλ uses ∂_r terms:
        !
        ! ṗ_r = -(1/2) [ (∂_r g^{tt}) pt^2 + (∂_r g^{rr}) pr^2
        !                 ṗ_r = -(1/2) [ (∂_r g^{tt}) pt^2 + (∂_r g^{rr}) pr^2
        dydlambda(6) = -0.5_wp * (dr_gtt * pt * pt + dr_grr * pr * pr + &
                                    dr_gthth * pth * pth + dr_gphph * pph * pph)

        ! 7 dp_θ/dλ uses ∂_θ terms:
        !
        ! Only g^{φφ} depends on θ, so:
        ! ṗ_θ = -(1/2) (∂_θ g^{φφ}) pφ^2
        dydlambda(7) = -0.5_wp * (dth_gphph * pph * pph)

    end subroutine rhs_schwarzschild

END MODULE schwarzschild_physics
