! -----------------------------------------------------------------------------
!  kerr_newman_physics.f95 / MODULE kerr_newman_physics
! -----------------------------------------------------------------------------

MODULE kerr_newman_physics

    use, intrinsic :: iso_fortran_env, only: wp => real64

    IMPLICIT NONE

    real(wp), parameter :: M = 1.0_wp                       !Mass
    real(wp), parameter :: a = 0.9_wp                       !Spin
    real(wp), parameter :: Q = 0.0_wp                       !Charge


CONTAINS

    pure real(wp) function B_of_r(r) result(B)
        real(wp), intent(in) :: r
        B = 2.0_wp * M * r - Q * Q
    end function B_of_r

    pure real(wp) function Delta_of_r(r) result(Delta)
        real(wp), intent(in) :: r
        Delta = r * r - 2.0_wp * M * r + a * a + Q * Q
    end function Delta_of_r

    pure real(wp) function Sigma_of_rth(r, th) result(Sigma)
        real(wp), intent(in) :: r, th
        real(wp) :: cth
        cth = cos(th)
        Sigma = r * r + a * a * cth * cth
    end function Sigma_of_rth


    pure real(wp) function hamiltonian_null(y) result(H)
        real(wp), intent(in) :: y(:)

        real(wp) :: r, th, pt, pr, pth, pph
        real(wp) :: B, Delta, Sigma, sth, sth2, gtt, grr, gthth, gphph, g_t_phi
        real(wp) :: N, D, P

        ! unpack what we need (ignore t, phi for H)
        r = y(2)
        th = y(3)
        pt = y(5)
        pr = y(6)
        pth = y(7)
        pph = y(8)

        sth = sin(th)
        sth2 = sth * sth
        !Guard at the poles.
        sth2 = max(sth2, 1.0e-30_wp)

        B = B_of_r(r)
        Delta = Delta_of_r(r)
        Sigma = Sigma_of_rth(r, th)
        P = (r * r + a * a) * (r * r + a * a) - a * a * Delta * sth2

        gtt = - P / (Delta * Sigma)
        grr = Delta / Sigma
        gthth = 1.0_wp / Sigma

        !Numerator for gphph
        N = Delta - a * a * sth2
        ! Denominator for gphph
        D = Delta * Sigma * sth2

        ! Guard agains sin(O/) = O/ at poles
        gphph = N / D

        g_t_phi = - (a * B) / (Delta * Sigma)

        ! Hamiltonian = 1/2 (g^{tt} pt^2 g^{rr} pr^2 g^{O/ O/} pO/^2 + g^{O|O|) pO|^2 + 2 g^{O/ O|} pO\ pO|)
        H = 0.5_wp * (gtt * pt * pt + grr * pr * pr + gthth * pth * pth + gphph * pph * pph + 2 * g_t_phi * pt * pph)
    end function hamiltonian_null

    subroutine rhs_kerr_newman(lambda, y, dydlambda)
        real(wp), intent(in) :: lambda
        real(wp), intent(in) :: y(:)
        real(wp), intent(out) :: dydlambda(size(y))

        real(wp) :: t, r, th, phi, pt, pr, pth, pph
        real(wp) :: sth, cth, sth2, B, Delta, Sigma

        real(wp) :: gtt, grr, gthth, gphph, g_t_phi
        real(wp) :: dr_gtt, dr_grr, dr_gthth, dr_gphph, dr_g_t_phi
        real(wp) :: dth_gtt, dth_grr, dth_gthth, dth_gphph, dth_g_t_phi

        real(wp) :: dr_Sigma, dth_Sigma, dr_Delta, dr_B, dr_D, dth_D, dr_N, dth_N

        real(wp) :: dr_Q, dth_Q, dr_P, dth_P

        real(wp) :: P, N, D, DS

        t = y(1)
        r = y(2)
        th = y(3)
        phi = y(4)
        pt = y(5)
        pr = y(6)
        pth = y(7)
        pph = y(8)

        sth = sin(th)
        sth2 = sth * sth
        !Guard at the poles.
        sth2 = max(sth2, 1.0e-30_wp)
        cth = cos(th)

        B = B_of_r(r)
        Delta = Delta_of_r(r)
        Sigma = Sigma_of_rth(r, th)

        dr_Sigma = 2.0_wp * r
        dth_Sigma = -2.0_wp * a * a * sth * cth

        dr_Delta = 2.0_wp * r - 2.0_wp * M
        dr_N = dr_Delta

        dth_N = -2.0_wp * a * a * sth * cth
        dth_D = Delta * (dth_Sigma * sth2 + Sigma * (2.0_wp * sth * cth))
        dr_D = (dr_Delta * Sigma + Delta * dr_Sigma) * sth2

        P = (r * r + a * a) * (r * r + a * a) - a * a * Delta * sth2
        DS = Delta * Sigma

        gtt = - P / DS
        grr = Delta / Sigma
        gthth = 1.0_wp / Sigma
        g_t_phi = - (a * B) / DS

        dr_B = 2.0_wp * M
        dr_g_t_phi = -a * (dr_B / DS - B * (dr_Delta * Sigma + Delta * dr_Sigma) / (DS * DS))
        dth_g_t_phi = a * B * dth_Sigma / (Delta * Sigma * Sigma)

        dr_gthth = - dr_Sigma / (Sigma * Sigma)
        dth_gthth = - dth_Sigma / (Sigma * Sigma)

        dr_grr = (dr_Delta * Sigma - Delta * dr_Sigma) / (Sigma * Sigma)
        dth_grr = -(Delta * dth_Sigma) / (Sigma * Sigma)

        dr_P = 4.0_wp * r * (r * r + a * a) - a * a * dr_Delta * sth2
        dth_P = -2.0_wp * a * a * Delta * sth * cth

        dr_Q = dr_Delta * Sigma + Delta * dr_Sigma
        dth_Q = Delta * dth_Sigma

        dr_gtt = -(dr_P * DS - P * dr_Q) / (DS * DS)
        dth_gtt = -(dth_P * DS - P * dth_Q) / (DS * DS )

        !Numerator for gphph
        N = Delta - a * a * sth2
        ! Denominator for gphph
        D = Delta * Sigma * sth2

        dr_gphph = (dr_N * D - N * dr_D) / (D * D)
        dth_gphph = (dth_N * D - N * dth_D) / (D * D)

        ! Guard agains sin(O/) = O/ at poles
        gphph = N / D

        dydlambda(1) = gtt * pt + g_t_phi * pph
        dydlambda(2) = grr * pr
        dydlambda(3) = gthth * pth
        dydlambda(4) = g_t_phi * pt + gphph * pph
        dydlambda(5) = 0.0_wp
        dydlambda(6) = -0.5_wp * (dr_gtt * pt * pt + 2 * dr_g_t_phi * pt * pph + &
                        dr_gphph * pph * pph + dr_grr * pr * pr + dr_gthth * pth * pth)
        dydlambda(7) = -0.5_wp * (dth_gtt * pt * pt + 2 * dth_g_t_phi * pt * pph + &
                        dth_gphph * pph * pph + dth_grr * pr * pr + dth_gthth * pth * pth)
        dydlambda(8) = 0.0_wp

    end subroutine rhs_kerr_newman


END MODULE kerr_newman_physics
