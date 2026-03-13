! -----------------------------------------------------------------------------
! camera.f95 / MODULE camera
!
! Builds initial null rays at the observer for a pinhole camera model.
!
! camera_make_ray:
!   Converts pixel coordinates (i,j) into a direction in the observer's local
!   orthonormal frame (static tetrad), then converts that photon momentum into
!   Schwarzschild coordinate covariant components p_mu.
!
! Output state:
!   y = [t, r, theta, phi, p_t, p_r, p_theta, p_phi]
! with p_mu being covariant.
! -----------------------------------------------------------------------------

MODULE camera

    use, intrinsic :: iso_fortran_env, only: wp => real64
    use kerr_newman_physics, only: M, a_bh => a, Q

    IMPLICIT NONE

    private
    public :: camera_make_ray

contains

! -----------------------------------------------------------------------------
! camera_make_ray(i, j, nx, ny, r_obs, theta, phi, fovx_deg, yt)
!
! Construct the initial ray at the observer position (r_obs, theta, phi).
!
! Steps:
!   1) Map pixel (i,j) to normalized screen coords (u,v) in [-1,1].
!   2) Create a camera-space direction using horizontal FOV.
!      Convention:
!        forward = -rhat (toward the black hole)
!        right   = +phihat
!        up      = -thetahat
!   3) Normalize the direction and build orthonormal photon momentum:
!        p^hat_t = 1, (p^hat_r, p^hat_theta, p^hat_phi) = unit direction
!   4) Convert orthonormal components to coordinate contravariant p^mu using
!      static-observer tetrad, then to covariant p_mu using the metric.
!
! Important conventions:
!   - The ray is "backward-traced": it starts at the camera and propagates inward.
!   - The stored momenta are covariant (p_t, p_r, p_theta, p_phi).
! -----------------------------------------------------------------------------

subroutine camera_make_ray(i, j, nx, ny, r_obs, theta, phi, fovx_deg, yt)
        ! Build initial ray state y at the observer for pixel (i,j),
        ! using a Kerr(-Newman) ZAMO (locally non-rotating) tetrad.
        !
        ! Output yt = [t, r, theta, phi, pt, pr, ptheta, pphi]  (COVARIANT momenta)
        !
        integer, intent(in) :: i, j, nx, ny
        real(wp), intent(in) :: r_obs, theta, phi, fovx_deg
        real(wp), intent(out) :: yt(8)

        real(wp), parameter :: pi = acos(-1.0_wp)
        real(wp) :: sth, cth, sth2
        real(wp) :: fovx, aspect, halfx, halfy
        real(wp) :: u, v, x, y
        real(wp) :: nxh, nyh, nzh, normn

        ! Kerr-Newman helpers
        real(wp) :: Sigma, Delta, B

        ! Covariant metric components at observer
        real(wp) :: g_tt, g_rr, g_thth, g_tphi, g_phiphi

        ! ZAMO quantities
        real(wp) :: alpha, omega
        real(wp) :: ut, uphi

        ! Orthonormal tetrad spatial scale factors
        real(wp) :: er, eth, ephi

        ! Photon momentum in orthonormal frame
        real(wp) :: pt_hat, pr_hat, pth_hat, pph_hat

        ! Photon momentum in coordinates (contravariant)
        real(wp) :: pt_con, pr_con, pth_con, pph_con

        ! Photon momentum in coordinates (covariant)
        real(wp) :: pt_cov, pr_cov, pth_cov, pph_cov

        ! 1) Position
        yt = 0.0_wp
        yt(1) = 0.0_wp   ! t
        yt(2) = r_obs    ! r
        yt(3) = theta    ! theta
        yt(4) = phi      ! phi

        sth  = sin(theta)
        cth  = cos(theta)
        sth2 = sth*sth
        sth2 = max(sth2, 1.0e-30_wp)   ! guard at poles

        ! 2) Kerr-Newman metric (covariant) at (r_obs, theta)
        Sigma = r_obs*r_obs + a_bh*a_bh*cth*cth
        Delta = r_obs*r_obs - 2.0_wp*M*r_obs + a_bh*a_bh + Q*Q
        B     = 2.0_wp*M*r_obs - Q*Q

        g_tt   = -(1.0_wp - B/Sigma)
        g_rr   =  Sigma / Delta
        g_thth =  Sigma
        g_tphi = -a_bh * B * sth2 / Sigma
        g_phiphi = sth2 * ( ((r_obs*r_obs + a_bh*a_bh)**2 - a_bh*a_bh*Delta*sth2) / Sigma )

        ! 3) ZAMO lapse and frame dragging angular velocity
        ! omega = -g_{tphi}/g_{phiphi}
        omega = -g_tphi / g_phiphi

        ! alpha = sqrt((g_{tphi}^2 - g_tt*g_phiphi)/g_phiphi)
        alpha = sqrt( max(0.0_wp, (g_tphi*g_tphi - g_tt*g_phiphi) / g_phiphi ) )

        ! ZAMO 4-velocity: u^mu = (1/alpha, 0, 0, omega/alpha)
        ut   = 1.0_wp / alpha
        uphi = omega / alpha

        ! 4) Pixel to camera screen coordinates u,v in [-1,1]
        u = ((real(i,wp) - 0.5_wp) / real(nx,wp)) * 2.0_wp - 1.0_wp
        v = 1.0_wp - ((real(j,wp) - 0.5_wp) / real(ny,wp)) * 2.0_wp

        fovx   = fovx_deg * pi / 180.0_wp
        aspect = real(ny,wp) / real(nx,wp)

        halfx = tan(0.5_wp * fovx)
        halfy = halfx * aspect

        x = u * halfx
        y = v * halfy

        ! 5) Build spatial direction in the local orthonormal basis (rhat, thetahat, phihat)
        ! Forward = -rhat (into scene), right = +phihat, up = -thetahat
        nxh = -1.0_wp
        nyh = -y
        nzh =  x

        normn = sqrt(nxh*nxh + nyh*nyh + nzh*nzh)
        nxh = nxh / normn
        nyh = nyh / normn
        nzh = nzh / normn

        ! 6) Null photon momentum in orthonormal frame
        pt_hat  = 1.0_wp
        pr_hat  = nxh
        pth_hat = nyh
        pph_hat = nzh

        ! 7) Spatial tetrad scale factors for ZAMO orthonormal basis:
        ! e_r^mu   = (0, 1/sqrt(g_rr),   0, 0)
        ! e_th^mu  = (0, 0, 1/sqrt(g_thth), 0)
        ! e_phi^mu = (0, 0, 0, 1/sqrt(g_phiphi))
        er   = 1.0_wp / sqrt(g_rr)
        eth  = 1.0_wp / sqrt(g_thth)
        ephi = 1.0_wp / sqrt(g_phiphi)

        ! 8) Convert tetrad components to coordinate contravariant components:
        ! p^mu = p^hat_t * u^mu + p^hat_r * e_r^mu + p^hat_th * e_th^mu + p^hat_phi * e_phi^mu
        pt_con  = pt_hat * ut
        pr_con  = pr_hat * er
        pth_con = pth_hat * eth
        pph_con = pt_hat * uphi + pph_hat * ephi

        ! 9) Convert to covariant momenta p_mu = g_{mu nu} p^nu
        pt_cov  = g_tt * pt_con + g_tphi * pph_con
        pr_cov  = g_rr * pr_con
        pth_cov = g_thth * pth_con
        pph_cov = g_tphi * pt_con + g_phiphi * pph_con

        yt(5) = pt_cov
        yt(6) = pr_cov
        yt(7) = pth_cov
        yt(8) = pph_cov

    end subroutine camera_make_ray

END MODULE camera
