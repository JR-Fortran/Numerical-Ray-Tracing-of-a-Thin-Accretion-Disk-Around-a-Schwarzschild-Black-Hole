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
    use schwarzschild_physics, only: A_of_r

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
        ! using horizontal FOV (fovx_deg).
        !
        ! Output y = [t, r, theta, phi, pt, pr, ptheta, pphi]
        !
        ! Backwards ray tracing: rays start at camera and go inwards.
        integer, intent(in) :: i, j, nx, ny
        real(wp), intent(in) :: r_obs, theta, phi, fovx_deg
        real(wp), intent(out) :: yt(8)

        real(wp) :: A, sth
        real(wp) :: fovx, aspect, halfx, halfy
        real(wp) :: u, v, x, y
        real(wp) :: nxh, nyh, nzh, normn
        real(wp) :: pt_hat, pr_hat, pth_hat, pph_hat
        real(wp) :: pt_con, pr_con, pth_con, pph_con
        real(wp) :: pt_cov, pr_cov, pth_cov, pph_cov

        ! 1 Position the observer in Schwarzschild coordinates.
        yt = 0.0_wp
        yt(1) = 0.0_wp   ! t
        yt(2) = r_obs    ! r
        yt(3) = theta    ! theta
        yt(4) = phi      ! phi

        A = A_of_r(r_obs)
        sth = sin(theta)

        ! 2 Convert pixel (i, j) to normalized screen coordinates (u, v) in [-1, 1]
        ! i runs 1 ... nx left -> right, j runs 1 ... ny top -> bottom
        u = ((real(i,wp) - 0.5_wp) / real(nx,wp)) * 2.0_wp - 1.0_wp
        v = 1.0_wp - ((real(j,wp) - 0.5_wp) / real(ny,wp)) * 2.0_wp

        ! 3 Horizontal FOV in radiants
        fovx = fovx_deg * acos(-1.0_wp) / 180.0_wp
        aspect = real(ny,wp) / real(nx,wp)

        ! Half width of the image plane at distance 1 in camera space
        halfx = tan(0.5_wp * fovx)
        halfy = halfx * aspect

        ! 4 Camera space direction before normalisation
        ! Forward = -rhat (into scene)
        ! right = +phihat
        ! up = - thetahat
        x = u * halfx
        y = v * halfy

        ! In the local orthonormal frame (rhat, thetahat, phihat):
        ! n_hat_r corresponds to forward/back (we choose -1 for forward)
        ! n_hat_thata corresponds to up/down (we use -thatahat for "up")
        ! n_hat_phi corresponds to left/right
        !
        ! We build n_hat = forward + x * right + y * up
        nxh = -1.0_wp   ! Along r-hat (forward = -rhat)
        nyh = -y        ! Along theta-hat (because up is -thatahat)
        nzh = x         ! along phi-hat

        normn = sqrt(nxh * nxh + nyh * nyh + nzh * nzh)
        nxh = nxh / normn
        nyh = nyh / normn
        nzh = nzh / normn

        ! 5 Local orthonormal photon momentum components
        ! For a null ray in orthonormal frame p^hat_t = 1, |p^hat_spatial| = 1
        pt_hat = 1.0_wp
        pr_hat = nxh
        pth_hat = nyh
        pph_hat = nzh

        ! 6 Convert from orthonormal p^hat to coordinates contravariant p^mu
        ! tetrad for static observer:
        ! p^t = (1/sqrt(A)) * p^hat_t
        ! p^r = sqrt(A) * p^hat_r
        ! p^theta = (1/r) * p^hat_theta
        ! p^phi = (1/(r sin(θ))) * p^hat_phi
        pt_con = (1.0_wp/sqrt(A)) * pt_hat
        pr_con = sqrt(A) * pr_hat
        pth_con = (1.0_wp/r_obs) * pth_hat
        pph_con = (1.0_wp/(r_obs * sth)) * pph_hat

        ! 7 Convert to covariant p_mu = g_muv p^v
        ! Schwarzschild covariant metric:
        ! g_tt = -A, g_rr = 1/A, g_tt = r^2, g_phiphi = r^2 sin^2(theta)
        pt_cov = (-A) * pt_con
        pr_cov = (1.0_wp/A) * pr_con
        pth_cov = (r_obs * r_obs) * pth_con
        pph_cov = (r_obs * r_obs * sth * sth) * pph_con

        ! 8 Store covariant momenta in state vector
        yt(5) = pt_cov
        yt(6) = pr_cov
        yt(7) = pth_cov
        yt(8) = pph_cov

    end subroutine camera_make_ray

END MODULE camera
