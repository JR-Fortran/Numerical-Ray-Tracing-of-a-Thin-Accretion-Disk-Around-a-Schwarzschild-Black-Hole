! -----------------------------------------------------------------------------
! render_shadow.f95 / MODULE render_shadow
!
! High-level renderer:
!   For each pixel, construct an initial null ray at the observer (camera_make_ray),
!   integrate the ray backwards in affine parameter using DP5(4),
!   test whether it crosses the equatorial plane within disk radii,
!   and shade the pixel based on a simple emissivity model + redshift factor g^4.
!
! Public API:
!   render_shadow_csv : writes a CSV brightness map (float per pixel)
!   render_shadow_ppm : writes a tone-mapped colour PPM (P6)
!
! Ray termination conditions:
!   - captured : r <= 2 + eps_hor
!   - escaped  : r >= r_escape  (scaled with r_obs)
!
! Disk model:
!   Disk is the equatorial plane theta = pi/2,
!   with inner radius r_in and outer radius r_out.
! -----------------------------------------------------------------------------

MODULE render_shadow
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use DP54, only: DP54_try_step
  use schwarzschild_physics, only: rhs_schwarzschild
  use camera, only: camera_make_ray
  use disk, only: disk_crossing
  implicit none

  private
  public :: render_shadow_ppm

contains

  ! -----------------------------------------------------------------------------
  ! render_shadow_ppm
  !
  ! Render a colour image and write binary PPM (P6).
  !
  ! Differences vs CSV path:
  !   - Stores two intermediate images:
  !       Limg(i,j): luminance proxy (I_em * g^4)
  !       Timg(i,j): "colour temperature" proxy
  !   - Uses a simple temperature-to-RGB gradient (temp_gradient).
  !   - Applies global tone mapping:
  !       lum  = exposure * (L/Lmax)
  !       lum  = clamp(lum, 0..1)
  !       lum  = lum^(1/gamma_out)
  !     then scales the RGB colour by lum.
  !
  ! Disk hit refinement:
  !   After a crossing is detected, disk_refine_hit(...) bisects the step interval
  !   to refine the equatorial crossing location, reducing speckle at the disk edge.
  !
  ! Optional inputs:
  !   exposure   : overall brightness multiplier (default 1)
  !   gamma_out  : output gamma for tone mapping (default 2.2)
  ! -----------------------------------------------------------------------------

  subroutine render_shadow_ppm(nx, ny, outname, r_obs, theta_deg, phi_deg, fovx_deg, rtol, atol, T0, exposure, gamma_out)
      use, intrinsic :: iso_fortran_env, only: wp => real64, int8
      use DP54, only: DP54_try_step
      use schwarzschild_physics, only: rhs_schwarzschild
      use kerr_newman_physics, only: rhs_kerr_newman, M, a, Q
      use camera, only: camera_make_ray
      use disk, only: disk_crossing, disk_refine_hit
      use Plancks_law_mod, only: plancks_law
      implicit none

      integer, intent(in) :: nx, ny
      character(len=*), intent(in) :: outname
      real(wp), intent(in) :: r_obs, theta_deg, phi_deg, fovx_deg
      real(wp), intent(in) :: rtol, atol, T0
      real(wp), intent(in), optional :: exposure, gamma_out


      integer, parameter :: nvar = 8
      integer :: i, j, unit, step
      real(wp), parameter :: pi = acos(-1.0_wp)

      ! Cam / integrator
      real(wp) :: theta0, phi0
      real(wp) :: y0(nvar), y(nvar), y_prev(nvar)
      real(wp) :: lam, h, h_new
      logical  :: accepted
      logical  :: have_fsal
      real(wp) :: k_fsal(nvar)
      real(wp) :: yt(nvar)
      real(wp) :: k1(nvar),k2(nvar),k3(nvar),k4(nvar),k5(nvar),k6(nvar),k7(nvar)
      real(wp) :: y5(nvar), y4(nvar), sc(nvar), e(nvar)
      real(wp) :: err
      real(wp) :: lam_hit
      real(wp) :: y_hit(8)
      real(wp) :: h_acc

      ! Disk hit
      logical :: hit_plane
      real(wp) :: r_hit, phi_hit, s_hit
      real(wp) :: r_in
      real(wp), parameter :: r_out = 300.0_wp
      logical, parameter :: prograde_disk = .false.

      ! Brightness / Doppler
      real(wp) :: I_em, g, Omega
      real(wp) :: A_schw, v, gamma, cospsi
      real(wp) :: pt_cov, pphi_cov
      real(wp) :: pt_con, pphi_con
      real(wp) :: p_that, p_phihat
      real(wp) :: g_tt, g_tphi, g_phiphi
      real(wp) :: Sigma_eq, Delta, B
      real(wp) :: ut, uphi, E_inf, E_em
      real(wp) :: Omega_disk
      real(wp) :: x

      ! Stop conditions
      real(wp), parameter :: eps_hor = 1e-6_wp
      real(wp) :: r_escape
      real(wp):: r_guard
      integer, parameter :: max_steps = 200000
      logical :: captured, escaped
      real(wp) :: r_hor, disc


      ! Images (store for global tone mapping)
      real(wp), allocatable :: Limg(:,:), gimg(:,:), Temg(:,:)
      real(wp) :: Lmax, L, lum, expo, gam
      real(wp) :: rcol, gcol, bcol
      integer(int8) :: rgb(3)
      character(len=:), allocatable :: header

      ! variables for Planck's law
      real(wp), parameter :: c_light = 2.99792458e8_wp
      real(wp), parameter :: lamR = 700.0e-9_wp, lamG = 546.0e-9_wp, lamB = 435.0e-9_wp
      real(wp), parameter :: nuR = c_light/lamR, nuG = c_light/lamG, nuB = c_light/lamB
      real(wp), parameter :: I_em_max = 2.622813863194299e-4_wp  ! max of your I_em shape

      real(wp) :: T_em
      integer :: ir, ig, ib
      real(wp) :: I_R, I_G, I_B, den
      real(wp) :: gp

      theta0 = theta_deg * pi/180.0_wp
      phi0   = phi_deg   * pi/180.0_wp
      r_escape = max(200.0_wp, 2.0_wp*r_obs)

      expo = 1.0_wp
      if (present(exposure)) expo = exposure
      gam = 2.2_wp
      if (present(gamma_out)) gam = gamma_out

      allocate(Limg(nx, ny), gimg(nx, ny), Temg(nx, ny))
      Limg = 0.0_wp
      gimg = 0.0_wp
      Temg = 0.0_wp

      disc = M * M - a * a - Q * Q              ! guard (avoid NaNs if over-extremal)
      if (disc < 0.0_wp) disc = 0.0_wp
      r_hor = M + sqrt(disc)                    ! outer horizon radius r_+
      r_guard = r_hor + 0.2_wp
      r_in = r_isco_kerr(a)

      ! Parallel loop, the variables below are private to each thread to avoid race conditions, and threads overwriting each other's values.
      !$omp parallel do collapse(2) schedule(dynamic) &
      !$omp private(i, j, y0, y, y_prev, lam, h, h_new, have_fsal, k_fsal, &
      !$omp        yt, k1, k2, k3, k4, k5, k6, k7, y4, y5, sc, e, err, lam_hit, y_hit, h_acc, &
      !$omp        hit_plane, r_hit, phi_hit, s_hit, I_em, g, Omega, A_schw, v, gamma, cospsi, &
      !$omp        pt_cov, pphi_cov, pt_con, pphi_con, p_that, p_phihat, &
      !$omp        captured, escaped, L, T_em, step, &
      !$omp        g_tt, g_tphi, g_phiphi, Sigma_eq, Delta, B, ut, uphi, E_inf, E_em, Omega_disk, x)
      do j = 1, ny
        do i = 1, nx

          call camera_make_ray(i, j, nx, ny, r_obs, theta0, phi0, fovx_deg, y0)

          y = y0
          y_prev = y0
          lam = 0.0_wp
          h = 0.5_wp
          have_fsal = .false.
          k_fsal = 0.0_wp

          captured = .false.
          escaped  = .false.
          L = 0.0_wp
          g = 0.0_wp
          T_em = 0.0_wp

          do step = 1, max_steps

            if (y(2) <= r_hor + eps_hor) then
              captured = .true.; exit
            end if
            if (y(2) >= r_escape) then
              escaped = .true.; exit
            end if

            if (y(2) < r_guard) h = min(h, 1e-3_wp)

            call DP54_try_step(rhs_kerr_newman, lam, y, h, rtol, atol, &
                               have_fsal, k_fsal, yt, &
                               k1,k2,k3,k4,k5,k6,k7, &
                               y5,y4, sc,e, err, h_new, accepted)

            if (accepted) then
              h_acc = h
              lam = lam + h_acc
              h = h_new

              call disk_crossing(y_prev, y, hit_plane, r_hit, phi_hit, s_hit)
              if (hit_plane) then
                if (r_hit >= r_in .and. r_hit <= r_out) then

                  call disk_refine_hit(rhs_kerr_newman, lam-h_acc, y_prev, lam, y, 12, lam_hit, y_hit)

                  r_hit   = y_hit(2)
                  phi_hit = y_hit(4)

                  x = max(0.0_wp, 1.0_wp - sqrt(r_in / r_hit))

                  pt_cov   = y_hit(5)
                  pphi_cov = y_hit(8)

                  I_em = x / (r_hit**3)

                  ! ------------------------------------------------------------
                  ! Kerr(-Newman) redshift factor: g = E_inf / E_em
                  ! using u^mu_em = u^t (1, 0, 0, Omega)
                  ! and E_em = -p_mu u^mu = -(p_t u^t + p_phi u^phi)
                  ! ------------------------------------------------------------

                  ! Covariant metric components at equator (theta = pi/2 so sin^2=1, cos=0)
                  Sigma_eq = r_hit*r_hit
                  B        = 2.0_wp*M*r_hit - Q*Q
                  Delta    = r_hit*r_hit - 2.0_wp*M*r_hit + a*a + Q*Q

                  g_tt     = -(1.0_wp - B/Sigma_eq)
                  g_tphi   = -a * B / Sigma_eq
                  g_phiphi = ((r_hit*r_hit + a*a)**2 - a*a*Delta) / Sigma_eq

                  ! Disk angular velocity (Kerr geodesic approximation; good first step even if Q=0)
                  if (prograde_disk) then
                    Omega_disk =  1.0_wp / (r_hit**1.5_wp + a)
                  else
                    Omega_disk = -1.0_wp / (r_hit**1.5_wp - a)
                  end if
                  if (.not. prograde_disk) then
                    if (abs(r_hit**1.5_wp - a) < 1e-6_wp) Omega_disk = 0.0_wp
                  end if

                  ! Normalise emitter 4-velocity
                  ut = 1.0_wp / sqrt( -(g_tt + 2.0_wp*g_tphi*Omega_disk + g_phiphi*Omega_disk*Omega_disk) )
                  uphi = Omega_disk * ut

                  ! Conserved energy at infinity and energy measured by emitter
                  E_inf = -pt_cov
                  E_em  = -(pt_cov*ut + pphi_cov*uphi)

                  if (E_em <= 0.0_wp) then
                    g = 0.0_wp
                  else
                    g = E_inf / E_em
                  end if

                  T_em = T0 * x**0.25_wp / sqrt(r_hit)   ! ~ r^{-1/2} scaling for T

                  L = max(0.0_wp, I_em * g**4)

                  exit
                end if
              end if

              y_prev = y
            else
              h = h_new
            end if

            h = min(h, 5.0_wp)
            if (h < 1e-12_wp) then
              captured = .true.; exit
            end if

          end do

          Temg(i,j) = T_em
          gimg(i,j) = g
          Limg(i,j) = L

        end do
      end do
      !$omp end parallel do

      Lmax = maxval(Limg); if (Lmax <= 0.0_wp) Lmax = 1.0_wp

      ! Write binary PPM (P6)
      open(newunit=unit, file=outname, status="replace", access="stream", form="unformatted", action="write")
      header = "P6"//achar(10)//itoa(nx)//" "//itoa(ny)//achar(10)//"255"//achar(10)
      write(unit) header

      do j = 1, ny
        do i = 1, nx

            if (Limg(i,j) <= 0.0_wp .or. gimg(i,j) <= 1.0e-12_wp .or. Temg(i,j) <= 0.0_wp) then
                rgb = 0_int8
            else
                gp = gimg(i,j)

                I_R = gp**3 * plancks_law(nuR/gp, Temg(i,j))
                I_G = gp**3 * plancks_law(nuG/gp, Temg(i,j))
                I_B = gp**3 * plancks_law(nuB/gp, Temg(i,j))

                den = max(I_R, max(I_G, I_B))
                if (den > 0.0_wp) then
                    rcol = I_R / den
                    gcol = I_G / den
                    bcol = I_B / den
                else
                    rcol = 0.0_wp
                    gcol = 0.0_wp
                    bcol = 0.0_wp
                end if

                lum = expo * (Limg(i,j) / Lmax)
                lum = max(0.0_wp, min(1.0_wp, lum))
                lum = lum**(1.0_wp/gam)

                rcol = rcol * lum
                gcol = gcol * lum
                bcol = bcol * lum

                ir = int(255.0_wp * max(0.0_wp, min(1.0_wp, rcol)) + 0.5_wp)
                ig = int(255.0_wp * max(0.0_wp, min(1.0_wp, gcol)) + 0.5_wp)
                ib = int(255.0_wp * max(0.0_wp, min(1.0_wp, bcol)) + 0.5_wp)

                rgb(1) = int(ir, int8)
                rgb(2) = int(ig, int8)
                rgb(3) = int(ib, int8)

            end if

          write(unit) rgb

        end do
      end do

      close(unit)
      print *, "Wrote ", trim(outname)

    contains

    pure function r_isco_kerr(a) result(risco)
      real(wp), intent(in) :: a
      real(wp) :: risco
      real(wp) :: aclip, z1, z2, term

      ! Clamp just in case of tiny numerical overshoot
      aclip = max(-1.0_wp, min(1.0_wp, a))

      z1 = 1.0_wp + (1.0_wp - aclip*aclip)**(1.0_wp/3.0_wp) * &
                  ( (1.0_wp + aclip)**(1.0_wp/3.0_wp) + (1.0_wp - aclip)**(1.0_wp/3.0_wp) )
      z2 = sqrt(3.0_wp*aclip*aclip + z1*z1)

      term = sqrt((3.0_wp - z1) * (3.0_wp + z1 + 2.0_wp * z2))

      if (prograde_disk) then
        risco = 3.0_wp + z2 - term
      else
        risco = 3.0_wp + z2 + term
      end if
    end function r_isco_kerr

    pure function itoa(n) result(s)
        integer, intent(in) :: n
        character(len=:), allocatable :: s
        character(len=32) :: buf
        write(buf,'(I0)') n
        s = trim(buf)
      end function itoa

    end subroutine render_shadow_ppm

END MODULE render_shadow
