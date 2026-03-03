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
  public :: render_shadow_csv, render_shadow_ppm

contains

! -----------------------------------------------------------------------------
! render_shadow_csv
!
! Render an image by ray tracing and write per-pixel brightness to CSV.
!
! For each pixel:
!   1) camera_make_ray(...) gives initial state y0 at the observer.
!   2) Integrate dy/dλ = rhs_schwarzschild with DP54_try_step.
!   3) Detect equatorial-plane crossing between accepted steps (disk_crossing).
!   4) If crossing radius r_hit is within [r_in, r_out], compute brightness:
!        I_em : local emissivity
!        g    : frequency-shift factor (Doppler + gravitational)
!        pix  : observed brightness proxy = I_em * g^4
!
! Notes on key variables:
!   lam        : affine parameter λ (integration parameter)
!   h, h_new   : current and suggested next step size in λ
!   y_prev     : last accepted state (used for plane-crossing interpolation)
!   hit_plane  : true if theta crossed pi/2 between y_prev and y
!   r_hit      : interpolated crossing radius
!   phi_hit    : interpolated crossing azimuth (with unwrap)
!
! Doppler / redshift terms:
!   A      = 1 - 2/r_hit
!   Omega  = -1 / r^(3/2)      (Keplerian angular velocity; sign sets rotation sense)
!   v      = (r*Omega)/sqrt(A) (orbital speed measured by static observer)
!   gamma  = 1/sqrt(1-v^2)
!   cospsi : photon direction cosine relative to disk motion in local static frame
!   g      = sqrt(A) / (gamma * (1 - v*cospsi))
!
! Output:
!   CSV with ny rows and nx columns (comma-separated floats).
! -----------------------------------------------------------------------------

  subroutine render_shadow_csv(nx, ny, outname, r_obs, theta_deg, phi_deg, fovx_deg, rtol, atol)
    integer, intent(in) :: nx, ny
    character(len=*), intent(in) :: outname
    real(wp), intent(in) :: r_obs, theta_deg, phi_deg, fovx_deg
    real(wp), intent(in) :: rtol, atol

    integer, parameter :: nvar = 8
    integer :: i, j, unit, step
    real(wp), parameter :: pi = acos(-1.0_wp)

    ! Camera angles (radians)
    real(wp) :: theta0, phi0

    ! Integrator
    real(wp) :: lam, h, h_new
    logical :: accepted, have_fsal

    ! Ray state and work arrays
    real(wp) :: y(nvar), y0(nvar), y_prev(nvar), y5(nvar), y4(nvar), yt(nvar)
    real(wp) :: k_fsal(nvar), k1(nvar),k2(nvar),k3(nvar),k4(nvar),k5(nvar),k6(nvar),k7(nvar)
    real(wp) :: sc(nvar), e(nvar), err

    ! Disk intersection
    logical :: hit_plane, hit_disk
    real(wp) :: r_hit, phi_hit, s_hit
    real(wp), parameter :: r_in  = 6.0_wp
    real(wp), parameter :: r_out = 100.0_wp

    ! Brightness
    real(wp) :: pix, I_em, g, nu_em, E_inf
    real(wp) :: u_t, u_phi, Omega

    ! Stop conditions
    real(wp), parameter :: eps_hor = 1e-6_wp
    real(wp) :: r_escape
    real(wp), parameter :: r_guard = 2.2_wp
    integer, parameter :: max_steps = 200000
    logical :: captured, escaped

    real(wp) :: A, v, gamma, cospsi
    real(wp) :: pt_cov, pphi_cov
    real(wp) :: pt_con, pphi_con
    real(wp) :: p_that, p_phihat

    theta0 = theta_deg * pi/180.0_wp
    phi0   = phi_deg   * pi/180.0_wp

    ! Make escape radius scale with observer distance (helps when you move r_obs around)
    r_escape = max(200.0_wp, 2.0_wp*r_obs)

    open(newunit=unit, file=outname, status="replace", action="write")

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
        escaped = .false.
        hit_disk = .false.
        pix = 0.0_wp

        do step = 1, max_steps

          if (y(2) <= 2.0_wp + eps_hor) then
            captured = .true.
            exit
          end if
          if (y(2) >= r_escape) then
            escaped = .true.
            exit
          end if

          if (y(2) < r_guard) h = min(h, 1e-3_wp)

          call DP54_try_step(rhs_schwarzschild, lam, y, h, rtol, atol, &
                             have_fsal, k_fsal, yt, &
                             k1,k2,k3,k4,k5,k6,k7, &
                             y5,y4, sc,e, err, h_new, accepted)

          if (accepted) then
            lam = lam + h
            h = h_new

            call disk_crossing(y_prev, y, hit_plane, r_hit, phi_hit, s_hit)
            if (hit_plane) then
              if (r_hit >= r_in .and. r_hit <= r_out) then
                ! emitted brightness model (placeholder)
                I_em = (1.0_wp/r_hit**3) * (1.0_wp - sqrt(6.0_wp/r_hit))

                A = 1.0_wp - 2.0_wp/r_hit

                pt_cov   = y(5)
                pphi_cov = y(8)

                ! contravariant components
                pt_con   = (-1.0_wp/A) * pt_cov
                pphi_con = (1.0_wp/(r_hit*r_hit)) * pphi_cov

                ! static orthonormal frame components
                p_that   = sqrt(A) * pt_con
                p_phihat = r_hit * pphi_con

                cospsi = p_phihat / p_that

                Omega = -1.0_wp / (r_hit**1.5_wp)
                v = (r_hit * Omega) / sqrt(A)

                if (v >= 1.0_wp) v = 0.999999_wp
                gamma = 1.0_wp / sqrt(1.0_wp - v*v)

                g = sqrt(A) / (gamma * (1.0_wp - v*cospsi))

                pix = I_em * g**4
                hit_disk = .true.
                exit
              end if
            end if

            y_prev = y
          else
            h = h_new
          end if

          h = min(h, 5.0_wp)
          if (h < 1e-12_wp) then
            captured = .true.
            exit
          end if

        end do

        ! Write brightness (float). No hit => 0.
        if (i < nx) then
          write(unit, "(ES16.8,A)", advance="no") pix, ","
        else
          write(unit, "(ES16.8)") pix
        end if

      end do
    end do

    close(unit)
    print *, "Wrote ", trim(outname)

  end subroutine render_shadow_csv

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

  subroutine render_shadow_ppm(nx, ny, outname, r_obs, theta_deg, phi_deg, fovx_deg, rtol, atol, exposure, gamma_out)
      use, intrinsic :: iso_fortran_env, only: wp => real64, int8
      use DP54, only: DP54_try_step
      use schwarzschild_physics, only: rhs_schwarzschild
      use camera, only: camera_make_ray
      use disk, only: disk_crossing, disk_refine_hit
      implicit none

      integer, intent(in) :: nx, ny
      character(len=*), intent(in) :: outname
      real(wp), intent(in) :: r_obs, theta_deg, phi_deg, fovx_deg
      real(wp), intent(in) :: rtol, atol
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
      real(wp), parameter :: r_in  = 6.0_wp
      real(wp), parameter :: r_out = 300.0_wp

      ! Brightness / Doppler
      real(wp) :: I_em, g, Omega
      real(wp) :: A, v, gamma, cospsi
      real(wp) :: pt_cov, pphi_cov
      real(wp) :: pt_con, pphi_con
      real(wp) :: p_that, p_phihat

      ! Stop conditions
      real(wp), parameter :: eps_hor = 1e-6_wp
      real(wp) :: r_escape
      real(wp), parameter :: r_guard = 2.2_wp
      integer, parameter :: max_steps = 200000
      logical :: captured, escaped

      ! Images (store for global tone mapping)
      real(wp), allocatable :: Limg(:,:), Timg(:,:)
      real(wp) :: Lmax, Tmax, L, T, tnorm, lum, expo, gam
      real(wp) :: rcol, gcol, bcol
      integer :: ir, ig, ib
      integer(int8) :: rgb(3)
      character(len=:), allocatable :: header

      theta0 = theta_deg * pi/180.0_wp
      phi0   = phi_deg   * pi/180.0_wp
      r_escape = max(200.0_wp, 2.0_wp*r_obs)

      expo = 1.0_wp
      if (present(exposure)) expo = exposure
      gam = 2.2_wp
      if (present(gamma_out)) gam = gamma_out

      allocate(Limg(nx,ny), Timg(nx,ny))
      Limg = 0.0_wp
      Timg = 0.0_wp

      ! Parallel loop, the variables below are private to each thread to avoid race conditions, and threads overwriting each other's values.
      !$omp parallel do collapse(2) schedule(dynamic) &
      !$omp private(i, j, y0, y, y_prev, lam, h, h_new, have_fsal, k_fsal, &
      !$omp        yt, k1, k2, k3, k4, k5, k6, k7, y4, y5, sc, e, err, lam_hit, y_hit, h_acc, &
      !$omp        hit_plane, r_hit, phi_hit, s_hit, I_em, g, Omega, A, v, gamma, cospsi, &
      !$omp        pt_cov, pphi_cov, pt_con, pphi_con, p_that, p_phihat, &
      !$omp        captured, escaped, L, T, step)
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
          T = 0.0_wp

          do step = 1, max_steps

            if (y(2) <= 2.0_wp + eps_hor) then
              captured = .true.; exit
            end if
            if (y(2) >= r_escape) then
              escaped = .true.; exit
            end if

            if (y(2) < r_guard) h = min(h, 1e-3_wp)

            call DP54_try_step(rhs_schwarzschild, lam, y, h, rtol, atol, &
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

                  call disk_refine_hit(rhs_schwarzschild, lam-h_acc, y_prev, lam, y, 12, lam_hit, y_hit)

                  r_hit   = y_hit(2)
                  phi_hit = y_hit(4)

                  pt_cov   = y_hit(5)
                  pphi_cov = y_hit(8)

                  I_em = (1.0_wp/r_hit**3) * (1.0_wp - sqrt(6.0_wp/r_hit))

                  A = 1.0_wp - 2.0_wp/r_hit

                  pt_con   = (-1.0_wp/A) * pt_cov
                  pphi_con = (1.0_wp/(r_hit*r_hit)) * pphi_cov

                  p_that   = sqrt(A) * pt_con
                  p_phihat = r_hit * pphi_con

                  cospsi = p_phihat / p_that

                  Omega = -1.0_wp / (r_hit**1.5_wp)
                  v = (r_hit * Omega) / sqrt(A)
                  if (v >= 1.0_wp) v = 0.999999_wp
                  gamma = 1.0_wp / sqrt(1.0_wp - v*v)

                  g = sqrt(A) / (gamma * (1.0_wp - v*cospsi))

                  L = max(0.0_wp, I_em * g**4)

                  ! local emission temperature proxy
                  T = max(0.0_wp, I_em)**0.25_wp

                  ! observed colour temperature (spectral shift)
                  T = g * T

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

          Limg(i,j) = L
          Timg(i,j) = T

        end do
      end do
      !$omp end parallel do

      Lmax = maxval(Limg); if (Lmax <= 0.0_wp) Lmax = 1.0_wp
      Tmax = maxval(Timg); if (Tmax <= 0.0_wp) Tmax = 1.0_wp

      ! Write binary PPM (P6)
      open(newunit=unit, file=outname, status="replace", access="stream", form="unformatted", action="write")
      header = "P6"//achar(10)//itoa(nx)//" "//itoa(ny)//achar(10)//"255"//achar(10)
      write(unit) header

      do j = 1, ny
        do i = 1, nx

          if (Limg(i,j) <= 0.0_wp) then
            rgb = 0_int8
          else
            tnorm = Timg(i,j) / Tmax
            tnorm = max(0.0_wp, min(1.0_wp, tnorm))

            call temp_gradient(tnorm, rcol, gcol, bcol)

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

      pure function itoa(n) result(s)
        integer, intent(in) :: n
        character(len=:), allocatable :: s
        character(len=32) :: buf
        write(buf,'(I0)') n
        s = trim(buf)
      end function itoa

      pure subroutine temp_gradient(t, r, g, b)
        real(wp), intent(in)  :: t
        real(wp), intent(out) :: r, g, b
        real(wp) :: u

        ! Cool (redshifted, outer) -> hot (blueshifted, inner)
        ! t=0: deep red, t=0.5: orange/yellow, t=1: blue-white
        if (t < 0.33_wp) then
          ! deep red -> orange
          u = t / 0.33_wp
          r = 0.6_wp + 0.4_wp*u
          g = 0.05_wp*u
          b = 0.0_wp
        else if (t < 0.66_wp) then
          ! orange -> yellow-white
          u = (t - 0.33_wp) / 0.33_wp
          r = 1.0_wp
          g = 0.05_wp + 0.75_wp*u
          b = 0.2_wp*u
        else
          ! yellow-white -> blue
          u = (t - 0.66_wp) / 0.34_wp
          r = (1.0_wp - u)*1.0_wp + u*0.3_wp
          g = (1.0_wp - u)*0.8_wp + u*0.5_wp
          b = (1.0_wp - u)*0.2_wp + u*1.0_wp
        end if
      end subroutine temp_gradient

    end subroutine render_shadow_ppm

END MODULE render_shadow
