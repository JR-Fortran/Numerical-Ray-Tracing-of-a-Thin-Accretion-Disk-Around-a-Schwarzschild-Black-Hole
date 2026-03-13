! -----------------------------------------------------------------------------
! main.f95
!
! Entry point. Configures camera + integrator tolerances and calls the renderer.
!
! Key parameters:
!   nx, ny      : output image resolution in pixels
!   r_obs       : observer radius (Schwarzschild r)
!   theta_deg   : observer polar angle in degrees
!   phi_deg     : observer azimuth in degrees
!   fovx_deg    : horizontal field of view in degrees
!   rtol, atol  : scalar relative/absolute tolerances for adaptive integrator
!
! Output:
!   disk.ppm    : rendered image (P6 binary PPM)
! -----------------------------------------------------------------------------

program main
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use render_shadow, only: render_shadow_ppm
    use utils, only: ui, settings_type
    use kerr_newman_physics, only: a, Q
    implicit none

    type(settings_type) :: settings
    real(wp) :: rtol, atol
    logical :: ask_about_os = .true.
    logical :: reset_settings = .true.

    10 continue

    call ui(settings, ask_about_os, reset_settings)

    ask_about_os = .false.
    reset_settings = .false.

    a = settings%spin_a
    Q = settings%charge_q

    rtol = 1e-9_wp
    atol = 1e-12_wp

    call render_shadow_ppm(settings%nx, settings%ny, settings%filename, settings%observer_radius, &
                         settings%observer_theta, settings%observer_phi, settings%fov, rtol, atol, &
                         settings%temperature_scale, settings%prograde, exposure=2.0_wp, gamma_out=2.2_wp)

    goto 10
end program main
