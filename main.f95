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

PROGRAM main
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use render_shadow, only: render_shadow_ppm
  use utils, only: read_user_input
  implicit none

  integer :: nx, ny
  !character(len=*), parameter :: filename = "disk.csv"
  character(len=*), parameter :: filename = "disk.ppm"

  real(wp) :: r_obs, theta_deg, phi_deg, fovx_deg, T0
  real(wp) :: rtol, atol


  call read_user_input(nx, ny, r_obs, theta_deg, phi_deg, fovx_deg, T0)

  rtol = 1e-9_wp
  atol = 1e-12_wp

  !call render_shadow_csv(nx, ny, filename, r_obs, theta_deg, phi_deg, fovx_deg, rtol, atol)
  call render_shadow_ppm(nx, ny, filename, r_obs, theta_deg, phi_deg, fovx_deg, rtol, atol, T0, &
                         exposure=2.0_wp, gamma_out=2.2_wp)
END PROGRAM main
