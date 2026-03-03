MODULE Plancks_law_mod
  use, intrinsic :: iso_fortran_env, only: wp => real64

  implicit none

  private
  public :: plancks_law

  ! Physics constants (SI units)
  real(wp), parameter :: h = 6.62607015e-34_wp
  real(wp), parameter :: k_B = 1.380649e-23_wp
  real(wp), parameter :: c = 2.99792458e8_wp

CONTAINS

  pure function plancks_law(nu, T) result(B)
    real(wp), intent(in) :: nu, T               ! Frequency and temperature
    real(wp) :: B                               ! B(nu, T)
    real(wp) :: x                               ! Intermediate variable

    if (nu <= 0.0_wp .or. T <= 0.0_wp) then
      B = 0.0_wp
      return
    end if

    x = (h * nu) / (k_B * T)
    B = (2.0_wp * h * (nu * nu * nu) / (c * c)) * (1 / (exp(x) - 1.0_wp))

  end function plancks_law

END MODULE Plancks_law_mod
