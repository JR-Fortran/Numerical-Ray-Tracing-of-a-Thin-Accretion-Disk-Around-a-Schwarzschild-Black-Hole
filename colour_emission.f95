MODULE colour_emission
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use Plancks_law_mod, only: plancks_law

    IMPLICIT NONE

    private

    public :: disk_emission_profile
    public :: observed_bolometric_luminance
    public :: observed_rgb_from_temperature
    public :: tone_map_rgb_to_int8

    ! Light
    real(wp), parameter :: c_light = 2.99792458e8_wp                 ! Speed of light in m/s
    real(wp), parameter :: lamR = 700.0e-9_wp                        ! Red 700 nm
    real(wp), parameter :: lamG = 546.0e-9_wp                        ! Green 546 nm
    real(wp), parameter :: lamB = 435.0e-9_wp                        ! Blue 435 nm

    ! Frequencies
    real(wp), parameter :: nuR = c_light / lamR
    real(wp), parameter :: nuG = c_light / lamG
    real(wp), parameter :: nuB = c_light / lamB

CONTAINS

    pure subroutine disk_emission_profile(r, r_in, T0, I_em, T_em)
        real(wp), intent(in) :: r, r_in, T0                         ! r is radius where the ray hits the disk, r_in is the inner edge of the disk, T0 is the temperature scale
        real(wp), intent(out) :: I_em, T_em                         ! I_em is emitted intesity proxy, T_em is the emitted temperature proxy

        real(wp) :: x                                               ! helper variable

        ! If the radius in nonsense, then emit nothing
        if (r <= 0.0_wp .or. r_in <= 0.0_wp .or. r < r_in) then
            I_em = 0.0_wp
            T_em = 0.0_wp
            return
        end if

        ! Inner edge factor, emission only turns on outside the inner edge
        x = max(0.0_wp, 1.0_wp - sqrt(r_in / r))

        ! Emitted intenisty proxy, brightness falls off as 1/r^3, modulated by inner edge factor
        I_em = x / (r * r * r)

        ! Temperature proxy,
        if (x > 0.0_wp) then
            T_em = T0 * (x / (r * r * r))**0.25_wp
        else
            T_em = 0.0_wp
        end if
    end subroutine disk_emission_profile

    pure function observed_bolometric_luminance(I_em, g) result(L)
        real(wp), intent(in) :: I_em, g                            ! I_em emitted brightness proxy, g redshift factor
        real(wp) :: L                                              ! L observed brightness proxy

        if (I_em <= 0.0_wp .or. g <= 0.0_wp) then
            L = 0.0_wp
        else
            ! Bolometric luminance scales as I_em * g^4
            ! one factor of g from photon energy
            ! one factor from arrival rate / time dilation
            ! two factors for frequency-bin / solid-angle
            L = I_em * g**4
        end if
    end function observed_bolometric_luminance

    pure subroutine observed_rgb_from_temperature(g, T_em, rcol, gcol, bcol)
        real(wp), intent(in) :: g, T_em                            ! g redshift factor, T_em temperature proxy
        real(wp), intent(out) :: rcol, gcol, bcol                  ! rcol, gcol, bcol RGB colour components (not yet pixel values)

        real(wp) :: I_R, I_G, I_B, den

        ! If ray contributes nothing useful, return to black.
        if (g <= 0.0_wp .or. T_em <= 0.0_wp) then
            rcol = 0.0_wp
            gcol = 0.0_wp
            bcol = 0.0_wp
            return
        end if

        ! For each observed channel frequency, look up the emitted spectrum nu/g, multply by g^3.
        I_R = (g * g * g) * plancks_law(nuR / g, T_em)
        I_G = (g * g * g) * plancks_law(nuG / g, T_em)
        I_B = (g * g * g) * plancks_law(nuB / g, T_em)

        ! Normalize colour ratios, this finds the strongest channel.
        den = max(I_R, max(I_G, I_B))

        ! This makes it so that the strongest channel contributes 1.0, and the others are fractions of it.
        if (den > 0.0_wp) then
            rcol = I_R / den
            gcol = I_G / den
            bcol = I_B / den
        else
            rcol = 0.0_wp
            gcol = 0.0_wp
            bcol = 0.0_wp
        end if
    end subroutine observed_rgb_from_temperature

    pure subroutine tone_map_rgb_to_int8(L, Lmax, exposure, gamma_out, &
                                        rcol_in, gcol_in, bcol_in, rgb)
        real(wp), intent(in) :: L, Lmax, exposure, gamma_out              ! L brightness of pixel, Lmax brightest pixel in image, exposure user scale, gamma_out display gamma
        real(wp), intent(in) :: rcol_in, gcol_in, bcol_in
        integer, intent(out) :: rgb(3)                             ! rgb(3) final 8-bit RGB triplet

        real(wp) :: lum
        real(wp) :: rcol, gcol, bcol
        integer :: ir, ig, ib

        ! If there is no light make the pixel black
        if (L <= 0.0_wp .or. Lmax <= 0.0_wp) then
            rgb = 0
            return
        end if

        ! Normalise image brightness by scaling the brightness relative to the brightest pixel.
        lum = exposure * (L / Lmax)
        ! Clamp, prevents values from going below 0 or above 1.
        lum = max(0.0_wp, min(1.0_wp, lum))
        ! Gamma correction, makes the image appear brighter, without this the image may appear too dark.
        lum = lum**(1.0_wp / gamma_out)

        ! Applying brightness to RGB channels, rcol_in, gcol_in, bcol_in carry hue, lum carries brightness. Multiplying them combines hue and brightness.
        rcol = max(0.0_wp, min(1.0_wp, rcol_in * lum))
        gcol = max(0.0_wp, min(1.0_wp, gcol_in * lum))
        bcol = max(0.0_wp, min(1.0_wp, bcol_in * lum))

        ! Convert to 0-255 range and round to nearest integer.
        ir = int(255.0_wp * rcol + 0.5_wp)
        ig = int(255.0_wp * gcol + 0.5_wp)
        ib = int(255.0_wp * bcol + 0.5_wp)

        ! Store the final RGB values as integers in the range 0-255.
        rgb(1) = ir
        rgb(2) = ig
        rgb(3) = ib
    end subroutine tone_map_rgb_to_int8

END MODULE colour_emission
