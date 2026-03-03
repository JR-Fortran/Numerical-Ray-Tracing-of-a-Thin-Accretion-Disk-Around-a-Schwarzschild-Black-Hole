! -----------------------------------------------------------------------------
! utils.f95
!
! Simple utility module for reading user input. So that the main file does not
! need to be edited for changing the most common parameters.
!
! If you are generating multiple images and dont want to input the same parameters
! repeatedly, you can modify the default values in this module.
!
! Key parameters:
!   nx, ny      : output image resolution in pixels
!   r_obs       : observer radius (Schwarzschild r)
!   theta_deg   : observer polar angle in degrees
!   phi_deg     : observer azimuth in degrees
!   fovx_deg    : horizontal field of view in degrees
!   answer      : character answer from the user (y/n)
! -----------------------------------------------------------------------------

MODULE utils

    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    private
    public :: read_user_input

CONTAINS

    ! subroutine for reading user input
    subroutine read_user_input(nx, ny, r_obs, theta_deg, phi_deg, fovx_deg, T0)
        ! Declaring output variables
        integer, intent(out) :: nx, ny
        real(wp), intent(out) :: r_obs, theta_deg, phi_deg, fovx_deg, T0

        ! Declaring local variables
        character(len=1) :: answer

        ! Print a welcome message
        write(*,'(A)') "-----------------------------------------------------------------"
        write(*,'(A)') "       Welcome to Schwarzschild accretion disk renderer!"
        write(*,'(A)') "-----------------------------------------------------------------"

        ! Asking the user if they want to use the default resolution
        write(*,'(A)', ADVANCE='NO') "Do you want to use the default resolution? (4k) (y/n): "
        read(*,*) answer
        if (answer == 'y') then
            nx = 4096                   ! Default width edit these if you want to change the default resolution
            ny = 2160                   ! Default height edit these if you want to change the default resolution
        else
            write(*,'(A)', ADVANCE='NO') "Enter the number of pixels along the x-axis: "
            read(*,*) nx
            write(*,'(A)', ADVANCE='NO') "Enter the number of pixels along the y-axis: "
            read(*,*) ny
        end if

        ! Asking the user if they want to use the default camera orientation
        write(*,'(A)') "-----------------------------------------------------------------"
        write(*,'(A)', ADVANCE='NO') "Do you want to use the default camera orientation? (theta = 80, phi = 0) (y/n): "
        read(*,*) answer
        if (answer == 'y') then
            theta_deg = 80.0_wp         ! Default polar angle edit these if you want to change the default
            phi_deg = 0.0_wp            ! Default azimuthal angle edit these if you want to change the default
        else
            write(*,'(A)', ADVANCE='NO') "Enter the polar angle in degrees: "
            read(*,*) theta_deg
            write(*,'(A)', ADVANCE='NO') "Enter the azimuthal angle in degrees: "
            read(*,*) phi_deg
        end if

        ! Asking the user if they want to use the default field of view and distance
        write(*,'(A)') "-----------------------------------------------------------------"
        write(*,'(A)', ADVANCE='NO') "Do you want to use the default field of view and distance?" //&
                                        "(fovx_deg = 40, r_obs = 150) (y/n): "
        read(*,*) answer
        if (answer == 'y') then
            fovx_deg = 40.0_wp            ! Default field of view edit these if you want to change the default
            r_obs = 150.0_wp              ! Default distance edit these if you want to change the default
        else
            write(*,'(A)', ADVANCE='NO') "Enter the field of view in degrees: "
            read(*,*) fovx_deg
            write(*,'(A)', ADVANCE='NO') "Enter the distance from the observer to the center: "
            read(*,*) r_obs
        end if

        ! Asking the user if they want to use the default temperature
        write(*,'(A)') "-----------------------------------------------------------------"
        write(*,'(A)', ADVANCE='NO') "Do you want to use the default temperature? (y/n): "
        read(*,*) answer
        if (answer == 'y') then
            T0 = 9000.0_wp            ! Default temperature edit this if you want to change the default
        else
            write(*,'(A)', ADVANCE='NO') "Enter the temperature: "
            read(*,*) T0
        end if

        ! Printing the settings
        write(*,'(A)') "================================================================="
        write(*,*) "Settings:"
        write(*,'(A)') "-----------------------------------------------------------------"
        write(*,*) "  nx = ", nx
        write(*,*) "  ny = ", ny
        write(*,'(A)') "-----------------------------------------------------------------"
        write(*,*) "  theta_deg = ", theta_deg
        write(*,*) "  phi_deg = ", phi_deg
        write(*,'(A)') "-----------------------------------------------------------------"
        write(*,*) "  fovx_deg = ", fovx_deg
        write(*,*) "  r_obs = ", r_obs
        write(*,'(A)') "-----------------------------------------------------------------"
        write(*,*) "  T0 = ", T0
        write(*,'(A)') "-----------------------------------------------------------------"

        ! Warn the user about the resolution
        if (nx > 4096 .or. ny > 2160) then
            write(*,'(A)') "-----------------------------------------------------------------"
            write(*,*) "Warning: it looks like you're using a resolution higher than 4k."
            write(*,*) "This may take a while to render."
            write(*,'(A)') "-----------------------------------------------------------------"
        end if

        ! Asking the user if they want to proceed with these settings
        write(*,'(A)', ADVANCE='NO') "Do you want to proceed with these settings? (y/n): "
        read(*,*) answer
        if (answer == 'y') then
            write(*,*) "Proceeding with settings..."
        else
            write(*,*) "Exiting..."
            stop
        end if

    end subroutine read_user_input

END MODULE utils
