module utils

    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    private

    public :: ui
    public :: settings_type
    public :: unix
    public :: progress_bar

    logical :: unix = .true.

    ! For progress bar
    logical :: calculate_cycels_frac = .false.
    real(wp) :: frac

    type :: settings_type
        integer :: nx = 4096
        integer :: ny = 2160
        real(wp) :: observer_radius   = 150.0_wp
        real(wp) :: observer_theta    = 85.0_wp
        real(wp) :: observer_phi      = 0.0_wp
        real(wp) :: fov               = 40.0_wp
        real(wp) :: spin_a            = 0.0_wp
        real(wp) :: charge_q          = 0.0_wp
        real(wp) :: temperature_scale = 50000.0_wp
        logical  :: prograde          = .false.
        character(len=256) :: filename = "BlackHole_image.ppm"
    end type settings_type

contains

    subroutine ui(settings, ask_about_os, reset_settings)
        type(settings_type), intent(inout) :: settings
        logical, intent(in) :: ask_about_os, reset_settings
        logical :: continue_loop
        integer :: menu_choice

        if (ask_about_os) then
            call choose_terminal_type()
        end if
        if (reset_settings) then
            call set_default_settings(settings)
        end if
        continue_loop = .true.
        do while (continue_loop)
            call main_menu(settings, menu_choice)
            call handle_main_menu_choice(settings, menu_choice, continue_loop)
        end do
    end subroutine ui


    subroutine choose_terminal_type()
        integer :: os_choice, ios
        logical :: valid_choice

        valid_choice = .false.

        do while (.not. valid_choice)
            call clear_screen()
            write(*,'(A)') "========================================"
            write(*,'(A)') "           Terminal Setup"
            write(*,'(A)') "========================================"
            write(*,'(A)') ""
            write(*,'(A)') "Choose what OS you are using:"
            write(*,'(A)') "1. macOS / Linux"
            write(*,'(A)') "2. Windows"
            write(*,'(A)') ""
            write(*,'(A)', advance='no') "Select an option: "

            read(*,*,iostat=ios) os_choice
            if (ios /= 0) then
                call flush_input_line()
                write(*,'(A)') ""
                write(*,'(A)') "Invalid input. Please enter 1 or 2."
                call wait_for_enter()
                cycle
            end if

            select case (os_choice)
            case (1)
                unix = .true.
                valid_choice = .true.
            case (2)
                unix = .false.
                valid_choice = .true.
            case default
                write(*,'(A)') ""
                write(*,'(A)') "Invalid option. Please enter 1 or 2."
                call wait_for_enter()
            end select
        end do
    end subroutine choose_terminal_type


    subroutine set_default_settings(settings)
        type(settings_type), intent(inout) :: settings

        settings%nx                = 4096
        settings%ny                = 2160
        settings%observer_radius   = 150.0_wp
        settings%observer_theta    = 85.0_wp
        settings%observer_phi      = 0.0_wp
        settings%fov               = 40.0_wp
        settings%spin_a            = 0.0_wp
        settings%charge_q          = 0.0_wp
        settings%temperature_scale = 50000.0_wp
        settings%prograde          = .false.
        settings%filename          = "BlackHole_image.ppm"
    end subroutine set_default_settings


    subroutine clear_screen()
        if (unix) then
            call execute_command_line("clear")
        else
            call execute_command_line("cls")
        end if
    end subroutine clear_screen


    subroutine main_menu(settings, menu_choice)
        type(settings_type), intent(in) :: settings
        integer, intent(out) :: menu_choice
        integer :: ios

        do
            call clear_screen()

            write(*,'(A)') "========================================"
            write(*,'(A)') "         Black Hole Ray Tracer"
            write(*,'(A)') "========================================"
            write(*,'(A)') ""

            call print_current_settings(settings)

            write(*,'(A)') ""
            write(*,'(A)') "Menu:"
            write(*,'(A)') "1. Change resolution"
            write(*,'(A)') "2. Change filename"
            write(*,'(A)') "3. Change observer settings"
            write(*,'(A)') "4. Change disk orientation"
            write(*,'(A)') "5. Change black hole parameters"
            write(*,'(A)') "6. Change temperature scale"
            write(*,'(A)') "7. Review settings"
            write(*,'(A)') "8. Start rendering"
            write(*,'(A)') "9. Reset to defaults"
            write(*,'(A)') "10. Quit"
            write(*,'(A)') ""
            write(*,'(A)', advance='no') "Enter choice: "

            read(*,*,iostat=ios) menu_choice
            if (ios == 0) exit

            call flush_input_line()
            write(*,'(A)') ""
            write(*,'(A)') "Invalid input. Please enter a number."
            call wait_for_enter()
        end do
    end subroutine main_menu


    subroutine print_current_settings(settings)
        type(settings_type), intent(in) :: settings

        write(*,'(A)') "Current settings:"
        write(*,'(A, I0, A, I0)')  "Resolution          : ", settings%nx, " x ", settings%ny
        write(*,'(A, A)')          "Filename            : ", trim(settings%filename)
        write(*,'(A, F10.3)')      "Observer radius     : ", settings%observer_radius
        write(*,'(A, F10.3, A)')   "Observer theta      : ", settings%observer_theta, " deg"
        write(*,'(A, F10.3, A)')   "Observer phi        : ", settings%observer_phi, " deg"
        write(*,'(A, F10.3, A)')   "FOV                 : ", settings%fov, " deg"

        if (settings%prograde) then
            write(*,'(A)') "Disk orientation    : Prograde"
        else
            write(*,'(A)') "Disk orientation    : Retrograde"
        end if

        write(*,'(A, F10.3)')      "Spin a              : ", settings%spin_a
        write(*,'(A, F10.3)')      "Charge Q            : ", settings%charge_q
        write(*,'(A, F10.3, A)')   "Temperature scale   : ", settings%temperature_scale, " K"
    end subroutine print_current_settings


    subroutine handle_main_menu_choice(settings, menu_choice, continue_loop)
        type(settings_type), intent(inout) :: settings
        integer, intent(in) :: menu_choice
        logical, intent(inout) :: continue_loop

        select case (menu_choice)
        case (1)
            call change_resolution(settings)
        case (2)
            call change_filename(settings)
        case (3)
            call change_camera_settings(settings)
        case (4)
            call change_disk_orientation(settings)
        case (5)
            call change_black_hole_parameters(settings)
        case (6)
            call change_temperature_scale(settings)
        case (7)
            call review_settings(settings)
            call wait_for_enter()
        case (8)
            if (confirm_render(settings)) then
                continue_loop = .false.
            end if
        case (9)
            call set_default_settings(settings)
            write(*,'(A)') ""
            write(*,'(A)') "Settings restored to defaults."
            call wait_for_enter()
        case (10)
            stop
        case default
            write(*,'(A)') ""
            write(*,'(A)') "Invalid menu choice."
            call wait_for_enter()
        end select
    end subroutine handle_main_menu_choice


    subroutine change_resolution(settings)
        type(settings_type), intent(inout) :: settings
        integer :: nx_new, ny_new, ios
        character(len=1) :: choice

        call clear_screen()
        write(*,'(A)') "========================================"
        write(*,'(A)') "          Change Resolution"
        write(*,'(A)') "========================================"
        write(*,'(A)') ""

        write(*,'(A)', advance='no') "Enter width  : "
        read(*,*,iostat=ios) nx_new
        if (ios /= 0) then
            call flush_input_line()
            write(*,'(A)') "Invalid width."
            call wait_for_enter()
            return
        end if

        write(*,'(A)', advance='no') "Enter height : "
        read(*,*,iostat=ios) ny_new
        if (ios /= 0) then
            call flush_input_line()
            write(*,'(A)') "Invalid height."
            call wait_for_enter()
            return
        end if

        if (nx_new <= 0 .or. ny_new <= 0) then
            write(*,'(A)') ""
            write(*,'(A)') "Resolution must be positive."
            call wait_for_enter()
            return
        end if

        if (nx_new > 4096 .or. ny_new > 2160) then
            write(*,'(A)') ""
            write(*,'(A)') "Warning: this is larger than 4K and may take a long time."
            write(*,'(A)', advance='no') "Continue anyway? (y/n): "
            read(*,'(A)') choice
            if (choice /= 'y' .and. choice /= 'Y') return
        end if

        settings%nx = nx_new
        settings%ny = ny_new
    end subroutine change_resolution


    subroutine change_filename(settings)
        type(settings_type), intent(inout) :: settings
        character(len=256) :: filename_new

        call clear_screen()
        write(*,'(A)') "========================================"
        write(*,'(A)') "           Change Filename"
        write(*,'(A)') "========================================"
        write(*,'(A)') ""
        write(*,'(A)', advance='no') "Enter output filename: "
        read(*,'(A)') filename_new

        if (len_trim(filename_new) > 0) then
            settings%filename = trim(filename_new)
        else
            write(*,'(A)') ""
            write(*,'(A)') "Filename unchanged."
            call wait_for_enter()
        end if
    end subroutine change_filename


    subroutine change_camera_settings(settings)
        type(settings_type), intent(inout) :: settings
        integer :: choice, ios
        logical :: done

        done = .false.
        do while (.not. done)
            call clear_screen()
            write(*,'(A)') "========================================"
            write(*,'(A)') "       Change observer Settings"
            write(*,'(A)') "========================================"
            write(*,'(A)') ""
            write(*,'(A)') "1. Change observer radius"
            write(*,'(A)') "2. Change observer theta"
            write(*,'(A)') "3. Change observer phi"
            write(*,'(A)') "4. Change FOV"
            write(*,'(A)') "5. Back to main menu"
            write(*,'(A)') ""
            write(*,'(A)', advance='no') "Enter choice: "

            read(*,*,iostat=ios) choice
            if (ios /= 0) then
                call flush_input_line()
                cycle
            end if

            select case (choice)
            case (1)
                call set_positive_real("Observer radius", settings%observer_radius)
            case (2)
                call set_real_in_range("Observer theta (deg)", settings%observer_theta, 0.0_wp, 180.0_wp)
            case (3)
                call set_any_real("Observer phi (deg)", settings%observer_phi)
            case (4)
                call set_positive_real("FOV (deg)", settings%fov)
            case (5)
                done = .true.
            case default
                write(*,'(A)') ""
                write(*,'(A)') "Invalid option."
                call wait_for_enter()
            end select
        end do
    end subroutine change_camera_settings


    subroutine change_disk_orientation(settings)
        type(settings_type), intent(inout) :: settings
        integer :: choice, ios

        call clear_screen()
        write(*,'(A)') "========================================"
        write(*,'(A)') "       Change Disk Orientation"
        write(*,'(A)') "========================================"
        write(*,'(A)') ""
        write(*,'(A)') "1. Prograde"
        write(*,'(A)') "2. Retrograde"
        write(*,'(A)') ""
        write(*,'(A)', advance='no') "Enter choice: "

        read(*,*,iostat=ios) choice
        if (ios /= 0) then
            call flush_input_line()
            write(*,'(A)') "Invalid input."
            call wait_for_enter()
            return
        end if

        select case (choice)
        case (1)
            settings%prograde = .true.
        case (2)
            settings%prograde = .false.
        case default
            write(*,'(A)') "Invalid option."
            call wait_for_enter()
        end select
    end subroutine change_disk_orientation


    subroutine change_black_hole_parameters(settings)
        type(settings_type), intent(inout) :: settings
        integer :: choice, ios
        logical :: done

        done = .false.
        do while (.not. done)
            call clear_screen()
            write(*,'(A)') "========================================"
            write(*,'(A)') "     Change Black Hole Parameters"
            write(*,'(A)') "========================================"
            write(*,'(A)') ""
            write(*,'(A)') "Note: a = 0 and q = 0 gives Schwarzschild."
            write(*,'(A)') ""
            write(*,'(A)') "1. Change spin a"
            write(*,'(A)') "2. Change charge Q"
            write(*,'(A)') "3. Back to main menu"
            write(*,'(A)') ""
            write(*,'(A)', advance='no') "Enter choice: "

            read(*,*,iostat=ios) choice
            if (ios /= 0) then
                call flush_input_line()
                cycle
            end if

            select case (choice)
            case (1)
                call set_real_in_range("Spin a", settings%spin_a, 0.0_wp, 0.999999_wp)
            case (2)
                call set_real_in_range("Charge Q", settings%charge_q, 0.0_wp, 0.999999_wp)
            case (3)
                done = .true.
            case default
                write(*,'(A)') ""
                write(*,'(A)') "Invalid option."
                call wait_for_enter()
            end select
        end do
    end subroutine change_black_hole_parameters


    subroutine change_temperature_scale(settings)
        type(settings_type), intent(inout) :: settings
        real(wp) :: t_new
        integer :: ios

        call clear_screen()
        write(*,'(A)') "========================================"
        write(*,'(A)') "       Change Temperature Scale"
        write(*,'(A)') "========================================"
        write(*,'(A)') ""
        write(*,'(A)', advance='no') "Enter temperature scale (K): "
        read(*,*,iostat=ios) t_new

        if (ios /= 0) then
            call flush_input_line()
            write(*,'(A)') "Invalid input."
            call wait_for_enter()
            return
        end if

        if (t_new <= 0.0_wp) then
            write(*,'(A)') ""
            write(*,'(A)') "Temperature scale must be positive."
            call wait_for_enter()
            return
        end if

        settings%temperature_scale = t_new
    end subroutine change_temperature_scale


    subroutine review_settings(settings)
        type(settings_type), intent(in) :: settings

        call clear_screen()
        write(*,'(A)') "========================================"
        write(*,'(A)') "           Review Settings"
        write(*,'(A)') "========================================"
        write(*,'(A)') ""
        call print_current_settings(settings)
    end subroutine review_settings


    logical function confirm_render(settings)
        type(settings_type), intent(in) :: settings
        character(len=1) :: answer

        call review_settings(settings)
        write(*,'(A)') ""
        write(*,'(A)', advance='no') "Proceed with render? (y/n): "
        read(*,'(A)') answer

        confirm_render = (answer == 'y' .or. answer == 'Y')
        calculate_cycels_frac = .true.
    end function confirm_render

    subroutine wait_for_enter()
        character(len=1) :: dummy
        write(*,'(A)') ""
        write(*,'(A)') "Press Enter to continue."
        read(*,'(A)') dummy
    end subroutine wait_for_enter


    subroutine flush_input_line()
        integer :: ios
        character(len=256) :: buffer
        read(*,'(A)',iostat=ios) buffer
    end subroutine flush_input_line


    subroutine set_positive_real(label, x)
        character(len=*), intent(in) :: label
        real(wp), intent(inout) :: x
        real(wp) :: temp
        integer :: ios

        call clear_screen()
        write(*,'(A)') trim(label)
        write(*,'(A)', advance='no') "Enter new value: "
        read(*,*,iostat=ios) temp

        if (ios /= 0) then
            call flush_input_line()
            write(*,'(A)') "Invalid input."
            call wait_for_enter()
            return
        end if

        if (temp <= 0.0_wp) then
            write(*,'(A)') "Value must be positive."
            call wait_for_enter()
            return
        end if

        x = temp
    end subroutine set_positive_real


    subroutine set_any_real(label, x)
        character(len=*), intent(in) :: label
        real(wp), intent(inout) :: x
        real(wp) :: temp
        integer :: ios

        call clear_screen()
        write(*,'(A)') trim(label)
        write(*,'(A)', advance='no') "Enter new value: "
        read(*,*,iostat=ios) temp

        if (ios /= 0) then
            call flush_input_line()
            write(*,'(A)') "Invalid input."
            call wait_for_enter()
            return
        end if

        x = temp
    end subroutine set_any_real


    subroutine set_real_in_range(label, x, xmin, xmax)
        character(len=*), intent(in) :: label
        real(wp), intent(inout) :: x
        real(wp), intent(in) :: xmin, xmax
        real(wp) :: temp
        integer :: ios

        call clear_screen()
        write(*,'(A)') trim(label)
        write(*,'(A, F8.3, A, F8.3, A)') "Allowed range: [", xmin, ", ", xmax, "]"
        write(*,'(A)', advance='no') "Enter new value: "
        read(*,*,iostat=ios) temp

        if (ios /= 0) then
            call flush_input_line()
            write(*,'(A)') "Invalid input."
            call wait_for_enter()
            return
        end if

        if (temp < xmin .or. temp > xmax) then
            write(*,'(A)') "Value out of range."
            call wait_for_enter()
            return
        end if

        x = temp
    end subroutine set_real_in_range

    subroutine progress_bar(count, nx, ny)
        implicit none
        integer, intent(in) :: count, nx, ny
        integer :: filled, percent
        character(len=30) :: bar

        integer, save :: number_of_pixels = 0
        integer, save :: last_filled = -1
        real(wp),    save :: pixels_per_step = 0.0

        if (calculate_cycels_frac) then
            number_of_pixels = nx * ny
            pixels_per_step = real(number_of_pixels) / 30.0
            calculate_cycels_frac = .false.
        end if

        filled = int(real(count) / pixels_per_step)
        percent = int(100.0 * real(count) / real(number_of_pixels))

        if (filled < 0)  filled = 0
        if (filled > 30) filled = 30

        if (filled /= last_filled) then
            bar = repeat('#', filled) // repeat('-', 30 - filled)
            write(*,'(A,"[",A,"] ",I3,"%")', advance='no') achar(13), bar, percent
            flush(6)
            last_filled = filled
        end if

    end subroutine progress_bar

end module utils
