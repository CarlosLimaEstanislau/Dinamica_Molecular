module ui_mod
	use, intrinsic :: iso_fortran_env, only: dp => real64
    use var_mod
    implicit none

    integer, parameter :: UI_DATA_LINE = 21

    contains

    subroutine initialize(sys, params, part, dt_event)
        implicit none
        real(dp), intent(in)        :: dt_event
        type(system_g), intent(in)  :: sys
        type(parameters), intent(in):: params
        type(particles), intent(in) :: part

        print *
        print *, '====================================================='
        print *, '        Molecular Dynamics Simulation Started        '
        print *, '====================================================='
        print *

        write(*,'(A,I10)')       ' Number of particles:', sys%num_particles
        write(*,'(A,F10.4)')     ' Density (rho):      ', sys%rho
        write(*,'(A,F10.4)')     ' Box length:         ', sys%box
        write(*,'(A,I10)')       ' Total steps:        ', params%max_steps
        write(*,'(A,F10.6)')     ' Time step (dt):     ', dt_event
        write(*,'(A,F10.4)')     ' Target temperature: ', sys%temp_target
        write(*,'(A,F10.4)')     ' Frac. radii type A: ', sys%frac_particles
        write(*,'(A,F10.4)')     ' Frac. charged:      ', sys%frac_charges
        write(*,'(A,F10.4)')     ' Charge magnitude Z: ', part%Z
        write(*,'(A,F10.4)')     ' Langevin gamma:     ', params%gamma

        print *
        print *, '-----------------------------------------------------'
        print *, ' Step        Energy           Temperature           '
        print *, '-----------------------------------------------------'
    end subroutine
    
    subroutine write_out(step, energy, temp)
        implicit none
        real(dp), intent(in) :: energy, temp
        integer, intent(in)  :: step
        character(len=1), parameter :: esc = achar(27)
        character(len=8) :: line_str

        write(line_str,'(I0)') UI_DATA_LINE
        write(*,'(A)', advance='no') esc//'['//trim(line_str)//';1H'
        write(*,'(A)', advance='no') esc//'[0G'

        write(*,'(I7, A, F14.4, A, F12.4, A)') &
             step, '  ', energy, '  ', temp, '         '
        print *, '-----------------------------------------------------'
        flush(6)
    end subroutine

end module