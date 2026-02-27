!---------------------------------- io_mod -----------------------------------------!
! This module handles all input/output routines of the simulation.                  !
!                                                                                    !
! read_config  -> Reads simulation parameters and initializes system variables.     !
! read_files   -> Reads positions, velocities and box length from init_conf.dat.    !
! write_pdb    -> Writes PDB trajectory files (positions, box, particle types).     !
! write_log    -> Writes simulation log (energy, temperature vs step).              !
! write_final_config -> Writes final configuration to disk.                          !
! last_config  -> Checks existence of a final configuration file.                   !
!-----------------------------------------------------------------------------------!

module io_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use var_mod
    implicit none

    contains
    
    subroutine read_nml(sys, part, params)
        implicit none

        type(system_g),  intent(inout)      :: sys
        type(particles), intent(inout)      :: part
        type(parameters), intent(inout)     :: params

        integer :: unit, ios

        ! Namelist variables 
        integer :: num_particles, max_steps
        real(dp) :: rho, frac_particles, frac_charges
        real(dp) :: radius1, radius2, Z

        namelist /config_nml/ num_particles, max_steps, rho, frac_particles, &
                              frac_charges, radius1, radius2, Z

        !Default values
        num_particles   = 100
        rho             = 0.01_dp
        frac_particles  = 1.0_dp
        frac_charges    = 1.0_dp
        radius1         = 0.5_dp
        radius2         = 0.0_dp
        Z               = 100.0_dp
        max_steps       = 1000000

        open(newunit=unit, file=params%file_nml, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "Error: Unable to open file: "//trim(params%file_nml)

        read(unit, nml=config_nml, iostat=ios)
        if (ios /= 0) error stop "Error: Unable to read file: "//trim(params%file_nml)

        close(unit)

        ! System_g
        sys%num_particles  = num_particles
        sys%rho            = rho
        sys%frac_particles = frac_particles
        sys%frac_charges   = frac_charges

        ! Particles
        part%radius1 = radius1
        part%radius2 = radius2
        part%Z       = Z

        ! Parameters
        params%max_steps = max_steps

    end subroutine

    subroutine read_init_config(sys, part, params)
        implicit none

        type(system_g),   intent(inout) :: sys
        type(particles),  intent(inout) :: part
        type(parameters), intent(in)    :: params

        integer :: i, ios, unit

        open(newunit=unit, file=params%filename, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "Error: Unable to open file: "//trim(params%filename)

        read(unit, *, iostat=ios) sys%box
        if (ios /= 0) error stop "Error: Unable to read box size!"

        do i = 1, sys%num_particles
            read(unit, *, iostat=ios) &
                part%positions(1,i), part%positions(2,i), part%positions(3,i), &
                part%velocities(1,i), part%velocities(2,i), part%velocities(3,i)

            if (ios /= 0) then
                error stop "Error: Unable to read positions and velocities!"
            end if
        end do

        close(unit)
    end subroutine

    subroutine write_pdb(sys, part, outdir)
        implicit none

        type(system_g),  intent(in) :: sys
        type(particles), intent(in) :: part
        character(len=*), intent(in):: outdir

        integer            :: i, fdi, ios, idx, ir, ic
        character(len=512) :: filepath
        character(len=2), parameter :: labelP(4) = ['C', 'H', 'N', 'O']

        ir = 0
        ic = 0
 
        write(filepath, '(A,"/pdb/traj_rho_",F4.2,"_Z_",F6.1,"_N_",I6,".pdb")') &
            trim(outdir), sys%rho, part%Z, sys%num_particles

        open(newunit=fdi, file=filepath, status="unknown", &
             action="write", position="append", iostat=ios)

        if (ios /= 0) then
            error stop "Error: Unable to open file: "//trim(filepath)
        end if

        write(fdi,'(A6,3F9.3,3F7.2)') 'CRYST1', &
            sys%box, sys%box, sys%box, 90.0, 90.0, 90.0
        
        do i=1, sys%num_particles
            
            if (part%radius(i) .eq. part%radius1) ir = 0
            if (part%radius(i) .eq. part%radius2) ir = 2
            if (part%charges(i) .eq. part%Z) ic = 1
            if (part%charges(i) .eq. -part%Z) ic = 2
            
            idx = ir + ic 

            write(fdi,'(A6,I5,1X,A4,1X,A3,1X,I4,4X,3F8.3,2F6.2,10X,A2)') &
                'ATOM  ', i, labelP(idx), 'ION', 1, &
                part%positions(1,i), part%positions(2,i), part%positions(3,i), &
                1.0, part%radius(i)
           
        end do

        write(fdi,'(A6)') 'ENDMDL'
        close(fdi)

    end subroutine write_pdb

    subroutine write_log(sys, part, outdir, step, temp, is_first)
        implicit none

        type(system_g),     intent(in) :: sys
        type(particles),    intent(in) :: part
        character(len=*),   intent(in) :: outdir
        integer, intent(in)            :: step
        logical, intent(in)            :: is_first
        real(dp), intent(in)           :: temp  

        integer            :: fdi, ios
        character(len=512) :: filepath
        character(len=32)  :: datetime
        character(len=8)   :: date
        character(len=10)  :: time

        call date_and_time(date=date, time=time)
        datetime = date//' '//time

        write(filepath, '(A,"/log/log_rho_",F4.2,"_Z_",F4.1,"_N_",I6,".log")') &
            trim(outdir), sys%rho, part%Z, sys%num_particles

        open(newunit=fdi, file=filepath, status="unknown", &
             action="write", position="append", iostat=ios)

        if (ios /= 0) then
            error stop "Error: Unable to open log file: "//trim(filepath)
        end if

        if (is_first) then
            write(fdi, '(A)') "# Simulation Log"
            write(fdi, '(A)') "# Date/time: "//trim(datetime)
            write(fdi, '(A)') "# step        total_energy            temp"
        end if

        write(fdi, '(I10, 3X, F14.6, 3X, F12.6)') &
             step, sys%total_energy, temp

        close(fdi)

    end subroutine write_log

    subroutine write_final_config(sys, part, outdir)
        implicit none

        type(system_g), intent(in)   :: sys
        type(particles), intent(in)  :: part
        character(len=*), intent(in) :: outdir
        integer :: fdi, l, ios
        character(len=512) :: filepath

        write(filepath, '(A,"/config/final_config_N",I0,"_Z_",F5.1,"_rho_",F5.3,".dat")') &
             trim(outdir), sys%num_particles, part%Z, sys%rho

        open(newunit=fdi, file=filepath, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            error stop "Error: unable to open final configuration file: "//trim(filepath)
        end if

        write(fdi, '(F20.15)') sys%box

        do l = 1, sys%num_particles
            write(fdi, '(6F20.15)') part%positions(:,l), part%velocities(:,l)
        end do

        close(fdi)
    end subroutine

    subroutine last_config(sys, part, outdir, exist)
        implicit none

        type(system_g),   intent(in)  :: sys
        type(particles),  intent(in)  :: part
        logical,          intent(out) :: exist
        character(len=*), intent(in)  :: outdir
        character(len=512)            :: filepath

        write(filepath,'(A,"/config/final_config_N",I0,"_Z_",F5.1,"_rho_",F5.3,".dat")') &
         trim(outdir), sys%num_particles, part%Z, sys%rho

        inquire(file=filepath, exist=exist)
    end subroutine
end module