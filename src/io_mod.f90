module io_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use var_mod, only: system_g, particles, parameters, a0
    implicit none
    private

    public :: read_nml
    public :: read_config, write_config
    public :: read_part_file, write_part_file
    public :: write_pdb, write_log, write_final_config

contains

    subroutine read_nml(sys, part, params)
        implicit none
        type(system_g),  intent(inout) :: sys
        type(particles), intent(inout) :: part
        type(parameters), intent(inout) :: params

        integer :: unit, ios

        integer  :: num_particles, max_steps
        real(dp) :: rho, frac_particles, frac_charges
        real(dp) :: radius1, radius2, Z, gamma

        namelist /config_nml/ num_particles, max_steps, rho, frac_particles, &
                              frac_charges, radius1, radius2, Z, gamma

        ! Defaults
        num_particles  = 100
        max_steps      = 1000000
        rho            = 0.01_dp
        frac_particles = 1.0_dp
        frac_charges   = 1.0_dp
        radius1        = 0.5_dp
        radius2        = 0.0_dp
        Z              = 100.0_dp
        gamma = 1.0_dp

        open(newunit=unit, file=params%file_nml, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Erro ao abrir namelist: ', trim(params%file_nml), ' IOSTAT=', ios
            error stop 'read_nml: falha ao abrir arquivo'
        end if

        read(unit, nml=config_nml, iostat=ios)

        if (num_particles < 2) error stop "read_nml: num_particles deve ser >= 2"
        if (radius1 <= 0.0_dp) error stop "read_nml: radius1 deve ser > 0"
        if (radius2 <= 0.0_dp) error stop "read_nml: radius2 deve ser > 0"

        if (ios /= 0) then
            close(unit)
            write(*,*) 'Erro ao ler namelist: ', trim(params%file_nml), ' IOSTAT=', ios
            error stop 'read_nml: falha ao ler arquivo'
        end if

        close(unit, iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Aviso: erro ao fechar namelist: ', trim(params%file_nml), ' IOSTAT=', ios
        end if

        ! Validação
        if (num_particles <= 0) error stop 'read_nml: num_particles deve ser > 0'
        if (max_steps <= 0)     error stop 'read_nml: max_steps deve ser > 0'
        if (rho <= 0.0_dp)      error stop 'read_nml: rho deve ser > 0'
        if (frac_particles < 0.0_dp .or. frac_particles > 1.0_dp) then
            error stop 'read_nml: frac_particles deve estar entre 0 e 1'
        end if
        if (frac_charges < 0.0_dp .or. frac_charges > 1.0_dp) then
            error stop 'read_nml: frac_charges deve estar entre 0 e 1'
        end if
        if (radius1 < 0.0_dp .or. radius2 < 0.0_dp) then
            error stop 'read_nml: raios devem ser >= 0'
        end if

        sys%num_particles  = num_particles
        sys%rho            = rho
        sys%frac_particles = frac_particles
        sys%frac_charges   = frac_charges
        sys%gamma          = gamma

        part%radius1 = radius1
        part%radius2 = radius2
        part%Z       = Z


        params%max_steps = max_steps
    end subroutine read_nml

    subroutine read_config(sys, part, filename)
        implicit none
        type(system_g),  intent(inout) :: sys
        type(particles), intent(inout) :: part
        character(len=*), intent(in)   :: filename

        integer :: i, ios, unit

        if (.not. allocated(part%positions))  error stop 'read_config: positions nao alocado'
        if (.not. allocated(part%velocities)) error stop 'read_config: velocities nao alocado'

        if (size(part%positions,1) /= 3)  error stop 'read_config: positions deve ser (3,N)'
        if (size(part%velocities,1) /= 3) error stop 'read_config: velocities deve ser (3,N)'
        if (size(part%positions,2) /= sys%num_particles) then
            error stop 'read_config: positions inconsistente com num_particles'
        end if
        if (size(part%velocities,2) /= sys%num_particles) then
            error stop 'read_config: velocities inconsistente com num_particles'
        end if

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Erro ao abrir config: ', trim(filename), ' IOSTAT=', ios
            error stop 'read_config: falha ao abrir arquivo'
        end if

        read(unit, *, iostat=ios) sys%box
        if (ios /= 0) then
            close(unit)
            error stop 'read_config: falha ao ler box'
        end if
        if (sys%box <= 0.0_dp) then
            close(unit)
            error stop 'read_config: box deve ser positivo'
        end if

        do i = 1, sys%num_particles
            read(unit, *, iostat=ios) &
                part%positions(1,i), part%positions(2,i), part%positions(3,i), &
                part%velocities(1,i), part%velocities(2,i), part%velocities(3,i)

            if (ios /= 0) then
                close(unit)
                write(*,*) 'Erro ao ler configuracao na particula ', i
                error stop 'read_config: falha ao ler posicoes/velocidades'
            end if
        end do

        close(unit, iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Aviso: erro ao fechar config: ', trim(filename), ' IOSTAT=', ios
        end if
    end subroutine read_config

    subroutine write_config(part, sys, filename)
        implicit none
        type(particles), intent(in)    :: part
        type(system_g),  intent(in)    :: sys
        character(len=*), intent(in)   :: filename

        integer :: fdi, ioerr, i

        if (.not. allocated(part%positions))  error stop 'write_config: positions nao alocado'
        if (.not. allocated(part%velocities)) error stop 'write_config: velocities nao alocado'

        if (size(part%positions,1) /= 3) error stop 'write_config: positions deve ser (3,N)'
        if (size(part%velocities,1) /= 3) error stop 'write_config: velocities deve ser (3,N)'
        if (size(part%positions,2) /= sys%num_particles) error stop 'write_config: positions inconsistente'
        if (size(part%velocities,2) /= sys%num_particles) error stop 'write_config: velocities inconsistente'
        if (sys%box <= 0.0_dp) error stop 'write_config: box deve ser positivo'

        open(newunit=fdi, file=filename, status='replace', action='write', iostat=ioerr)
        if (ioerr /= 0) then
            write(*,*) 'Erro ao abrir: ', trim(filename), ' IOSTAT=', ioerr
            error stop 'write_config: falha ao abrir arquivo'
        end if

        write(fdi,'(ES25.16)') sys%box

        do i = 1, sys%num_particles
            write(fdi,'(6ES25.16)') part%positions(:,i), part%velocities(:,i)
        end do

        close(fdi, iostat=ioerr)
        if (ioerr /= 0) then
            write(*,*) 'Aviso: erro ao fechar: ', trim(filename), ' IOSTAT=', ioerr
        end if
    end subroutine write_config

    subroutine write_part_file(part, sys, filename)
        implicit none
        type(particles), intent(in)    :: part
        type(system_g),  intent(in)    :: sys
        character(len=*), intent(in)   :: filename

        integer :: fdi, ioerr, i

        if (.not. allocated(part%radius))  error stop 'write_part_file: radius nao alocado'
        if (.not. allocated(part%charges)) error stop 'write_part_file: charges nao alocado'

        if (size(part%radius) /= sys%num_particles) then
            error stop 'write_part_file: radius inconsistente com num_particles'
        end if
        if (size(part%charges) /= sys%num_particles) then
            error stop 'write_part_file: charges inconsistente com num_particles'
        end if

        open(newunit=fdi, file=filename, status='replace', action='write', iostat=ioerr)
        if (ioerr /= 0) then
            write(*,*) 'Erro ao abrir: ', trim(filename), ' IOSTAT=', ioerr
            error stop 'write_part_file: falha ao abrir arquivo'
        end if

        write(fdi,'(I0)') sys%num_particles

        do i = 1, sys%num_particles
            write(fdi,'(2ES25.16)') part%radius(i), part%charges(i)
        end do

        close(fdi, iostat=ioerr)
        if (ioerr /= 0) then
            write(*,*) 'Aviso: erro ao fechar: ', trim(filename), ' IOSTAT=', ioerr
        end if
    end subroutine write_part_file

    subroutine read_part_file(part, sys, filename)
        implicit none
        type(particles), intent(inout) :: part
        type(system_g),  intent(in)    :: sys
        character(len=*), intent(in)   :: filename

        integer :: fdi, ioerr, i, nfile

        if (.not. allocated(part%radius))  error stop 'read_part_file: radius nao alocado'
        if (.not. allocated(part%charges)) error stop 'read_part_file: charges nao alocado'

        if (size(part%radius) /= sys%num_particles) then
            error stop 'read_part_file: radius inconsistente com num_particles'
        end if
        if (size(part%charges) /= sys%num_particles) then
            error stop 'read_part_file: charges inconsistente com num_particles'
        end if

        open(newunit=fdi, file=filename, status='old', action='read', iostat=ioerr)
        if (ioerr /= 0) then
            write(*,*) 'Erro ao abrir: ', trim(filename), ' IOSTAT=', ioerr
            error stop 'read_part_file: falha ao abrir arquivo'
        end if

        read(fdi,*,iostat=ioerr) nfile
        if (ioerr /= 0) then
            close(fdi)
            error stop 'read_part_file: erro ao ler numero de particulas'
        end if

        if (nfile /= sys%num_particles) then
            close(fdi)
            error stop 'read_part_file: numero de particulas inconsistente'
        end if

        do i = 1, sys%num_particles
            read(fdi,*,iostat=ioerr) part%radius(i), part%charges(i)
            if (ioerr /= 0) then
                close(fdi)
                write(*,*) 'Erro ao ler .part na particula ', i
                error stop 'read_part_file: erro ao ler radius/charges'
            end if
        end do

        close(fdi, iostat=ioerr)
        if (ioerr /= 0) then
            write(*,*) 'Aviso: erro ao fechar: ', trim(filename), ' IOSTAT=', ioerr
        end if
    end subroutine read_part_file

    subroutine write_pdb(sys, part, outdir)
        implicit none
        type(system_g),  intent(in) :: sys
        type(particles), intent(in) :: part
        character(len=*), intent(in):: outdir

        integer :: i, fdi, ios, idx, ir, ic
        real(dp), parameter :: tol = 1.0e-10_dp
        character(len=512) :: filepath
        character(len=2), parameter :: labelP(2) = ['C ', 'H ']

        if (.not. allocated(part%positions)) error stop 'write_pdb: positions nao alocado'
        if (.not. allocated(part%radius))    error stop 'write_pdb: radius nao alocado'
        if (.not. allocated(part%charges))   error stop 'write_pdb: charges nao alocado'

        write(filepath, '(A,"/pdb/traj_rho_",F4.2,"_Z_",F6.1,"_N_",I6,".pdb")') &
            trim(outdir), sys%rho, part%Z, sys%num_particles

        open(newunit=fdi, file=filepath, status='unknown', action='write', &
             position='append', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Erro ao abrir pdb: ', trim(filepath), ' IOSTAT=', ios
            error stop 'write_pdb: falha ao abrir arquivo'
        end if

        write(fdi,'(A6,3F9.3,3F7.2)') 'CRYST1', sys%box, sys%box, sys%box, 90.0, 90.0, 90.0

        do i = 1, sys%num_particles
            if (part%charges(i) == -part%Z) then
                idx = 1
            else if (part%charges(i) == part%Z) then
                idx = 2
            else
                close(fdi)
                write(*,*) 'Particula com raio nao identificado: i=', i, ' raio=', part%radius(i)
                error stop 'write_pdb: classificacao de particula falhou'
            end if

            write(fdi,'(A6,I5,1X,A4,1X,A3,1X,I4,4X,3F8.3,2F6.2,10X,A2)') &
                'ATOM  ', i, labelP(idx), 'ION', 1, &
                part%positions(1,i), part%positions(2,i), part%positions(3,i), &
                1.0_dp, part%radius(i)
        end do

        write(fdi,'(A6)') 'ENDMDL'

        close(fdi, iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Aviso: erro ao fechar pdb: ', trim(filepath), ' IOSTAT=', ios
        end if
    end subroutine write_pdb

    subroutine write_log(sys, part, outdir, step, temp, is_first)
        implicit none
        type(system_g),   intent(in) :: sys
        type(particles),  intent(in) :: part
        character(len=*), intent(in) :: outdir
        integer,          intent(in) :: step
        logical,          intent(in) :: is_first
        real(dp),         intent(in) :: temp

        integer :: fdi, ios
        character(len=512) :: filepath
        character(len=8)   :: date
        character(len=10)  :: time
        character(len=32)  :: datetime

        call date_and_time(date=date, time=time)
        datetime = date//' '//time

        write(filepath, '(A,"/log/log_rho_",F4.2,"_Z_",F4.1,"_N_",I6,".txt")') &
            trim(outdir), sys%rho, part%Z, sys%num_particles

        open(newunit=fdi, file=filepath, status='unknown', action='write', &
             position='append', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Erro ao abrir log: ', trim(filepath), ' IOSTAT=', ios
            error stop 'write_log: falha ao abrir arquivo'
        end if

        if (is_first) then
            write(fdi,'(A)') '# Simulation Log'
            write(fdi,'(A)') '# Date/time: '//trim(datetime)
            write(fdi,'(A)') '# step        total_energy            temp'
        end if

        write(fdi,'(I10,3X,F14.6,3X,F12.6)') step, sys%total_energy, temp

        close(fdi, iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Aviso: erro ao fechar log: ', trim(filepath), ' IOSTAT=', ios
        end if
    end subroutine write_log

    subroutine write_final_config(sys, part, outdir)
        implicit none
        type(system_g),  intent(in) :: sys
        type(particles), intent(in) :: part
        character(len=*), intent(in) :: outdir

        integer :: fdi, ios, i
        character(len=512) :: filepath

        if (.not. allocated(part%positions))  error stop 'write_final_config: positions nao alocado'
        if (.not. allocated(part%velocities)) error stop 'write_final_config: velocities nao alocado'
        if (sys%box <= 0.0_dp) error stop 'write_final_config: box deve ser positivo'

        write(filepath, '(A,"/config/final_config_N",I0,"_Z_",F5.1,"_rho_",F5.3,".dat")') &
             trim(outdir), sys%num_particles, part%Z, sys%rho

        open(newunit=fdi, file=filepath, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Erro ao abrir final config: ', trim(filepath), ' IOSTAT=', ios
            error stop 'write_final_config: falha ao abrir arquivo'
        end if

        write(fdi,'(ES25.16)') sys%box
        do i = 1, sys%num_particles
            write(fdi,'(6ES25.16)') part%positions(:,i), part%velocities(:,i)
        end do

        close(fdi, iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Aviso: erro ao fechar final config: ', trim(filepath), ' IOSTAT=', ios
        end if
    end subroutine write_final_config

end module io_mod