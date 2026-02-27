program init_conf
    use, intrinsic:: iso_fortran_env, only: dp => real64
    use io_mod
    use var_mod
    implicit none

    type(system_g)   :: sys
    type(particles)  :: part
    type(parameters) :: params
    
    call init_params(params)
    call get_command_argument(1, params%outdir)
    call read_nml(sys, part, params)

    call init_system(sys)
    call init_basic(part, sys)
    call last_config(sys, part, params%outdir, params%exists)


    if(params%exists) then
        call continue_config(sys, part, params%outdir)
    else
        call new_config(part, sys, params)
    end if

    call end_basic(part)
    call end_system(sys)

    contains

    subroutine lattice_positions(pos, box, n)
        implicit none

        real(dp), dimension(:,:), intent(out) :: pos
        real(dp), intent(in)                  :: box
        integer, intent(in)                   :: n
        integer                               :: nx, ny, nz, i, j, k, p
        integer                               :: n_side
        real(dp)                              :: space, delta_x, delta_y, delta_z

        n_side = ceiling(n**(1.0_dp/3.0_dp))
        space = box/real(n_side, dp)

        p = 0

        do i = 0, n_side-1
            do j = 0, n_side-1
                do k = 0, n_side-1
                    p = p+1
                    
                    if(p > n) exit
                    
                    ! Gera números aleatórios diferentes para cada coordenada
                    call random_number(delta_x)
                    call random_number(delta_y)
                    call random_number(delta_z)
                    
                    pos(1, p) = (i+0.5_dp)*space + (delta_x-0.5_dp)*0.01_dp*space
                    pos(2, p) = (j+0.5_dp)*space + (delta_y-0.5_dp)*0.01_dp*space
                    pos(3, p) = (k+0.5_dp)*space + (delta_z-0.5_dp)*0.01_dp*space
                end do
                if(p > n) exit
            end do
            if(p > n) exit
        end do

    end subroutine lattice_positions
    
    subroutine rand_velocities(vel, sig,n)
        implicit none

        real(dp), intent(in)                  :: sig
        integer, intent(in)                   :: n
        real(dp), dimension(:,:), intent(out) :: vel
        real(dp), dimension(3)                :: vcm
        real(dp)                              :: u1, u2, u3, u4, R1, R2, theta1, theta2
        integer                               :: i
        
        do i = 1, n
            ! Gera 4 números aleatórios independentes
            call random_number(u1)
            call random_number(u2)
            call random_number(u3)
            call random_number(u4)
                
            ! Evita log(0)
            if (u1 < 1.0e-12_dp) u1 = 1.0e-12_dp
            if (u3 < 1.0e-12_dp) u3 = 1.0e-12_dp
                
            ! Transformação de Box-Muller
            R1 = sqrt(-2.0_dp*log(u1))
            theta1 = 2.0_dp*pi*u2
            R2 = sqrt(-2.0_dp*log(u3))
            theta2 = 2.0_dp*pi*u4
                
            vel(1, i) = sig * R1 * cos(theta1)
            vel(2, i) = sig * R1 * sin(theta1)  
            vel(3, i) = sig * R2 * cos(theta2)
        end do
        
        vcm = 0.0_dp
        do i = 1, 3
            vcm(i) = sum(vel(i,:)) / real(n, dp)
        end do

        do i = 1, n
            vel(:, i) = vel(:, i) - vcm
        end do
        
    end subroutine rand_velocities

    subroutine new_config(part, sys, params)
        implicit none
        type(particles),  intent(inout) :: part
        type(system_g),   intent(inout) :: sys
        type(parameters), intent(in)    :: params
        
        real(dp), dimension(:), allocatable :: r     
        
        real(dp) :: sigma, Vt
        integer  :: fdi, ioerr, i, n

        Vt = 0.0_dp
    
        allocate(r(sys%num_particles))

        if (sys%frac_particles >= 0.0_dp .and. sys%frac_particles <= 1.0_dp) then
            n = max(1, nint(sys%frac_particles * sys%num_particles))
        else
            error stop "Error: frac_particles must be between 0 and 1"
        end if
        
        if(n .eq. sys%num_particles) then
            r(1:n) = part%radius1
        else
            r(1:n) = part%radius1
            r(n+1:sys%num_particles) = part%radius2 
        end if

        do i = 1, sys%num_particles
            Vt = Vt + (4.0_dp/3.0_dp)* pi * r(i)**3
        end do
        
        sys%box = (Vt/sys%rho)**(1.0_dp/3.0_dp)
        sigma = sqrt(sys%temp_target)

        call lattice_positions(part%positions, sys%box, sys%num_particles)
        call rand_velocities(part%velocities, sigma, sys%num_particles)
        
        open(newunit = fdi, file = params%filename, status = 'replace', iostat = ioerr)
        
        if(ioerr /= 0) then
            write(*,*) 'Error opening: ', params%filename, 'IOSTAT= ', ioerr
            stop 'Error in write initial configuration!'
        end if
        
        write(fdi, '(f20.15)') sys%box
        
        do i = 1, sys%num_particles
            write(fdi, '(6F20.15)') part%positions(:, i), part%velocities(:, i)
        end do

        close(fdi)

        deallocate(r)
    end subroutine new_config

    subroutine continue_config(sys, part, outdir)
    implicit none
        type(system_g), intent(inout)  :: sys
        type(particles), intent(inout) :: part

        integer :: fdi, ioerr, i
        character(len=256) :: filepath
        character(len=*) :: outdir
    
        write(filepath,'(A,"/config/final_config_N",I0,"_Z_",F5.1,"_rho_",F5.3,".dat")') &
         trim(outdir), sys%num_particles, part%Z, sys%rho


        open(newunit = fdi, file = filepath, status = 'old', action = 'read', iostat = ioerr)
        if (ioerr /= 0) then
            write(*,*) 'Error opening: ', trim(filepath)
            stop
        endif

        read(fdi, '(f20.15)') sys%box

        do i = 1, sys%num_particles
            read(fdi, '(6F20.15)') part%positions(:, i), part%velocities(:, i)
        end do

        close(fdi)
    end subroutine
end program init_conf
