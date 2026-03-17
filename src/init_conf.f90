program init_conf
    use, intrinsic:: iso_fortran_env, only: dp => real64
    use io_mod
    use var_mod
    implicit none

    type(system_g)   :: sys
    type(particles)  :: part
    type(parameters) :: params

    call random_seed()
    
    call init_params(params)
    call get_command_argument(1, params%outdir)
    call read_nml(sys, part, params)

    call init_system(sys)
    call init_particles(part, sys)
    call last_config(sys, part, params%outdir, params%exists)


    if(params%exists) then
        call continue_config(sys, part, params%outdir)
    else
        call new_config(part, sys, params)
    end if

    call end_basic(part)
    call end_system(sys)

    contains

    subroutine lattice_positions(pos, box, n, radii)
        implicit none

        real(dp), dimension(:,:), intent(out) :: pos
        real(dp), intent(in)                  :: box
        integer, intent(in)                   :: n
        real(dp), dimension(:), intent(in)    :: radii
        integer                               :: i, j, k, p, l, attempt
        integer                               :: n_side
        real(dp)                              :: space, delta(3), shift
        real(dp), dimension(3)                :: rij, pos_tent
        real(dp)                              :: rij2, min_dist2
        logical                               :: overlap
        integer, parameter                    :: max_attempts = 1000
        real(dp), parameter                    :: perturb_frac = 0.2_dp

        n_side = ceiling(n**(1.0_dp/3.0_dp))
        space = box / real(n_side, dp)
        shift = -box / 2.0_dp
        
        p = 0

        do i = 0, n_side-1
            do j = 0, n_side-1
                do k = 0, n_side-1
                    if (p >= n) exit
                    
                    p = p + 1
                    attempt = 0
                    overlap = .true.
                    
                    do while (overlap .and. attempt < max_attempts)
                        attempt = attempt + 1
                        
                        call random_number(delta)
                        delta = (delta - 0.5_dp) * 2.0_dp * perturb_frac * space
                        
                        pos_tent(1) = shift + (i + 0.5_dp) * space + delta(1)
                        pos_tent(2) = shift + (j + 0.5_dp) * space + delta(2)
                        pos_tent(3) = shift + (k + 0.5_dp) * space + delta(3)
                        
                        overlap = .false.
                        do l = 1, p-1
                            rij = pos_tent - pos(:, l)
                            rij = rij - box * nint(rij / box)
                            rij2 = dot_product(rij, rij)
                            
                            min_dist2 = (radii(p) + radii(l))**2
                            
                            if (rij2 < min_dist2) then
                                overlap = .true.
                                exit
                            end if
                        end do

                    end do
                    if (overlap) error stop "Could not place particle without overlap"

                    pos(:, p) = pos_tent
                end do
            end do
        end do

    end subroutine lattice_positions
    
    subroutine rand_velocities(vel, sig, n)
        real(dp), intent(in) :: sig
        integer, intent(in) :: n
        real(dp), dimension(:,:), intent(out) :: vel
        real(dp) :: u1, u2, u3, u4, r2
        real(dp), dimension(3) :: vcm 
        integer :: i
        
        do i = 1, n
            do
                call random_number(u1)  
                call random_number(u2)  
                call random_number(u3)
                call random_number(u4)

                if (u3 < 1.0e-12_dp) u3 = 1.0e-12_dp

                u1 = 2.0_dp*u1 - 1.0_dp
                u2 = 2.0_dp*u2 - 1.0_dp
                r2 = u1*u1 + u2*u2
                if (r2 < 1.0_dp .and. r2 > 1.0e-12_dp) exit
            end do
            
            vel(1, i) = sig * sqrt(-2.0_dp * log(r2) / r2) * u1
            vel(2, i) = sig * sqrt(-2.0_dp * log(r2) / r2) * u2
            vel(3, i) = sig * sqrt(-2.0_dp * log(u3)) * cos(2.0_dp * pi * u4)
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
        
        real(dp) :: sigma, Vt, a0
        integer  :: fdi, ioerr, i, n

        Vt = (4.0_dp*pi/3.0_dp) * sum((part%radius)**3)
        
        sys%box = (Vt/sys%rho)**(1.0_dp/3.0_dp)
        
        sigma = sqrt(sys%temp_target)

        call lattice_positions(part%positions, sys%box, sys%num_particles, part%radius)
        call rand_velocities(part%velocities, sigma, sys%num_particles)
        
        open(newunit = fdi, file = params%filename, status = 'replace', iostat = ioerr)
        
        if(ioerr /= 0) then
            write(*,*) 'Error opening: ', params%filename, 'IOSTAT= ', ioerr
            stop 'Error in write initial configuration!'
        end if
        
        write(fdi,'(ES25.16)') sys%box
        
        do i = 1, sys%num_particles
            write(fdi,'(6ES25.16)') part%positions(:,i), part%velocities(:,i)
        end do

        close(fdi)
    end subroutine new_config

    subroutine continue_config(sys, part, outdir)
    implicit none
        type(system_g), intent(inout)  :: sys
        type(particles), intent(inout) :: part

        integer :: fdi, ioerr, i
        character(len=256) :: filepath
        character(len=*) :: outdir
    
        write(filepath,'(A,"/config/final_config_N",I0,"_Z_",F5.1,"_rho_",F0.3,".dat")') &
         trim(outdir), sys%num_particles, part%Z, sys%rho


        open(newunit = fdi, file = filepath, status = 'old', action = 'read', iostat = ioerr)
        if (ioerr /= 0) then
            write(*,*) 'Error opening: ', trim(filepath)
            stop
        endif

        read(fdi, '(ES25.16)') sys%box

        do i = 1, sys%num_particles
            read(fdi, '(6ES25.16)') part%positions(:, i), part%velocities(:, i)
        end do

        close(fdi)
    end subroutine

end program init_conf