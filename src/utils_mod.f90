module utils_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use rng_mod, only: random_normal
    implicit none
    private

    public :: check_interruption, compute_observables, validate_configuration
    public :: u1_propagator, u2_propagator, langevin_step

contains

    subroutine check_interruption(outdir, interrupted)
        implicit none
        character(len=*), intent(in) :: outdir
        logical, intent(out) :: interrupted
        logical :: exists
        character(len=512) :: stopfile

        write(stopfile,'(A,"/STOP.txt")') trim(outdir)
        inquire(file=stopfile, exist=exists)

        interrupted = exists
    end subroutine check_interruption

    subroutine u1_propagator(t, r, v, box)
        implicit none
        real(dp), intent(in) :: t, box
        real(dp), intent(inout) :: r(:,:), v(:,:)

        r = r + t * v
        r = r - box * nint(r / box)
    end subroutine u1_propagator

    subroutine u2_propagator(t, velocities, forces, masses)
        implicit none
        real(dp), intent(in) :: t
        real(dp), intent(inout) :: velocities(:,:)
        real(dp), intent(in) :: forces(:,:), masses(:)
        integer :: i, n

        n = size(velocities, 2)
        do i = 1, n
            velocities(:,i) = velocities(:,i) + (t / masses(i)) * forces(:,i)
        end do
    end subroutine u2_propagator

    subroutine langevin_step(dt, velocities, masses, gamma, temp)
        implicit none
        real(dp), intent(in) :: dt, gamma, temp
        real(dp), intent(inout) :: velocities(:,:)
        real(dp), intent(in) :: masses(:)

        integer :: i, d, n
        real(dp) :: c1, c2, g

        if (any(masses <= 0.0_dp)) then
            error stop "compute_observables: todas as massas devem ser positivas"
        end if

        n  = size(velocities, 2)
        c1 = exp(-gamma * dt)

        do i = 1, n
            c2 = sqrt((1.0_dp - c1*c1) * temp / masses(i))
            do d = 1, size(velocities,1)
                g = random_normal(0.0_dp, 1.0_dp)
                velocities(d,i) = c1 * velocities(d,i) + c2 * g
            end do
        end do
    end subroutine langevin_step

    subroutine compute_observables(masses, velocities, g, total_potential, total_energy, temp)
        implicit none
        integer, intent(in) :: g
        real(dp), dimension(:,:), intent(in) :: velocities
        real(dp), dimension(:), intent(in)   :: masses
        real(dp), intent(in) :: total_potential
        real(dp), intent(out) :: total_energy, temp
        real(dp) :: kinetic
        
        if (size(velocities,1) /= 3) then
            error stop "compute_observables: velocities deve ter dimensão (3,N)"
        end if

        if (size(velocities,2) /= size(masses)) then
            error stop "compute_observables: tamanho de masses inconsistente com velocities"
        end if
        
        kinetic = 0.5_dp * sum(masses * sum(velocities**2, dim=1))
        total_energy = kinetic + total_potential
        
        if (g <= 0) error stop "compute_observables: g deve ser > 0"
        
        temp = (2.0_dp * kinetic) / real(g, dp)

    end subroutine compute_observables

    subroutine validate_configuration(positions, radius, box)
        implicit none
        real(dp), intent(in) :: positions(:,:), radius(:), box
        integer :: i, j, n
        real(dp) :: rij(3), r2, r, contact, tol

        if (box <= 0.0_dp) then
            error stop "validate_configuration: box deve ser positivo"
        end if

        if (size(positions,1) /= 3) then
            error stop "validate_configuration: positions deve ter dimensão (3,N)"
        end if

        if (size(positions,2) /= size(radius)) then
            error stop "validate_configuration: positions e radius incompatíveis"
        end if

        if (any(radius < 0.0_dp)) then
            error stop "validate_configuration: existem raios negativos"
        end if

        n = size(radius)
        tol = 1.0e-10_dp

        do i = 1, n-1
            do j = i+1, n
                rij = positions(:,i) - positions(:,j)
                rij = rij - box * nint(rij / box)

                r2 = sum(rij*rij)

                if (r2 <= tiny(1.0_dp)) then
                    write(*,*) "Particulas praticamente coincidentes: ", i, j
                    error stop "validate_configuration: distancia ~ 0"
                end if

                r = sqrt(r2)
                contact = radius(i) + radius(j)

                if (r < contact - tol) then
                    write(*,*) "Overlap detectado entre particulas: ", i, j
                    write(*,*) "distancia = ", r
                    write(*,*) "contato   = ", contact
                    error stop "validate_configuration: configuracao invalida"
                end if
            end do
        end do
    end subroutine validate_configuration
end module utils_mod