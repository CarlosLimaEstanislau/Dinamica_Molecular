module handler_mod
    use, intrinsic:: iso_fortran_env, only: dp => real64
    use var_mod, only: pi, dt
    implicit none

    contains
    !-----------------Utilidades --------------------------------------!
        function gaussian_random(mu, sigma) result(x)
            implicit none
            real(dp), intent(in):: mu, sigma
            real(dp):: x
            real(dp):: u1, u2

            call random_number(u1)
            call random_number(u2)
            if (u1 < 1.0e-12_dp) u1 = 1.0e-12_dp
            x = mu + sigma * sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*pi*u2)
        end function gaussian_random

        subroutine check_interruption(interrupted)
            implicit none
            logical, intent(out) :: interrupted
            logical :: exists

            inquire(file="\simulation\STOP.txt", exist=exists)

            interrupted = exists
        end subroutine

        subroutine random_gaussian(x)
            implicit none
            real(dp), intent(out) :: x
            real(dp) :: u1, u2

            call random_number(u1)
            call random_number(u2)

            x = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
        end subroutine random_gaussian

        function random_normal ( mean, std ) RESULT ( r )
            implicit none
            real(dp)             :: r    
            real(dp), intent(in) :: mean 
            real(dp), intent(in) :: std  

            real(dp), dimension(2) :: zeta
            REAL(dp), save :: s
            logical,  save :: saved = .false.

            if( saved )then     
                r = s              
                r = mean + std * r 
                saved = .false.    
            else                                              
                call random_number (zeta)                      
                r = sqrt(-2.0*log(zeta(1)))*cos(2*pi*zeta(2)) 
                s = sqrt(-2.0*log(zeta(1)))*sin(2*pi*zeta(2)) 
                r = mean + std * r                             
                saved = .true.                                
            end if
        end function random_normal

        subroutine random_normals ( mean, std, r )
            implicit none 
            real(dp),               intent(in)  :: mean 
            real(dp),               intent(in)  :: std  
            real(dp), dimension(:), intent(out) :: r    

            integer :: i

            do i = 1, size(r)
               r(i) = random_normal ( mean, std )
            end do
        end subroutine random_normals

        subroutine collision_hard(positions, velocities, radius, masses, box)
            implicit none
            real(dp), dimension(:,:), intent(inout) :: positions, velocities
            real(dp), dimension(:), intent(in) :: masses, radius
            real(dp), intent(in) :: box
            real(dp), parameter :: tiny = 1.0e-12_dp
            real(dp), dimension(3) :: rij_hat, rij, vij
            real(dp) :: coeff_i, coeff_j, overlap
            real(dp) ::rij_norm, vij_rel, delta

            integer :: n, i, j

            n = size(positions,2)

            do i = 1, n-1
                do j = i+1, n
                    delta = radius(i) + radius(j)
                    
                    rij = positions(:,i) - positions(:,j)
                    rij = rij - box * nint(rij / box)
                    rij_norm = sqrt(dot_product(rij, rij))
                    rij_hat = rij / rij_norm

                    vij = velocities(:,i) - velocities(:,j)
                    vij_rel = dot_product(vij, rij_hat)
                    
                    if (rij_norm < tiny) cycle     

                    if (rij_norm < delta .and. vij_rel < 0.0_dp) then

                        coeff_i = 2.0_dp * masses(j) / (masses(i) + masses(j)) * vij_rel
                        coeff_j = 2.0_dp * masses(i) / (masses(i) + masses(j)) * vij_rel    
                        
                        overlap = delta - rij_norm

                        positions(:,i) = positions(:,i) - 0.5_dp*overlap * rij_hat
                        positions(:,j) = positions(:,j) + 0.5_dp*overlap * rij_hat

                        velocities(:,i) = velocities(:,i) - coeff_i * rij_hat
                        velocities(:,j) = velocities(:,j) + coeff_j * rij_hat
                    end if
                end do
            end do

        end subroutine collision_hard

    !----------Integradores--------------------------------------------!

        subroutine u1_propagator(t, r, v, box)
            implicit none
            real(dp), intent(in):: t, box
            real(dp), dimension(:,:), intent(inout):: r, v
            r(:,:) = r(:,:) + t * v(:,:)
            r(:,:) = r(:,:) - box * nint(r(:,:)/box)
        end subroutine u1_propagator

        subroutine u2_propagator(t, velocities, forces)
            implicit none
            real(dp), intent(in)                    :: t
            real(dp), dimension(:,:), intent(inout) :: velocities
            real(dp), dimension(:,:), intent(in)    :: forces
            
            integer :: i, n
            
            n = size(velocities,2) 
            do i = 1, n
                velocities(:,i) = velocities(:,i) + t * forces(:,i)
            end do
        end subroutine u2_propagator
        
        subroutine langevin_step(dt, velocities, masses, gamma, temp)
            implicit none
            real(dp), intent(in) :: dt, gamma, temp
            real(dp), dimension(:,:), intent(inout) :: velocities
            real(dp), dimension(:), intent(in) :: masses

            integer :: i, d, n
            real(dp) :: c1, c2
            real(dp) :: gauss

            n = size(velocities,2)

            c1 = exp(-gamma * dt)

            do i = 1, n
                c2 = sqrt((1.0_dp - c1*c1) * temp / masses(i))

                do d = 1, 3
                    call random_gaussian(gauss)
                    velocities(d,i) = c1 * velocities(d,i) + c2 * gauss
                end do
            end do

        end subroutine langevin_step
end module handler_mod