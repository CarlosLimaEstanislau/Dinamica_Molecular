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

        subroutine shuffle_array(arr)
            implicit none
            real(dp), intent(inout) :: arr(:)
            integer :: i, j, n
            real(dp) :: r, temp
            
            n = size(arr)
            
            do i = n, 2, -1
                call random_number(r)
                j = floor(r * i) + 1
                temp = arr(i)
                arr(i) = arr(j)
                arr(j) = temp
            end do
        end subroutine shuffle_array

    !----------Integradores--------------------------------------------!

        subroutine u1_propagator(t, r, v, box)
            implicit none
            real(dp), intent(in):: t, box
            real(dp), dimension(:,:), intent(inout):: r, v
            r(:,:) = r(:,:) + t * v(:,:)
            r(:,:) = r(:,:) - box * nint(r(:,:)/box)
        end subroutine u1_propagator

        subroutine u2_propagator(t, velocities, forces, masses)
            implicit none
            real(dp), intent(in)                    :: t
            real(dp), dimension(:,:), intent(inout) :: velocities
            real(dp), dimension(:,:), intent(in)    :: forces
            real(dp), dimension(:), intent(in)      :: masses
            
            integer :: i, n
            
            n = size(velocities,2) 
            do i = 1, n
                velocities(:,i) = velocities(:,i) + t * (forces(:,i)/masses(i))
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