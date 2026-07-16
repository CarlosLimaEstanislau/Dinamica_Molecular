    module forces_mod
        use, intrinsic:: iso_fortran_env, only: dp => real64
        implicit none

        contains
            subroutine forces(positions, box, radius, p_tiu, delta, epsilon, bjl, kappa, r_cut, list, point, force, potential)
                implicit none

                real(dp), dimension(:,:), intent(in) :: positions
                real(dp), dimension(:), intent(in) :: p_tiu, radius
                integer, dimension(:), intent(in)  :: list, point
                real(dp), intent(in) :: kappa, box, epsilon, bjl, r_cut, delta
                real(dp), intent(out) :: potential
                real(dp), dimension(:,:), intent(out) :: force

                real(dp), dimension(3) :: rij, rij_hat
                integer :: i, j, n, k
                real(dp) :: rc, rc2, sigma, a1, a2, r2, r
                real(dp) :: dbeta_dr, dgamma_dr,dalpha_dr, surf_dist, sigmaLJ, gamma, beta, alpha
                real(dp) :: expfac, du_dx, du_dy, du_dz, du_dr, u
                real(dp) :: sr2, sr6, sr12, force_scalar_wca, pot_wca

                potential = 0.0_dp
                force = 0.0_dp

                n = size(positions, 2)

                do i = 1, n-1
                    do k = point(i), point(i+1) - 1

                        j = list(k)

                        rij = positions(:,i) - positions(:,j)
                        rij = rij - box * nint(rij / box)
                        r2 = dot_product(rij, rij)
                        r = sqrt(r2)
                        rij_hat = rij/r

                        sigma = radius(i) + radius(j)

                        surf_dist = r - sigma

                        sigmaLJ = sigma / 2.0_dp**(1.0_dp/6.0_dp)
                        
                        if (surf_dist > 0.0_dp .and. surf_dist <= r_cut) then

                            alpha = p_tiu(i)*p_tiu(j)*0.5_dp*cos(delta)*exp(-kappa*r)/(r/bjl)**3

                            beta = -(kappa*r + 1) + kappa*kappa*r*r + 3.0_dp*(kappa*r +1)

                            gamma = (rij(3)*rij(3))/(r*r)

                            dalpha_dr = -alpha * (kappa + 3.0_dp/r)

                            dbeta_dr  = 2.0_dp*(kappa*kappa*r + kappa)

                            dgamma_dr = -2.0_dp*(gamma/r) 

                            du_dr = dalpha_dr*beta*gamma + alpha*gamma*dbeta_dr + alpha*beta*dgamma_dr

                            du_dx = du_dr*(rij(1)/r)
                            du_dy = du_dr*(rij(2)/r)
                            du_dz = (du_dr + 2.0_dp*alpha*beta*(1.0_dp-gamma)/r)*(rij(3)/r)


                            u = alpha*beta*gamma

                            force(1,i) = force(1,i) - du_dx
                            force(2,i) = force(2,i) - du_dy
                            force(3,i) = force(3,i) - du_dz

                            force(1,j) = force(1,j) + du_dx
                            force(2,j) = force(2,j) + du_dy
                            force(3,j) = force(3,j) + du_dz
                            
                            potential = potential + u
                        end if

                        if (surf_dist <= 0.0_dp) then
                            sr2 = (sigmaLJ * sigmaLJ) / r2
                            sr6 = sr2 * sr2 * sr2
                            sr12 = sr6 * sr6
                            
                            pot_wca = 4.0_dp * epsilon * (sr12 - sr6) + epsilon
                            force_scalar_wca = 24.0_dp * epsilon * (2.0_dp * sr12 - sr6) / r2
                            
                            potential = potential + pot_wca
                            force(:,i) = force(:,i) + force_scalar_wca * rij
                            force(:,j) = force(:,j) - force_scalar_wca * rij
                        end if

                    end do
                end do
            end subroutine forces

            subroutine verlet_list(positions, box, r_list, list, point, r_save)
                implicit none

                real(dp), dimension(:,:), intent(in)    :: positions
                real(dp), intent(in)                    :: box, r_list
                integer, dimension(:), allocatable, intent(inout) :: list
                integer, dimension(:), intent(inout)    :: point
                real(dp), dimension(:,:), intent(inout) :: r_save

                integer :: i, j, k, n
                real(dp) :: r_list2, rij2
                real(dp), dimension(3) :: rij

                n = size(positions, 2)
                r_list2 = r_list * r_list

                k = 0
                point(1) = 1

                do i = 1, n - 1

                    do j = i + 1, n
                        rij = positions(:,i) - positions(:,j)
                        rij = rij - box * nint(rij / box)
                        rij2 = dot_product(rij, rij)

                        if (rij2 <= r_list2) then
                            k = k + 1
                            if (k > size(list)) then
                                stop "Erro: The Verlet list is too short!"
                            end if
                            list(k) = j
                        end if
                    end do
                    point(i+1) = k + 1
                end do
                point(n+1) = k + 1

                r_save = positions
            end subroutine verlet_list

            logical function check_list(positions, r_save, box, skin)
                implicit none

                real(dp), dimension(:,:), intent(in) :: positions, r_save
                real(dp), intent(in) :: box, skin

                integer :: n
                real(dp) :: dr2_max
                real(dp), dimension(3, size(positions,2)) :: dr

                n = size(positions, 2)

                dr = positions - r_save
                dr = dr - box * nint(dr / box)

                dr2_max = maxval(sum(dr**2, dim=1))

                check_list = (4.0_dp * dr2_max >= skin*skin)
            end function check_list
    end module forces_mod