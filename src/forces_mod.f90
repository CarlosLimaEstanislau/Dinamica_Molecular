    module forces_mod
        use, intrinsic:: iso_fortran_env, only: dp => real64
        implicit none

        contains
            subroutine forces(positions, box, radius, charges, epsilon, bjl, kappa, r_cut_dlvo, list, point, force, potential)
                implicit none

                real(dp), dimension(:,:), intent(in) :: positions
                real(dp), dimension(:), intent(in) :: radius, charges
                integer, dimension(:), intent(in)  :: list, point
                real(dp), intent(in) :: kappa, box, epsilon, bjl, r_cut_dlvo
                real(dp), intent(out) :: potential
                real(dp), dimension(:,:), intent(out) :: force

                real(dp), dimension(3) :: rij
                integer :: i, j, n, k
                real(dp) :: sigma, a1, a2, r2, r
                real(dp) :: Bi, Bj, A, surf_dist, sigmaLJ
                real(dp) :: expfac, force_scalar_dlvo, pot_dlvo
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

                        sigma = radius(i) + radius(j)
                        a1 = radius(i)
                        a2 = radius(j)
                        
                        surf_dist = r - sigma

                        sigmaLJ = sigma / 2.0_dp**(1.0_dp/6.0_dp)
                        
                        if (surf_dist > 0.0_dp .and. surf_dist <= r_cut_dlvo) then
                            Bi = exp(kappa * a1) / (1.0_dp + kappa * a1)
                            Bj = exp(kappa * a2) / (1.0_dp + kappa * a2)
                            A = charges(i) * charges(j) * Bi * Bj * bjl
                            expfac = exp(-kappa * r)

                            pot_dlvo = A * expfac / r
                            force_scalar_dlvo = A * expfac * (kappa * r + 1.0_dp) / (r2 * r)
                            
                            potential = potential + pot_dlvo
                            force(:,i) = force(:,i) + force_scalar_dlvo * rij
                            force(:,j) = force(:,j) - force_scalar_dlvo * rij
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