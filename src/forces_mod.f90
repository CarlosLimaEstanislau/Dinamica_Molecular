module forces_mod
        use, intrinsic:: iso_fortran_env, only: dp => real64
        implicit none

        contains
            subroutine DLVO(charges, positions, box, radius,a0, kappa, force, potential)
                implicit none

                real(dp), dimension(:,:), intent(in) :: positions
                real(dp), dimension(:), intent(in) :: radius, charges
                real(dp), intent(in) :: kappa, box, a0
                real(dp), intent(out) :: potential
                real(dp), dimension(:,:), intent(out) :: force

                integer :: i, j, n
                real(dp), dimension(3) :: rij
                real(dp) :: rc, sigma, a1, a2, rij2, r
                real(dp) :: Bi, Bj, A, bjerrum
                real(dp) :: expfac, force_scalar, pot

                potential = 0.0_dp
                force = 0.0_dp
                force_scalar = 0.0_dp
                bjerrum = 0.00072_dp/a0    !Com o comprimento de bjerrum em micrometros

                n = size(positions,2) 

                do i = 1, n-1
                    do j = i+1, n

                        rij = positions(:,i) - positions(:,j)
                        rij = rij - box * nint(rij / box)
                        rij2 = dot_product(rij, rij)
                        r = sqrt(rij2)

                        sigma = radius(i) + radius(j)
                        
                        a1 = radius(i)
                        a2 = radius(j)

                        rc = 2.0_dp**(1.0_dp/6.0_dp) * sigma
                        
                        if (r .gt. rc) then
                            Bi = exp(kappa*a1) / (1.0_dp + kappa*a1)
                            Bj = exp(kappa*a2) / (1.0_dp + kappa*a2)
                            A = charges(i) * charges(j) * Bi * Bj * bjerrum
                            expfac = exp(-kappa*r)

                            pot = A * expfac / r
                            force_scalar = A * expfac * (kappa*r + 1.0_dp)/ (rij2*r)

                            potential = potential + pot
                            
                            force(:,i) = force(:,i) + force_scalar * rij
                            force(:,j) = force(:,j) - force_scalar * rij
                        end if
                    end do
                end do
            end subroutine DLVO

            subroutine WCA(positions, box, radius, epsilon, force, potential)
                implicit none
                
                real(dp), intent(in) :: box, epsilon
                real(dp), dimension(:,:), intent(in) :: positions
                real(dp), dimension(:), intent(in) :: radius
                real(dp), dimension(:,:), intent(out) :: force
                real(dp), intent(out) :: potential

                integer :: i, j, n
                real(dp), dimension(3) :: rij
                real(dp) :: rij_norm, sigma, rc, rc2, rij_sq
                real(dp) :: sr2, sr6, sr12, force_scalar, pot, r

                force = 0.0_dp
                potential = 0.0_dp
                pot = 0.0_dp

                n = size(positions,2)

                do i = 1, n-1
                    do j = i+1, n

                        rij = positions(:,i) - positions(:,j)
                        rij = rij - box * nint(rij / box)
                        rij_sq = dot_product(rij, rij)
                        r = sqrt(rij_sq)

                        sigma = radius(i) + radius(j)

                        rc = 2.0_dp**(1.0_dp/6.0_dp) * sigma
                        rc2 = rc*rc

                        if (rij_sq < rc2) then
                            
                            sr2 = (sigma*sigma)/rij_sq
                            sr6 = sr2*sr2*sr2
                            sr12 = sr6*sr6
                            
                            pot = 4.0_dp * epsilon * (sr12 - sr6) + epsilon
                            force_scalar = 24.0_dp * epsilon * (2.0_dp*sr12 - sr6) / rij_sq
                            
                            potential = potential + pot 
                            
                            force(:,i) = force(:,i) + force_scalar * rij
                            force(:,j) = force(:,j) - force_scalar * rij
                        end if
                    end do
                end do
            end subroutine WCA
    end module forces_mod