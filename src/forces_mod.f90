module forces_mod
    use, intrinsic:: iso_fortran_env, only: dp => real64
    implicit none

    contains
        subroutine DLVO(charges, positions, box, radius, kappa, force, potential)
            implicit none

            real(dp), dimension(:,:), intent(in) :: positions
            real(dp), dimension(:), intent(in) :: radius, charges
            real(dp), intent(in) :: kappa, box
            real(dp), intent(out) :: potential
            real(dp), dimension(:,:), intent(out) :: force

            integer :: i, j, n
            real(dp), dimension(3) :: rij
            real(dp) :: rij2, rij_norm
            real(dp) :: rc, sigma
            real(dp) :: Bi, Bj, A, bjerrum
            real(dp) :: expfac, force_scalar

            potential = 0.0_dp
            force = 0.0_dp
            bjerrum = 0.0007

            n = size(positions,2)

            do i = 1, n-1
                do j = i+1, n

                    rij = positions(:,i) - positions(:,j)
                    rij = rij - box * nint(rij / box)
                    rij2 = dot_product(rij, rij)
                    rij_norm = sqrt(rij2)

                    sigma = radius(i) + radius(j)
                    rc = 2.0_dp**(1.0_dp/6.0_dp) * sigma
                    
                    if (rij_norm .gt. rc) then
                        Bi = exp(kappa * radius(i)) / (1.0_dp + kappa * radius(i))
                        Bj = exp(kappa * radius(j)) / (1.0_dp + kappa * radius(j))

                        A = charges(i) * charges(j) * Bi * Bj * bjerrum
                        expfac = exp(-kappa * rij_norm)

                        potential = potential + A * expfac / rij_norm

                        force_scalar = A * expfac * (kappa / rij2 + 1.0_dp / (rij2 * rij_norm))

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
            real(dp), dimension(:,:), intent(inout) :: force
            real(dp), intent(inout) :: potential

            integer :: i, j, n
            real(dp), dimension(3) :: rij
            real(dp) :: r2, r, sigma, rc, rc2
            real(dp) :: sr2, sr6, force_scalar, pot

            force = 0.0_dp
            potential = 0.0_dp

            n = size(positions,2)

            do i = 1, n-1
                do j = i+1, n

                    rij = positions(:,i) - positions(:,j)
                    rij = rij - box * nint(rij / box)
                    r2 = dot_product(rij, rij)

                    sigma = radius(i) + radius(j)
                    rc = 2.0_dp**(1.0_dp/6.0_dp) * sigma
                    rc2 = rc*rc

                    if (r2 < rc2) then

                        sr2 = (sigma*sigma) / r2
                        sr6 = sr2**3

                        force_scalar = 24.0_dp * epsilon * (2.0_dp*sr6*sr6 - sr6) / r2

                        force(:,i) = force(:,i) + force_scalar * rij
                        force(:,j) = force(:,j) - force_scalar * rij

                        pot = 4.0_dp * epsilon * (sr6*sr6 - sr6) + epsilon
                        potential = potential + pot
                    end if

                end do
            end do
        end subroutine WCA
end module forces_mod