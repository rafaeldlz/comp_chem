module force_energy_mod
      implicit none
      public
      contains
      subroutine force_energy(num_atoms,atoms,length,force,energy_p)
        real :: atoms(3,1000),force(3,1000),energy_p,length
        real :: rij(3), sr2,sr6,sr12,rij_sq,distance
        real, parameter :: sigma = 3.40E-10      ! meters
        real, parameter :: avogadro = 6.022E23
        real :: disp_e                           ! kcal/mol
        integer :: i,j,num_atoms
        disp_e = 0.997 * 4184 / avogadro         ! convert from kcal/mol to j/atom

        force = 0.0

        do i = 1, num_atoms
         do j = 1, num_atoms

         ! Only work on unique pairs of atoms
         if (i.eq.j) cycle

         ! Calculate the distance between the atoms, using PBC
         rij(:) = atoms(:,i) - atoms(:,j)
         rij(:) = rij(:) - (length*(anint(rij(:)/length)))
         rij_sq = sum(rij(:)**2)
         distance = sqrt(rij_sq)

         ! Calculate the forces between the atoms
         sr2 = (sigma**2)/rij_sq
         sr6 = sr2**3
         force(:,i) = force(:,i) + ((24*disp_e/rij_sq)*((2 * sr6**2) - sr6) * distance * (rij(:)/distance))

         ! Calculate the potential between the atoms
         sr12 = sr6*sr6
         energy_p = energy_p +( 4.0 * disp_e * (sr12 - sr6))

         end do
        end do

      end subroutine force_energy
end module force_energy_mod
