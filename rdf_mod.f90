module rdf_mod
      implicit none
      public
      contains

      subroutine rdf(num_atoms,atoms,histogram,length,t,sigma,g,num_steps)
              integer :: num_atoms, num_steps
              real :: atoms(:,:), length, sigma
              real, allocatable :: atoms_rdf(:,:)
              integer :: t, i, k, j, nk
              real :: histogram(:)                  ! Histogram for RDF
              real :: g(:)                          ! RDF
              real :: dr, r1, r2, density_sigma
              real :: rij(3)
              real :: rij_sq,r_hi,r_lo,h_id,const,rho
              real, parameter :: pi = 4.0*atan(1.0)

              allocate(atoms_rdf(3,num_atoms))

              dr = 0.05
              dr = dr * sigma / length
              nk = floor(0.5*length/dr)
              !atoms_rdf = atoms/length

              do i = 1, num_atoms
               do j = 1, num_atoms
                atoms_rdf(i,j) = atoms(i,j)/length
               end do
              end do

              rij = 0.0

              print *, "nk=" , nk

              do i = 1, (num_atoms-1)
               do j = (i+1), num_atoms
                rij(:) = atoms_rdf(:,i) - atoms_rdf(:,j)
                rij(:) = rij(:) - anint(rij(:))
                rij_sq = sum(rij**2)
                k = floor(sqrt(rij_sq)/dr) + 1
                if (k.le.nk) histogram(k) = histogram(k) + 2
               end do
              end do

              !check histogram
              open (17,file="hist-rdf.dat")
              do k = 1 , nk
               write (17,*) k, histogram(k)
              end do

              if (t.eq.num_steps) then
               print *, "Normalize rdf"
               rho = real(num_atoms)
               const = 4.0 * pi * rho / 3.0
               do k = 1, nk
                g(k) = histogram(k) / real(num_atoms*num_steps) ! average number
                r_lo = real(k-1) * dr
                r_hi = r_lo + dr
                h_id = const * (r_hi**3 - r_lo**3) ! ideal number
                g(k) = g(k) / h_id
               end do

               print *, "Done calculating RDF"
               print *, "Printing to output"

               dr = dr * length

               open (unit=16,file="rdf.dat")
               do k = 1, nk
                write (16, *) (real(k)-0.5)*dr, g(k)
               end do
               close (16)

              end if

              deallocate(atoms_rdf)
      end subroutine

end module rdf_mod
