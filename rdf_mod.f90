module rdf_mod
      implicit none

      private
      public rdf

      contains

      subroutine rdf(num_atoms,atoms,length,t,sigma,num_steps)
              integer, intent(in) :: num_atoms, num_steps
              real, intent(in) :: atoms(:,:), length, sigma
              real :: atoms_rdf(3,num_atoms)
              integer, intent(in) :: t
              integer :: i, j, k, nk
              real :: histogram(5000)        ! Histogram for RDF
              real :: g(5000)                ! RDF
              real :: dr, r1, r2, density_sigma
              real :: rij(3)
              real :: rij_sq,r_hi,r_lo,h_id,const,rho
              real, parameter :: pi = 4.0*atan(1.0)

              dr = 0.05
              dr = dr * sigma / length
              nk = floor(0.5*length/(sigma*dr))
              
              ! Rewrite the atomic coordinates in proportion to box length
              do i = 1, 3
               do j = 1, num_atoms
                atoms_rdf(i,j) = atoms(i,j)/length
               end do
              end do

              rij = 0.0
              histogram = 0.0
              g = 0.0

              ! Record the number of unique atom pairs
              ! whose vectors fall within each distance bin (x2)
              do i = 1, (num_atoms-1)
               do j = (i+1), num_atoms
                rij(:) = atoms_rdf(:,i) - atoms_rdf(:,j)
                rij(:) = rij(:) - anint(rij(:))
                rij_sq = sum(rij**2)
                k = floor(sqrt(rij_sq)/dr) + 1
                if (k.le.nk) histogram(k) = histogram(k) + 2
               end do
              end do

              ! Write the histogram
              open (17,file="hist-rdf.dat")
              do k = 1 , nk
               write (17,*) k, histogram(k)
              end do
                
              ! Calculate the rdf for each distance bin
              if (t.eq.num_steps) then
               print *, "Normalize rdf"
               rho = real(num_atoms)
               const = 4.0 * pi * rho / 3.0
               do k = 1, nk
                g(k) = histogram(k) / real(num_atoms*num_steps) ! Average number
                r_lo = real(k-1) * dr
                r_hi = r_lo + dr
                h_id = const * (r_hi**3 - r_lo**3) ! Ideal number
                g(k) = g(k) / h_id
               end do

               print *, "Done calculating RDF"
               print *, "Printing to output"

               dr = dr * length
                
               ! Write the rdf
               open (unit=16,file="rdf.dat")
               do k = 1, nk
                write (16, *) (real(k)-0.5)*dr, g(k)
               end do
               close (16)

              end if
      end subroutine

end module rdf_mod
