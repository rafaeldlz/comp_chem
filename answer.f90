program liquid
implicit none

double precision, parameter :: disp_e = 0.997          ! kcal/mol
double precision, parameter :: sigma = 3.40            ! angstroms
double precision, parameter :: density = 1.395         ! g/cm3
double precision, parameter :: avogadro = 6.022D23  
double precision, parameter :: mw = 40.0               ! g/mol
double precision, parameter :: temp = 85.0             ! Kelvin
real,             parameter :: dt = 1E-15              ! femtoseconds
double precision            :: length,r,potential
real                        :: atoms(3,1000) ! dim1:x,y,z ; dim2:atom index
real                        :: box(200,200,200)
integer,          parameter :: num_atoms = 1000
integer                     :: i,j,k,n,m,grid,nk,nstep,iter
real                        :: histogram(5000)
real                        :: g(5000) ! RDF
real                        :: dr,total_time,theta 
real                        :: rij(3)
real                        :: rij_sq,r_hi,r_lo,h_id,const,rho
real, parameter             :: pi = 4.0*atan(1.0)
double precision            :: density_sigma
real                        :: acceleration(1000), force(1000)
real                        :: total_energy, pot_energy, kin_energy 

! Output files
open (unit=11,file="liquid_no_pbc.xyz")
!open (unit=11,file="liquid_yes_pbc.xyz")
write (11,fmt="(A)") "1000"
open (unit=12,file="liquid_no_pbc_rdf.dat")
!open (unit=12,file="liquid_yes_pbc_rdf.dat")

!!-----------------------------------------------------------------------v
! Calculate the length of the box
length = (num_atoms)*(avogadro**-1)*(mw)*(density**-1)*(1E-6)
length = length**(1.0/3.0)
length = length*(1E10)
length = length/sigma
write (11,fmt="(A,1x,F8.2,1x,A)") "Length of box:",length,"*sigma"
print *, length

grid = int(length)
!!-----------------------------------------------------------------------^

!!-----------------------------------------------------------------------v
! Generate initial velocities



!!-----------------------------------------------------------------------^

!!-----------------------------------------------------------------------v
! Generate initial coordinates (units = sigma)
do i=1,num_atoms
 k=0.0
  do while (k.eq.0.0)
   call get_rand(r,n)
   atoms(1,i)=r*length
   call get_rand(r,n)
   atoms(2,i)=r*length
   call get_rand(r,n)
   atoms(3,i)=r*length
   if ((atoms(1,i).eq.atoms(1,i-1)).and.(atoms(2,i).eq.atoms(2,i-1)).and.(atoms(3,i).eq.atoms(3,i-1))) then
    k=k+0
   else
    k=k+1
   end if
  end do
 write (11,fmt=*) "Ar",atoms(1,i), atoms(2,i), atoms(3,i)
end do ! end entire loop
!!-----------------------------------------------------------------------^

!!-----------------------------------------------------------------------v
! Velocity Verlet Algorithm
total_time = 0.0
do iter = 1, 1000
 ! Calculate atom distances, forces, angles
 
 ! Sweep through the atoms and calculate forces and acceleration
 do i = 1, 1000
  acceleration(i) = force(i) / mw
  ! new velocities
  ! new positions
 end do

 ! Calculate atom distances, forces, angles
 
 ! Sweep through the atoms and calculate forces and acceleration
 do i = 1, 1000
  acceleration(i) = force(i) / mw
  ! new velocities
  ! new positions
 end do
 
! Check the kinetic, potential, and total energies.

total_energy = pot_energy + kin_energy
total_time = total_time + dt
end do ! end time step
!!-----------------------------------------------------------------------^

!!-----------------------------------------------------------------------v
! Calculate the potential distance r for potential
potential = 0
do i=1,(num_atoms-1)
 do j=(i+1),num_atoms
  r=sum((atoms(:,i)-atoms(:,j))**2)
  potential=potential + (4.0*(((1.0/r)**12)-((1.0/r)**6)))
 end do
end do
potential=potential*disp_e

print *, "The value of the potential is", potential
!!-----------------------------------------------------------------------^

!!-----------------------------------------------------------------------v
! Code from Comp. Sim. of Liquids, #8.1
!print *, "Calculating RDF"
dr = 0.05
dr = dr / length
nk = floor(0.5*length/dr)
histogram = 0
atoms = atoms/length
print *, "nk=" , nk
do i = 1, (num_atoms-1)
 do j = (i+1), num_atoms
!  print *, "calculating distance"
  rij(:) = atoms(:,i) - atoms(:,j) 
!  rij(:) = rij(:) - (length*(anint(rij(:)/length)))
  rij_sq = sum(rij**2)
  k = floor(sqrt(rij_sq)/dr) + 1
  if (k.le.nk) histogram(k) = histogram(k) + 2
!  print *, "histogram(k) = ", histogram(k)
 end do
end do

!rho = density/sigma
rho = real(num_atoms)
const = 4.0 * pi * rho / 3.0 
do k = 1, nk
 g(k) = histogram(k) / real(num_atoms) ! average number
 r_lo = real(k-1) * dr
 r_hi = r_lo + dr
 h_id = const * (r_hi**3 - r_lo**3) ! ideal number
 g(k) = g(k) / h_id
! print *, "final g(k):", g(k)
end do

!print *, "Done calculating RDF"
!print *, "Printing to output"

dr = dr * length
do k = 1, nk
 write (12, *) (real(k)-0.5)*dr, g(k)
end do
!!-----------------------------------------------------------------------^


close (11)
close (12)

end program liquid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_rand(r,n)
implicit none
integer, allocatable :: seed(:)
integer :: n
double precision ,intent(out) :: r

call random_seed(size=n)
allocate(seed(n))
call random_seed(get=seed)
call random_number(r)
deallocate(seed)
return
end subroutine get_rand
