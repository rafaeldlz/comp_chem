program liquid
implicit none

real :: disp_e                           ! kcal/mol
real, parameter :: sigma = 3.40E-10      ! meters
real ::  sigma_bk
real, parameter :: density = 1.395       ! g/cm3
real, parameter :: avogadro = 6.022E23
real, parameter :: mw = 40.0             ! g/mol
real :: length,r,potential
real, allocatable :: atoms(:,:)          ! positions
real, allocatable :: velocity(:,:)       ! velocities
real, allocatable :: velocity_sq(:)      ! sqrt of the sum of squares of velocity components
real, allocatable :: force(:,:),acceleration(:,:)
real :: box(200,200,200)
integer :: num_atoms = 0
integer :: i,j,k,n,m,t,grid,nk,nstep
real :: histogram(5000)                  ! Histogram for RDF
real :: g(5000)                          ! RDF
real :: dr, r1, r2, density_sigma
real :: rij(3)
real :: rij_sq,r_hi,r_lo,h_id,const,rho
real, parameter :: pi = 4.0*atan(1.0)
real, parameter :: temp = 85.0           ! Kelvin
real, parameter :: boltz_k = 1.380649E-23
real :: sr2, sr6, sr12
real :: energy_tot, energy_p, energy_k
real, parameter :: dt = 5.0E-15          ! femtoseconds
real :: total_time
real :: set_min, set_max, bin_width
integer :: n_bin, histo_velocity(100), histo_velx(100), bin_index
character, allocatable :: col1(:), col3(:), col4(:), col5(:)
integer, allocatable :: col2(:), col6(:)
real, allocatable :: col7(:), col8(:), col9(:)
character :: header
integer :: num = 3
integer, parameter :: num_steps = 100

! Input file
open (unit=10,file="water_box_1000.pdb")

! Output files
open (unit=11,file="liquid.xyz")
write (11,fmt="(A)") "1000"
open (unit=12,file="liquid_rdf.dat")
open (unit=13,file="energies.dat")
write (13, fmt="(4(A10,4x))") "Step","Total E","Potential","Kinetic"
open (unit=14, file="histo-velocity.dat")
open (unit=15, file="forces.dat")
open (unit=16, file="timestep-velocities.dat")

! Read the number of atoms
do
 read(10,*,iostat=i)
 if(i.ne.0) exit
 num_atoms = num_atoms + 1
end do
num_atoms = num_atoms - 6
print *, num_atoms
rewind(10)

allocate(atoms(3,num_atoms))
allocate(velocity(3,num_atoms))
allocate(velocity_sq(num_atoms))
allocate(force(3,num_atoms), acceleration(3,num_atoms))
allocate(col1(num_atoms), col2(num_atoms), col3(num_atoms))
allocate(col4(num_atoms), col5(num_atoms), col6(num_atoms))
allocate(col7(num_atoms), col8(num_atoms), col9(num_atoms))

disp_e = 0.997 * 4184 / avogadro

!!-----------------------------------------------------------------------v
! Calculate the length of the box
length = (num_atoms)*(avogadro**-1)*(mw)*(density**-1)*(1E-6)
length = length**(1.0/3.0)
!length = length*(1E10)
!length = length/sigma      ! Not using reduced units
write (11,fmt="(A,1x,E8.2,1x,A)") "Length of box:",length,"meters"
print *, length

grid = int(length)
!!-----------------------------------------------------------------------^

!!-----------------------------------------------------------------------v
! Generate initial coordinates (units = meters)
!n = 0
!do i = 1, 10
!  do j = 1, 10
!    do k = 1, 10
!      n = n + 1
!      atoms(1,n) = i*length*0.1
!      atoms(2,n) = j*length*0.1
!      atoms(3,n) = k*length*0.1
!      write (11,fmt=*) "Ar", atoms(:,n)
!    end do
!  end do
!end do
!!-----------------------------------------------------------------------^

! Generate initial coordinates (units = A)
! Skip header lines
do i = 1, 5
 read(10,*) header
end do
! Get the coordinates of the first num = 3 atoms
do i = 1, num_atoms
 read(10,*) col1(i), col2(i), col3(i), col4(i), col5(i), col6(i), col7(i), col8(i), col9(i)
 atoms(1,i) = col7(i)
 atoms(2,i) = col8(i)
 atoms(3,i) = col9(i)
! Shift the x axis by 10 A
! atoms(1,i) = atoms(1,i) + 10
! Write to output
 if (mod(i,3).eq.0) then
  write(11,*) "O", atoms(:,i)
 else
  write(11,*) "H", atoms(:,i)
 end if
end do

!!-----------------------------------------------------------------------v
! Initial velocity distribution
print *, "Initial velocity distribution"
sigma_bk = sqrt(boltz_k * temp / (mw*(1E-3)/avogadro))

do i=1,500 ! Generate 500 pairs each x y z
!print *, "getting velocities", i
 call get_rand(r,n)
 r1 = r
 call get_rand(r,n)
 r2 = r
 velocity(1,i)     = sigma_bk * sqrt(-2.0*LOG(r1)) * cos(2.0*pi*r2)
 velocity(1,i+500) = sigma_bk * sqrt(-2.0*LOG(r1)) * sin(2.0*pi*r2)

 call get_rand(r,n)
 r1 = r
 call get_rand(r,n)
 r2 = r
 velocity(2,i)     = sigma_bk * sqrt(-2.0*LOG(r1)) * cos(2.0*pi*r2)
 velocity(2,i+500) = sigma_bk * sqrt(-2.0*LOG(r1)) * sin(2.0*pi*r2)

 call get_rand(r,n)
 r1 = r
 call get_rand(r,n)
 r2 = r
 velocity(3,i)     = sigma_bk * sqrt(-2.0*LOG(r1)) * cos(2.0*pi*r2)
 velocity(3,i+500) = sigma_bk * sqrt(-2.0*LOG(r1)) * sin(2.0*pi*r2)

end do ! End generating initial velocities

velocity_sq(:) = sqrt(velocity(1,:)**2 + velocity(2,:)**2 + velocity(3,:)**2) 

! Make histogram
set_min = 800.0
set_max = -800.0
n_bin = 100
bin_width = (set_max - set_min)/n_bin
histo_velocity = 0
histo_velx = 0

! histogram for x-component
do i = 1, 1000
 bin_index = floor((velocity(1,i)-set_min)/bin_width) + 1
 if (bin_index .le. n_bin) histo_velx(bin_index) = histo_velx(bin_index) + 1
end do ! End one velocity component loop

! histogram for Velocities
do i = 1, 1000
 bin_index = floor((velocity_sq(i)-set_min)/bin_width) + 1
 if (bin_index .le. n_bin) histo_velocity(bin_index) = histo_velocity(bin_index) + 1
end do ! End one velocity component loop


do i = 1, n_bin
 write (14,fmt="(I6,2x,F8.2,2x,F8.2,2x,F8.2,2x,I6,2x,I6)") i,(set_min+(i-1)*bin_width),(set_min+i*bin_width), &
  &(set_min+(i-0.5)*bin_width),histo_velocity(i),histo_velx(i)
end do

!!-----------------------------------------------------------------------^

!!-----------------------------------------------------------------------v
! Calculate initial forces and energies
print *, "length", length
energy_p = 0.0
energy_tot = 0.0
energy_k = 0.0

call force_energy(num_atoms,atoms,length,force,energy_p)
print *, "start change in v", force(:,i) * dt * 0.5 /  (mw*(1E-3)/avogadro)

do i = 1, num_atoms
 energy_k = energy_k + (0.5*(mw*(1E-3)/avogadro)*sum(velocity(:,i)**2))
end do

energy_tot = energy_p + energy_k
write (13,fmt="(A4,5x,3(E10.4,3x))") "0", energy_tot, energy_p, energy_k

write (15,*) "Timestep = 0"
do i=1,num_atoms
 write (15,*) force(:,i)
end do
!!-----------------------------------------------------------------------^

!!-----------------------------------------------------------------------v
! Velocity Verlet

print *, "Starting Velocity Verlet"

total_time = 0.0
histogram = 0.0
g = 0.0

do t=1,num_steps
 print *, "Timestep #", t

 acceleration = 0.0

 ! Move the atoms
 do i = 1, num_atoms
  acceleration(:,i) = force(:,i) / (mw*(1E-3)/avogadro)
  velocity(:,i)   = velocity(:,i) + (0.5 * dt * acceleration(:,i))
  atoms(:,i)      = atoms(:,i) + (dt * velocity(:,i))
 end do

 ! PBC
 do i = 1, num_atoms
  do j = 1, 3
   if (atoms(j,i).gt.length) atoms(j,i) = atoms(j,i) - length
   if (atoms(j,i).lt.0.0) atoms(j,i) = atoms(j,i) + length
  end do
 end do

 ! Check velocities
 write (16,*) t
 do i = 1, num_atoms
  write (16,*) velocity(:,i)
 end do

 ! Call force subroutine
 energy_p = 0.0
 call force_energy(num_atoms,atoms,length,force,energy_p)

 ! Check forces
 write (15,*) "Timestep",t
 do i=1,num_atoms
  write (15,*) force(:,i)
 end do

! print *, "Update velocity. Full"
 do i = 1, num_atoms
  acceleration(:,i) = force(:,i) / (mw*(1E-3)/avogadro)
  velocity(:,i)   = velocity(:,i) + (0.5 * dt * acceleration(:,i))
 end do

 ! Check velocities
 write (16,*) t
 do i = 1, num_atoms
  write (16,*) velocity(:,i)
 end do

 ! Calculate total energy
 energy_k = 0.0
 do i = 1, num_atoms
  energy_k = energy_k + (0.5*(mw*(1E-3)/avogadro)*sum(velocity(:,i)**2))
 end do
 
 energy_tot = energy_p + energy_k
 write (13,fmt="(I4,5x,3(E10.4,3x))") t, energy_tot, energy_p, energy_k
 
 ! Write new positions to output
 write (11,fmt="(A)") "1000"
 write (11,*) "Timestep #", t
 do i = 1, num_atoms
  write (11,*) "Ar", atoms(:,i)
 end do
 
 ! Increment total time for reference
! total_time = total_time + dt
 
 ! Do RDF for the last 1000 steps
 if (t.ge.num_steps-1000) call rdf(num_atoms,atoms,length,t,sigma,g,num_steps)
 
end do ! End entire do loop

!!-----------------------------------------------------------------------v

! Close output files
do i=11,16
 close (i)
end do

deallocate(atoms)
deallocate(velocity)
deallocate(velocity_sq)
deallocate(force, acceleration)
deallocate(col1, col2, col3)
deallocate(col4, col5, col6)
deallocate(col7, col8, col9)

end program liquid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_rand(r,n)
implicit none
integer, allocatable :: seed(:)
integer :: n
real ,intent(out) :: r

call random_seed(size=n)
allocate(seed(n))
call random_seed(get=seed)
call random_number(r)
deallocate(seed)
return
end subroutine get_rand

!------------------------------------------------------------------------

subroutine force_energy(num_atoms,atoms,length,force,energy_p)
implicit none
real :: atoms(3,1000),force(3,1000),energy_p,length
real :: rij(3), sr2,sr6,sr12,rij_sq,distance
real, parameter :: sigma = 3.40E-10      ! meters
real, parameter :: avogadro = 6.022E23  
real :: disp_e        ! kcal/mol
integer :: i,j,num_atoms
disp_e = 0.997 * 4184 / avogadro ! convert from kcal/mol to j/atom

force = 0.0

do i = 1, num_atoms
 do j = 1, num_atoms
  
  if (i.eq.j) cycle

  ! Calculate distance between 2 atoms
  rij(:) = atoms(:,i) - atoms(:,j) 
  rij(:) = rij(:) - (length*(anint(rij(:)/length)))
  rij_sq = sum(rij(:)**2)
  distance = sqrt(rij_sq)

  ! Calculate forces
  sr2 = (sigma**2)/rij_sq
  sr6 = sr2**3
  force(:,i) = force(:,i) + ((24*disp_e/rij_sq)*((2 * sr6**2) - sr6) * distance * (rij(:)/distance))

! if (i == 1) then
!   print *, "in subroutine: ", rij(:), "atom", j, force(:,i)
! end if

  ! Calculate potential between 2 atoms
  sr12 = sr6*sr6
  energy_p = energy_p +( 4.0 * disp_e * (sr12 - sr6))

 end do
end do

!print *, "in subroutine: ", force(:,1)

end subroutine force_energy

!------------------------------------------------------------------------

subroutine rdf(num_atoms,atoms,length,t,sigma,g,num_steps)
implicit none
integer :: num_atoms, num_steps
real :: atoms(3,1000), length, sigma, atoms_rdf(3,1000)
integer :: t, i, k, j, nk
real :: histogram(5000)                  ! Histogram for RDF
real :: g(5000)                          ! RDF
real :: dr, r1, r2, density_sigma
real :: rij(3)
real :: rij_sq,r_hi,r_lo,h_id,const,rho
real, parameter :: pi = 4.0*atan(1.0)

dr = 0.05
dr = dr * sigma / length
nk = floor(0.5*length/dr)
!nk = 500
atoms_rdf = atoms/length
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
end subroutine
