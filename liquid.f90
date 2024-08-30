program liquid
use get_rand_mod
use make_histogram_mod
use force_energy_mod
use rdf_mod
implicit none

real :: disp_e                           ! kcal/mol
real, parameter :: sigma = 3.40E-10      ! meters
real ::  sigma_bk 
real, parameter :: density = 1.395       ! g/cm3
real, parameter :: avogadro = 6.022E23  
real, parameter :: mw = 40.0             ! g/mol
real :: length,r,potential
real :: atoms(3,1000)                    ! positions
real :: velocity(3,1000)                 ! velocities
real :: velocity_sq(1000)                ! sqrt of the sum of squares of velocity components
real :: force(3,1000),acceleration(3,1000),box(200,200,200)
integer, parameter :: num_atoms = 1000
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
integer :: num_steps = 100

disp_e = 0.997 * 4184 / avogadro

! Output files
open (unit=11,file="liquid.xyz")
write (11,fmt="(A)") "1000"
open (unit=12,file="liquid_rdf.dat")
open (unit=13,file="energies.dat")
write (13, fmt="(4(A10,4x))") "Step","Total E","Potential","Kinetic"
open (unit=14, file="histo-velocity.dat")
open (unit=15, file="forces.dat")
open (unit=16, file="timestep-velocities.dat")

! Calculate the length of the box
length = (num_atoms)*(avogadro**-1)*(mw)*(density**-1)*(1E-6)
length = length**(1.0/3.0)
write (11,fmt="(A,1x,E8.2,1x,A)") "Length of box:",length,"meters"
print *, length

grid = int(length)

! Generate initial coordinates (units = meters)
n = 0
do i = 1, 10
  do j = 1, 10
    do k = 1, 10
      n = n + 1
      atoms(1,n) = i*length*0.1
      atoms(2,n) = j*length*0.1
      atoms(3,n) = k*length*0.1
      write (11,fmt=*) "Ar", atoms(:,i)
    end do
  end do
end do

! Initial velocity distribution
print *, "Initial velocity distribution"
sigma_bk = sqrt(boltz_k * temp / (mw*(1E-3)/avogadro))

! Generate 500 pairs each x y z

do i=1,500
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

end do

velocity_sq(:) = sqrt(velocity(1,:)**2 + velocity(2,:)**2 + velocity(3,:)**2) 

! Make histogram
call make_histogram(velocity,velocity_sq,set_min,set_max,n_bin,bin_width,histo_velocity,histo_velx,bin_index)

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
 
 ! Do RDF for the last 1000 steps
 if (t.ge.9000) call rdf(num_atoms,atoms,length,t,sigma,g,num_steps)
 
end do

! Close output files
do i=11,16
 close (i)
end do

end program liquid
