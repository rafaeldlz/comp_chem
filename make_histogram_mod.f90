module make_histogram_mod
      implicit none

      private
      public make_histogram

      contains

      subroutine make_histogram(num_atoms,velocity,velocity_sq)
              integer, intent(in) :: num_atoms
              real, intent(in) :: velocity(:,:)
              real, intent(in) :: velocity_sq(:)
              real :: set_min, set_max, bin_width
              integer :: n_bin, histo_velocity(100), histo_velx(100), bin_index
              integer :: i

              ! Set range of velocities
              set_min = 800.0
              set_max = -800.0
              
              ! Set number of bins
              n_bin = 100
              
              ! Set bin width
              bin_width = (set_max - set_min)/n_bin
              histo_velocity = 0
              histo_velx = 0

              ! Histogram for x-component
              do i = 1, num_atoms
               bin_index = floor((velocity(1,i)-set_min)/bin_width) + 1
               if (bin_index .le. n_bin) histo_velx(bin_index) = histo_velx(bin_index) + 1
              end do

              ! Histogram for velocities
              do i = 1, num_atoms
               bin_index = floor((velocity_sq(i)-set_min)/bin_width) + 1
               if (bin_index .le. n_bin) histo_velocity(bin_index) = histo_velocity(bin_index) + 1
              end do

              ! For each bin, write the number, minimum, maximum, midpoint, number of velocities,
              ! and number of velocity x-components
              do i = 1, n_bin
               write (14,fmt="(I6,2x,F8.2,2x,F8.2,2x,F8.2,2x,I6,2x,I6)") i,(set_min+(i-1)*bin_width),(set_min+i*bin_width), &
                &(set_min+(i-0.5)*bin_width),histo_velocity(i),histo_velx(i)
              end do
      end subroutine

end module make_histogram_mod

