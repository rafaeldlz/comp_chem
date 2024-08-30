module make_histogram_mod
      implicit none
      public
      contains

      subroutine make_histogram(velocity,velocity_sq,set_min,set_max,n_bin,bin_width,histo_velocity,histo_velx,bin_index)
              real :: velocity(3,1000)
              real :: velocity_sq(1000)
              real :: set_min, set_max, bin_width
              integer :: n_bin, histo_velocity(100), histo_velx(100), bin_index, i

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
              end do

              ! histogram for Velocities
              do i = 1, 1000
               bin_index = floor((velocity_sq(i)-set_min)/bin_width) + 1
               if (bin_index .le. n_bin) histo_velocity(bin_index) = histo_velocity(bin_index) + 1
              end do

              do i = 1, n_bin
               write (14,fmt="(I6,2x,F8.2,2x,F8.2,2x,F8.2,2x,I6,2x,I6)") i,(set_min+(i-1)*bin_width),(set_min+i*bin_width), &
                &(set_min+(i-0.5)*bin_width),histo_velocity(i),histo_velx(i)
              end do
      end subroutine
end module make_histogram_mod

