module get_rand_mod
      implicit none
      public
      contains
      subroutine get_rand(r,n)
              integer, allocatable :: seed(:)
              integer :: n
              real ,intent(out) :: r

              ! Assign a random value to n, and set seed to size n
              call random_seed(size=n)
              allocate(seed(n))
              call random_seed(get=seed)
              
              ! Assign a random value to r
              call random_number(r)
              deallocate(seed)
              return  
      end subroutine get_rand
end module get_rand_mod
