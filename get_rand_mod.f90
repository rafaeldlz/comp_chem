module get_rand_mod
      implicit none
      public
      contains
      subroutine get_rand(r,n)
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
end module get_rand_mod
