module seed
    ! Seed field
    !
    use parameters
    use time_grid
    use physical_grid
    use make_a_grid
    !  
      implicit none
    !
      double precision, dimension(nx) :: xseed
    contains
      subroutine init_random_seed
    !
        call random_number(xseed)
    !
      end subroutine init_random_seed
end module seed