program main_program
implicit none

integer :: n2, n_grid, a, i, k, j, l, n
integer, allocatable, dimension(:) :: combination

write(*,*) 'm out of n?'
read(*,*) n2
read(*,*) n_grid

allocate ( combination(n2) )
open(unit=4, file='combination.txt', status='replace', action='write')
!
do i = 1, n2
   combination(i) = i !123
end do
write(4,*) combination
!
!--------------------------------------------------
!now, we generate all of the other combinations:
!--------------------------------------------------
!
!
a = 2
k = n2 + 1 !drives the increase of the last integer in the
!          combination, combination(m)
!
do while (combination(1) .lt. (n_grid - n2 + 1))
   if (combination(n2) .lt. n_grid) then
      combination(n2) = k !124 125 135 235
      k = k + 1
      !
      write(4,*) combination
      write(*,*) 'loop1', a, combination
      !   
   else if (combination(n2) .eq. n_grid) then
      do j = 1, n2
         l = n2 - j
         if ( combination(j) .eq. (n_grid - l) ) then
            combination(j-1) = combination(j-1) + 1
            do n = j, n2
               combination(n) = combination(n-1) + 1 !134 145 234 245 
                                                                !345
            end do
            k = combination(n2) + 1
            !
            write(4,*) combination
            write(*,*) 'loop2', a, combination
            !
            exit
         end if
      end do
   end if
   a = a + 1
end do 


end program main_program
