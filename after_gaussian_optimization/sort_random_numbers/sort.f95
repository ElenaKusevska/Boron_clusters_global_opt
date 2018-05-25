program sorting
use sorting_subroutines
implicit none

real(kind=8), allocatable, dimension(:) :: A1, A2
real(kind=8) :: B
integer :: n, i

write(*,*) 'size of array?'
read(*,*) n

allocate (A1(n), A2(n))

call random_seed()
do i = 1, n
   call random_number(B)
   A1(i) = B*500.00
   A2(i) = A1(i)
end do

call bubble_sort(n, A1)

call first_and_last(n, A2)

end program sorting
