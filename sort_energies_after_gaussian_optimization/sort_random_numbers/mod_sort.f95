module sorting_subroutines
implicit none
contains

!------------------------------------------------------------------
!Subroutine CPU:
!------------------------------------------------------------------

subroutine CPU(start, m, steps, A, method)
implicit none

real(kind=8), intent(in) :: start
integer, intent(in) :: m, steps
real, allocatable, dimension(:), intent(in) :: A
character(len=*), intent(in) :: method
real(kind=8) :: finish, time
integer :: b
character a
character c
character format1

call cpu_time(finish)
time = finish - start

open(unit=2, file='CPU_time.txt', status='old', action='write', &
   position='append')
write(2,*) 'method: ', metehod, 'steps:', steps
write(2,*) 'CPU TIME:', time
write(2,*)
close(2)

if (v .gt. 999999999) then !because of a limitation in the way that
!                          the format is defined (b must contain a single
!                          digit, i.e. b <= 10
   write(*,*) 'sorry, matrix is too large'
   stop
end if

!Let's say the array contains 10 numbers:

write(a,'(I9)') m ! m = '10-------'
b = len_trim(a) !b = 2
write(c,'(A4I1A3)') '(A1I', b, 'A5)' ! c = '(A1I2A5)' 
write(format1, c) '(', m, 'F7.3)' ! 'format1 = '(24F7.3)------------'
format1 = trim(format1) ! '(24F7.3)'

write(*,format) A

end subroutine CPU

!------------------------------------------------------------------
!Bubble sort:
!------------------------------------------------------------------

subroutine bubble_sort (n1, numbers1)
implicit none

integer, intent(in) :: n1
real(kind=8), allocatable, dimension(:), intent(inout) :: numbers1
integer :: i, j, test
real(kind=8) :: start1, a, b

call cpu_time(start1)

open(unit=3, file='bubble.txt', status='replace', action='write')

!-----------------------------------------------------
!test = 1: they are ordered in descending order
!test = 0: they are not ordered in descending order
!-----------------------------------------------------

test = 0 !initialize the loop
j = 1 !counter for the number of steps in the while loop
k = 0 !variable to avoid looping the entire array when it's not necessary
do while (test == 0)
   write(3,*) '---------------------------------------------------'
   write(3,*) 'pass-bubble', j
   write(3,'(10000F7.2)') (numbers1(i), i = 1, n1)
   write(3,*) '---------------------------------------------------'
   test = 1 !maybe this time they will all be ordered
   do i = 1, n-k
      if ( numbers1(i) .lt. numbers1(i+1) ) then !A(4) < A5
         a = numbers1(i)
         b = numbers1(i+1)
         numbers1(i + 1) = a
         numbers1(i) = b !A(4) > A(5)
         test = 0 !but if they weren't ordered then, change test to 0
      end if
   end do
   j = j + 1 !increase j by one for the next loop
end do     

close(3)

call CPU(start1)

end subroutine bubble_sort

!subroutine insertion
!
!do i = 1, n
!   do j = i+1, n
!      if (A1(j) .gt. A1(i)) then
!         a = A1(i)
!         b = A1(j)
!         A1(i) = b
!         A1(j) = a
!         write(*,*) i, j
!         write(*,'(100F7.2)') A1
!      end if
!   end do
!end do
!
!end subroutine insertion
!
!subroutine selection
!
!do i = 1, n-1
!   minimum = A1(i)
!   k = 0
!   do j = 1 + i, n
!      if (A1(j) .lt. minimum) then
!         minimum = A1(j)
!         k = j
!      end if
!   end do
!   if (k /= 0) then
!      A1(k) = A1(i)
!      A1(i) = minimum
!      write(*,*) i, k
!      write(*,'(100F7.2)') A1
!   end if
!end do
!
!
!end subroutine selection
!
!-----------------------------------------------------------------
!First-and-last:
!-----------------------------------------------------------------
!
!subroutine first_and_last (n2, numbers2)
!implicit none
!
!integer, intent(in) :: n2
!real(kind=8), allocatable, dimension(:), intent(inout) :: numbers2
!real(kind=8), allocatable, dimension(:) :: work
!integer, dimension(1) :: k !for maxlox; maxval
!integer :: dim_work, i, j
!real(kind=8) :: start2, a, b
!
!call cpu_time(start2)
!
!open(unit=4, file='half-and-half.txt', status='replace', action='write')
!
!write(4,*)'-------------------------------------------------------'
!write(4,'(A13,100000F7.2)') 'input_array: ', (numbers2(i), i = 1, n2)
!write(4,*) '------------------------------------------------------'
!write(4,*)
!
!j = 0 !number of currently sorted values on both ends of the series
!dim_work = n2 !dimension of first instance of "working array"
!
!do while (dim_work .gt. 1) !untill both ends meet
!   allocate ( work(dim_work) )
!   do i = 1, dim_work !current dimension of "working array"
!      work(i) = numbers2(j + i) !populate the working array
!   end do
!
!   !---------------------
!   !maximum value:
!   !--------------------
!
!   k = maxloc(work) !the location of the i-th to last maximum value
!   write(4,*) '------------------------------------------------------'
!   write(4,*) 'k-max', k
!   write(4,'(A15,100000F7.2)') 'working array: ', (work(i), &
!      i = 1, dim_work)
!   write(4,*)
!
!   a = numbers2(j + k(1)) !the i-the to first maximum is the (i+k)th element
!   b = numbers2(j + 1) !and it should be switched with the i-th to 
!                         first element:
!   numbers2(j + 1) = a
!   numbers2(j + k(1)) = b
!
!   !also, make the shift in the working array as well, so that there 
!   !are no errors when determining the minimum:
!
!   a = work(k(1))
!   b  = work(1)
!   work(k(1)) = b
!   work(1) = a
!   
!   !--------------------
!   !minimum value:
!   !---------------------
!   
!   k = minloc(work)
!   write(4,*) 'k-min', k
!   write(4,'(A15,100000F7.2)') 'working array: ', (work(i), &
!      i = 1, dim_work)
!   write(4,*)
!   
!   a = numbers2(j + k(1)) !the i-th to last minimimum is the 
!                                (i+k)th element
!   b = numbers2(j + dim_work) !and it should be switched with the 
!                                i+dim_work-th to last element:
!   numbers2(j + dim_work) = a
!   numbers2(j + k(1)) = b
!
!   j = j + 1 !one value from the left-hand side of the series has been
!!              sorted
!   write(4,*)
!   write(4,*) 'pass-first-and-last', j
!   write(4,'(A9,100000F7.2)') 'numbers: ', (numbers2(i), i = 1, n2)
!   write(4,*) '------------------------------------------------------'
!   write(4,*)
!
!   deallocate(work) !because the dimension will change in the next step:
!   dim_work = dim_work - 2
!end do
!
!close(4)
!
!call CPU(start2)
!
!end subroutine first_and_last
!
end module sorting_subroutines
