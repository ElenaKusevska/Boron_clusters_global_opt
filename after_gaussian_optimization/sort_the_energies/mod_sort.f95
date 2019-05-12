module sorting_subroutines
implicit none
contains

subroutine CPU(start1)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)

real(dp), intent(in) :: start1
real(dp) :: finish1, time1

call cpu_time(finish1)
time1 = finish1 - start1
write(3,*) 'CPU TIME:', time1

end subroutine CPU

subroutine bubble_sort(jobs2, energies2, njobs2)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)

integer, allocatable, dimension(:), intent(inout) :: jobs2
real(dp), allocatable, dimension(:), intent(inout) :: energies2
integer, intent(in) :: njobs2
integer :: i, test, x, y, j
real(dp) :: start2, a, b

call cpu_time(start2)

test = 0
j = 1
do while (test .eq. 0)
   test = 1
   do i = njobs2 - 1, 1, -1
      if ( energies2(i) .lt. energies2(i+1) ) then
         a = energies2(i)
         b = energies2(i+1)
         energies2(i + 1) = a
         energies2(i) = b
         x = jobs2(i)
         y = jobs2(i + 1)
         jobs2(i + 1) = x
         jobs2(i) = y
         test = 0
      end if
   end do
   write(3,*)
   write(3,*) 'pass-bubble', j
   write(3,*)
   do i = 1, njobs2
      write(3, *) 'energy', energies2(i), 'job_order', jobs2(i)
   end do
   j = j + 1
end do     

open(unit=30, file='bubble.txt', status='replace', action='write')

write(30,*) 'pass-bubble', j
write(30,*)
do i = 1, njobs2
   write(30, *) 'energy', energies2(i), 'job_order', jobs2(i)
end do

close(30)

call CPU(start2)

end subroutine bubble_sort

subroutine first_and_last(jobs3, energies3, njobs3)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
integer, allocatable, dimension(:), intent(inout) :: jobs3
real(dp), allocatable, dimension(:), intent(inout) :: energies3
integer, intent(in) :: njobs3
real(dp), allocatable, dimension(:) :: work
integer, dimension(1) :: k
integer :: dim_work,  x, y, i, j
real(dp) :: start3, a, b

call cpu_time(start3)

i = 0
dim_work = njobs3

do while (dim_work .gt. 1)
   allocate ( work(dim_work) )
   do j = 1, dim_work
      work(j) = energies3(i + j)
   end do
   
   !---------------------
   !maxima:
   !--------------------
   
   k = maxloc(work)
   write(3,*) 'k-max', k
   !
   !energies:
   !
   a = energies3(i + k(1))
   b = energies3(i + 1)
   energies3(i + 1) = a
   energies3(i + k(1)) = b
   !
   a = work(k(1))
   b = work(1)
   work(1) = a
   work(k(1)) = b
   !
   !jobs:
   !
   x = jobs3(i + k(1))
   y = jobs3(i + 1)
   jobs3(i + 1) = x
   jobs3(i + k(1)) = y
   
   !--------------------
   !minima:
   !---------------------
   
   k = minloc(work)
   write(3,*) 'k-min', k
   !
   !energies:
   !
   a = energies3(i + k(1))
   b = energies3(i + dim_work)
   energies3(i + dim_work) = a
   energies3(i + k(1)) = b
   !
   !jobs:
   !
   x = jobs3(i + k(1))
   y = jobs3(i + dim_work)
   jobs3(i + dim_work) = x
   jobs3(i + k(1)) = y
   
   deallocate(work)
   write(3,*)
   write(3,*) 'pass-first-and-last', i
   write(3,*)
   do j = 1, njobs3
      write(3, *) 'energy', energies3(j), 'job_order', jobs3(j)
   end do
   dim_work = dim_work - 2
   write(*,*) 'dim_work', dim_work
   i = i + 1
   write(*,*) 'i', i
end do

open(unit=40, file='half-and-half.txt', status='replace', action='write')

write(40,*) 'pass-first-and-last', i
write(40,*)
do j = 1, njobs3
   write(40, *) 'energy', energies3(j), 'job_order', jobs3(j)
end do

close(40)

call CPU(start3)

end subroutine first_and_last

end module sorting_subroutines
