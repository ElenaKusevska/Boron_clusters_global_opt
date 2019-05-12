program sort_energies
   use sorting_subroutines
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)

integer, allocatable, dimension(:) :: jobs, jobs1
real(dp), allocatable, dimension(:) :: energies, energies1
integer :: i, njobs, io, var1
real(dp) :: var2

!-------------------------------------------------
!count the number of jobs - energies
!-------------------------------------------------

open(unit=1, file='energies.txt', status='old', action='read')
njobs = 0
do
   read(1,*,iostat=io) var1, var2
   if (io /= 0) exit
   njobs = njobs + 1
end do
write(*,*) 'njobs', njobs 
close(1)

allocate( jobs(njobs), jobs1(njobs), energies(njobs), energies1(njobs) )

!-----------------------------------------------------------------
!read the job number and the energy
!-----------------------------------------------------------------

open(unit=2, file='energies.txt', status='old', action='read')

do i = 1, njobs
   read(2,*) jobs(i), energies(i)
end do

!--------------------------------------------------
!perform the sorting:
!-------------------------------------------------

open(unit=3, file='output.txt', status='replace', action='write')

energies1 = energies
jobs1 = jobs

call bubble_sort(jobs, energies, njobs)

call first_and_last(jobs1, energies1, njobs)

end program sort_energies
