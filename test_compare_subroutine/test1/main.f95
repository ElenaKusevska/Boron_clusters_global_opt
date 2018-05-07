
!------------------------------------------------------------------
! Program to test if my concept for checking that two structures 
! are equal works properly.
!------------------------------------------------------------------

program abcd
use global_module
use auxilary_subroutines
implicit none

real(kind=8), allocatable, dimension(:,:) :: matrixa, matrixb
real(kind=8), allocatable, dimension(:) :: mass
integer :: i, j, n, test_variable

n = 20
allocate( matrixa(n,3), matrixb(n,3), mass(n) )

! Define the coordinates of points on the grid. You can see the shape
! when the files are printed. It's kinda like an L:
do i = 1, 10
   matrixa(i,1) = i
   matrixa(i,2) = 0
   matrixa(i,3) = 0
end do

do i = 1, 5
   matrixa(10+i,1) = i
   matrixa(10+i,2) = 1
   matrixa(10+i,3) = 0
end do

do i = 1, 5
   matrixa(15+i,1) = 1
   matrixa(15+i,2) = i
   matrixa(15+i,3) = 0
end do

! Make these points be hydrogen atoms (needed for some of the
! subrutines that are called). This array will be used for matrixa and
! recycled/reused for matrixb:
do i = 1, n
   mass(i) = 1.0
end do

! Print the .xyz file of these "hydrogen" points on the grid:
open(unit=1, file='file1-1.xyz', action='write', status='replace')
write(1,'(A2)') '20'
write(1,'(A4)') 'abcd'
do i = 1, 20
   write(1,'(A19,3F11.7)') ' H                 ', (matrixa(i,j), j = 1, 3)
end do
close(1)

! Move then so that the center of mass is at (0,0,0):
call center_and_main_inertia_vector (mass,n,matrixa)

! And print the new coordinates:
open(unit=2, file='file1-2.xyz', action='write', status='replace')
write(2,'(A2)') '20'
write(2,'(A4)') 'abcd'
do i = 1, 20
   write(2,'(A19,3F11.7)') ' H                 ', (matrixa(i,j), j = 1, 3)
end do

close(2)

! Now, make another set of hydrogen atoms of the same shape, but
! rotated about itself:
do i = 1, 10
   matrixb(i,1) = i
   matrixb(i,2) = 0
   matrixb(i,3) = 0
end do

do i = 1, 5
   matrixb(10+i,1) = i
   matrixb(10+i,2) = -1
   matrixb(10+i,3) = 0
end do

do i = 1, 5
   matrixb(15+i,1) = 1
   matrixb(15+i,2) = -i
   matrixb(15+i,3) = 0
end do

! And write them to a file:
open(unit=3, file='file2-1.xyz', action='write', status='replace')
write(3,'(A2)') '20'
write(3,'(A4)') 'abcd'
do i = 1, 20
   write(3,'(A19,3F11.7)') ' H                 ', (matrixb(i,j), j = 1, 3)
end do
close(3)

! And move them so that the center of mass is at (0,0,0):
call center_and_main_inertia_vector (mass,n,matrixb)

! And write the new, translated coordinates to a file:
open(unit=4, file='file2-2.xyz', action='write', status='replace')
write(4,'(A2)') '20'
write(4,'(A4)') 'abcd'
do i = 1, 20
   write(4,'(A19,3F11.7)') ' H                 ', (matrixb(i,j), j = 1, 3)
end do

close(4)

! and now, compare matrixa(file1-2.xyz) and matrixb(file2-2.xyz):
call compare(matrixa,matrixb,n,test_variable)

end program abcd
