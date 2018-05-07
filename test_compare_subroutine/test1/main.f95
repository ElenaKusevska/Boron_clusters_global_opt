program abcd
   use global_module
   use auxilary_subroutines
implicit none

real(kind=8), allocatable, dimension(:,:) :: matrix
real(kind=8), allocatable, dimension(:) :: masa
integer :: i, j, n

n = 20
allocate( matrix(n,3), masa(n) )

! Add points that make the shape of a ge
do i = 1, 10
   matrix(i,1) = i
   matrix(i,2) = 0
   matrix(i,3) = 0
end do

do i = 1, 5
   matrix(10+i,1) = i
   matrix(10+i,2) = 1
   matrix(10+i,3) = 0
end do

do i = 1, 5
   matrix(15+i,1) = 1
   matrix(15+i,2) = i
   matrix(15+i,3) = 0
end do

do i = 1, n
   masa(i) = 1.0
end do

open(unit=1, file='file1-1.xyz', action='write', status='replace')
write(1,'(A2)') '20'
write(1,'(A4)') 'abcd'
do i = 1, 20
   write(1,'(A19,3F11.7)') ' H                 ', (matrix(i,j), j = 1, 3)
end do
close(1)

call center_and_main_inertia_tensor (masa,n,matrix)

open(unit=2, file='file1-2.xyz', action='write', status='replace')

write(2,'(A2)') '20'
write(2,'(A4)') 'abcd'
do i = 1, 20
   write(2,'(A19,3F11.7)') ' H                 ', (matrix(i,j), j = 1, 3)
end do

close(2)

do i = 1, 10
   matrix(i,1) = i
   matrix(i,2) = 0
   matrix(i,3) = 0
end do

do i = 1, 5
   matrix(10+i,1) = i
   matrix(10+i,2) = -1
   matrix(10+i,3) = 0
end do

do i = 1, 5
   matrix(15+i,1) = 1
   matrix(15+i,2) = -i
   matrix(15+i,3) = 0
end do

open(unit=3, file='file2-1.xyz', action='write', status='replace')
write(3,'(A2)') '20'
write(3,'(A4)') 'abcd'
do i = 1, 20
   write(3,'(A19,3F11.7)') ' H                 ', (matrix(i,j), j = 1, 3)
end do
close(3)

call center_and_main_inertia_tensor (masa,n,matrix)

open(unit=4, file='file2-2.xyz', action='write', status='replace')

write(4,'(A2)') '20'
write(4,'(A4)') 'abcd'
do i = 1, 20
   write(4,'(A19,3F11.7)') ' H                 ', (matrix(i,j), j = 1, 3)
end do

close(4)

end program abcd
