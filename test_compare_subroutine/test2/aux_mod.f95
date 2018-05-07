module auxilary_subroutines
implicit none
contains
!
!--------------------------------------------------------
! Combinations:
!-------------------------------------------------------
!
integer(kind=16) function combination_x_of_p(x,p)
implicit none
!
integer(kind=16), intent(in) :: x, p
integer(kind=16) :: diff, prod, x_fact, p_x_fact, i
!
diff = p - x
write(*,*) 'x', x, 'p', p
write(*,*) 'diff', diff
!
!--------------------
! p-x > x :
!--------------------
!
if (diff .ge. x) then
   write(*,*) 'p - x > x'
   prod = 1
   do i = 1, x
      prod = prod * (diff + i)
      write(*,*) 'prod', prod
   end do
   !
   if (x .ne. 0) then
      x_fact = 1
      do i = 1, x
         x_fact = x_fact * i
      end do
   end if
   if (x .eq. 0) then
      x_fact = 1
   end if
   !
   write(*,*) 'prod', prod
   write(*,*) 'x_fact', x_fact
   combination_x_of_p = prod/x_fact
   write(*,*) 'combination x of p', combination_x_of_p
!
!--------------------------
! p-x < x :
!--------------------------
!
else if (diff .lt. x) then
   write(*,*) 'x > p - x'
   prod = 1
   do i = 1, diff
      prod = prod * (x+i)
   end do
   !
   p_x_fact = 1
   do i = 1, p - x
      p_x_fact = p_x_fact * i
   end do
   !
   combination_x_of_p = prod/p_x_fact
end if
!
!-----------------------------------------------------------------
! example -> x .lt. p-x
!-----------------------------------------------------------------
!
!   p = 8
!   x = 2
!   p - x = 6
!
!   x < p-x
!
!   6! 7 8            p-x+1 p-x+2
! -------------- = --------------------
!   6! 2!               x!                      
!
!-----------------------------------------------------------------
! example -> x .gt. p=x
!-----------------------------------------------------------------
!
!   p = 8
!   x = 5
!   p - x = 3
!
!   5! 6 7 8        x+1 x+2 x+3
! ------------ = ------------------
!   5! 3!           p - x!
!
! 
end function combination_x_of_p
!
!-----------------------------------------------------------------
!Subroutine to generate the distance matrix for a structure:
!-----------------------------------------------------------------
!
subroutine distance_matrix (a_matrix, a_atoms, dist_matrix)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
real(dp), allocatable, dimension(:,:), intent(in) :: a_matrix
integer, intent(in) :: a_atoms
real(dp), allocatable, dimension(:,:), intent(out) :: dist_matrix
integer :: i, j
!
allocate(dist_matrix(a_atoms,a_atoms))
do i = 1, a_atoms
  do j = 1, a_atoms
     dist_matrix(i,j) = sqrt( (a_matrix(i,1)-a_matrix(j,1))**(2.) + &
       (a_matrix(i,2)-a_matrix(j,2))**(2.) + &
       (a_matrix(i,3)-a_matrix(j,3))**(2.) )
  end do
end do
!
write(*,*)'Distance matrix:'
write(*,*)'----------------------------------------'
!
do i = 1, a_atoms
  write(*,'(20F7.3)') (dist_matrix(i,j), j=1,a_atoms)
end do
!
end subroutine distance_matrix
!           
!----------------------------------------------------------------
!Periodic system - read the mass for every atom in the molecule:
!----------------------------------------------------------------
!
subroutine periodic (labels_in, a, masses_out)
Implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
character(len=1), allocatable, dimension(:), intent(in) :: labels_in
integer :: a, i, k !number of atoms
real(dp), dimension(2) :: masses
character(len=1), dimension(2) :: atoms
real(dp), allocatable, dimension(:), intent(out) :: masses_out
!
atoms = (/'B', 'V'/)     
masses = (/10.81, 50.9415/)
!
allocate ( masses_out(a) )
!
!-----------------------
!Assigning the masses:
!-----------------------
!
do i = 1, a
    do k = 1, 2
        if (labels_in(i) == atoms(k)) then
            masses_out(i) = masses(k)
        end if
    end do 
end do 
!
end subroutine periodic
!
!---------------------------------------------------------
!Calculate CPU time:
!---------------------------------------------------------
!
subroutine CPU (start, method15)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), intent(in) :: start
character(len=17), intent(in) :: method15
real(dp) :: finish, cputime
!
call CPU_time(finish)
cputime = finish - start
write(25,*)
write(25,'(A18)') method15, ':'
write(25,*) 'CPU time:', cputime
!
end subroutine CPU
!
end module auxilary_subroutines
