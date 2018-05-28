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
integer(kind=16) :: diff, aprod, a_fact, i, j, k, a, b, c
integer(kind=16), allocatable, dimension(:) :: prod, x_fact, p_x_fact
integer(kind=16), dimension(13) :: divisors
!
!primes = (/2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, &
!   59, 61, 67, 71, 73, 79, 83, 89, 97, 101/)
divisors = (/2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14/)
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
   !
   !determine the values of the arrays prod and x_fact:
   !
   allocate( prod(x), x_fact(x) )
   !
   do i = 1, x
      prod(i) = diff + i
      !
      x_fact(i) = i
   end do
   !
   !reduce to the smallest common divisor:
   !
   do i = x, 1, -1
      a = prod(i)
      do j = x, 1, -1
         b = x_fact(j)
         do k = 13, 1, -1
            c = divisors(k)
            if (b .ge. c) then
               if (a .ge. c) then
                  if (mod(a,c) == 0) then
                     if (mod(b,c) == 0) then
                        prod(i) = prod(i)/divisors(k)
                        a = prod(i)
                        x_fact(j) = x_fact(j)/divisors(k)
                        b = x_fact(j)
                        write(*,*) 'a, b, c', a, b, c
                        write(*,*) 'mod-a-c', mod(a,c)
                        write(*,*) 'mod-b-c', mod(b,c)
                        write(*,*) 'prod', prod
                        write(*,*) 'x_fact', x_fact
                     end if
                  end if
               end if
            end if
         end do
      end do
   end do
   !
   !determine the number of combinations:
   !
   aprod = 1
   do i = 1, x
      aprod = aprod * prod(i)
      write(*,*) 'prod', prod
   end do
   !
   if (x .ne. 0) then
      a_fact = 1
      do i = 1, x
         a_fact = a_fact * x_fact(i)
      end do
   end if
   if (x .eq. 0) then
      x_fact = 1
   end if
   !
   write(*,*) 'prod', aprod
   write(*,*) 'x_fact', a_fact
   combination_x_of_p = aprod/a_fact
   write(*,*) 'combination x of p', combination_x_of_p
!
!--------------------------
! p-x < x :
!--------------------------
!
else if (diff .lt. x) then
   write(*,*) 'x > p - x'
   !
   !determine the values of the arrays p_x_fact and prod
   !
   allocate( prod(diff), p_x_fact(diff) )
   !
   do i = 1, diff
      prod(i) = x + i
      !
      p_x_fact(i) = i
   end do
   !
   !reduce prod and p_x_fact to the smallest common divisor:
   !
   do i = diff, 1, -1
      a = prod(i)
      do j = diff, 1, -1
         b = p_x_fact(j)
         do k = 13, 1, -1
            c = divisors(k)
            if (b .ge. c) then
               if (a .ge. c) then
                  if (mod(a,c) == 0) then
                     if (mod(b,c) == 0) then
                        prod(i) = prod(i)/divisors(k)
                        a = prod(i)
                        p_x_fact(j) = p_x_fact(j)/divisors(k)
                        b = p_x_fact(j)
                     end if
                  end if
               end if
            end if
         end do
      end do
   end do
   !
   !determine the number of combinations:
   !
   aprod = 1
   do i = 1, diff
      aprod = aprod * prod(i)
   end do
   !
   a_fact = 1
   do i = 1, diff
      a_fact = a_fact * p_x_fact(i)
   end do
   !
   combination_x_of_p = aprod/a_fact
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
