module combinations_module
   use iso_fortran_env
implicit none
contains
!
!----------------------------------
!combinations - function:
!----------------------------------
!
integer function combination_x_of_p(x,p)
implicit none
!
integer(kind=int64), intent(in) :: x, p
integer(kind=int64) :: diff, prod, x_fact, p_x_fact, i
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
end function combination_x_of_p


end module combinations_module
