program main_program
   use iso_fortran_env
   Use combinations_module
implicit none
!
!---------------------------------------------------------------------
!Find the kth lexicographical combination
!in the set of combinations of m out of n
!---------------------------------------------------------------------
!
integer(kind=int64) :: k_comb, m, n, test, i, j, previous_step, m_comb, l
integer, allocatable, dimension(:) :: combination
logical :: finish
!
!--------------------------------------------------
!example:
!--------------------------------------------------
!
!number of combinations of 5 out of 10 = 252
!
!lets find the 100th in lexicographical order.
!
!----------------------------------------------
!The very first combination is:
!--------------------------------------------
!
![1 2 3 4 5]
!
!and, then the rest are combinations of the type:
!
!1 . . . . . . . . .
!
!But, how many of them are there?
!
!Well, there are nine elements left, and for the initial set, and four for
!the set of elements in the particular combination, so the number of
!1 x x x x combinations is:
!
!C 4 of 9 = 189
! 
!So, there are 189 more such combinations, and the 199th is of the type
!2 x x x x. More precisely, the 199th is 2 3 4 5 6.
!
!199th : . 2 . . . . . . . . 
!
!That is over 100, so we stick with the 1 . . . . type combinations, which are from 1 to 189, and thus contain the 100th
!combination.
!
!------------------------------------
!------------------------------------
!
!      [ 1 . . . . ]
!
!
!---------------------------------
!---------------------------------
!
!1st : 1 2 . . . . . . . . 
!
!---------------------------------
!---------------------------------
!
!c 3 of 8 = 56 + 1 = 57
!
!57th: 1 . 3 . . . . . . . 
!            
!--------------------------------------
!--------------------------------------
!
!c 3 of 7 = 35 + 57 = 92
!
!92nd: 1 . . 4 . . . . . .
!
!-------------------------------------
!-------------------------------------
!
!c 3 of 6 = 20
!
!112th : 1 . . . 5 . . . .
!
!---------------------------------------
!-----------------------------------------
!
!           [ 1 4 . . . ]
!
!----------------------------------------------
!--------------------------------------------
!
!92nd: 1 . . 4 5 . . . . .
!
!------------------------------------------------
!------------------------------------------------
!
!c 2 of 5 = 10 + 92 = 102
!
!102th 1 . . 4 . 6 . . . .
!
!------------------------------------------
!------------------------------------------
!
!           [1 4 5 . . ]
!
!92nd: 1 . . 4 5 6 . . . .
!
!-----------------------------------------
!------------------------------------------
!
!c 1 of 4 = 4
!
!96th 1 . . 4 5 . 7 . . .
!
!--------------------------------------
!----------------------------------------
!
!c 1 of 3 = 3
!
!99th: 1 . . 4 5 . . 8 . .
!
!------------------------------------------
!------------------------------------------
!
!c 1 of 2 = 2 + 99 = 101
!
!101st: 1 . . 4 5 . . . 9
!
!           [1 4 5 8 . ] 
!
!--------------------------------------------
!----------------------------------------------
!
!99th: 1 . . 4 5 . . 8 9 . 
!
!c 0 of 1  = 1
!
!--------------------------------------------------
!--------------------------------------------------
!
!99 + 1 = 100
!
!100th: 1 . . 4 5 . . 8 . 10
!
!-------------------------------------------------------------------
!---------------------------------------------------------------------
!
write(*,*) 'm out of n?'
read(*,*) m
read(*,*) n
!
allocate( combination(m) )
!
write(*,*) 'k_comb?'
read(*,*) k_comb
!
!------------------------
!First combination:
!------------------------
!
do i = 1, m
   combination(i) = i
end do
!
test = 1
i = 1 ! number of elements in the combination (m-x)
j = n !number of elements in the set from which the combination is
!           selected (n-p)
finish = .false.
m_comb = 1
!
!---------------------------------
!The other combinations:
!---------------------------------
!
do while ( (i .le. m) .and. (finish .eqv. .false.) )
   write(*,*) ' TEST !!!!! '
   do while (j .ge. 0)
      m_comb = m_comb + 1
      write(*,*) 'm_comb', m_comb
      combination(i) = m_comb
      previous_step = test
      test = test + combination_x_of_p(m-i,j-i)
      write(*,*) 'test', test
      if (test == k_comb) then
         do l = i + 1, m
            combination(l) = combination(l-1) + 1
         end do
         finish = .true.
         exit
      else if (test .gt. k_comb) then 
         m_comb = m_comb - 1
         combination(i) = m_comb
         test = previous_step
         m_comb = m_comb + 1
         exit
      else if (test .lt. k_comb) then
         j = j - 1
      end if
   end do 
   i = i + 1
end do
!
write(*,'(I5,A14,I5,A6,I5,A1)') k_comb, 'combination of', m, 'out of', n, ':'
write(*,*) combination
!
end program main_program
