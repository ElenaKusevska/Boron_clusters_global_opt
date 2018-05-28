program algorithm2
   use generate_structures
   use auxilary_subroutines
   use global_module
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
real(dp) :: beginning, d_molecule
integer :: vanadium, i, m0, c1, m1, n1, dim_grid
integer(kind=16) :: d_grid, n16
real(dp), allocatable, dimension(:,:) :: grid, coords1, first_combination
logical :: vanadium2_yes_or_no 
character(len=1) :: dim_2_3, many_combinations
character(len=17) :: method
!
call CPU_time(beginning)
call random_seed()
!
open(unit=25, file='output.txt', status='replace', action='write')
!
c1 = 150 !file names for priniting the coordinates
m1 = 300 !file names for priniting the molecule input files
!
write(*,*) 'how many vanadium atoms? [1/2]'
read(*,*) vanadium
write(*,*) 'how many atoms does the molecule contain?'
read(*,*) n1
!
if (vanadium == 1) then
  vanadium2_yes_or_no = .false.
else if (vanadium == 2) then
  vanadium2_yes_or_no = .true.
end if
!
!----------------------------------------------------
!1) generate coordinates using a grid
!----------------------------------------------------
method = 'combine_from_grid'
!
call generate_grid(dim_grid, grid)
d_grid = dim_grid
n16 = n1
call number_of_combinations(d_grid, n16, many_combinations)
if (many_combinations == 'N') then
   call grid_to_molecule_small (dim_grid, n1, method, &
      vanadium2_yes_or_no, grid, c1, m1, first_combination)
else if (many_combinations == 'Y') then
   call grid_to_molecule_large (dim_grid, n1, method, &
      vanadium2_yes_or_no, grid, c1, m1, first_combination)
 end if
!
call CPU(beginning, method)
!
!----------------------------------------------------
!2) generate coordinates using a random structure
!----------------------------------------------------
!
method = 'random_coordinate'
!
write(*,*) 'average size of the molecule'
read(*,*) d_molecule
!
i = 1
do while (i .le. 3)
   m0 = m1
   call generate_random3d(n1, d_molecule, coords1)
   call write_to_file (coords1, method, n1, vanadium2_yes_or_no, c1, m1)
   if (m1 .gt. m0) then
      i = i + 1
   end if
end do
!
i = 1
do while (i .le. 2)
   m0 = m1
   call generate_random2d(n1, d_molecule, coords1)
   call write_to_file (coords1, method, n1, vanadium2_yes_or_no, c1, m1)
   if (m1 .gt. m0) then
      i = i + 1
   end if
end do
!
call CPU (beginning, method)
!
!------------------------------------------------------------------
!3) Coordinate walking
!-----------------------------------------------------------------
!
method = 'coordinate_walkin'
!
!dim_2_3 = 'N'
!call coordinate_walk(first_combination, method, dim_2_3, vanadium2_yes_or_no, n1, c1, m1)
dim_2_3 = 'Y'
call coordinate_walk(first_combination, method, dim_2_3, vanadium2_yes_or_no, n1, c1, m1)
!
call CPU (beginning, method)
!
!---------------------------------------------------------------
!3) generate coordinates using a library of previous structures:
!---------------------------------------------------------------
!
!
end program algorithm2
