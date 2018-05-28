module global_module
   use auxilary_subroutines
implicit none
contains
!
!------------------------------------------------------------------
!subroutine in order to test if the molecular coordinates fulfill 
!certain criteria. 
!------------------------------------------------------------------
!
subroutine test1 (coords6, n6, logic6)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
real(dp), allocatable, dimension(:,:), intent(in) :: coords6
integer, intent(in) :: n6
real(dp), allocatable, dimension(:,:) :: d_matrix
real(dp), allocatable, dimension(:) :: distance_vector
integer :: i, j
character(len=1), intent(out) :: logic6
integer, allocatable, dimension(:) :: dv_test
integer :: dv_test_sum
!
!-------------------------------------------------------------------
!Build distance matrix in order to check if no two distances are less
!than 0.5 angstrom, and also, to make sure that there is no atom for
!which all of the distances are alrger than 2.5
!-------------------------------------------------------------
!
call distance_matrix (coords6, n6, d_matrix)
!
logic6 = 'Y'
!
!check if any two distances are less than 1.0 angstrom
!
do i = 1, n6-1
   do j = i+1, n6
      if (d_matrix(i,j) .lt. 0.8) then
         logic6 = 'N'
         exit
      end if 
   end do
   if (logic6 == 'N') then 
      deallocate(d_matrix)
      exit
   end if
end do
!
!check if at least n - 2 of the distances for every atom 
!are less than 2 angstrom this makes sure that the structure is compact 
!enough.
!
if (logic6 == 'Y') then
   allocate( dv_test(n6), distance_vector(n6) )
   do i = 1, n6
      do j = 1, n6
         distance_vector(j) = d_matrix(i,j) !take one row out of the 
!                                            distance matrix.
      end do
      !
      do j = 1, n6
         if (distance_vector(j) .ge. 2.0) then !record distances greater
!                                               than 2 in that row
            dv_test(j) = 1
         else if (distance_vector(j) .lt. 2.0)then
            dv_test(j) = 0
         end if
      end do
      !
      dv_test_sum = 0 !count number of distances greater than 2 in that
!                                                           row
      do j = 1, n6
         dv_test_sum = dv_test_sum + dv_test(j)
      end do
      !
      if (dv_test_sum .ge. n6 - 1 ) then
         logic6 = 'N'
         exit
         deallocate (d_matrix)
      end if
   end do
   deallocate (distance_vector)
   deallocate (dv_test)
end if      
!
if ( logic6 == 'Y' ) then
   deallocate(d_matrix)
end if
!
end subroutine test1
!
!------------------------------------------------------------------
!subroutine to place the center of mass of the molecule in the origin
!of the coordinate system, find the moment of inertia along which the
!maximum order of rotation is, and align the molecule along this axis
!-------------------------------------------------------------------
!
subroutine center_and_main_inertia_vector (mass,n3,coords)
implicit none 
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
real(dp), allocatable, dimension(:), intent(in) :: mass
integer, intent(in) :: n3
real(dp), allocatable, dimension(:,:), intent(inout) :: coords
real(dp), allocatable, dimension(:,:) :: coordsa
real(dp), dimension(3) :: center !center of mass
real(dp), dimension(3,3) :: inert !inertia tensor
real(dp) :: sumxa, sumxb, sumya, sumyb, sumza, sumzb
integer :: i, j
!
!For DSYEV:
!
integer :: n, lwork, lda, info
character :: jobz, uplo
real(dp) :: lwork1
real(dp), dimension(3) :: w, eigenv1, eigenv2, eigenv3, avector 
!                                               eigenvalues and
!                                                   eigenvectors
real(dp), allocatable, dimension(:,:) :: work
!
!For the maximum order of rotation:
!
real(dp), allocatable, dimension(:,:) :: mola, molb, molc
integer :: r, order, test_rot, z_axis
real(dp), dimension(3) :: z_axis_array
integer, dimension(1) :: max_axis 
integer, dimension(3) :: rot_order
real(dp) :: coss, sinn
real(dp), dimension(3,3) :: Rzz
!
!-------------------------------
!Determine center of mass:
!-------------------------------
!
!x-coordinate:
!
sumxa = 0.0
do i = 1, n3
    sumxa = sumxa + mass(i)*coords(i,1)
end do
!
sumxb = 0.0
do i =1, n3
    sumxb = sumxb + mass(i)
end do
!
center(1) = sumxa/sumxb
!
!y - coordinate:
!
sumya = 0.0
do i = 1, n3
    sumya = sumya + mass(i)*coords(i,2)
end do
!
sumyb = 0.0
do i = 1, n3
    sumyb = sumyb + mass(i)
end do
!
center (2) = sumya/sumyb
!
!z - coordinate:
!
sumza = 0.0
do i = 1, n3
    sumza = sumza + mass(i)*coords(i,3)
end do
!
sumzb = 0.0
do i = 1, n3
    sumzb = sumzb + mass(i)
end do 
!
center(3) = sumza/sumzb
!
!---------------------------------------
!moving the center of mass at (0,0,0)
!---------------------------------------
!
do i = 1, n3
    coords(i,1) = coords(i,1) - center(1)
end do
!
do i = 1, n3
    coords(i,2) = coords(i,2) - center(2)
end do
!
do i = 1, n3
  coords(i,3) = coords(i,3) - center(3)
end do
!
!-----------------------------------
!Calculating the inertia tensor: 
!-----------------------------------
!
inert(1,1) = 0.0
do i = 1, n3
    inert(1,1) = inert(1,1) + mass(i)*( (coords(i,2))**(2.) + &
      (coords(i,3))**(2.) )
end do
!
inert(1,2) = 0.0
do i = 1, n3
    inert(1,2) = inert(1,2) - mass(i)*coords(i,2)*coords(i,1)
end do
!
inert(2,1) = inert(1,2)
!
inert(1,3) = 0.0
do i = 1, n3
    inert(1,3) = inert(1,3) - mass(i)*coords(i,1)*coords(i,3)
end do
!
inert(3,1) = inert(1,3)
!
inert(2,2) = 0.0
do i = 1,n3
    inert(2,2) = inert(2,2) + mass(i)*( (coords(i,1))**(2.) + &
      (coords(i,3))**(2.) )
end do
!
inert(2,3) = 0.0
do i = 1, n3
    inert(2,3) = inert(2,3) - mass(i)*coords(i,2)*coords(i,3)
end do
!
inert(3,2) = inert(2,3)
!
inert(3,3) = 0.0
do i = 1, n3
    inert(3,3) = inert(3,3) + mass(i)*( (coords(i,1))**(2.) + &
      (coords(i,2))**(2.) )
end do 
!
!-----------------------------------------------------------------
!Diagonalize the inertia tensor - determine the three principal
!components of the inertia tensor. 
!-----------------------------------------------------------------
!
n = 3
lda = 3
jobz = 'N'
uplo = 'U'
lwork = -1
!
!--------------------------------------------
!determine the leading dimention of work:
!--------------------------------------------
!
allocate (work(1,3))
!
   call DSYEV(jobz, uplo, n, inert, lda, w, work, lwork, info) 
   write(*,*) 'info:', info, 'lwork1:', work(1,1)
!
   lwork1 = work(1,1)
   deallocate (work)
!
   lwork = int(lwork1)
!
allocate (work(1,lwork))
!
!---------------------------------
!diagonalize:
!---------------------------------
!
jobz = 'V'
!
    call DSYEV(jobz, uplo, n, inert, lda, w, work, lwork, info)
!
    write(*,*) 'eigenvalues of the inertia tensor:', w 
    write(*,*)  'info', info
!
deallocate (work)
!
!on input, the DSYEV subroutine takes the matrix inert, which represents
!the inertia tensor. On output it gives the same matrix (inert), but now
!it contains the eigenvectors of the inertia tensor in its columns. 
!
do i = 1, 3
    eigenv1(i) = inert(i,1)
    eigenv2(i) = inert(i,2)
    eigenv3(i) = inert(i,3)
end do
!
!----------------------------
!create molecule matrix:
!----------------------------
!
!the n - rows of this matrix represent the n - atoms
!
allocate (mola(n3,4))
!
do i = 1, n3
  mola(i,1) = mass(i) !column 1 contains the masses
  mola(i,2) = coords(i,1) !column 2 contains the x- coordinates
  mola(i,3) = coords(i,2) !column 3 contains the y - coordinates
  mola(i,4) = coords(i,3) !column 4 contains the z - coordinates
end do
!
!----------------------------------------------------------------
!Find the maximum order of rotation around every inertia vector
!----------------------------------------------------------------
!
allocate (molb(n3,4), molc(n3,4))
!
!---------------------------------------------------------------
!Find the maxium order of rotation around every inertia vector
!--------------------------------------------------------------
!
!mola - original molecule
!molb - the aligned molecule
!molc - the rotated molecule
!
write(*,*)
write(*,*) 'eigenv1'
write(*,*) eigenv1
write(*,*)
!
call align(eigenv1, eigenv2, mola, n3, molb) !align molecule with
                                             !inertia vector eigenv1
write(*,*) 'after align'
write(*,*) eigenv1
!
do i = 6, 1, -1
    r = i
    rot_order(1) = i
    call rotation(molb, n3, r, molc)
    call compare(molb, molc, n3, test_rot)
    if (test_rot == 1) exit !the maximum rotation order is recorded in
                        !the array rot_order
end do
write(*,*) 'rot_order:', rot_order(1)
!
write(*,*)
write(*,*) 'eigenv2'
write(*,*) eigenv2
write(*,*)
call align(eigenv2, eigenv1, mola, n3, molb)
write(*,*) 'after align'
write(*,*) eigenv2
!
do i = 6, 1, -1
    r = i
    rot_order(2) = i
    call rotation(molb, n3, r, molc)
    call compare(molb, molc, n3, test_rot)
    if (test_rot == 1) exit
end do
write(*,*) 'rot_order:', rot_order(2)
!
write(*,*)
write(*,*) 'eigenv3'
write(*,*) eigenv3
write(*,*)
call align(eigenv3, eigenv1, mola, n3, molb)
write(*,*) 'after align'
write(*,*) eigenv3
!
do i = 6, 1, -1
    r = i
    rot_order(3) = i
    call rotation(molb, n3, r, molc)
    call compare(molb, molc, n3, test_rot)
    if (test_rot == 1) exit
end do
write(*,*) 'rot_order:', rot_order(3)
!
order = maxval(rot_order) !order of maximum rotation axis
max_axis = maxloc(rot_order) !column of inertia matrix where the inertia
!                           vector along which is the maximum order of
!                           rotation, is situated. 
z_axis  = max_axis(1)
!
do i = 1, 3
  z_axis_array(i) =  inert(i,z_axis) !z - axis
end do
!
write(*,*) 'order', order, '# inertia tensor', max_axis
write(*,*) 'inertia tensor along maximum symmetry axis', z_axis_array
!
!----------------------------------------
!getting the final, aligned coordinates:
!----------------------------------------
!
do i = 1, 3
   if (i == z_axis) cycle
   do j = 1, 3
      avector = inert(j,i)
   end do
end do
call align(z_axis_array, avector, mola, n3, molb)
do i = 1, n3
   coords(i,1) = molb(i,2)
   coords(i,2) = molb(i,3)
   coords(i,3) = molb(i,4)
end do
!
do i = 1, 3
   if (avector(i) .lt. 0.0) then
      avector(i) = -avector(i)
      do j = 1, 3
         coords(j,i) = -coords(j,i)
      end do
   end if
end do
coss = avector(1) / ( sqrt( avector(1) + avector(2) ) )
sinn = avector(2) / ( sqrt( avector(1) + avector(2) ) )
call z_rotation_matrix(coss, sinn, Rzz)
allocate(coordsa(n3,3))
do i = 1, 3
   do j = 1, 3
      coordsa(i,j) = coords(i,j)
   end do
end do
!
coords = matmul(coordsa,Rzz)
!
deallocate(mola, molb, molc, coordsa)
!
end subroutine center_and_main_inertia_vector
!
!----------------------------------------------------------------------
!This subroutine rotates a molecule around z - axis by an angle 
!defined in the call of the subroutine
!-------------------------------------------------------------------
!
subroutine rotation(moly, n4, r, molz)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
real(dp), allocatable, dimension(:,:), intent(in) :: moly
real(dp), allocatable, dimension(:,:) :: coords3, rotated
real(dp), dimension(3,3) :: Rz ! rotation axis
integer, intent(in) :: n4 !dimension of molecule
integer, intent(in) :: r !order of rotation axis
real(dp), allocatable, dimension(:,:), intent(out) :: molz
real(dp), parameter :: pi = 3.141592653589793
real(dp) :: cosz, sinz, k ! = real(r)
integer :: i, j !counters
!
allocate ( coords3(n4,3), rotated(n4,3), molz(n4,4) )
!
do i = 1, n4
    do j = 1, 3
        coords3(i,j) = moly(i,j+1)
    end do
end do
!
k = real(r)
cosz = cos(((2.)*pi)/k)
sinz = sin(((2.)*pi)/k)
!
call z_rotation_matrix (cosz, sinz, Rz)
!
rotated = matmul(coords3, Rz) 
!
do i = 1, n4
    molz(i,1) = moly(i,1) ! just transferring the atomic masses
end do
! 
do i = 1, n4
    do j = 1, 3
        molz(i,j+1) = rotated(i,j)
    end do
end do
!
deallocate (coords3, rotated)
!
end subroutine rotation
!
!-----------------------------------------------------------
!subroutines to generate the rotation matrices:
!---------------------------------------------------------
!
subroutine x_rotation_matrix (cosx, sinx, Rx)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), intent(in) :: cosx, sinx
real(dp), dimension(3,3), intent(out) :: Rx
!
Rx(1,1)=1
Rx(1,2)=0
Rx(1,3)=0
Rx(2,1)=0
Rx(2,2)=cosx
Rx(2,3)=-sinx
Rx(3,1)=0
Rx(3,2)=sinx
Rx(3,3)=cosx
!
end subroutine x_rotation_matrix
!
subroutine y_rotation_matrix (cosy, siny, Ry)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), intent(in) :: cosy, siny
real(dp), dimension(3,3), intent(out) :: Ry
!
Ry(1,1)=cosy
Ry(1,2)=0
Ry(1,3)=-siny
Ry(2,1)=0
Ry(2,2)=1
Ry(2,3)=0
Ry(3,1)=siny
Ry(3,2)=0
Ry(3,3)=cosy
!
end subroutine y_rotation_matrix
!
subroutine z_rotation_matrix(cosz, sinz, Rz)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), intent(in) :: sinz, cosz
real(dp), dimension(3,3), intent(out) :: Rz
!
Rz(1,1)=cosz
Rz(1,2)=-sinz
Rz(1,3)=0
Rz(2,1)=sinz
Rz(2,2)=cosz
Rz(2,3)=0
Rz(3,1)=0
Rz(3,2)=0
Rz(3,3)=1
!
end subroutine z_rotation_matrix
!
!----------------------------------------------------------
!Align the molecule so that the selcted moment of inertia 
!is along the z - axis of the cartesian coordinate system:
!-----------------------------------------------------------
!
subroutine align (vector, vectorb, mol1, n5, mol2)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)

real(dp), dimension(3), intent(inout) :: vector, vectorb
real(dp), dimension(3,1) :: vector1, vectora
real(dp), dimension(3,1) :: vecxz, vec1xz !rotation around x-axis to 
                !move the vector
                    !to the xz plane
real(dp), dimension(3,1) :: vecz, vec1z ! the vector is aligned
                    !with z - axis after rotation around y - axis
integer, intent(in) :: n5 !number of atoms in the molecule
real(dp), allocatable, dimension(:,:), intent(in) :: mol1 !input molecule
real(dp), allocatable, dimension(:,:) :: matrix5, coords5, rotated1, rotated2
                    !working functions
real(dp) :: cosx, sinx, cosy, siny ! to rotate space about the x & y - axes
real(dp), dimension(3,3) :: Rx, Ry !Rotation about x & y - axes
real(dp), allocatable, dimension(:,:), intent(out) :: mol2
integer :: i, j
!
!-----------------------------------------------------------------
!This subroutine aligns the axis of interest with the z - axis.
!It does this by determining the sine and cosine of the angles between
!that axis and the z- and x- or y- axes, 
!based on its position relative to the coordinate system, 
!building the rotation matrices for rotation around the y - and x - 
!matrices and then applying them to every atom in the molecule. 
!
!This is equivalent to either rotating the entire
!molecule, or "moving space", i.e. rotating the coordinate system,
!to allow the axis of interest to be aligned with the z - axis
!--------------------------------------------------------------
!
do i = 1, 3 !convert to column vector
   vector1(i,1) = vector(i)
   vectora(i,1) = vectorb(i)
end do
!
allocate (coords5(3,n5), mol2(n5,4), matrix5(n5,3))
allocate (rotated1(3,n5), rotated2(3,n5))
!
do i = 1, n5
    do j = 1, 3
        matrix5(i,j) = mol1(i,j+1)
    end do
end do
!
!convert to column vectors:
!
do i = 1, n5
   do j = 1, 3
      coords5(j,i) = matrix5(i,j)
   end do
end do
!
if (dabs(vector(1)) .le. 0.000000000001) then
   if (dabs(vector(2)) .gt. 0.000000000001) then
      !
      !----------------------------
      ! [0,y,z]
      !--------------------------
      !
      !----------------------------------------------
      !find matrix for rotation about the x-axis
      !-----------------------------------------------
      !
      cosx=vector(3)/sqrt( vector(2)**(2.) + vector(3)**(2.) )
      sinx=vector(2)/sqrt( vector(2)**(2.) + vector(3)**(2.) )
      !
      !
      call x_rotation_matrix (cosx, sinx, Rx)
      !
      !--------------------------------------------------------
      !rotate space about the x - axis to bring the inertia vector in
      !the xz-plane:
      !----------------------------------------------------------
      !
      rotated2 = matmul(Rx,coords5)  
      vecz = matmul(Rx,vector1)
      vec1z = matmul(Rx,vectora)
      write(*,*) 'vecz', vecz
      !
   else if (dabs(vector(2)) .le. 0.000000000001) then
      !
      !-------------------------------------------
      ! [0,0,z]
      !-----------------------------------------------
      !
      vecz = vector1
      vec1z = vectora
      write(*,*) 'vecz', vecz
      do i = 1, 3
         do j = 1, n5
            rotated2(i,j) = coords5(i,j)
         end do
      end do
   end if
   !
   !
else if (dabs(vector(1)) .gt. 0.000000000001) then
   if (dabs(vector(2)) .le. 0.000000000001) then
      !
      !---------------------------------------------
      ! [x, 0, z]
      !----------------------------------------------------
      !
      !---------------------------------------------------
      !find matrix for rotation about the y-axis
      !--------------------------------------------------
      !
      cosy = vector(3)/sqrt(vector(1)**(2.) + vector(3)**(2.))
      siny = vector(1)/sqrt(vector(1)**(2.) + vector(3)**(2.))
      !
      call y_rotation_matrix (cosy, siny, Ry)
      !
      !-------------------------------------------------------------
      !rotate space around the y - axis, so that the rotation axis lies
      !along the positive z - axis
      !----------------------------------------------------------------
      !
      rotated2 = matmul(Ry,coords5)
      vecz = matmul(Ry,vector1) !should have the form [0,0,z]
      vec1z = matmul(Ry,vectora)
      write(*,*) 'vecz', vecz
   else if (dabs(vector(2)) .gt. 0.000000000001) then 
      !
      !-------------------------------------
      ! [x, y, z]
      !---------------------------------------
      !
      !----------------------------------------------
      !find matrix for rotation about the x-axis
      !-----------------------------------------------
      !
      cosx = vector(3)/sqrt( vector(2)**(2.) + vector(3)**(2.) )
      sinx = vector(2)/sqrt( vector(2)**(2.) + vector(3)**(2.) )
      !
      call x_rotation_matrix (cosx, sinx, Rx)
      !        
      !--------------------------------------------------------
      !rotate space about the x - axis to bring the inertia vector in
      !the xz-plane:
      !---------------------------------------------------------
      !
      rotated1 = matmul(Rx,coords5)
      vecxz = matmul(Rx,vector1)
      vec1xz = matmul(Rx,vectora)
      write(*,*) 'vecxz', vecxz
      !
      !---------------------------------------------------
      !find matrix for rotation about the y-axis
      !---------------------------------------------------
      !
      cosy = vecxz(3,1)/sqrt( vecxz(1,1)**(2.) + vecxz(3,1)**(2.) )
      siny = vecxz(1,1)/sqrt( vecxz(1,1)**(2.) + vecxz(3,1)**(2.) )
      !
      call y_rotation_matrix (cosy, siny, Ry)
      !
      !-------------------------------------------------------------
      !rotate space around the y - axis, so that the rotation axis lies
      !along the positive z - axis
      !----------------------------------------------------------------
      !
      rotated2 = matmul(Ry,rotated1)
      vecz = matmul(Ry,vecxz) !should have the form [0,0,z]
      vec1z = matmul(Ry,vec1xz)
      write(*,*) 'vecz', vecz
      !
   end if
end if
!
do i = 1, n5
   mol2(i,1) = mol1(i,1) ! just transferring the atomic masses
end do
!
do i = 1, n5
   do j = 1, 3
      matrix5(i,j) = rotated2(j,i) ! take the transpose
   end do
end do
!
do i = 1, n5
   do j = 1, 3
      mol2(i,j+1) = matrix5(i,j)
   end do
end do
!
deallocate(rotated1, rotated2)
!
end subroutine align
!
!------------------------------------------------------------------
!subroutine in order to test if the molecular coordinates are equal to 
!to any previously generated set of coordinates:
!------------------------------------------------------------------
!
subroutine test2 (mass8, coords8, n8, c8, logic8)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
real(dp), allocatable, dimension(:), intent(in) :: mass8
real(dp), allocatable, dimension(:,:), intent(in) :: coords8
real(dp), allocatable, dimension(:,:) :: mold, mole
integer, intent(in) :: n8, c8
integer :: i, j, k
character(len=1), intent(out) :: logic8
character(len=3) :: file8
integer :: logic5
!
logic8 = 'Y'
!
!check if it is equal to any previous set of generated coordinates
!
if (c8 .gt. 150 )then
   allocate( mold (n8,4), mole(n8,4) )
   do j = 1, n8 !this structure:
      mole(j,1) = mass8(j)
      mole(j,2) = coords8(j,1)
      mole(j,3) = coords8(j,2)
      mole(j,4) = coords8(j,3)
   end do
   do i = 150, c8 - 1
      !
      write(file8, "(I3)") i
      open(unit=i, file=file8, status='old', action='read')
      do k = 1, n8 !previous structure:
         read(i,*) (mold (k,j), j = 1, 4)
      end do
      close(i)
      !
      call compare(mole, mold, n8, logic5)
      if (logic5 == 1) then 
         logic8 = 'N'
         exit
      end if
   end do 
   deallocate (mold, mole)
end if
!
end subroutine test2
!
!-----------------------------------------------------------------
!Subroutine to compare two matrices (molecular structures), to see
!if they are trully identical.
!-----------------------------------------------------------------
!
subroutine compare (matrix1, matrix2, natoms7, prod)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
integer :: i, j
integer, intent(in) :: natoms7
real(dp), allocatable, dimension(:,:), intent(in) :: matrix1, matrix2
integer, allocatable, dimension (:) :: comparison, B, BB ! comparison array
integer, intent(out) :: prod ! for the product for the comparison
!
allocate (comparison(natoms7), B(natoms7), BB(natoms7))
!
do i = 1, natoms7 !set them to 0 at the start of the comparison.
    B(i) = 0
    comparison(i) = 0
    BB(i) = 0
end do
!
!--------------------------------------------------------------------
!check if for every atom in matrix1, there is an atom in matrix2 that
!has the same set of coordinates:
!--------------------------------------------------------------------
!
do i = 1, natoms7
   do j = 1, natoms7
      if (j == B(j)) cycle
      if (i == BB(i)) cycle
      if (dabs( matrix1(i,1) - matrix2(j,1)) .LE. 0.0001) then
         if (dabs(matrix1(i,2) - matrix2(j,2)) .LE. 0.015) then
            if (dabs(matrix1(i,3) - matrix2(j,3)) .LE. 0.015) then
               if (dabs(matrix1(i,4) - matrix2(j,4)) .LE. 0.015) then
                  comparison(i) = 1
                  B(j) = j
                  BB(i) = i
                  exit
               end if
            end if
         end if
      end if
   end do
end do
!
!---------------------------------------------------------------------
!and now to see if all the matrix elements have been found to be equal
!----------------------------------------------------------------------
!
prod = 1
do i = 1, natoms7
    prod = prod * comparison(i)
end do 
!
deallocate(comparison, B, BB)
!
end subroutine compare 
!
!-----------------------------------------------------------------
!To determine if there are symmetrically equivalent atoms:
!-----------------------------------------------------------------
!
subroutine SEA (b_matrix, b_atoms, SEA_array)
implicit none 
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
real(dp), allocatable, dimension(:,:), intent(in) :: b_matrix
integer, intent(in) :: b_atoms !number of atoms
real(dp), allocatable, dimension(:,:) :: distance_m
integer :: i, k, j, l, h, f, prod1
integer, allocatable, dimension(:) :: equivalent, A, B, C
real(dp), allocatable, dimension(:) :: D1, D2
logical, allocatable, dimension(:), intent(out) ::  SEA_array
!
!-------------------------------------------------------------------
!Build distance matrix in order to find two symetrically equivalent atoms
!-------------------------------------------------------------
!
call distance_matrix (b_matrix, b_atoms, distance_m)
!
!-----------------------------------------
!find two symetrically equivalent atoms:
!-----------------------------------------
!
allocate (D1(b_atoms), D2(b_atoms), A(b_atoms), equivalent(b_atoms), &
  B(b_atoms), C(b_atoms), SEA_array(b_atoms) )

do i = 1, b_atoms
    C(i) = 0
    SEA_array(i) = .false.
end do
!
!First, select a column (i) from the distance matrix. The set fo SEA
!to i will be searched for. Then, another one, (i+1), and the set if
!SEA for this one will be searched for as well. 
!
determine_SEA: do i = 1, b_atoms - 1
   if (i == C(i)) cycle !It has been found to be a SEA
   do k = i + 1, b_atoms
      if (k == C(k)) cycle !It has been found to be a SEA from
!                             another set
      do h = 1, b_atoms
         equivalent(h) = 0 !return the variable used for comparison to 0.
         B(h) = 0
         A(h) = 0
      end do
!     Now, fill the D1 and D2 vectors used for the comparison with elements
!     from the i and k rows:
      do j = 1, b_atoms 
         do l = 1, b_atoms  
            D1(j) = distance_m(i,j)
            D2(l) = distance_m(k,l)
         end do
      end do
!     Compare these two columns of the distance matrix
!     to see if they have the same elements
      do j = 1, b_atoms
         do l = 1 , b_atoms
            if (j == A(j)) cycle !to not check a distance that has alread
            if (l == B(l)) cycle !been found to be equal for both atoms
            if ( dabs(D1(j)-D2(l)) .LE. 0.0001 ) then
               equivalent(j) = 1
               A(j) = j
               B(l) = l
               exit
            end if
         end do
      end do
      write(*,*) 'equivalent?', equivalent
      prod1 = 1 !to see if for every element in D1 a corresponding equal 
!                element in D2 has been found.
     do f = 1, b_atoms
         prod1 = prod1 * equivalent(f)
         write(*,*) 'prod', prod1
      end do 
      if (prod1 == 1) then
!        fix one of the SEA atoms to not have it's atom changed from 
!        B to V
         SEA_array(k) = .true. !will not be included
         C(k) = k
      end if
  end do
end do determine_SEA
!
write(*,*) 'SEA Atoms', SEA_array
!
end subroutine SEA
!
!---------------------------------------------------------------------
!Find the kth lexicographical combination
!in the set of combinations of m out of n
!---------------------------------------------------------------------
!
subroutine lexicographical_order (k_comb, ng, n18, combination18)
implicit none
!
integer, intent(in) :: k_comb, ng, n18
integer(kind=16) :: n16, j, test, previous_step
integer :: i, m_comb, l
integer, allocatable, dimension(:), intent(out) :: combination18
logical :: finish
!
allocate( combination18(n18) )
!
!------------------------
!First combination:
!------------------------
!
do i = 1, n18
   combination18(i) = i
end do
!
test = 1
i = 1 ! number of elements in the combination array [1 2 3 ...] of n
j = ng !number of elements in the set from which the combination is
!           selected (n-p)
n16=n18 !integer(kind=16)
finish = .false.
m_comb = 1
!
!---------------------------------
!The other combinations:
!---------------------------------
!
do while ( (i .le. n18) .and. (finish .eqv. .false.) )
   do while (j .ge. 0)
      m_comb = m_comb + 1 
      combination18(i) = m_comb
      previous_step = test
      !
      test = test + combination_x_of_p(n16-i,j-i)
      !
      if (test == k_comb) then
         do l = i + 1, n18
            combination18(l) = combination18(l-1) + 1
         end do
         finish = .true.
         exit
      else if (test .gt. k_comb) then 
         m_comb = m_comb - 1
         combination18(i) = m_comb
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
write(*,'(I5,A14,I5,A6,I5,A1)') k_comb, 'combination of', n18, &
   'out of', ng, ':'
write(*,*) combination18
!
!--------------------------------------------------
!example:
!--------------------------------------------------
!
!number of combinations of 5 out of 10 = 252
!lets find the 100th in lexicographical order.
!
!The very first combination is:
!
!1 2 3 4 5 
!
!and, then the rest are combinations of the type:
!
!1 . . . . . . . . . .
!
!until we get to combinations of the type:
!
!. 2 . . . . . . . . .
!
!But, how many of them are there?
!
!Well, if we fix m(1) = 1, there are nine elements left in n to 
!make combinations out of, and four free spaces in
!the set of elements of the 100th combination, so the number of
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
!
!      [ 1 . . . . ]
!
!By the same logic, we can find the second element in 
!the 100th combination as
!
!---------------------------------
!
!1st : 1 2 . . . . . . . . 
!
!---------------------------------
!
!c 3 of 8 = 56 + 1 = 57
!
!57th: 1 . 3 . . . . . . . 
!            
!--------------------------------------
!
!c 3 of 7 = 35 + 57 = 92
!
!92nd: 1 . . 4 . . . . . .
!
!-------------------------------------
!
!c 3 of 6 = 20
!
!112th : 1 . . . 5 . . . .
!
!-----------------------------------------
!
!           [ 1 4 . . . ]
!
!--------------------------------------------
!
!92nd: 1 . . 4 5 . . . . .
!
!------------------------------------------------
!
!c 2 of 5 = 10 + 92 = 102
!
!102nd 1 . . 4 . 6 . . . .
!
!------------------------------------------
!
!           [1 4 5 . . ]
!
!92nd: 1 . . 4 5 6 . . . .
!
!-----------------------------------------
!
!c 1 of 4 = 4
!
!96th 1 . . 4 5 . 7 . . .
!
!--------------------------------------
!
!c 1 of 3 = 3
!
!99th: 1 . . 4 5 . . 8 . .
!
!------------------------------------------
!
!c 1 of 2 = 2 + 99 = 101
!
!101st: 1 . . 4 5 . . . 9
!
!           [1 4 5 8 . ] 
!
!--------------------------------------------
!
!99th: 1 . . 4 5 . . 8 9 . 
!
!c 0 of 1  = 1
!
!--------------------------------------------------
!
!99 + 1 = 100
!
!100th: 1 . . 4 5 . . 8 . 10
!
!--------------------------------------------------
!
end subroutine lexicographical_order
!
!----------------------------------------------------------------
!to determine the direction of the walk
!----------------------------------------------------------------
!
real(kind=8) function random_walk(x, change)
implicit none
!
real(kind=8), intent(in) :: x, change
real(kind=8), dimension(1) :: walk
!
call random_number(walk)
if (walk(1) .lt. 0.5) then
   random_walk = x + change
else
   random_walk = x - change
end if
!
end function random_walk
!
end module global_module
