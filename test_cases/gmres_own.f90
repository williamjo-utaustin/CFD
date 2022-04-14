program gmres_own
     implicit none

     !call test_1()
     
     call test_2()


end program gmres_own

subroutine test_2()

      implicit none

      integer, parameter :: n = 1000
      integer, parameter :: nz_num = (n * 3) - 2

      double precision, dimension(n,n) :: A
      double precision, dimension(n) :: b
      double precision, dimension(n) :: x_r, x_est
      integer :: i, j, k

      double precision, dimension(nz_num) :: a_spr 
      integer, dimension(nz_num) :: r_spr, c_spr
      
      integer :: itr_max, mr
      double precision :: tol_abs, tol_rel
      
      itr_max = 2000
      mr = n - 1
      tol_abs = 1.0D-08
      tol_rel = 1.0D-08


      call create_random_matrix(n,A)
      call create_test_x(n,x_r) 
      call create_b(n,A,x_r,b)

      !do i = 1, n
      !  write(6,*) A(i,:), x_r(i), b(i)
      !end do 

      write(6,*) "The test matrix has been created. Continue to solve the sparse matrix for Ax=b"
      pause


      ! Code block below determines which part of each matrix is not zero
      !nz_num = 0

      !do i = 1, n
      !  do j = 1, n
      !     if(A(i,j) > 0) then
      !          nz_num = nz_num + 1
      !     end if
      !  end do
      !end do 
        
      !write(6,*) nz_num
      
      k = 0
      do i = 1, n
        do j = 1, n
           if(A(i,j) > 0) then
                k = k + 1
                a_spr(k) = A(i,j)
                r_spr(k) = i
                c_spr(k) = j
           end if
        end do
      end do 
      
      ! Create the sparse matrix
      do k = 1, nz_num
         write(6,*) a_spr(k), r_spr(k), c_spr(k)
      end do

      x_est = 0

      call timestamp()   
      call mgmres_st(n, nz_num, r_spr, c_spr, a_spr, x_est, b, itr_max, mr, tol_abs, tol_rel)
      call timestamp()   
      
      do i = 1, n
        write (6, '(2x,i8,2x,g14.6)' ) i, x_est(i)
      end do

end subroutine test_2



subroutine create_b(n,A,x_r,b)
      implicit none
      integer, intent(in) :: n
      double precision, dimension(n), intent(in) :: x_r
      double precision, dimension(n,n), intent(in) :: A
      double precision, dimension(n), intent(out) :: b

      integer :: i, j

      b = 0
      do i = 1,n
        do j = 1,n
           b(i) = b(i) + A(i,j) * x_r(j)
        end do
      end do

end subroutine create_b


subroutine create_test_x(n,x_r)
      implicit none
      integer, intent(in) :: n
      double precision, dimension(n), intent(out) :: x_r
      integer :: i
      
      do i = 1,n
        x_r(i) = i
      end do

end subroutine create_test_x

subroutine create_random_matrix(n,A)

      implicit none

      integer, intent(in) :: n
      double precision, dimension(n,n), intent(out) :: A
      double precision :: rand
      integer :: i,j
      
      A = 0D0

      do i = 1, n

        do j = i-1,i+1
             if(j>=1.and.j<=n) then
                call random_number(rand)
                rand = rand * 50
                A(i,j) = rand
             end if

        end do
      end do

end subroutine create_random_matrix

subroutine test_1()
      implicit none

      ! this is the order of the linear system - 9x9
      integer, parameter :: n = 9

      ! this is the number of non-zero matrix values
      integer, parameter :: nz_num = 23

      ! this is the matrix vlaues excluding zeros
      double precision, dimension(nz_num) :: a = (/ &
         2.0D+00, -1.0D+00, &
         2.0D+00, -1.0D+00, &
        -1.0D+00,  2.0D+00, &
        -1.0D+00,  2.0D+00, -1.0D+00, &
        -1.0D+00,  2.0D+00, -1.0D+00, &
        -1.0D+00,  2.0D+00, -1.0D+00, &
        -1.0D+00,  2.0D+00, -1.0D+00, &
        -1.0D+00,  2.0D+00, -1.0D+00, &
        -1.0D+00,  2.0D+00 /)
      

      ! these are the row indices of the matrix values
      integer, dimension(nz_num) :: ia = (/ &
         1, 1, &
         2, 2, &
         3, 3, &
         4, 4, 4, &
         5, 5, 5, &
         6, 6, 6, &
         7, 7, 7, &
         8, 8, 8, &
         9, 9 /)

      ! these are the column indices of the matrix values
      integer, dimension(nz_num) :: ja = (/ &
         1, 4, &
         2, 3, &
         2, 3, &
         1, 4, 5, &
         4, 5, 6, &
         5, 6, 7, &
         6, 7, 8, &
         7, 8, 9, &
         8, 9 /)

      ! this is the b value in Ax = b
      double precision, dimension(n) :: rhs = (/ &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00 /)
 
      ! the x value you want in Ax = b
      double precision, dimension(n) :: x_estimate

      ! the maximum number of iterations to take  (itr_max)
      ! the maximum number of inner iterations to take
      integer :: itr_max, mr
      double precision :: tol_abs, tol_rel
      integer :: i


      itr_max = 20
      mr = n - 1
      tol_abs = 1.0D-08
      tol_rel = 1.0D-08
    
      call mgmres_st(n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr, tol_abs, tol_rel)

      do i = 1, n
        write (6, '(2x,i8,2x,g14.6)' ) i, x_estimate(i)
      end do


end subroutine test_1
