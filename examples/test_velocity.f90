module program_variables

      implicit none

      ! number of grids in x and y
      integer, parameter :: N_x = 4 ! columns - j
      integer, parameter :: N_y = 4 ! rows - i

      ! length of domain in x and y
      double precision, parameter :: L_x = 1D0
      double precision, parameter :: L_y = 1D0
 


end module program_variables


program test1_velocity

      use program_variables

      implicit none
      double precision :: h_x, h_y
      
      double precision, dimension(:), allocatable :: x_mesh_1
      double precision, dimension(:), allocatable :: y_mesh_1
      
      double precision, dimension(:), allocatable :: x_mesh_2
      double precision, dimension(:), allocatable :: y_mesh_2
      
      double precision, dimension(:,:), allocatable :: u, v
      double precision, dimension(:,:), allocatable :: u_h, v_h
      double precision, dimension(:,:), allocatable :: N_u, N_v
      integer :: i, j
      
      ! main grid will be size N, total grid will be size N+3 to account for boundaries
      allocate(x_mesh_1(0:N_x+2))
      allocate(y_mesh_1(0:N_y+2))
      
      ! staggered grid will be size N, total grid will be size N+2 to account for boundaries 
      allocate(x_mesh_2(0:N_x+1))
      allocate(y_mesh_2(0:N_y+1))
      
      allocate(u(0:N_y + 2, 0:N_x + 2))
      allocate(v(0:N_y + 2, 0:N_x + 2))
      allocate(u_h(0:N_y+1, 0:N_x+1))
      allocate(v_h(0:N_y+1, 0:N_x+1))
      allocate(N_u(0:N_y+1, 0:N_x+1))
      allocate(N_v(0:N_y+1, 0:N_x+1))
      
      u = 20D0
      v = 3D0

      
      call create_grid(h_x, h_y, x_mesh_1, y_mesh_1, x_mesh_2, y_mesh_2)

      write(6,*) h_x, h_y
      
      write(6,*) "x_mesh ", x_mesh_1
      write(6,*) "y_mesh ", y_mesh_1
      
      write(6,*) "x_mesh_hat", x_mesh_2
      write(6,*) "y_mesh_hat", y_mesh_2
        
      call convective_operator(u, v, h_x, h_y, N_u, N_v)
      
      do i = 0, N_y + 1
            write(6,*) N_u(i,:) 
      end do
            
      !write(6,*) N_v(i,:) 
      
      ! main grid will be size N, total grid will be size N+3 to account for boundaries
      deallocate(x_mesh_1)
      deallocate(y_mesh_1)
      
      ! staggered grid will be size N, total grid will be size N+2 to account for boundaries 
      deallocate(x_mesh_2)
      deallocate(y_mesh_2)
      
      deallocate(u)
      deallocate(v)
      deallocate(u_h)
      deallocate(v_h)
      deallocate(N_u)
      deallocate(N_v)

end program test1_velocity

subroutine create_grid(h_x, h_y, x_mesh_1, y_mesh_1, x_mesh_2, y_mesh_2)

          use program_variables

          ! length of a grid in x and y
          double precision, intent(out) :: h_x, h_y

          double precision, dimension(0:N_x+2), intent(out) :: x_mesh_1
          double precision, dimension(0:N_y+2), intent(out) :: y_mesh_1
          
          double precision, dimension(0:N_x+1), intent(out) :: x_mesh_2
          double precision, dimension(0:N_y+1), intent(out) :: y_mesh_2

          integer :: i, j

          h_x = L_x/dble(N_x)
          h_y = L_y/dble(N_y)


          do i = 0, N_x + 2
                x_mesh_1(i) = (i-1) * h_x      
          end do

          do j = 0, N_y + 2      
                y_mesh_1(j) = (j-1) * h_y      
          end do

          do i = 0, N_x + 1
                x_mesh_2(i) = (x_mesh_1(i) + x_mesh_1(i+1))/2
          end do
          
          do j = 0, N_y + 1
                y_mesh_2(j) = (y_mesh_1(j) + y_mesh_1(j+1))/2
          end do



end subroutine create_grid

subroutine vel_hat(u, v, u_h, v_h)

      use program_variables

      double precision, dimension(0:N_y + 2, 0:N_x + 2), intent(in) :: u, v
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(out) :: u_h, v_h

      integer :: i, j

      do j = 0, N_x + 1
            do i = 0, N_y + 1
                  u_h(i,j) = (u(i,j) + u(i+1,j))/2D0
                  v_h(i,j) = (v(i,j) + u(i,j+1))/2D0
            end do
      end do


end subroutine vel_hat



subroutine convective_operator(u, v, h_x, h_y, N_u, N_v)

      use program_variables

      implicit none
      double precision, dimension(0:N_y+2, 0:N_x+2), intent(in) :: u, v
      double precision, intent(in) :: h_x, h_y
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(out) :: N_u, N_v

      integer :: i,j
      double precision :: u_e, u_w, u_n, u_s
      double precision :: v_e, v_w, v_n, v_s

      do i = 1, N_y
            do j = 2, N_x
                  
                  u_e = 0.5 * (u(i,j+1) + u(i,j))
                  u_w = 0.5 * (u(i,j) + u(i,j-1))
                  u_n = 0.5 * (u(i+1,j) + u(i,j))
                  u_s = 0.5 * (u(i,j) + u(i-1,j))
           
                  v_n = 0.5 * (v(i+1,j-1) + v(i+1,j))
                  v_s = 0.5 * (v(i,j-1) + v(i,j))
           
                  N_u(i,j) = (u_e * u_e - u_w * u_w) * h_y + (v_n * u_n - v_s * u_s) * h_x
                  
                  write(6,*) "i = row ", "j = column"
                  write(6,*) "u(",i,",",j,")"
                  write(6,*) "u_e","u(",i,",",j+1,")"
                  write(6,*) "u_w","u(",i,",",j-1,")"
                  write(6,*) "u_n","u(",i+1,",",j,")"
                  write(6,*) "u_s","u(",i-1,",",j,")"
                  write(6,*) "v(",i+1,',',j-1,") ,v(",i+1,',',j,")"
                  write(6,*) "v(",i,',',j-1,") ,v(",i,',',j,")"

            end do
      end do


      do i = 2, N_y
            do j = 1, N_x

                  write(6,*) i,j
                  u_e = 0.5 * (u(i,j+1) + u(i-1,j+1))
                  u_w = 0.5 * (u(i-1,j) + u(i,j))
                  v_e = 0.5 * (v(i,j) + v(i,j+1))
                  v_w = 0.5 * (v(i,j) + v(i,j-1))
                  v_n = 0.5 * (v(i,j) + v(i+1,j))
                  v_s = 0.5 * (v(i,j) + v(i-1,j))
                  
                  write(6,*) "i = row ", "j = column"
                  write(6,*) "u(",i,",",j,")"
                  write(6,*) "v_e","v(",i,",",j+1,")"
                  write(6,*) "v_w","u(",i,",",j-1,")"
                  write(6,*) "v_n","v(",i+1,",",j,")"
                  write(6,*) "v_s","v(",i-1,",",j,")"
                  write(6,*) "u(",i,',',j+1,") ,u(",i-1,',',j+1,")"
                  write(6,*) "u(",i,',',j,") ,u(",i-1,',',j,")"

                  N_v(i,j) = (u_e * v_e - u_w * v_w) * h_y + (v_n * v_n - v_s * v_s) * h_x

            end do
      end do




end subroutine convective_operator



!subroutine convective_operator(u, v, u_h, v_h, N_u, N_v)
!
!      use program_variables
!
!      implicit none
!      double precision, dimension(0:N_x+2, 0:N_y+2), intent(in) :: u, v
!      double precision, dimension(0:N_x+1, 0:N_y+1), intent(in) :: u_h, v_h
!      double precision, dimension(0:N_x+1, 0:N_y+1), intent(out) :: N_u, N_v
!
!      integer :: i, j
!
!      N_u = 1
!      N_v = 1
!
!      do i = 1, N_x + 1
!            do j = 1, N_y + 1
!                  N_u(i,j) = -u(i,j) * (u(i+1,j) - u(i-1,j)) - ((v_h(i,j) + v_h(i-1,j))/2) * (u(i,j+1) - u(i,j-1)) 
!                  N_v(i,j) = -v(i,j) * (v(i,j+1) - u(i,j-1)) - ((u_h(i,j) + u_h(i,j-1))/2) * (v(i+1,j) - v(i-1,j)) 
!            end do
!      end do
!end subroutine convective_operator
!
