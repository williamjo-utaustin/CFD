module program_variables

      implicit none

      ! number of grids in x and y
      integer, parameter :: N_x = 5
      integer, parameter :: N_y = 5

      ! length of domain in x and y
      double precision, parameter :: L_x = 5D0
      double precision, parameter :: L_y = 5D0
 


end module program_variables


program test1_velocity

      use program_variables

      implicit none
      double precision :: h_x, h_y
      
      ! main grid will be size N, total grid will be size N+3 to account for boundaries
      double precision, dimension(0:N_x+2) :: x_mesh_1
      double precision, dimension(0:N_y+2) :: y_mesh_1
      
      ! staggered grid will be size N, total grid will be size N+2 to account for boundaries 
      double precision, dimension(0:N_x+1) :: x_mesh_2
      double precision, dimension(0:N_y+1) :: y_mesh_2
      
      double precision, dimension(0:N_x + 2, 0:N_y + 2) :: u, v
      double precision, dimension(0:N_x+1, 0:N_y+1) :: u_hat, v_hat
      double precision, dimension(0:N_x+1, 0:N_y+1) :: N_u, N_v
      integer :: i, j
      u = 20D0
      v = 3D0

      
      call create_grid(h_x, h_y, x_mesh_1, y_mesh_1, x_mesh_2, y_mesh_2)

      !write(6,*) h_x, h_y
      
      write(6,*) "x_mesh ", x_mesh_1
      write(6,*) "y_mesh ", y_mesh_1
      
      write(6,*) "x_mesh_hat", x_mesh_2
      write(6,*) "y_mesh_hat", y_mesh_2
        
      call vel_hat(u, v, u_hat, v_hat)
      
      !do i = 0, N_x + 1
      !    write(6,*) u_hat(i,:)
      !end do
      
      !do i = 0, N_x + 1
      !    write(6,*) v_hat(i,:)
      !end do

      call convective_operator(u, v, u_hat, v_hat, N_u, N_v)
      do i = 0, N_x + 1
            write(6,*) N_u(i,:) 
      end do
            
      !write(6,*) N_v(i,:) 

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

subroutine vel_hat(u, v, u_hat, v_hat)

      use program_variables

      double precision, dimension(0:N_x + 2, 0:N_y + 2), intent(in) :: u, v
      double precision, dimension(0:N_x+1, 0:N_y+1), intent(out) :: u_hat, v_hat

      integer :: i, j

      do i = 0, N_x + 1
            do j = 0, N_y + 1
                  u_hat(i,j) = (u(i,j) + u(i+1,j))/2D0
                  v_hat(i,j) = (v(i,j) + u(i,j+1))/2D0
            end do
      end do


end subroutine vel_hat


subroutine convective_operator(u, v, u_hat, v_hat, N_u, N_v)

      use program_variables

      implicit none
      double precision, dimension(0:N_x+2, 0:N_y+2), intent(in) :: u, v
      double precision, dimension(0:N_x+1, 0:N_y+1), intent(in) :: u_hat, v_hat
      double precision, dimension(0:N_x+1, 0:N_y+1), intent(out) :: N_u, N_v

      integer :: i, j

      N_u = 1
      N_v = 1

      do i = 1, N_x + 1
            do j = 1, N_y + 1
                  N_u(i,j) = -u(i,j) * (u(i+1,j) - u(i-1,j)) - ((v_hat(i,j) + v_hat(i-1,j))/2) * (u(i,j+1) - u(i,j-1)) 
                  N_v(i,j) = -v(i,j) * (v(i,j+1) - u(i,j-1)) - ((u_hat(i,j) + u_hat(i,j-1))/2) * (v(i+1,j) - v(i-1,j)) 
            end do
      end do
end subroutine convective_operator

