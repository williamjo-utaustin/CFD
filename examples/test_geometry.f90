module program_variables

        implicit none

          ! number of grids in x and y
          integer, parameter :: N_x = 10
          integer, parameter :: N_y = 15

          ! length of domain in x and y
          double precision, parameter :: L_x = 10D0
          double precision, parameter :: L_y = 15D0

end module program_variables


program test1_geometry

        use program_variables

        implicit none
        double precision :: h_x, h_y
        double precision, dimension(N_x + 1) :: x_mesh_1
        double precision, dimension(N_y + 1) :: y_mesh_1
        
        double precision, dimension(N_x) :: x_mesh_2
        double precision, dimension(N_y) :: y_mesh_2
        
        call create_grid(h_x, h_y, x_mesh_1, y_mesh_1, x_mesh_2, y_mesh_2)

        !write(6,*) h_x, h_y
        
        write(6,*) x_mesh_1
        write(6,*) y_mesh_1
        
        write(6,*) x_mesh_2
        write(6,*) y_mesh_2
        

end program test1_geometry

subroutine create_grid(h_x, h_y, x_mesh_1, y_mesh_1, x_mesh_2, y_mesh_2)

          use program_variables

          ! length of a grid in x and y
          double precision, intent(out) :: h_x, h_y

          double precision, dimension(N_x+1), intent(out) :: x_mesh_1
          double precision, dimension(N_y+1), intent(out) :: y_mesh_1
          
          double precision, dimension(N_x), intent(out) :: x_mesh_2
          double precision, dimension(N_y), intent(out) :: y_mesh_2

          integer :: i, j

          h_x = L_x/dble(N_x)
          h_y = L_y/dble(N_y)


          do i = 1, N_x + 1
                x_mesh_1(i) = (i-1) * h_x      
          end do

          do j = 1, N_y + 1      
                y_mesh_1(j) = (j-1) * h_y      
          end do

          do i = 1, N_x
                x_mesh_2(i) = (x_mesh_1(i) + x_mesh_1(i+1))/2
          end do
          
          do j = 1, N_y
                y_mesh_2(j) = (y_mesh_1(j) + y_mesh_1(j+1))/2
          end do



end subroutine create_grid

