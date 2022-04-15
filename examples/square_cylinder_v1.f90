module program_variables

      implicit none

      ! number of grids in x and y
      integer, parameter :: N_x = 4 ! columns - j
      integer, parameter :: N_y = 4 ! rows - i

      ! length of domain in x and y
      double precision, parameter :: L_x = 1D0
      double precision, parameter :: L_y = 1D0


      double precision, parameter :: nu = 1D0

      double precision, parameter :: delta_t = 0.01
 

      double precision, parameter :: u_inlet = 10
      double precision, parameter :: v_inlet = 0

end module program_variables


program test1_velocity

      use program_variables

      implicit none
      double precision :: h_x, h_y
      
      double precision, dimension(:), allocatable :: x_mesh_1
      double precision, dimension(:), allocatable :: y_mesh_1
      
      double precision, dimension(:), allocatable :: x_mesh_2
      double precision, dimension(:), allocatable :: y_mesh_2
      
      double precision, dimension(:,:), allocatable :: u, v, p
      double precision, dimension(:,:), allocatable :: u_temp, v_temp
      double precision, dimension(:,:), allocatable :: u_h, v_h
      double precision, dimension(:,:), allocatable :: N_u, N_v
      double precision, dimension(:,:), allocatable :: D_u, D_v
      integer :: i, j

      
      ! main grid will be size N, total grid will be size N+3 to account for boundaries
      allocate(x_mesh_1(0:N_x+2))
      allocate(y_mesh_1(0:N_y+2))
      
      ! staggered grid will be size N, total grid will be size N+2 to account for boundaries 
      allocate(x_mesh_2(0:N_x+1))
      allocate(y_mesh_2(0:N_y+1))
      
      allocate(u(0:N_y + 2, 0:N_x + 2))
      allocate(v(0:N_y + 2, 0:N_x + 2))
      allocate(p(0:N_y + 2, 0:N_x + 2))
      allocate(u_temp(0:N_y + 2, 0:N_x + 2))
      allocate(v_temp(0:N_y + 2, 0:N_x + 2))
      allocate(u_h(0:N_y+1, 0:N_x+1))
      allocate(v_h(0:N_y+1, 0:N_x+1))
      allocate(N_u(0:N_y+1, 0:N_x+1))
      allocate(N_v(0:N_y+1, 0:N_x+1))
      allocate(D_u(0:N_y+1, 0:N_x+1))
      allocate(D_v(0:N_y+1, 0:N_x+1))
     
      
      ! create initial conditions of the domain
      u = 0
      v = 0
      p = 0
      ! Loop here!

      ! maintain boundary conditions
      ! create the inlet velocity
      
      u(1:N_y,1) = u_inlet
      v(2:N_y,0) = 2*v_inlet - v(2:N_y,1)
      
      ! top and bottom wall
      v(1,1:N_x+1) = 0
      v(N_y+1,1:N_x+1) = 0
     

      !! velocities to test if boundary conditions work
      !! comment out when finished
      !! test for bottom wall
      !u(1,2:N_x+1) = 1
      !! test for top wall
      !u(N_y,2:N_x+1) = 2
      !! test for right outlet
      !u(2:N_y - 1,N_x) = 3     
      !v(2:N_y,1) = 4
      !! test for left boundary
      !v(2:N_y,N_x) = 5

      ! bottom wall ghost boundaries (u = 0 -> Dirichlet)
      u(0, 2:N_x) = -u(1, 2:N_x)
      
      ! top wall ghost boundary (u = 0 -> Dirichlet)
      u(N_y+1,2:N_x) = -u(N_y, 2:N_x)
      
      ! right outlet ghost boundary (dv/dx = 0 -> Neumann)
      v(2:N_y,N_x+1) = v(2:N_y,N_x)

      ! right outlet boundary (du/dx = 0 -> Neumann)
      u(2:N_y-1,N_x+1) = u(2:N_y - 1,N_x) 


      do i = 0, N_y + 2
            write(6,*) u(i,:) 
      end do
     
      write(6,*)
      do i = 0, N_y + 2
            write(6,*) v(i,:) 
      end do

      pause
      
      call create_grid(h_x, h_y, x_mesh_1, y_mesh_1, x_mesh_2, y_mesh_2)

      write(6,*) h_x, h_y
      
      write(6,*) "x_mesh ", x_mesh_1(0:N_x + 1)
      write(6,*) "y_mesh ", y_mesh_1(0:N_y +1 )
      
      write(6,*) "x_mesh_hat", x_mesh_2
      write(6,*) "y_mesh_hat", y_mesh_2
       

      ! call convective and diffusive operators for only the inner cells 
      ! create u* and v* for inner cells only
      call convective_operator(u, v, h_x, h_y, N_u, N_v)
      call diffusive_operator(u, v, h_x, h_y, D_u, D_v)
      call vel_temp(u, v, h_x, h_y, N_u, N_v, D_u, D_v, u_temp, v_temp)


      ! constrain boundary conditions
      ! use u* and v*, along with initial assumption of p to get u* and v* near or at boundaries
      ! create inlet u* (Dirichlet Condition)
      call u_temp_inlet(u_temp, h_x, h_y, p)



      ! create outlet u* (Neumann Condition)

      ! create top and bottom boundaries for u* (Dirichlet Condition)


      ! create inlet v*

      ! create outlet v* 

      ! create top and bottom boundaries for v*




      do i = 0, N_y + 1
            write(6,*) u_temp(i,0:N_x+1) 
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
      deallocate(p)
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
                  
                  !write(6,*) "i = row ", "j = column"
                  !write(6,*) "u(",i,",",j,")"
                  !write(6,*) "u_e","u(",i,",",j+1,")"
                  !write(6,*) "u_w","u(",i,",",j-1,")"
                  !write(6,*) "u_n","u(",i+1,",",j,")"
                  !write(6,*) "u_s","u(",i-1,",",j,")"
                  !write(6,*) "v(",i+1,',',j-1,") ,v(",i+1,',',j,")"
                  !write(6,*) "v(",i,',',j-1,") ,v(",i,',',j,")"

            end do
      end do


      do i = 2, N_y
            do j = 1, N_x

                  u_e = 0.5 * (u(i,j+1) + u(i-1,j+1))
                  u_w = 0.5 * (u(i-1,j) + u(i,j))
                  
                  v_e = 0.5 * (v(i,j) + v(i,j+1))
                  v_w = 0.5 * (v(i,j) + v(i,j-1))
                  v_n = 0.5 * (v(i,j) + v(i+1,j))
                  v_s = 0.5 * (v(i,j) + v(i-1,j))
                  
                  !write(6,*) "i = row ", "j = column"
                  !write(6,*) "u(",i,",",j,")"
                  !write(6,*) "v_e","v(",i,",",j+1,")"
                  !write(6,*) "v_w","u(",i,",",j-1,")"
                  !write(6,*) "v_n","v(",i+1,",",j,")"
                  !write(6,*) "v_s","v(",i-1,",",j,")"
                  !write(6,*) "u(",i,',',j+1,") ,u(",i-1,',',j+1,")"
                  !write(6,*) "u(",i,',',j,") ,u(",i-1,',',j,")"

                  N_v(i,j) = (u_e * v_e - u_w * v_w) * h_y + (v_n * v_n - v_s * v_s) * h_x

            end do
      end do

end subroutine convective_operator


subroutine diffusive_operator(u, v, h_x, h_y, D_u, D_v)

      use program_variables

      implicit none

      double precision, dimension(0:N_y+2, 0:N_x+2), intent(in) :: u, v
      double precision, intent(in) :: h_x, h_y
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(out) :: D_u, D_v
     
      integer :: i, j

      do i = 1, N_y
            do j = 2, N_x

                  D_u(i,j) = ((u(i,j+1) - 2 * u(i,j) + u(i,j-1))*h_y/h_x) + ((u(i+1,j) - 2*u(i,j) + u(i-1,j)) * h_x/h_y)
                  
                  !write(6,*) "------"
                  !write(6,*) i, j
                  !write(6,*) "------"
                  !write(6,*) i, j+1
                  !write(6,*) i, j
                  !write(6,*) i, j-1
                  !write(6,*) i+1, j
                  !write(6,*) i, j
                  !write(6,*) i-1, j
                  !write(6,*) "------"
            end do
      end do

      
      do i = 2, N_y
            do j = 1, N_x
                  D_v(i,j) = ((v(i,j+1) - 2 * v(i,j) + v(i,j-1))*h_y/h_x) + ((v(i+1,j) - 2*v(i,j) + v(i-1,j)) * h_x/h_y)

            end do
      end do

end subroutine diffusive_operator

subroutine vel_temp(u, v, h_x, h_y, N_u, N_v, D_u, D_v, u_temp, v_temp)

      use program_variables
      
      implicit none

      double precision, dimension(0:N_y+2, 0:N_x+2), intent(in) :: u, v
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(in) :: N_u, N_v
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(in) :: D_u, D_v
      double precision, intent(in) :: h_x, h_y
      
      double precision, dimension(0:N_y+2, 0:N_x+2), intent(out) :: u_temp, v_temp

      integer :: i, j
      
      do i = 1, N_y
            do j = 1, N_x
                 write(6,*) i,j 
                 u_temp(i,j) = u(i,j) + (delta_t/(h_x * h_y)) * (-N_u(i,j) + nu * D_u(i,j))
            end do
      end do

      do i = 1, N_y
            do j = 1, N_x
                 write(6,*) i,j 
                 v_temp(i,j) = v(i,j) + (delta_t/(h_x * h_y)) * (-N_v(i,j) + nu * D_v(i,j))
            end do
      end do

end subroutine vel_temp

subroutine u_temp_inlet(u_temp, h_x, h_y, p)
     
    use program_variables
    
    implicit none

    double precision, dimension(0:N_y+2,0:N_x+2), intent(inout) :: u_temp
    double precision, dimension(0:N_y+2,0:N_x+2), intent(in) :: p
    double precision, intent(in) :: h_x, h_y

    integer :: i, j
    double precision :: G

    do i = 1, N_y

        call g_operator(p, h_x, h_y, i, 1, 1, G)
        u_temp(i,1) = u_inlet + delta_t * G

    end do



end subroutine u_temp_inlet

subroutine g_operator(p, h_x, h_y, i, j, e_x, G)
   
    use program_variables

    implicit none

    double precision, dimension(0:N_y+2,0:N_x+2), intent(in) :: p
    integer, intent(in) :: i, j, e_x
    double precision, intent(out) :: G
    double precision, intent(in) :: h_x, h_y


    ! x-axis for G
    if(e_x.eq.1) then
           g = (p(i,j) - p(i,j-1))/h_x
    else
           g = (p(i,j) - p(i-1,j))/h_y
    end if
     
    
end subroutine g_operator



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
