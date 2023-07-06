module geometry
      implicit none


      contains

      subroutine create_grid(N_x, N_y, L_x, L_y, h_x, h_y)

                ! number of grids in x and y
                integer, intent(in) :: N_x, N_y 

                ! length of domain in x and y
                integer, intent(in) :: L_x, L_y

                ! length of a grid in x and y
                double precision, intent(out) :: h_x, h_y

                h_x = dble(L_x)/dble(N_x)
                h_y = dble(L_y)/dble(N_y)

      end subroutine create_grid

end module geometry
