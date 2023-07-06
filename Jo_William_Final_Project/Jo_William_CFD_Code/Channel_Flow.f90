! Channel Flow PDE Solver
!-----------------------------------------------------------------------
! This program computes channel flow in 2D.
! All subroutines and all variables are contained in this code aside from GMRES.f90
!
! By default the equations for x-velocity (U), y-velocity (V) and
! pressure correction (p') are solved. 
! To compile the program type "gfortran Channel_Flow.f90 gmres.f90 -o channel"
! Output files: - channel
! To run type: "./channel"
! Use Tecplot or Paraview for Visualization/Post-Processing
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Disclaimers and Notes
!-----------------------------------------------------------------------
! i is the index for rows, j is the index for columns
! Any array would read Array(Rows, Columns)
! Nx would be a count for columns 
! Ny would be a count for rows
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Default Values for Low Order Method and Kim and Moin
!-----------------------------------------------------------------------
! Default Channel Flow Variables
! L_x = 0.2D0
! L_y = 0.1D0
! delta_t = 0.00001
! Re = 10D0
!
! Default Values for Low Order:
! N_x = 25, N_y = 50, kim_moin = .False. 
!
! Default Values for KM:
! N_x = 50, N_y = 50, kim_moin = .True.
!-----------------------------------------------------------------------

module program_variables

      implicit none

      ! ------------------------------------------------------------------
      ! Modify these variables to change grid and flow geometry dimensions
      ! ------------------------------------------------------------------
      integer, parameter :: N_x = 50 ! columns - j
      integer, parameter :: N_y = 50 ! rows - i
      logical, parameter :: kim_moin = .True.

      double precision, parameter :: L_x = 0.2D0
      double precision, parameter :: L_y = 0.1D0
      
      double precision, parameter :: delta_t = 0.00001
      double precision, parameter :: Re = 10D0
      
      ! ------------------------------------------------------------------
      ! Do not touch these variables -------------------------------------
      ! ------------------------------------------------------------------
      
      double precision, parameter :: u_inlet = Re
      double precision, parameter :: v_inlet = 0      
      double precision, parameter :: t_step_trigger = 3
      integer, parameter :: N_eq = (N_x + 2) * (N_y + 2)
      integer, parameter :: N_eq_u_tmp = (N_x + 1) * (N_y + 2)
      integer, parameter :: N_eq_v_tmp = (N_x + 2) * (N_y + 1)
      double precision, parameter :: h_x = L_x/dble(N_x)
      double precision, parameter :: h_y = L_y/dble(N_y)
      double precision, parameter :: c_1 = 1D0/(h_x**2)
      double precision, parameter :: c_2 = 1D0/(h_y**2)
      double precision, parameter :: c_3 = 2 * (c_1 + c_2)
      double precision, parameter :: d_1 = 1D0/h_x
      double precision, parameter :: d_2 = 1D0/h_y
      double precision, parameter :: km = -(delta_t/2D0) * (1/Re)
      integer, parameter :: corner_eq = 4
      integer, parameter :: top_and_bottom_eq = 2 * (2 * N_x)
      integer, parameter :: inlet_eq = 2 * N_y
      integer, parameter :: outlet_eq = 2 * N_y
      integer, parameter :: int_eq = 5 * (N_x * N_y)
      integer, parameter :: nz_num = corner_eq + top_and_bottom_eq + inlet_eq + outlet_eq + int_eq
      integer, parameter :: inlet_eq_u_tmp = N_y
      integer, parameter :: outlet_eq_u_tmp = 3 * N_y
      integer, parameter :: t_and_b_eq_u_tmp = 2 * ((N_x-1) * 2)
      integer, parameter :: int_eq_u_tmp = 5 * (N_y * (N_x-1))
      integer, parameter :: nz_num_u_tmp = corner_eq + inlet_eq_u_tmp + outlet_eq_u_tmp + t_and_b_eq_u_tmp + int_eq_u_tmp
      integer, parameter :: inlet_eq_v_tmp = 2 * (N_y - 1)
      integer, parameter :: outlet_eq_v_tmp = 2 * (N_y - 1)
      integer, parameter :: t_and_b_eq_v_tmp = 2 * N_x
      integer, parameter :: int_eq_v_tmp = 5 * ((N_y - 1) * N_x)
      integer, parameter :: nz_num_v_tmp = corner_eq + inlet_eq_v_tmp + outlet_eq_v_tmp + t_and_b_eq_v_tmp + int_eq_v_tmp
      integer, parameter :: itr_max = 10000
      integer, parameter :: mr = N_eq - 1
      integer, parameter :: mr_u_tmp = N_eq_u_tmp - 1
      integer, parameter :: mr_v_tmp = N_eq_v_tmp - 1
      double precision, parameter :: tol_abs = 1.0D-03
      double precision, parameter :: tol_rel = 1.0D-03

end module program_variables

program channel_flow

      use program_variables

      implicit none
      
      double precision, dimension(:), allocatable :: x_mesh_1
      double precision, dimension(:), allocatable :: y_mesh_1
      double precision, dimension(:), allocatable :: x_h
      double precision, dimension(:), allocatable :: y_h
      double precision, dimension(:,:), allocatable :: u, v, p
      double precision, dimension(:,:), allocatable :: u_tmp, v_tmp
      double precision, dimension(:,:), allocatable :: u_h, v_h
      double precision, dimension(:,:), allocatable :: N_u, N_v
      double precision, dimension(:,:), allocatable :: D_u, D_v
      double precision, dimension(:,:), allocatable :: u_nm1, v_nm1
      double precision, dimension(:,:), allocatable :: u_h_nm1, v_h_nm1
      integer, dimension(:,:), allocatable :: spr_index, spr_index_rev
      integer, dimension(:,:), allocatable :: spr_km_u_tmp, spr_km_v_tmp, spr_rev_km_u_tmp, spr_rev_km_v_tmp
      integer, dimension(:), allocatable :: r_spr, c_spr, r_spr_u_tmp, r_spr_v_tmp, c_spr_u_tmp, c_spr_v_tmp
      double precision, dimension(:), allocatable :: a_spr, a_spr_u_tmp, a_spr_v_tmp
      double precision, dimension(:), allocatable :: P_est, u_tmp_spr, v_tmp_spr
      double precision, dimension(:), allocatable :: b_spr, b_spr_u_tmp, b_spr_v_tmp
      
      integer :: i, t_step
      double precision :: t_sim
      logical :: skip_to_k_m
      double precision :: rho
      integer :: itr_used

      ! main grid will be size N, total grid will be size N+3 to account for boundaries
      ! staggered grid will be size N, total grid will be size N+2 to account for boundaries 
      allocate(x_mesh_1(0:N_x+2))
      allocate(y_mesh_1(0:N_y+2))
      allocate(x_h(0:N_x+1))
      allocate(y_h(0:N_y+1))
      allocate(u(0:N_y + 2, 0:N_x + 2))
      allocate(v(0:N_y + 2, 0:N_x + 2))
      allocate(p(0:N_y + 2, 0:N_x + 2))
      allocate(spr_index(0:N_y + 1, 0:N_x + 1))
      allocate(spr_index_rev(N_eq, 2))
      allocate(u_tmp(0:N_y + 2, 0:N_x + 2))
      allocate(v_tmp(0:N_y + 2, 0:N_x + 2))
      allocate(u_h(0:N_y+1, 0:N_x+1))
      allocate(v_h(0:N_y+1, 0:N_x+1))
      allocate(N_u(0:N_y+1, 0:N_x+1))
      allocate(N_v(0:N_y+1, 0:N_x+1))
      allocate(D_u(0:N_y+1, 0:N_x+1))
      allocate(D_v(0:N_y+1, 0:N_x+1))
      allocate(b_spr(N_eq))
      allocate(P_est(N_eq))
      allocate(a_spr(nz_num))
      allocate(r_spr(nz_num))
      allocate(c_spr(nz_num))
      
      ! create sparse matrix for A => Aspr, rspr, cspr
      call create_sparse_index(spr_index, spr_index_rev)
      call sparse_vectors(spr_index, r_spr, c_spr, a_spr)
      
      if(kim_moin.eqv..True.) then
            
            allocate(spr_km_u_tmp(0:N_y + 1, 1:N_x + 1)) 
            allocate(spr_km_v_tmp(1:N_y + 1, 0:N_x + 1)) 
            
            allocate(spr_rev_km_u_tmp(N_eq_u_tmp, 2)) 
            allocate(spr_rev_km_v_tmp(N_eq_v_tmp, 2)) 
            
            allocate(a_spr_u_tmp(nz_num_u_tmp))
            allocate(r_spr_u_tmp(nz_num_u_tmp))
            allocate(c_spr_u_tmp(nz_num_u_tmp))
            allocate(b_spr_u_tmp(N_eq_u_tmp))
            
            allocate(a_spr_v_tmp(nz_num_v_tmp))
            allocate(r_spr_v_tmp(nz_num_v_tmp))
            allocate(c_spr_v_tmp(nz_num_v_tmp))
            allocate(b_spr_v_tmp(N_eq_v_tmp))

            ! Allocate previous timestep fields if kim and moin scheme is turned on
            allocate(u_nm1(0:N_y+2,0:N_x+2))
            allocate(v_nm1(0:N_y+2,0:N_x+2))
            
            allocate(u_h_nm1(0:N_y+1,0:N_x+1))
            allocate(v_h_nm1(0:N_y+1,0:N_x+1))
            
            allocate(u_tmp_spr(N_eq_u_tmp))
            allocate(v_tmp_spr(N_eq_u_tmp))
           
            ! Create the sparse coeff matrix for u* and v* (used to solve K-M Ax=b)
            call create_spr_index_km(spr_km_u_tmp, spr_km_v_tmp, spr_rev_km_u_tmp, spr_rev_km_v_tmp)
            call sparse_vectors_km(spr_km_u_tmp, spr_km_v_tmp, r_spr_u_tmp, r_spr_v_tmp, & 
                  c_spr_u_tmp, c_spr_v_tmp, a_spr_u_tmp, a_spr_v_tmp)
      
      end if
      
      ! create initial conditions of the domain
      u = 0
      v = 0
      p = 0
      
      skip_to_k_m = .False.
      
      ! create the inlet velocity
      u(1:N_y,1) = u_inlet
    
      call interp_vel(u, v, u_h, v_h)

      ! mesh the grid
      call create_grid(x_mesh_1, y_mesh_1, x_h, y_h)

      ! begin time stepping
      write(6,*) "Begin Time Steps"
      t_step = 0
      t_sim = 0D0
      
      do while (t_sim.lt.0.5)
            
            ! update the timestep
            t_step = t_step + 1
            t_sim = t_sim + delta_t
            
            ! need to set ghost cells/maintain BC's
            call ghost_operator(u,v)

            if (kim_moin.eqv..True.) then
                
                  ! switch to kim and moin method after 3 timesteps if method is turned on
                  if(t_step.gt.t_step_trigger) then
                  
                        skip_to_k_m = .True.
                       
                        if(t_step.eq.t_step_trigger) then
                  
                              write(6,*) "scheme switch"
                  
                        end if
                        
                        call sparse_b_km(u, v, u_h, v_h, p, u_nm1, v_nm1, u_h_nm1, v_h_nm1, b_spr_u_tmp, b_spr_v_tmp)
                        
                        call mgmres_st(N_eq_u_tmp, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp, a_spr_u_tmp, u_tmp_spr, b_spr_u_tmp, & 
                              itr_max, mr_u_tmp, tol_abs, tol_rel, rho, itr_used) 
                  
                        call update_u_tmp(u_tmp_spr, spr_rev_km_u_tmp, u_tmp)
                        
                        call mgmres_st(N_eq_v_tmp, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp, a_spr_v_tmp, v_tmp_spr, b_spr_v_tmp, & 
                              itr_max, mr_v_tmp, tol_abs, tol_rel, rho, itr_used) 
                        call update_v_tmp(v_tmp_spr, spr_rev_km_v_tmp, v_tmp)

                  end if
            end if

            if(skip_to_k_m.eqv..False.) then

                  ! call the convective and diffusive operators
                  call convective_operator(u, v, N_u, N_v)
                  call diffusive_operator(u, v, D_u, D_v)

                  ! create the temporary velocities for interior nodes
                  call vel_tmp(u, v, N_u, N_v, D_u, D_v, u_tmp, v_tmp)
                  
                  ! enforce u* and v* for border cells only
                  call u_tmp_boundaries(u_tmp, p)
                  call v_tmp_boundaries(v_tmp, p)

            end if
            

            ! create b vector for Ax = b
            call sparse_b(u_tmp, v_tmp, b_spr)
           
            ! only save previous step if kim and moin function is true
            ! u and v will be changed in this timestep, so the future "previous step" has to be updated to the current timestep
            if(kim_moin.eqv..True.) then
                  u_nm1 = u
                  v_nm1 = v
                  u_h_nm1 = u_h
                  v_h_nm1 = v_h
            end if

            ! if verbose is turned on then see matrix
            call verbose_poisson_matrix(.False., r_spr, c_spr, a_spr, b_spr)

            ! solve for P_est from Ax = b using sparse formulation
            call mgmres_st(N_eq, nz_num, r_spr, c_spr, a_spr, P_est, b_spr, itr_max, mr, tol_abs, tol_rel, rho, itr_used)

            ! Update the pressures to the next timestep
            call update_p(P_est, spr_index_rev, P)

            ! Update the velocity to the next timestep
            call update_vel(P, u_tmp, v_tmp, u, v)
         
            ! Interpolate the velocities
            call interp_vel(u, v, u_h, v_h)           

            ! Output to Terminal
            write(6,*) "Timestep = ", t_step, " Time simulation = ", t_sim 
            write(6,*) "Residual = ", rho, " Iterations Used = ", itr_used
            
            if(mod(t_step,1).eq.0) then
                  call write_output_to_file(t_step, x_h, y_h, u_h, v_h, P)
            end if 
      
      end do

      deallocate(a_spr)      
      deallocate(r_spr)      
      deallocate(c_spr)      
      deallocate(x_mesh_1)
      deallocate(y_mesh_1)
      deallocate(x_h)
      deallocate(y_h)
      deallocate(u)
      deallocate(v)
      deallocate(p)
      deallocate(spr_index)
      deallocate(u_h)
      deallocate(v_h)
      deallocate(N_u)
      deallocate(N_v)
      deallocate(b_spr)
      deallocate(P_est)
      
      if(kim_moin.eqv..True.) then
            deallocate(spr_km_u_tmp) 
            deallocate(spr_km_v_tmp) 
            
            deallocate(spr_rev_km_u_tmp) 
            deallocate(spr_rev_km_v_tmp) 
            
            deallocate(a_spr_u_tmp)
            deallocate(r_spr_u_tmp)
            deallocate(c_spr_u_tmp)
            deallocate(b_spr_u_tmp)
            
            deallocate(a_spr_v_tmp)
            deallocate(r_spr_v_tmp)
            deallocate(c_spr_v_tmp)
            deallocate(b_spr_v_tmp)
            
            deallocate(u_nm1)
            deallocate(v_nm1)
            deallocate(u_h_nm1)
            deallocate(v_h_nm1)
            
            deallocate(u_tmp_spr)
            deallocate(v_tmp_spr)
      end if

end program channel_flow

subroutine write_output_to_file(t_step, x_h, y_h, u_h, v_h, P)
     
      use program_variables

      implicit none
      integer, intent(in) :: t_step
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(in) :: u_h, v_h
      double precision, dimension(0:N_y + 2, 0:N_x + 2), intent(in) :: P
      double precision, dimension(0:N_x+1), intent(in) :: x_h
      double precision, dimension(0:N_y+1), intent(in) :: y_h
     
      character(len = 50) :: output_file_name, file_status
      logical :: file_exist
      integer :: i, j

      write(output_file_name, 201) t_step
      inquire(file = output_file_name, exist = file_exist)
      if(file_exist.eqv..True.) then
        file_status = "replace"
      else
        file_status = "new"
      end if
      open(unit = 1, file = output_file_name, form = 'formatted', status = file_status, action = 'write')

      write(1,101)
      write(1,102) N_x, N_y

      do i = 1, N_y
            do j = 1, N_x
                  write(1,103) x_h(j), y_h(i), u_h(i,j), v_h(i,j), P(i,j)
            end do
      end do
     
      close(1)

      101 format('VARIABLES = "X", "Y", "U", "V", "P"')
      102 format('ZONE I=', I6, ', J=', I6, ', K=1', ', F= POINT')
      103 format(5(es14.7,1x))
      201 FORMAT("flow2DIncomp_", I10.10,".dat")

end subroutine write_output_to_file

subroutine interp_vel(u, v, u_h, v_h)

      use program_variables
      implicit none

      double precision, dimension(0:N_y+2, 0:N_x+2), intent(in) :: u, v
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(out) :: u_h, v_h
      integer :: i, j

      do i = 1, N_y
            do j = 1, N_x
                  u_h(i,j) = (u(i,j) + u(i,j+1))/2D0
            end do
      end do
      
      do i = 1, N_y
            do j = 1, N_x
                  v_h(i,j) = (v(i+1,j) + v(i,j))/2D0
            end do
      end do

end subroutine interp_vel

subroutine update_vel(P, u_tmp, v_tmp, u, v)

      use program_variables
      implicit none
      double precision, dimension(0:N_y + 2, 0:N_x + 2), intent(in) :: P
      double precision, dimension(0:N_y+2, 0:N_x+2), intent(in) :: u_tmp, v_tmp
      double precision, dimension(0:N_y+2, 0:N_x+2), intent(inout) :: u, v
      integer :: i, j
      
      do i = 1, N_y
            do j = 1, N_x+1
                  u(i,j) = u_tmp(i,j) - delta_t * (p(i,j) - p(i,j-1))/h_x
            end do
      end do
      do i = 1, N_y + 1
            do j = 1, N_x
                  v(i,j) = v_tmp(i,j) - delta_t * (p(i,j) - p(i-1,j))/h_y
            end do
      end do

end subroutine update_vel

subroutine update_p(P_est, spr_index_rev, P)

      use program_variables
      implicit none

      integer, dimension(N_eq,2), intent(in) :: spr_index_rev
      double precision, dimension(N_eq), intent(in) :: P_est
      double precision, dimension(0:N_y + 2, 0:N_x + 2), intent(out) :: P
      integer :: i

      do i = 1, N_eq
            P(spr_index_rev(i,1), spr_index_rev(i,2)) = P_est(i)
      end do

end subroutine update_p

subroutine update_u_tmp(u_tmp_spr, spr_index_rev, u_tmp)

      use program_variables
      integer, dimension(N_eq_u_tmp,2), intent(in) :: spr_index_rev
      double precision, dimension(N_eq_u_tmp), intent(in) :: u_tmp_spr
      double precision, dimension(0:N_y + 2, 0:N_x + 2), intent(out) :: u_tmp
      integer :: i

      do i = 1, N_eq_u_tmp
            u_tmp(spr_index_rev(i,1), spr_index_rev(i,2)) = u_tmp_spr(i)
      end do

end subroutine update_u_tmp

subroutine update_v_tmp(v_tmp_spr, spr_index_rev, v_tmp)

      use program_variables
      integer, dimension(N_eq_v_tmp,2), intent(in) :: spr_index_rev
      double precision, dimension(N_eq_v_tmp), intent(in) :: v_tmp_spr
      double precision, dimension(0:N_y + 2, 0:N_x + 2), intent(out) :: v_tmp
      integer :: i

      do i = 1, N_eq_v_tmp
            v_tmp(spr_index_rev(i,1), spr_index_rev(i,2)) = v_tmp_spr(i)
      end do

end subroutine update_v_tmp

subroutine verbose_poisson_matrix(on, r_spr, c_spr, a_spr, b_spr)

      use program_variables
      
      implicit none

      logical, intent(in) :: on
      double precision, dimension(nz_num), INTENT(IN):: a_spr
      integer, dimension(nz_num), INTENT(IN) :: r_spr, c_spr
      double precision, dimension(N_eq, N_eq) :: A_full
      double precision, dimension(N_eq) :: b_spr
      integer :: i

      if(on.eqv..True.) then
            A_full = 0
            do i = 1, nz_num
                  write(6,*) i, r_spr(i), c_spr(i), a_spr(i)
                  A_full(r_spr(i),c_spr(i)) = a_spr(i)
            end do
            do i = 1, N_eq
                  write(6,*) A_full(i,:), b_spr(i), i, count(A_full(i,:)/=0)
            end do
      
      end if

end subroutine verbose_poisson_matrix

subroutine create_sparse_index(sparse_index, sparse_rev)
      
      use program_variables
     
      integer, dimension(0:N_y + 1, 0:N_x + 1), INTENT(OUT) :: sparse_index
      integer, dimension(n_eq, 2), INTENT(OUT) :: sparse_rev
      integer :: i, j, k
      
      k = 0
      do i = 0, N_y + 1
            do j = 0, N_x + 1
                  k = k + 1
                  sparse_index(i,j) = k
                  sparse_rev(k, 1) = i
                  sparse_rev(k, 2) = j
            end do
      end do

end subroutine create_sparse_index

subroutine create_spr_index_km(spr_km_u_tmp, spr_km_v_tmp, spr_rev_km_u_tmp, spr_rev_km_v_tmp)
      
      use program_variables
      
      integer, dimension(0:N_y + 1, 1:N_x + 1), INTENT(OUT) :: spr_km_u_tmp
      integer, dimension(1:N_y + 1, 0:N_x + 1), INTENT(OUT) :: spr_km_v_tmp
      integer, dimension(N_eq_u_tmp, 2), INTENT(OUT) :: spr_rev_km_u_tmp
      integer, dimension(N_eq_v_tmp, 2), INTENT(OUT) :: spr_rev_km_v_tmp
      integer :: i, j, k

      k = 0
      do i = 0, N_y+1
            do j = 1, N_x + 1
                  k = k + 1
                  spr_km_u_tmp(i,j) = k
                  spr_rev_km_u_tmp(k,1) = i
                  spr_rev_km_u_tmp(k,2) = j
                  write(6,*) i,j,spr_km_u_tmp(i,j)
            end do
      end do
      
      k = 0
      do i = 1, N_y+1
            do j = 0, N_x + 1
                  k = k + 1
                  spr_km_v_tmp(i,j) = k
                  spr_rev_km_v_tmp(k,1) = i
                  spr_rev_km_v_tmp(k,2) = j
                  write(6,*) i,j,spr_km_v_tmp(i,j)
            end do
      end do

      write(6,*) "u*"
      do k = 1, N_eq_u_tmp
            write(6,*) k, spr_rev_km_u_tmp(k,1), spr_rev_km_u_tmp(k,2)
      end do
     
      write(6,*) "v*"
      do k = 1, N_eq_v_tmp
            write(6,*) k, spr_rev_km_v_tmp(k,1), spr_rev_km_v_tmp(k,2)
      end do


end subroutine create_spr_index_km

subroutine sparse_vectors(spr_index, r_spr, c_spr, a_spr)
     
      use program_variables
      
      integer :: i, j, row_cnt, nz_index
      integer, intent(in), dimension(0:N_y + 1, 0:N_x + 1) :: spr_index
      integer, dimension(nz_num), intent(out) :: r_spr, c_spr
      double precision, dimension(nz_num), intent(out) :: a_spr

      row_cnt = 0
      nz_index = 0
      do i = 0, N_y + 1
            do j = 0, N_x + 1
           
                  row_cnt = row_cnt + 1

                  ! to access the corner cells set 1
                  if (i.eq.0.and.j.eq.0.or.i.eq.N_y+1.and.j.eq.N_x+1) then
                        ! set b in Ax = b to 0
                        !write(6,*) spr_index(i,j), "Corner"
                        call sparse_build(row_cnt, spr_index(i,j), nz_index, r_spr, c_spr)
                        a_spr(nz_index) = 1
                        cycle
                  end if

                  ! to access the corner cells set 2
                  if(i.eq.0.and.j.eq.N_x+1.or.i.eq.N_y+1.and.j.eq.0) then
                        ! set b in Ax = b to 0
                        !write(6,*) spr_index(i,j), "Corner"
                        call sparse_build(row_cnt, spr_index(i,j), nz_index, r_spr, c_spr) 
                        a_spr(nz_index) = 1
                        cycle
                  end if
           
                  !write(6,*) spr_index(i,j)
           
                  ! to access the bottom P ghost cells (grad P = 0)
                  if(j.ge.1.and.j.le.N_x.and.i.eq.0) then
                       !write(6,*) spr_index(i,j), spr_index(i+1,j), "Lower B"
                       call sparse_build(row_cnt, spr_index(i,j), nz_index, r_spr, c_spr)
                       !a_spr(nz_index) = -d_2 
                       a_spr(nz_index) = 1 
                       call sparse_build(row_cnt, spr_index(i+1,j), nz_index, r_spr, c_spr) 
                       a_spr(nz_index) = -1
                       cycle 
                  end if
                  
                  ! to access the top P ghost cells (grad P = 0)
                  if(j.ge.1.and.j.le.N_x.and.i.eq.N_y+1) then
                       !write(6,*) spr_index(i,j), spr_index(i-1,j), "Upper B"
                       call sparse_build(row_cnt, spr_index(i-1,j), nz_index, r_spr, c_spr)
                       !a_spr(nz_index) = -d_2 
                       a_spr(nz_index) = -1 
                       call sparse_build(row_cnt, spr_index(i,j), nz_index, r_spr, c_spr) 
                       a_spr(nz_index) = 1
                       cycle 
                  end if

                  ! to access the inlet P ghost cells (grad P = 0)
                  if(j.eq.0.and.i.ge.1.and.i.le.N_y) then
                        !write(6,*) spr_index(i,j), spr_index(i,j+1), "Left B"
                        call sparse_build(row_cnt, spr_index(i,j), nz_index, r_spr, c_spr)
                        !a_spr(nz_index) = -d_1 
                        a_spr(nz_index) = -1 
                        call sparse_build(row_cnt, spr_index(i,j+1), nz_index, r_spr, c_spr) 
                        !a_spr(nz_index) = d_1 
                        a_spr(nz_index) = 1 
                        cycle      
                  end if
                  
                  ! to access the outlet P ghost cells (P = 0)
                  if(j.eq.N_x+1.and.i.ge.1.and.i.le.N_y) then
                        !write(6,*) spr_index(i,j), spr_index(i,j-1), "Right B"
                        !call sparse_build(row_cnt, spr_index(i,j-3), nz_index, r_spr, c_spr) 
                        !a_spr(nz_index) = -d_1 
                        !call sparse_build(row_cnt, spr_index(i,j-2), nz_index, r_spr, c_spr) 
                        !a_spr(nz_index) = 5 * d_1 
                        !call sparse_build(row_cnt, spr_index(i,j-1), nz_index, r_spr, c_spr) 
                        !a_spr(nz_index) = -7 * d_1 
                        !call sparse_build(row_cnt, spr_index(i,j), nz_index, r_spr, c_spr) 
                        !a_spr(nz_index) = 3 * d_1 
                        
                        call sparse_build(row_cnt, spr_index(i,j-1), nz_index, r_spr, c_spr) 
                        a_spr(nz_index) = 1 
                        call sparse_build(row_cnt, spr_index(i,j), nz_index, r_spr, c_spr) 
                        a_spr(nz_index) = 1
                        cycle      
                  end if

                  !write(6,*) spr_index(i,j), spr_index(i-1,j), spr_index(i+1,j), spr_index(i,j-1), spr_index(i,j+1), "interior"
                  call sparse_build(row_cnt, spr_index(i-1,j), nz_index, r_spr, c_spr)
                  a_spr(nz_index) = c_2 
                  call sparse_build(row_cnt, spr_index(i,j-1), nz_index, r_spr, c_spr) 
                  a_spr(nz_index) = c_1
                  call sparse_build(row_cnt, spr_index(i,j), nz_index, r_spr, c_spr)
                  a_spr(nz_index) = -c_3
                  call sparse_build(row_cnt, spr_index(i,j+1), nz_index, r_spr, c_spr)
                  a_spr(nz_index) = c_1
                  call sparse_build(row_cnt, spr_index(i+1,j), nz_index, r_spr, c_spr)
                  a_spr(nz_index) = c_2


            end do
      end do

      !do i = 1, N_eq
      !      write(6,*) i, spr_index_rev(i, :)
      !end do

end subroutine sparse_vectors

subroutine sparse_vectors_km(spr_km_u_tmp, spr_km_v_tmp, r_spr_u_tmp, r_spr_v_tmp, & 
      c_spr_u_tmp, c_spr_v_tmp, a_spr_u_tmp, a_spr_v_tmp)

      use program_variables
      
      implicit none
      integer, dimension(0:n_y + 1, 1:n_x + 1), intent(in) :: spr_km_u_tmp
      integer, dimension(1:n_y + 1, 0:n_x + 1), intent(in) :: spr_km_v_tmp
      integer, dimension(nz_num_u_tmp), intent(out) :: r_spr_u_tmp, c_spr_u_tmp
      integer, dimension(nz_num_v_tmp), intent(out) :: r_spr_v_tmp, c_spr_v_tmp
      double precision, dimension(nz_num_u_tmp), intent(out) :: a_spr_u_tmp
      double precision, dimension(nz_num_v_tmp), intent(out) :: a_spr_v_tmp
      
      integer :: i, j, row_cnt, nz_index
      
      ! build sparse vectors for u*
      row_cnt = 0
      nz_index = 0
      
      do i = 0, N_y+1
            do j = 1, N_x + 1

                  row_cnt = row_cnt + 1

                  ! corner on one side
                  if (i.eq.0.and.j.eq.1.or.i.eq.0.and.j.eq.N_x+1) then
                        
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_u_tmp(i,j), "corner u*"
                        cycle
                  end if
                  
                  ! corner on other side
                  if (i.eq.N_y+1.and.j.eq.1.or.i.eq.N_y+1.and.j.eq.N_x+1) then
                        
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_u_tmp(i,j), "corner u*"
                        cycle
                  end if

                  ! lower wall
                  if(j.ge.2.and.j.le.N_x.and.i.eq.0) then
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = 1
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i+1,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_u_tmp(i,j), "lower wall u*"
                        cycle 
                   end if

                  ! upper wall
                  if(j.ge.2.and.j.le.N_x.and.i.eq.N_y+1) then
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i-1,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = 1
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_u_tmp(i,j), "upper wall u*"
                       cycle 
                  end if


                  ! inlet
                  if(i.ge.1.and.i.le.N_y.and.j.eq.1) then
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_u_tmp(i,j), "inlet u*"
                        cycle
                  end if

                  ! outlet
                  if(i.ge.1.and.i.le.N_y.and.j.eq.N_x+1) then
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j-2), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = 1
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j-1), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = -4
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = 3
                        write(6,*) i, j, spr_km_u_tmp(i,j), "outlet u*"
                        cycle
                  end if

                  ! interior
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i-1,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = km * c_2
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j-1), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = km * c_1
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = (1 - km * c_3)
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i,j+1), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = km * c_1
                        call sparse_build_km(row_cnt, spr_km_u_tmp(i+1,j), nz_index, nz_num_u_tmp, r_spr_u_tmp, c_spr_u_tmp)
                        a_spr_u_tmp(nz_index) = km * c_2
                  write(6,*) i, j, spr_km_u_tmp(i,j), "interior u*"

            end do
      end do


      write(6,*) " "
      ! build sparse vectors for v*

      row_cnt = 0
      nz_index = 0
      do i = 1, N_y+1
            do j = 0, N_x + 1
                  row_cnt = row_cnt + 1
                  ! corner
                  if (i.eq.1.and.j.eq.0.or.i.eq.1.and.j.eq.N_x+1) then
                        call sparse_build_km(row_cnt, spr_km_v_tmp(i,j), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                        a_spr_v_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_v_tmp(i,j), "corner v*"
                        cycle
                  end if
                  
                  ! corner on other side
                  if (i.eq.N_y+1.and.j.eq.0.or.i.eq.N_y+1.and.j.eq.N_x+1) then
                        call sparse_build_km(row_cnt, spr_km_v_tmp(i,j), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                        a_spr_v_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_v_tmp(i,j), "corner v*"
                        cycle
                  end if

                  ! lower wall
                  if(j.ge.1.and.j.le.N_x.and.i.eq.1) then
                        call sparse_build_km(row_cnt, spr_km_v_tmp(i,j), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                        a_spr_v_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_v_tmp(i,j), "lower wall v*"
                        cycle 
                   end if

                  ! upper wall
                  if(j.ge.1.and.j.le.N_x.and.i.eq.N_y+1) then
                        call sparse_build_km(row_cnt, spr_km_v_tmp(i,j), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                        a_spr_v_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_v_tmp(i,j), "upper wall v*"
                       cycle 
                  end if


                  ! inlet
                  if(i.ge.2.and.i.le.N_y.and.j.eq.0) then
                        call sparse_build_km(row_cnt, spr_km_v_tmp(i,j), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                        a_spr_v_tmp(nz_index) = 1
                        call sparse_build_km(row_cnt, spr_km_v_tmp(i,j+1), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                        a_spr_v_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_v_tmp(i,j), "inlet v*"
                        cycle
                  end if

                  ! outlet
                  if(i.ge.2.and.i.le.N_y.and.j.eq.N_x+1) then
                        call sparse_build_km(row_cnt, spr_km_v_tmp(i,j-1), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                        a_spr_v_tmp(nz_index) = -1
                        call sparse_build_km(row_cnt, spr_km_v_tmp(i,j), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                        a_spr_v_tmp(nz_index) = 1
                        write(6,*) i, j, spr_km_v_tmp(i,j), "outlet v*"
                        cycle
                  end if

                  ! interior
                  call sparse_build_km(row_cnt, spr_km_v_tmp(i-1,j), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                  a_spr_v_tmp(nz_index) = km * c_2
                  call sparse_build_km(row_cnt, spr_km_v_tmp(i,j-1), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                  a_spr_v_tmp(nz_index) = km * c_1
                  call sparse_build_km(row_cnt, spr_km_v_tmp(i,j), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                  a_spr_v_tmp(nz_index) = (1 - km * c_3)
                  call sparse_build_km(row_cnt, spr_km_v_tmp(i,j+1), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                  a_spr_v_tmp(nz_index) = km * c_1
                  call sparse_build_km(row_cnt, spr_km_v_tmp(i+1,j), nz_index, nz_num_v_tmp, r_spr_v_tmp, c_spr_v_tmp)
                  a_spr_v_tmp(nz_index) = km * c_2
                  write(6,*) i, j, spr_km_v_tmp(i,j), "interior v*"

            end do
      end do

end subroutine sparse_vectors_km

subroutine sparse_b(u_tmp, v_tmp, b_spr)
      use program_variables
     
      double precision, dimension(0:N_y+2, 0:N_x+2), intent(out) :: u_tmp, v_tmp
      double precision, dimension(N_eq), intent(out) :: b_spr
      integer :: i, j, row_cnt, nz_index
     
      row_cnt = 0
      nz_index = 0
      do i = 0, N_y + 1
            do j = 0, N_x + 1
           
                  row_cnt = row_cnt + 1
                  if (i.eq.0.and.j.eq.0.or.i.eq.N_y+1.and.j.eq.N_x+1) then
                        b_spr(row_cnt) = 0
                        cycle
                  end if

                  if(i.eq.0.and.j.eq.N_x+1.or.i.eq.N_y+1.and.j.eq.0) then
                        b_spr(row_cnt) = 0
                        cycle
                  end if
           
                  if(j.ge.1.and.j.le.N_x.and.i.eq.0) then
                        b_spr(row_cnt) = 0
                        
                       cycle 
                  end if
                  
                  if(j.ge.1.and.j.le.N_x.and.i.eq.N_y+1) then
                        b_spr(row_cnt) = 0
                       cycle 
                  end if

                  if(j.eq.0.and.i.ge.1.and.i.le.N_y) then
                        b_spr(row_cnt) = 0
                        cycle      
                  end if
                  
                  if(j.eq.N_x+1.and.i.ge.1.and.i.le.N_y) then
                        b_spr(row_cnt) = 0
                        cycle      
                  end if

                  b_spr(row_cnt) = (1D0/delta_t) * (((u_tmp(i,j+1) - u_tmp(i,j))/h_x) + ((v_tmp(i+1,j) - v_tmp(i,j))/h_y))


            end do
      end do


end subroutine sparse_b

subroutine sparse_b_km(u, v, u_h, v_h, p, u_nm1, v_nm1, u_h_nm1, v_h_nm1, b_spr_u_tmp, b_spr_v_tmp)
      use program_variables
      
      integer :: i, j, row_cnt, nz_index
      double precision :: G1, G2, G3, Nu_n, Nu_nm1, Nv_n, Nv_nm1, Lu, Lv
      double precision, dimension(0:N_y+2, 0:N_x+2), intent(in) :: u, v, p, u_nm1, v_nm1  
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(in) :: u_h, v_h, u_h_nm1, v_h_nm1  
     
      double precision, dimension(nz_num_u_tmp), intent(out) :: b_spr_u_tmp
      double precision, dimension(nz_num_v_tmp), intent(out) :: b_spr_v_tmp
      row_cnt = 0
      nz_index = 0
      
      do i = 0, N_y+1
            do j = 1, N_x + 1

                  row_cnt = row_cnt + 1

                  ! corner on one side
                  if (i.eq.0.and.j.eq.1.or.i.eq.0.and.j.eq.N_x+1) then
                        b_spr_u_tmp(row_cnt) = 0
                        cycle
                  end if
                  
                  ! corner on other side
                  if (i.eq.N_y+1.and.j.eq.1.or.i.eq.N_y+1.and.j.eq.N_x+1) then
                        b_spr_u_tmp(row_cnt) = 0
                        cycle
                  end if

                  ! lower wall
                  if(j.ge.2.and.j.le.N_x.and.i.eq.0) then
                        call g_operator(p, i, j, 1, G1)
                        call g_operator(p, i+1, j, 1, G2)
                        b_spr_u_tmp(row_cnt) = delta_t * (G1 + G2)
                        cycle 
                   end if

                  ! upper wall
                  if(j.ge.2.and.j.le.N_x.and.i.eq.N_y+1) then
                        call g_operator(p, i-1, j, 1, G2)
                        call g_operator(p, i, j, 1, G2)
                        b_spr_u_tmp(row_cnt) = delta_t * (G1 + G2)
                        cycle 
                  end if


                  ! inlet
                  if(i.ge.1.and.i.le.N_y.and.j.eq.1) then
                        b_spr_u_tmp(row_cnt) = u_inlet
                        cycle
                  end if

                  ! outlet
                  if(i.ge.1.and.i.le.N_y.and.j.eq.N_x+1) then
                        call g_operator(p, i, j, 1, G1)
                        call g_operator(p, i, j-1, 1, G2)
                        call g_operator(p, i, j-2, 1, G3)
                        b_spr_u_tmp(row_cnt) = delta_t * (3 * G1 - 4 * G2 + G3) 
                        cycle
                  end if

                  ! interior
                  call convective_operator_km(i, j, u, v, u_h, v_h, 1, Nu_n)
                  call convective_operator_km(i, j, u_nm1, v_nm1, u_h_nm1, v_h_nm1, 1, Nu_nm1)
                  call diffusive_operator_km(i, j, u, v, 1, Lu)
                  b_spr_u_tmp(row_cnt) = (delta_t/2)*(3*Nu_n - Nu_nm1) + u(i,j) - km * Lu

            end do
      end do


      row_cnt = 0
      nz_index = 0
      do i = 1, N_y+1
            do j = 0, N_x + 1
                  row_cnt = row_cnt + 1
                  ! corner
                  if (i.eq.1.and.j.eq.0.or.i.eq.1.and.j.eq.N_x+1) then
                        b_spr_v_tmp(row_cnt) = 0
                        cycle
                  end if
                  
                  ! corner on other side
                  if (i.eq.N_y+1.and.j.eq.0.or.i.eq.N_y+1.and.j.eq.N_x+1) then
                        b_spr_v_tmp(row_cnt) = 0
                        cycle
                  end if

                  ! lower wall
                  if(j.ge.1.and.j.le.N_x.and.i.eq.1) then
                        b_spr_v_tmp(row_cnt) = 0
                        cycle 
                   end if

                  ! upper wall
                  if(j.ge.1.and.j.le.N_x.and.i.eq.N_y+1) then
                        b_spr_v_tmp(row_cnt) = 0
                        cycle 
                  end if


                  ! inlet
                  if(i.ge.2.and.i.le.N_y.and.j.eq.0) then
                        call g_operator(p, i, j, 2, G1)
                        call g_operator(p, i, j+1, 2, G2)
                        b_spr_v_tmp(row_cnt) = 2 * v_inlet + delta_t * (G1 + G2)
                        cycle
                  end if

                  ! outlet
                  if(i.ge.2.and.i.le.N_y.and.j.eq.N_x+1) then
                        call g_operator(p, i, j-1, 2, G1)
                        call g_operator(p, i, j, 2, G2)
                        b_spr_v_tmp(row_cnt) = delta_t * (G2 - G1)
                        cycle
                  end if

                  ! interior
                  call convective_operator_km(i, j, u, v, u_h, v_h, 2, Nv_n)
                  call convective_operator_km(i, j, u_nm1, v_nm1, u_h_nm1, v_h_nm1, 2, Nv_nm1)
                  call diffusive_operator_km(i, j, u, v, 2, Lv)
                  b_spr_v_tmp(row_cnt) = (delta_t/2)*(3*Nv_n - Nv_nm1) + v(i,j) - km * Lv

            end do
      end do
      
end subroutine sparse_b_km

subroutine sparse_build(row_cnt, col_cnt, nz_index, r_spr, c_spr)

      use program_variables
      integer, intent(in)  :: row_cnt, col_cnt
      integer, intent(inout) :: nz_index
      integer, dimension(nz_num), intent(inout) :: r_spr, c_spr

      nz_index = nz_index + 1
      r_spr(nz_index) = row_cnt      
      c_spr(nz_index) = col_cnt

end subroutine sparse_build

subroutine sparse_build_km(row_cnt, col_cnt, nz_index, nz_num_uv, r_spr, c_spr)
      use program_variables
      integer, intent(in)  :: row_cnt, col_cnt, nz_num_uv
      integer, intent(inout) :: nz_index
      integer, dimension(nz_num_uv), intent(inout) :: r_spr, c_spr

      nz_index = nz_index + 1
      r_spr(nz_index) = row_cnt      
      c_spr(nz_index) = col_cnt
end subroutine sparse_build_km

subroutine create_grid(x_mesh_1, y_mesh_1, x_h, y_h)

          use program_variables

          ! length of a grid in x and y

          double precision, dimension(0:N_x+2), intent(out) :: x_mesh_1
          double precision, dimension(0:N_y+2), intent(out) :: y_mesh_1
          
          double precision, dimension(0:N_x+1), intent(out) :: x_h
          double precision, dimension(0:N_y+1), intent(out) :: y_h

          integer :: i, j


          do i = 0, N_x + 2
                x_mesh_1(i) = (i-1) * h_x      
          end do

          do j = 0, N_y + 2      
                y_mesh_1(j) = (j-1) * h_y      
          end do

          do i = 0, N_x + 1
                x_h(i) = (x_mesh_1(i) + x_mesh_1(i+1))/2
          end do
          
          do j = 0, N_y + 1
                y_h(j) = (y_mesh_1(j) + y_mesh_1(j+1))/2
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

subroutine convective_operator(u, v, N_u, N_v)

      use program_variables

      implicit none
      double precision, dimension(0:N_y+2, 0:N_x+2), intent(in) :: u, v
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(out) :: N_u, N_v

      integer :: i,j
      double precision :: u_e, u_w, u_n, u_s
      double precision :: v_e, v_w, v_n, v_s

      ! u nodes from Row 1 to Row N_y, Column 2 to N_x
      do i = 1, N_y
            do j = 2, N_x
                  u_e = 0.5 * (u(i,j+1) + u(i,j))
                  u_w = 0.5 * (u(i,j) + u(i,j-1))
                  u_n = 0.5 * (u(i+1,j) + u(i,j))
                  u_s = 0.5 * (u(i,j) + u(i-1,j))
                  
                  v_n = 0.5 * (v(i+1,j-1) + v(i+1,j))
                  v_s = 0.5 * (v(i,j-1) + v(i,j))
           
                  N_u(i,j) = (u_e * u_e - u_w * u_w) * h_y + (v_n * u_n - v_s * u_s) * h_x
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
                  
                  N_v(i,j) = (u_e * v_e - u_w * v_w) * h_y + (v_n * v_n - v_s * v_s) * h_x

            end do
      end do

end subroutine convective_operator

subroutine diffusive_operator(u, v, D_u, D_v)

      use program_variables

      implicit none

      double precision, dimension(0:N_y+2, 0:N_x+2), intent(in) :: u, v
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(out) :: D_u, D_v
     
      integer :: i, j

      do i = 1, N_y
            do j = 2, N_x
                  D_u(i,j) = ((u(i,j+1) - 2 * u(i,j) + u(i,j-1))*h_y/h_x) + ((u(i+1,j) - 2*u(i,j) + u(i-1,j)) * h_x/h_y)
            end do
      end do

      do i = 2, N_y
            do j = 1, N_x
                  D_v(i,j) = ((v(i,j+1) - 2 * v(i,j) + v(i,j-1))*h_y/h_x) + ((v(i+1,j) - 2*v(i,j) + v(i-1,j)) * h_x/h_y)
            end do
      end do

end subroutine diffusive_operator

subroutine ghost_operator(u, v)

      use program_variables

      implicit none
      double precision, dimension(0:N_y+2, 0:N_x+2), intent(inout) :: u, v
      integer :: i, j
      double precision :: u_wall

      u_wall = 0

      do j = 2, N_x
            
            ! inlet boundary conditions
            u(0,j) = 2 * u_wall - u(1,j)

            ! outlet boundary conditions
            u(N_y+1,j) = 2 * u_wall - u(N_y,j)

      end do

      do i = 2, N_y

            ! inlet boundary conditions (dirichlet v = 0)
            v(i,0) = 2*v_inlet - v(i,1)

            ! outlet boundary conditions (neumann dv/dx = 0)
            v(i,N_x + 1) = v(i,N_x)

      end do



end subroutine ghost_operator

subroutine vel_tmp(u, v, N_u, N_v, D_u, D_v, u_tmp, v_tmp)

      use program_variables
      
      implicit none

      double precision, dimension(0:N_y+2, 0:N_x+2), intent(in) :: u, v
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(in) :: N_u, N_v
      double precision, dimension(0:N_y+1, 0:N_x+1), intent(in) :: D_u, D_v
      double precision, dimension(0:N_y+2, 0:N_x+2), intent(out) :: u_tmp, v_tmp

      integer :: i, j
      
      do i = 1, N_y
            do j = 1, N_x
                 !write(6,*) i,j 
                 u_tmp(i,j) = u(i,j) + (delta_t/(h_x * h_y)) * (-N_u(i,j) + (1D0/Re) * D_u(i,j))
            end do
      end do

      do i = 1, N_y
            do j = 1, N_x
                 !write(6,*) i,j 
                 v_tmp(i,j) = v(i,j) + (delta_t/(h_x * h_y)) * (-N_v(i,j) + (1D0/Re) * D_v(i,j))
            end do
      end do

end subroutine vel_tmp

subroutine u_tmp_boundaries(u_tmp, p)
     
      use program_variables
    
      implicit none

      double precision, dimension(0:N_y+2,0:N_x+2), intent(inout) :: u_tmp
      double precision, dimension(0:N_y+2,0:N_x+2), intent(in) :: p

      integer :: i, j
      double precision :: G_Nxp1, G_Nx, G_Nxm1
      double precision :: G_0, G_1, G_Ny, G_Nyp1

      ! For u* (Row 0 (Bottom Wall) & Ny+1 (Top Wall) from Col 2 to Nx, Col 1 (inlet) & Nx+1 (outlet) from Row 1 to Ny)
      
      ! u* for the inlet and the outlet 
      do i = 1, N_y
      
        ! create u* for inlet boundaries (Dirichlet u = u_inlet)
        u_tmp(i,1) = u_inlet

        ! create u* for outlet boundaries (Neumann du/dx = 0)
        call g_operator(p, i, N_x+1, 1, G_Nxp1)
        call g_operator(p, i, N_x, 1, G_Nx)
        call g_operator(p, i, N_x-1, 1, G_Nxm1)
        u_tmp(i,N_x+1) = delta_t * G_Nxp1 + & 
            (4D0/3D0) * (u_tmp(i,N_x) - delta_t * G_Nx) - & 
            (1D0/3D0) * (u_tmp(i,N_x-1) - delta_t * G_Nxm1)

      end do

      ! u* for the top and bottom walls
      do j = 2, N_x

        ! create u* for the bottom wall
        call g_operator(p, 0, j, 1, G_0)
        call g_operator(p, 1, j, 1, G_1)
        u_tmp(0,j) = -u_tmp(1,j) + delta_t * (G_0 + G_1) 

        ! create u* for the top wall
        call g_operator(p, N_y, j, 1, G_Ny)
        call g_operator(p, N_y+1, j, 1, G_Nyp1)
        u_tmp(N_y+1,j) = -u_tmp(N_y,j) + delta_t * (G_Ny + G_Nyp1) 

      end do
      
end subroutine u_tmp_boundaries

subroutine v_tmp_boundaries(v_tmp, p)

      use program_variables
      implicit none

      double precision, dimension(0:N_y+2,0:N_x+2), intent(inout) :: v_tmp
      double precision, dimension(0:N_y+2,0:N_x+2), intent(in) :: p
      integer :: i, j
      double precision :: G_0, G_1, G_Nx, G_Nxp1

      ! For v* (Row 1 (Bottom Wall) & Ny+1 (Top Wall) from Col 1 to Nx, Col 0 (inlet) & Nx+1 (outlet) from Row 2 to Ny)
      do i = 2, N_y 

            ! v* for the inlet
            call g_operator(p, i, 0, 2, G_0)
            call g_operator(p, i, 1, 2, G_1)
            v_tmp(i,0) = 2 * v_inlet - v_tmp(i,1) + delta_t * (G_0 + G_1)


            ! v* for the outlet
            call g_operator(p, i, N_x, 2, G_Nx)
            call g_operator(p, i, N_x+1, 2, G_Nxp1)
            v_tmp(i,N_x+1) = v_tmp(i,N_x) - delta_t * (G_Nx - G_Nxp1)

      end do


      do j = 1, N_x
           
            v_tmp(1,j) = 0
            v_tmp(N_y+1,j) = 0

      end do

end subroutine v_tmp_boundaries

subroutine g_operator(p, i, j, e_x, G)
   
    use program_variables

    implicit none

    double precision, dimension(0:N_y+2,0:N_x+2), intent(in) :: p
    integer, intent(in) :: i, j, e_x ! if e_x = 1 (d/dx), else (d/dy)
    double precision, intent(out) :: G


    ! x-axis for G
    if(e_x.eq.1) then
           g = (p(i,j) - p(i,j-1))/h_x
    else
           g = (p(i,j) - p(i-1,j))/h_y
    end if
     
end subroutine g_operator

subroutine convective_operator_km(i, j, u, v, u_h, v_h, ex, N)

      use program_variables

      implicit none
      
      integer, intent(in) :: i, j, ex
      double precision, dimension(0:N_x+2, 0:N_y+2), intent(in) :: u, v
      double precision, dimension(0:N_x+1, 0:N_y+1), intent(in) :: u_h, v_h
      double precision, intent(out) :: N


      if (ex.eq.1) then
            N = ((-u(i,j) * (u(i+1,j) - u(i-1,j)))/(2*h_x)) - (((v_h(i,j) + v_h(i-1,j))/2) * ((u(i,j+1) - u(i,j-1))/(2*h_y))) 
      else
            N = ((-v(i,j) * (v(i,j+1) - v(i,j-1)))/(2*h_y)) - (((u_h(i,j) + u_h(i,j-1))/2) * ((v(i+1,j) - v(i-1,j))/(2*h_x))) 
      end if

end subroutine convective_operator_km

subroutine diffusive_operator_km(i, j, u, v, ex, L)

      use program_variables

      integer, intent(in) :: i, j, ex
      double precision, dimension(0:N_x+2, 0:N_y+2), intent(in) :: u, v
      double precision, intent(out) :: L
      if (ex.eq.1) then
            L =  ((u(i,j+1) - 2*u(i,j) + u(i,j-1))/(h_x**2)) + ((u(i+1,j) - 2*u(i,j) + u(i-1,j))/(h_y**2))
      else
            L =  ((v(i,j+1) - 2*v(i,j) + v(i,j-1))/(h_x**2)) + ((v(i+1,j) - 2*v(i,j) + v(i-1,j))/(h_y**2))
      end if


end subroutine diffusive_operator_km