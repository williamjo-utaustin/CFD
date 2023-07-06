!-----------------------------------------------------------------------
! Square Cylinder PDE Solver
!-----------------------------------------------------------------------
! This program computes square cylinder flow in 2D.
! All subroutines and all variables are contained in this code aside from GMRES.f90
!
! By default the equations for x-velocity (U), y-velocity (V) and
! pressure correction (p') are solved. 
! To compile the program type "gfortran Square_Cylinder_Flow.f90 gmres.f90 -o sq_cyl"
! Output files: - 
! To run type: "./sq_cyl"
! Use Tecplot or Paraview for Visualization/Post-Processing
!-----------------------------------------------------------------------
! Default Values for Low Order Method for Square Cylinder Flow
!-----------------------------------------------------------------------
! Default Channel Flow Variables
! L_x = 30D0
! L_y = 8D0
! delta_t = 0.000001
! Re = 1D0
!
! Default Values for Low Order:
! N_x = 150, N_y = 40, kim_moin = .False. 
!
! Default Values for KM:
! None
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
! Disclaimers and Notes
!-----------------------------------------------------------------------
! i is the index for rows, j is the index for columns
! Any array would read Array(Rows, Columns)
! Nx would be a count for columns 
! Ny would be a count for rows
!-----------------------------------------------------------------------











