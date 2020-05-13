subroutine constant_mass_flow_forcing
! 
! Purpose:
!  Correct driving force for a constant mass flow rate.
! 
  use types
  use parameters, only: magUbar, gradPcmf
  use variables, only: U
  use sparse_matrix, only: apu
  use geometry, only: numCells,Vol
  use fieldManipulation, only: volumeWeightedAverage

  implicit none

  ! integer :: inp
  real(dp):: magUbarStar, rUAw, gragPplus

   ! Extract the velocity in the flow direction
   magUbarStar = volumeWeightedAverage(U)

   ! Calculate the pressure gradient increment needed to
   ! adjust the average flow-rate to the correct value
   ! gragPplus = ( magUbar - magUbarStar ) / rUA.weightedAverage( mesh.V() )
   rUAw = volumeWeightedAverage(APU)
   gragPplus = ( magUbar - magUbarStar ) / rUAw

   ! Correction of velocity to satisfy mass flow
   U(1:numCells) =  U(1:numCells) + APU(1:numCells) * gragPplus * Vol(1:numCells)

   ! Pressure gradient force that will drive the flow - we use it in calcuvw.
   gradPcmf  = gradPcmf + gragPplus

   write(6,'(2(a,es13.6))') "  Uncorrected Ubar = ",magUbarStar," pressure gradient = ",gradPcmf

end subroutine