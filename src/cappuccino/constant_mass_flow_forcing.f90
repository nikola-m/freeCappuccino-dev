subroutine constant_mass_flow_forcing
! 
! Purpose:
!  Correct driving force for a constant mass flow rate.
! 
  use types
  use parameters, only: magUbar, gradPcmf
  use variables, only: U
  use sparse_matrix, only: apu
  use fieldManipulation, only: volumeWeightedAverage

  implicit none

  real(dp):: magUbarStar, rUAw, gragPplus, flowDirection

   ! # Extract the velocity in the flow direction
   ! magUbarStar = ( flowDirection & U ).weightedAverage( mesh.V() )
   magUbarStar = volumeWeightedAverage(U)

   ! # Calculate the pressure gradient increment needed to
   ! # adjust the average flow-rate to the correct value
   ! gragPplus = ( magUbar - magUbarStar ) / rUA.weightedAverage( mesh.V() )
   rUAw = volumeWeightedAverage(APU)
   gragPplus = ( magUbar - magUbarStar ) / rUAw

   ! #  Correction
   ! U.ext_assign( U + flowDirection * rUA * gragPplus )
   flowDirection = 1.0_dp
   U =  U + flowDirection * APU * gragPplus

   ! # Pressure gradient force that will drive the flow - we use it in calcuvw.
   gradPcmf  = gradPcmf + gragPplus

   write(6,'(2(a,es13.6))') "Uncorrected Ubar = ",magUbarStar," pressure gradient = ",gradPcmf

end subroutine