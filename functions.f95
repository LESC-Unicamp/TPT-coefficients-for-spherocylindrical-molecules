! ############################################################################################### !
!              Canonical Monte Carlo algorithm for ellipsoid-of-revolution molecules              !
!        This code contains all functions used in the main program and in the subroutines.        !
!                                                                                                 !
! Version number: 1.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Original Developer: Joyce Tavares Lopes                             !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                            June 20th                                            !
! ############################################################################################### !
! Main Reference:                     B. J. Berne, P. Pechukas                                    !
!                                  J. Chem. Phys. 56, 4213 (1972)                                 !
!                                      DOI: 10.1063/1.1677837                                     !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
! This function takes the relative orientations of two molecular ellipsoids of revolution i and j !
!    and the unit vector joining their centers of mass and calculates their contact distance.     !
!          See Berne and Pechukas, J. Chem. Phys. 56, 4213 (1972), for more information.          !
! *********************************************************************************************** !
double precision function sigmahgo(ei,ej,urij)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8, dimension (3)	:: ei	! Orientation of particle i
	real*8, dimension (3)	:: ej	! Orientation of particle j
	real*8, dimension (3)	:: urij	! Normalized vector distance between the centers of mass of particles i and j
	real*8			:: rei	! Dot product of normalized vector distance and orientation of particle i
	real*8			:: rej	! Dot product of normalized vector distance and orientation of particle j
	real*8			:: eiej	! Dot product of both orientations (particles i and j)

	! ***************************************************************************************
	! Dot product of the normalized vector distance and orientation of particle i
	! ***************************************************************************************
	rei  = ( urij(1) * ei(1) ) + ( urij(2) * ei(2) ) + ( urij(3) * ei(3) )
	! ***************************************************************************************
	! Dot product of the normalized vector distance and orientation of particle j
	! ***************************************************************************************
	rej  = ( urij(1) * ej(1) ) + ( urij(2) * ej(2) ) + ( urij(3) * ej(3) )
	! ***************************************************************************************
	! Dot product of both orientations (particles i and j)
	! ***************************************************************************************
	eiej = ( ei(1) * ej(1) ) + ( ei(2) * ej(2) ) + ( ei(3) * ej(3) )

	! ***************************************************************************************
	! HGO contact distance (reduced units)
	! ***************************************************************************************
	sigmahgo = ( 1.d0 - ( 0.5d0 * chi ) * ( ( ( rei + rej ) ** ( 2.d0 ) ) / ( 1.d0 + ( chi * eiej ) ) + &
	&	   ( ( rei - rej ) ** ( 2.d0 ) ) / ( 1.d0 - ( chi * eiej ) ) ) ) ** ( -0.5d0 )

	return

end function sigmahgo
