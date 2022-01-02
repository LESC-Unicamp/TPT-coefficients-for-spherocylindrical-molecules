! ############################################################################################### !
!                 Canonical Monte Carlo algorithm for spherocylindrical molecules                 !
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
!                                        January 2nd, 2022                                        !
! ############################################################################################### !
! Main References:                        C. Vega, S. Lago                                        !
!                                 Computers Chem. 18, 55-59 (1993)                                !
!                                DOI: 10.1016/0097-8485(94)80023-5                                !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!     This function takes the relative orientations of two molecular spherocylinders i and j      !
!   and the vector distance between their centers of mass and calculates the shortest distance    !
!                                       between their cores.                                      !
!          See Vega and Lago, Computers Chem. 18, 55-59 (1993), for more information.             !
! *********************************************************************************************** !
double precision function sigmavl(ei,ej,rij,rijsq)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8, dimension (3)	:: ei			! Orientation of particle i
	real*8, dimension (3)	:: ej			! Orientation of particle j
	real*8, dimension (3)	:: rij			! Vector distance between the centers of mass of particles i and j
	real*8                  :: rijsq		! Squared vector distance between the centers of mass of particles i and j
	real*8			:: rijei		! Dot product of vector distance and orientation of particle i
	real*8			:: rijej		! Dot product of vector distance and orientation of particle j
	real*8			:: eiej			! Dot product of both orientations (particles i and j)
	real*8			:: halfleng		! Half-length of spherocylinder
	real*8			:: dlambda, dmu		! Values that minimize r²
	real*8			:: cc, auxi, auxj	! Auxiliary variables

	! ***************************************************************************************
	! Identical spherocylinders
	! ***************************************************************************************
	halfleng = ( 0.5d0 * aspect_ratio )

	! ***************************************************************************************
	! Dot product of vector distance and orientation of particle i
	! ***************************************************************************************
	rijei = ( rij(1) * ei(1) ) + ( rij(2) * ei(2) ) + ( rij(3) * ei(3) )
	! ***************************************************************************************
	! Dot product of vector distance and orientation of particle j
	! ***************************************************************************************
	rijej = ( rij(1) * ej(1) ) + ( rij(2) * ej(2) ) + ( rij(3) * ej(3) )
	! ***************************************************************************************
	! Dot product of both orientations (particles i and j)
	! ***************************************************************************************
	eiej  = ( ei(1) * ej(1) ) + ( ei(2) * ej(2) ) + ( ei(3) * ej(3) )
	! ***************************************************************************************
	! Denominator 
	! ***************************************************************************************
	cc    = 1.d0 - ( eiej * eiej )

	! ***************************************************************************************
	! Parallel spherocylinders 
	! ***************************************************************************************
	if ( cc < 1.d-10 ) then
		if ( dabs( rijei ) > 1.d-10 ) then		! Parallel spherocylinders not perpendicular to the intermolecular axis
			dlambda = dsign( halfleng, rijei )	! Take the extreme side of the other particle
			dmu = ( dlambda * eiej ) - rijej	! Closest point between particle i and particle j
			if ( dabs( dmu ) > halfleng ) then
				dmu = dsign( halfleng, dmu )
			end if
			! ***********************************************************************
			! Shortest distance (squared)
			! ***********************************************************************
			sigmavl = rijsq + ( dlambda * dlambda ) + ( dmu * dmu ) - ( 2.d0 * dlambda * dmu * eiej ) + &
			&	  ( 2.d0 * dmu * rijej ) - ( 2.d0 * dlambda * rijei )
			return
		else						! Parallel spherocylinders "perfectly" orthogonal to the intermolecular axis
			dlambda = 0.d0
			dmu     = 0.d0
			! ***********************************************************************
			! Shortest distance (squared)
			! ***********************************************************************
			sigmavl = rijsq
			return
		end if
	end if

	! ***************************************************************************************
	! Non-parallel Spherocylinders
	! ***************************************************************************************

	! ***************************************************************************************
	! STEP 1: Evaluation of (λ’, μ’) according to Equations (3) and (4)
	! ***************************************************************************************
	!  See Vega and Lago, Computers Chem. (1994) for more information.
	! ***************************************************************************************
	dlambda = ( rijei - ( eiej * rijej ) ) / cc
	dmu     = ( -rijej + ( eiej * rijei ) ) / cc

	! ***************************************************************************************
	! STEP 2: Check whether the point (λ’, μ’) is in the plane (λ, μ)
	! ***************************************************************************************
	!  See Vega and Lago, Computers Chem. (1994) for more information.
	! ***************************************************************************************
	!   Point (λ’, μ’) in the plane (λ, μ)
	! ***************************************************************************************
	if ( ( dabs( dlambda ) <= halfleng ) .and. ( dabs( dmu ) <= halfleng ) ) then
		sigmavl = rijsq + ( dlambda * dlambda ) + ( dmu * dmu ) - ( 2.d0 * dlambda * dmu * eiej ) + &
		&	      ( 2.d0 * dmu * rijej ) - ( 2.d0 * dlambda * rijei )
		return
	end if
	! ***************************************************************************************
	!   Point (λ’, μ’) not in the plane (λ, μ)
	! ***************************************************************************************
	!   Determination of regions 1-4.
	! ***************************************************************************************
	!    See Fig. 2 of Vega and Lago, Computers Chem. (1994) for more information.
	! ***************************************************************************************
	auxi = dabs( dlambda ) - halfleng
	auxj = dabs( dmu ) - halfleng

	! ***************************************************************************************
	! STEP 3-7: Shortest distance from the side of the region to the line where the other rod 
	!  is contained
	! ***************************************************************************************
	!  See Vega and Lago, Computers Chem. (1994) for more information.
	! ***************************************************************************************
	if ( auxi > auxj ) then		! REGION 1 or 3
		dlambda = dsign( halfleng, dlambda )
		dmu = ( dlambda * eiej ) - rijej
		if ( dabs( dmu ) > halfleng ) then
			dmu = dsign( halfleng, dmu )
		end if
	else				! REGION 2 or 4
		dmu = dsign( halfleng, dmu )
		dlambda = ( dmu * eiej ) + rijei
		if ( dabs( dlambda ) > halfleng ) then
			dlambda = dsign( halfleng, dlambda )
		end if
	end if

	! ***************************************************************************************
	! STEP 8: Evaluate the shortest distance
	! ***************************************************************************************
	sigmavl = rijsq + ( dlambda * dlambda ) + ( dmu * dmu ) - ( 2.d0 * dlambda * dmu * eiej ) + &
	&	      ( 2.d0 * dmu * rijej ) - ( 2.d0 * dlambda * rijei )

	return

end function sigmavl
