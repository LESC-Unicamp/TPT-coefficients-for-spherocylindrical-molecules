! ############################################################################################### !
!                 Canonical Monte Carlo algorithm for spherocylindrical molecules                 !
!                  This code contains all subroutines used in the main program.                   !
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
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                                            T. Kihara                                            !
!                                    J. Phys. Soc. Japan (1951)                                   !
!                                     DOI: 10.1143/jpsj.6.289                                     !
!                             --------------------------------------                              !
!                                 C. R. A. Abreu, F. A. Escobedo                                  !
!                                J. Chem. Phys. 124, 054116 (2006)                                !
!                                     DOI: 10.1063/1.2165188                                      !
!                             --------------------------------------                              !
!                                           O. K. Smith                                           !
!                           Communications of the ACM, 4(4), 168 (1961)                           !
!                                    DOI: 10.1145/355578.366316                                   !
!                             --------------------------------------                              !
!                                          G. Marsaglia                                           !
!                            Ann. Math. Statist. 43(2), 645-646 (1972)                            !
!                                  DOI: 10.1214/aoms/1177692644                                   !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!                               Linear congruential generator (LCG)                               !
!    This subroutine generates a random number from the uniform distribution over the range       !
!   0 ≤ x < 1. It does not take any arguments. The number generator seed has an in/out intent,    !
!             i. e., its value is changed every time the ranf() subroutine is called.             !
!                  This seed is the heart of the pseudorandom number generator                   !
!                           and responsible for ensuring repeatability.                           !
!         See Allen and Tildesley, 2nd Edition (2017), Appendix E, for more information.          !
! *********************************************************************************************** !
subroutine ranf()

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Integer parameters
	! ***************************************************************************************
	integer*8, parameter	:: l = 1029
	integer*8, parameter	:: c = 221591
	integer*8, parameter	:: m = 1048576

	! ***************************************************************************************
	! Finite modulus arithmetic
	! ***************************************************************************************
	seed	 = mod( ( ( seed * l ) + c ), m )	! Number generator seed
	random_n = real( seed ) / real( m )		! Pseudorandom number

	return

end subroutine ranf

! *********************************************************************************************** !
!    This subroutine takes the body-fixed components (efixed) and the unit quaternion (qi) and    !
!           generates an active rotation/transformation (ei) via a 3D rotation matrix.            !
!       See Allen and Tildesley, 2nd Edition (2017), pages 106-111 for more information.          !
! *********************************************************************************************** !
subroutine active_transformation(efixed,qp,ep)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Integer variables
	! ***************************************************************************************
	integer*8		:: i, j		! Counters

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8, dimension (3)	:: efixed	! Body-fixed basis vector
	real*8, dimension (0:3)	:: qp		! Unit quaternion
	real*8, dimension (3)	:: ep		! Active rotation
	real*8, dimension (3,3)	:: A		! Rotation matrix
	real*8, dimension (3,3)	:: At		! Transpose of rotation matrix

	! ***************************************************************************************
	! Rotation Matrix A - Allen and Tildesley, 2nd edition, page 110
	! ***************************************************************************************

	! First row
	A(1,1) = ( qp(0) * qp(0) ) + ( qp(1) * qp(1) ) - ( qp(2) * qp(2) ) - ( qp(3) * qp(3) )
	A(1,2) = 2.d0 * ( ( qp(1) * qp(2) ) + ( qp(0) * qp(3) ) )
	A(1,3) = 2.d0 * ( ( qp(1) * qp(3) ) - ( qp(0) * qp(2) ) )

	! Second row
	A(2,1) = 2.d0 * ( ( qp(1) * qp(2) ) - ( qp(0) * qp(3) ) )
	A(2,2) = ( qp(0) * qp(0) ) - ( qp(1) * qp(1) ) + ( qp(2) * qp(2) ) - ( qp(3) * qp(3) )
	A(2,3) = 2.d0 * ( ( qp(2) * qp(3) ) + ( qp(0) * qp(1) ) )

	! Third row
	A(3,1) = 2.d0 * ( ( qp(1) * qp(3) ) + ( qp(0) * qp(2) ) )
	A(3,2) = 2.d0 * ( ( qp(2) * qp(3) ) - ( qp(0) * qp(1) ) )
	A(3,3) = ( qp(0) * qp(0) ) - ( qp(1) * qp(1) ) - ( qp(2) * qp(2) ) + ( qp(3) * qp(3) )

	! ***************************************************************************************
	! Transpose of rotation matrix
	! ***************************************************************************************
	do i = 1, 3
		do j = 1, 3
			At(i,j) = A(j,i)
		end do
	end do

	! ***************************************************************************************
	! Active tranformation (dot product of body-fixed vector and transpose of rotation matrix)
	! ***************************************************************************************
	ep(1) = ( efixed(1) * At(1,1) ) + ( efixed(2) * At(1,2) ) + ( efixed(3) * At(1,3) )
	ep(2) = ( efixed(1) * At(2,1) ) + ( efixed(2) * At(2,2) ) + ( efixed(3) * At(2,3) )
	ep(3) = ( efixed(1) * At(3,1) ) + ( efixed(2) * At(3,2) ) + ( efixed(3) * At(3,3) )

	return

end subroutine active_transformation

! *********************************************************************************************** !
!   This subroutine computes the total potential energy for the initial molecular configuration   !
! *********************************************************************************************** !
subroutine compute_total_energy()

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Integer variables
	! ***************************************************************************************
	integer*8				:: i, j		! Counters

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8				:: sigmavl	! Vega-Lago shortest distance (function)
	real*8				:: rvl		! Vega-Lago shortest distance (variable)
	real*8				:: rijsq	! Vector distance between particles i and j (squared)
	real*8				:: modrij	! Magnitude of the vector distance between particles i and j
	real*8, dimension (3)		:: rij		! Vector distance between particles i and j
	real*8, dimension (3)		:: ei, ej	! Orientation of particles i and j
	real*8, dimension (n_lambda)	:: vij_sw	! Pair potential energy (Square-Well potential)
	real*8, dimension (n_n)		:: vij_kh	! Pair potential energy (Kihara potential)

	! ***************************************************************************************
	! Initialization of array
	! ***************************************************************************************
	v(:) = 0.d0

	! First loop represents a particle with a fixed index i
	do i = 1, n_particles - 1

		! Second loop represents all other particles with indexes j
		do j = i + 1, n_particles

			! ***********************************************************************
			! Orientation of particles i and j
			! ***********************************************************************
			ei(:)   = e(:,i)
			ej(:)   = e(:,j)
			! ***********************************************************************
			! Vector distance between particles i and j
			! ***********************************************************************
			rij(:) = r(:,i) - r(:,j)
			! ***********************************************************************
			! Minimum Image Convention
			! (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
			! ***********************************************************************
			rij(:) = rij(:) - ( box_length(:) ) * anint( rij(:) / box_length(:) )
			! ***********************************************************************
			! Magnitude of the vector distance (squared)
			! ***********************************************************************
			rijsq  = ( rij(1) * rij(1) ) + ( rij(2) * rij(2) ) + ( rij(3) * rij(3) )
			! ***********************************************************************
			! Magnitude of the vector distance
			! ***********************************************************************
			modrij = dsqrt( rijsq )
			! ***********************************************************************
			! Vega-Lago shortest distance
			!  See 'functions' file for more information
			! ***********************************************************************
			rvl = sigmavl(ei,ej,rij,rijsq)
			! ***********************************************************************
			! Compute potential energy via the selected force field
			! ***********************************************************************
			if ( ff_selec(1) ) then
				! ***************************************************************
				! Compute the Square-Well potential
				! ***************************************************************
				call compute_potential_sw(modrij,vij_sw)
				! ***************************************************************
				! Total potential energy (iterative calculation)
				! ***************************************************************
				v(:) = v(:) + vij_sw(:)
			else if ( ff_selec(2) ) then
				! ***************************************************************
				! Compute the Kihara potential
				! ***************************************************************
				call compute_potential_kh(rvl,vij_kh)
				! ***************************************************************
				! Total potential energy (iterative calculation)
				! ***************************************************************
				v(:) = v(:) + vij_kh(:)
			end if

		end do

	end do

	return

end subroutine compute_total_energy

! *********************************************************************************************** !
!              This subroutine computes the potential energy of a random particle i               !
! *********************************************************************************************** !
subroutine compute_particle_energy(i,ei,ri)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Integer variables
	! ***************************************************************************************
	integer*8				:: i, j		! Counters

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8				:: sigmavl	! Vega-Lago shortest distance (function)
	real*8				:: rvl		! Vega-Lago shortest distance (variable)
	real*8				:: rijsq	! Vector distance between particles i and j (squared)
	real*8				:: modrij	! Magnitude of the vector distance between particles i and j
	real*8, dimension (3)		:: ri		! Position of particle i
	real*8, dimension (3)		:: rij		! Vector distance between particles i and j
	real*8, dimension (3)		:: ei		! Orientation of particle i
	real*8, dimension (n_lambda)	:: vij_sw	! Pair potential energy (Square-Well potential)
	real*8, dimension (n_n)		:: vij_kh	! Pair potential energy (Kihara potential)

	! ***************************************************************************************
	! Initialization of array
	! ***************************************************************************************
	vi(:) = 0.d0

	! ***************************************************************************************
	! Particle loops
	! ***************************************************************************************
	!  The loops below analyze the neighboring particles j around a fixed particle i
	!  and compute the pair potential between i and j.
	! ***************************************************************************************

	! First loop takes only particles whose j-indexes are below the i-index of the fixed particle
	do j = 1, i - 1

		! *******************************************************************************
		! Vector distance between particles i and j
		! *******************************************************************************
		rij(:) = ri(:) - rmc(:,j)
		! *******************************************************************************
		! Minimum Image Convention
		! (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
		! *******************************************************************************
		rij(:) = rij(:) - box_length(:) * anint( rij(:) / box_length(:) )
		! *******************************************************************************
		! Magnitude of the vector distance (squared)
		! *******************************************************************************
		rijsq  = ( rij(1) * rij(1) ) + ( rij(2) * rij(2) ) + ( rij(3) * rij(3) )
		! *******************************************************************************
		! Magnitude of the vector distance
		! *******************************************************************************
		modrij = dsqrt( rijsq )
		! ***********************************************************************
		! Vega-Lago shortest distance
		!  See 'functions' file for more information
		! ***********************************************************************
		rvl = sigmavl(ei,emc(:,j),rij,rijsq)
		! ***********************************************************************
		! Compute potential energy via the selected force field
		! ***********************************************************************
		if ( ff_selec(1) ) then
			! ***************************************************************
			! Compute the Square-Well potential
			! ***************************************************************
			call compute_potential_sw(modrij,vij_sw)
			! ***************************************************************
			! Total potential energy (iterative calculation)
			! ***************************************************************
			vi(:) = vi(:) + vij_sw(:)
		else if ( ff_selec(2) ) then
			! ***************************************************************
			! Compute the Kihara potential
			! ***************************************************************
			call compute_potential_kh(rvl,vij_kh)
			! ***************************************************************
			! Total potential energy (iterative calculation)
			! ***************************************************************
			vi(:) = vi(:) + vij_kh(:)
		end if

	end do

	! Second loop takes only particles whose j-indexes are above the i-index of the fixed particle
	do j = i + 1, n_particles

		! *******************************************************************************
		! Vector distance between particles i and j
		! *******************************************************************************
		rij(:) = ri(:) - rmc(:,j)
		! *******************************************************************************
		! Minimum Image Convention
		! (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
		! *******************************************************************************
		rij(:) = rij(:) - box_length(:) * anint( rij(:) / box_length(:) )
		! *******************************************************************************
		! Magnitude of the vector distance (squared)
		! *******************************************************************************
		rijsq  = ( rij(1) * rij(1) ) + ( rij(2) * rij(2) ) + ( rij(3) * rij(3) )
		! *******************************************************************************
		! Magnitude of the vector distance
		! *******************************************************************************
		modrij = dsqrt( rijsq )
		! ***********************************************************************
		! Vega-Lago shortest distance
		!  See 'functions' file for more information
		! ***********************************************************************
		rvl = sigmavl(ei,emc(:,j),rij,rijsq)
		! ***********************************************************************
		! Compute potential energy via the selected force field
		! ***********************************************************************
		if ( ff_selec(1) ) then
			! ***************************************************************
			! Compute the Square-Well potential
			! ***************************************************************
			call compute_potential_sw(modrij,vij_sw)
			! ***************************************************************
			! Total potential energy (iterative calculation)
			! ***************************************************************
			vi(:) = vi(:) + vij_sw(:)
		else if ( ff_selec(2) ) then
			! ***************************************************************
			! Compute the Kihara potential
			! ***************************************************************
			call compute_potential_kh(rvl,vij_kh)
			! ***************************************************************
			! Total potential energy (iterative calculation)
			! ***************************************************************
			vi(:) = vi(:) + vij_kh(:)
		end if

	end do

	return

end subroutine compute_particle_energy

! *********************************************************************************************** !
!              This subroutine computes the pair potential between particles i and j              !
!           It applies a discrete square-well potential to compute the pair potential.            !
!             In this case, we only consider the attractive part of the SW potential.             !
! Although we have previously determined the diameter of a sphere with the same molecular volume  !
!  as the ellipsoid of revolution (see 'main' program), the spherical geometry is non-physical,   !
!        that is, no physical overlap is generated when two spheres penetrate each other.         !
!     The spherical geometry only represents an attractive zone around the hard ellipsoid of      !
! revolution, which can be extended by the attractive range parameter (λ). In this case, two hard !
!      ellipsoids of revolution will only interact attractively when these 'extended zones'       !
!           overlap each other. This will thus represent a region of mutual attraction.           !
!  When there are no overlaps between these zones, the pair potential of particles i and j is 0.  !
!                         Otherwise, it is equal to -ϵ (the well depth).                          !
! *********************************************************************************************** !
subroutine compute_potential_sw(modrij,vij)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8				:: modrij	! Magnitude of the vector distance between particles i and j
	real*8, dimension (n_lambda)	:: vij		! Pair potential

	! ***************************************************************************************
	! Pair potential for every attractive range parameter (λ)
	! ***************************************************************************************
	!  These potentials can be computed simultaneously because the reference system remains
	!  unchanged. Please note that the reference system generates the non-overlapping
	!  configurations. Thus we only have to change the effective range of attraction,
	!  'swrange', for every λ, and compute all potentials simultaneously after every accepted
	!  trial move. Therefore, we do not need to run several simulations for different λ
	!  to achieve this.
	! ***************************************************************************************
	do counter_lambda = 1, n_lambda

		! *******************************************************************************
		! Overlap in the region of attraction
		! *******************************************************************************
		if ( modrij <= swrange(counter_lambda) ) then
			! ***********************************************************************
			! Pair potential (reduced units)
			! ***********************************************************************
			!  The fundamental unit of energy is ϵ, thus all energy quantities are 
			!  divided by ϵ.
			! ***********************************************************************
			vij(counter_lambda) = -1.d0
		else
			vij(counter_lambda) = 0.d0	! No overlap
		end if

	end do

	return

end subroutine compute_potential_sw

! *********************************************************************************************** !
!              This subroutine computes the pair potential between particles i and j              !
!       It applies a continuous anharmonic Kihara potential to compute the pair potential.        !
!     The Kihara potential is based on the shortest distance between the molecular hard cores.    !
!                 See T. Kihara, J. Phys. Soc. Japan (1951), for more information.                !
! *********************************************************************************************** !
subroutine compute_potential_kh(shortest_d,vij)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8			:: shortest_d	! Shortest distance between particles i and j
	real*8			:: mnum, nnum 	! Auxiliars
	real*8			:: m_attractive	! Attractive index
	real*8, dimension (n_n)	:: vij		! Pair potential

	m_attractive = 6.d0

	! ***************************************************************************************
	! Pair potential for every repulsive parameter (n)
	! ***************************************************************************************
	!  These potentials can be computed simultaneously because the reference system remains
	!  unchanged. Please note that the reference system generates the non-overlapping
	!  configurations.
	! ***************************************************************************************
	do counter_n = 1, n_n

		mnum = ( m_attractive ) / ( n_repulsive(counter_n) - m_attractive )
		nnum = ( n_repulsive(counter_n) ) / ( n_repulsive(counter_n) - m_attractive )

		vij(counter_n) = mnum * ( ( 1.d0 / shortest_d ) ** ( n_repulsive(counter_n) ) ) - &
		&		 nnum * ( ( 1.d0 / shortest_d ) ** ( m_attractive ) )

	end do

	return

end subroutine compute_potential_kh

! *********************************************************************************************** !
!     This subroutine takes a reference unit quaternion (qm) and combines it with a randomly      !
!        generated quaternion (qr) through quaternion multiplication, creating a randomly         !
!                                  composed unit quaternion (qn)                                  !
! *********************************************************************************************** !
subroutine composed_quaternion(qm,qn)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8, dimension (0:3) :: qm	! Reference unit quaternion
	real*8, dimension (0:3) :: qr	! Random quaternion
	real*8, dimension (0:3) :: qn	! Composed unit quaternion

	! ***************************************************************************************
	! Random quaternion generator
	! ***************************************************************************************
	call random_quaternion(qr)

	! ***************************************************************************************
	! Quaternion multiplication (composed rotation)
	! ***************************************************************************************
	call multiply_quaternions(qr,qm,qn)

	return

end subroutine composed_quaternion

! *********************************************************************************************** !
!        This subroutine generates a random quaternion from a random angle and random axis        !
! *********************************************************************************************** !
subroutine random_quaternion(qr)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8, dimension (0:3) :: qr	! Random quaternion
	real*8, dimension (3)   :: sr	! Random vector

	! ***************************************************************************************
	! Random rotation angle
	! ***************************************************************************************
	call ranf()
	random_angle = ( ( 2.d0 * random_n ) - 1.d0 ) * angmax	! Range [-angmax,angmax]

	! ***************************************************************************************
	! Random vector generator
	! ***************************************************************************************
	call random_vector(sr)

	! ***************************************************************************************
	! Quaternion algebra
	! ***************************************************************************************
	qr(0) = dcos( random_angle * 0.5d0 )		! Real part
	qr(1) = dsin( random_angle * 0.5d0 ) * sr(1)	! Imaginary part (Vector)
	qr(2) = dsin( random_angle * 0.5d0 ) * sr(2)	! Imaginary part (Vector)
	qr(3) = dsin( random_angle * 0.5d0 ) * sr(3)	! Imaginary part (Vector)

	return

end subroutine random_quaternion

! *********************************************************************************************** !
!            This subroutine generates a random vector on the surface of a unit sphere            !
!           See Allen and Tildesley, 2nd Edition (2017), page 514 for more information.           !
!        (Routine 'maths_module.f90' of Code A.1. in Marsaglia, Ann. Math. Statist., 1972)        !
! *********************************************************************************************** !   
subroutine random_vector(sr)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8, dimension (3)	:: sr			! Random vector
	real*8			:: csi_1, csi_2, zeta	! Random numbers

	! ***************************************************************************************
	! Marsaglia's routine
	! ***************************************************************************************
	zeta_loop: do
		! *******************************************************************************
		! Uniform random number, ξ₁
		! *******************************************************************************
		call ranf()
		csi_1 = ( 2.d0 * random_n ) - 1.d0
		! *******************************************************************************
		! Uniform random number, ξ₂
		! *******************************************************************************
		call ranf()
		csi_2 = ( 2.d0 * random_n ) - 1.d0
		! *******************************************************************************
		! Sum of squares
		! *******************************************************************************
		zeta  = ( csi_1 * csi_1 ) + ( csi_2 * csi_2 )
		! *******************************************************************************
		! Marseglia's criterion
		! *******************************************************************************
		if ( zeta < 1.d0 ) then
			exit zeta_loop
		end if
	end do zeta_loop

	! ***************************************************************************************
	! Random vector
	! ***************************************************************************************
	sr(1) = ( 2.d0 * csi_1 ) * dsqrt( 1.d0 - zeta )
	sr(2) = ( 2.d0 * csi_2 ) * dsqrt( 1.d0 - zeta )
	sr(3) = 1.d0 - ( 2.d0 * zeta )

	return

end subroutine random_vector

! *********************************************************************************************** !
!             This subroutine creates a composed rotation quaternion by multiplying a             !
!                          reference quaternion and a random quaternion                           !
! *********************************************************************************************** !
subroutine multiply_quaternions(qr,qm,qn)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8, dimension (0:3)	:: qr	! Random quaternion
	real*8, dimension (0:3)	:: qm	! Reference quaternion
	real*8, dimension (0:3)	:: qn	! Composed quaternion

	! ***************************************************************************************
	! Cross product of quaternions (qr × qm)
	! ***************************************************************************************
	qn(0) = ( qr(0) * qm(0) ) - ( qr(1) * qm(1) ) - ( qr(2) * qm(2) ) - ( qr(3) * qm(3) )
	qn(1) = ( qr(1) * qm(0) ) + ( qr(0) * qm(1) ) - ( qr(3) * qm(2) ) + ( qr(2) * qm(3) )
	qn(2) = ( qr(2) * qm(0) ) + ( qr(3) * qm(1) ) + ( qr(0) * qm(2) ) - ( qr(1) * qm(3) )
	qn(3) = ( qr(3) * qm(0) ) - ( qr(2) * qm(1) ) + ( qr(1) * qm(2) ) + ( qr(0) * qm(3) )

	return

end subroutine multiply_quaternions

! *********************************************************************************************** !
! This subroutine checks if a random displacement (translation or rotation) of a fixed particle i !
!                           causes any overlaps with other particles j                            !
! *********************************************************************************************** !
subroutine check_overlap(i,ei,ri,overlap)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Integer variables
	! ***************************************************************************************
	integer*8		:: i, j		! Counters

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8			:: sigmavl	! Vega-Lago shortest distance (function)
	real*8			:: rvl		! Vega-Lago shortest distance (variable)
	real*8			:: rijsq	! Vector distance between particles i and j (squared)
	real*8			:: cutoff	! Cutoff diameter
	real*8			:: cutoffsq	! Cutoff diameter (squared)
	real*8, dimension (3)	:: ri		! Position of particle i
	real*8, dimension (3)	:: ei		! Orientation of particle i
	real*8, dimension (3)	:: rij		! Vector distance between particles i and j

	! ***************************************************************************************
	! Logical variables
	! ***************************************************************************************
	logical			:: overlap	! Detects overlap between two particles (VL)	: TRUE = overlap detected;  FALSE = overlap not detected

	! ***************************************************************************************
	! Initialization of overlap detector
	! ***************************************************************************************
	overlap = .false.

	! ***************************************************************************************
	! Cutoff diameter
	! ***************************************************************************************
	!  A cutoff distance equivalent to the diameter of a spherical geometry circumscribing
	!  an spherocylinder.
	! ***************************************************************************************
	!  REMINDER: the fundamental unit of length is 'D'.
	! ***************************************************************************************
	cutoff   = ( 1.d0 + aspect_ratio )	! or 'L+D' for dimensional variables
	cutoffsq = cutoff * cutoff

	! ***************************************************************************************
	! Particle loops
	! ***************************************************************************************
	!  The loops below analyze the neighboring particles j around a fixed particle i
	!  and compute the contact distance (VL) between i and j.
	! ***************************************************************************************

	! First loop takes only particles whose j-indexes are below the i-index of the fixed particle
	do j = 1, i - 1
		! *******************************************************************************
		! Vector distance between particles i and j
		! *******************************************************************************
		rij(:) = ri(:) - rmc(:,j)
		! *******************************************************************************
		! Minimum Image Convention
		! (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
		! *******************************************************************************
		rij(:) = rij(:) - ( box_length(:) ) * anint( rij(:) / box_length(:) )
		! *******************************************************************************
		! Magnitude of the vector distance (squared)
		! *******************************************************************************
		rijsq  = ( rij(1) * rij(1) ) + ( rij(2) * rij(2) ) + ( rij(3) * rij(3) )
		! *******************************************************************************
		! Cutoff distance
		! *******************************************************************************
		if ( rijsq <= cutoffsq ) then
			! ***********************************************************************
			! Shortest distance - Vega-Lago Overlap Criterion
			!  See 'functions' file for more information
			! ***********************************************************************
			rvl = sigmavl(ei,emc(:,j),rij,rijsq)
			! ***********************************************************************
			! Overlap criterion
			! ***********************************************************************
			if ( rvl <= 1.d0 ) then
				overlap = .true.	! Overlap detected
				return			! Return immediately if an overlap is detected
			end if
		end if
    	end do

	! Second loop takes only particles whose j-indexes are above the i-index of the fixed particle
	do j = i + 1, n_particles
		! *******************************************************************************
		! Vector distance between particles i and j
		! *******************************************************************************
		rij(:) = ri(:) - rmc(:,j)
		! *******************************************************************************
		! Minimum Image Convention
		! (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
		! *******************************************************************************
		rij(:) = rij(:) - ( box_length(:) ) * anint( rij(:) / box_length(:) )
		! *******************************************************************************
		! Magnitude of the vector distance (squared)
		! *******************************************************************************
		rijsq  = ( rij(1) * rij(1) ) + ( rij(2) * rij(2) ) + ( rij(3) * rij(3) )
		! *******************************************************************************
		! Cutoff distance
		! *******************************************************************************
		if ( rijsq <= cutoffsq ) then
			! ***********************************************************************
			! Shortest distance - Vega-Lago Overlap Criterion
			!  See 'functions' file for more information
			! ***********************************************************************
			rvl = sigmavl(ei,emc(:,j),rij,rijsq)
			! ***********************************************************************
			! Overlap criterion
			! ***********************************************************************
			if ( rvl <= 1.d0 ) then
				overlap = .true.	! Overlap detected
				return			! Return immediately if an overlap is detected
			end if
		end if
	end do

	return	! No overlaps detected

end subroutine check_overlap

! *********************************************************************************************** !
!   This subroutine resets all Monte Carlo parameters, including acceptance ratios and logical    !
! variables, and returns the system to its initial configuration (the potential energy included)  !
! *********************************************************************************************** !
subroutine reset_system()

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Resetting of Monte Carlo parameters
	! ***************************************************************************************
	resetmc   = .false.	! Monte Carlo parameters			(initial value)
	mov_trans = .false.	! Translational move selector			(initial value)
	mov_rot   = .false.	! Rotational move selector			(initial value)
	drmax     = max_trans	! Maximum translational displacement		(initial value)
	angmax    = max_rot	! Maximum rotational displacement		(initial value)
	nacct     = 0		! Translational move acceptance counter		(initial value)
	naccr     = 0		! Rotational move acceptance counter		(initial value)
	movt      = 0		! Translational move counter			(initial value)
	movr      = 0		! Rotational move counter			(initial value)
	qmc(:,:)  = q(:,:)	! Quaternion algebra				(initial value)
	rmc(:,:)  = r(:,:)	! Position of particles				(initial value)
	emc(:,:)  = e(:,:)	! Orientation of particles			(initial value)
	vmc(:)    = v(:)	! Total potential energy			(initial value)

	! ***************************************************************************************
	! Resetting of random number generator seed
	! ***************************************************************************************
	seed = 123456

	return

end subroutine reset_system

! *********************************************************************************************** !
!    This subroutine calculates the order parameter of a nematic phase via the Q-tensor method    !
! *********************************************************************************************** !
subroutine order_parameter()

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Integer variables
	! ***************************************************************************************
	integer*8		:: i		! Counter
	integer*8		:: alpha, beta	! Unit vector specifiers (î, ĵ, k̂)

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8, dimension (3,3)	:: Qabi		! Order tensor Q of particle i (3 x 3 Matrix)
	real*8, dimension (3,3)	:: Qab		! Order tensor Q of particle i (Average)
	real*8, dimension (3)	:: eigenvalue	! Eigenvalues of order tensor Q
	real*8			:: kronecker 	! Kronecker delta
	real*8			:: m		! Trace (matrix)
	real*8			:: det		! Determinant (matrix)
	real*8			:: p		! Sum of squares of elements of a matrix
	real*8			:: hq, phi, pq	! Auxiliar variables

	! ***************************************************************************************
	! Initialization of Q matrix
	! ***************************************************************************************
	Qabi(:,:) = 0.d0	

	! ***************************************************************************************
	! Matrix construction
	! ***************************************************************************************
	do i = 1, n_particles
		do alpha = 1, 3
			do beta = 1, 3
				! ***************************************************************
				! Dyadic product of phase director and orientation
				! ***************************************************************
				!  0, if i ≠ j (perpendicular) and 1, if i = j (parallel)
				! ***************************************************************
				if ( alpha == beta ) then
					kronecker = 1.d0
				else
					kronecker = 0.d0
				end if
				! ***************************************************************
				! Legendre polynomial of the second order (P₂)
				! ***************************************************************
				Qabi(alpha,beta) = Qabi(alpha,beta) + ( 1.5d0 * emc(alpha,i) * emc(beta,i) ) - &
				&		   ( 0.5d0 * kronecker )
			end do
		end do
	end do

	! ***************************************************************
	! Arithmetic mean
	! ***************************************************************
	Qab(:,:) = Qabi(:,:) / dble( n_particles )

	! ***************************************************************
	! Eigenvalues of a Symmetric 3 x 3 Matrix
	! See Smith, Communications of the ACM, 1961, for more information
	! ***************************************************************

	! Trace of symmetric Q matrix
	m   = ( Qab(1,1) + Qab(2,2) + Qab(3,3) )

	! One-third of the trace of the symmetric Q matrix
	m   = m / 3.d0

	! Determinant of (Q - mI), where I is the identity matrix
	det = ( Qab(1,1) - m ) * ( Qab(2,2) - m ) * ( Qab(3,3) - m ) + ( Qab(1,2) * Qab(2,3) * Qab(3,1) ) + &
	&     ( Qab(1,3) * Qab(3,2) * Qab(2,1) ) - ( Qab(1,3) * ( Qab(2,2) - m ) * Qab(3,1) ) - &
	&     ( Qab(2,3) * Qab(3,2) * ( Qab(1,1) - m ) ) - ( ( Qab(3,3) - m ) * Qab(2,1) * Qab(1,2) )

	! Half of the determinant of (Q - mI)
	hq  = 0.5d0 * det

	! Sum of squares of elements of (Q - mI)
	p   = ( Qab(1,1) - m ) * ( Qab(1,1) - m ) + ( Qab(1,2) * Qab(1,2) ) + ( Qab(1,3) * Qab(1,3) ) + &
	&     ( Qab(2,2) - m ) * ( Qab(2,2) - m ) + ( Qab(2,1) * Qab(2,1) ) + ( Qab(2,3) * Qab(2,3) ) + &
	&     ( Qab(3,3) - m ) * ( Qab(3,3) - m ) + ( Qab(3,1) * Qab(3,1) ) + ( Qab(3,2) * Qab(3,2) )

	! One-sixth of the sum of squares of elements of (Q - mI)
	p   = p / 6.d0

	! Test condition
	pq  = ( p * p * p ) - ( hq * hq )

	! ***************************************************************
	! Real eigenvalues condition (p³ ≥ hq²)
	! ***************************************************************
	if ( pq >= 0.d0 ) then
		phi = datan( dsqrt( pq ) / hq ) / 3.d0	! 0 ≤ ϕ ≤ π
	else
		phi = 0.d0
	end if

	! ***************************************************************
	! Eigenvalues of Q
	! ***************************************************************
	eigenvalue(1) = m + ( 2.d0 * dsqrt( p ) * dcos( phi ) )
	eigenvalue(2) = m - ( dsqrt( p ) ) * ( dcos( phi ) + ( dsqrt( 3.d0 ) * dsin( phi ) ) )
	eigenvalue(3) = m - ( dsqrt( p ) ) * ( dcos( phi ) - ( dsqrt( 3.d0 ) * dsin( phi ) ) )

	! ***************************************************************
	! Nematic order parameter
	! ***************************************************************
	s = maxval( eigenvalue )	! Largest eigenvalue

	return

end subroutine order_parameter

! *********************************************************************************************** !
!     This subroutine uses a block-averaging method to calculate the first- and second-order      !
!         TPT coefficients, the perturbed Helmholtz free energy, and their uncertainties          !
! *********************************************************************************************** !
!                        Original developer: Luis Fernando Mercier Franco                         !
!                     University of Campinas, School of Chemical Engineering                      !
! *********************************************************************************************** !
!         See Allen and Tildesley, 2nd Edition (2017), pages 281-287 for more information         !
! *********************************************************************************************** !
subroutine block_averaging(flag,a1,a2,apert,dpa1,dpa2,dpapert)

	! Uses one module: global variables
	use globalvar

	implicit none

	! ***************************************************************************************
	! Integer variables
	! ***************************************************************************************
	integer*8, parameter	:: max_data = 1e4		! Maximum number of block data for linear regression
	integer*8		:: equil_save			! Equilibration data points in the potential file
	integer*8		:: n_run			! Production data points in the potential file
	integer*8		:: n_block			! Counter for the number of blocks
	integer*8		:: n_data			! Number of data points in a block
	integer*8		:: i, j, k, count_aux, step	! Auxiliary counters

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8			:: perturbed_moment1			! First moment (average),	<exp(-U*/T*)>
	real*8			:: perturbed_moment2			! Second moment,		<exp(-2U*/T*)>
	real*8			:: moment1				! First moment (average),	<U*>
	real*8			:: moment2				! Second moment,		<U*²>
	real*8			:: moment3				! Third moment,			<U*³>
	real*8			:: moment4				! Fourth moment,		<U*⁴>
	real*8			:: perturbed_var1_tot			! Variance of the 1st moment,	<exp(-2U*/T*)> - <exp(-U*/T*)>²
	real*8			:: var1_tot				! Variance of the 1st moment,	<U*²> - <U*>²
	real*8			:: var2_tot				! Variance of the 2nd moment,	<U*⁴> - <U*²>²
	real*8			:: perturbed_avg_blk, avg_blk		! Block average
	real*8			:: perturbed_var_block, var_block	! Block variance
	real*8			:: perturbed_var_sum, var_sum		! Auxiliary summation variable for the block variance
	real*8			:: step_blocks(max_data)		! Inverse of the number of data points in each block
	real*8			:: perturbed_sb(max_data), sb(max_data)	! Statistical inefficiency per number of blocks
	real*8			:: a, b					! Parameters of linear regression, s = f(step_blocks) = a*step_blocks + b
	real*8			:: r2, perturbed_r2			! Coefficient of determination
	real*8			:: perturbed_s_run, s_run		! Statistical inefficiency
	real*8			:: perturbed_sigsq1			! Variance of <exp(-U*/T*)>
	real*8			:: sigsq1				! Variance of <U*>
	real*8			:: sigsq2				! Variance of <U*²>
	real*8			:: cov					! Covariance between U* and U*²
	real*8			:: corr					! Correlation between U* and U*²
	real*8			:: expargument				! Boltzmann factor argument
	real*8			:: apert				! Perturbed Helmholtz free energy
	real*8			:: a1, a2				! TPT coefficients
	real*8			:: dpapert				! Perturbed Helmholtz free energy (standard deviation)
	real*8			:: dpa1, dpa2				! TPT coefficients (standard deviation)
	real*8			:: initial_time				! Initial time
	real*8			:: final_time				! Final time

	! ***************************************************************************************
	! Real variables (allocatable)
	! ***************************************************************************************
	real*8, dimension (:), allocatable	:: pot			! Potential energy (production only)

	! ***************************************************************************************
	! Logical variables
	! ***************************************************************************************
	logical			:: flag					! Generic true/false flag

	! ***************************************************************************************
	! Number of data points
	! ***************************************************************************************
	n_run 	   = int( ( dble( max_cycles ) - dble( n_equil ) ) / dble( n_save ) )	! Production points
	equil_save = int( dble( n_equil ) / dble( n_save ) )				! Equilibration points

	! ***************************************************************************************
	! Initialization of logical-type flag
	! ***************************************************************************************
	flag = .false.
	! Return to main program if number of block data exceeds max_data parameter
	if ( max_blocks - min_blocks >= max_data ) then
		flag = .true.
		return
	end if

	! ***************************************************************************************
	! Initialization of variables
	! ***************************************************************************************
	perturbed_moment1 = 0.d0
	perturbed_moment2 = 0.d0
	moment1		  = 0.d0
	moment2		  = 0.d0
	moment3		  = 0.d0
	moment4		  = 0.d0

	allocate ( pot(n_run) )

	! ***************************************************************************************
	! Start subroutine timer
	! ***************************************************************************************
	call cpu_time(initial_time)

	! ***************************************************************************************
	! Subroutine log
	! ***************************************************************************************
	if ( ff_selec(1) ) then
		open ( UNIT = 100, FILE = "Perturbed_Coefficient/SW/"//trim(descriptor_date)//"/Lambda_"//trim(descriptor_lamb)// &
		&			  "/logfile_ε"//trim(descriptor_file1)//"_ρ"//trim(descriptor_file2)//".dat" )
	else if ( ff_selec(2) ) then
		open ( UNIT = 100, FILE = "Perturbed_Coefficient/KH/"//trim(descriptor_date)//"/n_"//trim(descriptor_lamb)// &
		&			  "/logfile_ε"//trim(descriptor_file1)//"_ρ"//trim(descriptor_file2)//".dat" )
	end if
	write ( 100 , * ) "*******************************************************"
	write ( 100 , * ) "Block average calculation                              "
	write ( 100 , * ) "*******************************************************"
	write ( 100 , * ) "University of Campinas                                 "
	write ( 100 , * ) "School of Chemical Engineering                         "
	write ( 100 , * ) "Luis Fernando Mercier Franco                           "
	write ( 100 , * ) "March 26th, 2021                                       "
	write ( 100 , * ) "*******************************************************"
	write ( 100 , * ) " "
	write ( 100 , * ) "*******************************************************"
	write ( 100 , * ) "The author does not take any liability for the use     "
	write ( 100 , * ) "of this code!                                          "
	write ( 100 , * ) "*******************************************************"
	write ( 100 , * ) " "
	write ( 100 , * ) "*******************************************************"
	write ( 100 , * ) "Number of data points    = ", n_run 
	write ( 100 , * ) "Minimum number of blocks = ", min_blocks
	write ( 100 , * ) "Maximum number of blocks = ", max_blocks
	write ( 100 , * ) "*******************************************************"
	write ( 100 , * ) " "

	! ***************************************************************************************
	! Initialization of the Largest Exponential Argument
	! ***************************************************************************************
	!  This corresponds to a method to calculate the natural logarithm of large exponential 
	!  iterations, avoiding mathematical inaccuracies. This is applied to determine the 
	!  Helmholtz free energy of the perturbed system by taking the negative natural logarithm 
	!  of the canonical ensemble average over the reference fluid: -[ln<exp(-βU¹)>0].
	! ***************************************************************************************
	!  See Supplementary Material of Abreu and Escobedo (2006) for more information.
	! ***************************************************************************************
	maxarg = 0.d0

	! ***************************************************************************************
	! Largest Exponential Argument and Potential Energy
	! ***************************************************************************************
	if ( ff_selec(1) ) then
		open ( UNIT = 150, FILE = "Potential/SW/"//trim(descriptor_date)//"/Lambda_"//trim(descriptor_lamb)// &
		&			  "/thermo_ε"//trim(descriptor_file1)//"_ρ"//trim(descriptor_file2)//".dat" )
	else if ( ff_selec(2) ) then
		open ( UNIT = 150, FILE = "Potential/KH/"//trim(descriptor_date)//"/n_"//trim(descriptor_n)// &
		&			  "/thermo_ε"//trim(descriptor_file1)//"_ρ"//trim(descriptor_file2)//".dat" )
	end if
	! Skip equilibration data points (when necessary)
	if ( .not. potential_check ) then
		do i = 1, equil_save
			read ( 150, * ) anynumber, anynumber
		end do
	end if
	! Read production data points
	do i = 1, n_run
		read ( 150, * ) step, pot(i)
		beta_energy = -pot(i) / temp
		if ( beta_energy > maxarg ) then
			maxarg = beta_energy
		end if
	end do
	close ( 150 )

	! ***************************************************************************************
	! OBSERVATION: Boltzmann Factor Argument 
	! (calculation with reduced quantities)
	! ***************************************************************************************
	!  The Boltzmann factor argument is -βU. Since β = 1/(kB T), then:
	!
	!                                 -βU = -U/(kB T)
	!
	!  The fundamental unit of temperature is ε/kB, where kB is the Boltzmann
	!  constant. The reduced temperature, T*, will be:
	!
	!                      T* = T/(ε/kB), or (kB T)/ε, or 1/βε
	!
	!  In this case:
	!
	!                                   T = εT*/kB
	!
	!  Thus:
	!
	!                   -βU = -U/(kB T) = -(U kB)/(kB εT*) = -U/εT*
	!
	!  But, U* = U/ε is the reduced potential energy. Then:
	!
	!                                  -βU = -U*/T*
	!
	!  which is the Boltzmann factor argument calculated with reduced quantities.
	! ***************************************************************************************

	! ***************************************************************************************
	! Iterative process
	! ***************************************************************************************
	do i = 1, n_run
		expargument	  = ( -pot(i) / temp ) - ( maxarg )
		perturbed_moment1 = perturbed_moment1 + dexp( expargument )
		perturbed_moment2 = perturbed_moment2 + dexp( 2.d0 * expargument )
		moment1		  = moment1 + ( pot(i) )
		moment2		  = moment2 + ( pot(i) * pot(i) )
		moment3		  = moment3 + ( pot(i) * pot(i) * pot(i) )
		moment4		  = moment4 + ( pot(i) * pot(i) * pot(i) * pot(i) )
	end do

	! ***************************************************************************************
	! First moments (average)
	! (calculated using Eq. [8.11], p. 282)
	! ***************************************************************************************
	perturbed_moment1  = perturbed_moment1 / dble( n_run )
	moment1		   = moment1 / dble( n_run )

	! ***************************************************************************************
	! Second moments
	! (calculated using Eq. [8.11], p. 282)
	! ***************************************************************************************
	perturbed_moment2  = perturbed_moment2 / dble( n_run )
	moment2		   = moment2 / dble( n_run )

	! ***************************************************************************************
	! Third moment
	! (calculated using Eq. [8.11], p. 282)
	! ***************************************************************************************
	moment3		   = moment3 / dble( n_run )

	! ***************************************************************************************
	! Fourth moment
	! (calculated using Eq. [8.11], p. 282)
	! ***************************************************************************************
	moment4		   = moment4 / dble( n_run )

	! ***************************************************************************************
	! Variances of <X>
	! (calculated using Eq. [8.13], p. 282 for n_run >> 1)
	! ***************************************************************************************
	perturbed_var1_tot = perturbed_moment2 - ( perturbed_moment1 * perturbed_moment1 )
	var1_tot	   = moment2 - ( moment1 * moment1 )

	! ***************************************************************************************
	! Variance of <U*²>
	! (calculated using Eq. [8.13], p. 282 for n_run >> 1)
	! ***************************************************************************************
	var2_tot	   = moment4 - ( moment2 * moment2 )

	! ***************************************************************************************
	! Covariance
	! ***************************************************************************************
	cov		   = moment3 - ( moment1 * moment2 )

	write ( 100 , * ) "No. Blocks", "     ", "No. Data", "     ", "      var(TPT)      ", &
	&		  "      var(PERT)     ", "       s(TPT)       ", "       s(PERT)      "

 900	format ( I6, 9X, I6, 9X, E15.7, 5X, E15.7, 5X, E15.7, 5X, E15.7 )

	! ***************************************************************************************
	! Block Average
	! (calculated using Eq. [8.14], p. 282)
	! ***************************************************************************************
	do n_block = min_blocks, max_blocks

		n_data = int( dble( n_run ) / dble( n_block ) )	! Steps in each block for a certain number of blocks (n_block)

		! Initialization
		perturbed_var_sum = 0.d0
		var_sum		  = 0.d0
		j		  = 1

		do k = 1, n_block
			perturbed_avg_blk = 0.d0	! Reset
			avg_blk		  = 0.d0	! Reset
			do i = 1, n_data
				expargument	  = ( -pot(j) / temp ) - ( maxarg )
				perturbed_avg_blk = perturbed_avg_blk + dexp( expargument )
				avg_blk		  = avg_blk + pot(j)
				! Increment counter
				j		  = j + 1
			end do
			perturbed_avg_blk = perturbed_avg_blk / dble( n_data )
			perturbed_var_sum = perturbed_var_sum + ( perturbed_avg_blk - perturbed_moment1 ) ** ( 2.d0 )
			avg_blk		  = avg_blk / dble( n_data )
			var_sum		  = var_sum + ( avg_blk - moment1 ) ** ( 2.d0 )
		end do

		! *******************************************************************************
		! Variance of the block averages
		! (calculated using Eq. [8.15], p. 283)
		! *******************************************************************************
		perturbed_var_block = perturbed_var_sum / dble( n_block - 1 )
		var_block	    = var_sum / dble( n_block - 1 )

		! *******************************************************************************
		! Auxiliary counter
		! *******************************************************************************
		count_aux = ( n_block - min_blocks + 1 )

		! *******************************************************************************
		! Inverse of block size
		! *******************************************************************************
		step_blocks(count_aux) = ( 1.d0 ) / dble( n_data )

		! *******************************************************************************
		! Statistical inefficiency
		! (function of the number of blocks)
		! *******************************************************************************
		perturbed_sb(count_aux) = dble( n_data ) * ( perturbed_var_block / perturbed_var1_tot )
		sb(count_aux)	        = dble( n_data ) * ( var_block / var1_tot )
		write ( 100, 900 ) n_block, n_data, var_block, perturbed_var_block, &
		&		   sb(count_aux), perturbed_sb(count_aux)
	end do

	deallocate ( pot )

	! ***************************************************************************************
	! Linear regression (perturbed Helmholtz free energy)
	! ***************************************************************************************
	!  Fits the statistical inefficiency to a function of the inverse of the number of steps 
	!  in a certain block
	! ***************************************************************************************
	call linear_fit(step_blocks,perturbed_sb,a,b,perturbed_r2,count_aux)

	! ***************************************************************************************
	! Statistical inefficiency (perturbed Helmholtz free energy)
	! (calculated as defined by Eq. [8.16], p. 283)
	! ***************************************************************************************
	perturbed_s_run = a

	! ***************************************************************************************
	! Linear regression (TPT coefficients)
	! ***************************************************************************************
	!  Fits the statistical inefficiency to a function of the inverse of the number of steps 
	!  in a certain block
	! ***************************************************************************************
	call linear_fit(step_blocks,sb,a,b,r2,count_aux)

	! ***************************************************************************************
	! Statistical inefficiency (TPT coefficients)
	! (calculated as defined by Eq. [8.16], p. 283)
	! ***************************************************************************************
	s_run = a

	! ***************************************************************************************
	! Variances of <X>
	! (calculated using Eq. [8.17], p. 283)
	! ***************************************************************************************
	perturbed_sigsq1 = perturbed_var1_tot * ( perturbed_s_run / dble( n_run ) )
	sigsq1		 = var1_tot * ( s_run / dble( n_run ) )

	! ***************************************************************************************
	! Variance of <U*²>
	! ***************************************************************************************
	sigsq2 = var2_tot * ( s_run / dble( n_run ) )

	! ***************************************************************************************
	! Covariance between U* and U*²
	! ***************************************************************************************
	cov = cov * ( s_run / dble( n_run ) )

	! ***************************************************************************************
	! Correlation between U* and U*²
	! ***************************************************************************************
	corr = cov / ( dsqrt( sigsq1 ) ) / ( dsqrt( sigsq2 ) )

	! ***************************************************************************************
	! TPT Coefficients
	! ***************************************************************************************
	a1 = moment1 / n_particles
	a2 = ( -0.5d0 ) * ( moment2 - ( moment1 * moment1 ) ) / n_particles

	! ***************************************************************************************
	! TPT Coefficients (Uncertainty propagation)
	! ***************************************************************************************
	dpa1 = dsqrt( sigsq1 ) / n_particles
	dpa2 = ( 0.5d0 ) * ( dsqrt( ( 4.d0 * sigsq1 * moment1 * moment1 ) + ( sigsq2 ) - ( 4.d0 * moment1 * cov ) ) ) / &
	&      ( n_particles )

	! ***************************************************************************************
	! Perturbed Helmholtz free energy
	! ***************************************************************************************
	!  See Supplementary Material of Abreu and Escobedo (2006) for more information on how to
	!  calculate the natural logarithm of large exponential iterations.
	! ***************************************************************************************
	apert   = - ( 1.d0 / n_particles ) * ( dlog( perturbed_moment1 ) + maxarg )

	! ***************************************************************************************
	! Perturbed Helmholtz free energy (Uncertainty propagation)
	! ***************************************************************************************
	dpapert = ( 1.d0 / n_particles ) * ( dsqrt( ( perturbed_sigsq1 ) / ( perturbed_moment1 * perturbed_moment1 ) ) )

	! ***************************************************************************************
	! Subroutine results
	! ***************************************************************************************
	write ( 100 , * ) "****************************************************************"
	write ( 100 , "(A43,E15.7)" ) " Statistical inefficiency (Coefficients) = ", s_run
	write ( 100 , "(A43,E15.7)" ) " Statistical inefficiency (Helmholtz)    = ", perturbed_s_run
	write ( 100 , "(A43,E15.7)" ) " R² (Coefficients)                       = ", r2
	write ( 100 , "(A43,E15.7)" ) " R² (Helmholtz)                          = ", perturbed_r2
	write ( 100 , "(A43,E15.7)" ) " <U*>                                    = ", moment1
	write ( 100 , "(A43,E15.7)" ) " <exp(-U*/T*)>                           = ", perturbed_moment1
	write ( 100 , "(A43,E15.7)" ) " <U*²>                                   = ", moment2
	write ( 100 , "(A43,E15.7)" ) " <exp(-2U*/T*)>                          = ", perturbed_moment2
	write ( 100 , * ) " "
	write ( 100 , * ) "Covariance matrix = _                               _    "
	write ( 100 , "(A22,2E15.7,A3)" ) "                   | ", sigsq1, cov, "  |"
	write ( 100 , "(A22,2E15.7,A3)" ) "                   |_", cov, sigsq2, " _|"
	write ( 100 , * ) " "
	write ( 100 , "(A43,E15.7)" ) " Variance (Helmholtz)                    = ", perturbed_sigsq1
	write ( 100 , * ) " "
	write ( 100 , "(A43,e15.7)" ) " Correlation (Coefficients)              = ", corr
	write ( 100 , * ) "****************************************************************"
	write ( 100 , * ) " "

	! ***************************************************************************************
	! Finish subroutine timer
	! ***************************************************************************************
	call cpu_time(final_time)

	! ***************************************************************************************
	! Elapsed time
	! ***************************************************************************************
	write ( 100, * ) "Elapsed time = ", ( final_time - initial_time ), "s"

	close ( 100 )

	return

end subroutine block_averaging

! *********************************************************************************************** !
!                    This subroutine makes a linear regression of x and y data                    !
! *********************************************************************************************** !
subroutine linear_fit(x,y,a,b,r2,n)

	implicit none

	! ***************************************************************************************
	! Integer variables
	! ***************************************************************************************
	integer*8	:: i, n

	! ***************************************************************************************
	! Real variables
	! ***************************************************************************************
	real*8		:: sx		! Sum of X
	real*8		:: sy		! Sum of Y
	real*8		:: sxx		! Sum of X²
	real*8		:: sxy		! Sum of XY
	real*8		:: syy		! Sum of Y²
	real*8		:: a		! Y-intercept
	real*8		:: b		! Slope
	real*8		:: r2		! Coefficient of determination
	real*8		:: x(n)		! Independent variable
	real*8		:: y(n)		! Dependent variable

	! ***************************************************************************************
	! Initialization
	! ***************************************************************************************
	sx  = 0.d0
	sy  = 0.d0
	sxx = 0.d0
	sxy = 0.d0
	syy = 0.d0

	! ***************************************************************************************
	! Iteration
	! ***************************************************************************************
	do i = 1, n
	   sx  = sx + x(i)
	   sy  = sy + y(i)
	   sxx = sxx + ( x(i) * x(i) )
	   sxy = sxy + ( x(i) * y(i) )
	   syy = syy + ( y(i) * y(i) )
	end do

	! ***************************************************************************************
	! Linear regression (Y = a + bX)
	! ***************************************************************************************
	! Slope
	b  = ( ( dble( n ) * sxy ) - ( sx * sy ) ) / ( ( dble( n ) * sxx ) - ( sx * sx ) )
	! Y-intercept
	a  = ( sy - ( b * sx ) ) / dble( n )
	! Coefficient of correlation
	r2 = ( ( dble( n ) * sxy ) - ( sx * sy ) ) / &
	&    dsqrt( ( ( dble( n ) * sxx ) - ( sx * sx ) ) * ( ( dble( n ) * syy ) - ( sy * sy ) ) )
	! Coefficient of determination
	r2 = ( r2 * r2 )

	return

end subroutine linear_fit
