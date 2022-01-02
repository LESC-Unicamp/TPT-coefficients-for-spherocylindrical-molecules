! ############################################################################################### !
!              Canonical Monte Carlo algorithm for spherocylindrical molecules                    !
!  This module allows the user to choose one of the stretched/narrowed molecular configurations:  !
!             simple cube (SC), body-centered cube (BCC), or face-centered cube (FCC).            !
!        Molecules will be then allocated in accordance to the selected crystal structure.        !
! Please remember to enter a VALID number of particles, N, in the input file that corresponds to: !
!     an exact cube root of N for the SC configuration, an exact cube root of N/2 for the BCC     !
!             configuration, and an exact cube root of N/4 for the FCC configuration.             !
!                       See Macpherson et al. (2007) for more information.                        !
!     This module also writes out a file containing all particles' positions and quaternions.     !
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
! Main References:                 M. P. Allen, D. J. Tildesley                                   !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             doi: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                               G. Macpherson, M. K. Borg, J. Reese                               !
!              Molecular Simulation, Taylor & Francis, 2007, 33 (15), pp. 1199-1212.              !
!                                 doi: 10.1080/08927020701730724                                  !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

module initconfig

	! Uses one module: global variables
	use globalvar

	implicit none

	! *************************************************************************************** !
	!                                  INITIAL CONFIGURATION                                  !
	!                                    Simple cube = SC                                     !
	!                                Face-centered cube = FCC                                 !
	!                                Body-centered cube = BCC                                 !
	! *************************************************************************************** !

	contains

	! *************************************************************************************** !
	!      This subroutine allows the user to choose the initial molecular configuration      !
	! *************************************************************************************** !
	subroutine config_selection()

		implicit none

		! *******************************************************************************
		! Initialization of logical array
		!  (1) = SC
		!  (2) = BCC
		!  (3) = FCC
		! *******************************************************************************
		config_selec(:) = .false.

		! *******************************************************************************
		! Initial configuration inquiry
		! *******************************************************************************
		config_check_loop: do
			print *, "Select the initial configuration of the particles (SC, BCC, or FCC)."
			print *, "SC (Simple Cube); BCC (Body-Centered Cube); FCC (Face-Centered Cube)"
			read (*,*) config_inq
			if ( config_inq == "sc" .or. config_inq == "SC" ) then
				config_selec(1) = .true.
				exit config_check_loop
			else if ( config_inq == "bcc" .or. config_inq == "BCC" ) then
				config_selec(2) = .true.
				exit config_check_loop
			else if ( config_inq == "fcc" .or. config_inq == "FCC" ) then
				config_selec(3) = .true.
				exit config_check_loop
			! Subroutine returns error if an incorrect configuration is entered, restarting the loop
			else
				print *, "Invalid answer. Answer with: "
				print *, "SC (Simple cube), BCC (Body-centered cube), or FCC (Face-centered cube)."
			end if
		end do config_check_loop

	end subroutine config_selection

	! *************************************************************************************** !
	!     This subroutine allocates particles according to the SC molecular configuration     !
	! *************************************************************************************** !
	subroutine config_sc()

	implicit none

		integer*8	:: i, j, k, counter	! Counters

		allocate ( q (0:3,n_particles) )
		allocate ( r (3,n_particles) )
		allocate ( e (3,n_particles) )

		! *******************************************************************************
		! Unrotated reference orientation (Allen and Tildesley, 2nd Edition, pages 106-111)
		! *******************************************************************************
		!  Particles are shaped like spherocylinders and thus have only two 
		!  rotational degrees of freedom. Orientation of solids of revolution can be
		!  specified according to a special case for linear molecules described in
		!  page 111 (Allen and Tildesley, 2nd Edition, 2017).
		! *******************************************************************************
		fixed_axis(1) = 0.d0
		fixed_axis(2) = 0.d0
		fixed_axis(3) = 1.d0

		! *******************************************************************************
		! Atom ID (Required in some visualization and analysis software)
		! *******************************************************************************
		atom = "C"

		! *******************************************************************************
		! Convert degrees to radians
		! *******************************************************************************
		quaternion_angle = quaternion_angle * pi / 180.d0

		! *******************************************************************************
		! Quaternion Algebra
		! *******************************************************************************
		!  See 'Quaternion algebras (2021)' book by John Voight.
		!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.
		! *******************************************************************************
		q(0,:) = dcos( quaternion_angle * 0.5d0 )			! Real part
		q(1,:) = dsin( quaternion_angle * 0.5d0 ) * fixed_axis(1)	! Imaginary part (Vector)
		q(2,:) = dsin( quaternion_angle * 0.5d0 ) * fixed_axis(2)	! Imaginary part (Vector)
		q(3,:) = dsin( quaternion_angle * 0.5d0 ) * fixed_axis(3)	! Imaginary part (Vector)

		! *******************************************************************************
		! Number of unit cells per axis (Simple Cube)
		! *******************************************************************************
		n_cells = dnint( dble( n_particles ) ** ( 1.d0 / 3.d0 ) )

		! *******************************************************************************
		! Unit cell length (Simple Cube)
		! *******************************************************************************
		!  In other words, the cube root of the simulation box volume over the 
		!  number of unit cells.
		! *******************************************************************************
		cell_length = ( 1.d0 / rho ) ** ( 1.d0 / 3.d0 )

		! *******************************************************************************
		! Unit cell length (Stretched Simple Cube)
		! *******************************************************************************
		!  Makes the unit cell proportional to the molecular geometry.
		!  REMINDER: all length quantities are divided by 'D' (diameter).
		!
		!  Spherocylinders: the unit cell is only stretched along z-axis.
		!
		!  For spherocylinders, L/D > 0, so the cubic structure will be stretched.
		!  For spheres, L/D = 0, so the cubic structure will remain unchanged.
		! *******************************************************************************
		coordinate(1) = cell_length				! x-axis (does not change)
		coordinate(2) = cell_length				! y-axis (does not change)
		coordinate(3) = cell_length * ( 1.d0 + aspect_ratio )	! z-axis (stretching)

		! *******************************************************************************
		! Simulation box length
		! *******************************************************************************
		box_length(:) = coordinate(:) * n_cells

		! *******************************************************************************
		! Positioning of particles (centers of mass)
		! *******************************************************************************
		counter = 1
		do i = 1, n_cells
			do j = 1, n_cells
				do k = 1, n_cells
					! *******************************************************
					! Particles on the right vertex of unit cell
					! *******************************************************
					r(1,counter) = dble( i - 1 ) * coordinate(1)
					r(2,counter) = dble( j - 1 ) * coordinate(2)
					r(3,counter) = dble( k - 1 ) * coordinate(3)
					counter = counter + 1
				end do
			end do
		end do

		! *******************************************************************************
		! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
		! *******************************************************************************
		do i = 1, n_particles
			r(:,i) = r(:,i) - ( 0.5d0 * box_length(:) )
		end do

	end subroutine config_sc

	! *************************************************************************************** !
	!    This subroutine allocates particles according to the BCC molecular configuration     !
	! *************************************************************************************** !
	subroutine config_bcc()

        implicit none

		integer*8	:: i, j, k, counter	! Counters

		allocate ( q (0:3,n_particles) )
		allocate ( r (3,n_particles) )
		allocate ( e (3,n_particles) )

		! *******************************************************************************
		! Unrotated reference orientation (Allen and Tildesley, 2nd Edition, pages 106-111)
		! *******************************************************************************
		!  Particles are shaped like spherocylinders and thus have only two 
		!  rotational degrees of freedom. Orientation of solids of revolution can be
		!  specified according to a special case for linear molecules described in
		!  page 111 (Allen and Tildesley, 2nd Edition, 2017).
		! *******************************************************************************
		fixed_axis(1) = 0.d0
		fixed_axis(2) = 0.d0
		fixed_axis(3) = 1.d0

		! *******************************************************************************
		! Atom ID (Required in some visualization and analysis software)
		! *******************************************************************************
		atom = "C"

		! *******************************************************************************
		! Convert degrees to radians
		! *******************************************************************************
		quaternion_angle = quaternion_angle * pi / 180.d0

		! *******************************************************************************
		! Quaternion Algebra
		! *******************************************************************************
		!  See 'Quaternion algebras (2021)' book by John Voight.
		!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.
		! *******************************************************************************
		q(0,:) = dcos( quaternion_angle * 0.5d0 )			! Real part
		q(1,:) = dsin( quaternion_angle * 0.5d0 ) * fixed_axis(1)	! Imaginary part (Vector)
		q(2,:) = dsin( quaternion_angle * 0.5d0 ) * fixed_axis(2)	! Imaginary part (Vector)
		q(3,:) = dsin( quaternion_angle * 0.5d0 ) * fixed_axis(3)	! Imaginary part (Vector)

		! *******************************************************************************
		! Number of unit cells per axis (Body-Centered Cube)
		! *******************************************************************************
		n_cells = dnint( ( 0.5d0 * dble( n_particles ) ) ** ( 1.d0 / 3.d0 ) )

		! *******************************************************************************
		! Unit cell length (Body-Centered Cube)
		! *******************************************************************************
		!  In other words, the cube root of the simulation box volume over the 
		!  number of unit cells.
		! *******************************************************************************
		cell_length = ( 2.d0 / rho ) ** ( 1.d0 / 3.d0 )

		! *******************************************************************************
		! Unit cell length (Stretched Body-Centered Cube)
		! *******************************************************************************
		!  Makes the unit cell proportional to the molecular geometry.
		!  REMINDER: all length quantities are divided by 'D' (diameter).
		!
		!  Spherocylinders: the unit cell is only stretched along z-axis.
		!
		!  For spherocylinders, L/D > 0, so the cubic structure will be stretched.
		!  For spheres, L/D = 0, so the cubic structure will remain unchanged.
		! *******************************************************************************
		coordinate(1) = cell_length				! x-axis (does not change)
		coordinate(2) = cell_length				! y-axis (does not change)
		coordinate(3) = cell_length * ( 1.d0 + aspect_ratio )	! z-axis (stretching)

		! *******************************************************************************
		! Simulation box length
		! *******************************************************************************
		box_length(:) = coordinate(:) * n_cells

		! *******************************************************************************
		! Positioning of particles (centers of mass)
		! *******************************************************************************
		counter = 1
		do i = 1, n_cells
			do j = 1, n_cells
				do k = 1, n_cells
					! *******************************************************
					! Particles on the right vertex of unit cell
					! *******************************************************
					r(1,counter) = dble( i - 1 ) * coordinate(1)
					r(2,counter) = dble( j - 1 ) * coordinate(2)
					r(3,counter) = dble( k - 1 ) * coordinate(3)
					counter = counter + 1
					! *******************************************************
					! Particles on the center of unit cell
					! *******************************************************
					r(1,counter) = ( dble( i ) - 0.5d0 ) * coordinate(1)
					r(2,counter) = ( dble( j ) - 0.5d0 ) * coordinate(2)
					r(3,counter) = ( dble( k ) - 0.5d0 ) * coordinate(3)
					counter = counter + 1
				end do
			end do
		end do

		! *******************************************************************************
		! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
		! *******************************************************************************
		do i = 1, n_particles
			r(:,i) = r(:,i) - ( 0.5d0 * box_length(:) )
		end do

	end subroutine config_bcc

	! *************************************************************************************** !
	!    This subroutine allocates particles according to the FCC molecular configuration     !
	! *************************************************************************************** !
	subroutine config_fcc()

        implicit none

		integer*8	:: i, j, k, counter	! Counters

		allocate ( q (0:3,n_particles) )
		allocate ( r (3,n_particles) )
		allocate ( e (3,n_particles) )

		! *******************************************************************************
		! Unrotated reference orientation (Allen and Tildesley, 2nd Edition, pages 106-111)
		! *******************************************************************************
		!  Particles are shaped like spherocylinders and thus have only two 
		!  rotational degrees of freedom. Orientation of solids of revolution can be
		!  specified according to a special case for linear molecules described in
		!  page 111 (Allen and Tildesley, 2nd Edition, 2017).
		! *******************************************************************************
		fixed_axis(1) = 0.d0
		fixed_axis(2) = 0.d0
		fixed_axis(3) = 1.d0

		! *******************************************************************************
		! Atom ID (Required in some visualization and analysis software)
		! *******************************************************************************
		atom = "C"

		! *******************************************************************************
		! Convert degrees to radians
		! *******************************************************************************
		quaternion_angle = quaternion_angle * pi / 180.d0

		! *******************************************************************************
		! Quaternion Algebra
		! *******************************************************************************
		!  See 'Quaternion algebras (2021)' book by John Voight.
		!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.
		! *******************************************************************************
		q(0,:) = dcos( quaternion_angle * 0.5d0 )			! Real part
		q(1,:) = dsin( quaternion_angle * 0.5d0 ) * fixed_axis(1)	! Imaginary part (Vector)
		q(2,:) = dsin( quaternion_angle * 0.5d0 ) * fixed_axis(2)	! Imaginary part (Vector)
		q(3,:) = dsin( quaternion_angle * 0.5d0 ) * fixed_axis(3)	! Imaginary part (Vector)

		! *******************************************************************************
		! Number of unit cells per axis (Face-Centered Cube)
		! *******************************************************************************
		n_cells = dnint( ( 0.25d0 * dble( n_particles ) ) ** ( 1.d0 / 3.d0 ) )

		! *******************************************************************************
		! Unit cell length (Face-Centered Cube)
		! *******************************************************************************
		!  In other words, the cube root of the simulation box volume over the 
		!  number of unit cells.
		! *******************************************************************************
		cell_length = ( 4.d0 / rho ) ** ( 1.d0 / 3.d0 )

		! *******************************************************************************
		! Unit cell length (Stretched Face-Centered Cube)
		! *******************************************************************************
		!  Makes the unit cell proportional to the molecular geometry.
		!  REMINDER: all length quantities are divided by 'D' (diameter).
		!
		!  Spherocylinders: the unit cell is only stretched along z-axis.
		!
		!  For spherocylinders, L/D > 0, so the cubic structure will be stretched.
		!  For spheres, L/D = 0, so the cubic structure will remain unchanged.
		! *******************************************************************************
		coordinate(1) = cell_length				! x-axis (does not change)
		coordinate(2) = cell_length				! y-axis (does not change)
		coordinate(3) = cell_length * ( 1.d0 + aspect_ratio )	! z-axis (stretching)

		! *******************************************************************************
		! Simulation box length
		! *******************************************************************************
		box_length(:) = coordinate(:) * n_cells

		! *******************************************************************************
		! Positioning of particles (centers of mass)
		! *******************************************************************************
		counter = 1
		do i = 1, n_cells
			do j = 1, n_cells
				do k = 1, n_cells
					! *******************************************************
					! Particles on the right vertex of unit cell
					! *******************************************************
					r(1,counter) = dble( i - 1 ) * coordinate(1)
					r(2,counter) = dble( j - 1 ) * coordinate(2)
					r(3,counter) = dble( k - 1 ) * coordinate(3)
					counter = counter + 1
					! *******************************************************
					! Particles on the front face of unit cell
					! *******************************************************
					r(1,counter) = dble( i - 1 ) * coordinate(1)
					r(2,counter) = ( dble( j ) - 0.5d0 ) * coordinate(2)
					r(3,counter) = ( dble( k ) - 0.5d0 ) * coordinate(3)
					counter = counter + 1
					! *******************************************************
					! Particles on the left face of unit cell
					! *******************************************************
					r(1,counter) = ( dble( i ) - 0.5d0 ) * coordinate(1)
					r(2,counter) = dble( j - 1 ) * coordinate(2)
					r(3,counter) = ( dble( k ) - 0.5d0 ) * coordinate(3)
					counter = counter + 1
					! *******************************************************
					! Particles on the lower face of unit cell
					! *******************************************************
					r(1,counter) = ( dble( i ) - 0.5d0 ) * coordinate(1)
					r(2,counter) = ( dble( j ) - 0.5d0 ) * coordinate(2)
					r(3,counter) = dble( k - 1 ) * coordinate(3)
					counter = counter + 1
				end do
			end do
		end do

		! *******************************************************************************
		! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
		! *******************************************************************************
		do i = 1, n_particles
			r(:,i) = r(:,i) - ( 0.5d0 * box_length(:) )
		end do

	end subroutine config_fcc

	! *************************************************************************************** !
	!      This subroutine creates a file with all particles' positions and orientations      !
	! *************************************************************************************** !
	subroutine config_out()

	implicit none

		integer*8	:: i	! Counter

		! *******************************************************************************
		! Simple cubic structure
		! *******************************************************************************
		if ( config_selec(1) ) then
			open  ( UNIT = 10, FILE = "Initial_Configuration/OVITO/"//trim(descriptor_date)//"/initconf_ε" &
			&			  //trim(descriptor_file1)//"_ρ"//trim(descriptor_file2)//"_sc.xyz" )
			write ( 10, "(I4)" ) n_particles
			write ( 10, * ) " "
			do i = 1, n_particles
				! ***************************************************************
				! Initial configuration for OVITO (reduced units).
				! ***************************************************************
				!  The last two parameters are the spherocylindrical diameter 
				!  along x- and y-axes, and the spherocylindrical length along
				!  z-axis.
				! ***************************************************************
				write ( 10, * ) atom, r(1,i), r(2,i), r(3,i), q(0,i), q(1,i), q(2,i), q(3,i), 0.5d0, aspect_ratio
			end do
			close ( 10 )
		! *******************************************************************************
		! Body-centered cubic structure
		! *******************************************************************************
		else if ( config_selec(2) ) then
			open  ( UNIT = 10, FILE = "Initial_Configuration/OVITO/"//trim(descriptor_date)//"/initconf_ε" &
			&			  //trim(descriptor_file1)//"_ρ"//trim(descriptor_file2)//"_bcc.xyz" )
			write ( 10, "(I4)" ) n_particles
			write ( 10, * ) " "
			do i = 1, n_particles
				! ***************************************************************
				! Initial configuration for OVITO (reduced units).
				! ***************************************************************
				!  The last two parameters are the spherocylindrical diameter 
				!  along x- and y-axes, and the spherocylindrical length along
				!  z-axis.
				! ***************************************************************
				write ( 10, * ) atom, r(1,i), r(2,i), r(3,i), q(0,i), q(1,i), q(2,i), q(3,i), 0.5d0, aspect_ratio
			end do
			close ( 10 )
		! *******************************************************************************
		! Face-centered cubic structure
		! *******************************************************************************
		else if ( config_selec(3) ) then
			open  ( UNIT = 10, FILE = "Initial_Configuration/OVITO/"//trim(descriptor_date)//"/initconf_ε" &
			&			  //trim(descriptor_file1)//"_ρ"//trim(descriptor_file2)//"_fcc.xyz" )
			write ( 10, "(I4)" ) n_particles
			write ( 10, * ) " "
			do i = 1, n_particles
				! ***************************************************************
				! Initial configuration for OVITO (reduced units).
				! ***************************************************************
				!  The last two parameters are the spherocylindrical diameter 
				!  along x- and y-axes, and the spherocylindrical length along
				!  z-axis.
				! ***************************************************************
				write ( 10, * ) atom, r(1,i), r(2,i), r(3,i), q(0,i), q(1,i), q(2,i), q(3,i), 0.5d0, aspect_ratio
			end do
			close ( 10 )
		end if

	end subroutine config_out

end module initconfig
