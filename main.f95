! ############################################################################################### !
!                 Canonical Monte Carlo algorithm for spherocylindrical molecules                 !
!             The Vega-Lago Shortest Distance (T. Kihara, 1951) model is used to trial            !
!          moves, and establishes the molecular configurations (position and orientation)         !
!                                    of the Reference System.                                     !
!       The spherical square-well and Kihara potentials are used to model the intermolecular      !
!                              interactions in the Perturbed System.                              !
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
! Main References:                           T. Kihara                                            !
!                                    J. Phys. Soc. Japan (1951)                                   !
!                                     DOI: 10.1143/jpsj.6.289                                     !
!                             --------------------------------------                              !
!                                   M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                                 C. R. A. Abreu, F. A. Escobedo                                  !
!                                J. Chem. Phys. 124, 054116 (2006)                                !
!                                     DOI: 10.1063/1.2165188                                      !
!                             --------------------------------------                              !
!                                           O. K. Smith                                           !
!                           Communications of the ACM, 4(4), 168 (1961)                           !
!                                    DOI: 10.1145/355578.366316                                   !
!                             --------------------------------------                              !
!                                  J. T. Lopes, L. F. M. Franco                                   !
!                            Ind. Eng. Chem. Res. 58, 6850−6859 (2019)                            !
!                                  DOI: 10.1021/acs.iecr.9b00766                                  !
!                             --------------------------------------                              !
!                                      N. Metropolis et al.                                       !
!                                 J. Chem. Phys. 21, 1087 (1953)                                  !
!                                     DOI: 10.1063/1.1699114                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

program main

	! Uses four modules: global variables, variable initialization, initial configuration, and
	!                    directory creator
	use globalvar
	use initvar
	use initconfig
	use folders

	implicit none

	! ***************************************************************************************
	! Integer variables
	! ***************************************************************************************
	integer*8		:: i, j, k	! Loop counters
	integer*8		:: pot_param	! Auxiliar

	! ***************************************************************************************
	! Real variables/functions
	! ***************************************************************************************
	real*8			:: sigmavl	! Vega-Lago shortest distance (function)
	real*8			:: rvl	        ! Vega-Lago shortest distance (variable)
	real*8			:: rijsq	! Vector distance between particles i and j (squared)
	real*8			:: max_boxl	! Largest box length
	real*8			:: max_aspec	! Largest dimension of a spherocylindrical geometry
	real*8			:: msgaux	! Auxiliar
	real*8			:: cutoff	! Cutoff diameter
	real*8			:: cutoffsq	! Cutoff diameter (squared)
	real*8, dimension (3)	:: ei, ej	! Orientation of particles i and j
	real*8, dimension (3)	:: rij		! Vector distance between particles i and j
	real*8, dimension (3)	:: em, en	! Orientation (before/after a trial move)
	real*8, dimension (3)	:: rm, rn	! Position (before/after a trial move)
	real*8, dimension (0:3)	:: qm, qn	! Quaternion (before/after a trial move)

	! ***************************************************************************************
	! Real variables (allocatable)
	! ***************************************************************************************
	real*8, dimension (:), allocatable	:: vm, vn	! Perturbed potential energy (before/after a trial move)
	real*8, dimension (:), allocatable	:: dv		! Perturbed potential energy difference between two microstates m and n
	real*8, dimension (:), allocatable 	:: a1, a2	! First- and second-order TPT coefficients
	real*8, dimension (:), allocatable	:: apert	! Helmholtz free energy of the perturbed system
	real*8, dimension (:), allocatable	:: dpa1, dpa2	! First- and second-order TPT coefficients (standard deviation)
	real*8, dimension (:), allocatable	:: dpapert	! Helmholtz free energy of the perturbed system (standard deviation)

	! ***************************************************************************************
	! Character variables
	! ***************************************************************************************
	character		:: l_format*32	! Number of attractive/repulsive range parameters (format)

	! ***************************************************************************************
	! Logical variables
	! ***************************************************************************************
	logical			:: overlap	! Detects overlap between two particles (VL)   : TRUE = overlap detected;  FALSE = overlap not detected
	logical			:: flag		! Generic true/false flag

	! ***************************************************************************************
	! Random number generator seed (kept fixed for repeatability)
	! ***************************************************************************************
	seed = 123456

	! ***************************************************************************************
	! Molecular configuration selection (see 'initconfig' module)
	! ***************************************************************************************
	call config_selection()

	! ***************************************************************************************
	! Initialization of common variables (see 'initvar' module)
	! ***************************************************************************************
	call common_var()

	if ( ff_selec(1) ) then
		allocate ( a1(n_lambda) )
		allocate ( a2(n_lambda) )
		allocate ( apert(n_lambda) )
		allocate ( dpa1(n_lambda) )
		allocate ( dpa2(n_lambda) )
		allocate ( dpapert(n_lambda) )
		allocate ( v(n_lambda) )
		allocate ( vmc(n_lambda) )
		allocate ( vn(n_lambda) )
		allocate ( vm(n_lambda) )
		allocate ( swrange(n_lambda) )
		allocate ( dv(n_lambda) )
		allocate ( lfexist(n_lambda) )
	else if ( ff_selec(2) ) then
		allocate ( a1(n_n) )
		allocate ( a2(n_n) )
		allocate ( apert(n_n) )
		allocate ( dpa1(n_n) )
		allocate ( dpa2(n_n) )
		allocate ( dpapert(n_n) )
		allocate ( v(n_n) )
		allocate ( vmc(n_n) )
		allocate ( vn(n_n) )
		allocate ( vm(n_n) )
		allocate ( dv(n_lambda) )
		allocate ( lfexist(n_lambda) )
	end if
	allocate ( qmc (0:3,n_particles) )
	allocate ( rmc (3,n_particles) )
	allocate ( emc (3,n_particles) )

	! ***************************************************************************************
	! CPU Clock
	! ***************************************************************************************
	call date_and_time( VALUES = date_time )

	! ***************************************************************************************
	! Initial configuration (see 'initconfig' module)
	! ***************************************************************************************
	! Calls 'config_sc' subroutine if the user chooses a simple cubic structure
	if ( config_selec(1) ) then
		call config_sc()
	! Calls 'config_bcc' subroutine if the user chooses a body-centered cubic structure
	else if ( config_selec(2) ) then
		call config_bcc()
	! Calls 'config_fcc' subroutine if the user chooses a face-centered cubic structure
	else if ( config_selec(3) ) then
		call config_fcc()
	end if

	! ***************************************************************************************
	! Initialization of Monte Carlo parameters (see 'initvar' module)
	! ***************************************************************************************
	call montecarlo_var()

	! ***************************************************************************************
	! Initialization of Block-averaging parameters (see 'initvar' module)
	! ***************************************************************************************
	call blockav_var()

	! ***************************************************************************************
	! Initialization of Inquiry/Control variables (see 'initvar' module)
	! ***************************************************************************************
	call inquery_var()

	! ***************************************************************************************
	! Friendly Summary
	! ***************************************************************************************
	print *, " "
	print *, "Data read successfully!         "
	print *, "Number of particles:            ", n_particles
	print *, "Reduced number density (κρ*):   ", rho
	print *, "Reduced temperature:            ", temp
	print *, "Aspect ratio (L/D):             ", aspect_ratio
	print *, "Quaternion angle:               ", quaternion_angle
	if ( ff_selec(1) ) then
		print *, "No. of attractive parameters:   ", n_lambda
	else if ( ff_selec(2) ) then
		print *, "No. of repulsive parameters:    ", n_n
	end if
	print *, "Total number of cycles:         ", max_cycles
	print *, "Equilibration cycles:           ", n_equil
	print *, "Production cycles:              ", max_cycles-n_equil
	print *, "Adjustment frequency (cycles):  ", n_adjust
	print *, "Data saving frequency (cycles): ", n_save
	if ( coef_check ) then
		print *, "Minimum no. of blocks:          ", min_blocks
		print *, "Maximum no. of blocks:          ", max_blocks
	end if
	print *, "Write out trajectory?           ", traj_inq
	print *, "Calculate TPT coefficients?     ", coef_inq
	print *, " "
	print *, "Write out only production-related potentials? ", pot_inq
	print *, " "
	if ( ff_selec(1) ) then
		do i = 1, n_lambda
			write ( *, "(A28,I2.2,A2,F5.2)" ) "Attractive range parameter #", i, ": ", lambda(i)
		end do
	else if ( ff_selec(2) ) then
		do i = 1, n_n
			write ( *, "(A28,I2.2,A2,F5.2)" ) "Repulsive index parameter #", i, ":  ", n_repulsive(i)
		end do
	end if
	print *, " "
	! Wait for user
	print *, "Press and enter any KEY to continue..."
	read ( *, * ) dummy

	! Status
	print *, " "
	print *, "Creating folders and files..."
	print *, " "

	! ***************************************************************************************
	! Initial configuration folder (see 'folders' module)
	! ***************************************************************************************
	call initfolder()

	! ***************************************************************************************
	! Create output directories (see 'folders' module)
	! ***************************************************************************************
	call inquire_folders()

	! ***************************************************************************************
	! Create date subfolders (see 'folders' module)
	! ***************************************************************************************
	call date_folders()

	! ***************************************************************************************
	! Create attractive parameter subfolders (see 'folders' module)
	! ***************************************************************************************
	call lambda_folders()

	! ***************************************************************************************
	! Initial configuration file (see 'initconfig' module)
	! ***************************************************************************************
	call config_out()

	! ***************************************************************************************
	! Output file units
	! ***************************************************************************************

	! Trajectory file (depends on user's choice)
	if ( traj_check ) then
		open ( UNIT = 20, FILE = "Trajectories/"//trim(descriptor_date)//"/traj_ε"//trim(descriptor_file1)//"_ρ" &
		&			  //trim(descriptor_file2)//".xyz" )
		write ( 20, "(I4)" ) n_particles
		write ( 20, * ) " "
		do i = 1, n_particles
			write ( 20, * ) atom, r(1,i), r(2,i), r(3,i), q(0,i), q(1,i), q(2,i), q(3,i), &
			&		0.5d0, aspect_ratio
		end do
	end if

	! Ratio file (translation)
	open ( UNIT = 30, FILE = "Ratio/Translation/"//trim(descriptor_date)//"/ratio_ε"//trim(descriptor_file1)//"_ρ" &
	&			  //trim(descriptor_file2)//".dat" )

	! Ratio file (rotation)
	open ( UNIT = 40, FILE = "Ratio/Rotation/"//trim(descriptor_date)//"/ratio_ε"//trim(descriptor_file1)//"_ρ" &
	&			  //trim(descriptor_file2)//".dat" )

	! Order parameter file
	open ( UNIT = 50, FILE = "Order_Parameter/"//trim(descriptor_date)//"/order_ε"//trim(descriptor_file1)//"_ρ" &
	&			  //trim(descriptor_file2)//".dat" )

	! Potential parameters
	if ( ff_selec(1) ) then
		pot_param = n_lambda
		do counter_lambda = 1, n_lambda
			! File descriptor: attractive range parameter
			write ( descriptor_lamb, format_lamb ) lambda(counter_lambda)
			! Potential files
			open ( UNIT = ( 60 + counter_lambda ), FILE = "Potential/SW/"//trim(descriptor_date)//"/Lambda_" &
			&				       //trim(descriptor_lamb)//"/thermo_ε"//trim(descriptor_file1)//"_ρ" &
			&				       //trim(descriptor_file2)//".dat" )
		end do
	else if ( ff_selec(2) ) then
		pot_param = n_n
		do counter_n = 1, n_n
			! File descriptor: repulsive parameter
			write ( descriptor_n, format_n ) n_repulsive(counter_n)
			! Potential files
			open ( UNIT = ( 60 + counter_n ), FILE = "Potential/KH/"//trim(descriptor_date)//"/n_" &
			&				  //trim(descriptor_n)//"/thermo_ε"//trim(descriptor_file1)//"_ρ" &
			&				  //trim(descriptor_file2)//".dat" )
		end do
	end if

	! Status
	print *, "Done!"

	! ***************************************************************************************
	! Largest box length
	! ***************************************************************************************
	max_boxl = 0.d0
	do i = 1, 3
		if ( max_boxl <= box_length(i) ) then
			max_boxl = box_length(i)
		end if
	end do

	! ***************************************************************************************
	! Largest spherocylindrical dimension
	! ***************************************************************************************
	max_aspec = ( 1.d0 + aspect_ratio )

	! ***************************************************************************************
	! Hard-core volumetric relation (spherocylinders and spheres)
	! ***************************************************************************************
	!  The spherical and spherocylindrical volumes are related as follows:
	!
	!                          V_sphere = V_sc
	!      (π/6)*(sigma_sphere**3) = (π/6)*(D**3)+(π/4)*(D**2)*(L)
	!          sigma_sphere = {(D**3)+[(3/2)*(D**2)*(L)]}**(1/3)
	!
	!  The formula above is not in reduced units. The fundamental unit of length is 'D',
	!  thus all length quantities are divided by 'D'. Then:
	!
	!           sigma_sphere = {(1**3)+[(3/2)*(1**2)*(L/D)]}**(1/3)
	!                 sigma_sphere = {1+[(3/2)*(L/D)]}**(1/3)
	!                   sigma_sphere = [1+1.5*(L/D)]**(1/3)
	!
	! ***************************************************************************************
	sigsphere  = ( 1.d0 + ( 1.5d0 * aspect_ratio ) ) ** ( 1.d0 / 3.d0 )

	! ***************************************************************************************
	! Effective range of attraction
	! ***************************************************************************************
	swrange(:) = lambda(:) * sigsphere

	! ***************************************************************************************
	! Start simulation timer
	! ***************************************************************************************
	call cpu_time(start_timer)

	! Status
	print *, " "
	print *, "Verifying initial configuration..."
	print *, " "

	! ***************************************************************************************
	! Active transformation
	! ***************************************************************************************
	!  Converts the unit quaternion into an active rotation/transformation in the 3D Euclidean space
	!  See 'subroutines' code for more details.
	! ***************************************************************************************
	do i = 1, n_particles
		call active_transformation(fixed_axis,q(:,i),e(:,i))
	end do

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
	! Overlap check (initial configuration)
	! ***************************************************************************************
	!  Searches for possible molecular overlaps in the initial configuration
	! ***************************************************************************************

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
			rij(:)  = r(:,i) - r(:,j)
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
			! Cutoff distance
			! ***********************************************************************
			if ( rijsq <= cutoffsq ) then
				! ***********************************************************************
				! Vega-Lago Overlap Criterion - Shortest distance
				!  See 'functions' file for more information
				! ***********************************************************************
				rvl = sigmavl(ei,ej,rij,rijsq)
				if ( rvl <= 1.d0 ) then
					write ( *, 100 ) i, j, rho
					! Stop program if overlap is detected
					call exit()
				end if
			end if

		end do

	end do

	! Might be necessary to change if the number of decimal places is higher than 5.
 100	format ( "Overlap detected in initial configuration between particles ", I4.4, " and ", I4.4, &
	&	 " for a reduced number density of ", F0.5, "! Exiting..." )

	! Status
	print *, "No overlaps detected. Resuming..."

	! ***************************************************************************************
	! Monte Carlo Simulation
	! ***************************************************************************************
	print *, " "
	print *, "Now running Monte Carlo simulation. Computing total energy..."
	print *, " "

	! ***************************************************************************************
	! Computation of total potential energy (initial configuration)
	! See 'subroutines' code for more details.
	! ***************************************************************************************
	call compute_total_energy()

	! Status
	print *, "Done!"

	! ***************************************************************************************
	! Monte Carlo parameters
	! ***************************************************************************************
	!  Translation and rotation are independent moves.
	! ***************************************************************************************
	resetmc   = .false.	! Monte Carlo parameters resetting		(initial value)
	mov_trans = .false.	! Translational move selector		        (initial value)
	mov_rot   = .false.	! Rotational move selector			(initial value)
	stop_t    = .false.	! Translational ratio threshold modifier	(initial value)
	stop_r    = .false.	! Rotational ratio threshold modifier		(initial value)
	drmax     = max_trans	! Maximum translational displacement		(initial value)
	angmax    = max_rot	! Maximum rotational displacement		(initial value)
	reset_c   = 0		! Monte Carlo parameter resets counter		(initial value)
	nacct     = 0		! Translational move acceptance counter  	(initial value)
	naccr     = 0		! Rotational move acceptance counter	        (initial value)
	movt      = 0		! Translational move counter			(initial value)
	movr      = 0		! Rotational move counter			(initial value)
	r_acc_r   = 0.5d0	! Translational acceptance ratio threshold	(initial value)
	r_acc_t   = 0.5d0 	! Rotational acceptance ratio threshold		(initial value)
	qmc(:,:)  = q(:,:)	! Quaternion algebra				(initial value)
	rmc(:,:)  = r(:,:)	! Position of particles				(initial value)
	emc(:,:)  = e(:,:)	! Orientation of particles			(initial value)
	vmc(:)    = v(:)	! Total potential energy			(initial value)

	! ***************************************************************************************
	! Metropolis Non-sequential Algorithm - Importance Sampling
	! (see Metropolis et al., J. Chem. Phys. 21, 1087 (1953), for more information)
	! ***************************************************************************************
	print *, " "
	print *, "Running Metropolis algorithm..."
	print *, " "

	! ***************************************************************************************
	! Messages of parameter resetting
	! ***************************************************************************************
 200	format ( "Acceptance ratio threshold decreased at ", I10.10, " th cycle.", /, "Translational displacement reached ", &
 	&	 F7.3, "!", /, "Translational threshold changed from ", F6.4, " to ", F6.4, ".", /, &
	&	 "Resetting Monte Carlo parameters...", / )

 250	format ( "Acceptance ratio threshold increased at ", I10.10, " th cycle.", /, "Translational displacement reached ", &
 	&	 F7.3, "!", /, "Translational threshold changed from ", F6.4, " to ", F6.4, ".", /, &
	&	 "Resetting Monte Carlo parameters...", / )

 300	format ( "Acceptance ratio threshold decreased at ", I10.10, " th cycle.", /, "Rotational displacement reached ", &
 	&	 F5.3, "!", /, "Rotational threshold changed from ", F6.4, " to ", F6.4, ".", /, &
	&	 "Resetting Monte Carlo parameters...", / )

 350	format ( "Acceptance ratio threshold increased at ", I10.10, " th cycle.", /, "Rotational displacement reached ", &
 	&	 F5.3, "!", /, "Rotational threshold changed from ", F6.4, " to ", F6.4, ".", /, &
	&	 "Resetting Monte Carlo parameters...", / )

	! ***************************************************************************************
	! Simulation cycles
	! ***************************************************************************************
	!  A 'cycle' is characterized by a trial move (rotation or translation) of a randomized
	!  particle i.
	! ***************************************************************************************
	do cycles = 1, max_cycles

		! *******************************************************************************
		! Potential data
		! *******************************************************************************
		! Save all potential data (equilibration and production)
		if ( .not. potential_check ) then
			if ( mod( cycles, n_save ) == 0 ) then
			do i = 1, pot_param
				write ( 60 + i, *) cycles, vmc(i)
			end do
		end if
		!Save potential data (production-related only)
		else if ( potential_check ) then
			if ( cycles > n_equil .and. mod( cycles, n_save ) == 0 ) then
				do i = 1, pot_param
					write ( 60 + i, * ) cycles, vmc(i)
				end do
			end if
		end if

		! *******************************************************************************
		! Random Move Selection
		! *******************************************************************************
		!  Both moves (rotation and translation) have the same likelihood of being
		!  chosen. The variables 'movt' and 'movr' represent how many times translation
		!  and rotation are randomly chosen, respectively. These variables will be
		!  required during displacement adjustment.
		! *******************************************************************************

		if ( aspect_ratio /= 0.d0 ) then
			! ***********************************************************************
			! Pseudorandom number generator (uniform distribution)
			! (see 'subroutines' code for more details)
			! ***********************************************************************
			call ranf()
			! ***********************************************************************
			! Translation criterion
			! ***********************************************************************
			if ( random_n < 0.5d0 ) then
				mov_trans = .true.	! Enable translation
				mov_rot   = .false.	! Disable rotation
				movt      = movt + 1	! Increment move counter
			! ***********************************************************************
			! Rotation criterion
			! ***********************************************************************
			else if ( random_n >= 0.5d0 ) then
				mov_rot   = .true.	! Enable rotation
				mov_trans = .false.	! Disable translation
				movr      = movr + 1	! Increment move counter
			end if
		! *******************************************************************************
		! Spherical case (only translational DOFs)
		! *******************************************************************************
		!  Ignores random movement selection if the length-to-diameter aspect ratio is 0.
		! *******************************************************************************
		else
			mov_trans = .true.		! Enable translation (permanently)
			mov_rot	  = .false.		! Disable rotation (permanently)
			movt	  = movt+1		! Increment move counter
		end if

		! *******************************************************************************
		! Random selection of particles of i-index
		! *******************************************************************************
		!  Generates a random number and determines the i-index of the particle using the
		!  INT function. See <https://gcc.gnu.org/onlinedocs/gfortran/INT.html> for more 
		!  information.
		! *******************************************************************************
		call ranf()
		i = int( random_n * dble( n_particles ) ) + 1

		! *******************************************************************************
		! Assignment of previous configuration (Microstate m)
		! *******************************************************************************
		!  Assigns the previous configuration (position, quaternion, and orientation) of a 
		!  randomly selected particle i to a microstate m. A microstate m simply means 
		!  'a state of the system before a trial move'. The previous configuration may be 
		!  either the initial configuration or the lastest accepted trial move.
		! *******************************************************************************
		rm(:) = rmc(:,i)	! Position
		qm(:) = qmc(:,i)	! Quaternion
		em(:) = emc(:,i)	! Orientation

		! *******************************************************************************
		! Computation of potential energy of particle i (Microstate m)
		! See 'subroutines' code for more details.
		! *******************************************************************************
		call compute_particle_energy(i,em,rm,vm)

		! *******************************************************************************
		! Random Displacement (trial move of particle i)
		! *******************************************************************************
		!  Both moves, rotation and translation, are independent moves, which means
		!  a random particle i performs either a translational or rotational displacement
		!  at a time.
		! *******************************************************************************
		!  See Allen and Tildesley, 2nd Edition (2017), pages 169-173 for more information
		!  on random displacements.
		! *******************************************************************************

		! *******************************************************************************
		! TRANSLATION
		! *******************************************************************************
		if ( mov_trans ) then

			! ***********************************************************************
			! Random translation along x-axis
			! ***********************************************************************
			call ranf()
			rn(1) = rm(1) + ( ( 2.d0 * random_n ) - 1.d0 ) * drmax	! Range [-drmax,drmax]
			! ***********************************************************************
			! Random translation along y-axis
			! ***********************************************************************
			call ranf()
			rn(2) = rm(2) + ( ( 2.d0 * random_n ) - 1.d0 ) * drmax	! Range [-drmax,drmax]
			! ***********************************************************************
			! Random translation along z-axis
			! ***********************************************************************
			call ranf()
			rn(3) = rm(3) + ( ( 2.d0 * random_n ) - 1.d0 ) * drmax	! Range [-drmax,drmax]
			! ***********************************************************************
			! Minimum Image Convention
			! (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
			! ***********************************************************************
			rn(:) = rn(:) - ( box_length(:) ) * anint( rn(:) / box_length(:) )

		! *******************************************************************************
		! NO TRANSLATION
		! *******************************************************************************
		!  If translation is currently disabled for this cycle, assign all positions of the
		!  microstate m (before trial move) to the microstate n.
		! *******************************************************************************
		else if ( .not. mov_trans ) then
			rn(:) = rm(:)
		end if

		! *******************************************************************************
		! ROTATION
		! *******************************************************************************
		if ( mov_rot ) then

			! ***********************************************************************
			! Random Composed Unit Quaternion
			! (see 'subroutines' code for more details)
			! ***********************************************************************
			call composed_quaternion(qm,qn)

			! ***********************************************************************
			! Active transformation
			! ***********************************************************************
			!  Converts the unit quaternion into an active rotation/transformation
			!  in the 3D Euclidean space.
			! ***********************************************************************
			call active_transformation(fixed_axis,qn,en)

		! *******************************************************************************
		! NO ROTATION
		! *******************************************************************************
		!  If rotation is currently disabled for this cycle, assign all orientations and
		!  quaternions of the microstate m (before trial move) to the microstate n.
		! *******************************************************************************
		else if ( .not. mov_rot ) then
			qn(:) = qm(:)
			en(:) = em(:)
		end if

		! *******************************************************************************
		! Metropolis criterion (Boltzmann factor of the energy difference)
		! *******************************************************************************
		!  The Metropolis criterion establishes whether a trial move will be accepted or
		!  rejected, which is decided based on the potential energy difference of microstates 
		!  n and m, that is, δUₙₘ = (Uₙ - Uₘ).
		!
		!  If δUₙₘ < 0, then the trial move is accepted.
		!  If δUₙₘ ≥ 0, then we calculate exp(-βδUₙₘ) and compare the result to a random number:
		!
		!  	If exp(-βδUₙₘ) > random_number, then the trial move is accepted.
		!	If exp(-βδUₙₘ) ≤ random_number, then the trial move is rejected.
		!
		!  In case of hard-core fluids, δUₙₘ can only result in 0 (no overlap) or ∞ (overlap),
		!  which falls in the second condition above. Therefore, exp(-βδUₙₘ) can only assume 
		!  1 (maximum value) or 0 (minimum value). 1 will be greater than any random_number
		!  generated from a uniform distribution in the range of [0,1[, and 0 will be less 
		!  than or equal to any random_number.
		!
		!  In other words, for hard-core fluids, if there is no overlap of particles, then the
		!  trial move is immediately accepted. Otherwise, it is rejected.
		!
		!  See Allen and Tildesley, 2nd Edition (2017), pages 155-160 for more information
		! *******************************************************************************

		! *******************************************************************************
		! Overlap Verification (VL)
		! See 'subroutines' code for more details.
		! *******************************************************************************
		!  Checks if a random displacement (translation or rotation) of particle i causes 
		!  any overlaps with other particles.
		! *******************************************************************************
		call check_overlap(i,en,rn,overlap)

		! *******************************************************************************
		! Move Accepted
		! *******************************************************************************
		if ( .not. overlap ) then

			! ***********************************************************************
			! System configuration update
			! ***********************************************************************
			!  Assigns the configuration (position, quaternion, and orientation) of 
			!  a trial move of particle i to the system configuration.
			! ***********************************************************************
			rmc(:,i) = rn(:)	! Update position
			qmc(:,i) = qn(:)	! Update quaternion
			emc(:,i) = en(:)	! Update orientation

			! ***********************************************************************
			! Computation of potential energy of particle i (Microstate n)
			! ***********************************************************************
			call compute_particle_energy(i,en,rn,vn)

			! ***********************************************************************
			! Computation of energy difference of microstates n and m
			! ***********************************************************************
			dv(:)  = vn(:) - vm(:)

			! ***********************************************************************
			! System energy update
			! ***********************************************************************
			vmc(:) = vmc(:) + dv(:)

			! ***********************************************************************
			! Displacement counter update
			! ***********************************************************************
			!  Only updates if the move has been previously selected from a
			!  uniform distribution.
			! ***********************************************************************
			if ( mov_trans ) then
				nacct = nacct + 1	! Translational move counter
			else if ( mov_rot ) then
				naccr = naccr + 1	! Rotational move counter
			end if

		end if

		! *******************************************************************************
		! Adjustment of maximum displacement
		! *******************************************************************************
		if ( cycles <= n_equil ) then	! During equilibration only

			! ***********************************************************************
			! Adjustment of maximum translation
			! ***********************************************************************
			if ( mov_trans ) then

				! ***************************************************************
				! Adjustment frequency (translation)
				! ***************************************************************
				if ( mod( movt, n_adjust ) == 0 ) then

					! *******************************************************
					! Acceptance ratio 
					! (non-overlapping microstates over sampled microstates)
					! *******************************************************
					ratio = dble( nacct ) / dble( n_adjust )
					! *******************************************************
					! Translational adjustment
					! *******************************************************
					!  Decreases displacement if the acceptance ratio is 
					!  lower than or equal to a dynamic threshold.
					!  Increases displacement if the acceptance ratio is
					!  greater than a dynamic threshold.
					! 
					!  The adjustment is made to improve the sampling of
					!  microstates.
					! *******************************************************
					if ( ratio <= r_acc_t ) then
						drmax = 0.95d0 * drmax
					else
						drmax = 1.05d0 * drmax
					end if

					! *******************************************************
					! Dynamic threshold (translation) - Decrease value
					! *******************************************************
					if ( ( drmax < 1.d-5 ) .and. ( .not. stop_t ) ) then		! Arbitrary condition
						msgaux  = r_acc_t
						r_acc_t = r_acc_t * 0.95d0				! Arbitrary modification
						! Minimum threshold criterion (translation)
						if ( r_acc_t < ( 1.d0 / dble( n_adjust ) ) ) then
							r_acc_t  = 1.d0 / dble( n_adjust )
							stop_t   = .true.
						end if
						! Reset Monte Carlo parameters
						resetmc = .true.
						! Message
						write ( *, 200 ) cycles, drmax, msgaux, r_acc_t
						! Increment counter
						reset_c = reset_c + 1
					end if

					! *******************************************************
					! Dynamic threshold (translation) - Increase value
					! *******************************************************
					if ( ( drmax >= max_boxl ) .and. ( .not. stop_t ) ) then	! Arbitrary condition
						msgaux  = r_acc_t
						r_acc_t = r_acc_t * 1.05d0				! Arbitrary modification
						! Maximum threshold (translation)
						if ( r_acc_t > ( dble( n_adjust - 1 ) / dble( n_adjust ) ) ) then
							r_acc_t  = dble( n_adjust - 1 ) / dble( n_adjust )
							stop_t   = .true.
						end if
						! Reset Monte Carlo parameters
						resetmc = .true.
						! Message
						write ( *, 250 ) cycles, drmax, msgaux, r_acc_t
						! Increment counter
						reset_c = reset_c + 1
					end if

					! *******************************************************
					! Ratio data
					! *******************************************************
					write( 30, * ) cycles, ratio, drmax, r_acc_t

					! Reset counter
					nacct = 0

					! *******************************************************
					! Reset Monte Carlo parameters
					! See 'subroutines' code for more details.
					! *******************************************************
					!  OBS.: These parameters were not resetted in the 
					!  Monte Carlo simulations performed in our paper.
					! *******************************************************
					!if ( resetmc ) then
					!	call reset_system()
					!	! Advance loop
					!	cycle
					!end if

				end if

			end if

			! ***********************************************************************
			! Adjustment of maximum rotation
			! ***********************************************************************
			if ( mov_rot ) then

				! ***************************************************************
				! Adjustment frequency (rotation)
				! ***************************************************************
				if ( mod( movr, n_adjust ) == 0 ) then

					! *******************************************************
					! Acceptance ratio 
					! (non-overlapping microstates over sampled microstates)
					! *******************************************************
					ratio = dble( naccr ) / dble( n_adjust )
					! *******************************************************
					! Rotational adjustment
					! *******************************************************
					!  Decreases displacement if the acceptance ratio is 
					!  lower than or equal to a dynamic threshold.
					!  Increases displacement if the acceptance ratio is
					!  greater than a dynamic threshold.
					!
					!  The adjustment is made to improve the sampling of
					!  microstates.
					! *******************************************************
					if ( ratio <= r_acc_r ) then
						angmax = 0.95d0 * angmax
					else
						angmax = 1.05d0 * angmax
					end if

					! *******************************************************
					! Dynamic threshold (rotation) - Decrease value
					! *******************************************************
					if ( ( ( angmax * ( 0.5d0 * max_elong ) ) < 1.d-4 ) &		! Arbitrary condition
					&   .and. ( .not. stop_r ) ) then
						msgaux  = r_acc_r
						r_acc_r = 0.95d0 * r_acc_r				! Arbitrary modification
						! Minimum threshold (rotation)
						if ( r_acc_t < ( 1.d0 / dble( n_adjust ) ) ) then
							r_acc_t  = 1.d0 / dble( n_adjust )
							stop_r   = .true.
						end if
						! Reset Monte Carlo parameters
						resetmc = .true.
						! Message
						write ( *, 300 ) cycles, angmax, msgaux, r_acc_r
						! Increment counter
						reset_c = reset_c + 1
					end if

					! *******************************************************
					! Dynamic threshold (rotation) - Increase value
					! *******************************************************
					!  If θ = |s| > π, then the orientation will be identical
					!  the wrapped orientation s - 2πŝ, where s is a random
					!  vector (Allen & Tildesley, 2017).
					!
					!  See Allen and Tildesley, 2nd Edition (2017), page 172 
					!  for more information.
					! *******************************************************
					if ( ( angmax > pi ) .and. ( .not. stop_r ) ) then		! Half-turn condition
						msgaux  = r_acc_r
						r_acc_r = 1.05d0 * r_acc_r				! Arbitrary modification
						! Maximum threshold (rotation)
						if ( r_acc_r > ( dble( n_adjust - 1 ) / dble( n_adjust ) ) ) then
							r_acc_r  = dble( n_adjust - 1 ) / dble( n_adjust )
							stop_r   = .true.
						end if
						! Reset Monte Carlo parameters
						resetmc = .true.
						! Message
						write ( *, 350 ) cycles, angmax, msgaux, r_acc_r
						! Increment counter
						reset_c = reset_c + 1
					end if

					! *******************************************************
					! Ratio data
					! *******************************************************
					write( 40, * ) cycles, ratio, angmax, r_acc_r

					! Reset counter
					naccr = 0

					! *******************************************************
					! Reset Monte Carlo parameters
					! See 'subroutines' code for more details.
					! *******************************************************
					!  OBS.: These parameters were not resetted in the 
					!  Monte Carlo simulations performed in our paper.
					! *******************************************************
					!if ( resetmc ) then
					!	call reset_system()
					!	! Advance loop
					!	cycle
					!end if

				end if

			end if

		end if

		! *******************************************************************************
		! Order parameter data
		! *******************************************************************************
		if ( cycles > n_equil .and. mod( cycles, n_save ) == 0 ) then
			! ***********************************************************************
			! Nematic order parameter (Q-tensor method)
			! See 'subroutines' code for more details.
			! ***********************************************************************
			call order_parameter()
			! Only production-related data
			write ( 50, * ) cycles, s
		end if

		! *******************************************************************************
		! Trajectory data
		! *******************************************************************************
		if ( traj_check ) then
			if ( mod ( cycles, n_save ) == 0 ) then
				write ( 20, "(I4)" ) n_particles
				write ( 20, * ) " "
				do i = 1, n_particles
					! *******************************************************
					! OVITO format (reduced units)
					! *******************************************************
					write ( 20, * ) atom, rmc(1,i), rmc(2,i), rmc(3,i), &
					&		qmc(0,i), qmc(1,i), qmc(2,i), qmc(3,i), &
					&		0.5d0, aspect_ratio
				end do
			end if
		end if

		! *******************************************************************************
		! Simulation progress
		! *******************************************************************************
		if ( mod( cycles, ( 100 * n_save ) ) == 0 ) then
			write ( *, "(A10,F5.1,A1)" ) "Progress: ", ( 100.d0 ) * ( dble( cycles ) / dble( max_cycles ) ), "%"
			write ( *, "(A18,I10)" ) "Remaining cycles: ", ( max_cycles - cycles )
			if ( reset_c > 0 ) then
				write ( *, "(A19,I3)" ) "Threshold resets: ", reset_c
			end if
			write ( *, * ) " "
		end if

	end do

	! ***************************************************************************************
	! End of Metropolis algorithm
	! ***************************************************************************************
	print *, "Monte Carlo simulation finished successfully! See directories for results."
	print *, " "

	! ***************************************************************************************
	! Output units
	! ***************************************************************************************
	if ( traj_check ) then
		close ( 20 )
	end if
	close ( 30 )
	close ( 40 )
	close ( 50 )
	do i = 1, pot_param
		close ( 60 + i )
	end do

 400	format ( "Attractive/Repulsive parameter of ", F0.5, ":", F0.5 )
 
 	! ***************************************************************************************
	! Summary report
	! ***************************************************************************************
 	print *, "Final total potential energy per particle: "
	print *, " "
	do i = 1, pot_param
		write ( *, 400 ) lambda(i), vmc(i)/n_particles
	end do

	call order_parameter()

	write ( *, * ) " "
	write ( *, * ) "Final order parameter: ", s

	! ***************************************************************************************
	! Finish simulation timer
	! ***************************************************************************************
	call cpu_time(stop_timer)

	write ( *, * ) " "
	write ( *, * ) "Simulation Elapsed Time: ", ( stop_timer - start_timer ), " s"
	write ( *, * ) " "

	! ***************************************************************************************
	! Deallocation of arrays
	! ***************************************************************************************
	deallocate ( v, vmc, vn, vm, dv )
	deallocate ( swrange )
	deallocate ( qmc, rmc, emc )

	! ***************************************************************************************
	! Messages
	! ***************************************************************************************
 450	format ( "Now calculating TPT coefficients for an attractive/repulsive parameter of ", F0.5, "..." )

 500	format ( "Calculation of TPT coefficients for an attractive/repulsive parameter of ", F0.5, " finished successfully!" )
 
 550	format ( "The results are reported in a log file in the 'Perturbed_Coefficient' folder." )

	! ***************************************************************************************
	! Output file format
	! ***************************************************************************************
 600	format ( 1X, "'1st_Order_TPTCoefficient'", ",", "'1st_Order_TPTCoefficient_STD'", &
 	&	 ",", "'2nd_Order_TPTCoefficient'", ",", "'2nd_Order_TPTCoefficient_STD'", &
	&	 ",", "'Perturbed_Helmholtz_FEnergy'", ",", "'Perturbed_Helmholtz_FEnergy_STD'" )

 650	format ( 1X, F10.6, ",", F10.6, ",", F10.6, ",", F10.6, ",", F10.6, ",", F10.6 )

	! ***************************************************************************************
	! Thermodynamic Perturbation Theory Coefficients
	! ***************************************************************************************
	if ( coef_check ) then

		print *, "User enabled the calculation of the first- and second-order TPT coefficients. Resuming..."
		print *, " "

		! *******************************************************************************
		! Calculation of TPT coefficients and Perturbed Helmholtz free energy
		! *******************************************************************************
		block_av_loop: do i = 1, pot_param
			if ( ff_selec(1) ) then
				write ( *, 450 ) lambda(i)
				write ( descriptor_lamb, format_lamb ) lambda(i)
			else if ( ff_selec(2) ) then
				write ( *, 450 ) n_repulsive(i)
				write ( descriptor_n, format_n ) n_repulsive(i)
			end if
			write ( *, * ) " "
			! ***********************************************************************
			! Block-averaging method
			! See 'subroutines' code for more details.
			! ***********************************************************************
			call block_averaging(flag,a1(i),a2(i),apert(i),dpa1(i),dpa2(i),dpapert(i))
			! ***********************************************************************
			! Stop condition
			! ***********************************************************************
			if ( flag ) then
				exit block_av_loop	! Terminate calculation if min. block greater than max. blocks
			end if
			! Message
			if ( ff_selec(1) ) then
				write ( *, 500 ) lambda(i)
			else if ( ff_selec(2) ) then
				write ( *, 500 ) n_repulsive(i)
			end if
			write ( *, 550 )
			write ( *, * ) " "
		end do block_av_loop

		! *******************************************************************************
		! TPT Coefficients and Perturbed Helmholtz free energy results
		! *******************************************************************************
		if ( .not. flag ) then	! Subroutine returns no errors

			do i = 1, pot_param

				if ( ff_selec(1) ) then
					write ( descriptor_lamb, format_lamb ) lambda(i)
					open ( UNIT = 150, FILE = "Perturbed_Coefficient/SW/"//trim(descriptor_date)//"/Lambda_" &
					&    //trim(descriptor_lamb)//"/TPT_coefficients_ε"//trim(descriptor_file1)//"_ρ" &
					&    //trim(descriptor_file2)//".csv")
				else if ( ff_selec(2) ) then
					write ( descriptor_n, format_n ) n_repulsive(i)
					open ( UNIT = 150, FILE = "Perturbed_Coefficient/KH/"//trim(descriptor_date)//"/n_" &
					&    //trim(descriptor_n)//"/TPT_coefficients_ε"//trim(descriptor_file1)//"_ρ" &
					&    //trim(descriptor_file2)//".csv")
				end if

					write ( 150, "(A18,I4.4)" ) "No. of particles: ", n_particles
					write ( 150, "(A21,F0.5)" ) "Reduced temperature: ", temp
					write ( 150, "(A24,F0.5)" ) "Reduced number density: ", rho
					write ( 150, "(A12,F0.5)" ) "Aspect ratio (L/D): ", aspect_ratio
					if ( ff_selec(1) ) then
						write ( 150, "(A18,F0.5)" ) "Attractive range: ", lambda(i)
					else if ( ff_selec(2) ) then
						write ( 150, "(A18,F0.5)" ) "Repulsive parameter: ", n_repulsive(i)
					end if
					write ( 150, * ) " "
					write ( 150, 600 )
					write ( 150, 650 ) a1(i), dpa1(i), a2(i), dpa2(i), apert(i), dpapert(i)

				close ( 150 )

			end do

		else if ( flag ) then	! Subroutine returns errors

			write ( *, * ) "ERROR: the total number of blocks is larger than the allowed maximum (1e4)!"
			write ( *, * ) "TPT coefficients and Perturbed Helmholtz free energy not calculated!"
			write ( *, * ) " "

		end if

	end if

	! ***************************************************************************************
	! Write down results
	! ***************************************************************************************
	print *, "Writing simulation log..."
	print *, " "

	! ***************************************************************************************
	! Attractive range format
	! ***************************************************************************************
	format_file0 = "(I3)"
	write ( descriptor_file0, format_file0 ) pot_param
	l_format = "("//trim(descriptor_file0)//"F0.5)"

	! ***************************************************************************************
	! Simulation log descriptors
	! ***************************************************************************************
	write ( char_label(1),  "(I5)"     ) n_particles
	write ( char_label(2),  "(F0.5)"   ) aspect_ratio
	write ( char_label(3),  "(F0.5)"   ) rho
	write ( char_label(4),  "(F0.5)"   ) temp
	if ( coef_check ) then
		write ( char_label(5),  "(I4)" ) min_blocks
		write ( char_label(6),  "(I4)" ) max_blocks
	end if
	if ( ff_selec(1) ) then
		write ( char_label(7),  l_format   ) lambda
	else if ( ff_selec(2) ) then
		write ( char_label(7),  l_format   ) n_repulsive
	end if
	write ( char_label(8),  "(I12.12)" ) max_cycles
	write ( char_label(9),  "(I12.12)" ) n_equil
	write ( char_label(10), "(I12.12)" ) ( max_cycles - n_equil )
	write ( char_label(11), "(I5.5)"   ) n_save
	write ( char_label(12), "(I4.4)"   ) n_adjust
	write ( char_label(13), "(I4.4)"   ) reset_c
	write ( char_label(14), "(F9.2)"   ) ( stop_timer - start_timer )

	! ***************************************************************************************
	! File inquiry
	! ***************************************************************************************
	inquire ( FILE = "Simulation_Log.dat", EXIST = file_exist )

 700	format(I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)

	! ***************************************************************************************
	! Simulation log
	! ***************************************************************************************
	if ( .not. file_exist ) then
		open ( UNIT = 95, FILE = "Simulation_Log.dat" )
			write ( 95, * ) "Monte Carlo Simulation Log"
			write ( 95, * ) " "
			write ( 95, * ) "Canonical Monte Carlo algorithm for ellipsoid-of-revolution molecules"
			write ( 95, * ) " "
			write ( 95, * ) "Nathan Barros de Souza"
			write ( 95, * ) "University of Campinas"
			write ( 95, * ) "School of Chemical Engineering"
			write ( 95, * ) " "
			write ( 95, * ) "Original Developer: Joyce Tavares Lopes"
			write ( 95, * ) "Supervisor: Luís Fernando Mercier Franco"
			write ( 95, * ) "________________________________________________________"
			write ( 95, * ) " "
			write ( 95, "(A14)" ) "Execution Date"
			write ( 95, 700 ) date_time(1), date_time(2), date_time(3), &
			&		  date_time(5), date_time(6), date_time(7)
			write ( 95, * ) " "
			write ( 95, "(A27,5X,A20)" ) "Number of Particles:       ", adjustl ( char_label(1)  )
			write ( 95, "(A27,5X,A20)" ) "Aspect Ratio (L/D):        ", adjustl ( char_label(2)  )
			write ( 95, "(A27,5X,A20)" ) "Reduced Number Density     ", adjustl ( char_label(3)  )
			write ( 95, "(A27,5X,A20)" ) "Reduced Temperature:       ", adjustl ( char_label(4)  )
			if ( coef_check ) then
				write ( 95, "(A27,5X,A20)" ) "Minimum Number of Blocks :    ", adjustl ( char_label(5)  )
				write ( 95, "(A27,5X,A20)" ) "Maximum Number of Blocks :    ", adjustl ( char_label(6)  )
			end if
			if ( ff_selec(1) ) then
				write ( 95, "(A27,5X,A30)" ) "Attractive Range:          ", adjustl ( char_label(7)  )
			else if ( ff_selec(2) ) then
				write ( 95, "(A27,5X,A30)" ) "Repulsive Index:           ", adjustl ( char_label(7)  )
			end if
			write ( 95, "(A27,5X,A20)" ) "Total Number of Cycles:    ", adjustl ( char_label(8)  )
			write ( 95, "(A27,5X,A20)" ) "Equilibration Cycles:      ", adjustl ( char_label(9)  )
			write ( 95, "(A27,5X,A20)" ) "Production Cycles:         ", adjustl ( char_label(10) )
			write ( 95, "(A27,5X,A20)" ) "Data Saving Frequency :    ", adjustl ( char_label(11) )
			write ( 95, "(A27,5X,A20)" ) "Adjustment Frequency:      ", adjustl ( char_label(12) )
			write ( 95, "(A27,5X,A20)" ) "Threshold Resets:          ", adjustl ( char_label(13) )
			write ( 95, * ) " "
			if ( traj_check ) then
				write ( 95, "(A33)" ) "Trajectory of particles computed."
			else if ( .not. traj_check ) then
				write ( 95, "(A37)" ) "Trajectory of particles not computed."
			end if
			if ( coef_check ) then
				write ( 95, "(A62)" ) "TPT Coefficients and Perturbed Helmholtz free energy computed."
			else if ( .not. coef_check ) then
				write ( 95, "(A66)" ) "TPT Coefficients and Perturbed Helmholtz free energy not computed."
			end if
			if ( potential_check ) then
				write ( 95, "(A53)" ) "Potential energy computed only for production cycles."
			else if ( .not. potential_check ) then
				write ( 95, "(A46)" ) "Potential energy computed only for all cycles."
			end if
			if ( ff_selec(1) ) then
				write ( 95, "(A31)" ) "Square-Well potential selected."
			else if ( ff_selec(2) ) then
				write ( 95, "(A26)" ) "Kihara potential selected."
			end if
			write ( 95, * )  " "
			write ( 95, "(A27,5X,A20,A4)")  "Simulation length:            ", &
			&				 adjustl ( char_label(14) ), " sec"
		close ( 95 )
	! ***************************************************************************************
	! Simulation log (appending)
	! ***************************************************************************************
	else if ( file_exist ) then
		open ( UNIT = 95, FILE = "Simulation_Log.dat", POSITION = "append" )
			write ( 95, * ) "________________________________________________________"
			write ( 95, * ) " "
			write ( 95, "(A27,5X,A20)" ) "Number of Particles:       ", adjustl ( char_label(1)  )
			write ( 95, "(A27,5X,A20)" ) "Aspect Ratio (L/D):        ", adjustl ( char_label(2)  )
			write ( 95, "(A27,5X,A20)" ) "Reduced Number Density     ", adjustl ( char_label(3)  )
			write ( 95, "(A27,5X,A20)" ) "Reduced Temperature:       ", adjustl ( char_label(4)  )
			if ( coef_check ) then
				write ( 95, "(A27,5X,A20)" ) "Minimum Number of Blocks :    ", adjustl ( char_label(5)  )
				write ( 95, "(A27,5X,A20)" ) "Maximum Number of Blocks :    ", adjustl ( char_label(6)  )
			end if
			if ( ff_selec(1) ) then
				write ( 95, "(A27,5X,A30)" ) "Attractive Range:          ", adjustl ( char_label(7)  )
			else if ( ff_selec(2) ) then
				write ( 95, "(A27,5X,A30)" ) "Repulsive Index:           ", adjustl ( char_label(7)  )
			end if
			write ( 95, "(A27,5X,A20)" ) "Total Number of Cycles:    ", adjustl ( char_label(8)  )
			write ( 95, "(A27,5X,A20)" ) "Equilibration Cycles:      ", adjustl ( char_label(9)  )
			write ( 95, "(A27,5X,A20)" ) "Production Cycles:         ", adjustl ( char_label(10) )
			write ( 95, "(A27,5X,A20)" ) "Data Saving Frequency :    ", adjustl ( char_label(11) )
			write ( 95, "(A27,5X,A20)" ) "Adjustment Frequency:      ", adjustl ( char_label(12) )
			write ( 95, "(A27,5X,A20)" ) "Threshold Resets:          ", adjustl ( char_label(13) )
			write ( 95, * ) " "
			if ( traj_check ) then
				write ( 95, "(A33)" ) "Trajectory of particles computed."
			else if ( .not. traj_check ) then
				write ( 95, "(A37)" ) "Trajectory of particles not computed."
			end if
			if ( coef_check ) then
				write ( 95, "(A62)" ) "TPT Coefficients and Perturbed Helmholtz free energy computed."
			else if ( .not. coef_check ) then
				write ( 95, "(A66)" ) "TPT Coefficients and Perturbed Helmholtz free energy not computed."
			end if
			if ( potential_check ) then
				write ( 95, "(A53)" ) "Potential energy computed only for production cycles."
			else if ( .not. potential_check ) then
				write ( 95, "(A46)" ) "Potential energy computed only for all cycles."
			end if
			if ( ff_selec(1) ) then
				write ( 95, "(A31)" ) "Square-Well potential selected."
			else if ( ff_selec(2) ) then
				write ( 95, "(A26)" ) "Kihara potential selected."
			end if
			write ( 95, * )  " "
			write ( 95, "(A27,5X,A20,A4)")  "Simulation length:            ", &
			&				 adjustl ( char_label(14) ), " sec"
		close ( 95 )
	end if

	! Status
	print *, "Done!"

end program main
