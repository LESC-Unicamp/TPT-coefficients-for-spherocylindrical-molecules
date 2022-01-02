! ############################################################################################### !
!                 Canonical Monte Carlo algorithm for spherocylindrical molecules                 !
!           This module defines the variables used by the main program and most of the            !
!         subroutines and functions. A brief description is presented for each variable.          !
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
!                                 C. R. A. Abreu, F. A. Escobedo                                  !
!                                J. Chem. Phys. 124, 054116 (2006)                                !
!                                     DOI: 10.1063/1.2165188                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

module globalvar

	implicit none

	! ***************************************************************************************
	! Integer variables -*- THIS IS SINGLE PRECISION -*-
	! ***************************************************************************************
	integer, dimension (8)	:: date_time		! Computer clock (date and time)

	! ***************************************************************************************
	! Integer variables (General)
	! ***************************************************************************************
	integer*8		:: seed			! Random number generator seed
	integer*8		:: n_particles		! Number of particles
	integer*8		:: n_cells		! Number of unit cells
	integer*8		:: n_lambda		! Number of attractive range parameters
	integer*8		:: counter_lambda	! Counter (Square Well)
	integer*8		:: n_n		        ! Number of repulsive parameters
	integer*8		:: counter_n	        ! Counter (Kihara)

	! ***************************************************************************************
	! Integer variables (Monte Carlo parameters)
	! ***************************************************************************************
	integer*8		:: cycles		! Counter of cycles
	integer*8		:: reset_c		! Counter of Monte Carlo parameter resets
	integer*8		:: max_cycles		! Total number of cycles
	integer*8		:: n_equil		! Number of equilibration cycles
	integer*8		:: n_save		! Results saving frequency
	integer*8		:: n_adjust		! Maximum displacement adjustment frequency
	integer*8		:: nacct		! Move acceptance counter: Translation
	integer*8		:: naccr		! Move acceptance counter: Rotation
	integer*8		:: movt			! Move counter (Translation)
	integer*8		:: movr			! Move counter (Rotation)

	! ***************************************************************************************
	! Integer variables (Block Averaging)
	! ***************************************************************************************
	integer*8		:: min_blocks		! Minimum number of blocks
	integer*8		:: max_blocks		! Maximum number of blocks

	! ***************************************************************************************
	! Real variables (General)
	! ***************************************************************************************
	real*8			:: quaternion_angle	! Quaternion angle [real part, W] (for initial configuration only)
	real*8			:: cell_length		! Length of unit cell (cubic structure)
	real*8, dimension (3)	:: coordinate		! Length (x,y,z) of unit cell (stretched/narrowed cubic structure)
	real*8, dimension (3)	:: box_length		! Length (x,y,z) of simulation box
	real*8			:: rho			! Reduced number density (κρ*)
	real*8			:: temp			! Reduced temperature
	real*8			:: sigsphere		! Spherical diameter
	real*8			:: aspect_ratio		! Spherocylindrical aspect ratio parameter (L/D)
	real*8			:: s			! Nematic order parameter
	real*8			:: maxarg		! Largest exponential argument (see Supplementary Material of Abreu and Escobedo, J. Chem. Phys. 124)
	real*8			:: beta_energy		! Boltzmann factor argument, -βU
	real*8			:: check_root		! Checks cube root for the number of particles, N, in each structure (SC, BCC, and FCC)
	real*8			:: random_n		! Random number from a pseudorandom number generator subroutine
	real*8			:: random_angle		! Random angle
	real*8			:: start_timer		! Start timer of Monte Carlo simulation
	real*8			:: stop_timer		! Stop timer of Monte Carlo simulation
	real*8			:: anynumber		! Dummy (number)
	real*8, dimension (3)	:: fixed_axis		! Body-fixed axis of rotation (see Allen and Tildesley, page 111)

	! ***************************************************************************************
	! Real variables (Monte Carlo parameters)
	! ***************************************************************************************
	real*8			:: max_trans		! User maximum displacement [+/-] (Translation)
	real*8			:: max_rot		! User maximum displacement [+/-] (Rotation)
	real*8			:: drmax		! Maximum displacement [+/-] (Translation)
	real*8			:: angmax		! Maximum displacement [+/-] (Rotation)
	real*8			:: r_acc_t		! Acceptance ratio threshold (Translation)
	real*8			:: r_acc_r		! Acceptance ratio threshold (Rotation)
	real*8			:: ratio		! Acceptance ratio (Simulation)

	! ***************************************************************************************
	! Real variables (allocatable)
	! ***************************************************************************************
	real*8, dimension (:), allocatable	:: lambda	! Attractive range parameter
	real*8, dimension (:), allocatable	:: swrange	! Effective range of attraction
	real*8, dimension (:), allocatable	:: n_repulsive	! Repulsive parameter
	real*8, dimension (:), allocatable	:: v, vmc	! Potential energy array
	real*8, dimension (:,:), allocatable	:: q, qmc	! Quaternion array
	real*8, dimension (:,:), allocatable	:: r, rmc	! Position array
	real*8, dimension (:,:), allocatable	:: e, emc	! Orientation array

	! ***************************************************************************************
	! Real parameters
	! ***************************************************************************************
	real*8, parameter	:: pi = 4.d0 * datan ( 1.d0 )

	! ***************************************************************************************
	! Character variables
	! ***************************************************************************************
	character			:: atom*1		! Atom ID
	character			:: descriptor_file0*10	! Descriptor for output file
	character			:: descriptor_file1*10	! Descriptor for output file
	character			:: descriptor_file2*10	! Descriptor for output file
	character			:: descriptor_date*10	! Descriptor for output folder
	character			:: descriptor_lamb*10	! Descriptor for output folder
	character			:: descriptor_n*10	! Descriptor for output folder
	character			:: format_file0*32	! String format for output file
	character			:: format_file1*32	! String format for output file
	character			:: format_file2*32	! String format for output file
	character			:: format_date*32	! String format for output folder
	character			:: format_lamb*32	! String format for output folder
	character			:: format_n*32		! String format for output folder
	character			:: date_clock*8		! Computer clock (date)
	character			:: time_clock*10	! Computer clock (time)
	character			:: traj_inq*1		! Trajectory output inquiry
	character			:: pot_inq*1		! Potential output inquiry
	character			:: coef_inq*1		! Coefficient output inquiry
	character			:: config_inq*3		! Initial configuration inquiry
	character			:: ff_inq*3		! Force field inquiry
	character			:: dummy*1		! Dummy (character)
	character			:: get*100		! Variable names (.sci input file)
	character, dimension (14)	:: char_label*64	! Simulation log

	! ***************************************************************************************
	! Logical variables
	! ***************************************************************************************
	logical			:: file_exist		! Checks whether a file exists or not
	logical			:: traj_check		! Checks inquiry of trajectory output
	logical			:: potential_check	! Checks inquiry of potential output
	logical			:: coef_check		! Checks inquiry of coefficient output
	logical			:: mov_rot		! Rotation move selection				: TRUE = movement selected; FALSE = movement not selected
	logical			:: mov_trans		! Translation movement selection			: TRUE = movement selected; FALSE = movement not selected
	logical			:: stop_r		! Ratio threshold modifier (rotational move)		: TRUE = stop modification; FALSE = continue modification
	logical 		:: stop_t		! Ratio threshold modifier (translational move)		: TRUE = stop modification; FALSE = continue modification
	logical			:: resetmc		! Monte Carlo parameters reset				: TRUE = reset parameters;  FALSE = keep parameters
	logical, dimension (2)	:: ff_selec		! Checks the selected force field
	logical, dimension (3)	:: config_selec		! Checks the selected molecular configuration
	logical, dimension (6)	:: fexist		! Checks whether folder exists or not
	logical, dimension (7)	:: dfexist		! Checks whether date folders exist or not
	logical, dimension (3)	:: sfexist		! Checks whether subfolder exists or not

	! ***************************************************************************************
	! Logical variables (allocatable)
	! ***************************************************************************************
	logical, dimension (:), allocatable	:: lfexist	! Checks whether the attractive/repulsive parameter subfolders exist or not

end module globalvar
