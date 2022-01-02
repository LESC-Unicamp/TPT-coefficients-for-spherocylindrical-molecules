! ############################################################################################### !
!              Canonical Monte Carlo algorithm for spherocylindrical molecules                    !
!           This module creates folders and subfolders to organize simulation results.            !
! Directories are created by executing a shell command via an intrinsic function called 'system'. !
!               Please note that which shell is used to invoke the command line is                !
!                           system-dependent and environment-dependent.                           !
!         See <https://gcc.gnu.org/onlinedocs/gfortran/SYSTEM.html> for more information.         !
!                  The code below is meant for Ubuntu 20.04 operational systems.                  !
!                   We have not provided an alternative code for Windows users.                   !
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
!                                            January 2nd                                          !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

module folders

	! Uses one module: global variables
	use globalvar

	implicit none

	contains

	! *************************************************************************************** !
	!                          Initialization of parent directories                           !
	! *************************************************************************************** !
	subroutine initfolder()

		implicit none

		! Inquires whether a folder exists and stores the inquiry result in a logical variable
		inquire ( FILE = "Initial_Configuration", EXIST = fexist(1) )

		! *******************************************************************************
		! Initial configuration folder (holds information on the initial molecular structure)
		! *******************************************************************************
		if ( .not. fexist(1) ) then
			call system ( "mkdir Initial_Configuration" )
		end if

		! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
		!  The initial molecular structure at 'OVITO' subfolder is properly formatted to be analyzed by that software.
		inquire ( FILE = "Initial_Configuration/OVITO/", EXIST = sfexist(1) )

		! *******************************************************************************
		! Initial configuration subfolder
		! *******************************************************************************
		if ( .not. sfexist(1) ) then
			call system ( "mkdir Initial_Configuration/OVITO/" )
		end if

		! Date format (YYYY/MM/DD)
		format_date = "(I4,2I2.2)"

		write ( descriptor_date, format_date ) date_time(1), date_time(2), date_time(3)

		! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
		inquire ( FILE = "Initial_Configuration/OVITO/"//trim(descriptor_date)//"/", EXIST = dfexist(1) )

		! *******************************************************************************
		! Date subfolder
		! *******************************************************************************
		if ( .not. dfexist(1) ) then
			call system ( "mkdir Initial_Configuration/OVITO/"//trim(descriptor_date)//"/" )
		end if

	end subroutine initfolder

	! *************************************************************************************** !
	!                          Initialization of parent directories                           !
	! *************************************************************************************** !
	subroutine inquire_folders()

		implicit none

		! Inquires whether a folder exists and stores the inquiry result in a logical variable
		inquire ( FILE = "Trajectories", EXIST = fexist(2) )
		inquire ( FILE = "Potential", EXIST = fexist(3) )
		inquire ( FILE = "Ratio", EXIST = fexist(4) )
		inquire ( FILE = "Order_Parameter", EXIST = fexist(5) )
		inquire ( FILE = "Perturbed_Coefficient", EXIST = fexist(6) )

		! *******************************************************************************
		! Trajectory folder (holds information on orientation and position of particles)
		! *******************************************************************************
		if ( .not. fexist(2) ) then
			call system ( "mkdir Trajectories" )
		end if

		! *******************************************************************************
		! Potential folder (holds information on the potential energy of the system)
		! *******************************************************************************
		if ( .not. fexist(3) ) then
			call system ( "mkdir Potential" )
		end if

		! *******************************************************************************
		! Ratio folder (holds information on the equilibration cycles, like maximum displacement adjustment)
		! *******************************************************************************
		if ( .not. fexist(4) ) then
			call system ( "mkdir Ratio" )
		end if

		! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
		inquire ( FILE = "Ratio/Translation/", EXIST = sfexist(2) )
		inquire ( FILE = "Ratio/Rotation/", EXIST = sfexist(3) )

		! *******************************************************************************
		! Ratio subfolders
		! *******************************************************************************
		if ( .not. sfexist(2) ) then
			call system ( "mkdir Ratio/Translation/" )
		end if
		if ( .not. sfexist(3) ) then
			call system ( "mkdir Ratio/Rotation/" )
		end if

		! *******************************************************************************
		! Order parameter folder (holds information on the nematic order parameter)
		! *******************************************************************************
		if ( .not. fexist(5) ) then
			call system ( "mkdir Order_Parameter" )
		end if

		! *******************************************************************************
		! TPT Coefficients folder (holds information on the perturbation coefficients)
		! *******************************************************************************
		!  These are the terms of the high-temperature series expansion (HTSE) of the 
		!  perturbed Helmholtz free energy.
		! *******************************************************************************
		if ( .not. fexist(6) ) then
			call system ( "mkdir Perturbed_Coefficient" )
		end if

	end subroutine inquire_folders

	! *************************************************************************************** !
	!                            Initialization of date subfolders                            !
	! *************************************************************************************** !
	subroutine date_folders()

		implicit none

		! Date format (YYYY/MM/DD)
		format_date = "(I4,2I2.2)"

		write ( descriptor_date, format_date ) date_time(1), date_time(2), date_time(3)

		! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
		inquire ( FILE = "Trajectories/"//trim(descriptor_date)//"/", EXIST = dfexist(2) )
		inquire ( FILE = "Potential/"//trim(descriptor_date)//"/", EXIST = dfexist(3) )
		inquire ( FILE = "Ratio/Translation/"//trim(descriptor_date)//"/", EXIST = dfexist(4) )
		inquire ( FILE = "Ratio/Rotation/"//trim(descriptor_date)//"/", EXIST = dfexist(5) )
		inquire ( FILE = "Order_Parameter/"//trim(descriptor_date)//"/", EXIST = dfexist(6) )
		inquire ( FILE = "Perturbed_Coefficient/"//trim(descriptor_date)//"/", EXIST = dfexist(7) )

		! *******************************************************************************
		! Date subfolders
		! *******************************************************************************
		if ( .not. dfexist(2) ) then
			call system ( "mkdir Trajectories/"//trim(descriptor_date)//"/" )
		end if
		if ( .not. dfexist(3) ) then
			call system ( "mkdir Potential/"//trim(descriptor_date)//"/" )
		end if
		if ( .not. dfexist(4) ) then
			call system ( "mkdir Ratio/Translation/"//trim(descriptor_date)//"/" )
		end if
		if ( .not. dfexist(5) ) then
			call system ( "mkdir Ratio/Rotation/"//trim(descriptor_date)//"/" )
		end if
		if ( .not. dfexist(6) ) then
			call system ( "mkdir Order_Parameter/"//trim(descriptor_date)//"/" )
		end if
		if ( .not. dfexist(7) ) then
			call system ( "mkdir Perturbed_Coefficient/"//trim(descriptor_date)//"/" )
		end if

	end subroutine date_folders

	! *************************************************************************************** !
	!                    Initialization of attractive parameter subfolders                    !
	! *************************************************************************************** !
	subroutine lambda_folders()

		implicit none

		! *******************************************************************************
		! Attractive parameter format (free-width)
		!  Might be necessary to change if the number of decimal places is higher than 5.
		! *******************************************************************************
		format_lamb = "(F0.5)"

		! *******************************************************************************
		! Attractive parameter subfolders
		! *******************************************************************************
		do counter_lambda = 1, n_lambda
			write ( descriptor_lamb, format_lamb ) lambda(counter_lambda)
			! Parent folder: Potential
			inquire ( FILE = "Potential/"//trim(descriptor_date)//"/Lambda_"//trim(descriptor_lamb)//"/", &
			& 	  EXIST = lfexist(counter_lambda) )
			if ( .not. lfexist(counter_lambda) ) then
				call system ( "mkdir Potential/"//trim(descriptor_date)//"/Lambda_"//trim(descriptor_lamb)//"/" )
			end if
			! Parent folder: TPT Coefficients
			inquire ( FILE = "Perturbed_Coefficient/"//trim(descriptor_date)//"/Lambda_"//trim(descriptor_lamb)//"/", &
			&	  EXIST = lfexist(counter_lambda) )
			if ( .not. lfexist(counter_lambda) ) then
				call system ( "mkdir Perturbed_Coefficient/"//trim(descriptor_date)//"/Lambda_" &
			&		       //trim(descriptor_lamb)//"/" )
			end if
		end do

	end subroutine lambda_folders

	! *************************************************************************************** !
	!                    Initialization of repulsive parameter subfolders                     !
	! *************************************************************************************** !
	subroutine n_folders()

		implicit none

		! *******************************************************************************
		! Repulsive parameter format (free-width)
		!  Might be necessary to change if the number of decimal places is higher than 5.
		! *******************************************************************************
		format_n = "(F0.5)"

		! *******************************************************************************
		! Repulsive parameter subfolders
		! *******************************************************************************
		do counter_n = 1, n_n
			write ( descriptor_n, format_n ) n_repulsive(counter_n)
			! Parent folder: Potential
			inquire ( FILE = "Potential/"//trim(descriptor_date)//"/n_"//trim(descriptor_n)//"/", &
			& 	  EXIST = lfexist(counter_n) )
			if ( .not. lfexist(counter_n) ) then
				call system ( "mkdir Potential/"//trim(descriptor_date)//"/n_"//trim(descriptor_n)//"/" )
			end if
			! Parent folder: TPT Coefficients
			inquire ( FILE = "Perturbed_Coefficient/"//trim(descriptor_date)//"/n_"//trim(descriptor_n)//"/", &
			&	  EXIST = lfexist(counter_n) )
			if ( .not. lfexist(counter_n) ) then
				call system ( "mkdir Perturbed_Coefficient/"//trim(descriptor_date)//"/n_" &
			&		       //trim(descriptor_n)//"/" )
			end if
		end do

	end subroutine n_folders

end module folders
