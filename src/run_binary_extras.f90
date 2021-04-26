! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! *********************************************************************** 
      module run_binary_extras 

      use star_def
      use const_def
      use const_def
      use chem_def
      use binary_def
      use cc_def
      use ce_def
      use run_star_support
      
      use star_lib
      use math_lib
      use num_lib
      use cc_lib
      use ce_lib

      implicit none

      include 'test_suite_extras_def.inc'

      integer, parameter :: max_num_kicks = 100000

      ! b% lxtra(lx_ce_on) is true when ce is on
      integer, parameter :: lx_ce_on = 1
      ! b% lxtra(lx_ce_off) is true when ce is off
      integer, parameter :: lx_ce_off = 2

      ! b% ixtra(ix_ce_model_number) contains the model number of a ce start
      integer, parameter :: ix_ce_model_number = 1
      ! b% ixtra(ix_ce_num_star_1) contains the number of times star 1 has been donor in a ce
      integer, parameter :: ix_ce_num_star_1 = 2
      ! b% ixtra(ix_ce_num_star_1) contains the number of times star 2 has been donor in a ce
      integer, parameter :: ix_ce_num_star_2 = 3
      ! b% ixtra(ix_num_switches) contains the number of donor switches of a run
      integer, parameter :: ix_num_switches = 4

      ! b% xtra(x_cumulative_binding_energy) contains the cumulative Ebind removed
      integer, parameter :: x_cumulative_binding_energy = 1
      ! b% xtra(x_binding_energy_prev_step) contains Ebind on previous timestep
      integer, parameter :: x_binding_energy_prev_step = 2
      ! b% xtra(x_orbital_energy_prev_step) contains Eorb on previous timestep
      integer, parameter :: x_orbital_energy_prev_step = 3
      ! b% xtra(x_ce_duration) contains the duration of the CE
      integer, parameter :: x_ce_duration = 4
      ! b% xtra(x_ce_years_in_detach) contains time of binary in detach after CE
      integer, parameter :: x_ce_years_in_detach = 5

      ! bin2dco variables
      character(len=strlen) :: star_plus_star_filename, star_plus_pm_filename
      character(len=strlen) :: cc1_inlist_filename, cc2_inlist_filename
      character(len=strlen) :: ce1_inlist_filename, ce2_inlist_filename
      logical :: do_star_plus_star, stop_after_star_plus_star
      logical :: do_kicks, do_kicks_in_one_run
      character(len=strlen) :: natal_kicks_filename
      character(len=strlen) :: star_info_at_cc_filename, binary_info_at_cc_filename
      logical :: add_kick_id_as_suffix
      character(len=strlen) :: termination_codes_folder
      namelist /bin2dco_controls/ &
         do_star_plus_star, &
         star_plus_star_filename, star_plus_pm_filename, &
         cc1_inlist_filename, cc2_inlist_filename, &
         ce1_inlist_filename, ce2_inlist_filename, &
         stop_after_star_plus_star, &
         do_kicks, do_kicks_in_one_run, &
         natal_kicks_filename, &
         star_info_at_cc_filename, binary_info_at_cc_filename, &
         add_kick_id_as_suffix, &
         termination_codes_folder

      ! only used together with natal-kicks
      real(dp) :: mass_of_progenitor
      real(dp) :: mass_of_remnant, mass_of_companion
      real(dp) :: pre_cc_separation, pre_cc_period

      integer, parameter :: header_lines_in_natal_kicks_file = 9

      ! post-kick id (changed from run to run)
      character(len=strlen) :: kick_id
      real(dp) :: porb_kick, a_kick, ecc_kick

      ! utilities
      integer :: num_kicks
      integer :: num_switches
      logical :: second_collapse

      logical, parameter :: dbg = .false.
      
      contains


      include 'test_suite_extras.inc'


      subroutine do_run(extras_controls, extras_binary_controls, ierr)
         integer, intent(out) :: ierr
         character(len=strlen) :: inlist_fname, termination_code_fname
         character(len=strlen) :: str_index
         integer :: iounit, k, k_idx
         logical :: has_reached_cc
         logical :: already_done_star_plus_star

         interface

            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls

            subroutine extras_binary_controls(binary_id, ierr)
               integer :: binary_id
               integer, intent(out) :: ierr
            end subroutine extras_binary_controls

         end interface

         include 'formats'

         write(*,'(/,a)') 'do bin2dco run'

         ! avoid trying to get MESA_INLIST from the command-line arguments
         MESA_INLIST_RESOLVED = .true.

         ! get bin2dco controls
         inlist_fname = 'inlist_bin2dco'
         call read_bin2dco_controls(inlist_fname, ierr)
         if (ierr /= 0) return

         ! check for a hidden file created when the star + star has already been computed
         already_done_star_plus_star = .false.
         open(newunit=iounit, file='.skip_star_plus_star', status='old', action='read', iostat=ierr)
         if (ierr == 0) then
            already_done_star_plus_star = .true.
            close(iounit)
         end if

         ! star + star simulation (if chosen)
         if (do_star_plus_star .and. .not. already_done_star_plus_star) then
            write(*,'(a)') 'going to do star + star'
            call run1_binary(.true., extras_controls, extras_binary_controls, &
               ierr, &
               star_plus_star_filename)

            if (stop_after_star_plus_star) then
               write(*,'(a)') 'will not continue after star + star ended.'
               return
            end if

            if (do_kicks) then
               write(termination_code_fname,'(a)') trim(termination_codes_folder) // &
                  '/termination_code_star_plus_star'
               has_reached_cc = end_as_core_collapse(termination_code_fname)
               if (.not. has_reached_cc) then
                  write(*,'(a)') 'did not reach core-collapse.'
                  return
               end if
            end if
         end if

         ! two possible ways of doing star + point-mass evolution: either grabbing
         ! asymmetric kicks from a file after which binary parameters will be updated,
         ! or a typical MESAbinary evolution in which binary parameters needs to be
         ! set in the file `inlist_star_plus_pm`
         if (do_kicks) then

            ! get masses from files produced by core_collapse module
            call read_cc_parameters(star_info_at_cc_filename, binary_info_at_cc_filename, &
               mass_of_progenitor, mass_of_remnant, mass_of_companion, &
               pre_cc_separation, pre_cc_period, &
               ierr)

            if (dbg) then
               write(*,'(a)')
               write(*,'(a)') 'binary at core-collapse:'
               write(*,'(a)')
               write(*,1) 'progenitor_mass', mass_of_progenitor
               write(*,1) 'remnant_mass', mass_of_remnant
               write(*,1) 'companion_mass', mass_of_companion
               write(*,1) 'separation', pre_cc_separation
               write(*,1) 'period', pre_cc_period
               write(*,'(a)')
            end if

            if (do_kicks_in_one_run) then

               ! reset initilized_*_handles to false, else after 10 runs, MESA will not continue evolving binaries
               have_initialized_binary_handles = .false.
               have_initialized_star_handles = .false.

               ! read how many simulations are needed
               num_kicks = number_of_kicks(natal_kicks_filename, &
                  header_lines_in_natal_kicks_file, &
                  ierr)
               if (ierr /= 0) stop 'failed to get number of kicks in ' // trim(natal_kicks_filename)

               ! find out if this is a restart, in which case k must be upgraded to this value
               open(newunit=iounit, file='.last_kick_id', status='old', action='read', iostat=ierr)
               if (ierr == 0) then
                  read(iounit,*,iostat=ierr) k_idx
                  if (ierr /= 0) stop 'failed to read k-index for restart of kicks'
                  close(iounit)
               else
                  k_idx = 1
               endif

               do k=k_idx, num_kicks

                  ! get orbital parameters from file
                  call read_natal_kick(natal_kicks_filename, &
                     k+header_lines_in_natal_kicks_file, &
                     kick_id, porb_kick, a_kick, ecc_kick, &
                     ierr)
                  
                  if (dbg) then
                     write(*,'(a)')
                     write(*,'(a)') 'binary after core-collapse:'
                     write(*,'(a)')
                     write(*,11) 'k-index', k
                     write(*,'(a55,a26)') 'kick_id', trim(kick_id)
                     write(*,1) 'period', porb_kick
                     write(*,1) 'separation', a_kick
                     write(*,1) 'eccentricty', ecc_kick
                     write(*,'(a)')
                  end if

                  ! call MESAbinary run
                  write(*,'(a)') 'going to do star + point-mass after a natal kick'
                  call run1_binary(.true., extras_controls, extras_binary_controls, &
                     ierr, &
                     star_plus_pm_filename)

                  ! write k value into a file for restarts
                  open(newunit=iounit, file='.last_kick_id', action='write', iostat=ierr)
                  if (ierr /= 0) stop 'failed to write k-index for restart of kicks'
                  write(iounit,*) k
                  close(iounit)

               end do

            else

               ! get orbital parameters post kick from command-line argument
               call get_command_argument(1, str_index, status=ierr)
               if (ierr /= 0) stop 'failed to get index of ' // trim(natal_kicks_filename)
               ! convert str_index to integer
               read(str_index,*) k

               ! get orbital parameters from file
               call read_natal_kick(natal_kicks_filename, &
                  k, kick_id, porb_kick, a_kick, ecc_kick, &
                  ierr)

               ! call MESAbinary run
               write(*,'(a)') 'going to do star + point-mass after a natal kick'
               call run1_binary(.true., extras_controls, extras_binary_controls, &
                  ierr, &
                  star_plus_pm_filename)

            end if

         else
            
            write(*,'(a)') 'going to do star + point-mass'
            call run1_binary(.true., extras_controls, extras_binary_controls, &
               ierr, &
               star_plus_pm_filename)

         end if

         contains

         ! include bin2dco utility functions

         include 'bin2dco_misc.inc'

      end subroutine do_run
      
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! this is for having the correct values in the header of the log
         if (.not. b% job% evolve_both_stars .and. do_kicks) then
            b% m2 = mass_of_remnant
            b% initial_period_in_days = porb_kick
            b% initial_eccentricity = ecc_kick
         end if

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup => extras_binary_startup
         b% extras_binary_start_step => extras_binary_start_step
         b% extras_binary_check_model => extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve => extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls


      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items


      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id

         how_many_extra_binary_history_columns = 7

      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: mu, I1, I2

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! number of switches of donor star
         names(1) = 'num_switches'
         vals(1) = num_switches

         ! ce output
         names(2) = 'ce_phase'
         names(3) = 'ce_num_star_1'
         names(4) = 'ce_num_star_2' 
         names(5) = 'ce_cumulative_Ebind'
         if (ce_on) then
            vals(2) = 1d0
         else
            vals(2) = 0d0
         end if
         vals(3) = ce_num_star_1
         vals(4) = ce_num_star_2
         vals(5) = cumulative_removed_binding_energy

         ! outer lagrangian point equivalent radii
         names(6) = 'rl2_donor'
         vals(6) = eval_outer_roche_lobe(b% m(b% d_i), b% m(b% a_i), b% separation) / Rsun


         ! Darwin unstable separation
         names(7) = 'a_Darwin'
         
         mu = b% m(1) * b% m(2) / (b% m(1) + b% m(2))
         I1 = 0d0
         I2 = 0d0
         if (b% point_mass_i /= 1) &
            I1 = I1 + dot_product(b% s1% dm_bar(1:b% s1% nz), b% s1% i_rot(1:b% s1% nz))
         if (b% point_mass_i /= 2) &
            I2 = I2 + dot_product(b% s1% dm_bar(1:b% s1% nz), b% s1% i_rot(1:b% s1% nz))

         vals(7) = sqrt(3 * (I1 + I2) / mu) / Rsun
         
      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart

         include 'formats'

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
            
         ! set ce controls
         if (b% point_mass_i == 0) then
            call ce_set_controls(ce1_inlist_filename, ierr)
            if (ierr /= 0) return
         else
            call ce_set_controls(ce2_inlist_filename, ierr)
            if (ierr /= 0) return
         end if

         ! update log names to use kick_id
         if (b% point_mass_i /= 0 .and. do_kicks) then
            write(b% history_name, '(a)') 'binary_history_' // trim(kick_id) // '.data'
            write(b% s_donor% star_history_name, '(a)') 'history_' // trim(kick_id) // '.data'
            write(b% s1% Grid1_file_dir, '(a)') 'png/png_' // trim(kick_id)
            write(b% s1% Grid1_file_prefix, '(a)') trim('grid1_') // trim(kick_id) // '_'
         end if

         ! make different stuff for restart or new run
         if (.not. restart) then

            b% lxtra(lx_ce_on) = .false.
            b% lxtra(lx_ce_off) = .true.
            b% ixtra(ix_ce_model_number) = 0
            b% ixtra(ix_ce_num_star_1) = 0
            b% ixtra(ix_ce_num_star_2) = 0
            b% ixtra(ix_num_switches) = 0
            b% xtra(x_cumulative_binding_energy) = 0d0
            b% xtra(x_binding_energy_prev_step) = 0d0
            b% xtra(x_orbital_energy_prev_step) = 0d0
            b% xtra(x_ce_duration)= 0d0

            ce_on = .false.
            ce_off = .true.

         else

            num_switches = b% ixtra(ix_num_switches)

            if (b% lxtra(lx_ce_on)) then
               ! call this function as they will compute (wrongly) energies. use MESA extra vars
               ! to correct this
               call ce_store_initial_info(binary_id, ierr)
               call ce_init_binary_controls(binary_id, ierr)
               ! then set everything correctly
               ce_on = .true.
               ce_off = .false.
               ce_donor_id = b% d_i
               ce_accretor_id = b% a_i
               ce_initial_model_number = b% ixtra(ix_ce_model_number)
               ce_num_star_1 = b% ixtra(ix_ce_num_star_1)
               ce_num_star_2 = b% ixtra(ix_ce_num_star_2)
               cumulative_removed_binding_energy = b% xtra(x_cumulative_binding_energy)
               donor_bind_energy_prev_step = b% xtra(x_binding_energy_prev_step)
               orbital_energy_prev_step = b% xtra(x_orbital_energy_prev_step)
               ce_duration = b% xtra(x_ce_duration)
               ce_years_in_detach = b% xtra(x_ce_years_in_detach)
            end if

            if (b% lxtra(lx_ce_off)) then
               ce_on = .false.
               ce_off = .true.
            end if

         end if

         ! when doing kicks, change output name to contain kick id
         if (do_kicks .and. b% point_mass_i /= 0 .and. add_kick_id_as_suffix) then
            write(filename_donor_profile_pre_ce, '(a)') &
               trim(filename_donor_profile_pre_ce) // "_" // trim(kick_id)
            write(filename_donor_profile_after_ce, '(a)') &
               trim(filename_donor_profile_after_ce) // "_" // trim(kick_id)
            write(filename_donor_model_pre_ce, '(a)') &
               trim(filename_donor_model_pre_ce) // "_" // trim(kick_id)
            write(filename_donor_model_after_ce, '(a)') &
               trim(filename_donor_model_after_ce) // "_" // trim(kick_id)
         end if
         
         extras_binary_startup = keep_going

      end function extras_binary_startup
      

      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         
         extras_binary_start_step = keep_going

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
            
         if (do_kicks .and. b% point_mass_i /= 0 .and. add_kick_id_as_suffix) then
         end if

         ! turn pgstar flag on
         if (b% doing_first_model_of_run) then
            b% s_donor% job% pgstar_flag = .true.
            if (b% point_mass_i == 0) b% s_accretor% job% pgstar_flag = .true.
         end if

         ! check ce end
         if (ce_on) then
            call ce_check_state(binary_id, ierr)
            if (ierr /= 0) return
            if (ce_merge) then
               extras_binary_start_step = terminate
               return
            end if
         end if

      end function extras_binary_start_step
     

      !Return either keep_going, retry or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         
         include 'formats'
         
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         
         extras_binary_check_model = keep_going
         
         if (b% model_number /= b% s_donor% model_number) &
            b% model_number = b% s_donor% model_number
         
         ! check if ce
         if (b% r(b% d_i) > b% rl(b% d_i) .and. ce_off) then
            call ce_unstable_mt_phase(binary_id, ierr)
            if (ierr /= 0) return

            if (ce_on) then
               call ce_init(binary_id, ierr)
               if (ierr /= 0) return
            end if
         end if

      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         integer, intent(in) :: binary_id
         integer :: ierr
         logical :: do_switch
         integer :: star_cc_id
         integer :: number_io
         real(dp) :: mu, I1, I2, a_Darwin
         
         include 'formats'

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  

         extras_binary_finish_step = keep_going

         ! check for Darwin unstable binary
         mu = b% m(1) * b% m(2) / (b% m(1) + b% m(2))
         I1 = 0d0
         I2 = 0d0
         if (b% point_mass_i /= 1) &
            I1 = I1 + dot_product(b% s1% dm_bar(1:b% s1% nz), b% s1% i_rot(1:b% s1% nz))
         if (b% point_mass_i /= 2) &
            I2 = I2 + dot_product(b% s1% dm_bar(1:b% s1% nz), b% s1% i_rot(1:b% s1% nz))
         a_Darwin = sqrt(3 * (I1 + I2) / mu)

         if (b% doing_first_model_of_run .and. b% separation < a_Darwin) then
            write(*,'(a)') 'system is Darwin unstable'
            b% s_donor% termination_code = t_xtra1
            termination_code_str(t_xtra1) = 'Darwin unstable'
            extras_binary_finish_step = terminate
            return
         end if

         ! if mass-transfer is really high even though no RLOF, then terminate
         if (b% r(b% d_i) < b% rl(b% d_i) .and. ce_off &
            .and. abs(b% mtransfer_rate) * secyer/Msun > max_mdot_rlof) then
            write(*,'(a)') 'reach a really high MT rate without having RLOF'
            b% s_donor% termination_code = t_xtra1
            termination_code_str(t_xtra1) = 'mdot_atmospheric > max_mdot_rlof'
            extras_binary_finish_step = terminate
            return
         end if

         ! take success step on for ce vars
         if (ce_on) then
            call ce_step_success(binary_id, ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed in ce_step_success'
               return
            end if
            
            ! we need to store some stuff on binary pointer for restarts
            b% lxtra(lx_ce_on) = .true.
            b% lxtra(lx_ce_off) = .false.
            b% ixtra(ix_ce_model_number) = ce_initial_model_number
            b% ixtra(ix_ce_num_star_1) = ce_num_star_1
            b% ixtra(ix_ce_num_star_2) = ce_num_star_2
            b% xtra(x_cumulative_binding_energy) = cumulative_removed_binding_energy
            b% xtra(x_binding_energy_prev_step) = donor_bind_energy_prev_step
            b% xtra(x_orbital_energy_prev_step) = orbital_energy_prev_step
            b% xtra(x_ce_duration) = ce_duration
            b% xtra(x_ce_years_in_detach) = ce_years_in_detach

            ! also, save model if first ce step
            if (b% model_number == b% ixtra(ix_ce_model_number) .and. save_model_pre_ce) then
               write(*,11) 'save model at start of ce, model_number', b% model_number
               call do_saves_for_binary(b, ierr)
            end if
         end if

         ! first collapse
         star_cc_id = 0
         if (b% point_mass_i == 0) then
            if (b% s1% center_c12 < 1d-4 .and. b% s1% center_he4 < 1d-4) then
               star_cc_id = 1
               call star_write_model(2, 'companion_at_core_collapse.mod', ierr)
               b% s1% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'core-collapse'
               extras_binary_finish_step = terminate
            else if (b% s2% center_c12 < 1d-4 .and. b% s2% center_he4 < 1d-4) then
               star_cc_id = 2
               call star_write_model(1, 'companion_at_core_collapse.mod', ierr)
               b% s2% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'core-collapse'
               extras_binary_finish_step = terminate
            end if
         ! second collapse
         else if (b% point_mass_i == 1) then
            if (b% s2% center_c12 < 1d-4 .and. b% s2% center_he4 < 1d-4) then
               second_collapse = .true.
               star_cc_id = 2
               b% s2% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'core-collapse'
               extras_binary_finish_step = terminate
            end if
         else if (b% point_mass_i == 2) then
            if (b% s1% center_c12 < 1d-4 .and. b% s1% center_he4 < 1d-4) then
               second_collapse = .true.
               star_cc_id = 1
               b% s1% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'core-collapse'
               extras_binary_finish_step = terminate
            end if
         end if

         if (star_cc_id > 0 .and. .not. second_collapse) then

            call cc_set_controls(cc1_inlist_filename, ierr)
            if (ierr /= 0) return

            if (do_kicks) then
               write(filename_for_star_data, '(a)') trim(filename_for_star_data) // "_1"
               write(filename_for_binary_data, '(a)') trim(filename_for_binary_data) // "_1"
            end if

            call cc_compact_object_formation(star_cc_id, ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed in cc_compact_object_formation'
               return
            end if

            ! save profile of first collapsing star
            if (star_cc_id == 1) then
               s => b% s1
            else
               s => b% s2
            end if

            write(*,*) trim(s% log_directory) // '/profile_at_cc.data'

            call star_write_profile_info(star_cc_id, trim(s% log_directory) // '/profile_at_cc.data', ierr)
            if (ierr /= 0) return
            
            return

         else if (star_cc_id > 0 .and. second_collapse) then
               
            ! save profile of first collapsing star inside log directory
            if (star_cc_id == 1) then
               s => b% s1
            else
               s => b% s2
            end if

            if (do_kicks) then
               call cc_set_controls(cc2_inlist_filename, ierr)
               if (ierr /= 0) return

               ! replace filenames by adding natal-kick id
               write(filename_for_star_data, '(a)') trim(filename_for_star_data) // "_" // trim(kick_id)
               write(filename_for_binary_data, '(a)') trim(filename_for_binary_data) // "_" // trim(kick_id)
               
               call cc_compact_object_formation(star_cc_id, ierr)
               if (ierr /= 0) then
                  write(*,'(a)') 'failed in cc_compact_object_formation'
                  return
               end if

               call star_write_profile_info(star_cc_id, &
                  trim(s% log_directory) // '/profile_at_cc' // '_' // trim(kick_id) // '.data', ierr)
               if (ierr /= 0) return

            else

               call cc_set_controls(cc1_inlist_filename, ierr)
               if (ierr /= 0) return
               
               if (do_kicks) then
                  write(filename_for_star_data, '(a)') trim(filename_for_star_data) // "_1"
                  write(filename_for_binary_data, '(a)') trim(filename_for_binary_data) // "_1"
               end if

               call cc_compact_object_formation(star_cc_id, ierr)
               if (ierr /= 0) then
                  write(*,'(a)') 'failed in cc_compact_object_formation'
                  return
               end if
            
               ! save profile of first collapsing star
               call star_write_profile_info(star_cc_id, &
                  trim(s% log_directory) // '/profile_at_second_cc.data', ierr)
               if (ierr /= 0) return
            
            end if

            return
         end if

         ! check if we need to switch donor star only when ce is off
         if (ce_off) then
            do_switch = switch_donor_star(binary_id)
            if (do_switch) then
               write(*,'(a)') 'switching donor'
               num_switches = num_switches + 1
               b% ixtra(ix_num_switches) = num_switches
               if (b% d_i == 2) then
                  b% d_i = 1
                  b% d_i_old = 1
                  b% a_i = 2
                  b% a_i_old = 2
                  b% s_donor => b% s1
                  b% s_accretor => b% s2
               else
                  b% d_i = 2
                  b% d_i_old = 2
                  b% a_i = 1
                  b% a_i_old = 1
                  b% s_donor => b% s2
                  b% s_accretor => b% s1
               end if
            end if
         end if
         
      end function extras_binary_finish_step

      subroutine do_saves_for_binary(b, ierr)
         type(binary_info), pointer :: b
         integer, intent(out) :: ierr
         integer :: iounit, id
         character (len=strlen) :: str_photo, filename, iomsg, report_str

         call string_for_model_number('x', b% model_number, b% photo_digits, str_photo)

         filename = trim(trim(b% photo_directory) // '/b_' // str_photo)
         report_str = trim('save ' // filename)
         open(newunit=iounit, file=trim(filename), action='write', &
            status='replace', iostat=ierr, iomsg=iomsg, form='unformatted')
         if (ierr /= 0) then
            write(*,*) 'failed in do_saves_for_binary', trim(filename)
            return
         end if
         call binary_photo_write(b% binary_id, iounit)
         close(iounit)

         if (b% have_star_1) then
            filename = trim(trim(b% s1% photo_directory) // '/1_' // str_photo)
            call star_save_for_restart(b% s1% id, filename, ierr)
            report_str = trim(trim(report_str) // ', ' // filename)
         end if
         if (b% have_star_2) then
            filename = trim(trim(b% s2% photo_directory) // '/2_' // str_photo)
            call star_save_for_restart(b% s2% id, filename, ierr)
            report_str = trim(trim(report_str) // ', ' // filename)
         end if
         if (ierr /= 0) then
            write(*,*) 'failed in do_saves_for_binary'
            return
         end if

         write(*,*) trim(trim(report_str) // ' for model'), b% model_number

         contains

         subroutine binary_photo_write(binary_id, iounit)
            integer, intent(in) :: binary_id, iounit
            type(binary_info), pointer :: b

            integer :: ierr

            ierr = 0
            call binary_ptr(binary_id, b, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in binary_ptr'
               return
            end if

            write(iounit) star_def_version

            write(iounit, iostat=ierr) &
               b% binary_age, b% binary_age_old, &
               b% model_number, b% model_number_old, &
               b% mtransfer_rate, b% mtransfer_rate_old, &
               b% angular_momentum_j, b% angular_momentum_j_old, &
               b% separation, b% separation_old, &
               b% eccentricity, b% eccentricity_old, &
               b% rl_relative_gap(1), b% rl_relative_gap_old(1), &
               b% rl_relative_gap(2), b% rl_relative_gap_old(2), &
               b% r(1), b% r_old(1), &
               b% r(2), b% r_old(2), &
               b% rl(1), b% rl_old(1), &
               b% rl(2), b% rl_old(2), &
               b% m(1), b% m_old(1), &
               b% m(2), b% m_old(2), &
               b% dt, b% dt_old, &
               b% env, b% env_old, &
               b% eq_initial_bh_mass, &
               b% period, b% period_old, &
               b% max_timestep, b% max_timestep_old, &
               b% change_factor, b% change_factor_old, &
               b% min_binary_separation, &
               b% using_jdot_mb(1), b% using_jdot_mb_old(1), &
               b% using_jdot_mb(2), b% using_jdot_mb_old(2), &
               b% d_i, b% d_i_old, b% a_i, b% a_i_old, &
               b% point_mass_i, b% point_mass_i_old, &
               b% ignore_rlof_flag, b% ignore_rlof_flag_old, &
               b% model_twins_flag, b% model_twins_flag_old, &
               b% dt_why_reason, b% dt_why_reason_old, &
               b% have_star_1, b% have_star_2, &
               b% CE_flag, b% CE_flag_old, &
               b% CE_init, b% CE_init_old, &
               b% CE_nz, b% CE_initial_radius, b% CE_initial_separation, b% CE_initial_Mdonor, &
               b% CE_initial_Maccretor, b% CE_initial_age, b% CE_initial_model_number, &
               b% CE_b_initial_age, b% CE_b_initial_model_number, &
               b% CE_num1, b% CE_num1_old, &
               b% CE_num2, b% CE_num2_old, &
               b% CE_lambda1, b% CE_lambda1_old, &
               b% CE_lambda2, b% CE_lambda2_old, &
               b% CE_Ebind1, b% CE_Ebind1_old, &
               b% CE_Ebind2, b% CE_Ebind2_old, &
               b% ixtra(:), b% ixtra_old(:), &
               b% xtra(:), b% xtra_old(:), &
               b% lxtra(:), b% lxtra_old(:)

            if (b% CE_init) then
               write(iounit, iostat=ierr) &
                  b% CE_m(:), b% CE_entropy(:), b% CE_U_in(:), b% CE_U_out(:), b% CE_Omega_in(:), b% CE_Omega_out(:)
            end if

            if (ierr /= 0) stop "error in binary_photo_write"

         end subroutine binary_photo_write

      end subroutine do_saves_for_binary


      logical function switch_donor_star(binary_id) result(do_switch)
         ! HINT: Switching donor will only work by setting
         !       keep_donor_fixed = .false. in &binary_controls
         !       within "inlist"  and adding these lines in
         !       run_binary_extras.f90
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         real(dp) :: F1, q, rho, p, grav, hp, v_th, rl3, q_temp
         real(dp) :: mdot_thin_donor, mdot_thin_accretor

         do_switch = .false.

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         if (b% point_mass_i == 0  .and. .not. b% keep_donor_fixed .and. b% r(b% d_i) < b% rl(b% d_i)) then
            !--------------------- Optically thin MT rate --------------------------
            ! As described in H. Ritter 1988, A&A 202,93-100 and
            ! U. Kolb and H. Ritter 1990, A&A 236,385-392
            rho = b% s_donor% rho(1)  ! density at surface in g/cm^3
            p = b% s_donor% p(1)  ! pressure at surface in dynes/cm^2
            grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2  ! local gravitational
            ! acceleration
            hp = p/(grav*rho)  ! pressure scale height
            v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))
            q = b% m(b% a_i)/b% m(b% d_i)  ! Mass ratio, as defined in Ritter 1988

            ! (Kolb & Ritter 1990 use the opposite!)
            ! consider range of validity for F1, do not extrapolate Eq. A9 of Ritter 1988
            q_temp = min(max(q,0.5d0),10d0)
            F1 = (1.23d0  + 0.5d0* log10(q_temp))
            rl3 = (b% rl(b% d_i))*(b% rl(b% d_i))*(b% rl(b% d_i))
            mdot_thin_donor = (2.0d0*pi/exp(0.5d0)) * v_th*v_th*v_th * &
               rl3/(b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1

            rho = b% s_accretor% rho(1)  ! density at surface in g/cm^3
            p = b% s_accretor% p(1)  ! pressure at surface in dynes/cm^2
            grav = b% s_accretor% cgrav(1)*b% m(b% a_i)/(b% r(b% a_i))**2  ! local gravitational acceleration
            hp = p/(grav*rho) ! pressure scale height
            v_th = sqrt(kerg * b% s_accretor% T(1) / (mp * b% s_accretor% mu(1)))
            q = b% m(b% d_i)/b% m(b% a_i)  ! Mass ratio, as defined in Ritter 1988

            ! (Kolb & Ritter 1990 use the opposite!)
            ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
            q_temp = min(max(q,0.5d0),10d0)
            F1 = (1.23d0  + 0.5d0* log10(q_temp))
            rl3 = (b% rl(b% a_i))*(b% rl(b% a_i))*(b% rl(b% a_i))
            mdot_thin_accretor = (2.0d0*pi/exp(0.5d0)) * v_th*v_th*v_th * &
               rl3/(b% s_accretor% cgrav(1)*b% m(b% a_i)) * rho * F1

            if (abs(mdot_thin_accretor) < abs(mdot_thin_donor)) do_switch = .true.
         end if

      end function switch_donor_star

      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer :: iounit, ios
         character(len=strlen) :: fname, termination_code
         type(star_info), pointer :: s

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         if (b% point_mass_i == 0) then
            if (b% s1% termination_code > 0) then
               termination_code = trim(termination_code_str(b% s1% termination_code))
            else if (b% s2% termination_code > 0) then
               termination_code = trim(termination_code_str(b% s2% termination_code))
            else
               termination_code = 'unknown'
            end if
         else
            if (b% have_star_1) then
               s => b% s1
            else if (b% have_star_2) then
               s => b% s2
            else
               write(*,*) 'logic failure when trying to save termination code'
               return
            end if

            if (s% termination_code > 0) then
               termination_code = trim(termination_code_str(s% termination_code))
            else
               termination_code = 'unknown'
            end if
         end if

         ! run mkdir -p `termination_codes`
         call mkdir('termination_codes')

         if (b% point_mass_i == 0) then
            fname = trim(termination_codes_folder) // '/termination_code_star_plus_star'
         else
            if (do_kicks .and. b% point_mass_i /= 0 .and. add_kick_id_as_suffix) then
               fname = trim(termination_codes_folder) // '/termination_code_' // trim(kick_id)
            else
               fname = trim(termination_codes_folder) // '/termination_code_star_plus_point_mass'
            end if
         end if

         ! write termination code into a file
         open(newunit=iounit, file=trim(fname), iostat=ios)
         if (ios /= 0) return
         write(iounit, fmt='(a)', iostat=ios, advance='no') trim(termination_code)
         if (ios /= 0) return
         
      end subroutine extras_binary_after_evolve


      ! compute radii for outer lagrangian point
      real function eval_outer_roche_lobe(donor_mass, accretor_mass, separation) result(rl2)
         ! Based on Eggleton 2006 'Evolutionary processes in binary and multiple stars'
         real(dp), intent(in) :: donor_mass, accretor_mass, separation
         real(dp) :: q, q13

         q = donor_mass / accretor_mass
         q13 = pow(q,one_third)
         if (q > 1d0) then
            rl2 = 0.49d0 * q13*q13 + 0.15d0
         else
            rl2 = 0.49d0 * q13*q13 + 0.27d0 * q - 0.12d0 * q13*q13*q13*q13
         end if

         rl2 = (separation * rl2) / (0.6d0 * q13*q13 + log1p(q13))

      end function eval_outer_roche_lobe
      
      end module run_binary_extras
