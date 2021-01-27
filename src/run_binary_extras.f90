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

      logical :: high_mass_evolution, low_mass_evolution
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
         high_mass_evolution, low_mass_evolution, &
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

      ! post-kick id (changed from run to run)
      character(len=strlen) :: kick_id
      real(dp) :: v_kick, phi_kick, theta_kick
      real(dp) :: ecc_kick, porb_kick

      ! utilities
      integer :: num_switches
      logical :: second_collapse

      logical, parameter :: dbg = .true.
      
      contains


      include 'test_suite_extras.inc'


      subroutine do_run(extras_controls, extras_binary_controls, ierr)
         integer, intent(out) :: ierr
         character(len=strlen) :: inlist_fname, termination_code_fname
         logical :: has_reached_cc

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

         ! get bin2dco controls
         inlist_fname = 'inlist_bin2dco'
         call read_bin2dco_controls(inlist_fname, ierr)
         if (ierr /= 0) return

         ! just check that `*_mass_evolution` are not both false
         if (.not. high_mass_evolution .and. .not. low_mass_evolution) then
            write(*,'(/,a,/)') &
               'cannot have both `high_mass_evolution` and `low_mass_evolution` be .false.'
            return
         end if

         ! star + star simulation (if chosen)
         if (do_star_plus_star) then
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
         ! assymetric kicks from a file after which binary parameters will be updated,
         ! or a typical MESAbinary evolution in which binary parameters needs to be
         ! set in the file `inlist_star_plus_pm`
         if (do_kicks) then

            ! get masses from files produced by core_collapse module
            call read_cc_parameters(star_info_at_cc_filename, binary_info_at_cc_filename, &
               mass_of_progenitor, mass_of_remnant, mass_of_companion, &
               pre_cc_separation, pre_cc_period, &
               ierr)

            write(*,'(a)')
            write(*,'(a)') 'binary at core-collapse:'
            write(*,'(a)')
            write(*,1) 'progenitor_mass', mass_of_progenitor
            write(*,1) 'remnant_mass', mass_of_remnant
            write(*,1) 'companion_mass', mass_of_companion
            write(*,1) 'separation', pre_cc_separation
            write(*,1) 'period', pre_cc_period
            write(*,'(a)')

            if (do_kicks_in_one_run) then
               write(*,'(a)') 'TBD'
            else
               write(*,'(a)') 'TBD'
            end if

         else
            
            call run1_binary(.true., extras_controls, extras_binary_controls, &
               ierr, &
               star_plus_pm_filename)

         end if

         if (dbg) stop 'debugging'

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

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

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
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         if (.not. restart) then
            b% lxtra(lx_ce_on) = .false.
            b% lxtra(lx_ce_off) = .true.
            b% ixtra(ix_ce_model_number) = 0

            ! set initial ce variables, no need to check if it's a restart
            if (b% point_mass_i == 0) then
               call ce_variables_on_startup(ce1_inlist_filename, ierr)
               if (ierr /= 0) return
            else
               call ce_variables_on_startup(ce2_inlist_filename, ierr)
               if (ierr /= 0) return
            end if
         else
            if (b% lxtra(lx_ce_on)) then
               ce_on = .true.
               ce_off = .false.
            end if
            if (b% lxtra(lx_ce_off)) then
               ce_on = .false.
               ce_off = .true.
            end if

            ! set ce controls
            if (b% point_mass_i == 0) then
               call ce_set_controls(ce1_inlist_filename, ierr)
               if (ierr /= 0) return
            else
               call ce_set_controls(ce2_inlist_filename, ierr)
               if (ierr /= 0) return
            end if

            if (ce_on) then
               call ce_init_binary_controls(binary_id, ierr)
               call ce_donor_binding_energy_prev_step(binary_id, ierr)
               call ce_orbital_energy_prev_step(binary_id, ierr)
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

         ! check ce end
         if (ce_on) then

            !b% s_donor% use_gold2_tolerances = .false.
            !b% s_donor% use_gold_tolerances = .false.
            !b% s_donor% use_dedt_form_of_energy_eqn = .false.
            !b% s_donor% always_use_dedt_form_of_energy_eqn = .false.
            !b% s_donor% use_dedt_form_with_total_energy_conservation = .false.
            !b% s_donor% varcontrol_target = 1d-4
            !b% s_donor% Pextra_factor = 4
            !b% s_donor% mesh_delta_coeff = 0.5d0
            !b% s_donor% kap_rq% kap_blend_logT_lower_bdy = 4.1
            !b% s_donor% kap_rq% kap_blend_logT_upper_bdy = 4.2
            !b% s_donor% tau_base = 1d0
            !b% s_donor% limit_for_rel_error_in_energy_conservation = 1d-2
            !b% s_donor% hard_limit_for_rel_error_in_energy_conservation = 1d0

            !b% s_donor% mix_factor = 0d0
            !b% s_donor% dxdt_nuc_factor = 0d0

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
               b% lxtra(lx_ce_on) = .true.
               b% lxtra(lx_ce_off) = .false.
               b% ixtra(ix_ce_model_number) = ce_initial_model_number
               call ce_init(binary_id, ierr)
               if (ierr /= 0) return
               write(*,11) 'save model at start of ce, model_number', b% model_number
               call do_saves_for_binary(b, ierr)
            end if
         end if

      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         logical :: do_switch
         integer :: star_cc_id
         integer :: number_io
         real(dp), parameter :: chandra_mass = 1.4d0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  

         extras_binary_finish_step = keep_going

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
         end if

         ! for high-mass evolution
         if (high_mass_evolution) then
            ! first collapse
            star_cc_id = 0
            if (b% point_mass_i == 0) then
               if (b% s1% center_c12 < 1d-3 * b% s1% initial_z .and. b% s1% center_he4 < 1d-6) then
                  star_cc_id = 1
                  call star_write_model(2, 'companion_at_core_collapse.mod', ierr)
                  b% s1% termination_code = t_xtra1
                  termination_code_str(t_xtra1) = 'core-collapse'
                  extras_binary_finish_step = terminate
               else if (b% s2% center_c12 < 1d-3 * b% s2% initial_z .and. b% s2% center_he4 < 1d-6) then
                  star_cc_id = 2
                  call star_write_model(1, 'companion_at_core_collapse.mod', ierr)
                  b% s2% termination_code = t_xtra1
                  termination_code_str(t_xtra1) = 'core-collapse'
                  extras_binary_finish_step = terminate
               end if
            ! second collapse
            else if (b% point_mass_i == 1) then
               if (b% s2% center_c12 < 1d-3 * b% s2% initial_z .and. b% s2% center_he4 < 1d-6) then
                  second_collapse = .true.
                  star_cc_id = 2
                  b% s2% termination_code = t_xtra1
                  termination_code_str(t_xtra1) = 'core-collapse'
                  extras_binary_finish_step = terminate
               end if
            else if (b% point_mass_i == 2) then
               if (b% s1% center_c12 < 1d-3 * b% s1% initial_z .and. b% s1% center_he4 < 1d-6) then
                  second_collapse = .true.
                  star_cc_id = 1
                  b% s1% termination_code = t_xtra1
                  termination_code_str(t_xtra1) = 'core-collapse'
                  extras_binary_finish_step = terminate
               end if
            end if

            if (star_cc_id > 0 .and. .not. second_collapse) then
               call cc_set_controls(cc1_inlist_filename, ierr)
               if (do_kicks) then
                  write(filename_for_star_data, '(a)') &
                     trim(filename_for_star_data) // "_1"
                  write(filename_for_binary_data, '(a)') &
                     trim(filename_for_binary_data) // "_1"
               end if
               if (ierr /= 0) return
               call cc_compact_object_formation(star_cc_id, ierr)
               if (ierr /= 0) then
                  write(*,'(a)') 'failed in cc_compact_object_formation'
                  return
               end if
               return
            else if (star_cc_id > 0 .and. second_collapse) then
               write(*,'(a)') 'calling second collapse'
               if (do_kicks) then
                  call cc_set_controls(cc2_inlist_filename, ierr)
                  if (ierr /= 0) return
                  ! replace filenames by adding natal-kick id
                  write(filename_for_star_data, '(a)') &
                     trim(filename_for_star_data) // "_" // trim(kick_id)
                  write(filename_for_binary_data, '(a)') &
                     trim(filename_for_binary_data) // "_" // trim(kick_id)
                  call cc_compact_object_formation(star_cc_id, ierr)
                  if (dbg) then
                     write(*,'(a32, 2x, a)') &
                        'filename_for_star_data =', trim(filename_for_star_data)
                     write(*,'(a32, 2x, a)') &
                        'filename_for_binary_data =', trim(filename_for_binary_data)
                  end if
               else
                  call cc_set_controls(cc1_inlist_filename, ierr)
                  if (do_kicks) then
                     write(filename_for_star_data, '(a)') &
                        trim(filename_for_star_data) // "_1"
                     write(filename_for_binary_data, '(a)') &
                        trim(filename_for_binary_data) // "_1"
                  end if
                  if (ierr /= 0) return
                  call cc_compact_object_formation(star_cc_id, ierr)
               end if
               if (ierr /= 0) then
                  write(*,'(a)') 'failed in cc_compact_object_formation'
                  return
               end if
               return
            end if
         end if

         ! check if we need to switch donor star only when ce is off
         if (ce_off) then
            do_switch = switch_donor_star(binary_id)
            if (do_switch) then
               write(*,'(a)') 'switching donor'
               num_switches = num_switches + 1
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
            if (b% s1% termination_code > 0) then
               termination_code = trim(termination_code_str(b% s1% termination_code))
            else if (b% s2% termination_code > 0) then
               termination_code = trim(termination_code_str(b% s2% termination_code))
            else
               termination_code = 'unknown'
            end if
         end if

         ! run mkdir -p `termination_codes`
         call mkdir('termination_codes')

         if (b% point_mass_i == 0) then
            fname = trim(termination_codes_folder) // '/termination_code_star_plus_star'
         else
            fname = trim(termination_codes_folder) // '/termination_code_' // trim(kick_id)
         end if

         ! write termination code into a file
         open(newunit=iounit, file=trim(fname), iostat=ios)
         if (ios /= 0) return
         write(iounit, fmt='(a)', iostat=ios, advance='no') trim(termination_code)
         if (ios /= 0) return
         
      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
