
      module ce_lib
      
      use ce_def

      use star_def
      use star_lib
      use binary_def
      use binary_lib

      use math_lib
      
      implicit none

      contains


      ! load ce variables on startup
      subroutine ce_variables_on_startup(inlist_fname, ierr)
         use ce_support, only: variables_on_startup
         character(len=*), intent(in) :: inlist_fname
         integer, intent(out) :: ierr

         ierr = 0

         ! set ce flags
         call variables_on_startup

         ! set ce controls
         call ce_set_controls(inlist_fname, ierr)
         if (ierr /= 0) return

      end subroutine ce_variables_on_startup


      ! check condition for ce init
      subroutine ce_unstable_mt_phase(ce_id, ierr)
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         real(dp) :: mdot_edd_donor, mdot_edd_accretor
         character(len=strlen) :: ce_condition

         ierr = 0
         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         ! Eddington mass-changes rates
         mdot_edd_donor = (4d0 * pi * clight * b% s_donor% r(1)) &
            / b% s_donor% photosphere_opacity

         if (b% point_mass_i == 0) then
            mdot_edd_accretor = (4d0 * pi * clight * b% s_accretor% r(1)) &
               / b% s_accretor% photosphere_opacity
         else
            mdot_edd_accretor = 4.0d0 * pi * standard_cgrav * b% m(b% a_i) &
            / (clight * b% s_donor% opacity(1) * b% mdot_edd_eta)
         end if

         ! check for ce init
         ce_condition = 'stable mt'
         if (b% point_mass_i == 0) then
            if (b% m(b% d_i) > b% m(b% a_i) .and. abs(b% mtransfer_rate) > edd_scaling_factor * mdot_edd_donor) then
               ce_condition = 'mtransfer_rate > mdot_edd_donor'
            else if (b% m(b% d_i) > b% m(b% a_i) .and. abs(b% component_mdot(b% a_i)) > mdot_edd_accretor) then
               ce_condition = 'maccreted_rate > mdot_edd_accretor'
            else if (b% m(b% d_i) > b% m(b% a_i) .and. abs(b% mtransfer_rate) * secyer/Msun > max_mdot_rlof) then
               ce_condition = 'mtransfer_rate > max_mdot_rlof'
            end if
         else
            if (b% m(b% d_i) > b% m(b% a_i) .and. abs(b% mtransfer_rate) > edd_scaling_factor * mdot_edd_donor) then
               ce_condition = 'mtransfer_rate > mdot_edd_donor'
            else if (b% m(b% d_i) > b% m(b% a_i) .and. abs(b% mtransfer_rate) * secyer/Msun > max_mdot_rlof) then
               ce_condition = 'mtransfer_rate > max_mdot_rlof'
            end if
         end if

         if (report_mdot_for_ce) then
            call report_mdot_comparison
         end if

         if (ce_condition /= 'stable mt') then
            write(*,'(/,a)') 'turning ce on, condition: ' // trim(ce_condition)
            ce_on = .true.
            ce_off = .false.
         end if

         contains

         subroutine report_mdot_comparison
            real(dp) :: m_conv, f_conv
            type(star_info), pointer :: s
            integer :: k

            s => b% s_donor

            if (s% he_core_mass > 0d0 .and. s% center_h1 < 1d-6) then
               m_conv = 0d0
               if (s% he_core_k == 1 .or. (s% star_mass - s% he_core_mass) < 1d-3) then
                  f_conv = 0d0
               else
                  do k = 1, s% he_core_k
                     if (s% mixing_type(k) == convective_mixing) m_conv = m_conv + s% dm(k)
                  end do
                  f_conv = (m_conv / Msun) / (s% star_mass - s% he_core_mass)
               end if
            else
               m_conv = 0d0
               f_conv = 0d0
            end if

            if (b% point_mass_i == 0) then
               if (b% num_tries == 1) then
                  write(*,'(a)') '  binary_step  mtransfer_rate  mdot_edd_donor' // &
                     '      f_conv  mdot_edd_accretor  mdot_accretor  max_mdot_rlof                              result'
                  write(*,'(a)') '  -------------------------------------------' // &
                     '-------------------------------------------------------------------------------------------------'
               end if
               write(*,'(i13,2(1pe16.3,0p),1(1pe12.3,0p),1(1pe19.3,0p),2(1pe15.3,0p),a36)') b% s_donor% model_number, &
                  abs(b% mtransfer_rate) * secyer/Msun, &
                  mdot_edd_donor * edd_scaling_factor * secyer/Msun, &
                  f_conv, &
                  mdot_edd_accretor * secyer/Msun, &
                  abs(b% component_mdot(b% a_i)) * secyer/Msun, &
                  max_mdot_rlof, &
                  trim(ce_condition)
            else
               if (b% num_tries == 1) then
                  write(*,'(a)') '  binary_step  mtransfer_rate  mdot_edd_donor' // &
                     '          f_conv   max_mdot_rlof                              result'
                  write(*,'(a)') '  -------------------------------------------' // &
                     '--------------------------------------------------------------------'
               end if
               write(*,'(i13,4(1pe16.3,0p),a36)') b% s_donor% model_number, &
                  abs(b% mtransfer_rate) * secyer/Msun, &
                  mdot_edd_donor * edd_scaling_factor * secyer/Msun, &
                  f_conv, &
                  max_mdot_rlof, &
                  trim(ce_condition)
            end if
         end subroutine report_mdot_comparison


      end subroutine ce_unstable_mt_phase


      ! ce initialize
      subroutine ce_init(ce_id, ierr)
         use ce_ctrls_io, only: save_ce_profile, save_ce_model
         use utils_lib, only: mkdir
         integer, intent(in) :: ce_id

         integer, intent(out) :: ierr

         ierr = 0
         
         ! run `mkdir -p` over ce_data_directory if it does not exist
         if (save_profile_pre_ce .or. save_profile_after_ce) call mkdir(ce_data_directory)

         ! check type of ce
         call ce_type_of(ce_id, ierr)
         if (ierr /= 0) return

         ! get number of ce events
         call ce_num_events(ce_id, ierr)
         if (ierr /= 0) return

         ! store some info from start of ce
         call ce_store_initial_info(ce_id, ierr)
         if (ierr /= 0) return

         ! check if we need to save stuff
         if (save_profile_pre_ce) call save_ce_profile(ce_id, 'pre', ierr)
         if (ierr /= 0) return
         ! if (save_model_pre_ce) call save_ce_model(ce_id, 'pre', ierr)
         ! if (ierr /= 0) return

         ! eval energies from previous step
         call ce_donor_binding_energy_prev_step(ce_id, ierr)
         if (ierr /= 0) return
         call ce_orbital_energy_prev_step(ce_id, ierr)
         if (ierr /= 0) return

         ! modify binary controls
         call ce_init_binary_controls(ce_id, ierr)
         if (ierr /= 0) return

      end subroutine ce_init


      ! store some info from binary at start of ce
      subroutine ce_store_initial_info(ce_id, ierr)
         use ce_support, only: store_initial_info
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr

         ierr = 0

         call store_initial_info(ce_id, ierr)
         if (ierr /= 0) return

      end subroutine ce_store_initial_info

      
      ! ce controls
      subroutine ce_set_controls(inlist_fname, ierr)
         use ce_ctrls_io, only: set_default_controls, read_ce_controls
         character(len=*), intent(in) :: inlist_fname
         integer, intent(out) :: ierr

         ierr = 0

         call set_default_controls

         call read_ce_controls(inlist_fname, ierr)

      end subroutine ce_set_controls


      ! type of ce
      subroutine ce_type_of(ce_id, ierr)
         use ce_support, only: binary_type
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         logical :: val

         include 'formats'

         ierr = 0
         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         ce_type = binary_type(b, ierr)
         if (ierr /= 0) return

         if (ce_dbg) write(*,11) 'type of ce', ce_type

      end subroutine ce_type_of


      ! number of ce events
      subroutine ce_num_events(ce_id, ierr)
         use ce_utils, only: num_events
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr

         ierr = 0
         call num_events(ce_id, ierr)
         if (ierr /= 0) return

      end subroutine ce_num_events
         

      ! change controls of binary module
      subroutine ce_init_binary_controls(ce_id, ierr)
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         ierr = 0
         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         ! reduce timestep to be half a year max
         if (b% time_step > 0.5d0) then
            b% time_step = 0.5d0
         end if
         b% s_donor% dt = b% time_step * secyer
         if (ce_type == ce_two_stars) b% s_accretor% dt = b% s_donor% dt

         ! do not consider accretion during ce
         if (.not. add_accretion_on_ce) then
            b% mass_transfer_alpha = 0.0d0
            b% mass_transfer_beta = 1.0d0
            b% mass_transfer_delta = 0.0d0
            b% mass_transfer_gamma = 0.0d0
         else
            write(*,'(a)') 'this feature is not ready to use, so no accretion will be calculated'
            b% mass_transfer_alpha = 0.0d0
            b% mass_transfer_beta = 1.0d0
            b% mass_transfer_delta = 0.0d0
            b% mass_transfer_gamma = 0.0d0
         end if

         ! no photos during ce
         b% photo_interval = 100000

         ! use custom rlo_mdot and jdot
         b% use_other_rlo_mdot = .true.
         b% other_rlo_mdot => ce_rlo_mdot
         b% use_other_extra_jdot = .true.
         b% other_extra_jdot => ce_dot_j

         ! do not use default jdots
         b% do_jdot_gr = .false.
         b% do_jdot_ml = .false.
         b% do_jdot_ls = .false.
         b% do_jdot_missing_wind = .false.
         b% do_jdot_mb = .false.

         b% fj = 0.001d0
         b% fj_hard = 0.01d0

         b% max_explicit_abs_mdot = max_mdot_rlof
         b% max_tries_to_achieve = 0
         b% solver_type = 'bisect'

         !b% s_donor% delta_lgTeff_limit = 0.1d0
         !b% s_donor% delta_lgL_phot_limit_L_min = 1d99
         !b% s_donor% delta_lgL_limit_L_min = 1d99
         
         ! avoid burning info on donor
         b% s_donor% mix_factor = 0d0
         b% s_donor% dxdt_nuc_factor = 0d0

         ! only used in two stars to avoid accretor overflow
         b% accretor_overflow_terminate = 0d0

         ! instant circularization
         b% eccentricity = 0d0
         b% eccentricity_old = b% eccentricity

         ! do not change donor during ce
         if (ce_type == ce_two_stars) b% keep_donor_fixed = .true.

      end subroutine ce_init_binary_controls


      ! eval binding energy of star
      subroutine ce_donor_binding_energy_prev_step(ce_id, ierr)
         use ce_utils, only: eval_binding_energy
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         type(star_info), pointer :: s

         ierr = 0
         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         s => b% s_donor

         call eval_binding_energy(s, s% nz, 1, donor_grav_energy_prev_step, donor_int_energy_prev_step, ierr)
         if (ierr /= 0) return

         donor_bind_energy_prev_step = donor_grav_energy_prev_step + donor_int_energy_prev_step

      end subroutine ce_donor_binding_energy_prev_step


      ! orbital energy of binary
      subroutine ce_orbital_energy_prev_step(ce_id, ierr)
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         ierr = 0
         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         orbital_energy_prev_step = - 0.5d0 * standard_cgrav * b% m(b% d_i) * b% m(b% a_i) / b% separation

      end subroutine ce_orbital_energy_prev_step


      ! ce mdot
      subroutine ce_rlo_mdot(ce_id, mdot, ierr)
         use ce_mdot, only: eval_ce_mdot
         integer, intent(in) :: ce_id
         real(dp), intent(out) :: mdot
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         ierr = 0
         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         mdot = eval_ce_mdot(b, ierr)
         if (ierr /= 0) return

      end subroutine ce_rlo_mdot


      ! ce jdot
      subroutine ce_dot_j(ce_id, ierr)
         use ce_jdot, only: eval_ce_jdot
         integer, intent(in) :: ce_id
        integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         ierr = 0
         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         call eval_ce_jdot(b, ierr)
         if (ierr /= 0) return

      end subroutine ce_dot_j


      ! check condition of ce at end of each timestep
      subroutine ce_check_state(ce_id, ierr)
         use ce_support, only: is_detach, will_merge
         use ce_ctrls_io, only: save_ce_profile, save_ce_model
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         ierr = 0
         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         if (b% model_number <= ce_initial_model_number) then
            ce_detach = .false.
            ce_merge = .false.
            return
         end if

         ! check for detachment
         ce_detach = is_detach(b, ierr)
         if (ierr /= 0) return

         ! check for ce merger
         if (.not. ce_detach) ce_merge = will_merge(b, ierr)
         if (ierr /= 0) return

         if (ce_detach) then
            write(*,'(/,a,/)') 'reach ce detach'
            if (save_profile_after_ce) call save_ce_profile(ce_id, 'after', ierr)
           ! if (save_model_after_ce) call save_ce_model(ce_id, 'after', ierr)
            call ce_end(b, ierr)
            ce_on = .false.
            ce_off = .true.
         end if

         if (ce_merge) then
            b% s_donor% termination_code = t_xtra1
            termination_code_str(t_xtra1) = 'ce merge'
            return
         end if

      end subroutine ce_check_state


      ! call this when step is already accepted, i.e., inside the
      ! finish_step
      subroutine ce_step_success(ce_id, ierr)
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         include 'formats'

         ierr = 0
         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         ce_duration = ce_duration + b% time_step

         if (ce_duration < years_to_max_mdot_rlof .and. abs(b% mtransfer_rate) * secyer/Msun < max_mdot_rlof) then
            if (b% s_donor% dt_next > 0.5d0 * secyer) then
               write(*,1) 'reduce dt due to not reaching years_to_max_mdot_rlof', &
                  b% s_donor% dt_next / secyer, &
                  min(0.5d0 * secyer, b% s_donor% dt_next) / secyer
               b% s_donor% dt_next = min(0.5d0 * secyer, b% s_donor% dt_next)
               if (ce_type == ce_two_stars) b% s_accretor% dt_next = b% s_donor% dt_next
            end if
         end if

         ! start to count the time a binary is detached but only after reaching max_mdot_rlof
         ! after a certain amount of time, just exit CE
         if (b% r(b% d_i) < b% rl(b% d_i) .and. ce_duration > years_to_max_mdot_rlof) then
            ce_years_in_detach = ce_years_in_detach + b% time_step
            write(*,1) 'binary close ce end. years detached:', ce_years_in_detach
         end if

         ! add energy removed to cumulative
         cumulative_removed_binding_energy = cumulative_removed_binding_energy &
            + (donor_bind_energy_prev_step - donor_bind_energy)

         ! update prev step binding energy for next timestep to take
         donor_bind_energy_prev_step = donor_bind_energy
         ! update prev step orbital energy for next timestep to take
         call ce_orbital_energy_prev_step(ce_id, ierr)

      end subroutine ce_step_success


      ! ce end
      subroutine ce_end(b, ierr)
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         ierr = 0
         
         b% photo_interval = 100

         ! unset some binary controls
         b% use_other_rlo_mdot = .false.
         b% use_other_extra_jdot = .false.

         b% do_jdot_gr = .false.
         b% do_jdot_ml = .true.
         b% do_jdot_ls = .false.
         b% do_jdot_missing_wind = .false.
         b% do_jdot_mb = .false.

         b% mass_transfer_alpha = 0d0
         b% mass_transfer_beta = 0d0
         b% mass_transfer_delta = 0d0
         b% mass_transfer_gamma = 0d0

         b% fj = 0.001d0
         b% fj_hard = 0.01d0

         b% max_explicit_abs_mdot = 1d99
         b% max_tries_to_achieve = 200
         b% solver_type = 'both'
         
         b% s_donor% mix_factor = 1d0
         b% s_donor% dxdt_nuc_factor = 1d0

         b% accretor_overflow_terminate = 10d0

         if (ce_type == ce_two_stars) b% keep_donor_fixed = .false.

      end subroutine ce_end

      end module ce_lib

