
      module ce_support

      use ce_def

      use const_def
      use star_def
      use star_lib
      use binary_def
      use binary_lib
      use math_lib
      
      implicit none

      contains


      ! set ce vars on startup so they are always initialized
      subroutine variables_on_startup
         
         ce_on = .false.
         ce_off = .true.

         ce_num_star_1 = 0
         ce_num_star_2 = 0

         ce_num_two_stars = 0
         ce_num_xrb = 0

         cumulative_removed_binding_energy = 0d0

      end subroutine variables_on_startup


      ! store initial info
      subroutine store_initial_info(ce_id, ierr)
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         call binary_ptr(ce_id, b, ierr)
         if (ierr /= 0) return

         if (b% s_donor% star_age > b% binary_age) then
            ce_initial_age = b% s_donor% star_age
         else
            ce_initial_age = b% binary_age
         end if

         ce_initial_model_number = b% model_number

         ce_duration = 0d0
         ce_years_in_detach = 0d0

         cumulative_removed_binding_energy = 0d0

         ce_donor_id = b% d_i
         ce_accretor_id = b% a_i

         alpha_mt_start_ce = b% mass_transfer_alpha
         beta_mt_start_ce = b% mass_transfer_beta
         delta_mt_start_ce = b% mass_transfer_delta
         gamma_mt_start_ce = b% mass_transfer_gamma

         fj_start_ce = b% fj
         
         max_explicit_abs_mdot_ce = b% max_explicit_abs_mdot
         max_tries_to_achieve_ce = b% max_tries_to_achieve
         solver_type_start_ce = b% solver_type

         accretor_overflow_terminate_start_ce = b% accretor_overflow_terminate

         lg_mtransfer_rate_start_ce = safe_log10(abs(b% mtransfer_rate*secyer/Msun))

         donor_mass_start_ce = b% m(b% d_i)  ! in g
         donor_radius_start_ce = b% r(b% d_i)  ! in cm

         accretor_mass_start_ce = b% m(b% a_i)  ! in g
         accretor_mass_start_ce = b% r(b% a_i)  ! in g
 
      end subroutine store_initial_info

      
      ! search for type of ce
      integer function binary_type(b, ierr)
         type(binary_info), pointer :: b
         integer, intent(out) :: ierr

         if (b% point_mass_i == 0) then
            binary_type = ce_two_stars
         else
            binary_type = ce_xrb
         end if

      end function binary_type


      ! check ce detachment
      logical function is_detach(b, ierr)
         type(binary_info), pointer :: b
         integer, intent(out) :: ierr
         real(dp) :: tol

         is_detach = .false.
         ierr = 0
         ! avoid checking condition on first model
         if (b% model_number <= ce_initial_model_number) return

         ! get tolerance based on type of ce
         tol = tol_xrb
         if (ce_type == ce_two_stars) tol = tol_two_stars

         ! check detachment
         if (b% rl_relative_gap(ce_donor_id) < -tol) then
            if (ce_years_in_detach >= years_in_detachment) is_detach = .true.
         end if

      end function is_detach


      ! check for ce merger
      logical function will_merge(b, ierr)
         type(binary_info), pointer :: b
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: h_diff, he_diff
         real(dp) :: rlobe

         will_merge = .false.
         ierr = 0
         ! do not check if first timestep
         if (b% model_number <= ce_initial_model_number) return

         ! if core is overflowing, end as merge
         do k = 1, b% s_donor% nz
            h_diff = abs(b% s_donor% center_h1 - b% s_donor% X(k))
            he_diff = abs(b% s_donor% center_he4 - b% s_donor% Y(k))
            if (h_diff < 0.01 .and. he_diff < 0.01) then
               rlobe = binary_eval_rlobe(b% s_donor% m(k), b% m(b% a_i), b% separation)
               if (b% s_donor% r(k) > rlobe) then
                  write(*,'(a)') 'core of donor star is overflowing its Roche lobe'
                  will_merge = .true.
                  return
               end if
               exit
            end if
         end do

         ! check for many retries during ce
         if (b% s_donor% num_retries > max_number_retries_during_ce) then
            write(*,'(a)') 'reach max_number_retries_during_ce'
            will_merge = .true.
            return
         end if

         ! if exceeds max relative gap, also merge
         if (b% rl_relative_gap(ce_donor_id) > max_relative_gap) then
            write(*,'(a)') 'donor star exceeds max_relative_gap value'
            will_merge = .true.
            return
         end if

         ! if accretor inside donor star, merge
         if (b% separation < b% r(ce_accretor_id)) then
            write(*,'(a)') 'non-degenerate accretor is inside donor star'
            will_merge = .true.
            return
         end if

      end function will_merge

      end module ce_support

