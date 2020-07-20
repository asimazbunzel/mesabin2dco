

      module bin2dco_io

      implicit none

      contains


      ! write inlist pointing to new new post core-collapse binary
      subroutine create_inlist(inlist_fname, id, m2, period, e, ierr)
         use const_def, only: dp
         use utils_lib, only: alloc_iounit, free_iounit
         character(len=*) :: inlist_fname
         character(len=*) :: id
         real(dp), intent(in) :: m2
         real(dp), intent(in) :: period
         real(dp), intent(in) :: e
         integer, intent(out) :: ierr
         integer :: iounit

         1 format(3x, a, 1pd10.3)
         2 format(3x, a)

         ! allocate file unit
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit to create_inlist'
         open(unit=iounit, file=trim(inlist_fname), status='replace', action='write', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(inlist_fname)
            return
         end if

         ! replace inlist with new one
         write(iounit,'(a)')
         write(iounit,'(a)') '&binary_controls'
         write(iounit,'(a)')
         write(iounit,2) "history_name = 'binary_" // trim(id) // ".data'"
         write(iounit,'(a)')
         write(iounit,1) 'm2 =', m2
         write(iounit,'(a)')
         write(iounit,1) 'initial_period_in_days =', period
         write(iounit,'(a)')
         write(iounit, 1) 'initial_eccentricity =', e
         write(iounit,'(a)')
         write(iounit,'(a)') '/'
         write(iounit,'(a)')
         write(iounit,'(a)') '&controls'
         write(iounit,'(a)')
         write(iounit,2) "star_history_name = 'secondary_history_" // trim(id) // ".data'"
         write(iounit,'(a)')
         write(iounit,'(a)') '/'

         close(iounit)
         call free_iounit(iounit)

      end subroutine create_inlist


      ! modify inlist_cc for the second collapse using information from first collapse inlist
      subroutine modify_inlist_2nd_cc(inlist_fname_cc_1, inlist_fname_cc_2, id, ierr)
         use const_def, only: dp, strlen
         use utils_lib, only: alloc_iounit, free_iounit
         ! copy from core_collapse/private/cc_ctrls.io.f90
         character(len=strlen) :: model_name
         character (len=strlen) :: filename_for_binary_data
         character (len=strlen) :: filename_for_star_data
         real(dp) :: max_ns_mass
         logical :: continue_binary_evolution
         logical :: add_asymmetric_kick
         namelist /cc_controls/ &
            model_name, &
            max_ns_mass, &
            filename_for_star_data, &
            filename_for_binary_data, &
            continue_binary_evolution, &
            add_asymmetric_kick
         character(len=*) :: inlist_fname_cc_1, inlist_fname_cc_2
         character(len=*) :: id
         integer, intent(out) :: ierr
         integer :: iounit

         1 format(3x, a, 1pd10.3)
         2 format(3x, a)

         ! allocate file unit
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit to modify_inlist_2nd_cc'
         open(unit=iounit, file=trim(inlist_fname_cc_1), action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(inlist_fname_cc_1)
            call free_iounit(iounit)
            return
         end if
         
         ! read cc_controllers namelist
         read(iounit, nml=cc_controls, iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to read cc_controls'
            return
         end if
         close(iounit)
         call free_iounit(iounit)

         ! allocate file unit
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit to write inlist_2nd_cc'
         open(unit=iounit, file=trim(inlist_fname_cc_2), status='replace', action='write', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(inlist_fname_cc_2)
            return
         end if

         ! write cc controls to file
         write(iounit,'(a)')
         write(iounit,'(a)') '&cc_controls'
         write(iounit,'(a)')
         write(iounit,2) "model_name = '" // trim(model_name) // "'"
         write(iounit,'(a)')
         write(iounit,2) "filename_for_star_data = 'cc_data/core_collapse_" // trim(id) // ".data'"
         write(iounit,2) "filename_for_binary_data = 'cc_data/core_collapse_" // trim(id) // ".data'"
         write(iounit,'(a)')
         write(iounit,1) 'max_ns_mass = ', max_ns_mass
         write(iounit,'(a)')
         write(iounit,2) 'add_asymmetric_kick = .false.'
         write(iounit,'(a)')
         write(iounit,2) 'continue_binary_evolution = .false.'
         write(iounit,'(a)')
         write(iounit,'(a)') '/ ! end of cc_controls namelist'

         close(iounit)
         call free_iounit(iounit)

      end subroutine modify_inlist_2nd_cc


      ! modify inlist for ce of second part using info from first part ce inlist
      subroutine modify_inlist_2nd_ce(inlist_fname_ce_1, inlist_fname_ce_2, id, ierr)
         use const_def, only: dp, strlen
         use utils_lib, only: alloc_iounit, free_iounit
         ! copy from ce/private/ce_ctrls.io.f90
         real(dp) :: alpha_two_stars, alpha_xrb, edd_scaling_factor, tol_two_stars, tol_xrb
         real(dp) :: max_mdot_rlof, years_to_max_mdot_rlof, max_relative_gap
         integer :: max_number_retries_during_ce
         logical :: add_accretion_on_ce, save_profile_pre_ce, save_profile_after_ce
         logical :: save_model_pre_ce, save_model_after_ce
         character(len=strlen) :: filename_donor_profile_pre_ce, filename_donor_profile_after_ce
         character(len=strlen) :: filename_donor_model_pre_ce, filename_donor_model_after_ce
         character(len=strlen) :: filename_accretor_profile_pre_ce, filename_accretor_profile_after_ce
         character(len=strlen) :: filename_accretor_model_pre_ce, filename_accretor_model_after_ce
         namelist /ce_controls/ &
            alpha_two_stars, &
            alpha_xrb, &
            edd_scaling_factor, &
            tol_two_stars, &
            tol_xrb, &
            years_to_max_mdot_rlof, &
            max_mdot_rlof, &
            add_accretion_on_ce, &
            save_profile_pre_ce, &
            save_profile_after_ce, &
            save_model_pre_ce, &
            save_model_after_ce, &
            filename_donor_profile_pre_ce, &
            filename_donor_profile_after_ce, &
            filename_donor_model_pre_ce, &
            filename_donor_model_after_ce, &
            filename_accretor_profile_pre_ce, &
            filename_accretor_profile_after_ce, &
            filename_accretor_model_pre_ce, &
            filename_accretor_model_after_ce, &
            max_relative_gap, &
            max_number_retries_during_ce
         character(len=*) :: inlist_fname_ce_1, inlist_fname_ce_2
         character(len=*) :: id
         integer, intent(out) :: ierr
         integer :: iounit
         character(len=strlen) :: tmp1, tmp2, tmp3, tmp4

         1 format(3x, a, 1pd10.3)
         2 format(3x, a)
         3 format(3x, a, i3)

         ! open first ce inlist
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit to modify_inlist_2nd_ce'
         open(unit=iounit, file=trim(inlist_fname_ce_1), action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(inlist_fname_ce_1)
            call free_iounit(iounit)
            return
         end if

         ! read ce_controls namelist
         read(iounit, nml=ce_controls, iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to read ' // trim(inlist_fname_ce_1)
            return
         end if
         close(iounit)
         call free_iounit(iounit)

         ! write namelist into new file for second collapse (changes for each star + point_mass simulation)
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit to write inlist_2nd_ce'
         open(unit=iounit, file=trim(inlist_fname_ce_2), status='replace', action='write', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(inlist_fname_ce_2)
            call free_iounit(iounit)
            return
         end if

         write(iounit,'(a)')
         write(iounit,'(a)') '&ce_controls'
         write(iounit,'(a)')
         write(iounit,1) 'alpha_two_stars =', alpha_two_stars
         write(iounit,1) 'alpha_xrb =', alpha_xrb
         write(iounit,'(a)')
         write(iounit,1) 'edd_scaling_factor =', edd_scaling_factor
         write(iounit,'(a)')
         write(iounit,1) 'tol_two_stars =', tol_two_stars
         write(iounit,1) 'tol_xrb =', tol_xrb
         write(iounit,'(a)')
         write(iounit,1) 'years_to_max_mdot_rlof =', years_to_max_mdot_rlof
         write(iounit,'(a)')
         write(iounit,1) 'max_mdot_rlof =', max_mdot_rlof
         write(iounit,'(a)')
         write(iounit,2) 'add_accretion_on_ce = .false.'
         write(iounit,'(a)')
         write(iounit,2) 'save_profile_pre_ce = .true.'
         write(iounit,2) 'save_profile_after_ce = .true.'
         write(iounit,2) 'save_model_pre_ce = .false.'
         write(iounit,2) 'save_model_after_ce = .false.'
         write(iounit,'(a)')
         write(iounit,2) "filename_donor_profile_pre_ce = 'ce_data/profile_donor_pre_ce_" // trim(id) // "'"
         write(iounit,2) "filename_donor_profile_after_ce = 'ce_data/profile_donor_after_ce_" // trim(id) // "'"
         write(iounit,2) "filename_donor_model_pre_ce = 'ce_data/model_donor_pre_ce_" // trim(id) // "'"
         write(iounit,2) "filename_donor_model_after_ce = 'ce_data/model_donor_after_ce_" // trim(id) // "'"
         write(iounit,'(a)')
         write(iounit,1) 'max_relative_gap =', max_relative_gap
         write(iounit,'(a)')
         write(iounit,3) 'max_number_retries_during_ce = ', max_number_retries_during_ce
         write(iounit,'(a)')
         write(iounit,'(a)') '/ ! end of ce_controls namelist'

         close(iounit)
         call free_iounit(iounit)

      end subroutine modify_inlist_2nd_ce


      ! save termination_code to its final place
      subroutine store_termination_code(name_id)
         use const_def, only: strlen
         use utils_lib, only: alloc_iounit, free_iounit
         character(len=*) :: name_id
         character(len=strlen) :: second_part_end_code
         integer :: iounit, ierr

         second_part_end_code = ''
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit for termination code tmp_second_part'
         open(unit=iounit, file='tmp_second_part', action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            call free_iounit(iounit)
            write(*,*) 'failed to open tmp_second_part'
            call execute_command_line('rm -f tmp_second_part')
            iounit = alloc_iounit(ierr)
            if (ierr /= 0) stop 'could not alloc_iounit to read termination code of tmp_second_part'
            open(unit=iounit, file='termination_codes/termination_code_' // trim(name_id), action='write', iostat=ierr)
            if (ierr /= 0) stop 'could not open termination_code for ' // trim(name_id)
            write(iounit,'(a)') 'unknown -- cannot open tmp_second_part'
            close(iounit)
            call free_iounit(iounit)
         else
            read(iounit,'(a)') second_part_end_code
            close(iounit)
            call free_iounit(iounit)
            call execute_command_line('rm -f tmp_second_part')
            iounit = alloc_iounit(ierr)
            if (ierr /= 0) stop 'could not alloc_iounit to copy termination code to new file'
            open(unit=iounit, file='termination_codes/termination_code_' // trim(name_id), action='write', iostat=ierr)
            if (ierr /= 0) stop 'could not open termination_codes/termination_code_' // trim(name_id)
            write(iounit,'(a)') trim(second_part_end_code)
            close(iounit)
            call free_iounit(iounit)
         end if

      end subroutine store_termination_code

      end module bin2dco_io
