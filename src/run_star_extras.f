! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib
      use chem_def
      use binary_def
      
      implicit none

      ! find zams info
      logical :: found_zams

      real(dp) :: X_central_initial
      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% other_wind => brott_wind_prescription
         
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras =.false.       
            
      end subroutine extras_controls
      
      
      ! use Brott+ (2011) prescription for winds
      subroutine brott_wind_prescription(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr
         integer :: h1, he4
         real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa
         real(dp) :: L1, M1, R1, T1
         real(dp) :: vink_wind, nieu_wind, hamann_wind
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         L1 = Lsurf
         M1 = Msurf
         R1 = Rsurf
         T1 = Tsurf

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         Xs = s% xa(h1,1)
         Ys = s% xa(he4,1)
         ! Z=0.017 is Z from Grevesse et al. 1996
         Z_div_Z_solar = s% Zbase/0.017d0
         ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
         Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10_cr(Z_div_Z_solar)))

         vink_wind = 0d0
         nieu_wind = 0d0
         hamann_wind = 0d0
         w = 0d0

         call eval_Vink_wind(vink_wind)
         call eval_Nieuwenhuijzen_wind(nieu_wind)
         call eval_Hamann_wind(hamann_wind)

         ! use 1/10 hamann
         hamann_wind = hamann_wind/10d0

         if (T1 < Teff_jump) then
            ! low T wind
            w = max(vink_wind, nieu_wind)
         else
            ! high T wind
            alfa = 0d0
            if (Xs > 0.7d0) then
               alfa = 1d0
            else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
               alfa = (Xs - 0.4d0)/0.3d0
            end if
            w = alfa * vink_wind + (1d0-alfa) * hamann_wind
         end if

         ierr = 0

         contains

            subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = (T1 - (Teff_jump - dT)) / (2*dT)
               end if
            end if

            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow_cr(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10_cr(L1/Lsun/1d5) &
                  - 1.313d0*log10_cr(M1/Msun/30) &
                  - 1.226d0*log10_cr(vinf_div_vesc/2d0) &
                  + 0.933d0*log10_cr(T1/4d4) &
                  - 10.92d0*pow2(log10_cr(T1/4d4)) &
                  + 0.85d0*log10_cr(Z_div_Z_solar)
               w1 = exp10_cr(logMdot)
            else
               w1 = 0
            end if

            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow_cr(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10_cr(L1/Lsun/1d5) &
                  - 1.339d0*log10_cr(M1/Msun/30) &
                  - 1.601d0*log10_cr(vinf_div_vesc/2d0) &
                  + 1.07d0*log10_cr(T1/2d4) &
                  + 0.85d0*log10_cr(Z_div_Z_solar)
               w2 = exp10_cr(logMdot)
            else
               w2 = 0
            end if
            w = alfa*w1 + (1 - alfa)*w2
         end subroutine eval_Vink_wind

         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            log10w = -14.02d0 &
                     +1.24d0*log10_cr(L1/Lsun) &
                     +0.16d0*log10_cr(M1/Msun) &
                     +0.81d0*log10_cr(R1/Rsun) &
                     +0.85d0*log10_cr(Z_div_Z_solar)
            w = exp10_cr(log10w)
         end subroutine eval_Nieuwenhuijzen_wind

         subroutine eval_Hamann_wind(w)
            ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
            real(dp), intent(out) :: w
            real(dp) :: log10w
            log10w = -11.95d0 &
                     +1.5d0*log10_cr(L1/Lsun) &
                     -2.85d0*Xs &
                     + 0.85d0*log10_cr(Z_div_Z_solar)
            w = exp10_cr(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind_prescription


      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         extras_startup = 0
         
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if

         found_zams = .false.

         X_central_initial = s% center_h1

      end function extras_startup
      

      integer function extras_start_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: Lnuc_div_L
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         extras_start_step = 0

         ! Aaron Dotter's ZAMS definition in MIST is:
         !
         ! The ZAMS EEP is taken as the first point after the
         ! H-burning luminosity exceeds 99.9% of the total
         ! luminosity and before the central H mass fraction has
         ! fallen below its initial value by 0.0015.
         Lnuc_div_L = s% L_nuc_burn_total / s% L_phot
         if (Lnuc_div_L > 0.999d0 .and. .not. found_zams) then
            if (X_central_initial - s% center_h1 > 0.015d0) then
               write(*,'(a)') 'found ZAMS!'
               found_zams = .true.
               s% xtra30 = s% L_phot
            end if
         end if


      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         integer, parameter :: phase_carbon_burning = 7
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 3
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         integer :: i, j, h1, core_location, k
         real(dp) :: mdot_edd, E_grav, E_int
         real(dp) :: m_conv

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call binary_ptr(s% binary_id, b, ierr)
         
         names(1) = 'mdot_edd'
         mdot_edd = (pi4 * clight * s% r(1)) / s% photosphere_opacity  ! in g/s
         mdot_edd = (mdot_edd * secyer) / Msun  ! now, Msun/yr
         vals(1) = mdot_edd

         names(2) = 'conv_envelope_mass'
         names(3) = 'fract_conv_env'
         if (s% he_core_mass > 0d0 .and. s% center_h1 < 1d-6) then
            m_conv = 0d0
            if (s% he_core_k == 1 .or. (s% star_mass - s% he_core_mass) < 1d-3) then
               vals(2) = 0d0
               vals(3) = 0d0
            else
               do k = 1, s% he_core_k
                  if (s% mixing_type(k) == convective_mixing) m_conv = m_conv + s% dm(k)
               end do
               vals(2) = m_conv / Msun
               vals(3) = (m_conv / Msun) / (s% star_mass - s% he_core_mass)
            end if
         else
            vals(2) = 0d0
            vals(3) = 0d0
         end if

      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
      end subroutine data_for_extra_profile_columns


      subroutine how_many_extra_history_header_items(id, id_extra, num_cols)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: num_cols
         num_cols = 0
      end subroutine how_many_extra_history_header_items
     

      subroutine data_for_extra_history_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
         integer, intent(in) :: id, id_extra, num_extra_header_items
         character (len=*), pointer :: extra_header_item_names(:)
         real(dp), pointer :: extra_header_item_vals(:)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         !here is an example for adding an extra history header item
         !set num_cols=1 in how_many_extra_history_header_items and then unccomment these lines
         !extra_header_item_names(1) = 'mixing_length_alpha'
         !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items


      subroutine how_many_extra_profile_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols = 0
      end subroutine how_many_extra_profile_header_items
      
      subroutine data_for_extra_profile_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr

      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra profile header item
      !set num_cols=1 in how_many_extra_profile_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info


      end module run_star_extras
