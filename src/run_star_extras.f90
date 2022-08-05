! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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

      use const_def
      use chem_def
      use star_def
      use binary_def

      use math_lib
      use star_lib
      use utils_lib, only: is_bad

      implicit none

      integer :: time0, time1, clock_rate

      real(dp) :: delta_HR_limit
      real(dp) :: delta_HR_hard_limit

      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% other_wind => brott_wind

         delta_HR_limit = s% delta_HR_limit
         delta_HR_hard_limit = s% delta_HR_hard_limit

         ! other eps grav is used to restore the homologous vs non
         ! homologous eps grav calculation that was dropped in 15140
         s% use_other_eps_grav = .true.
         s% other_eps_grav => my_eps_grav

      end subroutine extras_controls

      subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         integer :: h1, he4
         real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
            vink_wind, nieu_wind, hamann_wind
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
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
         Z_div_Z_solar = s% kap_rq% Zbase / 0.017d0
         ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
         Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z_div_Z_solar)))

         vink_wind = 0d0
         nieu_wind = 0d0
         hamann_wind = 0d0
         w = 0

         if (b% ignore_rlof_flag) then
            return
         end if

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
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w1 = exp10(logMdot)
            else
               w1 = 0
            end if

            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10(L1/Lsun/1d5) &
                  - 1.339d0*log10(M1/Msun/30) &
                  - 1.601d0*log10(vinf_div_vesc/2d0) &
                  + 1.07d0*log10(T1/2d4) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w2 = exp10(logMdot)
            else
               w2 = 0
            end if

            w = alfa*w1 + (1 - alfa)*w2

         end subroutine eval_Vink_wind

         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 &
                     +1.24d0*log10(L1/Lsun) &
                     +0.16d0*log10(M1/Msun) &
                     +0.81d0*log10(R1/Rsun) &
                     +0.85d0*log10(Z_div_Z_solar)
            w = exp10(log10w)
         end subroutine eval_Nieuwenhuijzen_wind

         subroutine eval_Hamann_wind(w)
            ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -11.95d0 &
                     +1.5d0*log10(L1/Lsun) &
                     -2.85d0*Xs &
                     + 0.85d0*log10(Z_div_Z_solar)
            w = exp10(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (.not. restart) then
            call system_clock(time0, clock_rate)
         end if

      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_start_step = 0

      end function extras_start_step


      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         include 'formats'

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_check_model = keep_going

      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id

         how_many_extra_history_columns = 13

      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: dt
         type(star_info), pointer :: s
         type(binary_info), pointer :: b
         real(dp) :: m_conv, Ebind
         integer :: k

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return

         names(1) = 'Mconv_env'
         names(2) = 'fconv_env'

         if (s% he_core_mass > 0d0 .and. s% center_h1 < 1d-6) then
            m_conv = 0d0
            if (s% he_core_k == 1 .or. (s% star_mass - s% he_core_mass) < 1d-3) then
               vals(1) = 0d0
               vals(2) = 0d0
            else
               do k = 1, s% he_core_k
                  if (s% mixing_type(k) == convective_mixing) m_conv = m_conv + s% dm(k)
               end do
               vals(1) = m_conv / Msun
               vals(2) = (m_conv / Msun) / (s% star_mass - s% he_core_mass)
            end if
         else
            vals(1) = 0d0
            vals(2) = 0d0
         end if

         names(3) = 'Ebind_abs'

         Ebind = 0d0
         do k=1, s% nz
            Ebind = Ebind + s% dm(k) * (-standard_cgrav * s% m(k) / s% r(k) + s% energy(k))
         end do
         vals(3) = abs(Ebind)

         names(4) = 'spin_magnitude'
         vals(4) = clight * s% total_angular_momentum / (standard_cgrav * s% m(1)**2)

         names(5) = 'log_mtransfer_rate'
         vals(5) = safe_log10(abs(b% step_mtransfer_rate)/Msun*secyer)

         names(6) = 'eff_xfer_fraction'
         if (b% component_mdot(b% d_i) == 0d0) then
            vals(6) = 1d0
         else
            vals(6) = (-b% component_mdot(b% a_i))/(b% component_mdot(b% d_i))
         end if

         names(7) = 'orbital_period'
         vals(7) = b% period/(60d0*60d0*24d0)

         names(8) = 'separation'
         vals(8) = b% separation/Rsun

         names(9) = 'Prot'
         vals(9) = 2 * pi / s% omega_avg_surf /(60d0*60d0*24d0)

         names(10) = 'Prot_div_Porb'
         vals(10) = 2 * pi / s% omega_avg_surf / b% period

         names(11) = 'log_accretion_luminosity'
         vals(11) = safe_log10(b% accretion_luminosity / Lsun)

         names(12) = 'accretion_mode'
         vals(12) = b% accretion_mode

         dt = dble(time1 - time0) / clock_rate / 60
         names(13) = 'runtime_minutes'
         vals(13) = dt

      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id

         how_many_extra_profile_columns = 1

      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         type(star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'Ebind'
         vals(:,1) = 0d0
         vals(1,1) = s% dm(1) * (-standard_cgrav * s% m(1) / s% r(1) + s% energy(1))
         do k=2, s% nz
            vals(k,1) = vals(k-1,1) + s% dm(k) * (-standard_cgrav * s% m(k) / s% r(k) + s% energy(k))
         end do

      end subroutine data_for_extra_profile_columns


      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call system_clock(time1, clock_rate)

         extras_finish_step = keep_going

      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

      end subroutine extras_after_evolve


      subroutine my_eps_grav(id, k, dt, ierr)
         integer, intent(in) :: id, k
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call eval1_eps_grav_and_partials(s, k, ierr)
         if (ierr /= 0) return

      end subroutine my_eps_grav

      subroutine eval1_eps_grav_and_partials(s, k, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr

         include 'formats'
         ierr = 0

         call do_eps_grav_with_lnd(s, k, ierr)

         if (ierr /= 0 .or. is_bad(s% eps_grav(k))) then
            ierr = -1
            if (s% report_ierr) then
               write(*,2) &
                  'failed in eval_eps_grav_and_partials', k, s% eps_grav(k)
            end if
            if (s% stop_for_bad_nums) then
               stop 'eval1_eps_grav_and_partials'
            end if
            return
         end if

      end subroutine eval1_eps_grav_and_partials

      subroutine do_eps_grav_with_lnd(s, k, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: alfa
         include 'formats'
         alfa = get_Eulerian_fraction_for_eps_grav(s, k, ierr)
         if (ierr /= 0) return
         call combine_two_eps_gravs( &
            s, k, alfa, 1d0 - alfa, &
            do_eps_grav_with_lnd_Eulerian, &
            do_eps_grav_with_lnd_Lagrangian, ierr)
         if (ierr /= 0) return
         !call include_composition_in_eps_grav(s, k)
      end subroutine do_eps_grav_with_lnd

      subroutine do_eps_grav_with_lnd_Eulerian(s, k, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         include 'formats'
         call combine_two_eps_gravs( &
            s, k, 1d0, 1d0, do_spatial_term, do_dt_at_const_q_with_lnd, ierr)
      end subroutine do_eps_grav_with_lnd_Eulerian

      subroutine do_eps_grav_with_lnd_Lagrangian(s, k, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         call do_lnd_eps_grav(s, k, s% dlnd_dt(k), s% dlnT_dt(k), ierr)
      end subroutine do_eps_grav_with_lnd_Lagrangian

      subroutine do_dt_at_const_q_with_lnd(s, k, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         call do_lnd_eps_grav( &
            s, k, s% dlnd_dt_const_q(k), s% dlnT_dt_const_q(k), ierr)
      end subroutine do_dt_at_const_q_with_lnd

      subroutine do_spatial_term(s, k, ierr)
         use eos_def, only: i_Cp, i_grad_ad
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr

         ! Townsley & Bildsten, The Astrophysical Journal, 600:390–403, 2004 January 1
         ! ds/dm = (ds/dP)*(dP/dm)
         ! dP/dm = -Gm/(4 pi r^4)
         ! ds/dP = (gradT - grada)*Cp/P
         ! gradT = (dlnT/dlnP)_actual
         ! grada = (dlnT/dlnP)_adiabatic
         ! in mesa/star m = M_center + q*xmstar,
         ! so (dm/dt)_q = q*Mdot since M_center is constant during timesteps
         ! and d(xmstar)/dt = Mdot.
         ! Putting it all together,
         ! we get eps_grav_h = (-G m q Mdot Cp T)*(gradT - grada)/(4 pi r^4 P)

         real(dp) :: &
            P2, rmid3, rmid4, inv_rmid4, inv_4pi_r4, P, inv_4pi_r4_P, &
            gradT_mid, dg, mmid, qmid, G_m_q_mdot, CpT, &
            d_inv_rmid4_dlnR00, d_inv_rmid4_dlnRp1, &
            d_inv_4pi_r4_P_dlnR00, d_inv_4pi_r4_P_dlnRp1, &
            d_dg_dlnR00, d_dg_dlnRp1, d_dg_dL00, d_dg_dLp1, &
            dgradT_mid_dlndm1, dgradT_mid_dlnd00, dgradT_mid_dlndp1, &
            d_dg_dlndm1, d_dg_dlnd00, d_dg_dlndp1, &
            dgradT_mid_dlnTm1, dgradT_mid_dlnT00, dgradT_mid_dlnTp1, &
            d_dg_dlnTm1, d_dg_dlnT00, d_dg_dlnTp1, d_CpT_dlnd, d_CpT_dlnT, &
            d_inv_4pi_r4_P_dlnd, d_inv_4pi_r4_P_dlnT, &
            dP_dlnPgas_const_T, dP_dlnT_const_Pgas, &
            dPinv_dlnPgas_const_T, dPinv_dlnT_const_Pgas, &
            d_inv_4pi_r4_P_dlnPgas_const_T, d_inv_4pi_r4_P_dlnT_const_Pgas, &
            d_CpT_dlnPgas_const_T, d_CpT_dlnT_const_Pgas, &
            d_dg_dlnPgas00_const_T, d_dg_dlnT00_const_Pgas, &
            d_dg_dlnPgasp1_const_T, d_dg_dlnTp1_const_Pgas, &
            d_dg_dlnPgasm1_const_T, d_dg_dlnTm1_const_Pgas

         include 'formats'
         ierr = 0

         call zero_eps_grav_and_partials(s, k)
         if (k == s% nz) return

         rmid3 = s% rmid(k)*s% rmid(k)*s% rmid(k)
         rmid4 = rmid3*s% rmid(k)
         inv_rmid4 = 1d0/rmid4
         inv_4pi_r4 = inv_rmid4/(4*pi)
         P = s% P(k)
         inv_4pi_r4_P = inv_4pi_r4/P
         gradT_mid = 0.5d0*(s% gradT(k) + s% gradT(k+1))
         dg = s% grada(k) - gradT_mid
         mmid = 0.5d0*(s% m(k) + s% m(k+1))
         qmid = 0.5d0*(s% q(k) + s% q(k+1))
         G_m_q_mdot = s% cgrav(k)*mmid*qmid*s% mstar_dot
         CpT = s% Cp(k)*s% T(k)

         s% eps_grav(k) = G_m_q_mdot*CpT*dg*inv_4pi_r4_P

         if (is_bad(s% eps_grav(k))) then
            ierr = -1
            if (s% report_ierr) &
               write(*,2) 'do_spatial_term -- bad value for eps_grav', k, s% eps_grav(k)
            if (s% stop_for_bad_nums) stop 'do_spatial_term'
            return
         end if

         d_inv_rmid4_dlnR00 = -2*inv_rmid4*s% r(k)*s% r(k)*s% r(k)/rmid3
         d_inv_rmid4_dlnRp1 = -2*inv_rmid4*s% r(k+1)*s% r(k+1)*s% r(k+1)/rmid3
         d_inv_4pi_r4_P_dlnR00 = d_inv_rmid4_dlnR00/(4*pi*P)
         d_inv_4pi_r4_P_dlnRp1 = d_inv_rmid4_dlnRp1/(4*pi*P)

         d_dg_dlnR00 = -0.5d0*s% d_gradT_dlnR(k)
         d_dg_dlnRp1 = -0.5d0*s% d_gradT_dlnR(k+1)

         s% d_eps_grav_dlnR00(k) = G_m_q_mdot*CpT* &
            (d_dg_dlnR00*inv_4pi_r4_P + dg*d_inv_4pi_r4_P_dlnR00)
         s% d_eps_grav_dlnRp1(k) = G_m_q_mdot*CpT* &
            (d_dg_dlnRp1*inv_4pi_r4_P + dg*d_inv_4pi_r4_P_dlnRp1)

         d_dg_dL00 = -0.5d0*s% d_gradT_dL(k)
         d_dg_dLp1 = -0.5d0*s% d_gradT_dL(k+1)

         s% d_eps_grav_dL00(k) = d_dg_dL00*G_m_q_mdot*CpT*inv_4pi_r4_P
         s% d_eps_grav_dLp1(k) = d_dg_dLp1*G_m_q_mdot*CpT*inv_4pi_r4_P

         dgradT_mid_dlndm1 = 0.5d0*s% d_gradT_dlndm1(k)
         dgradT_mid_dlnd00 = 0.5d0*(s% d_gradT_dlnd00(k) + s% d_gradT_dlndm1(k+1))
         dgradT_mid_dlndp1 = 0.5d0*s% d_gradT_dlnd00(k+1)

         d_dg_dlndm1 = -dgradT_mid_dlndm1
         d_dg_dlnd00 = s% d_eos_dlnd(i_grad_ad,k) - dgradT_mid_dlnd00
         d_dg_dlndp1 = -dgradT_mid_dlndp1

         dgradT_mid_dlnTm1 = 0.5d0*s% d_gradT_dlnTm1(k)
         dgradT_mid_dlnT00 = 0.5d0*(s% d_gradT_dlnT00(k) + s% d_gradT_dlnTm1(k+1))
         dgradT_mid_dlnTp1 = 0.5d0*s% d_gradT_dlnT00(k+1)

         d_dg_dlnTm1 = -dgradT_mid_dlnTm1
         d_dg_dlnT00 = s% d_eos_dlnT(i_grad_ad,k) - dgradT_mid_dlnT00
         d_dg_dlnTp1 = -dgradT_mid_dlnTp1

         d_CpT_dlnd = s% d_eos_dlnd(i_Cp,k)*s% T(k)
         d_CpT_dlnT = s% d_eos_dlnT(i_Cp,k)*s% T(k) + CpT

         d_inv_4pi_r4_P_dlnd = -s% chiRho_for_partials(k)*inv_4pi_r4_P
         d_inv_4pi_r4_P_dlnT = -s% chiT_for_partials(k)*inv_4pi_r4_P

         s% d_eps_grav_dlndm1(k) = d_dg_dlndm1*G_m_q_mdot*CpT*inv_4pi_r4_P
         s% d_eps_grav_dlndp1(k) = d_dg_dlndp1*G_m_q_mdot*CpT*inv_4pi_r4_P
         s% d_eps_grav_dlnd00(k) = G_m_q_mdot*( &
            d_dg_dlnd00*CpT*inv_4pi_r4_P + &
            dg*d_CpT_dlnd*inv_4pi_r4_P + &
            dg*CpT*d_inv_4pi_r4_P_dlnd)

         s% d_eps_grav_dlnTm1(k) = d_dg_dlnTm1*G_m_q_mdot*CpT*inv_4pi_r4_P
         s% d_eps_grav_dlnTp1(k) = d_dg_dlnTp1*G_m_q_mdot*CpT*inv_4pi_r4_P
         s% d_eps_grav_dlnT00(k) = G_m_q_mdot*( &
            d_dg_dlnT00*CpT*inv_4pi_r4_P + &
            dg*d_CpT_dlnT*inv_4pi_r4_P + &
            dg*CpT*d_inv_4pi_r4_P_dlnT)

      end subroutine do_spatial_term

      ! this uses the given args to calculate -T*ds/dt
      subroutine do_lnd_eps_grav(s, k, dlnd_dt, dlnT_dt, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         use chem_def, only: ihe4
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dlnd_dt, dlnT_dt
         integer, intent(out) :: ierr

         real(dp) :: dT1_dlnTdot, a1, da1_dlnd, da1_dlnT, &
            T1, dT1_dlnd, dT1_dlnT00, dT1_dlnd_dt, dT1_d_dlnTdt, &
            a2, da2_dlnd, da2_dlnT, &
            T2, dT2_dlnT, dT2_dlnd00, &
            T3, dT3_dlnd, dT3_dlnT, center_he4
         integer :: he4
         logical :: test_partials

         include 'formats'
         ierr = 0

         !test_partials = (k == s% hydro_test_partials_k)
         test_partials = .false.

         call zero_eps_grav_and_partials(s, k)

         he4 = s% net_iso(ihe4)
         if (he4 > 0) then
            center_he4 = s% xa(he4,s% nz)
         else
            center_he4 = 0d0
         end if
         !if ((s% power_he_burn > s% alt_eps_grav_he_burn_limit .and. center_he4 > 0.8d0) .or. &
         !      s% T(s% nz) > s% alt_eps_grav_T_center_limit) then
         !   s% eps_grav(k) = -s% T(k)*s% cp_start(k)* &
         !         ((1-s% grada_start(k)*s% chiT_start(k))*dlnT_dt &
         !           - s% grada_start(k)*s% chiRho_start(k)*dlnd_dt)
         !   s% d_eps_grav_dlnd00(k) = &
         !      s% T(k)*s% cp_start(k)*s% grada_start(k)*s% chiRho_start(k)*s% dVARDOT_dVAR
         !   s% d_eps_grav_dlnT00(k) = s% eps_grav(k) - &
         !      s% T(k)*s% cp_start(k)*(1-s% grada_start(k)*s% chiT_start(k))*s% dVARDOT_dVAR
         !   return
         !end if

         a1 = 1 - s% grada(k)*s% chiT(k)
         da1_dlnd = -(s% d_eos_dlnd(i_grad_ad,k)*s% chiT(k) + s% grada(k)*s% d_eos_dlnd(i_chiT,k))
         da1_dlnT = -(s% d_eos_dlnT(i_grad_ad,k)*s% chiT(k) + s% grada(k)*s% d_eos_dlnT(i_chiT,k))

         T1 = dlnT_dt*a1
         dT1_dlnd = dlnT_dt*da1_dlnd

         dT1_d_dlnTdt = a1
         dT1_dlnT00 = s% dVARDOT_dVAR*a1 + dlnT_dt*da1_dlnT

         a2 = s% grada(k)*s% chiRho(k)
         da2_dlnd = s% d_eos_dlnd(i_grad_ad,k)*s% chiRho(k) + s% grada(k)*s% d_eos_dlnd(i_chiRho,k)
         da2_dlnT = s% d_eos_dlnT(i_grad_ad,k)*s% chiRho(k) + s% grada(k)*s% d_eos_dlnT(i_chiRho,k)

         T2 = dlnd_dt*a2
         dT2_dlnT = dlnd_dt*da2_dlnT

         dT2_dlnd00 = s% dVARDOT_dVAR*a2 + dlnd_dt*da2_dlnd

         T3 = -s% T(k)*s% cp(k)
         dT3_dlnd = -s% T(k)*s% d_eos_dlnd(i_Cp,k)
         dT3_dlnT = -s% T(k)*(s% cp(k) + s% d_eos_dlnT(i_Cp,k))

         ! eps_grav = T3*(T1-T2)
         s% eps_grav(k) = T3*(T1-T2)

         if (is_bad(s% eps_grav(k))) then
            ierr = -1
            if (s% report_ierr) &
               write(*,2) 'do_lnd_eps_grav -- bad value for eps_grav', k, s% eps_grav(k)
            if (s% stop_for_bad_nums) stop 'do_lnd_eps_grav'
            return
         end if

         s% d_eps_grav_dlndm1(k) = 0
         s% d_eps_grav_dlndp1(k) = 0
         s% d_eps_grav_dlnd00(k) = (T3*(dT1_dlnd - dT2_dlnd00) + dT3_dlnd*(T1-T2))

         s% d_eps_grav_dlnTm1(k) = 0
         s% d_eps_grav_dlnTp1(k) = 0
         s% d_eps_grav_dlnT00(k) = (T3*(dT1_dlnT00 - dT2_dlnT) + dT3_dlnT*(T1-T2))

      end subroutine do_lnd_eps_grav

      recursive subroutine combine_two_eps_gravs( &
            s, k, alfa, beta, eps_grav_proc1, eps_grav_proc2, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: alfa, beta
         interface
            subroutine eps_grav_proc1(s, k, ierr)
               use star_def, only: star_info
               use const_def, only: dp
               type (star_info), pointer :: s
               integer, intent(in) :: k
                     integer, intent(out) :: ierr
            end subroutine eps_grav_proc1
            subroutine eps_grav_proc2(s, k, ierr)
               use star_def, only: star_info
               use const_def, only: dp
               type (star_info), pointer :: s
               integer, intent(in) :: k
                     integer, intent(out) :: ierr
            end subroutine eps_grav_proc2
         end interface
         integer, intent(out) :: ierr

         real(dp) :: &
            eps_grav, d_eps_grav_dlndm1, d_eps_grav_dlnd00, d_eps_grav_dlndp1, &
            d_eps_grav_dlnTm1, d_eps_grav_dlnT00, d_eps_grav_dlnTp1, &
            d_eps_grav_dlnR00, d_eps_grav_dlnRp1, d_eps_grav_dL00, d_eps_grav_dLp1, &
            d_eps_grav_dlnPgas00_const_T, &
            d_eps_grav_dlnPgasm1_const_T, d_eps_grav_dlnPgasp1_const_T, &
            d_eps_grav_dlnTm1_const_Pgas, d_eps_grav_dlnT00_const_Pgas, &
            d_eps_grav_dlnTp1_const_Pgas, d_eps_grav_dv00, d_eps_grav_dvp1

         include 'formats'
         ierr = 0

         ! alfa is multiplier of result from calling eps_grav_proc1
         ! beta is multiplier of result from calling eps_grav_proc2
         ! i.e., eps_grav = alfa*eps_grav1 + beta*eps_grav2

         if (alfa > 1d0 .or. alfa < 0d0 .or. beta > 1d0 .or. beta < 0d0) then
            if (s% report_ierr) &
               write(*,2) 'combine_two_eps_gravs: alfa beta', k, alfa, beta
            ierr = -1
            return
         end if

         ! result is alfa*eps_grav_proc1 + beta*eps_grav_proc2

         if (alfa > 0d0) then
            call eps_grav_proc1(s, k, ierr)
            if (ierr /= 0) return
            if (beta == 0d0) return
            ! save results
            eps_grav = s% eps_grav(k)
            d_eps_grav_dlndm1 = s% d_eps_grav_dlndm1(k)
            d_eps_grav_dlnd00 = s% d_eps_grav_dlnd00(k)
            d_eps_grav_dlndp1 = s% d_eps_grav_dlndp1(k)
            d_eps_grav_dlnTm1 = s% d_eps_grav_dlnTm1(k)
            d_eps_grav_dlnT00 = s% d_eps_grav_dlnT00(k)
            d_eps_grav_dlnTp1 = s% d_eps_grav_dlnTp1(k)
            d_eps_grav_dlnR00 = s% d_eps_grav_dlnR00(k)
            d_eps_grav_dlnRp1 = s% d_eps_grav_dlnRp1(k)
            d_eps_grav_dL00 = s% d_eps_grav_dL00(k)
            d_eps_grav_dLp1 = s% d_eps_grav_dLp1(k)
            d_eps_grav_dlnPgas00_const_T = s% d_eps_grav_dlnPgas00_const_T(k)
            d_eps_grav_dlnPgasm1_const_T = s% d_eps_grav_dlnPgasm1_const_T(k)
            d_eps_grav_dlnPgasp1_const_T = s% d_eps_grav_dlnPgasp1_const_T(k)
            d_eps_grav_dlnTm1_const_Pgas = s% d_eps_grav_dlnTm1_const_Pgas(k)
            d_eps_grav_dlnT00_const_Pgas = s% d_eps_grav_dlnT00_const_Pgas(k)
            d_eps_grav_dlnTp1_const_Pgas = s% d_eps_grav_dlnTp1_const_Pgas(k)
            d_eps_grav_dv00 = s% d_eps_grav_dv00(k)
            d_eps_grav_dvp1 = s% d_eps_grav_dvp1(k)
         else ! not needed, but to keep the compiler happy we set these to 0
            eps_grav = 0
            d_eps_grav_dlndm1 = 0
            d_eps_grav_dlnd00 = 0
            d_eps_grav_dlndp1 = 0
            d_eps_grav_dlnTm1 = 0
            d_eps_grav_dlnT00 = 0
            d_eps_grav_dlnTp1 = 0
            d_eps_grav_dlnR00 = 0
            d_eps_grav_dlnRp1 = 0
            d_eps_grav_dL00 = 0
            d_eps_grav_dLp1 = 0
            d_eps_grav_dlnPgas00_const_T = 0
            d_eps_grav_dlnPgasm1_const_T = 0
            d_eps_grav_dlnPgasp1_const_T = 0
            d_eps_grav_dlnTm1_const_Pgas = 0
            d_eps_grav_dlnT00_const_Pgas = 0
            d_eps_grav_dlnTp1_const_Pgas = 0
            d_eps_grav_dv00 = 0
            d_eps_grav_dvp1 = 0
         end if

         call eps_grav_proc2(s, k, ierr)
         if (ierr /= 0) return
         if (alfa == 0d0) return

         ! combine results
         s% eps_grav(k) = alfa*eps_grav + beta*s% eps_grav(k)

         s% d_eps_grav_dlndm1(k) = alfa*d_eps_grav_dlndm1 + beta*s% d_eps_grav_dlndm1(k)
         s% d_eps_grav_dlnd00(k) = alfa*d_eps_grav_dlnd00 + beta*s% d_eps_grav_dlnd00(k)
         s% d_eps_grav_dlndp1(k) = alfa*d_eps_grav_dlndp1 + beta*s% d_eps_grav_dlndp1(k)

         s% d_eps_grav_dlnTm1(k) = alfa*d_eps_grav_dlnTm1 + beta*s% d_eps_grav_dlnTm1(k)
         s% d_eps_grav_dlnT00(k) = alfa*d_eps_grav_dlnT00 + beta*s% d_eps_grav_dlnT00(k)
         s% d_eps_grav_dlnTp1(k) = alfa*d_eps_grav_dlnTp1 + beta*s% d_eps_grav_dlnTp1(k)

         s% d_eps_grav_dlnPgas00_const_T(k) = &
            alfa*d_eps_grav_dlnPgas00_const_T + beta*s% d_eps_grav_dlnPgas00_const_T(k)
         s% d_eps_grav_dlnPgasm1_const_T(k) = &
            alfa*d_eps_grav_dlnPgasm1_const_T + beta*s% d_eps_grav_dlnPgasm1_const_T(k)
         s% d_eps_grav_dlnPgasp1_const_T(k) = &
            alfa*d_eps_grav_dlnPgasp1_const_T + beta*s% d_eps_grav_dlnPgasp1_const_T(k)

         s% d_eps_grav_dlnTm1_const_Pgas(k) = &
            alfa*d_eps_grav_dlnTm1_const_Pgas + beta*s% d_eps_grav_dlnTm1_const_Pgas(k)
         s% d_eps_grav_dlnT00_const_Pgas(k) = &
            alfa*d_eps_grav_dlnT00_const_Pgas + beta*s% d_eps_grav_dlnT00_const_Pgas(k)
         s% d_eps_grav_dlnTp1_const_Pgas(k) = &
            alfa*d_eps_grav_dlnTp1_const_Pgas + beta*s% d_eps_grav_dlnTp1_const_Pgas(k)

         s% d_eps_grav_dlnR00(k) = alfa*d_eps_grav_dlnR00 + beta*s% d_eps_grav_dlnR00(k)
         s% d_eps_grav_dlnRp1(k) = alfa*d_eps_grav_dlnRp1 + beta*s% d_eps_grav_dlnRp1(k)

         s% d_eps_grav_dL00(k) = alfa*d_eps_grav_dL00 + beta*s% d_eps_grav_dL00(k)
         s% d_eps_grav_dLp1(k) = alfa*d_eps_grav_dLp1 + beta*s% d_eps_grav_dLp1(k)

         s% d_eps_grav_dv00(k) = alfa*d_eps_grav_dv00 + beta*s% d_eps_grav_dv00(k)
         s% d_eps_grav_dvp1(k) = alfa*d_eps_grav_dvp1 + beta*s% d_eps_grav_dvp1(k)

      end subroutine combine_two_eps_gravs

      real(dp) function get_Eulerian_fraction_for_eps_grav(s, k, ierr) result(alfa)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: q_Eulerian, q_Lagrangian
         include 'formats'
         ! alfa is fraction of Eulerian form in result

         ierr = 0

         if (k < s% k_below_const_q .or. &
               s% k_below_const_q > s% nz) then
            alfa = 1d0 ! pure Eulerian
         else if (k >= s% k_const_mass) then
            alfa = 0d0 ! pure Lagrangian
         else
            if (s% k_below_const_q == 1) then
               q_Eulerian = 1d0
            else
               q_Eulerian = s% q(s% k_below_const_q-1)
            end if
            if (s% k_const_mass > s% nz) then
               q_Lagrangian = 0d0
            else
               q_Lagrangian = s% q(s% k_const_mass)
            end if
            alfa = max(0d0, min(1d0, (s% q(k) - q_Lagrangian)/(q_Eulerian - q_Lagrangian)))
            if (is_bad(alfa)) then
               write(*,2) 's% k_below_Eulerian_eps_grav', s% k_below_const_q, q_Eulerian
               write(*,2) 's% q(k)', k, s% q(k)
               write(*,2) 's% k_Lagrangian_eps_grav', s% k_const_mass, q_Lagrangian
               write(*,2) 'alfa', k, alfa
               write(*,*) 'failed in get_Eulerian_fraction_for_eps_grav'
            end if
         end if

      end function get_Eulerian_fraction_for_eps_grav

      subroutine zero_eps_grav_and_partials(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% eps_grav(k) = 0
         s% d_eps_grav_dlndm1(k) = 0
         s% d_eps_grav_dlnd00(k) = 0
         s% d_eps_grav_dlndp1(k) = 0
         s% d_eps_grav_dlnTm1(k) = 0
         s% d_eps_grav_dlnT00(k) = 0
         s% d_eps_grav_dlnTp1(k) = 0
         s% d_eps_grav_dlnPgasm1_const_T(k) = 0
         s% d_eps_grav_dlnPgas00_const_T(k) = 0
         s% d_eps_grav_dlnPgasp1_const_T(k) = 0
         s% d_eps_grav_dlnTm1_const_Pgas(k) = 0
         s% d_eps_grav_dlnT00_const_Pgas(k) = 0
         s% d_eps_grav_dlnTp1_const_Pgas(k) = 0
         s% d_eps_grav_dlnR00(k) = 0
         s% d_eps_grav_dlnRp1(k) = 0
         s% d_eps_grav_dL00(k) = 0
         s% d_eps_grav_dLp1(k) = 0
         s% d_eps_grav_dv00(k) = 0
         s% d_eps_grav_dvp1(k) = 0
         s% d_eps_grav_dx(:,k) = 0
         s% eps_grav_composition_term(k) = 0
      end subroutine zero_eps_grav_and_partials

      end module run_star_extras
