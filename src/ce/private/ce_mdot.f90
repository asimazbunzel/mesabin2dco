
      module ce_mdot

      use ce_def
      use const_def
      use star_def
      use star_lib
      use binary_def
      use binary_lib
      use math_lib
      
      implicit none

      contains

      
      ! evaluate ce mdot
      ! NOTE: rlo_mdot is the variable that must be set here, and it must be in g/s
      real(dp) function eval_ce_mdot(b, ierr) result(mdot)
         type(binary_info), pointer :: b
         integer, intent(out) :: ierr
         real(dp) :: rl_rel_limit
         real(dp) :: x0, x1, y0, y1, m

         ierr = 0

         if (ce_duration <= years_to_max_mdot_rlof .and. &
             abs(b% mtransfer_rate) * secyer/Msun < max_mdot_rlof) then
            x0 = 0d0; x1 = years_to_max_mdot_rlof
            y0 = lg_mtransfer_rate_start_ce; y1 = safe_log10(max_mdot_rlof)
            m = (y1 - y0) / (x1 - x0)
            mdot = - exp10(m*(ce_duration-x0)+y0) * Msun/secyer
            return
         end if

         if (ce_type == ce_two_stars) then
            rl_rel_limit = tol_two_stars
         else
            rl_rel_limit = tol_xrb
         end if

         if (b% rl_relative_gap(b% d_i) > 0d0) then
            mdot = - max_mdot_rlof * Msun/secyer
         else if (b% rl_relative_gap(b% d_i) < -rl_rel_limit) then
            mdot = - mdot_kh() * Msun/secyer
         else
            x0 = 0d0; x1 = -rl_rel_limit
            y0 = safe_log10(max_mdot_rlof)
            y1 = safe_log10(mdot_kh()) - 2d0
            m = (y1 - y0) / (x1 - x0)
            mdot = - exp10(m*(b% rl_relative_gap(b% d_i)-x0)+y0) * Msun/secyer
            if (ce_dbg) write(*,*) 'reducing mdot_rlof. mdot, target', -mdot * secyer/Msun, exp10(y1)
         end if

         contains

         real(dp) function mdot_kh()
            type(star_info), pointer :: s
            real(dp) :: t_kh

            s => b% s_donor

            t_kh =  0.75d0 * s% cgrav(1) * s% mstar * s% mstar / (s% r(1) * s% L(1))
            mdot_kh = (s% mstar / Msun) / (t_kh / secyer)

         end function mdot_kh
            
      end function eval_ce_mdot

      end module ce_mdot

