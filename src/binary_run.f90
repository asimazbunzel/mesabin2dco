      program binary_run
      use run_star_extras
      use run_binary_extras

      integer :: ierr

      call do_run(extras_controls, extras_binary_controls, ierr)
      if (ierr /= 0) stop

      end program
