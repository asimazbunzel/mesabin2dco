
      module cc_ctrls_io

      use cc_def
      use const_def

      implicit none

      namelist /cc_controls/ &
         model_name, &
         max_ns_mass, &
         cc_data_directory, &
         filename_for_star_data, &
         filename_for_binary_data

      contains

      subroutine set_default_controls

         model_name = 'rapid'
         max_ns_mass = 2.5d0
         cc_data_directory = 'cc_data'
         filename_for_star_data = 'star_at_core_collapse'
         filename_for_binary_data = 'binary_at_core_collapse'

      end subroutine set_default_controls

      subroutine read_cc_controls(fname, ierr)
         use utils_lib, only: mkdir
         character (len=*), intent(in) :: fname
         integer, intent(out) :: ierr
         integer :: iounit

         ierr = 0

         open(newunit=iounit, file=trim(fname), action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname)
            return
         end if
         read(iounit, nml=cc_controls, iostat=ierr)
         close(iounit)
         if (ierr /= 0) then
            write(*,'(a)')
            write(*,'(a)') 'failed to read cc_controls'
            write(*,'(a)') 'the following runtime error can help you find the problem'
            write(*,'(a)')
            open(newunit=iounit, file=trim(fname), action='read', delim='quote', status='old', iostat=ierr)
            read(iounit, nml=cc_controls)
            close(iounit)
            return
         end if

         ! after reading cc_data_directory, run `mkdir -p` on it
         call mkdir(cc_data_directory)

      end subroutine read_cc_controls

      end module cc_ctrls_io
