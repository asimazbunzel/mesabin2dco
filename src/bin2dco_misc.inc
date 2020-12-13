
         subroutine read_bin2dco_controls(filename, ierr)
            character(len=*) :: filename
            integer, intent(out) :: ierr
            integer :: unit

            open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'Failed to open control namelist file ', trim(filename)
            else
               read(unit, nml=bin2dco_controls, iostat=ierr)
               close(unit)
               if (ierr /= 0) then
                  write(*,'(a)') 'Failed while trying to read bin2dco controls namelist file ', trim(filename)
                  write(*,'(a)') &
                     'The following runtime error message might help you find the problem'
                  write(*, *)
                  open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
                  read(unit, nml=bin2dco_controls)
                  close(unit)
               end if
            end if

         end subroutine read_bin2dco_controls

         logical function end_as_core_collapse(filename)
            use const_def, only: strlen
            character(len=*) :: filename
            character(len=strlen) :: string
            integer :: iounit, ierr

            end_as_core_collapse = .false.

            open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed to open ' // trim(filename)
               return
            else
               read(iounit, '(a)') string
               close(iounit)
               if (string == 'core-collapse') end_as_core_collapse = .true.
            end if

         end function end_as_core_collapse

         ! search for string in file and return value associated to it
         subroutine read_cc_parameters(star_filename, binary_filename, &
            mass_of_progenitor, mass_of_remnant, mass_of_companion, &
            pre_cc_separation, pre_cc_period, &
            ierr)
            character(len=*) :: star_filename, binary_filename
            real(dp), intent(out) :: mass_of_progenitor, mass_of_remnant, mass_of_companion
            real(dp), intent(out) :: pre_cc_separation, pre_cc_period
            integer, intent(out) :: ierr

            call read_parameter('remnant_mass', star_filename, mass_of_remnant,  ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed to get remnant_mass with read_parameter'
               return
            end if
            call read_parameter('progenitor_mass', binary_filename, mass_of_progenitor, ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed to get progenitor_mass with read_parameter'
               return
            end if
            call read_parameter('companion_mass', binary_filename, mass_of_companion, ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed to get companion_mass with read_parameter'
               return
            end if
            call read_parameter('separation_pre_cc', binary_filename, pre_cc_separation, ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed to get separation_pre_cc with read_parameter'
               return
            end if
            call read_parameter('period_pre_cc', binary_filename, pre_cc_period, ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed to get period_pre_cc with read_parameter'
               return
            end if

         end subroutine read_cc_parameters

         subroutine read_parameter(string, filename, val, ierr)
            use const_def, only: dp, strlen
            character(len=*) :: string
            character(len=*) :: filename
            real(dp), intent(out) :: val
            integer, intent(out) :: ierr
            character(len=strlen) :: buffer, trimmed_buffer, label
            integer :: pos
            logical :: dbg = .false.
            integer :: iounit

            ! allocate file unit
            open(newunit=iounit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(filename)
               return
            end if

            do while (ierr == 0)
               read(iounit, '(A)', iostat=ierr) buffer
               trimmed_buffer = trim(adjustl(buffer))

               ! Find the first instance of whitespace.  Split label and data.
               pos = index(trimmed_buffer,' ')
               label = trim(adjustl(trimmed_buffer(1:pos)))
               buffer = trim(adjustl(trimmed_buffer(pos+1:)))

               if (label == trim(string)) then
                  read(buffer, *, iostat=ierr) val
                  if (dbg) write(*, '(a32, 2x, 1pes40.16e3)') 'Read ' // trim(label) // ': ', val
                  exit
               end if
            end do

            close(iounit)

         end subroutine read_parameter


         ! get number of natal kicks saved in the natal-kick file
         integer function number_of_kicks(filename, ierr)
            character(len=*) :: filename
            integer, intent(out) :: ierr
            integer :: iounit
            integer :: Nsim

            Nsim = 0
            open(newunit=iounit, file=trim(filename), status='old', action='read', iostat=ierr)
            if (ierr /= 0) stop 'could not open ' // trim(filename)
            do
               if (Nsim > max_num_kicks) then
                  write(*,'(a, 2x, i12)') 'reach max number of kicks that can be simulated of', Nsim
                  stop 'modify value of max_num_kicks on run_binary_extras.f and re-compile to have more kicks'
               end if

               read(iounit,*,iostat=ierr)
               if (ierr /= 0) exit
               Nsim = Nsim + 1
            end do
            ierr = 0
            close(iounit)

            number_of_kicks = Nsim

         end function number_of_kicks
