

         subroutine set_default_bin2dco_controls
            include 'bin2dco_controls.defaults'
         end subroutine set_default_bin2dco_controls


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
         integer function number_of_kicks(filename, skip_rows, ierr)
            character(len=*) :: filename
            integer, intent(in) :: skip_rows
            integer, intent(out) :: ierr
            integer :: iounit
            integer :: Nsim

            Nsim = 0
            open(newunit=iounit, file=trim(filename), status='old', action='read', iostat=ierr)
            if (ierr /= 0) stop 'could not open ' // trim(filename)

            do k = 1, skip_rows
               read(iounit,*,iostat=ierr)
               if (ierr /= 0) return
            end do

            do
               if (Nsim > max_num_kicks) then
                  write(*,'(a, 2x, i12)') 'reach max number of kicks', Nsim
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


         ! read from file natal kick information: id, velocity, angles
         subroutine read_natal_kick(filename, line_number, name_id, porb, a, ecc, ierr)
            use const_def, only: dp, strlen
            character(len=*) :: filename
            integer, intent(in) :: line_number
            character(len=*), intent(out) :: name_id
            real(dp), intent(out) :: porb, a, ecc
            integer, intent(out) :: ierr

            integer :: iounit
            integer :: i, j
            character(len=strlen), dimension(4) :: str_entry
            character(len=strlen) :: buffer, trimmed_buffer, string
            logical :: dbg = .false.

            include 'formats'

            ! open up natal kicks file
            open(newunit=iounit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(filename)
               return
            end if

            i = 1
            do while (ierr == 0)

               ! skip first (line_number - 1) lines
               do
                  if (i > line_number - 1) exit
                  read(iounit, *, iostat=ierr)
                  if (ierr > 0) then
                     if (dbg) write(*,'(a)') 'problem in reading ' // trim(filename)
                     return
                  else if (ierr < 0) then
                     if (dbg) write(*,'(a)') 'end of file ' // trim(filename) // ' reached'
                     return
                  else
                     i = i + 1
                  end if
               end do

               read(iounit, '(A)', iostat=ierr) buffer
               read(buffer,*) (str_entry(i), i=1,4)

               do j = 1, 4
                  string = trim(adjustl(str_entry(j)))
                  select case (j)
                  case (1)
                     read(string, *, iostat=ierr) name_id
                     if (dbg) write(*,'(a55,99(a55))') 'Read name_id:', name_id
                  case (2)
                     read(string, *, iostat=ierr) porb
                     if (dbg) write(*,1) 'Read porb:', porb
                  case (3)
                     read(string, *, iostat=ierr) a
                     if (dbg) write(*,1) 'Read a:', a
                  case (4)
                     read(string, *, iostat=ierr) ecc
                     if (dbg) write(*,1) 'Read ecc:', ecc
                  end select
               end do
               ierr = -1
            end do

            ierr = 0
            close(iounit)

         end subroutine read_natal_kick
