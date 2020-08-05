

      module bin2dco_io

      implicit none

      contains

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
