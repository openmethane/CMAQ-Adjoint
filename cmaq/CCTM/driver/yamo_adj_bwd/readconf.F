!-----------------------------------------------------------------------
! Function:
!    The function reads the namelist file and sets the environmental variables
!
! Preconditions:
!    The name of the namelist file is read from environmental variable CTM_NAMELIST.
!    Default value is namelist.cmaq.
!
! Revision history:
!   Aug 2011 Jaroslav Resler, ICS Prague: initial version
!-----------------------------------------------------------------------


      subroutine readconf(namelist)

        character(*), intent(in)   :: namelist

        character(16)               :: pname = 'readconf'
        integer                     :: fnum
        integer                     :: ierr
        character(255)              :: line, vname, value
        integer                     :: maxfiles = 1000
        logical                     :: opened
        integer                     :: i


        do fnum = 10, maxfiles
          inquire ( unit=fnum, opened=opened, iostat=ierr )
          if ( ierr .ne. 0 ) then
            call m3exit(pname, 0, 0, 'Error testing files', ierr)
            return
          end if
          if ( .not. opened ) then
            goto 11
          end if
        end do
        call m3exit(pname, 0, 0, 'Max. number of opened files exceeded.', -1)
        return

11      continue
        open(unit=fnum, file=namelist, status='old', action='read', iostat=ierr)
        if (ierr .ne. 0) then
          line = 'Cannot open namelist file: '//trim(namelist)
          call m3warn( pname, 0, 0, line)
          return
        end if

        do while (.true.)
          read(fnum,'(a)', err=111, end=112) line
          i = index(line,'!')
          if ( i.gt.0 ) then
            line = trim(adjustl(line(:i-1)))
          end if
          i = index(line,'=')
          if ( i.gt.1 ) then
            vname = adjustl(line(:i-1))
            value = adjustl(line(i+1:))
            if (len_trim(vname).gt.0) then
               call setenvvar( trim(vname), trim(value) )
            end if
          end if
        end do

111     continue
        call m3exit(pname, 0, 0, 'Error reading namelist file', -1)
        return

112     continue
        close(fnum)
        return


      end subroutine readconf
