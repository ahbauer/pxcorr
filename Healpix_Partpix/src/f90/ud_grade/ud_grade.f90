!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2010 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
module udgmod

  USE healpix_types
  USE fitstools, ONLY : getsize_fits, write_bintab, input_map
  USE pix_tools, ONLY : nside2npix, npix2nside
  USE misc_utils
  USE head_fits, ONLY : add_card, get_card, write_minimal_header
  USE udgrade_nr, ONLY : udgrade_ring, udgrade_nest
  USE extension, ONLY : getArgument, nArguments
  USE paramfile_io
  IMPLICIT none

  character(len=*), parameter :: CODE = "UD_GRADE"
  character(len=FILENAMELEN)  :: lcode

  private
  public :: udg_sub_s, udg_sub_d, CODE, lcode

  contains

    subroutine udg_sub_s(parafile)
      ! single precision implementation
      integer(I4B),               parameter  :: KMAP  = SP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'udg_sub_inc.f90'
    end subroutine udg_sub_s
    !------------------------------------


    subroutine udg_sub_d(parafile)
      ! double precision implementation
      integer(I4B),               parameter  :: KMAP  = DP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'udg_sub_inc.f90'
    end subroutine udg_sub_d

end module udgmod
!====================================================================================

program ud_grade

  use healpix_types
  use extension, only: getArgument, nArguments
  use udgmod,    only: udg_sub_s, udg_sub_d, CODE, lcode
  use misc_utils,only: assert, fatal_error, strlowcase

  integer(i4b)                :: n_args, i
  character(len=FILENAMELEN)  :: arg
  logical(lgt)                :: do_double
  character(len=FILENAMELEN)  :: parafile

  ! count arguments, should be 0, 1 or 2
  n_args = nArguments()
  lcode = trim(strlowcase(CODE))
  call assert(n_args <= 2,' Usage: '//trim(lcode)//' [-s|--single|-d|--double] [parameter_file_name]')

  parafile = ''
  do_double = .false.

  ! parse arguments
  do i=1, n_args
     call getArgument(i, arg)
     arg = trim(adjustl(arg))
     if (arg(1:1) == '-') then
        select case (trim(arg))
        case ('-d', '--double')
           do_double = .true.
        case ('-s', '--single')
           do_double = .false.
        case default
           call fatal_error(CODE//': Invalid argument: '//trim(arg))
        end select
     else
        parafile = arg
     endif
  enddo

  ! start calculations
  if (do_double) then
     call udg_sub_d(parafile)
  else
     call udg_sub_s(parafile)
  endif

  stop
end program ud_grade


