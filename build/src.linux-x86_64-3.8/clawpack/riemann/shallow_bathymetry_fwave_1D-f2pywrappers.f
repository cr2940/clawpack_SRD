C     -*- fortran -*-
C     This file is autogenerated with f2py (version:2)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pyinitcparam(setupfunc)
      external setupfunc
      real(kind=8) grav
      real(kind=8) dry_tolerance
      real(kind=8) sea_level
      common /cparam/ grav,dry_tolerance,sea_level
      call setupfunc(grav,dry_tolerance,sea_level)
      end

