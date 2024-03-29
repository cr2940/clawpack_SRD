C     -*- fortran -*-
C     This file is autogenerated with f2py (version:2)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pywrapphilim (philimf2pywrap, a, b, meth)
      external philim
      double precision a
      double precision b
      integer meth
      double precision philimf2pywrap, philim
      philimf2pywrap = philim(a, b, meth)
      end


      subroutine f2pyinitcomxyt(setupfunc)
      external setupfunc
      double precision dtcom
      double precision dxcom
      double precision dycom
      double precision tcom
      integer icom
      integer jcom
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
      call setupfunc(dtcom,dxcom,dycom,tcom,icom,jcom)
      end


