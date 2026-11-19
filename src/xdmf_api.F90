!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
!----------------------------------------------------------------------
!                M O D U L E     X D M F   A P I
!----------------------------------------------------------------------
! Single point of access to the XDMF3 Fortran API. XDMF ships only the
! Xdmf.f include file (API constants and the result types of its
! functions), with no Fortran interfaces, so the API routines are
! declared external here. Procedures that use XDMF import what they
! need from this module instead of including Xdmf.f themselves.
!----------------------------------------------------------------------
module mod_xdmf_api
   implicit none

   public

   include 'Xdmf.f'

   ! Functions typed by Xdmf.f
   ! allow(C091) -- XDMF provides no Fortran interfaces to declare instead
   external :: xdmfaddattribute, xdmfaddcoordinate, xdmfaddfunctionasattribute, &
      xdmfaddfunctionascoordinate, xdmfaddfunctionasset, xdmfaddinformation, &
      xdmfaddset, xdmfaddsubsetasattribute, xdmfaddsubsetascoordinate, &
      xdmfaddsubsetasset, xdmfgetmaxopenedfiles, xdmfsetbrick, &
      xdmfsetdimensions, xdmfsetfunctionasgeometry, xdmfsetfunctionastopology, &
      xdmfsetgeometry, xdmfsetorigin, xdmfsetsubsetasgeometry, &
      xdmfsetsubsetastopology, xdmfsettopology, xdmfstoremap

   ! Subroutines
   ! allow(C091) -- XDMF provides no Fortran interfaces to declare instead
   external :: xdmfinit, xdmfinithdf5, xdmfclose, xdmfwrite, xdmfaddgrid, &
      xdmfaddgridcollection, xdmfclosegridcollection, xdmfsettime, &
      xdmfaddpreviousattribute, xdmfaddpreviousinformation, &
      xdmfaddinformationarray, xdmfretrievenuminformation, &
      xdmfretrievegridcollectionnumgrids, xdmfsetpreviousgeometry, &
      xdmfsetprevioustopology

   ! Subroutines used to read XDMF input (mesh.F, nodalattr.F)
   ! allow(C091) -- XDMF provides no Fortran interfaces to declare instead
   external :: xdmfopendomaingrid, xdmfopenset, xdmfread, xdmfretrieveattributecenter, &
      xdmfretrieveattributename, xdmfretrieveattributesize, &
      xdmfretrieveattributetype, xdmfretrieveattributevalues, &
      xdmfretrieveattributevaluetype, xdmfretrievedomaingridname, &
      xdmfretrievegeometrynumpoints, xdmfretrievegeometrysize, &
      xdmfretrievegeometrytype, xdmfretrievegeometryvalues, &
      xdmfretrievegeometryvaluetype, xdmfretrieveinformation, &
      xdmfretrievenumattributes, xdmfretrievenumsets, &
      xdmfretrievesetnumproperties, xdmfretrievesetproperty, &
      xdmfretrievesetsize, xdmfretrievesetvalues, xdmfretrievesetvaluetype, &
      xdmfretrievetopologynumelements, xdmfretrievetopologysize, &
      xdmfretrievetopologytype, xdmfretrievetopologyvalues, &
      xdmfretrievetopologyvaluetype

end module mod_xdmf_api
