c     $Id: wmultibody.F,v 1.12 2008/10/14 14:09:46 marcus_martin Exp $
#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2008 Marcus G. Martin                       *
c     * see the file license.gpl for the full license information      *
c     *                                                                *
c     * This program is free software; you can redistribute it and/or  *
c     * modify it under the terms of the GNU General Public License    *
c     * as published by the Free Software Foundation; either version 2 *
c     * of the License, or (at your option) any later version.         *
c     *                                                                *
c     * This program is distributed in the hope that it will be useful,*
c     * but WITHOUT ANY WARRANTY; without even the implied warranty of *
c     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
c     * GNU General Public License for more details.                   *
c     *                                                                *
c     * You should have received a copy of the GNU General Public      *
c     * License along with this program; if not, write to the Free     *
c     * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
c     * MA  02111-1307, USA.                                           *
c     ******************************************************************
      function twh_wmultibody( rijsq,ntii,ntjj,ichain,jchain)
c     ******************************************************************
c     * Computes the intermolecular pair virial function divided by r^2*
c     * as defined by equation 2.60 in Allen and Tildesley             *
c     * w(r) = r dv(r)/dr                                              *
c     * w(r)/r^2 = (dv(r)/dr) /r                                       *
c     * does this for the embedded atom multibody portion using the    *
c     * chain rule d/dr F(rho(r)) = dF/drho * drho/dr                  *
c     * have to do this twice, once for i and once for -j              *
c     *                                                                *
c     * originally written 03-25-2005 by M.G. Martin                   *
c     * last update 10-13-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_WMULTIBODY
#include "functions.h"
c     --- variables passed to/from the function
      integer ntii,ntjj,ichain,jchain
      double precision rijsq
c     --- local variables
c     --- logical scalars
      logical lcontinue
c     --- character string scalars
      character*20 eamembedstyle,eamdensstyle
c     --- integer scalars
      integer itrial,icompute,icoeff,trialchain,iother,cpnumber
c     --- double precision scalars
      double precision x,rij,density,eamdensone,eamdenstwo
c     --- double precision arrays
      double precision embedone,embedtwo,embedthr
      double precision embedprime,rhoprime
      dimension embedprime(2),rhoprime(2)

c      write(6,*) 'starting WMULTIBODY'

c     --- initialize wmultibody
      twh_wmultibody = 0.0d0

      call twh_cp_number(GLB_GET,cpnumber)
      if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c        --- compute rij
         rij = dsqrt(rijsq)
c        --- the embedding and density parts
         do itrial = 1,2
            rhoprime(itrial) = 0.0d0
            embedprime(itrial) = 0.0d0
            if ( itrial .eq. 1 ) then
               icompute = ntii
               iother = ntjj
               trialchain = ichain
            else
               icompute = ntjj
               iother = ntii
               trialchain = jchain
            endif
c           --- density derivative
            call twh_eam_dens_style(GLB_GET,iother,icompute
     &           ,eamdensstyle)
            if ( eamdensstyle .eq. 'Ackland cubic sum' )then
c              --- cubic sum with heaviside
               call twh_eam_ndens(GLB_GET,iother,icompute,icoeff)
               if ( icoeff .gt. 0 ) then
                  lcontinue = .true.
               else
                  lcontinue = .false.
               endif
               do while ( lcontinue )
                  call twh_eam_dens(GLB_GET,1,iother,icompute,icoeff
     &                 ,eamdensone)
                  x = eamdensone - rij
                  call twh_eam_dens(GLB_GET,2,iother,icompute,icoeff
     &                 ,eamdenstwo)
                  rhoprime(itrial) = rhoprime(itrial) +
     &                 eamdenstwo*(-3.0d0)*x*x
                  icoeff = icoeff - 1
                  if ( icoeff .eq. 0 ) then
                     lcontinue = .false.
                  else
                     call twh_eam_dens(GLB_GET,1,iother,icompute,icoeff
     &                    ,eamdensone)
                     if ( eamdensone .le. rij ) lcontinue = .false.
                  endif
               enddo
            elseif ( eamdensstyle .eq. 'none' ) then
c              --- no density contribution from this pair interaction
               rhoprime(itrial) = 0.0d0
            endif
c           --- embedding derivative
c           --- set the density
            call twh_eam_rho(GLB_GET,CRD_REAL,trialchain,density)
            call twh_eam_embed_style(GLB_GET,icompute,eamembedstyle)
            if ( eamembedstyle .eq. 'power 0.5 and 2' ) then
c              --- power series with 0.5, and 2
               if ( density .ne. 0.0d0 ) then
                  call twh_eam_embed(GLB_GET,2,icompute,1,embedone)
                  call twh_eam_embed(GLB_GET,2,icompute,2,embedtwo)
                  embedprime(itrial) = (0.5d0)*embedone/dsqrt(density) 
     &                 +(2.0d0)*embedtwo*density
               endif
            elseif ( eamembedstyle .eq. 'power 0.5, 2, and 4' ) then
c              --- power series with 0.5, 2, and 4 powers
               if ( density .ne. 0.0d0 ) then
                  call twh_eam_embed(GLB_GET,2,icompute,1,embedone)
                  call twh_eam_embed(GLB_GET,2,icompute,2,embedtwo)
                  call twh_eam_embed(GLB_GET,2,icompute,3,embedthr)
                  embedprime(itrial) =
     &                 (0.5d0)*embedone/dsqrt(density)
     &                 +(2.0d0)*embedtwo*density
     &                 +(4.0d0)*embedthr*(density**3)
               endif
            endif
         enddo
c        --- compute wmultibody (i terms + jterms)/r
         twh_wmultibody = ( rhoprime(1)*embedprime(1) 
     &        + rhoprime(2)*embedprime(2) )/rij
      endif

      return
      end
