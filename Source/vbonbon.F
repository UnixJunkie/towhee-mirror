#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *c     * Copyright (C) 2000-2020 Marcus G. Martin                       *
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
      function twh_vbonbon(imolty,iunit,index,length)
c     ******************************************************************
c     * computes the IJ-KL bond-bond cross term for a given bond-bond  *
c     * type and the two bond lengths                                  *
c     * the bond-bond cross type is determined from the molecule type, *
c     * the unit number and the bond-bond index                        *
c     *                                                                *
c     * originally written 02-28-2000 by M.G. Martin                   *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_VBONBON
#include "functions.h"
c     --- variables passed to/from the function
      integer imolty,iunit,index
      double precision length
      dimension length(2)
c     --- local variables
      integer style,type,iorder,itemp
      double precision vtone,vttwo,benprefact

      call twh_itbnbn(GLB_GET,imolty,iunit,index,type)
      if ( type .gt. 0 ) then
         iorder = 1
      else
         type = -type
         iorder = -1
      endif
      call twh_anglestyle(GLB_GET,type,style)

      if ( style .eq. ANG_STYLE_COMPASS_AUTO ) then
c        --- Compass bond-bond cross term
c        --- energy = k*(bond(1) - bond0(1))*(bond(2) - bond0(2))
         call twh_ibtbnbn1(GLB_GET,imolty,iunit,index,itemp)
         call twh_vibcoeff(GLB_GET,itemp,0,vtone)
         call twh_ibtbnbn2(GLB_GET,imolty,iunit,index,itemp)
         call twh_vibcoeff(GLB_GET,itemp,0,vttwo)
         call twh_bencoeff(GLB_GET,type,6,benprefact)
         twh_vbonbon = benprefact*(length(1)-vtone)*(length(2)-vttwo)
      elseif ( style .eq. ANG_STYLE_COMPASS_EX ) then
c        --- Compass bond-bond cross term with explicit distances
c        --- energy = k*(bond(1) - bond0(1))*(bond(2) - bond0(2))
         if ( iorder .eq. 1 ) then
            call twh_bencoeff(GLB_GET,type,9,vtone)
            call twh_bencoeff(GLB_GET,type,10,vttwo)
         else
            call twh_bencoeff(GLB_GET,type,10,vtone)
            call twh_bencoeff(GLB_GET,type,9,vttwo)
         endif
         call twh_bencoeff(GLB_GET,type,8,benprefact)
         twh_vbonbon = benprefact*(length(1)-vtone)*(length(2)-vttwo)
      else
c        --- these angle types should not have a bond-bond term
         stop 'bad angle type in vbonbon'
      endif

      return
      end
