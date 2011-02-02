c     $Id: ffcheck.F,v 1.7 2006/02/09 20:47:33 mattwycz Exp $
#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2005 Marcus G. Martin                       *
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
      subroutine ffcheck(style,type)
c     ******************************************************************
c     * checks to see if we have exceeded any arrays                   *
c     *                                                                *
c     * originally written 12-19-2002 by M.G. Martin                   *
c     * last modified 10-18-2005 M.G. Martin                           *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      integer style,type

      if ( style .eq. 0 ) then
c        --- nonbond
         write(6,*) 'Nonbonded types:',type
         if ( type .gt. NNTYPE ) then
            write(6,*) 'type exceeds NNTYPE'
            stop
         endif
      elseif ( style .eq. 1 ) then
c        --- vibration
         write(6,*) 'Bond types:',type
         if ( type .gt. TVIBMAX ) then
            write(6,*) 'type exceeds TVIBTYPE'
            stop
         endif
      elseif ( style .eq. 2 ) then
c        --- angles
         write(6,*) 'Angle types:',type
         if ( type .gt. TBENMAX ) then
            write(6,*) 'type exceeds TBENMAX'
            stop
         endif
      elseif ( style .eq. 3 ) then
c        --- torsions
         write(6,*) 'Torsion types:',type
         if ( type .gt. TTORMAX ) then
            write(6,*) 'type exceeds TTORMAX'
            stop
         endif
      elseif ( style .eq. 4 ) then
c        --- impropers
         write(6,*) 'Improper types:',type
         if ( type .gt. TIMPMAX ) then
            write(6,*) 'type exceeds TIMPMAX'
            stop
         endif
      elseif ( style .eq. 5 ) then
         write(6,*) 'Angle-angle types:',type
         if ( type .gt. TAAMAX ) then
            write(6,*) 'type exceeds TAAMAX'
            stop
         endif
      elseif ( style .eq. 6 ) then
c        --- bond increment
         write(6,*) 'Bond Increment types:',type
         if ( type .gt. TBIMAX ) then
            write(6,*) 'type exceeds TBITYPE'
            stop
         endif
      elseif ( style .eq. 7 ) then
c        --- one-five
         write(6,*) 'One-Five types:',type
         if ( type .gt. TOFMAX ) then
            write(6,*) 'type exceeds TOFMAX'
            stop
         endif
      elseif ( style .eq. 8 ) then
         write(6,*) 'H-bond types:',type
         if ( type .gt. THBONDMAX ) then
            write(6,*) 'type exceeds THBONDMAX'
            stop
         endif
      else
         write(6,*) 'FFCHECK: unknown style'
         stop
      endif

c     --- reset type
      type = 0

      return
      end
