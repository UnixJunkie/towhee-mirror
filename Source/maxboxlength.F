c     $Id: maxboxlength.F,v 1.7 2008/10/14 14:09:45 marcus_martin Exp $
#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2008 Marcus G. Martin                       *
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
      function twh_maxboxlength(ibox)
c     ******************************************************************
c     * finds the maximum dimension of a simulation box                *
c     *                                                                *
c     * originally written 08-13-2004 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_MAXBOXLENGTH
#include "functions.h"

c     --- variables passed to/from the subroutine
      integer ibox

c     --- local variables
      integer idim,jdim
      double precision boxlength,dvalue

      if ( ibox .eq. 0 ) then
c        --- ideal gas box, infinite length, set to zero
         twh_maxboxlength = 0.0d0
      else
         do jdim = 1,3
            boxlength = 0.0d0
            do idim = 1,3
               call twh_hmatrix(GLB_GET,ibox,idim,jdim,dvalue)
               boxlength = boxlength + dabs(dvalue)
            enddo
            if ( jdim .eq. 1 ) then
               twh_maxboxlength = boxlength
            elseif (boxlength .gt. twh_maxboxlength ) then
               twh_maxboxlength = boxlength
            endif
         enddo
      endif

      return
      end
