#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2011 Marcus G. Martin                       *
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
      function twh_vembed(lfailure,ntii,density)
c     ******************************************************************
c     * computes the embedding function for the EAM potential          *
c     *                                                                *
c     * originally written 11-19-2003 by  M.G. Martin                  *
c     * last update 08-09-2011 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_VEMBED
#define FUNCTION_EXTRACTEMBED
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure
      integer ntii
      double precision density
c     --- local variables
      character*20 eamembedstyle
      double precision dsq,embedone,embedtwo,embedthr

      call twh_eam_embed_style(GLB_GET,ntii,eamembedstyle)
      if ( eamembedstyle .eq. 'table' ) then
         twh_vembed = twh_extractembed(lfailure,ntii,density)
         if ( lfailure ) return
      elseif ( eamembedstyle .eq. 'logarithmic') then
c        --- functional form from JJ Hoyt and J Zimmerman
c        --- f(density) = density*(a*ln(density)+b)
         call twh_eam_embed(GLB_GET,2,ntii,1,embedone)
         call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
         twh_vembed = density*((embedone*log(density))+embedtwo)
      elseif ( eamembedstyle .eq. 'squareroot') then
c        --- functional form from JJ Hoyt and J Zimmerman
c        --- f(density) = b*density - a*(sqrt(density))
         call twh_eam_embed(GLB_GET,2,ntii,1,embedone)
         call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
         twh_vembed = (embedtwo*density) - (embedone*dsqrt(density))
      elseif ( eamembedstyle .eq. 'power 0.5 and 2' ) then
c        --- Ackland two terms of square root and square
         call twh_eam_embed(GLB_GET,2,ntii,1,embedone)
         call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
         twh_vembed = embedone * dsqrt(density) 
     &        + embedtwo*density*density 
      elseif ( eamembedstyle .eq. 'power 0.5, 2, and 4' ) then
c        --- power series with 0.5, 2, and 4 powers
         dsq = density*density
         call twh_eam_embed(GLB_GET,2,ntii,1,embedone)
         call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
         call twh_eam_embed(GLB_GET,2,ntii,3,embedthr)
         twh_vembed = embedone * dsqrt(density)
     &        + embedtwo*dsq + embedthr*dsq*dsq
      else
         write(6,*) 'VEMBED: invalid eam_embed_style'
         write(6,*) 'ntii',ntii
     &        ,'eam_embed+style=',eamembedstyle,'='
         stop
      endif

c      write(6,*) 'finished VEMBED'

      return
      end
