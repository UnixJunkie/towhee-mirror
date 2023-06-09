#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2018 Marcus G. Martin                       *
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
      subroutine twh_recip(lfailure,ibox,type)
c     ******************************************************************
c     * updates the reciprocal ewald-sum.                              *
c     *                                                                *
c     * rewritten on 06-25-1999 by Bin Chen.                           *
c     * last modified 04-03-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to and from the subroutine
      logical lfailure
      integer ibox,type
c     --- local variables
      logical lewald
      integer ncount
      double precision kx,ky,kz,prefact,ssumr,ssumi

c     --- bail out if not using ewald sum
      call twh_lewald(GLB_GET,lewald)
      if ( .not. lewald ) return
c     --- bail out if this is the ideal chain box
      if ( ibox .eq. 0 ) return
c     --- get ncount
      call twh_ewald_numvect(GLB_GET,ibox,ncount)
      if ( type .eq. RCP_ACCEPT_NEW_VOL )  then
c        --- update for a global change
         call twh_ewald_ssumr(GLB_STORE,CRD_N_TO_R,ncount,ibox,ssumr)
         call twh_ewald_ssumi(GLB_STORE,CRD_N_TO_R,ncount,ibox,ssumi)
         call twh_ewald_kvector(GLB_STORE,CRD_T_TO_R,ncount,ibox
     &        ,kx,ky,kz)
         call twh_ewald_prefact(GLB_STORE,CRD_T_TO_R,ncount,ibox
     &        ,prefact)
      elseif ( type .eq. RCP_ACCEPT_NEW ) then
c        --- update the reciprocal space k vectors, volume didn't change
         call twh_ewald_ssumr(GLB_STORE,CRD_N_TO_R,ncount,ibox,ssumr)
         call twh_ewald_ssumi(GLB_STORE,CRD_N_TO_R,ncount,ibox,ssumi)
      elseif ( type .eq. RCP_STORE ) then
c        --- store the reciprocal space k vectors         
         call twh_ewald_kvector(GLB_STORE,CRD_R_TO_T,ncount,ibox
     &        ,kx,ky,kz)
         call twh_ewald_prefact(GLB_STORE,CRD_R_TO_T,ncount,ibox
     &        ,prefact)
         call twh_ewald_ssumr(GLB_STORE,CRD_R_TO_O,ncount,ibox,ssumr)
         call twh_ewald_ssumi(GLB_STORE,CRD_R_TO_O,ncount,ibox,ssumi)

      elseif ( type .eq. RCP_RESTORE ) then
c        --- restore the reciprocal space k vectors         
         call twh_ewald_kvector(GLB_STORE,CRD_T_TO_R,ncount,ibox
     &        ,kx,ky,kz)
         call twh_ewald_prefact(GLB_STORE,CRD_T_TO_R,ncount,ibox
     &        ,prefact)
         call twh_ewald_ssumr(GLB_STORE,CRD_O_TO_R,ncount,ibox,ssumr)
         call twh_ewald_ssumi(GLB_STORE,CRD_O_TO_R,ncount,ibox,ssumi)
      else
c        --- invalid type
         write(6,*) 'RECIP: invalid type flag:',type
         lfailure = .true.
         return
      endif

      return 
      end
