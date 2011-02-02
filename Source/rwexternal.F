#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2007-2008 Marcus G. Martin                       *
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
      subroutine twh_rwexternal(lfailure,rwflag,ifile)
c     ******************************************************************
c     * reads or writes all input parameters from towhee_input that are*
c     * related to the external potential style                        *
c     * rwflag is either READ_FLAG or WRITE_FLAG                       *
c     *                                                                *
c     * rewritten from readquantum 01-12-2007 by M.G. Martin           *
c     * last modifed 11-10-2008 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ifile,rwflag
c     --- local variables
      character*15 idname,rwformat
      character*20 externalcode
      character*50 goodstring

c     --- check rwflag
      if ( rwflag .ne. READ_FLAG .and. rwflag .ne. WRITE_FLAG ) then
         write(6,*) 'RWEXTERNAL: unknown rwflag'
         write(6,*) rwflag
         lfailure = .true.
         return
      endif
c     --- set subroutine idname
      idname = 'rwexternal'

      goodstring = 'external_code'
      if ( rwflag .eq. READ_FLAG ) then
         call twh_read_labeled_string(ifile,goodstring,.true.,idname
     &        ,externalcode)
         call twh_external_code(GLB_SET,externalcode)
      else
         rwformat = '(a11)'
         write(ifile,rwformat) goodstring
         call twh_external_code(GLB_GET,externalcode)
         write(ifile,'(a1,a20,a1)') char(39),externalcode,char(39)
      endif

      if ( externalcode .eq. 'lcao' ) then
         call twh_readquantum(lfailure,rwflag,ifile)
         if ( lfailure ) return
      else
         write(6,*) 'RWEXTERNAL: invalid external_code'
         lfailure = .true.
         return
      endif

      return
      end
