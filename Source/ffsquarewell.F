#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2020 Marcus G. Martin                       *
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
      subroutine twh_ffsquarewell(lfailure,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes)
c     ******************************************************************
c     * sets up generic force field parameters for the square well     *
c     * potential                                                      *
c     * has now been modified so adding a new nonbonded type is enough *
c     * to set up all of the other interactions                        *
c     * bonds are set to the arithmetic average of the hard sphere of  *
c     *       the nonbonded types                                      *
c     * angles and torsions are strictly nonbonded interactions        *
c     *                                                                *
c     * originally written 04-08-2004 by M.G. Martin                   *
c     * last modified 03-31-2020 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
c     --- character scalars
      character*(FFNAMELEN) ffname
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
c     --- character arrays
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
c     --- integer scalars
      integer ii,jj,itype,jtype,type
      integer nvslocal,nbslocal,ntslocal
c     --- local double precision scalars
      double precision masslocal,nbc_ii,nbc_jj
c     --- local double precision arrays
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)

      ffname = 'SquareWell'
c     --- square well potential
      call twh_classical_potential(GLB_SET,'Square Well')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

c     --- initialize types
      type = 0

c     --- type a
      type = type + 1
      masslocal = 1.0d0
      elelocal = 'Sw'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='type_a'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
      nbclocal(2) = 1.5d0
      nbclocal(3) = 1.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- type b
      type = type + 1
      masslocal = 1.0d0
      elelocal = 'Sw'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='type_b'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.0d0
      nbclocal(2) = 3.0d0
      nbclocal(3) = 1.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     --- bond potentials
c     --- fixed bond length that is the average of the inner hard 
c     --- sphere radii

      do itype = 1,ntypes
         ii = twh_getnbtype(itype,itype)
         do jtype = itype,ntypes
            jj = twh_getnbtype(jtype,jtype)
c           --- increment bond type
            type = type + 1
            call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c           --- set to the arithmetic average of the hard radii
            call twh_nbcoeff(GLB_GET,1,ii,nbc_ii)
            call twh_nbcoeff(GLB_GET,1,jj,nbc_jj)
            vclocal(0) = 0.5d0*(nbc_ii + nbc_jj)
            call twh_set_vibcoeffs(type,vclocal)
            call twh_vibff(GLB_SET,type,ffname)
            nvslocal = 1
c           --- use the bond names of those types
            call twh_cbname(GLB_GET,itype,vnlocal(1,1))
            call twh_cbname(GLB_GET,jtype,vnlocal(1,2))
            call twh_set_vibnames(type,nvslocal,vnlocal)
         enddo
      enddo

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     --- angle potentials nonbond interaction for all
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_NONBOND)
      call twh_bendff(GLB_SET,type,ffname)
c     wild - wild - wild
      nbslocal = 1
      bnlocal(1,1) = 'wild'
      bnlocal(1,2) = 'wild'
      bnlocal(1,3) = 'wild'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     --- torsion potentials nonbond interaction for all
      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_NB)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- wild - wild - wild - wild
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'wild'
      tnlocal(1,3) = 'wild'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
