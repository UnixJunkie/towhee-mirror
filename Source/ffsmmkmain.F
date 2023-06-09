#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2014 Marcus G. Martin                       *
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
      subroutine twh_ffsmmkmain(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,noftypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the mund1996     *
c     * flexible varient of the SMMK Main force field                  *
c     *                                                                *
c     * originally written 12-14-2004 by M.G. Martin                   *
c     * last modified 08-13-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,noftypes
c     --- local variables
      integer ii,type,isame
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
      character*(FFNAMELEN) ofnlocal
      dimension ofnlocal(MAXOFNAME)

      ffname = 'SMMKmain'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Geometric')
      type = 0

c     ******************************************************************
c     * SMMK Main                                                      *
c     ******************************************************************
c     References:
c     [mp1997]
c     J.I. Siepmann; M.G. Martin; C.J. Mundy; M.L. Klein;
c     "Intermolecular potentials for branched alkanes and the vapour-
c     liquid phase equilibria of n-heptane, 2-methylhexane, and 
c     3-ethylpentane", Mol. Phys. 90, 687-693 (1997)
c     [faraday1996]
c     C.J. Mundy; S. Balasubramanian; K. Bagchi; J.I. Siepmann;
c     M.L. Klein;
c     "Equilibrium and non-equilibrium simulation studies of fluid
c     alkanes in bulk and at interfaces", Faraday Discuss. 104,
c     17-36 (1996)

c     --- we are not using the 1-4 nonbonded terms so set them to zero
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0

c     (CH3long) CH3 methyl group on an n-alkane or alkyl side chain
c     of propyl or longer [mp1997]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.93d0
      nbclocal(2) = 114.0d0
      masslocal = 12.001d0 + 3.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3long'
      cbnlocal = 'C'
      canlocal = 'CH3'
      ctnlocal = 'CH3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CH3ethyl) CH3 methyl group on propane or an ethyl side chain
c     [mp1997]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.93d0
      nbclocal(2) = 97.0d0
      masslocal = 12.001d0 + 3.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3ethyl'
      cbnlocal = 'C'
      canlocal = 'CH3'
      ctnlocal = 'CH3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CH3methyl) CH3 methyl group on ethane or a methyl side chain
c     [mp1997]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.93d0
      nbclocal(2) = 78.0d0
      masslocal = 12.001d0 + 3.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3methyl'
      cbnlocal = 'C'
      canlocal = 'CH3'
      ctnlocal = 'CH3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CH2) CH2 methylene group [mp1997]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.93d0
      nbclocal(2) = 47.0d0
      masslocal = 12.001d0 + 2.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2'
      cbnlocal = 'C'
      canlocal = 'CH2'
      ctnlocal = 'CH2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CH) CH methine group [mp1997]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.85d0
      nbclocal(2) = 32.0d0
      masslocal = 12.001d0 + 1.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH'
      cbnlocal = 'C'
      canlocal = 'CH'
      ctnlocal = 'CH'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (C) C quaternary carbon [faraday1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      masslocal = 12.001d0 
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C'
      cbnlocal = 'C'
      canlocal = 'C'
      ctnlocal = 'C'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * harmonic bond lengths                                          *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.54d0
      vclocal(1) = (452900.0d0)*(0.5d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     [faraday1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * simple harmonic angles                                         *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.0d0
      bclocal(1) = 0.5d0*(62500.0d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     angle centered on a CH2 [faraday1996]
      bnlocal(1,1) = 'wild'
      bnlocal(1,2) = 'CH2'
      bnlocal(1,3) = 'wild'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.0d0
      bclocal(1) = 0.5d0*(62500.0d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     angle centered on a CH [faraday1996]
      isame = isame + 1
      bnlocal(isame,1) = 'wild'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'wild'
c     angle centered on a C [faraday1996]
      isame = isame + 1
      bnlocal(isame,1) = 'wild'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'wild'
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * opls style torsions                                            *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_O_OPLS_C)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 0.0d0
      tclocal(1) = 355.03d0
      tclocal(2) = -68.19d0
      tclocal(3) = 791.32d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     -CH2-CH2- [faraday1996]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CH2'
      tnlocal(isame,3) = 'CH2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_O_OPLS_C)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = -251.06d0
      tclocal(1) = 428.73d0
      tclocal(2) = -111.85d0
      tclocal(3) = 441.27d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     -CH2-CH- [faraday1996]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CH2'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_O_OPLS_C)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 0.0d0
      tclocal(1) = 0.0d0
      tclocal(2) = 0.0d0
      tclocal(3) = 461.29d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     -CH2-C- [faraday1996]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CH2'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * special 1-5 interactions                                       *
c     ******************************************************************

      type = type + 1
      call twh_ofstyle(GLB_SET,type,1)
      call twh_offf(GLB_SET,type,ffname)
c     --- different Lennard-Jones parameters for this interaction
      call twh_ofcoeff(GLB_SET,type,1,3.62d0)
      call twh_ofcoeff(GLB_SET,type,2,40.0d0)
c     different nonbonded parameters for a 1-5 interaction between
c     CH3methyl groups along a molecule like 
c     2,2,4,4,6,8,8-heptamethylnonane
      ofnlocal(1) = 'CH3methyl'
      ofnlocal(2) = 'wild'
      ofnlocal(3) = 'wild'
      ofnlocal(4) = 'wild'
      ofnlocal(5) = 'CH3methyl'
      call twh_set_ofnames(type,ofnlocal)

      noftypes = type
      call twh_ffcheck(lfailure,7,type)
      if ( lfailure ) return

      return
      end
