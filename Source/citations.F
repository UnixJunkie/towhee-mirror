#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2011-2011 Marcus G. Martin                       *
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
      subroutine twh_citations(lfailure)
c     ******************************************************************
c     * writes out the suggested citations for a Towhee simulation     *
c     *                                                                *
c     * originally written 08-09-2001 by M.G. Martin                   *
c     * last modified 08-16-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure

c     --- local variables
      character*3 ensstring
      character*7 version
      character*20 filename
      integer citefile,citecount,numboxes,ensemble,rng_code,luxlevel

c     --- output message to the standard out (or redirect)
      write(6,'(a,1x,a)')
     &     'Please see towhee_citations for a list of suggested'
     &     ,'citations for this simulation'

c     --- open file
      citefile = 19
      filename = 'towhee_citations'
      call twh_openfile(lfailure,citefile,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')

c     --- initialize citecount
      citecount = 0

c     --- disclaimer
      write(citefile,'(a,a)') 'Suggested citations for the methods,'
     &     ,' algorithms and force fields used in this simulation'

c     --- The code
      write(citefile,*)
      citecount = citecount + 1
      call twh_version(GLB_GET,version)
      write(citefile,'(a40,1x,a7,1x,a17,i2,a1)') 
     &     'Simulations were performed using version',version
     &     ,'of MCCCS Towhee [',citecount,']'
      write(citefile,'(a,i2,a)')
     &     '[',citecount,'] http://towhee.sourceforge.net'

c     --- random number generators
      write(citefile,*)
      call twh_random_code(GLB_GET,rng_code)
      if ( rng_code .eq. RNG_RANLUX ) then
         call twh_random_luxlevel(GLB_GET,luxlevel)
c        --- the original description of the RANLUX algorithm
         citecount = citecount + 1
         write(citefile,'(a,i2,a1,i2,a1)')
     &        'pseudo-random numbers generated using RANLUX ['
     &        ,citecount,',',citecount+1,']'
         write(citefile,'(a1,i2,a1,1x,a)')
     &        '[',citecount,']'
     &        ,'M. Luscher;'
         write(citefile,'(a1,a,1x,a,a1)')
     &        char(34)
     &        ,'A portable high-quality random number generator'
     &        ,'for lattice field theory simulations'
     &        ,char(34)
         write(citefile,'(a)') 'Comp. Phys. Comm. 79 100-110 (1994)'
c        --- the implementation of the RANLUX algorithm
         citecount = citecount + 1
         write(citefile,'(a1,i2,a1,1x,a)')
     &        '[',citecount,']'
     &        ,'F. James;'
         write(citefile,'(a1,a,1x,a,a1)')
     &        char(34)
     &        ,'RANLUX: A Fortran implementation of the high-quality'
     &        ,'pseudorandom number generator of Luscher'
     &        ,char(34)
         write(citefile,'(a)') 'Comp. Phys. Comm. 79 111-114 (1994)'
      elseif ( rng_code .eq. RNG_DX_1597_2_7 ) then
         citecount = citecount + 1
         write(citefile,'(a,i2,a1)')
     &        'pseudo-random numbers generated using DX-1597-2-7 ['
     &        ,citecount,']'
         write(citefile,'(a1,i2,a1,1x,a)')
     &        '[',citecount,']'
     &        ,'L.-Y. Deng;'
         write(citefile,'(a1,a,1x,a,a1)')
     &        char(34)
     &        ,'Efficient and Portable Multiple Recursive Generators'
     &        ,'of Large Order'
     &        ,char(34)
         write(citefile,'(a,1x,a)') 'ACM Transactions on Modeling and'
     &        ,'Computer Simulation 15 1-13 (2005)'
      elseif ( rng_code .eq. RNG_MRG32K3A ) then
         citecount = citecount + 1
         write(citefile,'(a,i2,a1)')
     &        'pseudo-random numbers generated using MRG32k3a ['
     &        ,citecount,']'
         write(citefile,'(a1,i2,a1,1x,a,a1,a)')
     &        '[',citecount,']'
     &        ,'P. l',char(39),'Ecuyer;'
         write(citefile,'(a1,a,1x,a,a1)')
     &        char(34)
     &        ,'Good Parameters and Implementations for Combined'
     &        ,'Multiple Recursive Random Number Generators'
     &        ,char(34)
         write(citefile,'(a)') 'Operations Research 47 159-164 (1999)'
      elseif ( rng_code .eq. RNG_KISS99 ) then
         citecount = citecount + 1
         write(citefile,'(a,i2,a1)')
     &        'pseudo-random numbers generated using KISS ['
     &        ,citecount,']'
         write(citefile,'(a1,i2,a1,1x,a,1x,a)')
     &        '[',citecount,']'
     &        ,'No literature references for KISS as it was published'
     &        ,'on internet newsgroups by G. Marsaglia'
      else
         write(citefile,*) 'Unknown random number generator'
      endif

c     --- Ensemble
      call twh_numboxes(GLB_GET,numboxes)
      call twh_ensemble(GLB_GET,ensemble)
      write(citefile,*)
      if ( ensemble .eq. ENS_NVT .and. numboxes .eq. 1 ) then
         citecount = citecount + 1
         write(citefile,'(a,i2,a1)')
     &        '... in the Canonical ensemble [',citecount,']'
         write(citefile,'(a1,i2,a1,1x,a,1x,a)')
     &        '[',citecount,']'
     &        ,'N. Metropolis; A. W. Rosenbluth; M. N. Rosenbluth'
     &        ,'A. H. Teller; E. Teller;'
         write(citefile,'(a1,a38,1x,a18,a1)') 
     &        char(34),'Equation of State Calculations by Fast'
     &        ,'Computing Machines',char(34)
         write(citefile,'(a)')
     &        ,'J. Chem. Phys. 21 1088-1092 (1953)'

      elseif ( ensemble .eq. ENS_NPT .and. numboxes .eq. 1 ) then
         citecount = citecount + 1
         write(citefile,'(a1,i2,a1,a,i2,a1)')
     &        '[',citecount,']'
     &        ,'... in the isobaric-isothermal ensemble [',citecount,']'
         write(citefile,'(a)') 'I. R. McDonald;'
         write(citefile,'(a1,a41,1x,a22,a1)')
     &        char(34),'NpT-ensemble Monte Carlo calculations for'
     &        ,'binary liquid mixtures',char(34)
         write(citefile,'(a)') 'Mol. Phys. 23 41-58 (1972)'

      elseif ( ensemble .eq. ENS_muVT ) then
         citecount = citecount + 1
         write(citefile,'(a,i2,a1)')
     &        '... in the Grand Canonical ensemble [',citecount,']'
         write(citefile,'(a1,i2,a1,1x,a)')
     &        '[',citecount,']'
     &        ,'G. E. Norman; V. S. Filinov;'
         write(citefile,'(a1,a,1x,a,a1)') 
     &        char(34)
     &        ,'Investigations of Phase Transitions by a Monte-Carlo'
     &        ,'Method'
     &        ,char(34)
         write(citefile,'(a)')
     &        ,'High Temperature 7 216-222 (1969)'

      elseif ( numboxes .gt. 1 .and.
     &        (ensemble .eq. ENS_NVT .or. ensemble .eq. ENS_NPT) ) then
         citecount = citecount + 1
         if ( ensemble .eq. ENS_NVT ) then
            ensstring = 'NVT'
         else
            ensstring = 'NpT'
         endif
         write(citefile,'(a10,1x,a3,a17,i2,a1,i2,a1)')
     &        '... in the',ensstring,'-Gibbs ensemble [',citecount
     &        ,char(44),citecount+1,']'
c        --- first Gibbs ensemble reference
         write(citefile,'(a1,i2,a1,1x,a)') 
     &        '[',citecount,']'
     &        , 'A. Z. Panagiotopoulos;'
         write(citefile,'(a1,a,1x,a,1x,a,a1)')
     &        char(34),'Direct determination of phase coexistence'
     &        ,'properties of fluids by Monte Carlo simulation in'
     &        ,'a new ensemble',char(34)
         write(citefile,'(a)')
     &        'Mol. Phys. 61 813-826 (1987)'
c        --- update citecount for the extra citation
         citecount = citecount + 1
c        --- second Gibbs ensemble reference
         write(citefile,'(a1,i2,a1,1x,a,1x,a)') 
     &        '[',citecount,']'
     &        , 'A. Z. Panagiotopoulos; N. Quirke; M. Stapleton;'
     &        ,'D. J. Tildesley;'
         write(citefile,'(a1,a,1x,a,1x,a,a1)')
     &        char(34)
     &        ,'Phase equilibria by simulation in the Gibbs ensemble:'
     &        ,'Alternative derivation, generalization and'
     &        ,'application to mixture and membrane equilibria'
     &        ,char(34)
         write(citefile,'(a)')
     &        'Mol. Phys. 63 527-545 (1988)'

      else
         write(citefile,*) 'citations: unknown ensemble'
      endif


c     --- close file
      close(citefile)

      return
      end
