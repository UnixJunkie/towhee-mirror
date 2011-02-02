#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000-2010 Marcus G. Martin                       *
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
      subroutine twh_atd_cbmc_setup(lfailure)
c     ******************************************************************
c     * sets up all of the arbitrary trial distribution information for*
c     * use in the configurational-bias routines                       *
c     *                                                                *
c     * split from checkstruc 05-18-2006 by M.G. Martin                *
c     * last modified 07-25-2010 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- debug logical
#define DEBUG_ATD_CBMC_SETUP 0
c     --- function declarations
#if DEBUG_ATD_CBMC_SETUP
#define FUNCTION_GAUSSPROB
#endif
#define FUNCTION_EXPON
#define FUNCTION_GET_ANGLESTYLE
#define FUNCTION_GET_BENCOEFF
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_ITBEN
#define FUNCTION_GET_INTOR
#define FUNCTION_LAWOFCOSINE
#define FUNCTION_VANGLE
#define FUNCTION_VBOND
#define FUNCTION_VTORSION
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
      character*30 peakstyle
      character*50 twobbias_style,twobbias_init,threebbias_init
      character*50 threebbias_style
      logical lovrlap,lfound,lrigid
c     --- local integer parameters
      integer maxextrema
      parameter (maxextrema = 2*(CBMC_MAXTORPEAKS+CBMC_MAXBENDPEAKS))
c     --- local integer scalars
      integer imolty,iunit,ivib,junit,itest,itype,ibend,kunit,itor
      integer nextrema,nmaxima,iextrema,lunit,ipeak,npeak
      integer imax,junitb,kunitb,icrossbend,jbend,kbend,ibondmax
      integer bij_type,bjk_type,iscan_bij,iscan_bjk
      integer type,nmolty,nunit,ibox,cbmc_bondgen_code,cbmc_bendgen_code
      integer cbmc_dihedgen_code,style,bonddist
#if DEBUG_ATD_CBMC_SETUP
      integer count,style_code,imgm,ivalue
#endif
c     --- local dp scalars
      double precision rmax,sum_x,sum_xx,length,theta,truenorm,arbnorm
      double precision energy,weight,sum_w,distone,disttwo,angcos
      double precision phione,phitwo,thetaone,thetatwo,lolimit,hilimit
      double precision xub,yub,zub,gamma,phic,xxx,yyy,zzz
      double precision bij_equil,bjk_equil,costheta
      double precision bij_energy,bjk_energy,aijk_energy,bik_distance
      double precision rik,riksq,onepi,dvalue,twopi,beta
      double precision twobvibone,twobvibtwo,bond_sdev_multiplier
      double precision benda_mult,bendb_mult,dihed_mult
      double precision viblo,vibhi,onethird,mean,std_dev
#if DEBUG_ATD_CBMC_SETUP
      double precision ideal_frac
      double precision phidisp,pdf_bias,lodim,hidim
#endif
c     --- local dp arrays
      double precision xbacord,ybacord,zbacord
      dimension xbacord(3),ybacord(3),zbacord(3)
      double precision xtcord,ytcord,ztcord
      dimension xtcord(4),ytcord(4),ztcord(4)
      double precision extrema
      dimension extrema(maxextrema,2)
      double precision prevweight,prevloc,prevenergy
      dimension prevweight(2),prevloc(2),prevenergy(2)
#if DEBUG_ATD_CBMC_SETUP
      double precision delta_dist
      dimension delta_dist(0:CBMCMAXTOR)
#endif

c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_constant_onethird(GLB_GET,onethird)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_beta(GLB_GET,beta)

c     --- set ibox to zero so we use the full coulombic energies instead
c     --- of one scaled by Ewald parameters
      ibox = 0

c     --- set up information needed to generation trials in a cbmc move
      call twh_cbmc_bondgen_code(GLB_GET,cbmc_bondgen_code)
      call twh_max_bond_length(GLB_GET,rmax)

      write(6,'(a35)') 'Determining CBMC bond distributions'
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            do ivib = 1,twh_get_invib(imolty,iunit)
               call twh_itvib(GLB_GET,imolty,iunit,ivib,itype)
               call twh_vibstyle(GLB_GET,itype,style)
c              --- check the style to see if is one of the special cases
               if ( style .eq. 1 ) then
c                 --- fixed bond length, everything else infinite
                  bonddist = DIST_DELTA
c                 --- only one value, normalization is 1.0 for both
                  truenorm = 1.0d0
                  arbnorm = 1.0d0
               elseif ( style .eq. 7 .or. style .eq. 8 
     &                 .or. cbmc_bondgen_code .eq. CBMC_IDEAL) then
c                 --- either no interaction, or just nonbonded
                  bonddist = DIST_R_SQ_IDEAL
c                 --- store the maximum bond length
                  call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib,3
     &                 ,rmax)
c                 --- normalization is ideal based only on max: 3/rmax^3
                  truenorm = 3.0d0/rmax**3
                  arbnorm = 3.0d0/rmax**3
               elseif ( style .eq. 10 ) then
c                 --- infinite square well with defined upper and lower
c                 --- bounds defined by the potential
                  bonddist = DIST_R_SQ_WITH_BOUNDS
                  call twh_vibcoeff(GLB_GET,itype,0,length)
                  call twh_vibcoeff(GLB_GET,itype,1,lolimit)
                  call twh_vibcoeff(GLB_GET,itype,2,hilimit)
c                 --- transform for use generating potentials
                  hilimit = hilimit**3 - lolimit**3
                  lolimit = lolimit**3
c                 --- store these values for use generating
                  call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                 ,3,lolimit)
                  call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                 ,4,hilimit)
c                 --- renorm to the full continuous distribution as
c                 --- the "penalty" for a small bond range
                  truenorm = 3.0d0/rmax**3
                  arbnorm = 3.0d0/hilimit
               else
c                 --- not a special case based on the potential
                  if ( cbmc_bondgen_code 
     &                 .eq. CBMC_R_SQ_WITH_BOUNDS ) then
                     bonddist = DIST_R_SQ_WITH_BOUNDS
c                    --- reset vibrang values for use in geometry
                     call twh_vibcoeff(GLB_GET,itype,0,length)
                     call twh_vibrang(GLB_GET,1,viblo)
                     call twh_vibrang(GLB_GET,2,vibhi)
                     lolimit = viblo*length
                     hilimit = vibhi*length
c                    --- transform for use generating potentials
                     hilimit = hilimit**3 - lolimit**3
                     lolimit = lolimit**3
c                    --- store these values for use generating
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,3,lolimit)
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,4,hilimit)
c                    --- even though our generation is bound, we use
c                    --- the ideal distribution to normalize
                     truenorm = 3.0d0/rmax**3
                     arbnorm = 3.0d0/hilimit
                  elseif ( cbmc_bondgen_code 
     &                    .eq. CBMC_GLOBAL_GAUSSIAN ) then
                     bonddist = DIST_GAUSSIAN
c                    --- uses a gaussian with parameters based on the
c                    --- equilibrium bond length and the default
c                    --- standard deviation
                     call twh_vibcoeff(GLB_GET,itype,0,mean)
                     call twh_sdevvib(GLB_GET,std_dev)
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,3,mean)
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,4,std_dev)
c                    --- normalization is determined in the gaussian
c                    --- selection process, true is the usual continuous
                     truenorm = 3.0d0/rmax**3
                     arbnorm = 1.0d0
                  elseif ( cbmc_bondgen_code 
     &                    .eq. CBMC_AUTOFIT_GAUSSIAN ) then
                     bonddist = DIST_GAUSSIAN
c                    --- determine gaussian parameters
                     junit = twh_get_ijvib(imolty,iunit,ivib)
                     imax = 1000
                     sum_w = 0.0d0
                     sum_x = 0.0d0
                     sum_xx = 0.0d0
                     do itest = 0,imax
                        length = dble(itest)*rmax/dble(imax)
                        lovrlap = .false.
                        energy = twh_vbond(lfailure,itype,length,imolty
     &                       ,iunit,junit,lovrlap)
                        if ( lfailure ) return
                        if ( lovrlap ) then
                           weight = 0.0d0
                        else
                           weight = twh_expon(-beta*energy)
     &                          *length*length
                        endif
                        sum_w = sum_w + weight
                        sum_x = sum_x + weight*length
                        sum_xx = sum_xx + weight*length*length
                     enddo
c                    --- average both quantities
                     sum_x = sum_x / sum_w
                     sum_xx = sum_xx / sum_w
c                    --- assign the mean
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,3,sum_x)
c                    --- assign the standard deviation
                     call twh_bond_sdev_multiplier(GLB_GET
     &                    ,bond_sdev_multiplier)
                     dvalue = bond_sdev_multiplier
     &                    *dsqrt(sum_xx - sum_x*sum_x)
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,4,dvalue)
#if DEBUG_ATD_CBMC_SETUP
                     write(6,'(a7,1x,i5,1x,a6,1x,2(i5,1x),a5,1x,f10.5
     &,1x,a8,1x,f10.7)') 
     &                    'Moltyp:',imolty,'Units:',iunit,junit
     &                    ,'Mean:',sum_x,'Std Dev:',dvalue
#endif
c                    --- normalization is ideal based only on rmax
                     truenorm = 3.0d0/rmax**3
                     arbnorm = 1.0d0
                  else
                     write(6,*) 'ATD_CBMC_SETUP: unknown cbmc_bondgen'
                     write(6,*) cbmc_bondgen_code
                     lfailure = .true.
                     return
                  endif
               endif
c              --- store the bonddist
               call twh_cbbonddist(GLB_SET,itype,bonddist)
c              --- store the bondnorms
               call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &              ,1,truenorm)
               call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &              ,2,arbnorm)
            enddo
         enddo
      enddo

      call twh_cbmc_bendgen_code(GLB_GET,cbmc_bendgen_code)
      if ( cbmc_bendgen_code .eq. CBMC_AUTOFIT_GAUSSIAN .or.
     &     cbmc_bendgen_code .eq. CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN
     &     ) then
         write(6,'(a46)') 
     &        'Determining Autofit Gaussian bend A parameters'
         call twh_bend_a_sdev_multiplier(GLB_GET,benda_mult)
         call twh_bend_b_sdev_multiplier(GLB_GET,bendb_mult)
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do ibend = 1,twh_get_inben(imolty,iunit)
                  junit = twh_get_ijben2(imolty,iunit,ibend)
                  kunit = twh_get_ijben3(imolty,iunit,ibend)
c                 --- take the equilibrium bond lengths for the angle
c                 --- computation
                  do ivib = 1,twh_get_invib(imolty,iunit)
                     if ( twh_get_ijvib(imolty,iunit,ivib) .eq. junit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,iunit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,distone)
                     endif
                  enddo
                  do ivib = 1,twh_get_invib(imolty,junit)
                     if ( twh_get_ijvib(imolty,junit,ivib) .eq. kunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,junit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,disttwo)
                     endif
                  enddo
                  imax = 1000
                  sum_w = 0.0d0
                  sum_x = 0.0d0
                  sum_xx = 0.0d0
                  do itest = 0,imax
                     theta = dble(itest)*onepi/dble(imax)
                     angcos = dcos(theta)
                     lovrlap = .false.
                     energy = twh_vangle(lfailure,lovrlap,1,imolty,iunit
     &                    ,ibend,0,angcos,distone,disttwo,xbacord
     &                    ,ybacord,zbacord)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        weight = 0.0d0
                     else
                        weight = twh_expon(-beta*energy)*dsin(theta)
                     endif
                     sum_w = sum_w + weight
                     sum_x = sum_x + weight*theta
                     sum_xx = sum_xx + weight*theta*theta
                  enddo
c                 --- average both quantities
                  sum_x = sum_x / sum_w
                  sum_xx = sum_xx / sum_w
c                 --- increment count
                  call twh_cbmcbendcount(GLB_INCR,imolty,iunit,ibend
     &                 ,ibend,1)
                  ipeak = 1
c                 --- assign the mean
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,1,sum_x)
c                 --- assign the standard deviation
                  dvalue = benda_mult*dsqrt(sum_xx - sum_x*sum_x)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,2,dvalue)
c                 --- assign limits
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,3,0.0d0)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,4,onepi)
#if DEBUG_ATD_CBMC_SETUP
                  write(6,'(a7,1x,i5,1x,a6,1x,3(i5,1x),a5,1x,f10.5
     &                 ,1x,a8,1x,f10.7)') 
     &                 'Moltyp:',imolty,'Units:',iunit,junit,kunit
     &                 ,'Mean:'
     &                 ,sum_x*180.0d0/onepi,'Std Dev:'
     &                 ,dvalue*180.0d0/onepi
#endif
               enddo
            enddo
         enddo
         write(6,'(a46)') 
     &        'Determining Autofit Gaussian bend B parameters'
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do ibend = 1,twh_get_inben(imolty,iunit)
                  junit = twh_get_ijben2(imolty,iunit,ibend)
                  kunit = twh_get_ijben3(imolty,iunit,ibend)
c                 --- take the equilibrium bond lengths for the angle
c                 --- computation
                  do ivib = 1,twh_get_invib(imolty,junit)
                     if ( twh_get_ijvib(imolty,junit,ivib) .eq. kunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,junit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,distone)
                     endif
                  enddo
                  itype = twh_get_itben(imolty,iunit,ibend)
                  thetaone = twh_get_bencoeff(itype,0)
                  phione = 0.0d0
                  do jbend = 1,twh_get_inben(imolty,iunit)
                     junitb = twh_get_ijben2(imolty,iunit,jbend)
                     kunitb = twh_get_ijben3(imolty,iunit,jbend)
                     if ( ibend .ne. jbend .and.
     &                    junit .eq. junitb ) then
                        do ivib = 1,twh_get_invib(imolty,junit)
                           if ( twh_get_ijvib(imolty,junit,ivib) 
     &                          .eq. kunitb ) then
                              call twh_itvib(GLB_GET,imolty,junit,ivib
     &                             ,itype)
                              call twh_vibcoeff(GLB_GET,itype,0,disttwo)
                           endif
                        enddo
                        itype = twh_get_itben(imolty,iunit,jbend)
                        thetatwo = twh_get_bencoeff(itype,0)
c                       --- get the values for the combo-angle
                        icrossbend = 0
                        do kbend = 1,twh_get_inben(imolty,kunit)
                           if ( twh_get_ijben2(imolty,kunit,kbend) 
     &                          .eq. junit
     &                          .and. twh_get_ijben3(imolty,kunit,kbend)
     &                          .eq. kunitb ) then
                              icrossbend = kbend
                           endif
                        enddo
                        if ( icrossbend .eq. 0 ) then
                           write(6,*)
     &                          'ATD_CBMC_SETUP: problem finding bend'
                           lfailure = .true.
                           return
                        endif
c                       --- initialize extrema
                        nextrema = 0
c                       --- check to see if we are sampling a rigid bend
                        type = twh_get_itben(imolty,kunit,icrossbend)
                        if ( twh_get_anglestyle(type) .eq. 0 ) then
c                          --- rigid angle, handled in CBMC move
                           lrigid = .true.
                        else
                           lrigid = .false.
                        endif

                        if ( .not. lrigid ) then
                           imax = 1000
                           do itest = -1,imax
                              phitwo = dble(itest)*twopi/dble(imax)
     &                             -onepi
                              lovrlap = .false.
                              call twh_coneangle(thetaone,phione
     &                             ,thetatwo,phitwo,angcos)
                              energy = twh_vangle(lfailure,lovrlap,1
     &                             ,imolty,kunit,icrossbend,0,angcos
     &                             ,distone,disttwo,xbacord,ybacord
     &                             ,zbacord)
                              if ( lfailure ) return
                              if ( lovrlap ) then
                                 weight = 0.0d0
                              else
                                 weight = twh_expon(-beta*energy)
                              endif
                              if ( itest .gt. 0 ) then
c                                --- see if the previous location was an
c                                --- extrema
                                 lfound = .false.
                                 if ( (prevweight(1) .gt. weight .and.
     &                                prevweight(1) .gt. prevweight(2))
     &                                .or.
     &                                (prevweight(1) .lt. weight .and.
     &                                prevweight(1) .lt. prevweight(2))
     &                                ) then
                                    lfound = .true.
                                 elseif ( prevweight(1) .eq. 0.0d0 .and.
     &                                   prevenergy(2) .lt.prevenergy(1)
     &                                   .and. energy .lt. prevenergy(1)
     &                                   ) then
                                    lfound = .true.
                                 endif
                                 if ( lfound ) then
                                    nextrema = nextrema + 1
                                    if ( nextrema .gt. maxextrema ) then
                                       write(6,*)
     &                                      'ATD_CBMC_SETUP:maxextrema'
     &                                      ,' exceeded'
                                       lfailure = .true.
                                       return
                                    endif
                                    extrema(nextrema,1) = prevloc(1)
                                    extrema(nextrema,2) = prevweight(1)
                                 endif
                              endif
c                             --- shift previous weights and locations
                              prevweight(2) = prevweight(1)
                              prevweight(1) = weight
                              prevloc(2) = prevloc(1)
                              prevloc(1) = phitwo
c                             --- shift previous energies
                              prevenergy(2) = prevenergy(1)
                              prevenergy(1) = energy
                           enddo
                           nmaxima = nextrema/2
                           if ( nmaxima .eq. 0 ) then
                              write(6,*)
     &                             'ATD_CBMC_SETUP: insufficient '
     &                             ,'extrema found for bends'
                              write(6,*) 'imolty,iunit,ibend',imolty
     &                             ,iunit,ibend
                              write(6,*) 'nextrema:',nextrema
                              do iextrema = 1,nextrema
                                 write(6,*) 'location, weight'
     &                                ,extrema(iextrema,1)
     &                                ,extrema(iextrema,2)
                              enddo
                              lfailure = .true.
                              return
                           endif
                           if ( nmaxima .gt. CBMC_MAXBENDPEAKS )then
                              write(6,*) 'ATD_CBMC_SETUP: nmaxima >'
     &                             ,' CBMC_MAXBENDPEAKS'
                              write(6,*) 'nmaxima',nmaxima
                              write(6,*) 'CBMC_MAXBENDPEAKS'
     &                             ,CBMC_MAXBENDPEAKS
                              lfailure = .true.
                              return
                           endif
                           iextrema = 0
                           do ipeak = 1,nmaxima
                              iextrema = iextrema + 1
                              lfound = .false.
                              do while ( .not. lfound ) 
                                 if ( iextrema .gt. 1 ) then
c                                   --- compare with previous extrema to
c                                   --- see if this is a maxima
                                    if ( extrema(iextrema,2) .gt. 
     &                                   extrema(iextrema-1,2) ) then
                                       lfound = .true.
                                    endif
                                 else
c                                   --- compare with next extrema to see
c                                   --- if this is a maxima
                                    if ( extrema(iextrema,2) .gt.
     &                                   extrema(iextrema+1,2) ) then
                                       lfound = .true.
                                    endif
                                 endif
                                 if ( .not. lfound ) then
                                    iextrema = iextrema + 1
                                 endif
                              enddo
c                             --- get the lo and hi bounding values
                              if ( iextrema .eq. 1 ) then
                                 lolimit = extrema(nextrema,1) - twopi
                                 hilimit = extrema(iextrema+1,1)
                              elseif ( iextrema .eq. nextrema ) then
                                 lolimit = extrema(nextrema-1,1)
                                 hilimit = extrema(1,1) + twopi
                              else
                                 lolimit = extrema(iextrema-1,1)
                                 hilimit = extrema(iextrema+1,1)
                              endif
                              sum_w = 0.0d0
                              sum_x = 0.0d0
                              sum_xx = 0.0d0
                              do itest = 0,imax
                                 phitwo = lolimit + dble(itest)
     &                                *(hilimit-lolimit)/dble(imax)
                                 lovrlap = .false.
                                 call twh_coneangle(thetaone,phione
     &                                ,thetatwo,phitwo,angcos)
                                 energy = twh_vangle(lfailure,lovrlap,1
     &                                ,imolty,kunit,icrossbend,0,angcos
     &                                ,distone,disttwo,xbacord,ybacord
     &                                ,zbacord)
                                 if ( lfailure ) return
                                 if ( lovrlap ) then
                                    weight = 0.0d0
                                 else
                                    weight = twh_expon(-beta*energy)
                                 endif
                                 sum_w = sum_w + weight
                                 sum_x = sum_x + weight*phitwo
                                 sum_xx = sum_xx + weight*phitwo*phitwo
                              enddo
c                             --- average both quantities
                              sum_x = sum_x / sum_w
                              sum_xx = sum_xx / sum_w
c                             --- increment count
                              call twh_cbmcbendcount(GLB_INCR,imolty
     &                             ,iunit,ibend,jbend,1)
c                             --- assign the mean
                              call twh_cbmcbendinfo(GLB_SET,imolty
     &                             ,iunit,ibend,jbend,ipeak,1,sum_x)
c                             --- assign the standard deviation
                              dvalue = bendb_mult
     &                             *dsqrt(sum_xx - sum_x*sum_x)
                              call twh_cbmcbendinfo(GLB_SET,imolty
     &                             ,iunit,ibend,jbend,ipeak,2,dvalue)
c                             --- assign the limits
                              call twh_cbmcbendinfo(GLB_SET,imolty
     &                             ,iunit,ibend,jbend,ipeak,3,lolimit)
                              call twh_cbmcbendinfo(GLB_SET,imolty
     &                             ,iunit,ibend,jbend,ipeak,4,hilimit)
#if DEBUG_ATD_CBMC_SETUP
                           call twh_cbmcbendcount(GLB_GET,imolty,iunit
     &                             ,ibend,jbend,count)
                           write(6,'(a7,1x,i5,1x,a6,1x,4(i5,1x),a6,1x
     &                          ,i1,1x,a5,1x,f7.2,1x,a8,1x,f7.2
     &                          ,1x,a10,1x,f7.2,1x,a10,1x,f7.2)') 
     &                          'Moltyp:',imolty,'Units:',iunit,junit
     &                          ,kunit,kunitb,'Count:',count
     &                          ,'Mean:',sum_x*180.0d0/onepi
     &                          ,'Std Dev:',dvalue*180.0d0/onepi
     &                          ,'Low Bound:',lolimit*180.0d0/onepi
     &                          ,'High Bound:',hilimit*180.0d0/onepi
#endif
                           enddo
                        endif
                     endif
                  enddo
               enddo
            enddo
         enddo
      endif

      call twh_cbmc_dihedgen_code(GLB_GET,cbmc_dihedgen_code)
      if ( cbmc_dihedgen_code .eq. CBMC_AUTOFIT_GAUSSIAN .or.
     &     cbmc_dihedgen_code .eq. CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN
     &     ) then
c        --- get the peak weighting style
         call twh_dihedral_peak_weight_style(GLB_GET,peakstyle)
         call twh_dihedral_sdev_multiplier(GLB_GET,dihed_mult)
         write(6,'(a48)') 
     &        'Determining Autofit Gaussian dihedral parameters'
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do itor = 1,twh_get_intor(imolty,iunit)
                  junit = twh_get_ijtor2(imolty,iunit,itor)
                  kunit = twh_get_ijtor3(imolty,iunit,itor)
                  lunit = twh_get_ijtor4(imolty,iunit,itor)
c                 --- set the iunit position
                  xtcord(1) = 0.0d0
                  ytcord(1) = 0.0d0
                  ztcord(1) = 0.0d0
c                 --- use the iunit-junit bond length to get the
c                 --- junit coordinates
                  do ivib = 1,twh_get_invib(imolty,iunit)
                     if ( twh_get_ijvib(imolty,iunit,ivib) .eq. junit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,iunit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  xtcord(2) = length
                  ytcord(2) = 0.0d0
                  ztcord(2) = 0.0d0
c                 --- use equilibrium junit-kunit length and the
c                 --- equilibrium iunit-junit-kunit angle to get
c                 --- the kunit coordinates
                  do ivib = 1,twh_get_invib(imolty,junit)
                     if ( twh_get_ijvib(imolty,junit,ivib) .eq. kunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,junit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  do ibend = 1,twh_get_inben(imolty,iunit)
                     if ( twh_get_ijben2(imolty,iunit,ibend) .eq. junit
     &                    .and.
     &                    twh_get_ijben3(imolty,iunit,ibend) .eq. kunit 
     &                    ) then
                        itype = twh_get_itben(imolty,iunit,ibend)
                        theta = twh_get_bencoeff(itype,0)
                     endif
                  enddo
c                 --- setup the cone
                  xxx = 1.0d0
                  yyy = 0.0d0
                  zzz = 0.0d0
                  gamma = 0.0d0
                  call twh_cone(lfailure,1,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
c                 --- set the angle
                  call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
                  xtcord(3) = xtcord(2) + length*xub
                  ytcord(3) = ytcord(2) + length*yub
                  ztcord(3) = ztcord(2) + length*zub
c                 --- reset the cone to get ready for the next atom
                  xxx = xub
                  yyy = yub
                  zzz = zub
                  call twh_cone(lfailure,1,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
                  do ivib = 1,twh_get_invib(imolty,kunit)
                     if ( twh_get_ijvib(imolty,kunit,ivib) .eq. lunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,kunit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  do ibend = 1,twh_get_inben(imolty,junit)
                     if ( twh_get_ijben2(imolty,junit,ibend) .eq. kunit
     &                    .and.
     &                    twh_get_ijben3(imolty,junit,ibend) .eq. lunit
     &                    ) then
                        itype = twh_get_itben(imolty,junit,ibend)
                        theta = twh_get_bencoeff(itype,0)
                     endif
                  enddo
c                 --- initialize extrema
                  nextrema = 0
                  imax = 1000
c                 --- initialize prevweights
                  prevweight(1) = 0.0d0
                  prevweight(2) = 0.0d0
                  do itest = -1,imax
                     gamma = dble(itest)*twopi/dble(imax) - onepi
                     call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                    ,xub,yub,zub)
                     if ( lfailure ) return
                     xtcord(4) = xtcord(3) + length*xub
                     ytcord(4) = ytcord(3) + length*yub
                     ztcord(4) = ztcord(3) + length*zub
                     lovrlap = .false.
                     energy = twh_vtorsion(lfailure,ibox,1,imolty,iunit
     &                    ,itor,xtcord,ytcord,ztcord,phic,lovrlap)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        weight = 0.0d0
                     else
                        weight = twh_expon(-beta * energy)
                     endif
                     if ( itest .gt. 0 ) then
c                       --- see if the previous location was an extrema
                        if ( (prevweight(1) .gt. weight .and.
     &                       prevweight(1) .gt. prevweight(2))
     &                       .or.
     &                       (prevweight(1) .lt. weight .and.
     &                       prevweight(1) .lt. prevweight(2))
     &                       ) then
                           nextrema = nextrema + 1
                           if ( nextrema .gt. maxextrema ) then
                              write(6,*) 'ATD_CBMC_SETUP:maxextrema'
     &                             ,' exceeded'
                              lfailure = .true.
                              return
                           endif
                           extrema(nextrema,1) = prevloc(1)
                           extrema(nextrema,2) = prevweight(1)
#if DEBUG_ATD_CBMC_SETUP
                           write(6,*) 'extrema count,loc,weight:'
     &                          ,nextrema,prevloc(1),prevweight(1)
#endif
                        endif
                     endif
c                    --- shift previous weights and locations
                     prevweight(2) = prevweight(1)
                     prevweight(1) = weight
                     prevloc(2) = prevloc(1)
                     prevloc(1) = gamma
                  enddo
                  nmaxima = nextrema/2
#if DEBUG_ATD_CBMC_SETUP
                  if ( nmaxima .eq. 0 ) then
                     call twh_cbmctorcount(GLB_GET,imolty,iunit,itor
     &                    ,ivalue)
                     write(6,'(a7,1x,i5,1x,a6,1x,4(i5,1x),a6,1x
     &                    ,i1,1x,a5,1x,f7.2,1x,a8,1x,f7.2
     &                    ,1x,a10,1x,f7.2,1x,a10,1x,f7.2)') 
     &                    'Moltyp:',imolty,'Units:',iunit,junit
     &                    ,kunit,lunit,'Count:',ivalue
                  endif
                  write(6,*) 'nmaxima:',nmaxima
                  write(6,*) 'nextrema:',nextrema
#endif
                  if ( nmaxima .gt. CBMC_MAXTORPEAKS )then
                     write(6,*) 'ATD_CBMC_SETUP: nmaxima exceeds'
     &                    ,' CBMC_MAXTORPEAKS'
                     write(6,*) 'nmaxima',nmaxima
                     write(6,*) 'CBMC_MAXTORPEAKS',CBMC_MAXTORPEAKS
                     lfailure = .true.
                     return
                  endif
                  iextrema = 0
                  do ipeak = 1,nmaxima
                     iextrema = iextrema + 1
                     lfound = .false.
                     do while ( .not. lfound ) 
                        if ( iextrema .gt. 1 ) then
c                          --- compare with previous extrema to
c                          --- see if this is a maxima
                           if ( extrema(iextrema,2) .gt. 
     &                          extrema(iextrema-1,2) ) then
                              lfound = .true.
                           endif
                        elseif ( iextrema .lt. nextrema ) then
c                          --- compare with next extrema to see if
c                          --- this is a maxima
                           if ( extrema(iextrema,2) .gt.
     &                          extrema(iextrema+1,2) ) then
                              lfound = .true.
                           endif
                        endif
                        if ( .not. lfound ) then
                           iextrema = iextrema + 1
                        endif
                     enddo
c                    --- get the lo and hi bounding values
                     if ( iextrema .eq. 1 ) then
                        lolimit = extrema(nextrema,1) - twopi
                        hilimit = extrema(iextrema+1,1)
                     elseif ( iextrema .eq. nextrema ) then
                        lolimit = extrema(nextrema-1,1)
                        hilimit = extrema(1,1) + twopi
                     else
                        lolimit = extrema(iextrema-1,1)
                        hilimit = extrema(iextrema+1,1)
                     endif
                     sum_w = 0.0d0
                     sum_x = 0.0d0
                     sum_xx = 0.0d0
                     do itest = 0,imax
                        gamma = lolimit + dble(itest)
     &                       *(hilimit-lolimit)/dble(imax)
                        call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                       ,xub,yub,zub)
                        if ( lfailure ) return
                        xtcord(4) = xtcord(3) + length*xub
                        ytcord(4) = ytcord(3) + length*yub
                        ztcord(4) = ztcord(3) + length*zub
                        lovrlap = .false.
                        energy = twh_vtorsion(lfailure,ibox,1,imolty
     &                       ,iunit,itor,xtcord,ytcord,ztcord,phic
     &                       ,lovrlap)
                        if ( lfailure ) return
                        if ( lovrlap ) then
                           weight = 0.0d0
                        else
                           weight = twh_expon(-beta*energy)
                        endif
                        sum_w = sum_w + weight
                        sum_x = sum_x + weight*gamma
                        sum_xx = sum_xx + weight*gamma*gamma
                     enddo
c                    --- average both quantities
                     sum_x = sum_x / sum_w
                     sum_xx = sum_xx / sum_w
c                    --- increment count
                     call twh_cbmctorcount(GLB_INCR,imolty,iunit,itor,1)
c                    --- assign the mean
                     call twh_cbmctorinfo(GLB_SET,imolty,iunit,itor
     &                    ,ipeak,1,sum_x)
c                    --- assign the standard deviation
                     dvalue = dihed_mult*dsqrt(sum_xx - sum_x*sum_x)
                     call twh_cbmctorinfo(GLB_SET,imolty,iunit,itor
     &                    ,ipeak,2,dvalue)
c                     --- assign the limits
                     call twh_cbmctorinfo(GLB_SET,imolty,iunit,itor
     &                    ,ipeak,3,lolimit)
                     call twh_cbmctorinfo(GLB_SET,imolty,iunit,itor
     &                    ,ipeak,4,hilimit)
c                    --- assign the raw weight
                     dvalue = sum_w / (hilimit - lolimit)
                     call twh_dihedral_peak_weights(GLB_SET,imolty
     &                    ,iunit,itor,ipeak,dvalue)

#if DEBUG_ATD_CBMC_SETUP
      call twh_cbmctorcount(GLB_GET,imolty,iunit,itor,ivalue)
      write(6,'(a7,1x,i5,1x,a6,1x,4(i5,1x),a6,1x
     &     ,i1,1x,a5,1x,f7.2,1x,a8,1x,f7.2
     &     ,1x,a10,1x,f7.2,1x,a10,1x,f7.2)') 
     &     'Moltyp:',imolty,'Units:',iunit,junit
     &     ,kunit,lunit,'Count:',ivalue
     &     ,'Mean:',sum_x*180.0d0/onepi
     &     ,'Std Dev:',dvalue*180.0d0/onepi
     &     ,'Low Bound:',lolimit*180.0d0/onepi
     &     ,'High Bound:',hilimit*180.0d0/onepi
#endif
                  enddo
c                 --- reset the peak weights
                  call twh_cbmctorcount(GLB_GET,imolty,iunit,itor
     &                 ,npeak)
                  if ( peakstyle .eq. 'uniform' ) then
c                    --- weight all peaks evenly
                     do ipeak = 1,npeak
                        dvalue = 1.0d0/dble(npeak)
                        call twh_dihedral_peak_weights(GLB_SET,imolty
     &                       ,iunit,itor,ipeak,dvalue)
                     enddo
                  elseif ( peakstyle .eq. 'isolated Boltzmann' ) then
c                    --- compute the total weight
                     sum_w = 0.0d0
                     do ipeak = 1,npeak
                        call twh_dihedral_peak_weights(GLB_GET,imolty
     &                       ,iunit,itor,ipeak,dvalue)
                        sum_w = sum_w + dvalue
                     enddo
c                    --- normalize
                     do ipeak = 1,npeak
                        call twh_dihedral_peak_weights(GLB_GET,imolty
     &                       ,iunit,itor,ipeak,dvalue)
                        dvalue = dvalue / sum_w
                        call twh_dihedral_peak_weights(GLB_SET,imolty
     &                       ,iunit,itor,ipeak,dvalue)
                     enddo
                  else
                     write(6,*) 'ATD_CBMC_SETUP unknown '
     &                    ,'dihedral_peak_weight_style'
                     write(6,*) peakstyle
                     lfailure = .true.
                     return
                  endif
#if DEBUG_ATD_CBMC_SETUP
      call twh_openfile(lfailure,99,'towhee_dist','DEFAULT','UNKNOWN'
     &                 ,'DEFAULT')
      if ( lfailure ) return
      do ipeak = 1,npeak
         call twh_cbmctorinfo(GLB_GET,imolty,iunit,itor,ipeak,3,lodim)
         call twh_cbmctorinfo(GLB_GET,imolty,iunit,itor,ipeak,4,hidim)
         write(6,*) 'set ipeak,lodim,hidim',ipeak,lodim,hidim
         do imgm = 0,100
            phidisp = (dble(imgm)/dble(100))
     &           *(hidim-lodim) + lodim
            write(6,*) 'phidisp:',phidisp
            call twh_dihedral_distribution(lfailure,.true.,.false.
     &           ,style_code
     &           ,npeak,delta_dist,lolimit,hilimit,mean,std_dev
     &           ,phidisp,dvalue,.true.,imolty,iunit,itor)
            if ( lfailure ) return
            write(6,*) 'lolimit,phidisp,hilimit',lolimit,phidisp,hilimit
            call twh_dihedral_ideal_fraction(GLB_GET,ideal_frac)
            pdf_bias = ideal_frac/twopi + (1.0d0-ideal_frac)
     &           *twh_gaussprob(phidisp,mean,std_dev
     &           ,lolimit,hilimit)*dvalue
            if ( phidisp .lt. 0.0d0 ) phidisp = phidisp + twopi
            write(99,*) phidisp*180.0d0/onepi
     &           ,pdf_bias*onepi/180.0d0
         enddo
      enddo
      close(99)
#endif
               enddo
            enddo
         enddo
      endif

c     --- two bond (and three bond that used two bond) fixed endpoint 
c     --- biasing setup
      call twh_two_bond_fe_bias_style(GLB_GET,twobbias_style)
      call twh_two_bond_bias_initial_value(GLB_GET
     &     ,twobbias_init)
      call twh_three_bond_fe_bias_style(GLB_GET,threebbias_style)
      if ( twobbias_style .eq. 'autofit gaussian'
     &     .or. threebbias_style .eq.
     &     'autofit gaussian using max and min 2-4 distance'
     &     .or. twobbias_init .eq. 'autofit gaussian' ) then
         write(6,'(a52,1x,a10)') 
     &        'Determining Autofit Gaussian two bond fixed endpoint'
     &        ,'parameters'
c        --- get the vibrange values
         call twh_two_bond_bias_vibrange(GLB_GET,1,twobvibone)
         call twh_two_bond_bias_vibrange(GLB_GET,2,twobvibtwo)
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do ibend = 1,twh_get_inben(imolty,iunit)
                  junit = twh_get_ijben2(imolty,iunit,ibend)
                  kunit = twh_get_ijben3(imolty,iunit,ibend)
c                 --- take the equilibrium bond lengths for the angle
c                 --- computation
                  ivib = 0
                  lfound = .false.
                  bij_equil = 1.0d0
                  do while ( ivib .lt. twh_get_invib(imolty,iunit) .and.
     &                 .not. lfound )
                     ivib = ivib + 1
                     if ( twh_get_ijvib(imolty,iunit,ivib) .eq. junit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,iunit,ivib
     &                       ,bij_type)
                        call twh_vibcoeff(GLB_GET,itype,0,bij_equil)
                        lfound = .true.
                     endif
                  enddo
                  if ( .not. lfound ) then
                     write(6,*) 'ATD_CBMC_SETUP: bij vib not found'
                     lfailure = .true.
                     return
                  endif
                  ivib = 0
                  lfound = .false.
                  bjk_equil = 1.0d0
                  do while ( ivib .lt. twh_get_invib(imolty,junit) .and.
     &                 .not. lfound )
                     ivib = ivib + 1
                     if ( twh_get_ijvib(imolty,junit,ivib) .eq. kunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,junit,ivib
     &                       ,bjk_type)
                        call twh_vibcoeff(GLB_GET,itype,0,bjk_equil)
                        lfound = .true.
                     endif
                  enddo
                  if ( .not. lfound ) then
                     write(6,*) 'ATD_CBMC_SETUP: bjk vib not found'
                     lfailure = .true.
                     return
                  endif
                  imax = 180
                  ibondmax = 100
                  sum_w = 0.0d0
                  sum_x = 0.0d0
                  sum_xx = 0.0d0
c                 --- scan the angle on 0,Pi
                  do itest = 0,imax
                     theta = dble(itest)*onepi/dble(imax)
                     angcos = dcos(theta)
                     lovrlap = .false.
c                    --- scan bij on two_bond_bias_range
                     do iscan_bij = 0,ibondmax
                        distone = bij_equil * ( twobvibone
     &                       + (dble(iscan_bij)/dble(ibondmax))
     &                       *(twobvibtwo - twobvibone) )
c                       --- compute the bij energy
                        bij_energy = twh_vbond(lfailure,bij_type,distone
     &                       ,imolty,iunit,junit,lovrlap)
                        if ( lfailure ) return
c                       --- scan bjk on two_bond_bias_range
                        do iscan_bjk = 0,ibondmax
                           disttwo = bjk_equil * ( twobvibone
     &                          + (dble(iscan_bjk)/dble(ibondmax))
     &                          *(twobvibtwo - twobvibone) )
c                          --- compute the bjk energy
                           bjk_energy = twh_vbond(lfailure,bjk_type
     &                          ,disttwo,imolty,junit,kunit,lovrlap)
                           if ( lfailure ) return
c                          --- compute the aijk energy
                           aijk_energy = twh_vangle(lfailure,lovrlap,1
     &                          ,imolty,iunit,ibend,0,angcos,distone
     &                          ,disttwo,xbacord,ybacord,zbacord)
                           if ( lfailure ) return
                           if ( lovrlap ) then
                              weight = 0.0d0
                           else
                              energy = bij_energy + bjk_energy 
     &                             + aijk_energy
                              weight = twh_expon(-beta*energy)
     &                             *dsin(theta)*distone*distone
     &                             *disttwo*disttwo
                           endif
c                          --- compute the ik distance
                           costheta = dcos(theta)
                           bik_distance = twh_lawofcosine(distone
     &                          ,disttwo,costheta)
                           sum_w = sum_w + weight
                           sum_x = sum_x + weight*bik_distance
                           sum_xx = sum_xx + weight*bik_distance**2
                        enddo
                     enddo
                  enddo
c                 --- average both quantities
                  sum_x = sum_x / sum_w
                  sum_xx = sum_xx / sum_w
c                 --- assign the mean
                  call twh_two_bond_bias_info(GLB_SET,imolty,iunit
     &                 ,ibend,1,sum_x)
c                 --- assign the standard deviation
                  dvalue = dsqrt(sum_xx - sum_x*sum_x)
                  call twh_two_bond_bias_info(GLB_SET,imolty,iunit
     &                 ,ibend,2,dvalue) 
#if DEBUG_ATD_CBMC_SETUP
                  write(6,'(a7,1x,i5,1x,a6,1x,3(i5,1x),a5,1x,f10.5
     &                 ,1x,a8,1x,f10.7)') 
     &                 'Moltyp:',imolty,'Units:',iunit,junit,kunit
     &                 ,'Mean:',sum_x
     &                 ,'Std Dev:',dvalue
#endif
               enddo
            enddo
         enddo
      endif

      call twh_three_bond_bias_initial_value(GLB_GET,threebbias_init)
      if ( threebbias_init  .eq. 'autofit gaussian' ) then
         write(6,'(a50)') 
     &        'Determining Autofit Gaussian Three Bond parameters'
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do itor = 1,twh_get_intor(imolty,iunit)
                  junit = twh_get_ijtor2(imolty,iunit,itor)
                  kunit = twh_get_ijtor3(imolty,iunit,itor)
                  lunit = twh_get_ijtor4(imolty,iunit,itor)
c                 --- set the iunit position
                  xtcord(1) = 0.0d0
                  ytcord(1) = 0.0d0
                  ztcord(1) = 0.0d0
c                 --- use the iunit-junit bond length to get the
c                 --- junit coordinates
                  do ivib = 1,twh_get_invib(imolty,iunit)
                     if ( twh_get_ijvib(imolty,iunit,ivib) .eq. junit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,iunit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  xtcord(2) = length
                  ytcord(2) = 0.0d0
                  ztcord(2) = 0.0d0
c                 --- use equilibrium junit-kunit length and the
c                 --- equilibrium iunit-junit-kunit angle to get
c                 --- the kunit coordinates
                  do ivib = 1,twh_get_invib(imolty,junit)
                     if ( twh_get_ijvib(imolty,junit,ivib) .eq. kunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,junit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  do ibend = 1,twh_get_inben(imolty,iunit)
                     if ( twh_get_ijben2(imolty,iunit,ibend) .eq. junit
     &                    .and.
     &                    twh_get_ijben3(imolty,iunit,ibend) .eq. kunit
     &                    ) then
                        itype = twh_get_itben(imolty,iunit,ibend)
                        theta = twh_get_bencoeff(itype,0)
                     endif
                  enddo
c                 --- setup the cone
                  xxx = 1.0d0
                  yyy = 0.0d0
                  zzz = 0.0d0
                  gamma = 0.0d0
                  call twh_cone(lfailure,1,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
c                 --- set the angle
                  call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
                  xtcord(3) = xtcord(2) + length*xub
                  ytcord(3) = ytcord(2) + length*yub
                  ztcord(3) = ztcord(2) + length*zub
c                 --- reset the cone to get ready for the next atom
                  xxx = xub
                  yyy = yub
                  zzz = zub
                  call twh_cone(lfailure,1,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
                  do ivib = 1,twh_get_invib(imolty,kunit)
                     if ( twh_get_ijvib(imolty,kunit,ivib) .eq. lunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,kunit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  do ibend = 1,twh_get_inben(imolty,junit)
                     if ( twh_get_ijben2(imolty,junit,ibend) .eq. kunit
     &                    .and.
     &                    twh_get_ijben3(imolty,junit,ibend) .eq. lunit
     &                    ) then
                        itype = twh_get_itben(imolty,junit,ibend)
                        theta = twh_get_bencoeff(itype,0)
                     endif
                  enddo
                  imax = 1000
                  sum_w = 0.0d0
                  sum_x = 0.0d0
                  sum_xx = 0.0d0
                  lolimit = -onepi
                  hilimit = onepi
                  do itest = 0,imax
                     gamma = lolimit + dble(itest)
     &                    *(hilimit-lolimit)/dble(imax)
                     call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                    ,xub,yub,zub)
                     if ( lfailure ) return
                     xtcord(4) = xtcord(3) + length*xub
                     ytcord(4) = ytcord(3) + length*yub
                     ztcord(4) = ztcord(3) + length*zub
c                    --- compute the 1-4 distance
                     riksq = (xtcord(4)-xtcord(1))**2 
     &                    + (ytcord(4)-ytcord(1))**2
     &                    + (ztcord(4)-ztcord(1))**2
                     rik = dsqrt(riksq)
                     lovrlap = .false.
                     energy = twh_vtorsion(lfailure,ibox,1,imolty,iunit
     &                    ,itor,xtcord,ytcord,ztcord,phic,lovrlap)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        weight = 0.0d0
                     else
                        weight = twh_expon(-beta*energy)
                     endif
                     sum_w = sum_w + weight
                     sum_x = sum_x + weight*rik
                     sum_xx = sum_xx + weight*riksq
                  enddo
c                 --- average both quantities
                  sum_x = sum_x / sum_w
                  sum_xx = sum_xx / sum_w
c                 --- assign the mean
                  call twh_three_bond_bias_info(GLB_SET,imolty,iunit
     &                 ,itor,1,sum_x)
c                 --- assign the standard deviation
                  dvalue = dsqrt(sum_xx - sum_x*sum_x)
                  call twh_three_bond_bias_info(GLB_SET,imolty,iunit
     &                 ,itor,2,dvalue)
#if DEBUG_ATD_CBMC_SETUP
                     write(6,'(a7,1x,i5,1x,a6,1x,4(i5,1x)
     &                    ,a5,1x,f7.2,1x,a8,1x,f7.2)') 
     &                    'Moltyp:',imolty,'Units:',iunit,junit
     &                    ,kunit,lunit,'Mean:',sum_x,'Std Dev:',dvalue
#endif
               enddo
            enddo
         enddo
      endif

      return
      end

