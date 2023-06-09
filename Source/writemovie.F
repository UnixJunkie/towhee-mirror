#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2020 Marcus G. Martin                       *
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
      subroutine twh_writemovie(lfailure,step)
c     ******************************************************************
c     * this routine outputs a single frame of information to the      *
c     * towhee_movie file for post-processing analysis                 *
c     *                                                                *
c     * rewritten 07-09-2003 by M.G. Martin from pupdate               *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer step
c     --- local variables
c     --- character string scalars
      character*30 filename
c     --- integer scalars
      integer ibox,imolty,idim,ichain,iunit,jdim,ifile,nstep,numboxes
      integer nframe,movieversion,nmolty,nunit,moviefreq,nboxi,ntype
      integer ivib,itor,nchain,ijtor2,ijtor3,ijtor4,invib,intorsion
c     --- integer arrays
      integer ijvibdata
      dimension ijvibdata(NNBOND)
      integer localncmt
      dimension localncmt(MAXBOX)
c     --- double precision scalars
      double precision xtemp,ytemp,ztemp,rcut,qqatom
c     --- double precision arrays
      double precision hmtemp
      dimension hmtemp(3)

c     --- retrieve constants
      call twh_nmolty(GLB_GET,nmolty)
      call twh_numboxes(GLB_GET,numboxes)
      call twh_rcut(GLB_GET,rcut)

c     --- set the movie file number
      ifile = FILE_MOVIE

      call twh_nchain(GLB_GET,nchain)
      if ( step .eq. 0 ) then
c        --- initializing movie file output
c        --- calculate number of frames in the movie file
         call twh_nstep(GLB_GET,nstep)
         call twh_moviefreq(GLB_GET,moviefreq)
         nframe = nstep / moviefreq
 
         if ( moviefreq .le. nstep ) then
c           --- open the movie file
            filename = 'towhee_movie'
            call twh_openfile(lfailure,ifile,filename,'DEFAULT'
     &           ,'UNKNOWN','DEFAULT')
            if ( lfailure ) return
c           --- set current movie version
            movieversion = 2
c           --- write out movie-header 
            write(ifile,*) movieversion
            write(ifile,*) nframe,nchain,nmolty,numboxes
     &           ,(rcut,ibox=1,numboxes)
            
            do imolty = 1,nmolty
               call twh_nunit(GLB_GET,imolty,nunit)
               write(ifile,*) nunit
c              --- output bond connectivity information
               do iunit = 1,nunit
                  call twh_invib(GLB_GET,imolty,iunit,invib)
                  do ivib = 1,invib
                     call twh_ijvib(GLB_GET,imolty,iunit,ivib
     &                    ,ijvibdata(ivib))
                  enddo
                  write(ifile,*) invib,(ijvibdata(ivib),ivib=1,invib)
               enddo
            
c              --- output torsional connectivity information
               do iunit = 1,nunit
                  call twh_intor(GLB_GET,imolty,iunit,intorsion)
                  write(ifile,*) intorsion
                  do itor = 1,intorsion
                     call twh_ijtor2(GLB_GET,imolty,iunit,itor,ijtor2)
                     call twh_ijtor3(GLB_GET,imolty,iunit,itor,ijtor3)
                     call twh_ijtor4(GLB_GET,imolty,iunit,itor,ijtor4)
                     write(ifile,*) ijtor2,ijtor3,ijtor4
                  enddo
               enddo
            enddo
         endif
      else
c        --- regular output of a movie frame
         write(ifile,*) step
         do ibox = 1,numboxes
            do imolty = 1,nmolty
               call twh_ncmt(GLB_GET,ibox,imolty,localncmt(imolty))
            enddo
            write(ifile,*) (localncmt(imolty),imolty=1,nmolty)
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_GET,BXD_REAL,ibox,idim,jdim
     &                 ,hmtemp(jdim))
               enddo
               write(ifile,*) (hmtemp(jdim),jdim=1,3)
            enddo
         enddo
         
         do ichain = 1, nchain
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
            call twh_nboxi(GLB_GET,ichain,nboxi)
            write(ifile,'(4(1x,i10))') ichain,imolty,nunit,nboxi
            do iunit = 1, nunit
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               call twh_ntype(GLB_GET,imolty,iunit,ntype)
               call twh_qqatom(GLB_GET,imolty,iunit,qqatom)
               write(ifile,'(4(1x,e15.6),i10)') xtemp,ytemp,ztemp
     &              ,qqatom,ntype
            enddo
         enddo
      endif

      return
      end
