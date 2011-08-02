#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2011 Marcus G. Martin                       *
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
      subroutine twh_readinput(lfinish)
c     ******************************************************************
c     * opens towhee_input and reads in the inputformat                *
c     * then it calls the appropriate input reading routine            *
c     *                                                                *
c     * rewritten 05-15-2002 by M.G. Martin                            *
c     * last modified 02-08-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "wrapperF.h"
#define FUNCTION_CHECK_LABEL
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_READ_LABELED_LOGICAL
#include "functions.h"

c     --- variables passed to/from the subroutine
      logical lfinish
c     --- local variables
      logical lallow,loutput
      integer*4 luxvalue,iseed
      character*15 inputformat,idname
      character*(MAXDIRLENGTH) io_dir

      idname = 'readinput'
      lfinish = .false.

c     --- retrieve constants
      call twh_loutput(GLB_GET,loutput)

c     --- open the input file for towhee
c     --- This is quite klunky right now: for parstyle = jobfarm,
c     --- io_directory is set; for other parstyles, 
c     --- warp_towhee_input_file may be set from the command line to
c     --- allow reading from different input files.  It
c     --- would be good to simplify this...
      call twh_openfile(lfinish,4,wrap_towhee_input_file,'DEFAULT'
     &     ,'OLD','DEFAULT')
      if ( lfinish ) return
      if ( loutput ) then
         write(6,'(1x,a,1x,a)') "Reading from towhee_input file:", 
     &        wrap_towhee_input_file
         call twh_io_directory(GLB_GET,io_dir)
         if (io_dir .eq. '') then
            write(6,'(1x,a)') "in directory: current directory"
         else 
            write(6,'(1x,a,1x,a)') "in directory:",io_dir
         endif
      endif

c     --- open towhee_altinp
      call twh_openfile(lfinish,15,'towhee_altinp','DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfinish ) return
c     --- read in inputformat
      call twh_read_labeled_string(4,'inputformat',loutput
     &     ,idname,inputformat)
c     --- optional read in the random number seed
      if (twh_check_label(4,'randomseed')) then
         call twh_get_randomseed(4, iseed, .false., idname)
         if ( iseed .lt. 0 ) then
            write(6,*) 'READINPUT: randomseed must be postitive'
            lfinish = .true.
            return
         endif
c        --- set the global randomseed
         call twh_randomseed(GLB_SET,iseed)
      endif

c     --- optional get the random number generator luxury level
      if ( twh_check_label(4,'random_luxlevel') ) then
         luxvalue = twh_read_labeled_integer(4,'random_luxlevel'
     &        ,.true.,idname)
c        --- error check
         if ( luxvalue .lt. 0 .or. luxvalue .gt. 8 ) then
            write(6,*) 'READINPUT: random_luxlevel must be between '
     &           ,'0 and 8 inclusive (levels 5-8 choose various RNGs)'
            lfinish = .true.
            return
         endif
         call twh_random_luxlevel(GLB_SET,luxvalue)
      endif

c      --- optional read random_allow_restart
      if ( twh_check_label(4,'random_allow_restart') ) then
         lallow = twh_read_labeled_logical(4,'random_allow_restart'
     &        ,.true.,idname)
         call twh_random_allow_restart(GLB_SET,lallow)
      else
c        --- use the default value
         call twh_random_allow_restart(GLB_GET,lallow)
      endif
c     --- see if we are forcing initialiation here
      if ( .not. lallow ) then
c        --- initialize random number generator now
         call twh_init_random_single()
      endif

c     --- head off to the various read routines
      if ( inputformat .eq. 'LAMMPS' ) then
         call twh_readlammps(lfinish)
         if ( lfinish ) return
         write(6,*) 'In order to run Towhee with the LAMMPS generated'
         write(6,*) 'input files use the following commands and then'
         write(6,*) 'run towhee again'
         write(6,*) 'cp towhee_altinp towhee_input'
         write(6,*) 'cat lammps_coords_* > towhee_coords'
         lfinish = .true.
      elseif ( inputformat .eq. 'Towhee' ) then
c        --- read in the rest of the file according to everything 
c        --- below
         call twh_readtowhee(lfinish)
      else
         write(6,*) 'READINPUT: invalid inputformat:',inputformat
         lfinish = .true.
      endif

c     --- close towhee_input and towhee_altinp
      close(4)
      close(15)

      return
      end



c     ******************************************************************
c     * Returns randomseed.  This value may come from one of two       *
c     * sources: from the towhee_input file or,                        *
c     * if lreadwrapper is set (in a call to towheemainloop) from      *
c     * variables set (by the C wrapper) in the wrapper common block.  *
c     *                                                                *
c     * Variables:                                                     *
c     * io_unit specifies the file number to read from randseed is set *
c     * by this subroutine.  This local variable will typically be     *
c     * assigned to the appropriate global variable by the caller.     *
c     * print_value specifies whether the label and value output to    *
c     * screen.                                                        *
c     * idname is the calling subroutine, printed only in case of error*
c     *                                                                *
c     * originally written prior to 2006                               *
c     * last modified 06-25-2008 by M.G. Martin                        *
c     ******************************************************************
      subroutine twh_get_randomseed(io_unit, randseed, print_value, 
     &   idname)
      implicit none
#include "wrapperF.h"
#define FUNCTION_CHECK_LABEL
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"

c     --- variables passed to/from the subroutine
      integer randseed, io_unit
      logical print_value
      character*(*) idname
c     --- local variables
      logical lreadwrapper

      call twh_lreadwrapper(GLB_GET,lreadwrapper)
      if (lreadwrapper) then
c     --- if randomseed comes from wrapper this parameter is allowed,
c     --- but optional, in towhee_input.  See if it exists,  read it
c     --- so that its skipped over.
         if (twh_check_label(io_unit,'randomseed')) then
            randseed = twh_read_labeled_integer(io_unit,'randomseed'
     &           ,.false.,idname)
         endif

c     --- actual value comes from wrapper
         randseed = wrap_randomseed
         if (print_value) then
            write(6,*) 'randomseed (from wrapper): ',randseed
         endif
      else
         randseed = twh_read_labeled_integer(io_unit,'randomseed', 
     &      print_value, idname)
      endif
      return
      end

      
      subroutine twh_init_random_single()
c     ******************************************************************
c     * initializes the random number generator using a single integer *
c     * seed                                                           *
c     *                                                                *
c     * originally written 11-13-2006 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#include "functions.h"
c     --- local variables
      logical loutput
      integer*4 luxvalue,iseed
      integer itest
      double precision rtest
      dimension rtest(5)

c     --- retrieve constants
      call twh_loutput(GLB_GET,loutput)
      call twh_random_luxlevel(GLB_GET,luxvalue)
      call twh_randomseed(GLB_GET,iseed)
      call rnginit(luxvalue,iseed)
      do itest=1,5
         rtest(itest) = twh_random()
      enddo
      if ( loutput ) then
c        --- print random numbers
         write(6,'(1x,a)') 'Testing random number generator'
         write(6,'(2x,5f10.6)') (rtest(itest),itest=1,5)
      endif

      return
      end
