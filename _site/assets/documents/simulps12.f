c23456789a123456789b123456789c123456789d123456789e123456789f123456789g12
c  Program SIMULPS12 (November 27, 1993)
c   This version inverts for Vp and Vp/Vs (previous versions inverted
c   for Vp and Vs).  The input data are P travel-times and S-P times
c   (previous versions used P travel-times and S travel-times.  See
c   Thurber (ref 1 below) for discussion of solution for Vp/Vs
c   using S-P data.
c   This version also allows the user to vary the weighting, in the 
c   hypocenter solution, of the S-P data relative to the P data. (see wtsp)
c   Vp and Vp/Vs are input and are used to compute Vs at each iteration.
c   The S velocities are on the same grid as the p velocities and are
c   stored in additional nz layer positions.  
c   Uses up to 10000 velocity nodes (and station parameters), but only invert for up to 700.
c   Has the option of outputing the raypaths.
c   Has Thurber's psuedo-bending ray-tracing from his SIMUL3M version
c   Allows less curvature for initial arcuate rays below "moho"
c   (see i3d input parameter).
c
c        Donna Eberhart-Phillips, Ph.D.
c        U. S. Geological Survey
c        Office of Earthquakes, Volcanoes & Engineering
c        525 S. Wilson Ave.
c        Pasadena, CA  91106
c        (818) 405-7240, usgs (818) 405-7823
c        fax: (818) 405-7827
c        email:  eberhart@bombay.gps.caltech.edu
c
c
c  If you are using this program please read and reference the following
c     chapters by myself & Cliff Thurber:
c
c  1)  Thurber, C. H., Local earthquake tomography: velocities and Vp/Vs - theory,
c            in Seismic Tomography: Theory and Practice, edited by
c            H. M. Iyer and K. Hirahara, 1993. 
c
c  2)  Eberhart-Phillips, D., Local earthquake tomography: earthquake source regions,
c            in Seismic Tomography: Theory and Practice, edited by
c            H. M. Iyer and K. Hirahara, 1993. 
c
c  3)  Thurber, C. H., Earthquake locations and three-dimensional crustal structure
c             in the Coyote Lake area, central California, J. Geophys. Res., 
c             v. 88, p. 8226-8236, 1983.
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c  USGS disclaimer statement:
c  ALTHOUGH THIS PROGRAM HAS BEEN USED BY THE USGS, NO WARRANTY,
c  EXPRESSED OR IMPLIED, IS MADE BY THE USGS OR THE UNITED
c  STATES GOVERNMENT AS TO THE ACCURACY AND FUNCTIONING OF THE
c  PROGRAM AND RELATED PROGRAM MATERIAL NOR SHALL THE FACT OF
c  DISTRIBUTION CONSTITUTE ANY SUCH WARRANTY, AND NO 
c  RESPONSIBILITY IS ASSUMED BY THE USGS IN CONNECTION THEREWITH.
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c 
c  Record of changes I have made in program (DMEP).
c
c  30-jun-94 Also save ndof1,wndof1 in OUTADJ
c  22-jun-94 Save var1,varw1 in OUTADJ instead of recomputing
c     in DECIDE
c  10-jun-94 Clarified ftest variables.
c 21-apr-94 Changed 'a' to 'ac' for hypocenter error ellipse
c   calculation in LOCEQK
c 15-mar-94 Added output file 36 that has summary output.
c 11-mar-94 Changes throughout program to have common blocks
c    in a separate file that is included: simulps_common.inc
c 7-mar-94 Change to TTMDER.  Improved check for writing to
c    boundary nodes.
c 27-nov-93  Several changes to take care of "uninitialized variables".
c    Most of these were in HYPOINVERSE write statements in subroutine OUT
c 25-nov-93  This version has been modified by J. Evans so that
c    it will compile better on sun machines.  (changes to do loops
c  4-oct-93 sdl was not being used (in s-p version).  Change to TTMDER so that sdl is
c      used.
c  12-jul-93 Minor change to OUTEND in format for writing nparvi to nodes output file
c  08-jul-93 Do not allow damping to become smaller (with idmp)
c     Also do not calculate new damp if zero nodes observed
c  07-jul-93 Change npar dimension to 10000
c  03-jul-93 Change to OUTEND for writing station output
c  23-jun-93 Minor change to OUTEND for writing number of stations
c  22-jun-93 Change number of stations allowed to 1800.
c  21-jun-93 Change to INPUT4 to note observed stations that are not in station list
c            if kout2 eq 0 or 1
c  16-jun-93 Change to BEND to avoid divide by zero for constant velocity (zero gradient)
c  03-jun-93 Change to FORWRD.  Corrected way synthetic data was computed for
c      S-P observations.
c  27-apr-93 Change to MEDDER.  Deleted old wr= ... dres that should have
c       been taken out when res3 was added to residual downweighting.
c  23-apr-93 Changes to VELADJ so that solution statistics will be correct
c      for the case where there are relatively few velocity nodes relative
c      to the number of station corrections.
c  19-mar-93 Minor change in INPUT4. 
c  23-feb-93 Change to subroutine OUTEND so that ttobs for S-P is correct
c      for file024.
c  10-feb-93 Change to subroutine OUTRES so that ttobs for S-P is correct
c  1-feb-93 Various minor changes to make more compatible with 
c      sun compiler
c      (Note compile on sun with f77 -lV77 simulps12.for -o simulps12)
c  7-jan-93 Also changed VELBKU for vp/vs
c  6-jan-93 Changes of OUTADJ, OUTEND so do not print out grids that
c      have dws=0.0
c      Changes to VELADJ for case of vp and sta corr only (no vp/vs).
c  5-jan-93 Changes to VELADJ to output sol norm and damping
c      for both vp and vp/vs.
c      Added "wtsp" to allow user to change relative weighting of
c      S-P observations in hypocentral solution.  Changes to WTHYP.
c  30-dec-92 Made changes so that vpvs array is used directly. (Cliff
c       had solved for vp/vs but then used the vp/vs perturbation in
c       veladj to perturb the associated vs element.)  Now vp/vs is 
c       input instead of vs, and vs is calculated from vp and vp/vs
c       in input3 and at end of veladj.
c  5-nov-92 Now have option to create synthetic travel-time data.
c           Use nitmax= -1, use file007 for input hypocenters and
c           stations.  Calculated travel-times will be output to
c           file028.  Made changes to Main and Outend.
c  3-nov-92 Made various changes suggested by Cliff Thurber in order to
c     use S-P times to invert for Vp/Vs instead of using S times to
c     invert for Vs.
c  12-mar-92 Changes to RAYWEB to have different curvature for rays below "moho"
c  24-feb-92 as noted by L. Hutchings, added nbls to variance calculation in RESCOV
c  24-feb-92 Added "bld" to input3.  Bld is factor for setting up velocity
c      interpolation arrays (1.0 or 0.1). Now velocity nodes can be defined to 1/10th km
c      if desired.  Changes to BLDMAP, INPUT3, INTMAP, OUTEND.
c 15-mar-91 In VELADJ, DECIDE, RESCOV, added khit=0 check to all hitct checks 
c    (and deleted hit=0 checks).  Parsep uses khit=0 to setup solution matrix so it
c    should also be used here.
c 22-feb-91 Changes from Egill Hauksson. Corrected variable name in 
c    zeroing out section of RESCOV.  Added hitct check to rhm summation
c    in VELADJ. Also minor printout change in INPUT3
c 9-oct-90 In DECIDE use rmsw instead of rms for comparison to rmstop
c 4-oct-90 In OUTEND use ssqrw instead of ssqr in calculating 'rat' for
c      use in stderr.
c  7-sep-90 Write out 1/covariance(diag) if kout2.eq.5 and ires.gt.0
c      to for016 and for045; changes for OUTEND, RESCOV.
c  3-aug-90 Made changes suggested by Goran Ekstrom for convergence 
c      check in LOCEQK.
c  12-jul-90 Put in a second linear residual weighting, 98% of downweighting
c      is done res1 to res2, last 2% of downweighting is done 
c      res2 to res3.  This is useful in case you happen to get a poor
c      hypocenter with mostly very high residuals.  If res3=res2,
c      100% of downweighting is done res1 to res2, as before.
c      Remove "inew" from MAIN since it's never used.
c      In MAIN, if nitmax.eq.0, skip parsep for blasts & shots also.
c  11-jul-90 Use hit array instead of khit to do cutoff in VELADJ,
c      also changing "nhitct" to "hitct".  A peripheral node can have
c      a high enough khit even though it has no rays near it.
c  27-jun-90 Change input4 so it can also read Alberto's format with
c      millisec. travel-times
c  15-jun-90 Output data variance as rms**2 also.  Correct weighted 
c      variance which should use wnobt in calculating wnodf.  Changes
c      to DECIDE.
c  11-apr-90 Allow ires=3 to only compute resolution on 1st iteration.
c      Saves cpu and is appropriate when idmp=1(damping increases,
c      resolution decreases on later iterations.)
c      Changes to MAIN.
c  6-oct-89 Add "snrmct" to stop iterating if solution norm is
c      less than the specified cutoff value.
c      Changes to MAIN,INPUT1,VELADJ.
c  4-oct-89 Can now invert for station corrections at selected stations
c      by fixing delays (to 0) at other stations.  Changes to common
c      OBSERVE, routines INPUT2, INPUT3, TTMDER, VELADJ.
c  6-jun-89 Fixed hypocentral error (diagonal elements were incorrectly
c      coded as noted by Goran Ekstrom), and error ellipse (for 
c      coordinate rotation) in LOCEQK
c     Have LOCEQK calculate error ellipses and print out
c      HYPOINVERSE-style summary cards, if kout=5.  OUT now called
c      from LOCEQK.
c     In LOCEQK, set ster to 0 before calc. so correct for blasts
c     In SETUP,TTMDER calculation of number of raypath segments, changed
c      from IFIX to NINT, should make it more accurate for short paths
c      where few segments are used
c  2-mar-89 Minor changes to RESCOV so more appropriate for ires=1
c     Changes to TTMDER so that it checks for writing outside defined
c     dtm-array.  This may have caused errors before (numerous 'poorly
c     constrained depth') if deepest grid was not far enough away.
c  22-sep-87 Made a variety of changes to make the program work properly
c     when station delays are included in the inversion (also for the
c     case when all velocity gridpoints are fixed).  Removed line in 
c     Parsep that was wrong for station corrections in inversion
c     Set index arrays properly for station delays, in Input3.
c     Skip sections of ttmder, veladj, velbku, outadj, outend if all velocity nodes are fixed
c     Don't calculate station delay partial derivatives when Ttmder 
c     is called from Loceqk.
c     Fixed array sizes in Main, Strt.
c     Have 3 damping parameters (3rd for station delays), changes to Input1
c  20-feb-87 Changed program so that when earthquake locations are fixed on the
c     first iterations, they are not weighted like shots or blasts.
c  12-jul-86  Fixed dr5 calculation in Path.
c  10-jul-86  Fixed initial value of xmov, pl(no) in LOCEQK.
c  08-jul-86  Also tally station observations by P and S.
c     In Loceqk have last 3d raypath remembered and used as the initial
c     path in Minima, if small change in hypocenter location on last
c     iteration.
c  05-jul-86  Fixed jndex array (which is used for p-or-s decisions).  
c     It had not been updated after poorly sampled nodes were deleted.
c  30-jun-86  Fixed a counter error in MEDDER.  Added write
c     statements to parsep to where resolution problem is
c     arising.
c
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
c
c  (In addition to these notes, see also the handout written by
c   Cliff)
c
c  Written by Cliff Thurber  as part of his PhD thesis.
c  Modified by W. Prothero to include station delays.
c  Parts of this were programmed by Steve Taylor of LLL.
c  Obtained from Prothero&Thurber in 1983 and subsequently
c  modified by Donna Eberhart-Phillips, U.S. Geological Survey.
c
c
c  Input data file list:
c     File 01 - Control Parameters
c        line 1- (free format)
c        neqs - number of earthquakes
c        nsht  - number of shots
c        nbls  - number of blasts with known location, but unknown
c                origin time
c        wtsht - weighting of shots relative to quake weighting
c        kout - output control parameter
c               value    files created
c                 0      16,36
c                 1      16,13,36
c                 2      16,13,22,23,24,36
c                 3      16,13,22,23,24,34,36
c                 4      16,13,22,23,24,25,34,36
c                 5      16,12,13,22,23,24,25,34,36
c        kout2 - Printout control parameter
c                0 = full printout including station residuals and
c                    location steps
c                1 = printout station residuals
c                2 = printout location steps
c                3 = don't printout location steps or station residuals
c                4 = as 3 above, also don't printout stations in input2
c                5 = as 4. Also printout 1/(diag. covariance) to for016 
c                    and for045, if ires.gt.0.
c        kout3 - Yet another output control parameter
c                0 = Don't output raypath points or tt differences
c                1 = Output raypath points to file 15, for all raypaths
c                    Print out travel-time differences between ART and
c                    psuedo-bending to File 19.
c                    This is useful for making plots of raypaths.
c                    (For instance, if you want to test a range of psuedo-
c                    bending parameters (xfac,tlim,nitpb) )
c                    **NOTE that this option should not be used regularly,
c                    ** but only to check a few selected events, since
c                    ** it creates a lot of output.
c
c        line 2 - (free format)
c        nitloc - max of iterations for hypocenter location.
c        wtsp   - for hypocenter solution, weight of S-P residual 
c                 relative to P residual (ie:wtsp=1.0 gives equal wt,
c                 wtsp<1 downweights S-P)
c        eigtol - SVD cutoff in hypocentral adjustments. If smallest
c                 eigenvalue in geiger's matrix is < eigtol, the depth
c                 is not adjusted, and a message is printed.
c        rmscut - value for rms residual below which hypocentral adjustments
c                 are terminated.
c        zmin - minimum hypocenter depth
c        dxmax - maximum horizontal hypocentral adjustment allowed in
c                each hypocenter iteration.
c
c        rderr - estimate of reading error, used to estimate hypocenter error
c        ercof - for hypoinverse-like error calculations. Set > 0 and
c                < 1 if you want to include rms.res in hypocenter error
c                estimate. (sigsq=rderr**2 + ercof * rms**2)
c        line 3 - (free format)
c        nhitct - of observations for a parameter to be included
c                 in the inversion. This uses the variable khit.
c        dvpmx - maximum P-velocity adjustment allowed per iteration. 
c        dvsmx - maximum S-velocity adjustment allowed per iteration.
c        idmp -  set to 1 to recalculate the damping value for succeeding
c                iterations.  set to 0 to have constant damping
c        vdamp - damping parameter used in velocity inversion.
c                 vdamp(1)=damping for p-velocity
c                 vdamp(2)=damping for Vp/Vs
c                 vdamp(3)=damping for station delays
c        stepl - (km) used for calculation of partial derivatives along
c                the raypath.
c
c        line 4 - (free format)
c        ires - set to 1 to compute the resolution and print diagonal elements,  
c             also prints out 1/(diagonal covariance)
c           2 to print full resolution to file 17 (recomputed on each
c              iteration.)
c           3 to calculate full resolution on 1st iteration only.
c           0 no resolution calculations.
c        i3d - flag for using psuedo-bending:
c              0=no psuedo-bending
c              1=use in forward problem to compute velocity partial derivatives
c              2=also use in earthquake location subroutine
c              3=psuedo-bending; also use less curvature for initial arcuate
c                rays below "moho".  Assumes last z grid (k=nz-1) is "moho".
c        nitmax - max of iterations of the velocity inversion-hypocenter
c                 relocation loop.
c                 For locations only, set nitmax=0.
c                 For creating synthetic data, set nitmax= -1, and
c                   use file007 for hypocenters and stations.
c                   ***NOTE that synthetic data option has not been ***
c                   ***changed for S-P.  It computes S travel-times. ***
c        snrmct - cutoff value for solution norm. Program will stop
c                 iterating if the solution norm is less than snrmct.
c        ihomo - flag for using 1-d starting model (1=yes,0=no)
c               if flagged, on first ihomo iterations do 2-d art
c        rmstop - overall rms residual for termination of program.
c         
c        ifixl - number of velocity inversion steps to keep hypocenters
c                fixed at start of progressive inversion
c
c        line 5 - (free format)
c        delt1, delt2 - distance weighting factors. The weight is 1 for
c               x<delt1, but tapers linearly to 0 between delt1 and delt2.
c        res1,res2 - same pattern as above, but for residual weighting.
c
c        line 6 - (free format)
c        ndip - of rotation angles of the plane of the ray, which will
c               be computed in the exhaustive search for the fastest time.
c        iskip - of rotation angles which will be skipped.
c                ndip=9, iskip=3 will give a vertical plane, and 2 swung
c                at angles of 22.5 degrees on each side.
c                ndip=9, iskip=4 will give only the vertical plane, and should
c                be used for 1 dim models, to save computer time.
c           *** Most of the computer time is spent in the raytracing. Careful
c               selection of the starting model and use of iskip can save
c               a lot.
c         scale1  - set scale1 to the step length for the travel
c                           time computation. Set no larger than the grid
c                           spacing.
c         scale2 - scale for the number of paths tried in the raytracing.
c                  Cliff uses a value of 1, and this seems ok, but would
c                  need to be tested in detail. If scale2 is smaller, the
c                  number of paths increases, and the computation time also
c                  goes up.
c
c
c        line 7 - (free format)
c        xfac - Convergence enhancement factor for psuedo-bending
c               Cliff suggests 1.2 to 1.5
c        tlim - Cutoff value for travel time difference to terminate
c               iteration. (0.0005 to 0.002 s)
c        nitpb - maximum permitted iterations for psuedo-bending (5 to 10)
c                nitpb(1) for shorter raypaths
c                nitpb(2) for raypaths > delt1
c              
c        line 8 - (free format)
c        iusep - flag to tell whether or not to use P arrivals (and invert for
c                P velocities).  0 = no, 1 = yes.
c        iuses - flag to tell whether or not to use S arrivals (and invert for
c                S velocities).  0 = no, 1 = yes.
c        invdel - flag to control inclusion of station delays in the inversion.
c                 =0 to not include stn delays.
c                 =1 to include stn delays.
c
c     File 02 - Station Data
c        line 1 - (free format)
c         ltdo, oltm, lndo, olnm, rota - sets origin and rotation of the
c              coordinate system. Choose the origin at the lower right corner
c              of the region.
c              Y points in North direction.
c              X points West.
c              rota - angle of rotation counterclockwise (degrees). This is
c                     used to rotate the entire coordinate system.
c
c         line 2 - (free format)
c         nsts - number of stations in station list to follow.
c
c         station list: see documentation.
c
c     File 03 - velocity model
c         This contains the number of nodes in the x, y , and z directions.
c         Then the velocity model. See the Thurber's documentation for a complete
c         description.  In simulps12 this has Vp, followed by Vp/Vs.
c         Specify which nodes (if any) you want to fix velocity for:
c           2  3  4   means x2(2nd column), y3(3rd row), z4(4th layer)
c
c     File 04 - travel-time data for earthquakes. Note that for simulps12
c         S data should be made into S-P data.
c         Use program 'convert6' to convert hypo71 summary and phase data
c         to this format.
c         (should still work with old 'convert3' files also, as well as
c         Michelini's format)
c
c     File 07 - travel-time data for shots
c
c     File 08 - travel-time data for blasts
c
c  OUTPUT FILES:
c     File 16 - printed output, send to line printer
c
c     File 12 - hypocenters from each iteration, hypoinverse format with error ellipse
c
c     File 13 - Final hypocenters in hypo71 summary card format
c     File 15 - Written when kout3=1, A separate file for each event
c         containing all the raypath points, useful for plotting
c
c     File 17 - Full resolution matrix.  Created if ires.ge.2.
c
c     File 18 - Pointer from full nodes to inversion nodes.
c          Output when fixed nodes are used.
c     File 19 - Written when kout3=1, Contains the difference in travel 
c          time between ART and psuedo-bending
c     File 22 - Station data with new P and S delays.  This
c          can be used as input to future runs. NOT created if
c          invdel=0 (no station delays inverted for).
c
c     File 20 - Station residual output (created for kout2=0 or 1)
c
c     File 23 - Final velocity model.  This can be used as input
c          to future runs.
c
c     File 24 - Earthquake travel-time data for new hypocenters.
c          This can be used as input in future runs.
c
c     File 25 - New velocities in station format so can be plotted 
c          with qplot.
c
c     File 26 - List of observations that used maximum allowed 
c          psuedo-bending iterations (nitpb)
c
c     File 28 - Blast travel-time data with new origin times.
c
c     File 34 - Output file similar to HYPO71 listing file, which can
c          be used as input to FPFIT fault plane solution program.
c          Only for earthquakes and shots since called from Loceqk.
c          (Note that DIST is the hypocentral distance, whereas HYPO71
c          outputs the epicentral distance.)  Written on last (nitmax) iteration.
c
c     File 36 - Summary output file that contains key solution
c          statistics.  This is useful when doing numerous
c          damping runs and wanting to compare variances.
c
c     File 45 - 1/diagonal elements of covariance matrix.  In same
c          format as velocity model input.  Created if ires.ge.2.
c   HINTS:
c
c       1. Invert for a one-dimensional model first. Use nx=3, ny=3, nz=
c          whatever you want in the final model. Set ndip=9, iskip=4.
c
c          Or, better yet, use VELEST, locate your events with the VELEST
c          model, then fix the earthquake locations on the 1st iteration.
c
c       2. Run the inversion on calculated data. This will tell you what
c          kind of averaging is inherent in the inversion process, since
c          the result will probably be different from your original input
c          model.
c
c       3. When the program crashes due to divide by 0, etc, the cause can
c          most often be traced to errors in the setup of the velocity
c          model. No part of a ray must reach over halfway to an outer node.
c          Try increasing the distance of the outer node.
c
c
c  ARRAY DIMENSIONING:
c      11-mar-94 These notes are now in the simulps_common.inc file.
c      PLEASE READ that file to understand the array dimensions
c      and how many parameters you are allowed to invert for.
c
c
c  declaration statements:
      character*26 dash
c
c  common block variables:
      include 'simulps_common.inc'
      character*1 rmk(maxobs, maxev)
c
      dash='--------------------------'
c
c
c   open output files
        open(unit=16,file='fort.16',status='new',recl=132)
        rewind (16)
        open(unit=36,file='fort.36',status='new')
        rewind (36)
c    file to check psuedo-bending
        open(unit=26,file='fort.26',status='new',recl=132)
        rewind (20)
        rewind (26)
c
c  input routines
c
c  input control parameters
      open(unit=01,status='old',file='control.dat',form='formatted')
      rewind (01)
      call input1
      close(01)
c  open resolution file if ires=2
      if(ires.ge.2) open(unit=17,file='fort.17',status='new')
c  open residual output file if kout2=0,1
      if(kout2.le.1) open(unit=20,file='fort.20',status='new')
c  open files to write raypath points, tt differences to if kout3=1
      if(kout3.eq.0) goto 70
      open(unit=19,file='fort.19',status='new')
      write(19,1901)
 1901 format('  ne  stn  delta  fstime  ttime   tdif')
c  initializations
   70 call strt(0)
      nit=0
c  input station list, set up center of coordinates, calculate
c  cartesian coordinates
      open(unit=02,status='old',file='stations.dat',form='formatted')
      rewind (02)
      call input2
      close(02)
c  input medium model
      open(unit=03,status='old',file='velocity.dat',form='formatted')
      rewind (03)
      call input3
      close(03)
      if(neqs.eq.0) goto 71
      open(unit=04,status='old',file='earthq.dat',form='formatted')
   71 if(nsht.eq.0) goto 72
      open(unit=07,status='old',file='shots.dat',form='formatted')
   72 if(nbls.eq.0) goto 73
      open(unit=08,status='old',file='blasts.dat',form='formatted')
c
   73 if(kout.lt.3) goto 74
      open(unit=34,file='fort.34',status='new')
      rewind(34)
      if(kout.lt.4) goto 74
c       open file for hypocenter summary card output on each iteration
      open(unit=12,file='fort.22',status='new')
        rewind (12)
c
   74 istop=0
      istop1=0
c
      neb=neqs+nbls
      nevt=neb+nsht
c
c  iterative inversion loop
c
    1 continue
c
      netemp=neqs
      nbtemp=nbls
c  if flagged, on first ifixl iterations, treat all quakes as blasts
      if (ifixl.le.nit) go to 110
      nbls=nbls+neqs
      neqs=0
  110 continue
c
      istemp=iskip
      ndtemp=ndip
c  if flagged, on first ihomo iterations do 2-d art
      if (ihomo.le.nit) go to 111
      iskip = 0
      ndip = 1
  111 continue
c
c
      write(16,1005) nit
 1005 format(///,' iteration step',i3,'; hypocenter adjustments')
      if (neqs.eq.0) go to 11
c
c  loop over all earthquakes
      ne=1
    9 continue
c  input observations of event
      if((nit.eq.0).or.(kout2.eq.0)) write(16,130)
     2 dash,dash,dash,dash,dash
  130 format(1x,5a26)
      if (nit.eq.0) call input4(ne,4,rmk)
      if(ne.gt.neqs) goto 10
c  locate individual earthquake to reduce residuals
      call loceqk(ne,nit,nwr)
      if((nit.eq.0).and.(nwr.lt.4)) goto 9
      if((nit.gt.0).and.(nwr.lt.4)) goto 9998
c  skip if only locating
      if((nitmax.le.0).and.(kout3.eq.0).and.(kout2.gt.1)) goto 208
c  do forward problem and calculate partial derivatives
      call forwrd(ne)
      if((nitmax.le.0).and.(kout2.gt.1)) goto 208
c  perform parameter separation
      call parsep(ne,nwr)
c  if last iteration, write out station residuals
  208 if((nit.eq.nitmax).and.(kout2.lt.2)) call outres(ne,rmk)
      if((nit.eq.nitmax).and.(kout.ge.3)) call outlis(ne,rmk)
c  skip this event if not enough good readings,
c        stop if not first iteration
      if((nit.gt.0).and.(nwr.lt.4)) goto 9998
      if(nwr.ge.4) ne=ne+1
      if(ne.le.neqs) goto 9
   10 continue
c
c  loop over all blasts
   11 if(nbls.eq.0) goto 15
      write(16,1613) nbls
 1613 format(/,2x,'Following ',i5,' Events are Blasts with Unknown',
     2 ' Origin Time')
      nb=1
   12 ne=nb+neqs
   13 infile=8
      if((nit.eq.0).or.(kout2.eq.0)) write(16,130)
     2 dash,dash,dash,dash,dash
      if((ifixl.gt.nit).and.(ne.le.netemp)) infile=4
      if(nit.eq.0) then
        call input4(ne,infile,rmk)
        if(ne.gt.neb) goto 14
        if((ne.gt.netemp).and.(infile.eq.4)) goto 13
      endif
      call loceqk(ne,nit,nwr)
      if((nitmax.le.0).and.(kout3.eq.0).and.(kout2.gt.1)) goto 308
      call forwrd(ne)
      if((nitmax.le.0).and.(kout2.gt.1)) goto 308
      call parsep(ne,nwr)
  308 if((nit.eq.nitmax).and.(kout2.lt.2)) call outres(ne,rmk)
      if((nit.eq.nitmax).and.(kout.ge.3)) call outlis(ne,rmk)
      if(nwr.ge.4) nb=nb+1
      if(nb.le.nbls) goto 12
   14 continue
c
c  loop over all shots
   15 continue
      if (nsht.eq.0) go to 25
      write(16,1614) nsht
 1614 format(/,2x,'Following ',i5,' Events are Shots with Known',
     2 ' Origin Time')
      ns=1
c  input observations of shot
   16 ne=ns+neqs+nbls
      if((nit.eq.0).or.(kout2.eq.0)) write(16,130)
     2 dash,dash,dash,dash,dash
      if (nit.gt.0) goto 17
      call input4(ne,7,rmk)
c  do forward problem and calculate partial derivatives
      if(nit.eq.0) then
        jfl=0
      else
        jfl=2
        if(ihomo.eq.nit) jfl=1
      end if
   17 call forwrd(ne)
c  add medium derivatives from shot to medium matrix
      call medder(ne)
c** Change for synthetic data, nitmax= -1
      if((nitmax.le.0).and.(kout2.gt.1)) goto 408
c  put derivatives into g matrix
      call parsep(ne,nwr)
  408 ns=ns+1
      if((nit.eq.nitmax).and.(kout2.lt.2)) call outres(ne,rmk)
      if(ns.le.nsht) goto 16
   20 continue
c
   25 continue
c** Change for synthetic data, nitmax= -1
      if(nitmax.le.0) goto 9999
c  continue iterating or terminate?
      call decide(istop,nit)
c     if (nit.gt.0) call decide(istop)
        if(istop1.gt.2) goto 9999  !If backed-up twice, don't adjust again,end
        if(istop.eq.2) goto 28     !If variance ratio bad, backup
        if(istop1.eq.2) goto 9999  !If backed-up once and now okay, end
      if(nit.ge.nitmax) goto 9999  !As in CT version, do hyp last
c
 900  continue
      nit=nit+1
      write(16,1015) nit
 1015 format(///,' iteration step',i3,'; simultaneous inversion')
c  invert for velocity model adjustments
      call veladj(nit)
      if(nit.eq.99) goto 9999
c  output results of iteration
      call outadj(nit,istop,istop1)
c  compute resolution?
      if(ires.eq.0) goto 901
      if((ires.ne.3).or.(nit.eq.1)) call rescov
  901 if(istop.eq.1) go to 9999
      if(istop.eq.0) goto 34
c
   28   write(16,902)
  902   format(' ******* Variance Ratio is less than Critical Ratio',
     2 /,12x,' ******* Backup Parameters Halfway *******',/)
        call velbku(istop1)
        istop1=istop1+istop
        call outadj(nit,istop,istop1)
c
   34 continue
c
c  restore neqs,nbls,iskip,ndip
      neqs=netemp
      nbls=nbtemp
      iskip=istemp
      ndip=ndtemp
c
      call strt(nit)
      rewind(26)
c
      go to 1
c
 9998 continue
      write(16,1620)
 1620 format(//,'  ********** STOP **********')
 9999 continue
c
      call outend(rmk)
c
      close(04)
      close(07)
      close(12)
      close(16)
      close(26)
      close(36)
      if(ires.eq.2) close(17)
      if(kout3.eq.0) stop
      close(19)
      stop
c***** end of main program *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine avsd(cnull,x,nx,sd,av,devtot)
c  program to find the average and standard deviation of a 
c  list of numbers.  Those with value cnull are not included.
c
      real x(5000)
      sum=0
      i=0
      do 50 ix=1,nx
         if(x(ix).eq.cnull) goto 50
         sum=sum+x(ix)
         i=i+1
   50 continue
      nx1=nx
      nx=i
      if(nx.eq.0) goto 800
      av=sum/nx
      devtot=0
      do 260 i=1,nx1
         if(x(i).eq.cnull) goto 260 
         dev=x(i)-av
         devtot=devtot+dev*dev
  260 continue
      sd=sqrt(devtot/nx)
c     write(6,620) nx,av,sd
  620 format(' nx=',i5,', average=',e11.4,', sd=',e11.4)
      return
  800 continue
      sd=0.00
      devtot=0.00
      av=0.00
      nx=0.00
      return
c ***** end of subroutine avsd *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine aztoa(x1,x2,y1,y2,z1,z2,xr,yr,azim,tkofan)
c  this subroutine computes the azimuth and take-off-angle
c  for an individual observation
c  (called from Loceqk)
c  pr is station, p1 and p2 are hypocenter and adjoining point
c  on the raypath
c
      common/shortd/ xltkm,xlnkm,rota,xlt,xln,snr,csr
      parameter (drad=1.7453292d-02)
      parameter (pi=3.1415926536)
      parameter (twopi=6.2831853072)
c
c  Azimuth
      xd=xr-x1
      yd=yr-y1
      xda=abs(xd)
      yda=abs(yd)
      phi=atan(xda/yda)
c  compute correct azimuth depending on quadrant
      if(xd.ge.0.0) then
        if(yd.ge.0.0) then
          theta=twopi-phi
        else
         theta=pi+phi
        endif
      else
        if(yd.ge.0.0) then
          theta=phi
        else
          theta=pi-phi
        endif
      endif
c  rotate back to real north, convert to degrees
      azim=(theta-rota)/drad
      if(azim.gt.360.0) azim=azim-360.0
      if(azim.lt.0.0) azim=azim+360.0
c
c  Take-off-angle
      xd=x2-x1
      yd=y2-y1
      r=sqrt(xd*xd+yd*yd)
      zd=z2-z1
      zda=abs(zd)
      phi=atan(r/zda)
      if(zd.lt.0.0) phi=pi-phi
      tkofan=phi/drad
c
      return
c ***** end of subroutine aztoa *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine bend(isp,xfac)
c*****this routine perturbs the initial path in the direction
c      of the normal to the ray path tangent at each point
c      by the optimal distance r
c
      common/pathm/x(130),y(130),z(130),v(130),tra,n,nn
      common/temp/xtemp(130),ytemp(130),ztemp(130),rtemp(130),ttemp(130)
c
c ***
      xtemp(1)=x(1)
      ytemp(1)=y(1)
      ztemp(1)=z(1)
c ***
      do 200 k=2,nn
c
         kk=k-1
         kkk=k+1
c
c*****compute the normal direction of maximum gradient of velocity
c
         dx=x(kkk)-xtemp(kk)
         dy=y(kkk)-ytemp(kk)
         dz=z(kkk)-ztemp(kk)
         dn=dx*dx+dy*dy+dz*dz
         ddn=sqrt(dn)
         rdx=dx/ddn
         rdy=dy/ddn
         rdz=dz/ddn
c
         xk=0.5*dx+xtemp(kk)
         yk=0.5*dy+ytemp(kk)
         zk=0.5*dz+ztemp(kk)
c ***
         call vel3(isp,xk,yk,zk,vk)
         call veld(isp,xk,yk,zk,vx,vy,vz)
c
c ***
         vrd=vx*rdx+vy*rdy+vz*rdz
         rvx=vx-vrd*rdx
         rvy=vy-vrd*rdy
         rvz=vz-vrd*rdz
c
         rvs=sqrt(rvx*rvx+rvy*rvy+rvz*rvz)
         if(rvs.eq.0.0) goto 200
         rvx=rvx/rvs
         rvy=rvy/rvs
         rvz=rvz/rvs
c
c*****compute the optimal distance r
         rcur=vk/rvs
         rtemp(k)=rcur-sqrt(rcur*rcur-0.25*dn)
c
c*****compute the new points and distance of perturbations
c
         xxk=xk+rvx*rtemp(k)
         yyk=yk+rvy*rtemp(k)
         zzk=zk+rvz*rtemp(k)
c
c  convergence enhancement
         xxk=xfac*(xxk-x(k))+x(k)
         yyk=xfac*(yyk-y(k))+y(k)
         zzk=xfac*(zzk-z(k))+z(k)
c
         ttemp(k)=sqrt((x(k)-xxk)**2+(y(k)-yyk)**2+(z(k)-zzk)**2)
         xtemp(k)=xxk
         ytemp(k)=yyk
         ztemp(k)=zzk
         call vel3(isp,xxk,yyk,zzk,vk)
         v(k)=vk
c ***
200   continue
c
      return
c ***** end of subroutine bend *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine bldmap
c  common block variables:
      include 'simulps_common.inc'
c
c     array size limits
c     ixkms=iykms=izkms=1500
c
c     write(6,400)
c 400 format(' subroutine bldmap')
      xl=bld-xn(1)
      ixmax=(xn(nx)+xl)/bld
      yl=bld-yn(1)
      iymax=(yn(ny)+yl)/bld
      zl=bld-zn(1)
      izmax=(zn(nz)+zl)/bld
c     write(6,402)ixmax,iymax,izmax
c 402 format(' array sizes: ',3i5)
c
c  Check for array size overflow
      if(ixmax.gt.ixkms.or.iymax.gt.iykms.or.izmax.gt.izkms)goto 330
      ix=1
      do 10 i=1,ixmax
c
         ix1=ix+1
c
         xnow=float(i)*bld-xl
         if (xnow.ge.xn(ix1)) ix=ix1
c
         ixloc(i)=ix
   10 continue
c  Fill remainder of array with zeroes.
      do 12 i=ixmax,ixkms
         ixloc(i)=0
   12 continue
c
c
      iy=1
      do 15 i=1,iymax
c
         iy1=iy+1
c
         ynow=float(i)*bld-yl
         if (ynow.ge.yn(iy1)) iy=iy1
c
         iyloc(i)=iy
   15 continue
c
c  Fill rest of array with zeroes.
      do 17 i=iymax,iykms
         iyloc(i)=0
 17   continue
c
      iz=1
      do 20 i=1,izmax
c
         iz1=iz+1
c
         znow=float(i)*bld-zl
         if (znow.ge.zn(iz1)) iz=iz1
c
         izloc(i)=iz
   20 continue
c
c  Fill remainder of array with zeroes.
      do 22 i=izmax,izkms
         izloc(i)=0
  22  continue
      return
 330   continue
      write(16,331)ixkms,iykms,izkms
 331  format(' ***** error in array size in common/locate/',/,
     *' maximum map dimensions (km)=',/,' x=',i5,' y=',i5,' z=',i5)
      write(16,332)ixmax,iymax,izmax
  332 format(' Actual map size (km): ',/,' x=',i5,' y=',i5,' z=',i5)
      stop
c***** end of subroutine bldmap *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine cmpdpv(xe,ye,ze,xr,yr,zr,scale2,ndip,dipvec)
c
c  parameters
      real xe,ye,ze,xr,yr,zr,scale2,dipvec(3,9)
c
      integer ndip
c  local variables
      real dx,dy,dz,xh1,yh1,zh1,xh2,yh2,zh2,size,xv,yv,zv,
     *     rescal,x451,y451,z451,x452,y452,z452
c
      integer nv
c
      dx=xr-xe
      dy=yr-ye
      dz=zr-ze
c
c  near-vertical vector
      xv=-dx*dz
      yv=-dy*dz
      zv=dx*dx+dy*dy
c  rescale vector to length scale2
      size=sqrt(xv*xv+yv*yv+zv*zv)
      rescal=scale2/size
c
      xv=xv*rescal
      yv=yv*rescal
      zv=zv*rescal
c
c  store this vector
      nv=(ndip+1)/2
      dipvec(1,nv)=xv
      dipvec(2,nv)=yv
      dipvec(3,nv)=zv
c
      if (ndip.eq.1) return
c
c  horizontal vectors
      xh1=dy
      yh1=-dx
      zh1=0.0
      xh2=-dy
      yh2=dx
      zh2=0.0
c  rescale the vectors to length scale2
      size=sqrt(xh1*xh1+yh1*yh1)
      rescal=scale2/size
c
      xh1=xh1*rescal
      yh1=yh1*rescal
      xh2=xh2*rescal
      yh2=yh2*rescal
c
c  store these two vectors
      dipvec(1,1)=xh1
      dipvec(2,1)=yh1
      dipvec(3,1)=zh1
c
      dipvec(1,ndip)=xh2
      dipvec(2,ndip)=yh2
      dipvec(3,ndip)=zh2
c
      if (ndip.eq.3) return
c
c  determine two 45 degree dip vectors
      rescal=0.7071068
c
      n1=(1+nv)/2
      n2=(nv+ndip)/2
c
      x451=(xh1+xv)*rescal
      y451=(yh1+yv)*rescal
      z451=(zh1+zv)*rescal
c
      x452=(xh2+xv)*rescal
      y452=(yh2+yv)*rescal
      z452=(zh2+zv)*rescal
c
      dipvec(1,n1)=x451
      dipvec(2,n1)=y451
      dipvec(3,n1)=z451
c
      dipvec(1,n2)=x452
      dipvec(2,n2)=y452
      dipvec(3,n2)=z452
c
      if (ndip.eq.5) return
c
c  determine four 22.5 degree dip vectors
      rescal=0.5411961
c
      dipvec(1,2)=(xh1+x451)*rescal
      dipvec(2,2)=(yh1+y451)*rescal
      dipvec(3,2)=(zh1+z451)*rescal
c
      dipvec(1,4)=(x451+xv)*rescal
      dipvec(2,4)=(y451+yv)*rescal
      dipvec(3,4)=(z451+zv)*rescal
c
      dipvec(1,6)=(xv+x452)*rescal
      dipvec(2,6)=(yv+y452)*rescal
      dipvec(3,6)=(zv+z452)*rescal
c
      dipvec(1,8)=(x452+xh2)*rescal
      dipvec(2,8)=(y452+yh2)*rescal
      dipvec(3,8)=(z452+zh2)*rescal
c
c***** end of subroutine cmpdpv *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine cmpdsv(ndip,iskip,ns,dipvec,disvec)
c
c  parameters
      real dipvec(3,9),disvec(390,9)
c
      integer ndip,iskip,ns
c  local variables
      real darc(129)
c
      integer inc,narc,ndp,np,n,nd
c  coefficients of standard arc
c
       data darc/0.,0.0342346,0.0676860,0.1003707,0.1323047,0.1635029,
     *.1939796,.2237485,.2528226,.2812141,.3089350,.3359963,
     *.3624091,.3881833,.4133289,.4378553,.4617713,.4850857,
     *.5078065,.5299417,.5514988,.5724850,.5929070,.6127717,
     *.6320853,.6508538,.6690831,.6867787,.7039459,.7205898,
     *.7367154,.7523272,.7674298,.7820274,.7961241,.8097238,
     *.8228301,.8354468,.8475771,.8592244,.8703916,.8810817,
     *.8912975,.9010416,.9103164,.9191245,.9274679,.9353487,
     *.9427691,.9497307,.9562353,.9622845,.9678797,.9730224,
     *.9777138,.9819550,.9857470,.9890908,.9919872,.9944367,.9964401,
     *.9979979,.9991102,.9997776,1.0000000,.9997776,.9991102,.9979979,
     *.9964401,.9944367,.9919872,.9890908,.9857470,.9819550,.9777138,
     *.9730224,.9678797,.9622845,.9562353,.9497307,.9427691,.9353487,
     *.9274679,.9191245,.9103164,.9010416,.8912975,.8810817,.8703916,
     *.8592244,.8475771,.8354468,.8228301,.8097238,.7961241,.7820274,
     *.7674298,.7523272,.7367154,.7205898,.7039459,.6867787,.6690831,
     *.6508538,.6320853,.6127717,.5929070,.5724850,.5514988,.5299417,
     *.5078065,.4850857,.4617713,.4378553,.4133289,.3881833,.3624091,
     *.3359963,.3089350,.2812141,.2528226,.2237485,.1939796,.1635029,
     *.1323047,.1003707,.0676860,.0342346,0.0/
      inc=128/ns
        ndip1=1+iskip
      ndip2=ndip-iskip
c
c  loop over dips
      do 30 ndp=ndip1,ndip2
         narc=1
         nd=3
c  loop over points on the path (skip first and last)
         do 20 np=2,ns
            narc=narc+inc
c
c  loop over x,y,z
            do 10 n=1,3
               nd=nd+1
               disvec(nd,ndp)=darc(narc)*dipvec(n,ndp)
c
   10       continue
   20    continue
   30 continue
c
c***** end of subroutine cmpdsv *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine cmpsep(path,pthsep,ns)
c
c  parameters
      real path(390),pthsep(130)
c
      integer ns
c  local variables
      integer nx,ny,nz,nx1,ny1,nz1
c
      nx=-2
c  loop over pairs of points in one set of stored vectors
      do 10 n=1,ns
         nx=nx+3
         ny=nx+1
         nz=nx+2
         nx1=nx+3
         ny1=nx+4
         nz1=nx+5
c
         pthsep(n)=sqrt((path(nx1)-path(nx))**2
     *              +(path(ny1)-path(ny))**2
     *              +(path(nz1)-path(nz))**2)
c
   10 continue
c
c***** end of subroutine cmpsep *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine curvdr(isp,nc,ndp,dtt,tt1,tt2,disvec,pthsep,
     *  strpth,trpth1)
c  This subroutine computes the difference in time between the nc
c  and the nc+1 curve with dip = ndp. dtt is the tt difference,
c  tt1 is the tt of the nc curve, and tt2 is the tt of the nc+1 curve.
      dimension disvec(390,9),pthsep(130),strpth(390),trpth1(390)
      common/raytr/trpath(390,9),npt,ns
      ncv=nc
      call curvtm(isp,ncv,ndp,tt,disvec,pthsep,strpth,trpth1)
      tt1=tt
      nc1=nc+1
      call curvtm(isp,nc1,ndp,tt,disvec,pthsep,strpth,trpth1)
      tt2=tt
      dtt=tt2-tt1
      return
c***** end of subroutine curvdr *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine curvtm(isp,nc,ndp,ttm,disvec,pthsep,strpth,trpth1)
c  This computes the travel time for curve nc at dip ndp. It is
c  used in the faster search programmed by Prothero.
      dimension disvec(390,9),pthsep(130),strpth(390),trpth1(390)
      common/raytr/trpath(390,9),npt,ns
c  loop to determine points along one path
      npt2=npt-2
      do 45 np=1,npt2
         n1=3*np+1
         n3=n1+2
         do 44 nn=n1,n3
            trpath(nn,ndp)=nc*disvec(nn,ndp)+strpth(nn)
            trpth1(nn)=trpath(nn,ndp)
  44     continue
  45  continue
c  set up pthsep array for travel time calculations
      call cmpsep(trpth1,pthsep,ns)
c  compute travel time along the path
      call ttime(isp,ns,npt,trpth1,pthsep,tt)
      ttm=tt
      return
c***** end of subroutine curvtm *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine decide(istop,nit)
c  common block variables:
      include 'simulps_common.inc'
c
      istop=0
      rms=sqrt(ssqr/float(nobt))
      dvar=ssqrw/wnobt
      rmsw=sqrt(dvar)
      write(16,1610) rms,rmsw,dvar
 1610 format(//,' unweighted rms=',f8.5,'; weighted rms=',f8.5,
     2 ' data var.(ssqrw/wnobt ie:rms**2)=',f8.5)
      dvarp=ssqrwp/wnobtp
      if(iuses.eq.1) then
        write(16,1611) dvarp
 1611   format(50x,'P data var.=',f8.5)
      else
        dvars=ssqrws/wnobts
        write(16,1612) dvarp,dvars
 1612   format(50x,'P data var.=',f8.5,'  S-P data var.=',f8.5)
      endif
      if (rmsw.lt.rmstop) istop=1
c
c  f-test
      mbl0=0
c  mbl1 was wrong, did not include stations, 1-april-1983, dmep
      do 10 n=1,npar
         if(khit(n).eq.0) goto 10
         if ((hit(n).lt.hitct).or.(nfix(n).eq.1)) go to 10
         mbl0=mbl0+1
   10 continue
      ndof=nobt-4*neqs-mbl0-nbls
      wndof=wnobt-float(4*neqs+mbl0+nbls)
        write(16,1601) mbl0,mbl1
 1601 format(' subroutine decide, mbl0 = ',i6,' mbl1 = ',i6)
      var=ssqr/ndof
      varw=ssqrw/wndof
      write(16,1606) ssqr,var,ssqrw,varw
 1606 format(/,'  Unweighted:  ssqr =',f12.3,'    var. (ssqr/ndof) =',
     2  f8.5,/,
     3  '  Weighted:   ssqrw =',f12.3,'  varw.(ssqrw/wndof) =',f8.5,/)
c
      if(nit.eq.0) then 
        write(36,3606) nit,ssqrw,varw
 3606   format(' Iteration:',i3,', ssqrw =',f12.3,
     *    '  varw.(ssqrw/wndof) =',f8.5)
        return
      endif
c
      call ftest(ndof,ndof1,ratio)
c
      rat=var1/var
      ratw=varw1/varw
c      write(16,6601) nit, var,ndof,var1,ndof1,rat,ratio
 6601 format(//,' f-test iteration',i3,/,
     *  ' new variance and ndof =',f10.5,i6,
     *  /,' old variance and ndof =',f10.5,i6,
     *  /,' variance ratio and critical ratio =',2f10.3)
      write(16,1605) 
 1605 format (//,'*****WEIGHTED****** use for f-test since weighted ',
     2  'throughout inversion',/)
      iwndof=nint(wndof)
      iwndf1=nint(wndof1)
      call ftest(iwndof,iwndf1,ratio)
      write(16,6602) nit,varw,wndof,varw1,wndof1,ratw,ratio
 6602 format(//,' f-test iteration',i3,/,
     *  ' new variance and ndof =',f10.5,f7.0,
     *  /,' old variance and ndof =',f10.5,f7.0,
     *  /,' variance ratio and critical ratio =',2f10.3)
      write(36,6601) nit,varw,ndof,varw1,ndof1,ratw,ratio
c
      if (ratw.lt.ratio) istop=2
c
c***** end of subroutine decide *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine disto(j,n,x,y)
c  this routine calculates distance of station or event
c  from given coordinate origin in terms of (possibly
c  rotated) cartesian coords x and y
c  uses short distance conversion factors from setorg
c
c  declaration statements:
      double precision drad,drlt
      parameter (drad=1.7453292d-02)
      parameter (drlt=9.9330647d-01)
c
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,xlt,xln,snr,csr
      include 'simulps_common.inc'
c
c  if j=1, calculate coords of station n w.r.t.
c  center of coords, doing rotation if any
c  otherwise, calculate coords of event n
      if (j.ne.1) go to 15
      plt=60.*ltds(n)+sltm(n)
      pln=60.*lnds(n)+slnm(n)
      go to 25
   15 continue
      plt=60.*ltde(n)+eltm(n)
      pln=60.*lnde(n)+elnm(n)
   25 continue
c  now convert lat and lon differences to km
      x=pln-xln
      y=plt-xlt
      xlt1=atan(drlt*tan(drad*(plt+xlt)/120.))
      x=x*xlnkm*cos(xlt1)
      y=y*xltkm
c  now do rotation
      if(rota.eq.0.0) return
      ty=csr*y+snr*x
      x=csr*x-snr*y
      y=ty
      return
c***** end of subroutine disto *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine fksvd (a, s, v, mmax, nmax, m, n, p, withu, withv)
c
c  common block variables:
      include 'simulps_common.inc'
c
      common/machin/ eta,tol
c
      integer    mmax, nmax, m, n, p
      real       r, w, cs, sn, tol, f, x, eps, gf, t, y
      real       eta, h, q, z
      integer    i, j, k, l, l1, n1, np
      logical    withu, withv
      double precision a(maxobs,maxobs)
      dimension        s(4), v(4,4)
      dimension as1(maxobs),as2(maxobs)
      dimension  t(maxobs)
c
c     ------------------------------------------------------------------
c
c     this is a translation of a cdc 6600 fortran program to ibm 360
c     fortran iv.  this subroutine uses short precision arithmetic.
c     a long precision version is available under the name 'dsvd'.
c
c     this subroutine replaces earlier subroutines with the same name,
c    689   6       &   &0&  s of a complex ar&thmetic program, published
c     as algorithm 358.  this current program is faster, more accurate
c     and less obscure in describing its capabilities.
c
c     original programmer=  r. c. singleton
c     360 version by=       j. g. lewis
c     last revision of this subroutine=  4 december 1973
c
c     ------------------------------------------------------------------
c
c     additional subroutine needed=  rotate
c
c     ------------------------------------------------------------------
c
c
c     this subroutine computes the singular value decomposition
c     of a real m*n matrix a, i.e. it computes matrices u, s, and v
c     such that
c
c                  a = u * s * vt ,
c     where
c              u is an m*n matrix and ut*u = i, (ut=transpose
c                                                    of u),
c              v is an n*n matrix and vt*v = i, (vt=transpose
c                                                    of v),
c        and   s is an n*n diagonal matrix.
c
c     description of parameters=
c
c     a = real array. a contains the matrix to be decomposed.
c         the original data are lost.  if withv=.true., then
c         the matrix u is computed and stored in the array a.
c
c     mmax = integer variable.  the number of rows in the
c            array a.
c
c     nmax = integer variable.  the number of rows in the
c            array v.
c
c     m,n = integer variables.  the number of rows and columns
c           in the matrix stored in a.  (ng=mg=100.  if it is
c           necessary to solve a larger problem, then the
c           amount of storage allocated to the array t must
c           be increased accordingly.)  if mlt n , then either
c           transpose the matrix a or add rows of zeros to
c           increase m to n.
c
c     p = integer variable.  if p'0, then columns n+1, . . . ,
c         n+p of a are assumed to contain the columns of an m*p
c         matrix b.  this matrix is multiplied by ut, and upon
c         exit, a contains in these same columns the n*p matrix
c         ut*b. (p'=0)
c
c     withu, withv = logical variables.  if withu=.true., then
c         the matrix u is computed and stored in the array a.
c         if withv=.true., then the matrix v is computed and
c         stored in the array v.
c
c     s = real array.  s(1), . . . , s(n) contain the diagonal
c         elements of the matrix s ordered so than s(i)>=s(i+1),
c         i=1, . . . , n-1.
c
c     v = real array.  v contains the matrix v.  if withu
c         and withv are not both =.true., then the actual
c         parameter corresponding to a and v may be the same.
c
c     this subroutine is a real version of a fortran subroutine
c     by businger and golub, algorithm 358=  singular value
c     decomposition of a complex matrix, comm. acm, v. 12,
c     no. 10, pp. 564-565 (oct. 1969).
c     with revisions by rc singleton, may 1972.
c     ------------------------------------------------------------------
c
c
c     VAX version:  machine constants calculated internally in strt
c     eta is the machine epsilon (relative accuracy)
c     tol is the smallest representable real divided by eta.
c
      np = n + p
      n1 = n + 1
c
c     householder reduction to bidiagonal form
      gf = 0.0
      eps = 0.0
      l = 1
   10 t(l) = gf
      k = l
      l = l + 1
c
c     elimination of a(i,k), i=k+1, . . . , m
      s(k) = 0.0
      z = 0.0
      do 20 i = k,m
         z = z + a(i,k)**2
   20 continue
      if (z.lt.tol) goto 50
      gf = sqrt(z)
      f = a(k,k)
      if (f.ge.0.0) gf = - gf
      s(k) = gf
      h = gf * (f - gf)
      a(k,k) = f - gf
      if (k.eq.np) goto 50
      do 45 j = l,np
         f = 0
         do 30 i = k,m
            f = f + a(i,k)*a(i,j)
   30    continue
         f = f/h
         do 40 i = k,m
            a(i,j) = a(i,j) + f*a(i,k)
   40    continue
   45 continue
c
c     elimination of a(k,j), j=k+2, . . . , n
   50 eps = amax1(eps,abs(s(k)) + abs(t(k)))
      if (k.eq.n) goto 100
      gf = 0.0
      z = 0.0
      do 60 j = l,n
         z = z + a(k,j)**2
   60 continue
      if (z.lt.tol) goto 10
      gf = sqrt(z)
      f = a(k,l)
      if (f.ge.0.0) gf = - gf
      h = gf * (f - gf)
      a(k,l) = f - gf
      do 70 j = l,n
         t(j) = a(k,j)/h
   70 continue
      do 95 i = l,m
         f = 0
         do 80 j = l,n
            f = f + a(k,j)*a(i,j)
   80    continue
      if (abs(f).lt.1.0e-25) f=0.
         do 90 j = l,n
            a(i,j) = a(i,j) + f*t(j)
   90    continue
   95 continue
c
      goto 10
c
c     tolerance for negligible elements
  100 eps = eps*eta
c
c     accumulation of transformations
      if (.not.withv) goto 160
      k = n
      goto 140
  110    if (t(l).eq.0.0) goto 140
         h = a(k,l)*t(l)
         do 135 j = l,n
            q = 0
            do 120 i = l,n
               q = q + a(k,i)*v(i,j)
  120       continue
            q = q/h
            do 130 i = l,n
               v(i,j) = v(i,j) + q*a(k,i)
  130       continue
  135    continue
  140    do 151 j = 1,n
            v(k,j) = 0
  151    continue
         v(k,k) = 1.0
         l = k
         k = k - 1
         if (k.ne.0) goto 110
c
  160 k = n
      if (.not.withu) goto 230
      gf = s(n)
      if (gf.ne.0.0) gf = 1.0/gf
      go to 210
  170    do 180 j = l,n
            a(k,j) = 0
  180    continue
         gf = s(k)
         if (gf.eq.0.0) goto 210
         h = a(k,k)*gf
         do 205 j = l,n
            q = 0
            do 190 i = l,m
               q = q + a(i,k)*a(i,j)
  190       continue
            q = q/h
            do 200 i = k,m
               a(i,j) = a(i,j) + q*a(i,k)
  200       continue
  205    continue
         gf = 1.0/gf
  210    do 220 j = k,m
            a(j,k) = a(j,k)*gf
  220    continue
         a(k,k) = a(k,k) + 1.0
         l = k
         k = k - 1
         if (k.ne.0) goto 170
c
c     qr diagonalization
      k = n
c
c     test for split
  230    l = k
  240       if (abs(t(l)).le.eps) goto 290
            l = l - 1
            if (abs(s(l)).gt.eps) goto 240
c
c     cancellation
         cs = 0.0
         sn = 1.0
         l1 = l
         l = l + 1
         do 280 i = l,k
            f = sn*t(i)
            t(i) = cs*t(i)
            if (abs(f).le.eps) goto 290
            h = s(i)
            w = sqrt(f*f + h*h)
            s(i) = w
            cs = h/w
            sn = - f/w
            do 221 ia=1,m
               as1(ia) = sngl(a(ia,l1))
               as2(ia) = sngl(a(ia,i))
  221       continue
            if (withu) call hyrot(mmax,as1, as2, cs, sn, m)
            do 222 ia=1,m
               a(ia,l1)=dble(as1(ia))
               a(ia,i)=dble(as2(ia))
  222       continue
            if (np.eq.n) goto 280
            do 270 j = n1,np
               q = a(l1,j)
               r = a(i,j)
               a(l1,j) = q*cs + r*sn
               a(i,j) = r*cs - q*sn
  270       continue
  280    continue
c
c     test for convergence
  290    w = s(k)
         if (l.eq.k) goto 360
c
c     origin shift
         x = s(l)
         y = s(k-1)
         gf = t(k-1)
         h = t(k)
         f = ((y - w)*(y + w) + (gf - h)*(gf + h))/(2.0*h*y)
         gf = sqrt(f*f + 1.0)
         if (f.lt.0.0) gf = - gf
         f = ((x - w)*(x + w) + (y/(f + gf) - h)*h)/x
c
c     qr step
         cs = 1.0
         sn = 1.0
         l1 = l + 1
         do 350 i = l1,k
            gf = t(i)
            y = s(i)
            h = sn*gf
            gf = cs*gf
            w = sqrt(h*h + f*f)
            t(i-1) = w
            cs = f/w
            sn = h/w
            f = x*cs + gf*sn
            gf = gf*cs - x*sn
            h = y*sn
            y = y*cs
            if (withv) call hyrot(nmax,v(1,i-1),v(1,i),cs,sn,n)
            w = sqrt(h*h + f*f)
            s(i-1) = w
            cs = f/w
            sn = h/w
            f = cs*gf + sn*y
            x = cs*y - sn*gf
            do 338 ia=1,m
               as1(ia) = sngl(a(ia,i-1))
               as2(ia) = sngl(a(ia,i))
  338       continue
            if (withu) call hyrot(mmax,as1, as2, cs, sn, m)
            do 339 ia=1,m
               a(ia,i-1)=dble(as1(ia))
               a(ia,i)=dble(as2(ia))
  339       continue
            if (n.eq.np) goto 350
            do 340 j = n1,np
               q = a(i-1,j)
               r = a(i,j)
               a(i-1,j) = q*cs + r*sn
               a(i,j) = r*cs - q*sn
  340       continue
  350    continue
c
         t(l) = 0.0
         t(k) = f
         s(k) = x
         goto 230
c
c     convergence
  360    if (w.ge.0.0) goto 380
         s(k) = - w
         if (.not.withv) goto 380
         do 370 j = 1,n
            v(j,k) = - v(j,k)
  370    continue
  380    k = k - 1
         if (k.ne.0) go to 230
c
c     sort singular values
      do 450 k = 1,n
         gf = -1.0
         do 390 i = k,n
            if (s(i).lt.gf) goto 390
            gf = s(i)
            j = i
  390    continue
         if (j .eq. k) go to 450
         s(j) = s(k)
         s(k) = gf
         if (.not.withv) goto 410
         do 400 i = 1,n
            q = v(i,j)
            v(i,j) = v(i,k)
            v(i,k) = q
  400    continue
  410    if (.not.withu) goto 430
         do 420 i = 1,m
            q = a(i,j)
            a(i,j) = a(i,k)
            a(i,k) = q
  420    continue
  430    if (n.eq.np) goto 450
         do 440 i = n1,np
            q = a(j,i)
            a(j,i) = a(k,i)
            a(k,i) = q
  440    continue
  450 continue
c
      return
c***** end of subroutine fksvd *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine forwrd(ne)
c  this routine calculates the forward problem
c  given the stations, and initial hypocenters and
c  velocity model, the routine calculates theoretical
c  travel times in a way appropriate to the assumed
c  form of the velocity model.  three-d ray tracing
c  is performed. travel time derivatives with respect to
c  hypocentral and medium parameters are calculated
c
c  declaration statements:
      character*8 fil15
      character*1 phs(2)
c
c  common block variables:
      include 'simulps_common.inc'
c
      data phs/'P','S'/
      zoff=99.0
c  event coordinates
      xe=evc(1,ne)
      ye=evc(2,ne)
      ze=evc(3,ne)
c  loop over all observations of this event
      nobs=kobs(ne)
      write(26,2601) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),
     2 seco(ne),ltde(ne),eltm(ne),lnde(ne),elnm(ne)
 2601 format(3h **,1x,i3,1x,a4,a2,1x,a2,i2,f6.2,i3,'n',f5.2,i4,'w',f5.2,
     2 2f7.2,3f6.2,3x,3i3)
      if(kout3.eq.0) goto 20
      write(fil15(1:8),1500) ne
 1500 format('ev',i3,'.rp')
      open(unit=15,file=fil15,status='new')
      write(15,1800) zoff,zoff,zoff
 1800 format('RAYPATH POINTS:',t68,f10.5,/,
     2 ' lat lon z x y -z',t68,f10.5,/,
     3 'format:17x,6f10.5',t68,f10.5)
   20 do 77 no=1,nobs
         call path(ne,no,xe,ye,ze,ttime)
c** Change for synthetic data, nitmax= -1
         if(nitmax.eq.-1) then
            secp(no,ne)=ttime+seco(ne)
c * * Different for S-P synthetic data * *
            isp=intsp(no,ne)
            if(isp.eq.1) then
               intsp(no,ne)=0
               call path(ne,no,xe,ye,ze,ptime)
               intsp(no,ne)=1
               smptime=ttime-ptime
               secp(no,ne)=smptime
            endif
         else
            call ttmder(ne,no,1,ttime)
         endif
         if(kout3.eq.0) goto 77
c  write out raypath points
         write(15,1801) ne,no,zoff
 1801    format(' ev=',i4,' obs=',i4,t68,f10.5)
         write(15,1803) iyrmo(ne),iday(ne),ihr(ne),mino(ne),
     2   seco(ne),zoff
 1803    format(a4,a2,1x,a2,i2,1x,f5.2,t68,f10.5)
         write(15,1804) stn(isto(no,ne)),phs(intsp(no,ne)+1),
     2   dlta(no,ne),zoff
 1804    format(a4,a1,f6.2,t68,f10.5)
         do 50 i=1,nrp(no)
            call latlon(rp(1,i,no),rp(2,i,no),lat,xlat,lon,xlon)
            xlat=float(lat)+xlat/60.0
            xlon=float(lon)+xlon/60.0
c  also print out negative z for easy plotting
            zd= -1.0 * rp(3,i,no)
            write(15,1802) xlat,xlon,rp(3,i,no),rp(1,i,no),
     2         rp(2,i,no),zd
   50    continue
 1802    format(17x,6f10.5)
   77 continue
      if(kout3.eq.1) close(15)
      return
c***** end of subroutine forwrd *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ftest(ndf,ndi,ratio)
c
c  interpolate in ftest ratio table
c  ndi is the initial number of degrees of freedom
c  ndf is the final number of degrees of freedom
c
      dimension rattab(4,4),valu(4)
      data rattab/1.69,1.64,1.58,1.51,1.59,1.53,1.47,1.39,
     *            1.50,1.43,1.35,1.25,1.39,1.32,1.22,1.00/
      data valu/0.025,0.016667,0.008333,0.00/
c
      if (ndi.ge.40.and.ndf.ge.40) go to 10
c
c  number of degrees of freedom outside table range
      ratio=2.0
      return
c
c  determine points in table for interploating
   10 continue
      index=ndi/30.
      index1=ndf/30.
      if (index.gt.4) index=4
      if (index1.gt.4) index1=4
c
      go to (11,12,12,14), index
c
   11 continue
      m=1
      m1=2
      go to 20
c
   12 continue
      m=2
      m1=3
      go to 20
c
   14 continue
      m=3
      m1=4
c
   20 continue
c
      go to (21,22,22,24), index1
c
   21 continue
      n=1
      n1=2
      go to 30
c
   22 continue
      n=2
      n1=3
      go to 30
c
   24 continue
      n=3
      n1=4
c
   30 continue
c
c  compute interpolated value from table
      fm=1.0/ndi
      fn=1.0/ndf
c
      vm=valu(m)
      vm1=valu(m1)
c
c  interpolate along m
      vf=rattab(m,n)+(rattab(m1,n)-rattab(m,n))*
     *   (fm-vm)/(vm1-vm)
      vf1=rattab(m,n1)+(rattab(m1,n1)-rattab(m,n1))*
     *    (fm-vm)/(vm1-vm)
c
c
      vn=valu(n)
      vn1=valu(n1)
c  interpolate across n
      ratio=vf+(vf1-vf)*(fn-vn)/(vn1-vn)
c***** end of subroutine ftest *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine h12(mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv,
     2               c1,ic1e,ic1v,nc1v)
c  c.l. lawson and r.j. hanson, jpl
c  from "solving least squares problems"
c  modified by c.h. thurber, spring 1980
c  construction and application of a single
c  householder transformation...  q =  i + u*(u**t)/b
      dimension u(m),c(*),c1(*)
      double precision sm,b
      one=1.
c
      if (0.ge.lpivot.or.lpivot.ge.l1.or.l1.gt.m) return
      cl=abs(u(lpivot))
      if (mode.eq.2) go to 60
c
c  construct the transformation
c
      do 10 j=l1,m
         cl=amax1(abs(u(j)),cl)
   10 continue
      if (cl) 130,130,20
   20 clinv=one/cl
      sm=(dble(u(lpivot))*clinv)**2
      do 30 j=l1,m
         sm=sm+(dble(u(j))*clinv)**2
   30 continue
c
c  convert dble prec sm to sngl prec sm1
c
      sm1=sm
      cl=cl*sqrt(sm1)
      if (u(lpivot)) 50,50,40
   40 cl=-cl
   50 up=u(lpivot)-cl
      u(lpivot)=cl
      go to 70
c
c  apply the transformation i+u*(u**t)/b to matrices c & c1
c
   60 if (cl) 130,130,70
   70 if (ncv.le.0) return
      b=dble(up)*u(lpivot)
c
c  b must be nonpositive here.  if b=0 return
c
      if (b) 80,130,130
   80 b=one/b
      i2=1-icv+ice*(lpivot-1)
      incr=ice*(l1-lpivot)
      do 120 j=1,ncv
         i2=i2+icv
         i3=i2+incr
         i4=i3
         sm=c(i2)*dble(up)
         do 90 i=l1,m
            sm=sm+c(i3)*dble(u(i))
            i3=i3+ice
   90    continue
         if (sm) 100,120,100
  100    sm=sm*b
         c(i2)=c(i2)+sm*dble(up)
         do 110 i=l1,m
            c(i4)=c(i4)+sm*dble(u(i))
            i4=i4+ice
  110    continue
  120 continue
      i2=1-ic1v+ic1e*(lpivot-1)
      incr=ic1e*(l1-lpivot)
      do 220 j=1,nc1v
         i2=i2+ic1v
         i3=i2+incr
         i4=i3
         sm=c1(i2)*dble(up)
         do 190 i=l1,m
            sm=sm+c1(i3)*dble(u(i))
            i3=i3+ic1e
  190    continue
         if (sm) 200,220,200
  200    sm=sm*b
         c1(i2)=c1(i2)+sm*dble(up)
         do 210 i=l1,m
            c1(i4)=c1(i4)+sm*dble(u(i))
            i4=i4+ic1e
  210    continue
  220 continue
  130 return
c***** end of subroutine h12 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine hyrot  (nmax,x, y, cs, sn, n)
      integer n
      real    x(nmax), y(nmax), cs, sn
c
c
      real    xx
      integer j
c
c
      do 10 j = 1, n
         xx = x(j)
         x(j) = xx*cs + y(j)*sn
         y(j) = y(j)*cs - xx*sn
   10 continue
      return
c***** end of subroutine hyrot *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine input1
c  this routine reads in control parameters and number of eq's
c
c  declaration statements:
      character*9 day
      character*8 tm
c
c  common block variables:
      include 'simulps_common.inc'
c
      call date_and_time(DATE=day)
      call date_and_time(TIME=tm)
      write(16,1605) day,tm
 1605 format(' Computation began at ',a9,1x,a8)
      write(16,9937)maxpar,mxpari,maxev,maxsta,maxobs
 9937  format(' Program Simulps12 (30-Jun-94 DMEP) Solves for Vp and '
     2  ,'Vp/Vs; Input data is P travel-time and S-P time.',/,
     3 '       Can vary relative weighting of S-P times in '
     4  , 'hypocenter location.',/,
     5 '       Allows fixed nodes (up to ',i5,' parameters,'
     *  , ' up to ',i5,' solution parameters);'/,
     6 '       up to ',i4,' events, ',i5,' stations, ',i4,
     *  ' observations per event',
     7 /,'       Psuedo-bending; Allows less curvature "moho" for ',
     8 'initial arcuate paths.')
      read(1,*,err=999) neqs,nsht,nbls,wtsht,kout,kout2,kout3
 1033 format(2i3,f4.1)
        if(neqs+nsht+nbls.le.maxev)go to 20
        write(16,22)
 22     format('0too many events for program arrays.')
        stop
  20    continue
c     write(16,9938)
 9938 format(/,'  program simul3m (jan 1985)',/,'  * fast art *',
     * /,' * pseudo-bending *')
      read(1,*) nitloc,wtsp,eigtol,rmscut,zmin,dxmax,rderr,ercof
      read(1,*) hitct,dvpmx,dvsmx,idmp,(vdamp(j),j=1,3),stepl
      read(1,*) ires,i3d,nitmax,snrmct,ihomo,rmstop,ifixl
      read(1,*) delt1,delt2,res1,res2,res3
      read(1,*) ndip,iskip,scale1,scale2
      read(1,*) xfac,tlim,nitpb(1),nitpb(2)
      read(1,*) iusep,iuses,invdel
 1051 format(2f7.2,2f5.2)
 1011 format(4i3,f5.3)
 1021 format(i3,3f5.3)
      write(16,1030)
      write(16,1040) kout,kout2,kout3
      write(16,1031)
      write(16,1041) neqs,nsht,nbls,wtsht,nitloc,wtsp,zmin,eigtol,
     2 rmscut,hitct,dvpmx,dvsmx
      write(16,1032)
      write(16,1042) idmp,(vdamp(j),j=1,3),ires,
     2 nitmax,snrmct,dxmax,rderr,ercof,ihomo,rmstop,ifixl,delt1,delt2,
     3 res1,res2,res3
      write(16,1061) stepl
      write(16,1155)
 1155 format(' parameters for approximate ray tracer',
     * /,6x,'ndip iskip scale1 scale2')
      write(16,1151) ndip,iskip,scale1,scale2
 1151 format(5x,i3,i6,2f7.2)
 1061 format(' step length for integration:',f6.3)
 1030 format(/, ' control parameters',/,' kout kout2 kout3')
 1040 format(1x,i3,i5,i6)
 1031 format(/,' neqs nsht nbls wtsht',
     2 ' nitloc wtsp  zmin eigtol  rmscut hitct',
     3 ' dvpmx dvpvsmx')
 1032 format(/,' idmp  vpdamp vpvsdmp stadamp ires nitmax snrmct ',
     2 'dxmax rderr ercof')
 1041 format(1x,3i4,4x,f4.1,1x,i3,f8.2,f6.2,f6.3,f9.3,f6.0,2f6.2)
 1042 format(i3,f10.2,2f8.2,i3,i7,f9.5,3f6.2,
     * /,' # its. for 1-d vel. model = ',i3,/,' rms for term. = ',f5.3
     * ,/,' fix locations for iterations = ',i3
     * ,/,' distance weighting:',2f7.2,'; residual weighting:',3f5.2)
      write(16,1166) i3d,xfac,tlim,nitpb(1),nitpb(2)
 1166 format('  parameters for pseudo-bending:',/,
     * '  i3d xfac   tlim  nitpb(1) nitpb(2)',/,i4,f6.2,f7.4,i6,i10)
      if(iuses.gt.0) then
        if(iusep.eq.1) then
          write(16,1676) iusep,iuses
 1676     format(/,' Both P-velocity and Vp/Vs are included in ',
     2    'inversion  (iusep=',i2,', iuses=',i2,')')
        else
          write(16,1677) iusep,iuses
 1677     format(/,' Only S velocities are included in ',
     2    'inversion  (iusep=',i2,', iuses=',i2,')')
        endif
      else
       write(16,1678) iusep,iuses
 1678  format(/,' Only P velocities are included in ',
     2    'inversion  (iusep=',i2,', iuses=',i2,')')
c  put station damping in 2nd position if no Vs used
      vdamp(2)=vdamp(3)
      endif
c  increment iuses so it can be an index to arrays
      iuses=iuses+1
      if(invdel.ne.0)write(16,1080)
 1080 format(/,' Station P and S-P delays included in inversion',/)
      ddlt=1.0/(delt2-delt1)
      if(res3.gt.res2) then
        dres12=0.98/(res2-res1)
        dres23=0.02/(res3-res2)
      else
        dres12=1.0/(res2-res1)
        dres23=0.0
      endif
      nevt=neqs+nsht+nbls
      return
  999 write(16,1699)
 1699 format('error in input1, probably forgot kout3')
c***** end of subroutine input1 *****
      stop
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine input2
c  this routine reads in the station list, sets up the
c  coordinate system, and calculates the stations' cartesian coordinates.
c  (reads from file02 )
c  subroutines required: setorg; disto;
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,xlt,xln,snr,csr
      include 'simulps_common.inc'
c
      write(16,2000)
 2000 format(/,'  origin :  latitude   longitude   rotation')
c  read in center of coordinates and rotation angle
      read(2,*) ltdo,oltm,lndo,olnm,rota
c  If inverting station corrections, setup output file22 to
c  be used as input file in future runs.
      if((invdel.eq.0).or.(kout.lt.2)) goto 5
      open(unit=22,file='fort.22',status='new')
	rewind 22
      write(22,2001) ltdo,oltm,lndo,olnm,rota
 2001 format(i3,1x,f5.2,i4,1x,f5.2,f8.2)
    5 write(16,2002) ltdo,oltm,lndo,olnm,rota
 2002 format(12x,i3,1x,f5.2,2x,i4,1x,f5.2,3x,f7.2)
c  set up short-distance conversion factors, given center of coords
      call setorg(ltdo,oltm,lndo,olnm)
      write(16,2004)
 2004 format('   station   latitude   longitude   elev',
     2 '     x      y      z   pdl s-pdl nfixst')
c  read in number of stations
      read(2,*) nsts
        if(nsts.le.maxsta)go to 40
        write(16,41)
 41     format('0Too many stations for input arrays.')
        stop
 40     continue
c  read in station list
      do 20 j=1,nsts
         read(2,2007,end=99) stn(j),ltds(j),sltm(j),lnds(j),
     *   slnm(j),ielev,pdl(j),sdl(j),nfixst(j)
 2007    format(2x,a4,i2,1x,f5.2,i4,1x,f5.2,i5,2f5.2,i3)
         z=-ielev*1.0e-3
c  calculate cartesian coordinates of station
         call disto(1,j,x,y)
         if((kout2.lt.4).or.(j.eq.1)) 
     2      write(16,2009) j,stn(j),ltds(j),sltm(j),lnds(j),
     3      slnm(j),ielev,x,y,z,pdl(j),sdl(j),nfixst(j)
 2009    format(1X,i4,3x,a4,1x,i3,1x,f5.2,2x,i4,1x,f5.2,2x,i5,
     *      1x,3f7.2,2f5.2,i3)
c        call latlon(x,y,latsta,xltsta,lonsta,xlnsta)
c        write(16,1666) latsta,xltsta,lonsta,xlnsta
c1666    format(2x,'checking latlon ',2(i4,f6.2))
c  store station coordinates
         stc(1,j)=x
         stc(2,j)=y
         stc(3,j)=z
   20 continue
      if(kout2.lt.4) return
c
c  kout2=4, condensed printout
   50 continue
      write(16,1605) nsts
 1605 format('  number of stations read in Input2 =',i5,
     2 ' (print 1st and last only)')
      j=nsts
      write(16,2009) j,stn(j),ltds(j),sltm(j),lnds(j),slnm(j),
     2 ielev,x,y,z,pdl(j),sdl(j),nfixst(j)
      return
c
c  not enough stations in list - error
   99 nsts=j-1
      write(16,2099) nsts
 2099 format('/ *** too few stations in list *** continue with ',i5)
      if(kout2.eq.4) goto 50
      return
c***** end of subroutine input2 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine input3
c  this routine reads in the initial velocity model in the
c  form of velocity specified on a uniform but not
c  necessarily evenly spaced grid of points
c  (reads from file03 )
c
c  common block variables:
      include 'simulps_common.inc'
c
c  declaration statements:
      integer ixf(maxpar),iyf(maxpar),izf(maxpar)
      character*1 vtype(2)
      parameter(zero=0.0,izero=0)
c
      vtype(1)='P'
      vtype(2)='S'
c
c  for this version the gridpoints can be unevenly spaced
c  the origin of the coordinate system is at (x,y,z)=(0,0,0)
c  which will not in general correspond to the point
c  xn(1),yn(1),zn(1).
c  xn,yn,zn should be factors of bld (ie: a.0 for bld=1.0 or a.b for bld=0.1)
c
c input the number of gridpoints in x, y and z directions
c  and bld factor (1.0 or 0.1 km) used to set up velocity interpolation grid
      read(3,3002) bld,nx,ny,nz
 3002 format(f4.1,3i3)
      if((bld.ne.1.0).and.(bld.ne.0.1)) then
        write(16,1625) bld
 1625   format(/, '******** STOP *********, bld must be 1.0 or 0.1,
     2   not ',f6.2)
      endif
        atemp=iuses*(nx-2)*(ny-2)*(nz-2)
        if(atemp.le.maxpar)goto 40
        write(16,42)
 42     format('0Too many nodes for program array sizes.')
        stop
 40     continue
c
c  input the x grid, y grid, and z grid
        read(3,3004) (xn(i),i=1,nx)
        read(3,3004) (yn(i),i=1,ny)
        read(3,3004) (zn(i),i=1,nz)
 3003 format(3i3)
 3004 format(20f6.1)
c
      write(16,3005) bld,nx,ny,nz
 3005 format(//,' velocity grid size:',/,
     * 'bld =',f4.1,5x,' nx =',i3,5x,'ny =',i3,5x,'nz =',i3)
c
      write(16,3006) (xn(i),i=1,nx)
 3006 format(/,' xgrid',/,3x,20f6.1)
      write(16,3007) (yn(i),i=1,ny)
 3007 format(/,' ygrid',/,3x,20f6.1)
      write(16,3008) (zn(i),i=1,nz)
 3008 format(/,' zgrid',/,3x,20f6.1,/)
c
c  read in which nodes to have fixed velocity
c  end with blank line
      i=1
   50 read(3,3003) ixf(i),iyf(i),izf(i)
      if(ixf(i).le.0) goto 60
      i=i+1
      goto 50
   60 continue
      inf=i-1
c
c  now read in the velocity values
   65 write(16,3101)
c     do 38 kv=1,iuses
         kv=1
         do 37 k=1,nz
            k2=k + (kv-1)*nz
            write(16,3015) k,vtype(kv),zn(k)
            do 36 j=1,ny
               read(3,3011) (vel(i,j,k2),i=1,nx)
               write(16,3013) (vel(i,j,k2),i=1,nx)
   36       continue
   37    continue
c  38 continue
c CHANGE FOR VP/VS INVERSION
      if(iuses.eq.2) then
        do 100 k=1,nz
           write(16,3016) k,zn(k)
           do 99 j=1,ny
              read(3,3011) (vpvs(i,j,k),i=1,nx)
              write(16,3013) (vpvs(i,j,k),i=1,nx)
   99     continue
  100  continue
c  compute Vs from Vp and Vp/Vs
        kv=2
        do 120 k=1,nz
           ks=k+nz
           write(16,3015) k,vtype(kv),zn(k)  
           do 115 j=1,ny
              do 110 i=1,nx
                 vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
  110         continue
              write(16,3013) (vel(i,j,ks),i=1,nx)
  115      continue
  120   continue
      endif
c
 3013 format(20f6.2)
 3015 format(/,' layer',i3,5x,a1,' velocity',10x,'z =',f7.1)
 3016 format(/,' layer',i3,5x,'Vp/Vs',10x,'z =',f7.1)
 3011 format(20f5.2)
 3101 format(//,' velocity values on three-dimensional grid')
c  compute total number of gridpoints (nodes)
      nodes=nx*ny*nz
      nxy=nx*ny
      nx2=nx-2             ! number non-edge nodes in row
      nxy2=nx2*(ny-2)      ! number non-edge nodes in layer
      nz2=nz-2
      nodes2=nz2*nxy2
c  peripheral nodes
      nx1=nx-1
      ny1=ny-1
      nz1=nz-1
c  Number of medium parameters to invert for
      npar=nodes2*iuses
      nparv=npar
      if(invdel.ne.0)npar=(npar+nsts*iuses)
c
c  Check to see whether medium parameters fit within array sizes
      if(nparv.gt.maxpar) goto 980
c
c  fix specified nodes by setting nfix(k)=1, else=0
      if(inf.eq.0) goto 496
        do 70 i=1,inf
           iizf=izf(i)-2
c  if s velocity node
           if(izf(i).gt.nz) iizf=izf(i)-4
           k=iizf*nxy2 + (iyf(i)-2)*nx2 + (ixf(i)-1)
           nfix(k)=1
   70   continue
c
       write(16,1610)
 1610 format(/,' velocity FIXED at the following nodes(1):')
  311 do 495 kv=1,iuses
         nz1=nz-1
         ny1=ny-1
         do 320 k=2,nz1
            if(kv.eq.1) write(16,1009) k,vtype(kv),zn(k)
 1009       format(/,' layer',i3,5x,a1,'-velocity nodes',
     2         10x,'z =',f7.1)
            if(kv.eq.2) write(16,3016) k,zn(k)
            kk=k+(kv-1)*nz2
            do 310 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1005) (nfix(i),i=n1,n2)
  310       continue
  320    continue
 1005 format('    ',18i6)
  495 continue
  496 continue
c
c  ndexfx: index from full nodes to nodes reduced by fixed (invert nodes)
c  mdexfx: index from inversion solution nodes to full velocity nodes
      in=0
      do 80 i=1,nparv
         if(nfix(i).eq.1) goto 80
         in=in+1
         ndexfx(i)=in
         mdexfx(in)=i
   80 continue
      inf2=nparv-in
      if(inf2.eq.inf) goto 85
        write(16,1615) inf,nparv,in,inf2
 1615   format(/,' **** number of fixed nodes input,',i4,
     2  ', does not equal velocity nodes,',i4,', minus invert',
     3  ' nodes,',i4,'.  Continue with inf=',i5,' ****',/)
        inf=inf2
   85 continue
      nparvi=nparv-inf
      npari=npar-inf
      if(invdel.eq.0) goto 95
c  also set indices if station delays are included in inversion
      i1=nparv+1
      do 90 i=i1,npar
         is=i-nparv
c s-delay
         if(is.gt.nsts) is=is-nsts
         if(nfixst(is).eq.1) goto 90
         in=in+1
         ndexfx(i)=in
         mdexfx(in)=i
   90 continue
      npari=in
   95 continue
      write(16,1620) npar,nparv,npari,nparvi
 1620 format(' INPUT3:npar,nparv,npari,nparvi',4i6)
c  Check to see whether medium parameters fit within array sizes
      if(npari.gt.mxpari) goto 990
      if(npar.gt.10000) goto 995
c
c  Set up an array which is used to point to node indices, for any x,y,z
      call bldmap
c
      return
c
  980 continue
      write(16,1698) nparv,maxpar
 1698 format(/,'  ****** STOP ******',/,i8,' velocity nodes, program',
     2 ' arrays only allow',i6)
      stop
  990 continue
      write(16,1699) npari,mxpari
 1699 format(/,'  ****** STOP ******',/,i8,' parameters to invert for',
     2 ', program arrays only allow',i6)
      stop
  995 continue
      write(16,1695) npar
 1695 format(/,'  ****** STOP ******',/,i8,' parameters, arrays are ',
     2 'for 10000')
      stop
c
c***** end of subroutine input3 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine input4(n,nfile,rmk)
c  this routine reads in the p travel times for all stations observing
c  the current event.  first trial hypocentral parameters are read, then
c  station names and weights and travel times.  arrivals at stations not
c  in the station list are discarded; unnormalized weights are calculated.
c  (reads earthquakes from file04, reads shots from file07 )
c  subroutine required: disto;
c
c  declaration statements:
c  local variables:
      real dep,x,y,pwt
      dimension nwav(2),ip(6),tt(6),sta(6)
      character*85 line
c
c  common block variables:
      include 'simulps_common.inc'

      character*1 rmk(maxobs, maxev), rmki(6), is(6)
c
      data  blank,iblank,blank0 /4h    ,4h    ,4h0   /
c
      if (n.gt.1) go to 1
c  print out heading for event cards
      write(16,4000)
 4000 format(//,' trial event locations:',/,
     2 5x,'n    origin time     latitude longitude  depth',
     3 '   mag',5x,'x     y     z',4x,'nob np ns')
      nobt=0
      nobtp=0
      nobts=0
    1 nobs=0
      nwav(1)=0
      nwav(2)=0
      wsum=0
c  check for extra travel time card
c  read in event card


c   10 read(nfile,4111,end=99) nchar,line,chk
   10  read(nfile, 4111, end=99)line  
 4111  format(a85)
    5  if(line(1:2).eq.' '.or.line(1:2).eq.'0 ')goto 10 

c 4111 format(q,a85,t1,a4)
c     write(16,4012) nchar,line,chk
4012   format(2x,i5,' characters in line',/,2x,a85,/,2x,'first 4 ',
     2 'char = ',a4)
c    5 if((chk.eq.blank).or.(chk.eq.blank0)) goto 10


   15 read(line,4001,err=5098) iyrmo(n),iday(n),ihr(n),
     2 mino(n),seco(n),ltde(n),eltm(n),lnde(n),elnm(n),dep,rmag(n)

 
4001  format(a4,a2,1x,a2,i2,1x,f5.2,i3,1x,f5.2,1x,i3,1x,f5.2,2f7.2)
      if (iyrmo(n).eq.iblank) go to 10
c  calculate event location in cartesian coordinate system
      call disto(2,n,x,y)
c  store event coordinates
      evc(1,n)=x
      evc(2,n)=y
      evc(3,n)=dep
c  read in travel times to stations in groups of six
   20 continue
c  format including phase-remarks
      read(nfile,4007) line
 4007 format(a85)
      read(line,4006,err=22) (sta(j),rmki(j),tt(j),j=1,6),
     2  (is(j),ip(j),j=1,6)
 4006 format(6(2a4,f6.2),t1,6(5x,a1,1x,i1,6x))
      goto 25
c Alberto Michelini (LBL) format
   22 read(line,4008,err=23) (sta(j),rmki(j),tt(j),j=1,5),
     2  (is(j),ip(j),j=1,5)
 4008 format(5(2a4,f7.4),t1,5(5x,a1,1x,i1,7x))
      goto 25
c  old format
   23 read(line,4005,err=5099) (sta(j),is(j),ip(j),tt(j),j=1,6)
 4005 format(6(a4,a1,i1,f6.2))
c
c  loop over the six readings
c  terminate if station is a blank (end of observations for event)
c  skip if weight is incorrect in some way
c  search list for current station, and index it if in list
c  also calculate p-arrival time and unnormalized weight
c  and increment number of observations for event
c
   25 continue
c  blank line separates the events
      if ((sta(1).eq.blank).or.(sta(1).eq.blank0)) go to 50
      do 30 j=1,6
         if(ip(j).ge.4) goto 30
         pwt=1.0/(2.**(ip(j)))
         if (pwt.le.0.0) go to 30
c search station list
         do 33 k=1,nsts
            if (sta(j).eq.stn(k)) go to 35
   33    continue
         if(kout2.gt.1) goto 30
         if(sta(j).eq.blank .or. sta(j).eq.blank0) then
            write(16,1609)
 1609       format(' Blank-station observation ignored')
         else
            write(16,1610) sta(j),iyrmo(n),iday(n),ihr(n),
     *      mino(n),seco(n)
 1610       format(' ** WARNING:  Observed station "',a4,
     *      '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *      ') is not in station list--observation ignored **')
         endif
         go to 30
c  add to set of observing stations if on list
   35    continue
c  don't use arrivals further away than delt2
         dx=evc(1,n)-stc(1,k)
         dy=evc(2,n)-stc(2,k)
         dz=evc(3,n)-stc(3,k)
         delta=sqrt((dx*dx)+(dy*dy)+(dz*dz))
         if(delta.lt.delt2) goto 36
c        write(16,1602) delta,stn(k),n
 1602    format(' **** delta=',f7.2,2x,a4,' not used for event',i4)
         goto 30
c        set iuses to 0 (changed to iuses+1 in input1) to disable
c           S-wave processing.
   36    if(iuses.eq.1.and.is(j).eq.'S') goto 30
         if((iusep.eq.0).and.((is(j).eq.'P').or.(is(j).eq.'p')))
     *      goto 30
         nobs=nobs+1
c check for too many observations 
         if(nobs.gt.maxobs) then
            write(6,1603) iyrmo(n),iday(n),ihr(n),mino(n),maxobs
 1603       format('  INPUT ERROR, Too many readings for event:',
     2      a4,a2,1x,a2,i2,/,'    CONTINUE with first',i5,
     3      'observations')
            write(16,1603) iyrmo(n),iday(n),ihr(n),mino(n),maxobs
            nobs=maxobs
            goto 50
         endif
         dlta(nobs,n)=delta
         iw(nobs,n)=ip(j)
         isto(nobs,n)=k
         secp(nobs,n)=tt(j)+seco(n)
         rmk(nobs,n)=rmki(j)
c  here, except for the addition of pdl(k) in the above line
c  are additions by wap to read in s times. If the time is
c  an s time, is(j) will equal 'S'. If this is true, the
c  proper adjustment will be made for the station delay, and
c  a 1 will be put in intsp(nobs,n).
         intsp(nobs,n)=0
         if(rmk(nobs,n).eq.' ') rmk(nobs,n)=' P  '
         if((is(j).ne.'S').and.(is(j).ne.'s')) goto 37
c  intsp(nobs,n) tells whether the nobs'th observation of the
c  arrival time is a P (0), or an S (1).
         intsp(nobs,n)=1
c** S-P CHANGE
c   For S-P we do not want to add origin time
         secp(nobs,n)=tt(j)
         if(rmk(nobs,n).eq.' ') rmk(nobs,n)=' S  '
  37     continue
         kv=intsp(nobs,n)+1
         nrd(k,kv)=nrd(k,kv)+1
         nwav(kv)=nwav(kv)+1
         wt(nobs,n)=pwt
c  19:july-1983 Change in normalizing factor
         wsum=wsum+pwt
c        wsum=wsum+pwt*pwt
   30 continue
c  continue reading stations and travel times
      go to 20
   50 continue
c  end of travel time readings for current event
c  check that there are at least four readings for this event
c  (does not matter for shots)
      if(nobs.eq.0) goto 90
      if((nobs.lt.4).and.(n.le.neb)) go to 90
      write(16,4004) n,iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2 ltde(n),eltm(n),lnde(n),elnm(n),dep,rmag(n),x,y,dep,
     3 nobs,nwav(1),nwav(2)
 4004 format(3h **,1x,i3,1x,a4,a2,1x,a2,i2,f6.2,i3,'n',f5.2,i4,'w',f5.2,
     2 2f7.2,3f6.2,3x,3i3)
c  normalize reading weights
c   change in normalizing factor 19-jul-1983
      wfac=nobs/wsum
c     wfac=sqrt(nobs/wsum)
c  check for too many readings
      if(nobs.le.maxobs) goto 55
        write(16,1698) n,nobs,maxobs,maxobs
 1698   format(' **** ERROR **** in event ',i4,', TOO MANY ',
     2  'OBSERVATIONS, nobs=',i6,' array size allows',i6,/,
     2  '   continuing with',i6,' observations',/)
        nobs=maxobs
   55 kobs(n)=nobs
      nobt=nobt+nobs
      nobtp=nobtp+nwav(1)
      nobts=nobts+nwav(2)
      if(n.gt.netemp) then
        nobtex=nobtex+nobs
        nobt=nobt+(wtsht-1)*nobs
      else
        nobteq=nobteq+nobs
      endif
      do 60 j=1,nobs
      wt(j,n)=wt(j,n)*wfac
   60 continue
      return
c  not enough readings for event n
   90 continue
      write(16,4009)  iyrmo(n),iday(n),ihr(n),mino(n)
 4009 format(' *** event ',a4,a2,1x,a2,i2,' should be discarded - ',
     2 'too few observations ***',/,' *** SKIPPING TO NEXT EVENT ***')
      if(nfile.eq.7) then       ! shot datafile
        nsht=nsht-1
      else
        if(nfile.eq.4) then     !earthquake datafile
          if(ifixl.le.0) then
            neqs=neqs-1
          else
            nbls=nbls-1
          endif
          netemp=netemp-1
        else                    ! blast datafile
          nbls=nbls-1
          nbtemp=nbtemp-1
        endif
        neb=neb-1
      endif
      nevt=nevt-1
      goto 1
c  ran out of event cards - error]
   99 write(16,4099) nfile
 4099 format('**** end of data file',i2,' - neqs too large ****')
      if(nfile.ne.4) goto 199
c  assume input error if neqs was less than 10% off, and continue with smaller neqs
c  assume other error,possibly wrong file, if more than 10% off, then stop
        if((abs(netemp-n)).gt.(0.10*netemp)) goto 199
        netemp=n-1
        neb=nbtemp+netemp
        nevt=neb+nsht
          if(ifixl.le.0) then
            neqs=n-1
          else
            nbls=neb
          endif
      write(16,1699) netemp
 1699 format(' continue with neqs= ',i4)
      return
  199 continue
      write(16,1690)
 1690 format(' stop since neqs more than 10% off, or problem '
     2 ,'was nbls or nsht')
      stop
 5098 write(16,5050) n
      write(16,4012) nchar,line,chk
      write(16,4011) iyrmo(n),iday(n),ihr(n),
     2 mino(n),seco(n),ltde(n),eltm(n),lnde(n),elnm(n),dep,rmag(n)
 4011 format(1x,a4,a2,a3,i2,1x,f5.2,i3,1x,f5.2,1x,i3,1x,f5.2,2f7.2)
      goto 10
c     stop
 5099   write(16,5050) n
      write(16,4015) (sta(j),is(j),ip(j),tt(j),j=1,6)
 4015 format(1x,6(a4,a1,i1,f6.2))
 5050   format('0error in data in input4. event= ',i8)
        stop
c****** end of subroutine input4 ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine intmap(x,y,z,ip,jp,kp)
c  Modified by W. Prothero so a single call can get the indices
c  common block variables:
      include 'simulps_common.inc'
c
      ip=int((x+xl)/bld)
      ip=ixloc(ip)
      jp=int((yl+y)/bld)
      jp=iyloc(jp)
      kp=int((z+zl)/bld)
      kp=izloc(kp)
c  If an array element=0, the position is off the map.
      return
c***** end of subroutine intmap *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
       subroutine latlon(x,y,lat,xlat,lon,xlon)
c  Subroutine to convert from Cartesian coordinates back to
c  latitude and longitude.
c
        common/shortd/xltkm,xlnkm,rota,xlt,xln,snr,csr
c
      rad=1.7453292e-2
      rlt=9.9330647e-1
c
      fy=csr*y-snr*x
      fx=snr*y+csr*x
c
      fy=fy/xltkm
      plt=xlt+fy
c
      xlt1=atan(rlt*tan(rad*(plt+xlt)/120.))
      fx=fx/(xlnkm*cos(xlt1))
      pln=xln+fx
c
      lat=plt/60.
      xlat=plt-lat*60.
      lon=pln/60.
      xlon=pln-lon*60.
c
      return
c
c***** end of subroutine latlon ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine loceqk(ne,niti,nwr)
c  this routine locates the given event in the initial velocity
c  structure using approximate ray tracing
c
c  common block variables:
      common/shortd/xltkm,xlnkm,rota,xlt,xln,snr,csr
      include 'simulps_common.inc'
c
c  declaration statements:
      real dr5(3),dr5i(3),s(4),v(4,4),adj(4),ster(4),cov(4)
      real w(maxobs)
      double precision ah(maxobs,maxobs),ac(3,3)
c  variables for error ellipse calculations
      integer iaz(3),idip(3)
      real covar(4,4),v3(4,4),serr(3)
      parameter (drad=1.7453292d-02)
c
c  Define small change in hypocenter, depending on bld factor
      small=0.5
      if(bld.eq.0.1) small=0.1
c
      xmov=900.0
c  nitl=counter for location iterations
      nitl=0
c  niti=counter for inversion iterations, Thurber starts with niti=1,
c       Eberhart-Phillips starts with niti=0
      if(niti.eq.0) then
        jfl=0
      else
        jfl=2
        if(ihomo.eq.niti) jfl=1
      end if
      jflag=jfl
      if (nitloc.eq.0) go to 9999
      rmswt=0.0
      rmslst=10.0
      rmsbfl=20.0
      nobs=kobs(ne)
c  location iteration loop
    1 continue
c  event coordinates
      xe=evc(1,ne)
      ye=evc(2,ne)
      ze=evc(3,ne)
c  loop over all observations of this event
      do 10 no=1,nobs
c  check for P or S reading 
      isp=intsp(no,ne)
      ns=isto(no,ne)
      xr=stc(1,ns)
      yr=stc(2,ns)
      zr=stc(3,ns)
c  If small change in hypocenter location, start with
c  old pb 3D path for additional hypocenter only iterations
      if(i3d.lt.2) goto 120
      if(pl(no).eq.0.0) pl(no)=dlta(no,ne)
c     if(xmov.gt.(0.05*pl(no))) goto 120
c  small change is defined at beginning of loceqk
      if(xmov.gt.small) goto 120
c  Move end of raypath
      dr5(1)=xe-rp(1,1,no)
      dr5(2)=ye-rp(2,1,no)
      dr5(3)=ze-rp(3,1,no)
      do 105 k=1,3
        dr5i(k)=dr5(k)/5.0
  105 continue
      do 115 j=1,5
        dj=float(5-(j-1))
        do 110 k=1,3
          rp(k,j,no)=rp(k,j,no)+dr5i(k)*dj
  110   continue
  115 continue
      ttime=ttc(no)
      goto 125
c  determine 3-d path using circular raypaths.
c  determine approximate path
  120 ncrold=nco(no,ne)
      ndpold=ndo(no,ne)
      call rayweb(no,isp,xe,ye,ze,xr,yr,zr,
     * fstime,jflag,ncrold,ndpold)
      nco(no,ne)=ncrold
      ndo(no,ne)=ndpold
c
c  calculate residual and hypocentral partial derivatives
      ttime=fstime
c
      if (i3d.lt.2) go to 12
  125 continue
c  number of pb iter depends on distance
      nitpbu=nitpb(1)
      if(dlta(no,ne).gt.delt1) nitpbu=nitpb(2)
      call minima(no,isp,ttime,nitpbu,jpb)
      if(jpb.lt.nitpbu) goto 130
        write(26,2601) no,stn(ns),jpb,nitpbu
 2601   format(' Minima: no=',i4,2x,a4,', used maximum ',
     2  'number PB iter.: j=',i3,', nitpb=',i3)
  130 continue
   12 continue
c
      call ttmder(ne,no,0,ttime)
c        write(16,305)ns,ttime
c305    format(' for station ',i8,' the ttime= ',f7.3)
      ttc(no)=ttime
      if((niti.lt.nitmax).or.(kout.lt.3)) goto 10
        x1=rp(1,1,no)
        x2=rp(1,2,no)
        y1=rp(2,1,no)
        y2=rp(2,2,no)
        z1=rp(3,1,no)
        z2=rp(3,2,no)
        call aztoa(x1,x2,y1,y2,z1,z2,xr,yr,azim,tkofan)
        az(no)=azim
        toa(no)=tkofan
   10 continue
c
      jfl=2
      jflag=2
c  calculate weights and apply to hypocenter matrix;
c  calculate rms residual
c  place derivs + resids into a-matrix for inversion
      call wthyp(ne,nobs,ah,rmswt,nwr,w)
      if((kout2.eq.0).or.(kout2.eq.2)) write(16,1001) nitl,xe,ye,ze,
     2 seco(ne),rmswt
 1001 format('  location, o.t., rmswt on iteration',i3,' are -',3f7.2,
     2 ',',f7.2,',',f6.3)
      if(nwr.lt.4) return
c  terminate or continue iterating?
c  check for worsened location, over 2 iterations
      if(rmswt.gt.rmslst) then
        if(rmswt.gt.rmsbfl) then
          evc(1,ne)=xbflst
          evc(2,ne)=ybflst
          evc(3,ne)=zbflst
          seco(ne)=otbflt
          rmswt=rmsbfl
          goto 999
        endif
      endif
      if (nitl.ge.nitloc) go to 999
c ************* wap *****************************
c  if rapid convergence, continue despite rmscut
c     write(16,1666) rmswt,rmslst,rmscut
 1666 format(' rmswt,rmslst,rmscut ',3f6.3)
      if((rmslst-rmswt).ge.0.02) goto 15
c  quit if 2 iterations better than rmscut
      if (rmswt.le.rmscut) then
        if(rmslst.le.rmscut) then
          if(rmsbfl.le.rmscut) goto 999
        endif
      endif
c  WP added this to end iteration when convergence in reached.
c  GE terminate when no improvement of rms in 2 iterations
      if((rmslst-rmswt).lt.0.005) then
        if((rmsbfl-rmslst).lt.0.005) goto 999
      endif
      if(ne.gt.neqs) goto 15
      if(xmov.lt.0.05) goto 999
c  find new hypocenter
   15 continue
      xbflst=xlst
      ybflst=ylst
      zbflst=zlst
      otbflt=otlst
      rmsbfl=rmslst
      xlst=evc(1,ne)
      ylst=evc(2,ne)
      zlst=evc(3,ne)
      otlst=seco(ne)
      rmslst=rmswt
      nitl=nitl+1
c  determine singular value decomposition of hypocenter matrix
      call fksvd(ah,s,v,maxobs,4,nobs,4,1,.false.,.true.)
c  calculate adjustments and apply to hypocentral parameters
c  determine number of non-zero singular values
      nfre=0
      do 20 i=1,4
      if (s(i).gt.eigtol) nfre=nfre+1
   20 continue
c  calculate adjustments
c   only adjust origin time for blast
      if(ne.gt.neqs) then
c    calculate adjustment
        adj(1)=0.0
        do 25 j=1,nfre
        adj(1)=adj(1)+v(1,j)*ah(j,5)/s(j)
   25   continue
c    apply adjustment
        seco(ne)=seco(ne)+adj(1)
c  calculate adjustments for earthquakes
      else
        do 30 i=1,4
           adj(i)=0.0
           do 29 j=1,nfre
              adj(i)=adj(i)+v(i,j)*ah(j,5)/s(j)
   29      continue
   30   continue
c  check for poor depth constraint
        adj(4)=adj(4)*0.75    !Cliff uses 0.75, Bill uses 0.50
        if (s(4).le.eigtol) adj(4)=0.0
        if(s(4).le.eigtol) write(16,9876)
 9876 format(' *** poorly constrained depth ***')
c  apply adjustments to hypocentral parameters
        seco(ne)=seco(ne)+adj(1)
        do 40 i=1,3
           i1=i+1
c  check to see that hyp. adjustment is less than dxmax
           if(adj(i1)) 32,38,34
   32      if(adj(i1).lt.(-dxmax)) adj(i1)= -dxmax
           goto 38
   34      if(adj(i1).gt.dxmax) adj(i1)=dxmax
   38      evc(i,ne)=evc(i,ne)+adj(i1)
   40   continue
c  check for depth less than zmin
      if(evc(3,ne).gt.zmin) goto 45
        write(16,1610) zmin,evc(3,ne)
 1610 format(5x,'***** Depth Fixed at ZMIN *****',' zmin, zadj= ',2f7.2)
        evc(3,ne)=zmin
   45 end if
      dx=xlst-evc(1,ne)
      dy=ylst-evc(2,ne)
      dz=zlst-evc(3,ne)
      dot=otlst-seco(ne)
      xmov=sqrt(dx*dx + dy*dy + dz*dz + dot*dot)
      go to 1
  999  continue
c  residual reduction achieved - proceed to parameter separation
c  if location only, calculate ssqr here instead of parsep
      if(nitmax.gt.0) goto 48
c  compute contribution to ssqr
      sqwtsh=sqrt(wtsht)
      totrms(ne)=0.0
      do 18 i=1,nobs
         resi=res(i)
         resi2=resi*resi
         totrms(ne)=totrms(ne)+resi2
         rrw=resi2*w(i)
         ssqrw=ssqrw+rrw
         wnobt=wnobt+w(i)
         if(intsp(i,ne).eq.0) then
c  P arrival
            ssqrwp=ssqrwp+rrw
            wnobtp=wnobtp+w(i)
         else
c  S arrival
            ssqrws=ssqrws+rrw
            wnobts=wnobts+w(i)
         endif
         if (ne.gt.netemp) resi=resi*sqwtsh
         ssqr=ssqr+resi*resi
   18 continue
      totrms(ne)=sqrt(totrms(ne)/nobs)
c
c  Covariance calculations for hypocenter parameters
c  Compute diagonal elements of covariance matrix
c  set reading error (sec).
   48 sig=rderr
      sig2=sig*sig
c  initially set st.er. to zero (so correct for blasts)
      do 47 i=1,4
         ster(i)=0.0
   47 continue
      ni=4
      if (ne.gt.neqs) ni=1
      do 50 i=1,ni
         sum=0.0
         do 60 j=1,4
            sss=s(j)
c fix s so won't get arithmetic fault
            if(sss.eq.0.0) sss=0.00000001
            sum=sum+v(i,j)*v(i,j)/(sss*sss)
   60    continue
         cov(i)=sig2*sum
         ster(i)=sqrt(cov(i))
   50 continue
      wrmsr(ne)=rmswt
      write(16,1099)ne,nitl,rmswt,seco(ne),(evc(j,ne),j=1,3),ster
 1099 format(2x,'event ',i3,'; nit=',i3,'; wtd rms res=',f6.3,
     2 '; ot,x,y,z =',4f6.2,
     3 '; St.Er.=',f8.4,'(ot)',f7.3,'(x)',f7.3,'(y)',
     4 f9.3,'(z)')
c  this is from Fred Klein's program HYPOINVERSE (6mar1989)
c  set Fred's "n" to 4
      n=4
C--SKIP THE REMAINING CALCS IF THEY ARE NOT NEEDED FOR FINAL OUTPUT
      if((kout.lt.5).or.(ne.gt.neqs)) goto 9999
C--CALCULATE COVARIANCE MATRIX AS SIGMA**2 * V * EIGVAL**-2 * VT
C--ESTIMATED ARRIVAL TIME ERROR
	SIGSQ=RDERR*rderr+ERCOF*rmswt*rmswt
	TEMP=EIGTOL**2
	DO 252 I=1,4
	DO 250 J=1,I
	COVAR(I,J)=0.
	IF (I.GT.N .OR. J.GT.N) THEN
	  IF (I.EQ.J) COVAR(I,J)=999.
	  GO TO 250
	END IF
	DO 245 L=1,N
c 245 COVAR(I,J)=COVAR(I,J)+V(I,L)*V(J,L)/(s(L)**2+TEMP)
245	COVAR(I,J)=COVAR(I,J)+V(I,L)*V(J,L)/(s(L)*s(L))
	COVAR(I,J)=SIGSQ*COVAR(I,J)
250	COVAR(J,I)=COVAR(I,J)
252	CONTINUE

C--EVALUATE THE HYPOCENTER ERROR ELLIPSE BY DIAGONALIZING
C--THE SPATIAL PART OF THE COVARIANCE MATRIX
C--USE ac AS TEMPORARY STORAGE
	DO 257 I=1,3
	DO 255 J=1,3
255	ac(I,J)=COVAR(I+1,J+1)
257	CONTINUE
	CALL fkSVD (ac,SERR,V3,MMAX,3,3,3,0,.FALSE.,.TRUE.)
	DO 260 I=1,3
	SERR(I)=SQRT(SERR(I))
	IF (SERR(I).GT.99.) SERR(I)=99.
260	CONTINUE

C--COMPUTE ERH AND ERZ AS THE LARGEST OF THE HORIZ AND VERTICAL
C--PROJECTIONS OF THE PRINCIPAL STANDARD ERRORS
	ERH=0.
	ERZ=0.
	DO 265 I=1,3
	TEMP=SERR(I)*SQRT(V3(1,I)**2+V3(2,I)**2)
	IF (TEMP.GT.ERH) ERH=TEMP
	TEMP=SERR(I)*ABS(V3(3,I))
	IF (TEMP.GT.ERZ) ERZ=TEMP
265	CONTINUE
	IF (ERZ.GT.99.) ERZ=99.
	IF (ERH.GT.99.) ERH=99.
c     write(16,1098)ne,nitl,rmswt,seco(ne),(evc(j,ne),j=1,3),
c    2 erh,erz
 1098 format(2x,'event ',i3,'; nit=',i3,'; wtd rms res=',f6.3,
     2 '; ot,x,y,z =',4f6.2,
     3 '; Hypinv.code St.Er.=',f10.3,' ERH',6x,f9.3,'ERZ')

C--NOW CALC THE ORIENTATIONS OF THE PRINCIPAL STD ERRORS
	DO 290 J=1,3
	IAZ(J)=0
	IDIP(J)=90
	TEMP=SQRT(V3(1,J)**2+V3(2,J)**2)
	IF (TEMP.EQ.0.) GO TO 290
c  Change indices because different coordinates. Hypoinverse array has 1=lat, 2=long(west=positive)
c     IAZ(J)=ATAN2(-V3(2,J),V3(1,J))/drad
c   Simul has evc(1)=west, evc(2)=north
c   and allows rotation
      iaz(j)=(atan2(-v3(1,j),v3(2,j))-rota)/drad
	IDIP(J)=ATAN2(V3(3,J),TEMP)/drad
	IF (IDIP(J).LT.0) THEN
	  IDIP(J)=-IDIP(J)
	  IAZ(J)=IAZ(J)+180
	END IF
	IF (IAZ(J).LT.0) IAZ(J)=IAZ(J)+360
        if(iaz(j).gt.360) iaz(j)=iaz(j)-360
290	CONTINUE
      call out(ne,niti,nwr,serr,iaz,idip,erh,erz)
 9999 continue
      return
c***** end of subroutine loceqk *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ludecp(a,ul,n,ier)
c  routine to perform cholesky decomposition
c
c  common block variables:
      include 'simulps_common.inc'
c
      dimension a(mgsol),ul(mgsol)
      data zero,one,four,sixtn,sixth/0.0,1.,4.,16.,.0625/
      d1=one
      d2=zero
      rn=one/(n*sixtn)
      ip=1
      ier=0
      do 45 i=1,n
         iq=ip
         ir=1
            do 40 j=1,i
            x=a(ip)
            if(j.eq.1) go to 10
               do 5 k=iq,ip1
               x=x-ul(k)*ul(ir)
               ir=ir+1
    5       continue
   10       if (i.ne.j) go to 30
            d1=d1*x
            if (a(ip)+x*rn.le.a(ip)) go to 50
   15       if (abs(d1).le.one) go to 20
            d1=d1*sixth
            d2=d2+four
            go to 15
   20       if (abs(d1).ge.sixth) go to 25
            d1=d1*sixtn
            d2=d2-four
            go to 20
   25       ul(ip)=one/sqrt(x)
            go to 35
   30       ul(ip)=x*ul(ir)
   35       ip1=ip
            ip=ip+1
            ir=ir+1
   40    continue
   45 continue
      go to 9005
   50 ier=129
 9000 continue
 9005 return
c* * * end of subroutine ludecp * * *
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine luelmp(a,b,n,x)
c  routine to perform elimination part of solution of ax=b
c
c  common block variables:
      include 'simulps_common.inc'
c
      dimension a(mgsol),b(mxpri1),x(mxpri1)
      data zero/0./
c  solution of ly=b
      ip=1
      kw=0
      do 15 i=1,n
         t=b(i)
         im1=i-1
         if (kw.eq.0) go to 9
         ip=ip+kw-1
         do 5 k=kw,im1
            t=t-a(ip)*x(k)
            ip=ip+1
    5    continue
         go to 10
    9    if (t.ne.zero) kw=i
         ip=ip+im1
   10    x(i)=t*a(ip)
         ip=ip+1
   15 continue
c  solution of ux=y
      n1=n+1
      do 30 i=1,n
         ii=n1-i
         ip=ip-1
         is=ip
         iq=ii+1
         t=x(ii)
         if (n.lt.iq) go to 25
         kk=n
         do 20 k=iq,n
            t=t-a(is)*x(kk)
            kk=kk-1
            is=is-kk
   20    continue
   25    x(ii)=t*a(is)
   30 continue
      return
c***** end of subroutine luelmp *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine medder(ne)
c routine to add medium derivatives from shot data
c  to the medium matrix dtmp for inversion
c
c  common block variables:
      include 'simulps_common.inc'
      common/wtpars/w(maxobs)
c
      nobs=kobs(ne)
      rms=0.0
c  compute weighting
c     as in parsep
c
      nwr=0
      wnorm=0.0
      do 13 j=1,nobs
c  reading weight
         w(j)=wt(j,ne)
c  residual weighting
c  downweighting(linear) 0 to 98% res1 to res2, 98 to 100% res2 to res3
         ares=abs(res(j))
         if(ares.le.res2) then
            wr=1.0-(ares-res1)*dres12
            if (wr.gt.1.0) wr=1.0
         else
            if(res3.gt.res2) then
               wr=0.02-(ares-res2)*dres23
               if (wr.lt.0.0) wr=0.0
            else
               wr=0.0
            endif
         endif
c  distance weighting
         wd=1.0-(dlta(j,ne)-delt1)*ddlt
         if (wd.gt.1.0) wd=1.0
         if (wd.lt.0.0) wd=0.0
c  unnormalized weight
         w(j)=w(j)*wr*wd
c  19-jul-1983 Change normalizing factor
         wnorm=wnorm+w(j)
c        wnorm=wnorm+w(j)*w(j)
         if (w(j).gt.0.0) nwr=nwr+1
  13  continue
c  check to be sure 4 or more readings with nonzero weight
c    dmep 02dec83
      if(nwr.ge.4) goto 12
      write(16,1612)
 1612 format(' ****** less than 4 readings with nonzero weights')
c  normalize weights
c    19-jul-1983 Change in normalizing factor
c     wfac=sqrt(nwr/wnorm)
   12 wfac=nwr/wnorm
      if(ne.gt.netemp) wfac=wfac*wtsht
      nswrt=nswrt+nwr
      do 14 j=1,nobs
         w(j)=w(j)*wfac
   14 continue
c  store residuals and compute rms residual
c  compute weighted rms residual as in wthyp
      wnorm=0.0
      rmswt=0.0
      do 10 ii=1,nobs
         resp(ii)= res(ii)*w(ii)
         w2=w(ii)*w(ii)
         rmswt= rmswt + w2*res(ii)* res(ii)
         wnorm=wnorm+w2
   10 continue
      rmswt=sqrt(rmswt/wnorm)
      wrmsr(ne)=rmswt
      write(16,1001) ne,rmswt
 1001 format(' explosion ',i3,'; weighted rms residual =',f7.3)
c  store partial derivatives
      nono=npari*nobs
      do 20 j=1,nono
         i=(j-1)/npari+1
         dtmp(j)=dtm(j)*w(i)
c check station corrections
c        im=j-(i-1)*npari
c        is=mdexfx(im)-nparv
c        if((is.le.0).or.(dtm(j).ne.1.0000)) goto 20
c        write(16,1620) is,j,dtm(j),dtmp(j)
c1620    format(' MEDDER:is,j,dtm,dtmp',i5,i10,2f14.6)
   20 continue
c  zero out dtm matrix
      do 45 m=1,nono
         dtm(m)=0.0
   45 continue
      return
c***** end of subroutine medder *****
      end
c
c-------------------------------------------------------
c
      subroutine minima(no,isp,fstime,npbmax,jpb)
c
c*****this routine finds the minimum path using pseudo-bending
c
c  common block variables:
      common/pathm/x(130),y(130),z(130),v(130),tra,n,nn
      common/temp/xtemp(130),ytemp(130),ztemp(130),rtemp(130),ttemp(130)
      common/pat/xt(130),yt(130),zt(130),rt(130),rm,rs,rtm,rts,nt,j
      common/upb/ upbtot,nupb
      include 'simulps_common.inc'
c
      tra=fstime
      n=nrp(no)
c
      do 20 i=1,n
         x(i)=rp(1,i,no)
         y(i)=rp(2,i,no)
         z(i)=rp(3,i,no)
   20 continue
c
      do 21 i=1,n
         xp=x(i)
         yp=y(i)
         zp=z(i)
         xtemp(i)=xp
         ytemp(i)=yp
         ztemp(i)=zp
         call vel3(isp,xp,yp,zp,vp)
         v(i)=vp
   21 continue
      call travel
c
      nn=n-1
      nupb=nupb+1
c     write(6,6000) nupb
 6000 format(' nubp=',i8)
      do 100 j=1,npbmax
         ta=tra
         call bend(isp,xfac)
         call travel
         deltat=ta-tra
c
         if (deltat.lt.0.0) go to 102
         do 22 i=1,n
            x(i)=xtemp(i)
            y(i)=ytemp(i)
            z(i)=ztemp(i)
   22    continue
         if(deltat.le.tlim) go to 102
100   continue
c
c
  102 continue


      if(j.gt.npbmax) j=npbmax
      upbtot=upbtot + float(j)
      jpb=j
  105 do 300 i=1,n
         rp(1,i,no)=x(i)
         rp(2,i,no)=y(i)
         rp(3,i,no)=z(i)
300   continue
      fstime=tra
c
      return
c ***** end of subroutine minima *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine out(ne,nit,nwr,serr,iaz,idip,erh,erz)
c
c  declaration statements:
      integer iaz(3),idip(3)
      real serr(3)
	DIMENSION KSIG(3)
	CHARACTER KSFL*1
      character*1 is,ie
      character*9 day
      character*8 tm
c
c  common block variables:
      include 'simulps_common.inc'
c
      call date_and_time(DATE=day)
      call date_and_time(TIME=tm)
c
c  Write new hypocenters to file12 in hypoinverse format
      if(ne.eq.1) write(12,1201) nit,day,tm
 1201 format(' iteration step=',i3,'  computed',1x,a9,1x,a8,
     1 ' hyp error ellipse',
     2 /,' Date HrMn Sec   Lat    Long DepthMgNwr',7x,'Rms',
     3 'Az1D1Ser1Az2D2Ser2Mg',3x,'Ser3',
     4 5x,'Erh Erz')
      i=ne
      call latlon(evc(1,i),evc(2,i),lat,xlat,lon,xlon)
c       convert origin time into gmt
        sec= seco(i)
        nin=mino(i)
 210    if(sec.lt.0) goto 230
 220    if(sec.lt.60) goto 240
        sec=sec-60
        nin=nin+1
        goto 220
  230   sec=sec+60
        nin=nin-1
        goto 210
  240 continue
c       write(12,260) iyrmo(i),iday(i),ihr(i),nin,sec,lat,xlat,lon,
c    * xlon,evc(3,i),rmag(i)
c260  format(a4,a2,1x,a2,i2,1x,f5.2,i3,1x,f5.2,1x,i3,1x,f5.2,2f7.2)
c
c  use Fred Klein's HYPOINVERSE format from his routine HYSUM
C--CALLED BY HYPOINVERSE TO OUTPUT SUMMARY DATA
c  set latitude to be north and longitude to be west
      is='n'
      ie='w'
c  use hypoinverse format
c  NOTE: I've left all the variables in the HYPOINVERSE output
c  for possible later use.
c  However those that are not calculated will be zero.
      ih71s=0
      iunit=12
c  Find minimum source-receiver distance
      do 5001 j=1,kobs(i)
      if(dlta(j,i).lt.dmin) dmin=dlta(j,i)
 5001 continue
C--CONVERT SOME DATA TO INTEGER FOR OUTPUT
	KLTM=xlat*100.+.5
	KLNM=xlon*100.+.5
	KQ=sec*100.+.5
	KZ=evc(3,i)*100.+.5
	LFMAG=rmag(i)*10.+.5
	LXMAG=rmag(i)*10.+.5
	KDMIN=DMIN+.5
	IF (KDMIN.GT.999) KDMIN=999
	KRMS=wrmsr(i)*100.+.5
	KERH=ERH*100.+.5
	KERZ=ERZ*100.+.5
c  nfm not used
c	NFM=NFRM
c	IF (NFM.GT.99) NFM=99
c use rmag(i) as input from for004 for both mxmag and mfmag
        mxmag=nint(10.0*rmag(i))
        mfmag=nint(10.0*rmag(i))
	IXMAG=.1*MXMAG+.5
	IF (IXMAG.GT.999) IXMAG=999
	IFMAG=.1*MFMAG+.5
	IF (IFMAG.GT.999) IFMAG=999
	DO 10 k=1,3
10	KSIG(k)=SERR(k)*100.
C--WRITE A SUMMARY RECORD
C--HYPO71 FORMAT
c  maxgap is not used although it could be calculated from az 
c  array (in AZTOA and LOCEQK)
	IF (IH71S.EQ.2 .AND. IUNIT.EQ.12) THEN
	  KSFL=' '
	  IF (NWS.GT.0) KSFL='S'
	  WRITE (12,1002) iyrmo(i),iday(i),ihr(i),nin,sec,LAT,IS,
c	2 xlat,LON,IE,xlon,evc(3,i),rmag(i),nwr,MAXGAP,DMIN,
	2 xlat,LON,IE,xlon,evc(3,i),rmag(i),nwr,DMIN,
     3 wrmsr(i),erh,erz,ksfl
1002	FORMAT (a4,a2,1x,a2,I2,F6.2,I3,A1,F5.2,I4,A1,F5.2,2F7.2,
c     2 I3,I4,F5.1,F5.2,2F5.1,1X,A1)
	2 I3,4x,F5.1,F5.2,2F5.1,1X,A1)
	ELSE

C--HYPOINVERSE FORMAT
	  WRITE (IUNIT,1001) iyrmo(i),iday(i),ihr(i),nin, KQ,LAT,IS,
c    2 KLTM,LON,IE,KLNM,KZ, LXMAG,nwr,MAXGAP, KDMIN,KRMS,
     2 KLTM,LON,IE,KLNM,KZ, LXMAG,nwr,KDMIN,KRMS,
c    3 (IAZ(k),IDIP(k),KSIG(k),k=1,2), LFMAG,REMK,KSIG(3),
     3 (IAZ(k),IDIP(k),KSIG(k),k=1,2), LFMAG,KSIG(3),
c    4 RMK1,RMK2,NWS, KERH,KERZ,NFM, IXMAG,IFMAG,KXMMAD,KFMMAD,
     4 KERH,KERZ, IXMAG,IFMAG
1001	  FORMAT (a4,a2,a2,i2, I4,I2,A1,
c     2 I4,I3,A1,I4,I5, I2,3I3,I4,
	2 I4,I3,A1,I4,I5, I2,i3,3x,i3,I4,
c    3 2(I3,I2,I4), I2,A3,I4,
	3 2(I3,I2,I4), I2,3x,I4,
c    4 2A1,I2, 2I4,I2, 4I3,
	4 2x,2x,2I4,2x, 2i3)
c    5 A3,A1, 3A1, I1,I3)
	END IF
c***** end of subroutine out *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outadj(nit,istop,istop1)
c
c  declaration statements:
      real fillin(20)
      real x(5000)
      parameter(zero=0.0,izero=0)
      character*5 vtype(2)
c
c  common block variables:
      include 'simulps_common.inc'
c
      vtype(1)='P-Vel'
      vtype(2)='Vp/Vs'
c
      write(16,1000) ssqr,nobt,nobtp,nobts,nobteq,nobtex,wtsht,
     2 nwrt,nswrt
       write(36,1000) ssqr,nobt,nobtp,nobts,nobteq,nobtex,wtsht,
     2 nwrt,nswrt
 1000 format(/,' sum of squared residuals =',f10.3,
     2 '; total number of observations =',i6,
     3 ', P obs =',i6,', S obs =',i6,
     4 /,53x,'earthquake obs.=',i6,', explosion obs.=',i6,
     5 ' (wtsht=',f5.2,')',
     6 /,20x,'with non-zero wt:',t65,'nwrt=',i6,',',10x,
     7 'nswrt=',i6)
      rmsres=sqrt(ssqr/float(nobt))
      rmssep=0.0
      do 88 ne=1,nevt
         rmssep=rmssep+seprms(ne)
   88 continue
      rmssep=sqrt(rmssep/nobt)
      write(16,1999) rmsres,rmssep
 1999 format(' rms residual =',f8.5,'   model rms =',f8.5)
      rmsw=sqrt(ssqrw/wnobt)
      write(16,1600) ssqrw,wnobt,rmsw
 1600 format(/,' weighted sum of squared residuals =',f10.3,
     * '; weighted total number of obs =',f10.3,
     * '; weighted rms res = ',f10.5)
c  store ssqr and number of velocity parameters for f-test
c  unless last step was backup
      if(istop1.ge.2) goto  5
      ssqrw1=ssqrw
      ssqr1=ssqr
      mbl1=mbl
      var1=var
      varw1=varw
      ndof1=ndof
      wndof1=wndof
    5 if(invdel.eq.0)go to 123
c=================================break================================
c  Output station corrections
      write(16,111)
  111 format(/,6x,'STATION CORRECTIONS AND ADJUSTMENTS',
     2 /,' station  Pdelay Adjusted S-Pdelay Adjusted')
      do 40 n=1,nsts
         nv=n+nparv
         nv1=nv+nsts
         if(hit(nv).lt.hitct.and.hit(nv1).lt.hitct) go to 40
         if(hit(nv).lt.hitct) vadj(nv)=0.
         if(hit(nv1).lt.hitct) vadj(nv1)=0.
c        write(16,120)stn(n),vadj(nv),vadj(nv1)
         write(16,130) stn(n),pdl(n),vadj(nv),sdl(n),vadj(nv1)
  130    format(1x,a4,2(f10.3,f9.3))
  120    format(3(1x,a4,1x,2f10.3))
   40 continue
c     write(16,110)
  110 format(/,' station corrections')
c     write(16,120)(stn(i),pdl(i),sdl(i),i=1,nsts)
 123  continue
c
      if(nparvi.eq.0) goto 999
c
      do 919 i=1,nx
         fillin(i)=0.0
  919 continue
      nx1=nx-1
      ny1=ny-1
      nz1=nz-1
c
      do 25 kv=1,iuses
         if((kv.eq.1).and.(iusep.eq.0)) goto 25
c
c        write(16,1688) nit
c1688    format(' outadj, nit=',i5)
         if(nit.gt.1) goto 488
         write(16,1633)
 1633    format(/,' DERIVATIVE WEIGHT SUM ')
         do 486 k=2,nz1
            k2=k + (kv-1)*nz
            write(16,1009) k,vtype(kv),zn(k)
            kk=k+(kv-1)*nz2
            do 485 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1635)(hit(i),i=n1,n2),zero
 1635          format(' 0.',18f7.0)
               do 484 i=n1,n2
                  sumhit(k2)=sumhit(k2)+hit(i)
  484          continue
  485       continue
  486    continue
  488    continue
c
         write(16,1001)
 1001    format(/,'  velocity model changes')
         do 12 k=2,nz1
            k2=k + (kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 12
            write(16,1009) k,vtype(kv),zn(k)
            write(16,6002) (fillin(i),i=1,nx)
            kk=k + (kv-1)*nz2
            do 10 j=2,ny1
               j1=(kk-2)*nxy2+(j-2)*nx2+1
               j2=j1+nx2-1
               write(16,1003) (vadj(i),i=j1,j2),fillin(1)
   10       continue
            write(16,6002) (fillin(i),i=1,nx)
   12    continue
 6002    format(20f6.2)
 1003    format('  0.00',20f6.2)
         write(16,2001)
 2001    format(/,' corrected velocity model',/)
         do 23 k=2,nz1
            k2=k + (kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 23
            write(16,1009) k,vtype(kv),zn(k)
 1009       format(/,' layer',i3,5x,a5,' nodes',10x,'z =',f7.1)
            do 22 j=1,ny
               if(kv.eq.1) write(16,6002) (vel(i,j,k2),i=1,nx)
               if(kv.eq.2) write(16,6002) (vpvs(i,j,k),i=1,nx)
  22        continue
  23     continue
c  Find average and variance (section taken from COMPVEL.FOR)
         write(16,3030)
 3030    format(/,'  STATISTICS, by Depth, excluding edge nodes')
         iv=0
         varv=0.0
         do 200 k=2,nz1
            write(16,3015) zn(k),vtype(kv)
            k2=k+(kv-1)*nz
            ix=0
            do 190 j=2,ny1
               do 180 i=2,nx1
                  ix=ix+1
                  if(kv.eq.1) then
                     x(ix)=vel(i,j,k2)
                  else
                     x(ix)=vpvs(i,j,k)
                  endif
  180          continue
  190       continue
            call avsd(0.00,x,ix,sd,av,devtot)
            iv=iv+ix 
            varv=varv+devtot
            write(16,3017) ix,av,sd
  200    continue
         varv=varv/iv
         sdv=sqrt(varv)
         write(16,3031)
 3031    format(/,' VARIANCE OF THIS VELOCITY MODEL')
         write(36,3018) iv,vtype(kv),varv,sdv
         write(16,3018) iv,vtype(kv),varv,sdv
 3015    format(/,' grid z=',f6.2,5x,a5)
 3017    format('  For ', i5,' gridpts, av=',f7.2,', sd=',f7.3)
 3018    format('  For all',i5,1x,a5,' gridpts, variance=',f12.8,
     2    ', sd=',f10.6)
c
   25 continue
      goto 999     !skip to make output shorter
c     if (nitloc.eq.0.or.neqs.eq.0) go to 999
c     write(16,1004)
c1004 format(/,' new hypocenters',/,4x,'number',5x,'ot',7x,'x',8x,
c    2 'y',8x,'z')
c     do 20 i=1,neqs
c        write(16,1008) I,SECo(i),(evc(j,i), j=1,3)
c1008    format(I9,4F9.2)
c  20 continue
  999 continue
c***** end of subroutine outadj *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outend(rmk)
c  Routine for final output.  See description of output files
c  and "kout" output control parameter at the beginning of
c  the program.
c
c  common block variables:
      include 'simulps_common.inc'
      common/machin/ eta,tol
        common/shortd/xltkm,xlnkm,rota,xlt,xln,snr,csr
      common/upb/ upbtot,nupb

      character*1 rmk(maxobs, maxev)
c
c  declaration statements:
      integer lat,lon
      real xlat,xlon
      real zeroi(maxnx)
      real type(2),ttime(maxobs),sekms(maxnx)
      data type(1),type(2) /1hP,1hS/
      parameter(zero=0.0,izero=0)
      character*50 formvl,formrs,line,frmrs1,frmrs2,frmvl1,frmvl2
      character*2 ch2nx2
      character*1 ch2nx3,ch2nx4
      character*5 vtype(2)
      character*9 day
      character*8 tm
c
      vtype(1)='P-Vel'
      vtype(2)='Vp/Vs'
      call date_and_time(DATE=day)
      call date_and_time(TIME=tm)
      formrs='(1x,3hres,5x,10(1h(,f5.2,1h),5x),3hres)'
      frmrs1='(1x,3hres,5x,10(1h(,f5.2,1h),5x))'
      frmrs2='(9x,2(12x),8(1h(,f5.2,1h),5x),3hres)'
      formvl='(1h0,f4.2,10(f6.2,1h+,f5.4),f6.2)'
      frmvl1='(1h0,f4.2,10(f6.2,1h+,f5.4))'
      frmvl2='(1h0,4x,2(12x),8(f6.2,1h+,f5.4),f6.2)'
      do 2 i=1,20
         zeroi(i)=zero
    2 continue
c
c
      write(16,100) ssqr,nobt,nobtp,nobts,nobteq,nobtex,wtsht,
     2 nwrt,nswrt
  100 format(/,' sum of squared residuals =',f10.3,
     2 '; total number of observations =',i6,
     3 ', P obs =',i6,', S obs =',i6,
     4 /,53x,'earthquake obs.=',i6,', explosion obs.=',i6,
     5 ' (wtsht=',f5.2,')',
     6 /,20x,'with non-zero wt:',t65,'nwrt=',i6,',',10x,
     7 'nswrt=',i6)
c** change for creating synthetic data, nitmax= -1
      if(nitmax.eq.-1) goto 76
      rmsres=sqrt(ssqr/float(nobt))
      write(16,1999) rmsres
 1999 format(' rms residual =',f8.5)
      rmsw=sqrt(ssqrw/wnobt)
      write(16,1600) ssqrw,wnobt,rmsw
 1600 format(/,' weighted sum of squared residuals =',f10.3,
     * '; weighted total number of obs =',f10.3,
     * '; weighted rms res = ',f10.5)
c
c  Write new hypocenters to output(file16)
   76 write(16,1650)
 1650 format(//,15x,'FINAL LOCATIONS',//,5x,'YRMODY HRMN  SEC',
     2 2x,'LATITUDE LONGITUDE  DEPTH    MAG  NO RMSRES',5x,'x',6x,
     3 'y',6x,'z')
c  Map 0,1,2,3,4,5 to -1,0,1,2,2,3
      kkout=nint(kout*0.85)-1
      if(kkout) 301,295,290
c  Write new travel time file24 to be used as input in future runs
  290 open(unit=24,file='fort.24',status='new')
      if(nbls.gt.0) open(unit=28,file='fort.28',status='new')
c  Write new hypocenters to file13 in hypo71 format
  295 open(unit=13,file='fort.13',status='new')
	rewind(13)
      write(13,1301) neqs,day,tm
 1301 format(i5,75x,'computed',2x,a9,1x,a8)
  301 neqbl=neqs+nbls
      iunit=24
      do 5 i=1,nevt
         if(i.gt.neqs) iunit=28
c  convert origin time into gmt
         sec= seco(i)
         nin=mino(i)
 210     if(sec.lt.0) goto 230
 220     if(sec.lt.60) goto 240
         sec=sec-60
         nin=nin+1
         goto 220
 230     sec=sec+60
         nin=nin-1
         goto 210
  240    call latlon(evc(1,i),evc(2,i),lat,xlat,lon,xlon)
         if(kkout) 304,303,302
c** Change for synthetic data, nitmax= -1
  302    if((i.gt.neqbl).and.(nitmax.ne.-1)) goto 303
         if(i.eq.1) then
         write(iunit,2404) iyrmo(i),iday(i),ihr(i),
     2      nin,sec,lat,xlat,lon,xlon,evc(3,i),rmag(i),
     3      day,tm
         else
         write(iunit,2401) iyrmo(i),iday(i),ihr(i),
     2      nin,sec,lat,xlat,lon,xlon,evc(3,i),rmag(i)
 2401    format(a4,a2,1x,a2,i2,1x,f5.2,i3,1x,f5.2,1x,i3,1x,f5.2,2f7.2)
 2404    format(a4,a2,1x,a2,i2,1x,f5.2,i3,1x,f5.2,1x,i3,1x,f5.2,2f7.2,
     2      10x,'computed',2x,a9,1x,a8)
         end if
c  write travel time observations for each event
         do 30 jobs=1,kobs(i)
            ttime(jobs)=secp(jobs,i)
            if(intsp(jobs,i).eq.0) ttime(jobs)=secp(jobs,i)-seco(i)
   30    continue
         write(iunit,2402) (stn(isto(j,i)),rmk(j,i),
     2      ttime(j),j=1,kobs(i))
         write(iunit,2403)
  303    write(13,1302) iyrmo(i),iday(i),ihr(i),nin,sec,
     2      lat,xlat,lon,xlon,evc(3,i),rmag(i),kobs(i),wrmsr(i)
  304    write(16,1602) i,iyrmo(i),iday(i),ihr(i),nin,sec,
     2      lat,xlat,lon,xlon,evc(3,i),rmag(i),kobs(i),wrmsr(i),
     3      (evc(j,i),j=1,3)
    5 continue
c
c Official version:
c2402 format(6(a4,a4,f6.2))
c Evans/Miller version:
 2402 format(6(a4,a4,f6.3))
c
 2403 format(4h0   )
 1302 format(a4,a2,1x,a2,i2,f6.2,1x,i2,'n',f5.2,1x,i3,'w',f5.2,f7.2,
     2 2x,f5.2,i3,9x,f5.2)
 1602 format(1x,i3,1x,a4,a2,1x,a2,i2,f6.2,1x,i2,'n',f5.2,1x,i3,'w',
     2 f5.2,f7.2,2x,f5.2,i4,f6.2,3f7.2)
    8 if(kkout) 310,306,305
  305 close(unit=24)
  306 close(unit=13)
c
  310 continue
c  output number of observations by station
      write(16,1640)
 1640 format(//,15x,'TALLY OF OBSERVATIONS',//,1x,
     2 10('station obs  '),/,10(7x,'P  S-P'),/,
     3 10(' ==== === ==='))
c     write(16,1645) (stn(i),(nrd(i,kv),kv=1,2),i=1,nsts)
c the above statement gave nonsense in the sun compiled version
      write(16,1645) (stn(i),nrd(i,1),nrd(i,2),i=1,nsts)
 1645 format(10(1x,a4,2i4))
c
c** Change for synthetic data, nitmax= -1
      if(nitmax.le.0) goto 900
      if(invdel.eq.0)go to 353
c       output final station corrections
        write(16,110)
 110  format(//,' FINAL STATION CORRECTIONS',/,' stn  typ Delay  StErr',
     2 ' Resol   nobs  stn  stn grdpt  inv soln',/,29x,
     3         'nrd hit name num   num arry arry',/,51x,
     4                           '  num  num',/,
     5 ' ____ ___ ______ _____ _____ ___ ___',5(' ____'))
        rat=ssqrw/ssqrw1
c  also write new station corrections to file22 to be used as input
c  in future runs.
      if(kout.ge.2) write(22,2201) nsts,day,tm
 2201 format(i4,24x,'computed',2x,a9,1x,a8)
      do 120 i=1,nsts
         ielev= stc(3,i)*(-1.0e+03)
         write(22,2207) stn(i),ltds(i),sltm(i),lnds(i),
     2      slnm(i),ielev,pdl(i),sdl(i),nfixst(i)
 2207    format(2x,a4,i2,1x,f5.2,i4,1x,f5.2,i5,2f5.2,i3)
         if(nfixst(i).eq.1) goto 118
         nvp=nparv+i
         if(hit(nvp).eq.zero) goto 115
         stderr(nvp)=rat*stderr(nvp)
         nvpi=ndexfx(nvp)
         nvpa=index(nvpi)
         ihit=nint(hit(nvp))
         write(16,125) stn(i),pdl(i),stderr(nvp),drm(nvp),
     2      nrd(i,1),ihit,stn(i),i,nvp,nvpi,nvpa
  125    format(1x,a4,' P  ',f7.3,2f6.3,2i4,1x,a4,2i5,2i5)
  115    continue
         if(iuses.eq.2) then
            nvs=nvp+nsts
            if(khit(nvs).eq.izero) goto 120
            stderr(nvs)=rat*stderr(nvs)
            nvsi=ndexfx(nvs)
            nvsa=index(nvsi)
            ihit=nint(hit(nvs))
            write(16,195) sdl(i),stderr(nvs),drm(nvs),
     2         nrd(i,2),ihit,stn(i),i,nvs,nvsi,nvsa
  195       format(5x,' S-P',f7.3,2f6.3,2i4,1x,a4,2i5,2i5)
         endif
  118    if(kout.lt.2) goto 120
  120 continue
      if(kout.ge.2) close(22)
c
 353  continue
c  skip velocity output if all fixed
      if(nparvi.eq.0) goto 900
c  compute vp/vs
c     if(iuses.eq.1) goto 365
c     do 360 k=1,nz
c        ks=k+nz
c        do 358 j=1,ny
c           do 356 i=1,nx
c              vpvs(i,j,k)=vel(i,j,k)/vel(i,j,ks)
c 356       continue
c 358    continue
  360 continue
c
  365 if(kout.lt.2) goto 311
c  Write the final velocity model to file23 in a format that
c  can be used as input for another run.                     
c  and to file25 in station format for plotting.
      open(unit=23,file='fort.23',status='new')
      if(kout.gt.3)
     2 open(unit=25,file='fort.25',status='new')
      write(23,2302) bld,nx,ny,nz,iuses,day,tm
 2302 format(f4.1,4i3,7x,'computed',2x,a9,1x,a8)
      write(23,2304) (xn(i),i=1,nx)
      write(23,2304) (yn(i),i=1,ny)
      write(23,2304) (zn(i),i=1,nz)
 2304 format(20f6.1)
      write(23,2305) izero,izero,izero
 2305 format(3i3)
      do 38 kv=1,iuses
         do 36 k=1,nz
            if((kout.gt.3).and.(k.gt.1).and.(k.lt.nz))
     2         write(25,2501) k,zn(k)
 2501       format(1x,'LAYR',1x,i2,8x,f7.2,'k')
            k2=k+(kv-1)*nz
            write(23,2311) (vel(i,1,k2),i=1,nx)
            do 34 j=2,ny-1
               write(23,2311) (vel(i,j,k2),i=1,nx)
               if((k.eq.1).or.(k.eq.nz)) goto 34
               if(kout.le.3) goto 34
               do 33 i=2,nx-1
                  call latlon(xn(i),yn(j),lat,xlat,lon,xlon)
                  write(25,2502) vel(i,j,k2),lat,xlat,lon,xlon,
     2            xn(i),yn(j),zn(k)
 2502             format(2x,f4.2,i2,f5.2,1x,i3,f5.2,10x,3f7.2)
   33          continue
   34       continue
            write(23,2311) (vel(i,ny,k2),i=1,nx)
   36    continue
   38 continue
c  write Vp/Vs at end of new velocity file
      if(iuses.eq.1) goto 50
      do 42 k=1,nz
         write(23,2311) (vpvs(i,1,k),i=1,nx)
         do 134 j=2,ny-1
            write(23,2311) (vpvs(i,j,k),i=1,nx)
            if((k.eq.1).or.(k.eq.nz)) goto 134
            if(kout.le.3) goto 134
            do 133 i=2,nx-1
               call latlon(xn(i),yn(j),lat,xlat,lon,xlon)
               write(25,2502) vpvs(i,j,k),lat,xlat,lon,xlon,
     2         xn(i),yn(j),zn(k)
  133       continue
  134    continue
         write(23,2311) (vpvs(i,ny,k),i=1,nx)
   42 continue
c
 2311 format(20f5.2)
   50 close(unit=23)
      if(kout.gt.3) close(unit=25)
c
  311 do 495 kv=1,iuses
         if((kv.eq.1).and.(iusep.eq.0)) goto 495
  312    write(16,1000) type(kv)
 1000    format(//,' FINAL ',a1,'-VELOCITY MODEL')
         do 10 k=1,nz
            k2=k+(kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 10
c write out vp/vs
            if(kv.eq.2) then
               write(16,1009) k,vtype(kv),zn(k)
               do 9 j=1,ny
                  write(16,1001) (vpvs(i,j,k),i=1,nx)
    9          continue
            endif
            write(16,1010) k,type(kv),zn(k)
 1009       format(/,' layer',i3,5x,a5,' nodes',10x,'z =',f7.1)
 1010       format(/,' layer',i3,5x,a1,'-velocity nodes',10x,
     2         'z =',f7.1)
            do 11 j=1,ny
               write(16,1001) (vel(i,j,k2),i=1,nx)
   11       continue
   10    continue
c
         write(16,1003) 
         nz1=nz-1
         ny1=ny-1
         do 22 k=2,nz1
            k2=k+(kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 22
            write(16,1009) k,vtype(kv),zn(k)
            kk=k+(kv-1)*nz2
            do 20 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1005) (khit(i),i=n1,n2),izero
   20       continue
   22    continue
c
 1003    format(/,' OBSERVATION MATRIX - KHIT - (will be 0 for',
     2      ' fixed nodes)')
 1005    format('     0',20i6)
 1001    format(20f6.2)
c
         write(16,1633)
 1633    format(/,' DERIVATIVE WEIGHT SUM ')
         do 487 k=2,nz1
            k2=k+(kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 487
            write(16,1009) k,vtype(kv),zn(k)
            kk=k+(kv-1)*nz2
            do 486 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1635)  (hit(i),i=n1,n2),zero
 1635          format(' 0.',18f7.0)
  486       continue
  487    continue
c
         if(ires.eq.0) goto 495
c
c  assign -1.0 to diagonal resolution for fixed gridpoint
         do 490 n=1,nparv
            if(nfix(n).eq.1) drm(n)=-1.0
  490    continue
c
         write(16,1634)
 1634    format(/,' RESOLUTION : GRIDOINT NUMBER, DIAGONAL RESOLUTION',
     2      ' ELEMENT (-1.0 indicates fixed velocity gridpoint)')
         do 489 k=2,nz1
            k2=k+(kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 489
            write(16,1009) k,vtype(kv),zn(k)
            kk=k+(kv-1)*nz2
            do 488 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1636) (i,drm(i),i=n1,n2)
 1636          format(1x,10(i5,':',f7.4))
  488       continue
  489    continue
  495 continue
c
      rat=ssqrw/ssqrw1
      if(ires.eq.0) goto 900
c
c  write out pointer from full nodes to inversion nodes
      if(inf.le.0) goto 545
      open(unit=18,file='fort.18',status='new')
      write(18,1800) npar,nparv,npari,nparvi
 1800 format(4i7,'   0.0001')
      do 540 l=1,npar
         write(18,1800) l,ndexfx(l)
  540 continue
      close(18)
c
  545 do 550 l=1,npar
         if(khit(l).eq.0) goto 550
         if ((hit(l).lt.hitct).or.(nfix(l).eq.1)) go to 550
         stderr(l)=rat*stderr(l)
c        write(16,2004) l,drm(l),stderr(l)
 2004    format('  parameter number ',i3,': resolution:',f8.4,
     *      '  error (%):',f8.4)
c
  550 continue
c  Also write resolution in more readable format
      write(16,1026)
 1026 format(/,10x,'VELOCITY WITH STANDARD ERROR(km/s) AND RESOLUTION'
     2 /,'  where error and res are both 0.00, the node was not',
     3 ' inverted for',
     4 ', where resol.=-1.00, gridpoint velocity was fixed')
      nx2=nx-2
      nxy2=nx2*(ny-2)
c  create proper sized format for nx nodes
      if(nx2.le.10) then
        write(line,27) nx2
        read(line,28) ch2nx2
   27   format(1x,i2)
   28   format(1x,a2)
        formvl(11:12)=ch2nx2
        formrs(14:15)=ch2nx2
      else
        nx3=nx2-10
        write(line,25) nx3
        read(line,26) ch2nx3
   25   format(1x,i1)
   26   format(1x,a1)
        nx4=10-nx3
        write(line,25) nx4
        read(line,26) ch2nx4
        frmvl2(9:9)=ch2nx4
        frmvl2(16:16)=ch2nx3
        frmrs2(5:5)=ch2nx4
        frmrs2(12:12)=ch2nx3
      endif
      do 75 kv=1,iuses
         do 70 k=2,nz-1
            k2=k+(kv-1)*nz 
            if(sumhit(k2).eq.0.0) goto 70
            kk=k+(kv-1)*nz2
            write(16,1009) k,vtype(kv),zn(k)
            if(kv.eq.1) then
               write(16,1021) (vel(i,1,k2),i=1,nx)
            else
               write(16,1021) (vpvs(i,1,k),i=1,nx)
            endif
            kp=(kk-2)*nxy2
            do 60 j=2,ny-1
               jp=(j-2)*nx2
               kjp=kp+jp
c  convert standard error from % to km/s
               do 55 i=1,nx2
                  v=vel(i+1,j,k2)
                  if(kv.eq.2) v=vpvs(i+1,j,k)
                  sekms(i)=stderr(kjp+i)*v
                  if(stderr(kjp+i).gt.tol) then
                     covdi(kjp+i)=1.0/(sekms(i)*sekms(i))
                  endif
   55          continue
               if(nx2.le.10) then
                  if(kv.eq.1) then
                     write(16,formvl) vel(1,j,k2),(vel(i+1,j,k2),
     2                  sekms(i),i=1,nx2),vel(nx,j,k2)
                  else
                     write(16,formvl) vpvs(1,j,k),(vpvs(i+1,j,k),
     2                  sekms(i),i=1,nx2),vpvs(nx,j,k)
                  endif
                  write(16,formrs) (drm(kjp+i),i=1,nx2)
               else
                  if(kv.eq.1) then
                     write(16,frmvl1) vel(1,j,k2),(vel(i+1,j,k2),
     2                  sekms(i),i=1,nx2),vel(nx,j,k2)
                  else
                     write(16,frmvl1) vpvs(1,j,k),(vpvs(i+1,j,k),
     2                  sekms(i),i=1,nx2),vpvs(nx,j,k)
                  endif
                  write(16,frmrs1) (drm(kjp+i),i=1,10)
                  if(kv.eq.1) then
                     write(16,frmvl2) (vel(i+1,j,k2),sekms(i),
     2                  i=11,nx2),vel(nx,j,k2)
                  else
                     write(16,frmvl2) (vpvs(i+1,j,k),sekms(i),
     2                  i=11,nx2),vpvs(nx,j,k)
                  endif
                  write(16,frmrs2) (drm(kjp+i),i=11,nx2)
               endif
   60       continue
            if(kv.eq.1) then
               write(16,1024) (vel(i,ny,k2),i=1,nx)
            else
               write(16,1024) (vpvs(i,ny,k),i=1,nx)
            endif
   70    continue
   75 continue
   16 format('0',f4.2,12(f7.2,'+',f4.2))
 1021 format(f5.2,f7.2,10f12.2)
 1023 format(1x,'res',6x,10('(',f4.2,')',6x))
 1024 format('0',f4.2,f7.2,10f12.2)
c
c  write out 1/covariance(diag) file in format like velocity file
c  so can be plotted
      if((kout2.eq.5).and.(ires.gt.0)) then
        open(unit=45,file='fort.45',status='new')
        write(45,4502) nx,ny,nz,day,tm
 4502   format(3i3,6x,'1/Diag. Covariance, Computed',2x,a9,
     2   1x,a8)
        write(45,2304) (xn(i),i=1,nx)
        write(45,2304) (yn(i),i=1,ny)
        write(45,2304) (zn(i),i=1,nz)
        write(45,2305) izero,izero,izero
      write(16,1739)
 1739 format(/,' 1/(Covariance Diag. )')
      do 746 kv=1,iuses
c  write peripheral grids to file45 also
         do 725 j=1,ny
            write(45,1737) (zeroi(i),i=1,nx)
  725    continue
         do 736 k=2,nz1
            k2=k+(kv-1)*nz
            if(sumhit(k2).gt.0.0)
     2         write(16,1009) k,vtype(kv),zn(k)
            kk=k+(kv-1)*nz2
            write(45,1737) (zeroi(i),i=1,nx)
            do 735 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               if(sumhit(k2).gt.0.0)
     2            write(16,1737) zero,(covdi(i),i=n1,n2),zero
 1737          format(10e10.3)
               write(45,1737) zero,(covdi(i),i=n1,n2),zero
  735       continue
  736    continue
         write(45,1737) (zeroi(i),i=1,nx)
         do 745 j=1,ny
            write(45,1737) (zeroi(i),i=1,nx)
  745    continue
  746 continue
      close(45)
      endif
c
  900 continue
      if(i3d.eq.0) goto 950
      upbavg=upbtot/(float(nupb))
      write(16,1665) nupb, upbavg
 1665 format(/,' For',i10,' calls to MINIMA, average number',
     2 ' psuedo-bending iter. =',f7.2)
  950 call date_and_time(TIME=tm)
      write(16,1660) day,tm
 1660 format(/,' Computation finished at ',a9,1x,a8)
c***** end of subroutine outend *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outlis(ne,rmk)
c  write out each event in a format similar to HYPO71 listing format
c  so that it can be used as input to FPFIT fault-plane solution
c  program
c
c  declaration statements:
      parameter (one=1.00)
      character*1 phs(2)
      character*9 day
      character*8 tm
c
c  common block variables:
      include 'simulps_common.inc'
      common/wtpars/w(maxobs)

      character*1 rmk(maxobs, maxev) 
c
      data phs/'P','S'/
      if (ne.gt.1) goto 5
c  print header line
      call date_and_time(DATE=day)
      call date_and_time(TIME=tm)
      write(34,3407) day,tm
 3407 format(/,' AZIM and TOA calculated with 3D velocity model ',
     2 '(SIMULPS9) ',a9,1x,a8,/)
    5 write(34,3401)
 3401 format(2x,'DATE',4x,'ORIGIN',3X,'LATITUDE LONGITUDE  DEPTH',4X,
     2 'MAG NO',11X,'RMS')
      call latlon(evc(1,ne),evc(2,ne),lat,xlat,lon,xlon)
  303 write(34,3402) iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 lat,xlat,lon,xlon,evc(3,ne),rmag(ne),kobs(ne),wrmsr(ne)
 3402 format(1x,a4,a2,1x,a2,i2,f6.2,1x,i2,'n',f5.2,1x,i3,'w',f5.2,f7.2,
     2 2x,f5.2,i3,9x,f5.2)
      write(34,3403)
 3403 format(/,2x,'STN  DIST  AZ TOA PRMK HRMN  PSEC TPOBS',
     2 14X,'PRES  PWT')
      nobs=kobs(ne)
      do 100 i=1,nobs
c  don't use S arrivals for fault-plane solution
         if(intsp(i,ne).eq.1)  goto 100
         iaz=nint(az(i))
         itoa=nint(toa(i))
         write(34,3405) stn(isto(i,ne)),dlta(i,ne),iaz,
     2      itoa,rmk(i,ne),ihr(ne),mino(ne),secp(i,ne),
     3      (secp(i,ne)-seco(ne)),res(i),w(i)
  100 continue
 3405 format(1x,a4,1x,f5.1,1x,i3,1x,i3,1x,a4,
     2 1x,a2,i2,1x,f5.2,1x,f5.2,
     3 13x,f5.2,1x,f5.2)
      write(34,3406)
 3406 format(/)
      return
c ****** end of subroutine outlis ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outres(ne,rmk)
c  write out station residuals and weights ( from parsep) for each event
c
c  common block variables:
      include 'simulps_common.inc'
      common/wtpars/w(maxobs)

      character*1 rmk(maxobs, maxev)
c
c  declaration statements:
      real ttobs(maxobs)
      character*1 phs(2)
      parameter (one=1.00)
c
      data phs/'P','S'/
c
      do 10 j=1,kobs(ne)
         ttobs(j)=secp(j,ne)
         if(intsp(j,ne).eq.0) ttobs(j)=secp(j,ne)-seco(ne)
   10 continue
      if(ne.gt.(neqs+nbls)) goto 200
      write(16,3) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
      write(20,3) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
 3    format(1h0,' residuals and weights(parsep) for event=',i4, 
     2 3x,a4,a2,1x,a2,i2,1x,f5.2,'  mag ',f4.2)
      write(16,51)
51    format(1x,4('sta  ph  wt  res:O-C ttobs delta',1x))
      write(20,55)
   55 format(1x,'sta  ph  wt  res:O-C ttobs ttcal  delta',
     2 4x,'x ev',4x,'y ev',4x,'z ev',4x,'x st',4x,'y st',
     3 4x,'z st')
      write(16,53) (stn(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2 ttobs(j),dlta(j,ne),j=1,kobs(ne))
53    format(4(1x,a4,a4,f5.2,f7.3,2f6.2))
      do 20 j=1,kobs(ne)
         ttcal=ttobs(j)-res(j)
         write(20,54) stn(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2      ttobs(j),ttcal,dlta(j,ne),
     3      (evc(ie,ne),ie=1,3),(stc(is,isto(j,ne)),is=1,3)
   20 continue
   54 format(1x,a4,a4,f5.2,f7.3,2f6.2,f7.2,6f8.2)
      write(16,1601)
 1601 format(/)
      return
  200 continue
      write(16,33) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
      write(20,33) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
 33   format(1h0,' residuals and weights(medder) for event=',i4, 
     2 3x,a4,a2,1x,a2,i2,1x,f5.2,'  mag ',f4.2)
      write(16,51)
      write(20,55)
      write(16,53) (stn(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2 ttobs(j),dlta(j,ne),j=1,kobs(ne))
      do 30 j=1,kobs(ne)
         ttcal=ttobs(j)-res(j)
         write(20,54) stn(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2      ttobs(j),ttcal,dlta(j,ne),
     3      (evc(ie,ne),ie=1,3),(stc(is,isto(j,ne)),is=1,3)
   30 continue
      write(16,1601)
      return
c ****** end of subroutine outres ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine parsep(ne,nwr)
c  routine to perform parameter separation using qr
c
c  common block variables:
      include 'simulps_common.inc'
      common/wtpars/w(maxobs)

      character*1 rmk(maxobs, maxev)
c
c  declaration statements:
      real c1(maxobs,maxobs),u(maxobs),up
      parameter(zero=0.0,one=1.0)
c
      nobs=kobs(ne)
c
c  for shots, don't compute dtmp
      if (ne.gt.(neqs+nbls)) go to 100
c
c  compute weighting
c
      nwr=0
      wnorm=0.0
      do 13 j=1,nobs
c  reading weight
         w(j)=wt(j,ne)
c  residual weighting
c  downweighting(linear) 0 to 98% res1 to res2, 98 to 100% res2 to res3
         ares=abs(res(j))
         if(ares.le.res2) then
            wr=1.0-(ares-res1)*dres12
            if (wr.gt.1.0) wr=1.0
         else
            if(res3.gt.res2) then
               wr=0.02-(ares-res2)*dres23
               if (wr.lt.0.0) wr=0.0
            else
               wr=0.0
            endif
         endif
c  distance weighting
         wd=1.0-(dlta(j,ne)-delt1)*ddlt
         if (wd.gt.1.0) wd=1.0
        if (wd.lt.0.0) wd=0.0
c  unnormalized weight
         w(j)=w(j)*wr*wd
c  19-jul-1983 Change normalizing factor
         wnorm=wnorm+w(j)
c        wnorm=wnorm+w(j)*w(j)
         if (w(j).gt.0.0) nwr=nwr+1
   13 continue
c  check to be sure 4 or more readings with nonzero weight
c    dmep 02dec83
      if(nwr.ge.4) goto 12
      write(16,1612)
 1612 format(' ****** less than 4 readings with nonzero weights',
     2 /,20x,'SKIP TO NEXT EVENT *********')
      call outres(ne,rmk)
      return
c  normalize weights
c    19-jul-1983 Change in normalizing factor
c     wfac=sqrt(nwr/wnorm)
   12 wfac=nwr/wnorm
      if(ne.gt.netemp) wfac=wfac*wtsht
      resmax=zero
      do 14 j=1,nobs
         w(j)=w(j)*wfac
         ares=abs(res(j))
         if(ares.gt.resmax) resmax=ares
   14 continue
      if(resmax.gt.(3.0*res2)) call outres(ne,rmk)
      nwrt=nwrt+nwr
      if(nitmax.le.0) return
c
      do 10 i=1,nobs
         do 15 j=1,nobs
            c1(j,i)=zero
   15    continue
         c1(i,i)=w(i)
   10 continue
c  apply weighting to hypocenter matrix (1=o.t., 2,3,4=hypocenter)
      nsep=4
c  for blast, only do origin time
      if(ne.gt.neqs) nsep=1
      do 17 i=1,nsep
         do 16 j=1,nobs
            dth(j,i)=dth(j,i)*w(j)
   16    continue
   17 continue
c
c  perform qr decomposition
c  prepare parameters for input to h12 routine
c  loop over the four columns of dth for decomposition
      do 20 i=1,nsep
         irow=i
         irow1=i+1
c  put column of dth into u-vector
         do 22 j=1,nobs
            u(j)=dth(j,i)
   22    continue
c  compute single transformation
         call h12(1,irow,irow1,nobs,u,1,up,dth,1,maxobs,4,
     2      c1,1,maxobs,nobs)
   20 continue
c
c
c  qr decomposition complete - multiply residuals and medium matrix
c  by u0 matrix (q matrix minus its first four columns)
      nobs4=nobs-nsep
      do 35 j=1,nobs
         resj=res(j)
         do 34 i=1,nobs4
c  28 jan84 Changed from 4 to nsep to match CT. dmep
            i4=i+nsep
            resp(i)=resp(i)+c1(i4,j)*resj
   34    continue
   35 continue
c  compute separated medium matrix and store
      do 30 n=1,nobs
         nn=(n-1)*npari
         do 37 i=1,nobs4
            ii=(i-1)*npari
            i4=i+nsep
            c4n=c1(i4,n)
            do 36 m=1,npari
               im=ii+m
               nm=nn+m
               dtmp(im)=dtmp(im)+c4n*dtm(nm)
   36       continue
   37    continue
   30 continue
c
c  add to elements of g matrix
c
  100 continue
c
c  compute contribution to ssqr
      sqwtsh=sqrt(wtsht)
      totrms(ne)=0.0
      do 18 i=1,nobs
         resi=res(i)
         resi2=resi*resi
         totrms(ne)=totrms(ne)+resi2
         rrw=resi2*w(i)
         ssqrw=ssqrw+rrw
         wnobt=wnobt+w(i)
         if(intsp(i,ne).eq.0) then
c  P arrival
            ssqrwp=ssqrwp+rrw
            wnobtp=wnobtp+w(i)
         else
c  S arrival
            ssqrws=ssqrws+rrw
            wnobts=wnobts+w(i)
         endif
         if (ne.gt.netemp) resi=resi*sqwtsh
         ssqr=ssqr+resi*resi
   18 continue
      totrms(ne)=sqrt(totrms(ne)/nobs)
c
c  loop over all observations (in separated form)
      n1=1
      n2=nobs4
      if (ne.gt.(neqs+nbls)) n2=nobs
c     write(16,359)n1,n2,ne,neqs
c359  format(' parsep: n1,n2,ne,neqs= ',4i8)
c  loop over observations
      do 140 i=n1,n2
c        write(16,364)i,n1,n2
c364     format(' parsep: i,n1,n2= ',3i8)
         ii=(i-1)*npari
c  for a given observation, loop over all nodes
         do 130 j=1,npari
            jj=j+ii
            if (dtmp(jj).eq.zero) go to 130
c  add unknown for first observation only
            if (khit(mdexfx(j)).gt.0) go to 125
            mbl=mbl+1
            index(j)=mbl
c  zero next row of g matrix
            mbtot=nbtot+1
            nbtot=nbtot+mbl
            do 120 l=mbtot,nbtot
               g(l)=0.0
  120       continue
  125       khit(mdexfx(j))=khit(mdexfx(j))+1
c CHECK (only for shots since parameter separation messes up indices)
c           if (ne.le.(neqs+nbls)) go to 126
c           if((intsp(i,ne).eq.1).and.(j.le.nodes2))
c    2         write(36,3601) dtmp(jj),jj,j,i,khit(j),index(j)
c3601       format(' Parsep: dtmp(jj)=',f7.2,'jj=',i6,',j=',i4,
c    2         ',i=',i4,',khit(j)=',i5,',index(j)=',i5)
  126       k=index(j)
            kk2=((k-1)*k)/2
c  build rhs
            rhs(k)=rhs(k)+resp(i)*dtmp(jj)
c  build g matrix
            do 128 l1=1,j
               lj=l1+ii
               if (dtmp(lj).eq.zero) go to 128
               l=index(l1)
               m=l+kk2
               if (k.lt.l) m=k+((l-1)*l)/2
               g(m)=g(m)+dtmp(jj)*dtmp(lj)
c              if(g(m).lt.0.0) write(16,1620) g(m),m,dtmp(jj),
c    2            jj,dtmp(lj),lj,k,l,l1,i,ii,j
c1620          format(' g(m)=',f6.3,'m=',i6,',dtmp(jj)=',f6.3,
c    2            'jj=',i6,',dtmp(lj)=',f6.3,'lj=',i6,',k=',i4,
c    3            ',l=',i4,',l1=',i4,',i=',i4,',ii=',i6,',j=',i4)
  128       continue
  130    continue
  140 continue
c  zero out dtm and dtmp matrices for next event
      nono=nobs*npari
      seprms(ne)=zero
      do 45 m=1,nono
         dtm(m)=zero
         dtmp(m)=zero
   45 continue
      do 945 i=1,nobs
         seprms(ne)=seprms(ne)+resp(i)*resp(i)
         resp(i)=zero
  945 continue
      rmssep=sqrt(seprms(ne)/nobs)
      if(ne.gt.(neqs+nbls)) then
         write(16,1698) totrms(ne),rmssep
 1698    format('  ** total rms =',f8.5,
     2      '  **  model rms =',f8.5,' **')
      else
         write(16,1699) nwr,totrms(ne),rmssep
 1699    format(5x,'nwr=',i3,', ','** total rms =',f8.5,
     2      '  **  model rms =',f8.5,' **')
      endif
      if(totrms(ne).ge.(2.0*res1)) call outres(ne,rmk)
      return
c***** end of subroutine parsep *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine path(ne,no,xe,ye,ze,ttime)
c  this routine determines the minimum-time ray path
c  in two steps:  first, an approximate path is
c  determined using approximate ray tracing
c  then the approximate path is used as a starting point in
c  shooting ray tracing routine to determine the path in 3-d
c  ***  note - current version does not contain full 3-d ray
c  ***  tracing - routines are under development
c
c  declaration statements:
c
c  common block variables:
      include 'simulps_common.inc'
c
c  set up parameters to call ray tracing routine
c
c  S or P reading
      isp=intsp(no,ne)
c  receiver coordinates
      ns=isto(no,ne)
      xr=stc(1,ns)
      yr=stc(2,ns)
      zr=stc(3,ns)
c  determine 3-d path using circular raypaths.
c  determine approximate path
  120 jflag=jfl
      ncrold=nco(no,ne)
      ndpold=ndo(no,ne)
c
      call rayweb(no,isp,xe,ye,ze,xr,yr,zr,
     * fstime,jflag,ncrold,ndpold)
      ttime=fstime
      ttc(no)=ttime
      nco(no,ne)=ncrold
      ndo(no,ne)=ndpold
c
c  do psuedo-bending if i3d>0
      if (i3d.eq.0) return
c
  125 continue
c  number of pb iter depends on distance
      nitpbu=nitpb(1)
      if(dlta(no,ne).gt.delt1) nitpbu=nitpb(2)
      call minima(no,isp,ttime,nitpbu,jpb)
      if(jpb.lt.nitpbu) goto 130
        write(26,2601) no,stn(ns),jpb,nitpbu
 2601   format(' Minima: no=',i4,2x,a4,', used maximum ',
     2  'number PB iter.: j=',i3,', nitpb=',i3)
  130 continue
      ttc(no)=ttime
      if(kout3.eq.0) return
c  write out travel-time differences to file 19
      tdif=fstime-ttime
      write(19,1900) ne,stn(isto(no,ne)),dlta(no,ne),fstime,ttime,tdif
 1900 format(i4,1x,a4,f7.2,3f8.4)
c
      return
c***** end of subroutine path *****
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine rayweb(no,isp,xe,ye,ze,xr,yr,zr,
     *            fstime,jflag,ncrold,ndpold)
c  approximate ray tracing package art2
c    with fast ray tracing code
c     by Cliff Thurber (from his simul3l version)
c
c  common block variables:
      include 'simulps_common.inc'
c
c  declaration statements:
c  parameters
      real xe,ye,ze,xr,yr,zr,fstime
c  local variables
      real delx,dely,delz,sep,delsep,pthsep(130),strpth(390),
     * fstpth(390),dipvec(3,9),disvec(390,9),trpath(390,9),
     * trtime(9),tmin,tt,trpth1(390)
      integer nd,ns,npt,ncr,i,ic,ndp,nc,n1,n2,n3,nn,np,ndpfst
c
c  compute source-receiver separation
      delx=xr-xe
      dely=yr-ye
      delz=zr-ze
      sep=sqrt(delx*delx+dely*dely+delz*delz)
c  determine integer parameters for set of curves to be constructed
      call setup(sep,scale1,scale2,nd,ns,npt,ncr)
c
c  set up pthsep array for straight-line travel time calculation
      sn1=1.0/ns
      delsep=sep*sn1
      do 20 i=1,ns
         pthsep(i)=delsep
   20 continue
c
c  determine points along straight-line path
      xstep=delx*sn1
      ystep=dely*sn1
      zstep=delz*sn1
c
      ic=0
      ns1=ns+1
      do 25 ii=1,ns1
c
         i=ii-1
c
         ic=ic+1
         strpth(ic)=xe+xstep*i
         fstpth(ic)=strpth(ic)
         ic=ic+1
         strpth(ic)=ye+ystep*i
         fstpth(ic)=strpth(ic)
         ic=ic+1
         strpth(ic)=ze+zstep*i
         fstpth(ic)=strpth(ic)
   25 continue
c
c  compute travel time along straight-line path
      call ttime(isp,ns,npt,strpth,pthsep,fstime)
c
      if (ncr.eq.1) go to 65
c
c  compute the dip vectors of length scale2
      call cmpdpv(xe,ye,ze,xr,yr,zr,scale2,ndip,dipvec)
c
c  compute the basic set of displacement vectors
      call cmpdsv(ndip,iskip,ns,dipvec,disvec)
c
c  set first and last points of all trial paths to source and receiver
      n1=3*npt-2
      n2=n1+1
      n3=n1+2
      ndip1=1+iskip
      ndip2=ndip-iskip
c
c  fast ray tracing code
c
      nz1=nz-1
      ncr0=1
      ncr1=ncr-1
      if (jflag.eq.0) go to 28
      ncr0=ncrold-1
      if (ncr0.lt.1) ncr0=1
      ncr1=ncrold+1
      if (ncr1.gt.ncr) ncr1=ncr
c  ndip was changed (in main) for ihomo iterations., make ndpold=vertical plane.
c     if(jfl.eq.1) ndpold=(ndip+1)/2     ! 21-feb-86, this is done in strt, so unnecessary here
      ndip1=ndpold-1
      if (ndip1.lt.1+iskip) ndip1=1+iskip
      ndip2=ndpold+1
      if (ndip2.gt.ndip-iskip) ndip2=ndip-iskip
   28 continue
c  set "old" values to straight line
      ncrold=0
      ndpold=(ndip+1)/2
c
c
      do 30 ndp=ndip1,ndip2
         trpath(1,ndp)=xe
         trpath(2,ndp)=ye
         trpath(3,ndp)=ze
         trpath(n1,ndp)=xr
         trpath(n2,ndp)=yr
         trpath(n3,ndp)=zr
   30 continue
      trpth1(1)=xe
      trpth1(2)=ye
      trpth1(3)=ze
      trpth1(n1)=xr
      trpth1(n2)=yr
      trpth1(n3)=zr
c
c  loop over the curve sets
      do 40 nc=ncr0,ncr1
         iz0=0
c
c  loop over different dips for one set of curves
         do 42 ndp=ndip1,ndip2
c
            npt2=npt-2
c  loop to determine points along one path
            do 44 np=1,npt2
               n1=3*np+1
               n3=n1+2
               do 43 nn=n1,n3
                  trpath(nn,ndp)=nc*disvec(nn,ndp)+strpth(nn)
                  trpth1(nn)=trpath(nn,ndp)
   43          continue
               if((i3d.lt.3).or.(ze.gt.zn(nz1))) goto 44
c  Use less curvature if below "moho"
               if(trpth1(n3).gt.zn(nz1)) then
                  if(iz0.eq.0) then
                     z0=trpth1(n3)-disvec(n3,ndp)
                     iz0=1
                  endif
                  trpath(n3,ndp)=disvec(n3,ndp)+z0
c                 write(6,600) trpth1(n3),zn(nz1),trpath(n3,ndp)
  600             format('trpth1(n3)',f7.2,',zn(nz1)',f7.2,
     2               'new trpth1(n3)',f7.2)
                  trpth1(n3)=trpath(n3,ndp)
               endif
   44       continue
c
c  set up pthsep array for travel time calculations
            if (ndp.eq.ndip1) call cmpsep(trpth1,pthsep,ns)
c  compute travel time along one path
            call ttime(isp,ns,npt,trpth1,pthsep,tt)
            trtime(ndp)=tt
   42    continue
c
c  sort through trtime to find fastest path from current set
         tmin=1.0e15
         do 50 ndp=ndip1,ndip2
            if (trtime(ndp).gt.tmin) go to 50
            tmin=trtime(ndp)
            ndpfst=ndp
   50    continue
c
c  compare fastest trtime to current value of fstime
c  replace fstime and fstpth if needed
         if (tmin.ge.fstime) go to 40
         fstime=tmin
c  reset "old" values
         ncrold=nc
         ndpold=ndpfst
c
         npt3=3*npt
         do 52 np=1,npt3
            fstpth(np)=trpath(np,ndpfst)
   52    continue
c
   40 continue
c
c  put fstpth into rp array
   65 continue
      do 60 np=1,npt
         n3=3*np
         n1=n3-2
         n2=n1+1
         rp(1,np,no)=fstpth(n1)
         rp(2,np,no)=fstpth(n2)
         rp(3,np,no)=fstpth(n3)
   60 continue
      nrp(no)=npt
c
c***** end of subroutine rayweb *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine rescov
c  compute resolution and covariance
c  common block variables:
      include 'simulps_common.inc'
c
      varnce=ssqrw/float(nobt-mbl-4*neqs)
c     write(16,1001) varnce
 1001 format(/,' RESCOV: data variance is ',f10.5)
c  resolution and error calculations
      if(ires.eq.2) rewind (17)
c  avoid printing zeroth element of rhs array
c  index point from nodes(input) to solution arrays(rhs)
      do 5 kc=1,npari
         if(khit(mdexfx(kc)).eq.0) index(kc)=npari+1
         if(hit(mdexfx(kc)).lt.hitct) index(kc)=npari+1
         if(index(kc).eq.0) index(kc)=npari+1
    5 continue
c  zero out drm and stderr elements
      do 300 l=1,npar
         drm(l)=0.0
         stderr(l)=0.0
  300 continue
c
   10 do 500 l=1,npari
         l1=l
         l1n=mdexfx(l1)
         if(khit(l1n).eq.0) goto 500
         if(hit(l1n).lt.hitct) goto 500
         k=index(l1)
         jj=0
         ii=0
c put appropriate elements of g1 into rhs
         do 435 i=1,mbl
            do 430 j=1,i
               jj=jj+1
               if (i.eq.k.or.j.eq.k) go to 420
               go to 430
  420          ii=ii+1
               rhs(ii)=g1(jj)
c CHECK
c              if(g1(jj).lt.-5.0) write(16,1620)
c    2            g1(jj),jj,ii,i,j,mbl,l,k,npar
c1620          format(' g1(jj)=',f9.3,'jj=',i6,',ii=',i4,',i=',i4,
c    2            ',j=',i4,',mbl=',i4,',l=',i4,',k=',i4,',npar=',i4)
  430       continue
  435    continue
c  compute vector of resolution matrix for this parameter
         call luelmp(g,rhs,mbl,rhs)
c print out full resolution matrix if ires=2 or 3
         if(ires.lt.2) goto 440
         write(17,1720) l1,mdexfx(l1)
 1720    format(/,'row=',i5,', gridpoint no.=',i5,/)
         write(17,1721) (rhs(index(kc)),kc=1,npari)
 1721    format(18f7.4)
c  store diagonal element
  440    drm(l1n)=rhs(k)
c  standard error and covariance calculation
         do 450 i=1,mbl
            rhs(i)=rhs(i)*varnce
  450    continue
c  compute vector of covariance matrix for this parameter
         call luelmp(g,rhs,mbl,rhs)
c  store standard error of slowness perturbation
         rhs(k)=abs(rhs(k))
         stderr(l1n)=sqrt(rhs(k))
  500 continue
  900 return
c***** end of subroutine rescov *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine setorg(ltdo,oltm,lndo,olnm)
c  this routine establishes the short distance conversion factors
c  given the origin of coordinates
c  the rotation angle is converted to radians also
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,xlt,xln,snr,csr
c  local variables:
      double precision dlt1,dlt2,dxlt,drad,drlt
      parameter ( re=6378.163, ell=298.26)
      parameter (drad=1.7453292d-2)
      parameter (drlt=9.9330647d-1)
      rota=rota*drad       !convert from degrees to radians
c
      dxlt=dble(60.*ltdo+oltm)
      xln=60.*lndo+olnm
      xlt=sngl(dxlt)
c  conversion factor for latitude
      dlt1=datan(drlt*dtan(dxlt*drad/60.d0))
      dlt2=datan(drlt*dtan((dxlt+1.d0)*drad/60.d0))
      del=sngl(dlt2-dlt1)
      r=re*(1.0-sngl(dsin(dlt1)**2)/ell)
      xltkm=r*del
c  conversion factor for longitude
      del=sngl(dacos(1.0d0-(1.0d0-dcos(drad/60.d0))*dcos(dlt1)**2))
      bc=r*del
      xlnkm=bc/sngl(dcos(dlt1))
      write(16,3001) xltkm,bc
 3001 format(5x,'short distance conversion factors',/,
     2       10x,'one min lat  ',f7.4,' km',/,
     3       10x,'one min lon  ',f7.4,' km',/)
c
c  convert coordinates with rotation cosines
      snr=sin(rota)
      csr=cos(rota)
      return
c***** end of subroutine setorg *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine setup(sep,scale1,scale2,nd,ns,npt,ncr)
c
c  parameters
      real sep,scale1,scale2
c
      integer nd,ns,npt,ncr
c
c  determine the number of path divisions - use scale1
      nd=1+nint(3.32193*log10(sep/scale1))
      if (nd.gt.7) nd=7
c  number of segments along path
      ns=2**nd
c  number of points on path
      npt=ns+1
c
c  determine the number of curves - use scale2
      ncr=1+0.5*sep/scale2
c
      if (sep.gt.scale1) return
      nd=0
      ns=1
      npt=2
      ncr=1
c
c***** end of subroutine setup *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine strt(nit)
c
c  declaration statements:
      parameter(zero=0.0)
c
c  common block variables:
      common/machin/ eta,tol
      include 'simulps_common.inc'
c
      iarr=maxpar
      iarri=mxpari
c
      nwrt=zero
      nswrt=zero
      ssqr=zero
      ssqrw=zero
      ssqrwp=zero
      ssqrws=zero
      wnobt=zero
      wnobtp=zero
      wnobts=zero
      do 77 n=1,iarri
         index(n)=0
         jndex(n)=0
         rhs(n)=zero
   77 continue
c
c
c  reset fast art arrays?
      if (nit.gt.ihomo) go to 80
      do 790 m=1,maxev
         do 79 n=1,maxobs
            ndo(n,m)=(ndip+1)/2
   79    continue
  790 continue
c
      if (nit.ge.1) go to 80
      do 780 m=1,maxev
         do 78 n=1,maxobs
            nco(n,m)=1
   78    continue
  780 continue
c
c
   80 mbl=0
      nbtot=0
c  zero hit counter arrays
      do 25 j=1,iarr
         hit(j)=zero
         khit(j)=0
   25 continue
c
c  zero out partial derivative arrays
      mxobs1=maxobs-1
      do 30 j=1,iarri
         jj2=maxobs*j
         jj1=jj2-mxobs1
         do 35 i=jj1,jj2
            dtm(i)=zero
            dtmp(i)=zero
   35    continue
   30 continue
      do 39 i=1,maxobs
         resp(i)=zero
   39 continue
c
      if (nit.gt.0) return
c
c  machine constants
      halfu=0.5
   50 temp1=1.0+halfu
      if (temp1.le.1.0) go to 100
      halfu=0.5*halfu
      go to 50
  100 eta=2.0*halfu
c
      temp2=0.5
  150 continue
      if (temp2.le.0.0) go to 200
      tol=temp2
      temp2=0.5*temp2
      go to 150
  200 continue
c
c     tol=tol/eta
c for unix make tol bigger
      tol=tol/(eta*eta)
c
      write(16,1055) eta,tol
 1055 format(/,'  computed machine constants eta and tol:',2e16.6)
c
      return
c***** end of subroutine strt *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine travel
c
      common/temp/xtemp(130),ytemp(130),ztemp(130),rtemp(130),ttemp(130)
c
      common/pathm/x(130),y(130),z(130),v(130),tra,n,nn
c
      tra=0
      do 60 i=2,n
         i1=i-1
         xd=xtemp(i)-xtemp(i1)
         yd=ytemp(i)-ytemp(i1)
         zd=ztemp(i)-ztemp(i1)
         ds=sqrt(xd*xd+yd*yd+zd*zd)
         tv=ds*(1.0/v(i)+1.0/v(i1))
         tra=tra+tv
  60  continue
      tra=0.5*tra
c
      return
c ***** end of subroutine travel *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ttime(isp,ns,npt,pathr,pthsep,tt)
c
c  travel time along path via trapezoidal rule integration
c
c  parameters
      real pathr(390),pthsep(130),tt
c
      integer ns,npt
c  local variables
      real vpt(130),x,y,z,v
c
      integer ip,np
c
c  loop over points along path to determine velocities
      ip=0
      do 10 np=1,npt
         ip=ip+1
         x=pathr(ip)
         ip=ip+1
         y=pathr(ip)
         ip=ip+1
         z=pathr(ip)
         call vel3(isp,x,y,z,v)
   10 vpt(np)=v
c
      tt=0.0
c  sum up travel time - use trapezoidal rule
      do 20 np=1,ns
         np1=np+1
c  Check for value outside defined area
         if((vpt(np).le.0.0).or.(vpt(np1).le.0.0)) goto 99
         tt=tt+pthsep(np)/(vpt(np)+vpt(np1))
   20 continue
      tt=2.0*tt
c
      return
   99 write(16,100) np,vpt(np),np1,vpt(np1)
  100 format(' **** ERROR IN TTIME ****, velocity le 0.0',
     2 /,'    np   vpt(np)   np1   vpt(np1)',/,1x,2(i5,f10.2))
      stop
c***** end of subroutine ttime *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ttmder(ne,no,kopt,ttime)
c
c  declaration statements:
      dimension in(8)
c
c  common block variables:
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'simulps_common.inc'
c
c  resegment ray path
c  calculate travel time derivatives with respect to hypocentral
c  parameters - use coords of first two points on raypath
c  determine slowness at source
      xe=rp(1,1,no)
      ye=rp(2,1,no)
      ze=rp(3,1,no)
c  check for P or S reading
      isp=intsp(no,ne)
      call vel3(isp,xe,ye,ze,v)
      us=1.0/v
c**  The program has been modified to include S-P times for the inversion.
c**  The corresponding travel-time derivatives are calculated.(06/29/92)
      if (isp.eq.1) then
              isp=0
              call vel3(isp,xe,ye,ze,vp)
              up=1.0/vp
              us=us-up
        isp=1
      endif
c**
c  determine cartesian derivatives from direction cosines
      dx=rp(1,2,no)-rp(1,1,no)
      dy=rp(2,2,no)-rp(2,1,no)
      dz=rp(3,2,no)-rp(3,1,no)
      ds=sqrt(dx*dx+dy*dy+dz*dz)
      uds=-us/ds
c  hypocentral derivatives
      dth(no,2)=uds*dx
      dth(no,3)=uds*dy
      dth(no,4)=uds*dz
c  origin time derivative
      dth(no,1)=1.0
c** S-P CHANGE
c   zero the origin time derivative for S-P
      if (isp.eq.1) dth(no,1)=0.0
c**
c  zero cartesian derivatives for quarry blasts
      if((ne.le.neqs).or.(ne.gt.(neqs+nbls))) goto 20
      dth(no,2)=0.0
      dth(no,3)=0.0
      dth(no,4)=0.0
   20 continue
c
c  skip next section if all velocity nodes are fixed (nparvi=0)
      if(nparvi.eq.0) goto 88
c  skip next section if only doing earthquake location
        if(kopt.eq.0) go to 88
c
c  travel time and velocity partial derivatives
      tt=0.0
      half=0.5
c  loop over segments comprising the ray path
      nrp1=nrp(no)-1
      pl(no)=0.0
      do 50 i=1,nrp1
         i1=i+1
         rx=rp(1,i,no)
         ry=rp(2,i,no)
         rz=rp(3,i,no)
         dx=rp(1,i1,no)-rx
         dy=rp(2,i1,no)-ry
         dz=rp(3,i1,no)-rz
c  compute segment length
         sl=sqrt(dx*dx+dy*dy+dz*dz)
         pl(no)=pl(no)+sl
c  decide on number of subsegments and compute length
         nseg=nint(sl/stepl)+1
         fnsegi=1.0/float(nseg)
         ssl=sl*fnsegi
         dxs=dx*fnsegi
         dys=dy*fnsegi
         dzs=dz*fnsegi
c
         xp=rx-half*dxs
         yp=ry-half*dys
         zp=rz-half*dzs
c  loop over subsegments
         do 55 is=1,nseg
            xp=xp+dxs
            yp=yp+dys
            zp=zp+dzs
c
            call vel3(isp,xp,yp,zp,v)
            dt=ssl/v
c
c
c  The next section is a change from 'block' to 'linear'
c   partial derivatives, by C.Thurber,may10,1983.
c  Nodes with non-zero weight
            in(1)=ip-1+nx2*(jp-2)+nxy2*(kp-2)-nxy2*(2*isp)
            in(2)=in(1)+1
            in(3)=in(1)+nx2
            in(4)=in(3)+1
            in(5)=in(1)+nxy2
            in(6)=in(5)+1
            in(7)=in(5)+nx2
            in(8)=in(7)+1
c
c  Assign zero weight to boundary nodes (these nodes are not 
c  included in the inversion, but are in the velocity array,
c  thus we want to avoid writing to negative or incorrect 
c  elements of the partial derivative matrix)
c
            if(ip.eq.1) then
c              write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(3)=0.0
               wv(5)=0.0
               wv(7)=0.0
            else
               if(ip.eq.nx1) then
c                 write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
                  wv(2)=0.0
                  wv(4)=0.0
                  wv(6)=0.0
                  wv(8)=0.0
               end if
            endif
c
            if(jp.eq.1) then
c              write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(2)=0.0
               wv(5)=0.0
               wv(6)=0.0
            else
               if(jp.eq.ny1) then
c                 write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
                  wv(3)=0.0
                  wv(4)=0.0
                  wv(7)=0.0
                  wv(8)=0.0
               endif
            endif
c
            if((kpg.eq.1).or.(kpg.eq.(nz1+1))) then
c              write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
               do 30 izg=1,4
                  wv(izg)=0.0
   30          continue
            else
               if((kpg.eq.nz1).or.(kpg.eq.(2*nz1))) then
c                 write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
                  do 35 izg=5,8
                     wv(izg)=0.0
   35             continue
               endif
            endif
 1610       format(' ASSIGNING ZERO WEIGHTS IN TTMDER',
     2         ' no=',i3,',xp=',f7.2,',yp=',f7.2,',zp=',
     3         f7.2,',v=',f5.3,',ip=',i2,',jp=',i2,',kpg=',i2)
c
c  Accumulate model partial derivatives
            do 48 kk=1,2
               kk1=kk-1
               do 47 jj=1,2
                  jj1=jj-1
                  do 46 ii=1,2
                     ii1=ii-1
                     ijk=ii+2*jj1+4*kk1
c skip boundary nodes
                     if(wv(ijk).lt.0.05) goto 46
c  skip fixed nodes
                     if(nfix(in(ijk)).eq.1) goto 46
                     ini=ndexfx(in(ijk))
c check for writing to an array element that is outside of the inversion
c  solution array
                     if((ini.lt.1).or.(ini.gt.nparvi)) then
                        write(16,1606) ini,ijk
 1606                   format(' *** Error in TTMDER, accessing',
     2                     ' gridpoint outside of velocity inversion',
     3                     ' gridpoints, ini=',i5,', ijk=',i5,/,
     4                     22x,'Probably boundary gridpoints are',
     5                     ' too close (TTMDER tries to write DTM',
     6                     ' elements with wv >= 0.05)')
                        write(16,1603) ne,no,xp,yp,zp,v,ip,jp,kp,kpg
 1603                   format(' ne=',i5,', no=',i5,', xp=',f8.2,
     2                     ', yp=',f8.2,', zp=',f8.2,', v=',f8.3,/,
     3                     21x,'ip=',i6,',   jp=',i6,',   kp=',i6,
     4                     ',   kpg=',i6)
                        write(16,1607) (j,in(j),j,wv(j),j=1,8)
 1607                   format(' in(',i1,')=',i6,' wv(',i1,')=',e15.5)
                        write(16,1608)
 1608                   format(' * * * * STOP * * * * (to avoid',
     2                     ' writing outside of defined DTM array)')
                        stop
                     end if
                     inp=ini+(no-1)*npari
                     hit(in(ijk))=hit(in(ijk))+wv(ijk)
c** S-P CHANGE
c   set up equation to solve for delta-Vp/Vs
                     if (isp.eq.0) dtm(inp)=dtm(inp)+dt*wv(ijk)*
     &                  vel(ip+ii1,jp+jj1,kp+kk1)/v
                     if (isp.eq.1) then
                        isp=0
                        call vel3(isp,xp,yp,zp,vp)
                        dtm(inp)=dtm(inp)+ssl*wv(ijk)/vp
                        isp=1
                     endif
c***
   46             continue
   47          continue
   48       continue
c
   55    continue
   50 continue
c  arrival time residual
  88    continue
      is=isto(no,ne)
      res(no)=secp(no,ne)-seco(ne)-ttime-pdl(is)
c** S-P CHANGE
c**     Estimate the calculated P and S times and then the S-P time
c       Then estimate the residuals between the observed and calculated.
      if (isp.eq.1) then
      intsp(no,ne)=0
      call path(ne,no,xe,ye,ze,ptime)
      intsp(no,ne)=1
      smptime=ttime-ptime
      res(no)=secp(no,ne)-smptime - sdl(is)
      endif
c**
c
  93  continue
      if(invdel.eq.0)return
c  Calculate station delay partial derivatives (unless doing location only)
      if(kopt.eq.0) return
      if(nfixst(is).eq.1) return
      ncor=nparv+is
c  If the observation is an S-wave, put the derivative nsts further up.
      if(intsp(no,ne).ne.0)ncor=ncor+nsts
      ini=ndexfx(ncor)
c   The derivative array, dtm(), is a long, linear array with all
c   of the observations and their derivatives end to end.
c   hit is an array which is used for this observation only. Khit
c   is set in subroutine parsep.
c   Each observation occupies npar locations, which = 1+(no-1)*npar to
c   no*npar
c   npar=nodes2, or nodes2+2*nsts when stn delays are being inverted, too.
      inp=(no-1)*npari+ini
      dtm(inp)=1.0
c checking station corrections
c     write(16,1620) is,inp,dtm(inp)
c1620 format(' TTMDER:is,inp,dtm',i5,i10,f14.6)
      hit(ncor)=hit(ncor)+1.0
      return
c***** end of subroutine ttmder *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine vel3(isp,x,y,z,v)
c  This routine is Cliff Thurber's
c  common block variables:
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'simulps_common.inc'
c
c  use Prothero's intmap here
      call intmap(x,y,z,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
c	write(16,100)x,xl,ip,y,yl,jp,z,zl,kp
c100	format(3(2f7.3,i3))
      xf=(x-xn(ip))/(xn(ip1)-xn(ip))
      yf=(y-yn(jp))/(yn(jp1)-yn(jp))
      zf=(z-zn(kp))/(zn(kp1)-zn(kp))
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
      wv(1)=xf1*yf1*zf1
      wv(2)=xf*yf1*zf1
      wv(3)=xf1*yf*zf1
      wv(4)=xf*yf*zf1
      wv(5)=xf1*yf1*zf
      wv(6)=xf*yf1*zf
      wv(7)=xf1*yf*zf
      wv(8)=xf*yf*zf
c  calculate velocity
c  S-velocity is stored after P-velocity
      kpg=kp
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
      v=wv(1)*vel(ip,jp,kp)+wv(2)*vel(ip1,jp,kp)
     2 +wv(3)*vel(ip,jp1,kp)+wv(4)*vel(ip1,jp1,kp)
     * +wv(5)*vel(ip,jp,kp1)+wv(6)*vel(ip1,jp,kp1)
     * +wv(7)*vel(ip,jp1,kp1)+wv(8)*vel(ip1,jp1,kp1)
      return
c***** end of subroutine vel3 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine veladj(nit)
c  routine to set up matrices for velocity inversion
c
c  declaration statements:
      integer nm(4),mbla(3)
      real r2(4),rhm(4),rhsvar(4),sqnorm(4),snrm(4),c1(2),dampc1(2)
      character*11 solnam(4)
      parameter(zero=0.0)
c
c  common block variables:
      include 'simulps_common.inc'
c
      data solnam/'Velocity   ','P-Velocity ','Vp/Vs ratio',
     2 'Sta. Corr. '/
c
c  compute adjustments to velocity model
c
        write(16,385)
 385    format(/,' Compute velocity adjustments to model-Veladj')
c  remove unknowns with fewer than nmin observations
c  mbl was defined in PARSEP as no. rows (ie: no. velocity gridpts.
c    observed), in VELADJ, it's adjusted by the <hitct gridpts.
      kbl=mbl
      if (hitct.eq.0.0) go to 230  !skip next section if hitcutoff not used
c
      do 10 i=1,3
         mbla(i)=0
   10 continue
      mbl=0
      ii=0
      jj=0
c  construct inverse mapping for index
c  index points from nodes(input) to solution arrays (rhs,g), while 
c  jndex points in the other direction.
c       loop over all sampleable nodes.
c  twice as many nodes if both P and S velocity
      do 170 i=1,npari
c  g has no equations for Unsampled nodes
         if (khit(mdexfx(i)).eq.0) go to 170
         jndex(index(i))=i
  170 continue
c       now remove equations from g for poorly sampled nodes
c  jj=counter for old g array, ii=counter for new g array
      do 220 i=1,kbl
         if(khit(mdexfx(jndex(i))).eq.0) goto 180
         if(hit(mdexfx(jndex(i))).ge.hitct) goto 190
  180    jj=jj+i
         go to 220
  190    mbl=mbl+1
c  reorder rhs and index
         rhs(mbl)=rhs(i)
         index(jndex(i))=mbl
         do 210 j=1,i
            jj=jj+1
            if (khit(mdexfx(jndex(j))).eq.0) goto 210
            if(hit(mdexfx(jndex(j))).ge.hitct) goto 200
            go to 210
  200       ii=ii+1
            g(ii)=g(jj)
  210    continue
  220 continue
c  revise jndex array for reduced g array
      do 225 i=1,npari
         if(khit(mdexfx(i)).eq.0) goto 225
         if(hit(mdexfx(i)).lt.hitct) goto 225
         jndex(index(i))=i
  225 continue
c  add damping (different for Vp and Vs) to diagonal elements of g
  230 j=0
      do 240 i=1,mbl
         j=j+i
         kv=(mdexfx(jndex(i))-1)/nodes2+1
c  Correction for vp & sta corr only (w/out vp/vs)
         if((kv.eq.2).and.(iuses.eq.1)) kv=3
         if(kv.gt.3) kv=3
         g(j)=g(j)+vdamp(kv)
  240 continue
c  if desired store g and rhs for resolution and error calculation
      if (ires.eq.0) go to 290
c  store .rhs
      do 260 n=1,mbl
         rhs1(n)=rhs(n)
  260 continue
      ng=mbl*(mbl+1)/2
      do 270 n=1,ng
         g1(n)=g(n)
  270 continue
c  subtract off vdamp from diagonal of g1
      j=0
      do 280 i=1,mbl
         j=j+i
         kv=(mdexfx(jndex(i))-1)/nodes2+1
         if((kv.eq.2).and.(iuses.eq.1)) kv=3
         if(kv.gt.3) kv=3
         mbla(kv)=mbla(kv)+1
         g1(j)=g1(j)-vdamp(kv)
  280 continue
c
  290 ier=0
      write(16,20) kbl,mbl,(solnam(i+1),mbla(i),i=1,3)
  20  format(' total nodes observed=',i4,', nodes inverted for=',i4,
     2 3(',',1x,a11,':',i4))
c  perform lu decomposition
      call ludecp(g,g,mbl,ier)
      if (ier.ne.0) write(16,1001)
 1001 format('  *** error in ludecp ***')
      if(ier.eq.129) write(16,1610)
 1610 format('  ier=129, matrix A is algorithmically not positive ',
     2 'definite')
c  perform lu elimination
      call luelmp(g,rhs,mbl,rhs)
c  inversion complete
c
c  calculate solution vector statistics
c  mean of rhs
      ivadj=0
      if(nparvi.eq.0) goto 36
      do 22 i=1,4
         rhm(i)=0.0
         nm(i)=0
   22 continue
      do 23 n=1,npari
         if(khit(mdexfx(n)).eq.0) goto 23
         if(hit(mdexfx(n)).lt.hitct) goto 23
         k=(mdexfx(n)-1)/nxy2 + 2
         ia=k/nz + 2
         if((ia.eq.3).and.(iuses.eq.1)) ia=4
         if(ia.gt.4) ia=4
         rhm(ia)=rhm(ia)+rhs(index(n))
         nm(ia)=nm(ia)+1
         if(ia.gt.3) goto 23
         rhm(1)=rhm(1)+rhs(index(n))
         nm(1)=nm(1)+1
   23 continue
      do 24 i=1,4
         if(nm(i).gt.0.0) rhm(i)=rhm(i)/float(nm(i))
         rhsvar(i)=0.0
         sqnorm(i)=0.0
   24 continue
      do 25 n=1,npari
         if(khit(mdexfx(n)).eq.0) goto 25
         if(hit(mdexfx(n)).lt.hitct) goto 25
         k=(mdexfx(n)-1)/nxy2 + 2
         ia=k/nz + 2
         if((ia.eq.3).and.(iuses.eq.1)) ia=4
         if(ia.gt.4) ia=4
         sqnorm(ia)=sqnorm(ia)+rhs(index(n))*rhs(index(n))
         r2(ia)=rhs(index(n))-rhm(ia)
         rhsvar(ia)=rhsvar(ia)+r2(ia)*r2(ia)
         if(ia.gt.3) goto 25
         sqnorm(1)=sqnorm(1)+rhs(index(n))*rhs(index(n))
         r2(1)=rhs(index(n))-rhm(1)
         rhsvar(1)=rhsvar(1)+r2(1)*r2(1)
   25 continue
      write(16,1624) nit
      write(36,1624) nit
 1624 format(' iteration',i3,' rhs solution vector:')
      do 27 i=1,4
         if(nm(i).le.0) goto 27 
         snrm(i)=sqrt(sqnorm(i))
         rhsvar(i)=rhsvar(i)/float(nm(i))
         write(16,1625) solnam(i),rhm(i),rhsvar(i),sqnorm(i),snrm(i)
         write(36,1625) solnam(i),rhm(i),rhsvar(i),sqnorm(i),snrm(i)
   27 continue
 1625 format(9x,a11,': mean',f13.9,', variance',
     2 f13.9,', norm squared',f13.9,', norm',f13.9)
c check whether velocity solution norm is below cutoff value
      if(snrm(1).lt.snrmct) then
        nit=99
        write(16,1630) snrm(1),snrmct
 1630   format(/,' **** STOP since snrm',f13.9,' is below',
     2  ' snrmct',f13.9,' ****',/)
        return
      endif
c
c  apply adjustments to velocity model
      do 30 n=1,nparvi
c  gridpoint number in array of all nodes
         nn=mdexfx(n)
         vadj(nn)=zero
         if(khit(nn).eq.0) goto 30
         if(hit(nn).lt.hitct) goto 30
         rh=rhs(index(n))
c  calculate x and z indices of velocity grid
         k=(nn-1)/nxy2+2
         j=2+(nn-1+(2-k)*nxy2)/nx2
         i=1+nn+nx2*(2-j)+nxy2*(2-k)
         if(k.ge.nz) k=k+2      ! if s velocity node
         delm=-vel(i,j,k)*rh/(1.0+rh)
c**   Modification for S-P times
c    equations solve for delta-Vp/Vs directly
         if(k.ge.nz) delm=rh
         delma=abs(delm)
c  place upper bound on velocity adjustment
         if(k.gt.nz) then
            dvmax=dvsmx
         else
            dvmax=dvpmx
         end if
         if (delma.gt.dvmax) delm=dvmax*delm/delma
         vadj(nn)=delm
c  apply adjustment to model
c** S-P CHANGE
         if (k.le.nz) vel(i,j,k)=vel(i,j,k)+delm
         if(k.ge.nz) then
            vpvs(i,j,k-nz)=vpvs(i,j,k-nz)+delm
         endif
         ivadj=ivadj+1
   30 continue
c
      if(iuses.eq.1) goto 130
c  Compute Vs with new Vp and new Vp/Vs
        do 120 k=1,nz
           ks=k+nz
           do 115 j=1,ny
              do 110 i=1,nx
                 vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
  110         continue
  115      continue
  120   continue
c**
c
c Recalculate damping
c consider c1 to give relative size of 2 terms being minimized
c  min [ sum(data resid squared) + damping*(solution norm squared)]
  130 do 34 i=1,iuses
         if(i.eq.1) ssq=ssqrwp
         if(i.eq.2) ssq=ssqrws
         if(ssq.le.0) goto 34
         if(mbla(i).eq.0) goto 34
         if(nit.eq.1) c1(i)=vdamp(i)*sqnorm(i+1)/ssq
         if(nit.gt.1) dampc1(i)=c1(i)*ssq/sqnorm(i+1)
         write(16,1626) solnam(i+1),vdamp(i),c1(i),dampc1(i)
 1626    format('Damping: ',a11,' damp =',f8.2,', c1 =',f8.4,
     2      ', damp(c1) =',f8.2)
c Replace damping if idmp is turned on
c Do not change damping to be smaller
         if(dampc1(i).lt.vdamp(i)) goto 34
         if((idmp.eq.1).and.(nit.gt.1)) vdamp(i)=dampc1(i)
   34 continue
   36 continue
      isadj=0
      if(invdel.eq.0)goto 99
c  Station correction adjustments
      do 40 n=1,nsts
         nvn=n+nparv
         if(khit(nvn).eq.0) goto 40
         if((hit(nvn).lt.hitct).or.(nfixst(n).eq.1)) goto 40
         nv=ndexfx(nvn)
         vadj(nvn)=zero
         nin=index(nv)
         vadj(nvn)=rhs(nin)
         pdl(n)=pdl(n)+vadj(nvn)
         isadj=isadj+1
   40 continue
      if(iuses.eq.1) goto 99
c  S delay adjustments
      do 45 n=1,nsts
         nvn=n+nsts+nparv
         if(khit(nvn).eq.0) goto 45
         if((hit(nvn).lt.hitct).or.(nfixst(n).eq.1)) goto 45
         nv=ndexfx(nvn)
         vadj(nvn)=zero
         nin=index(nv)
         vadj(nvn)=rhs(nin)
         sdl(n)=sdl(n)+vadj(nvn)
         isadj=isadj+1
   45 continue
   99 continue
        write(16,1602) ivadj,isadj
 1602 format(i10,' velocity adjustments,',i10,' station corr. '
     2 , 'adjustments')
c
      return
c***** end of subroutine veladj *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine velbku(istop1)
c  This routine will backup the velocity and station parameters halfway
c  common block variables:
      include 'simulps_common.inc'
c
c  apply adjustments to velocity model
        ivadj=0
      if(nparvi.eq.0) goto 36
      do 30 n=1,nparvi
c  gridpoint number in array of all nodes
         nn=mdexfx(n)
         if(vadj(nn).eq.0.0) goto 30
c  calculate x and z indices of velocity grid
         k=(nn-1)/nxy2+2
         j=2+(nn-1+(2-k)*nxy2)/nx2
         i=1+nn+nx2*(2-j)+nxy2*(2-k)
         if(k.ge.nz) k=k+2   ! if S velocity node
c  compute half adjustment if first backup
         if(istop1.eq.0) vadj(nn)= -0.5*vadj(nn)
c  apply adjustment to model
c** S-P CHANGE
         if (k.le.nz) vel(i,j,k)=vel(i,j,k)+vadj(nn)
         if(k.ge.nz) then
            vpvs(i,j,k-nz)=vpvs(i,j,k-nz)+vadj(nn)
         endif
         ivadj=ivadj+1
   30 continue
c
      if(iuses.eq.1) goto 36
c  Compute Vs with new Vp and new Vp/Vs
        do 120 k=1,nz
           ks=k+nz
           do 115 j=1,ny
              do 110 i=1,nx
                 vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
  110         continue
  115      continue
  120   continue
c**
c
   36 continue
        isadj=0
      if(invdel.eq.0)goto 99
c  Station correction adjustments
      do 40 n=1,nsts
         nvn=n+nparv
         if(khit(nvn).eq.0) goto 40
         if((hit(nvn).lt.hitct).or.(nfixst(n).eq.1)) goto 40
         if(istop1.eq.0) vadj(nvn)= -0.5*vadj(nvn)
         pdl(n)=pdl(n)+vadj(nvn)
         isadj=isadj+1
   40 continue
      if(iuses.eq.0) goto 99
c  S delay adjustments
      do 45 n=1,nsts
         nvn=n+nsts+nparv
         if(khit(nvn).eq.0) goto 45
         if((hit(nvn).lt.hitct).or.(nfixst(n).eq.1)) goto 45
         if(istop1.eq.0) vadj(nvn)= -0.5*vadj(nvn)
         sdl(n)=sdl(n)+vadj(nvn)
         isadj=isadj+1
   45 continue
   99 continue
        write(16,1602) ivadj,isadj
 1602 format(i10,' velocity adjustments,',i10,'station corr. '
     2  , 'adjustments')
c
      return
c***** end of subroutine velbku *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine veld(isp,xx,yy,zz,vx,vy,vz)
c
c*****this routine computes the derivatives of velocity
c     in x, y, and z directions
c
c  common block variables:
      include 'simulps_common.inc'
c
c  use Prothero's intmap here
      call intmap(xx,yy,zz,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
      xd=xn(ip1)-xn(ip)
      yd=yn(jp1)-yn(jp)
      zd=zn(kp1)-zn(kp)
c
      xf=(xx-xn(ip))/xd
      yf=(yy-yn(jp))/yd
      zf=(zz-zn(kp))/zd
c
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
c  S-velocity is stored in the 2nd half of the velocity array
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
c
c*****calculate derivatives of velocity
c
      vx=(yf1*zf1*(vel(ip1,jp,kp)-vel(ip,jp,kp))
     *+yf*zf1*(vel(ip1,jp1,kp)-vel(ip,jp1,kp))
     *+yf1*zf*(vel(ip1,jp,kp1)-vel(ip,jp,kp1))
     *+yf*zf*(vel(ip1,jp1,kp1)-vel(ip,jp1,kp1)))/xd
c
      vy=(xf1*zf1*(vel(ip,jp1,kp)-vel(ip,jp,kp))
     *+xf*zf1*(vel(ip1,jp1,kp)-vel(ip1,jp,kp))
     *+xf1*zf*(vel(ip,jp1,kp1)-vel(ip,jp,kp1))
     *+xf*zf*(vel(ip1,jp1,kp1)-vel(ip1,jp,kp1)))/yd
c
      vz=(xf1*yf1*(vel(ip,jp,kp1)-vel(ip,jp,kp))
     *+xf*yf1*(vel(ip1,jp,kp1)-vel(ip1,jp,kp))
     *+xf1*yf*(vel(ip,jp1,kp1)-vel(ip,jp1,kp))
     *+xf*yf*(vel(ip1,jp1,kp1)-vel(ip1,jp1,kp)))/zd
c
      return
c ***** end of subroutine veld *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine wthyp(ne,nobs,ah,rmswt,nwr,w)
c  routine to apply weighting to hypocenter matrix + residuals
c
c  common block variables:
      include 'simulps_common.inc'
c
c  declaration statements:
      real w(maxobs)
      double precision ah(maxobs,maxobs)
      character*1 phs(2)
      parameter (pi=3.1415927)
c
      data phs/'P','S'/
c
c  compute weighting
c
      nwr=0
      wnorm=0.0
      do 13 j=1,nobs
c  reading weight
         w(j)=wt(j,ne)
c  change relative weighting of S-P observations
         if(intsp(j,ne).eq.1) w(j)=w(j)*wtsp
c  residual weighting
c  downweighting(linear) 0 to 98% res1 to res2, 98 to 100% res2 to res3
         ares=abs(res(j))
         if(ares.le.res2) then
            wr=1.0-(ares-res1)*dres12
            if (wr.gt.1.0) wr=1.0
         else
            if(res3.gt.res2) then
               wr=0.02-(ares-res2)*dres23
               if (wr.lt.0.0) wr=0.0
            else
               wr=0.0
            endif
         endif
c  distance weighting
         wd=1.0-(dlta(j,ne)-delt1)*ddlt
         if (wd.gt.1.0) wd=1.0
         if (wd.lt.0.0) wd=0.0
c  unnormalized weight
         w(j)=w(j)*wr*wd
         wnorm=wnorm+w(j)
         if (w(j).gt.0.0) nwr=nwr+1
  13  continue
c  check to be sure 4 or more readings with nonzero weight
      if(nwr.lt.4) goto 900
c  normalize weights
   12 wfac=nwr/wnorm
c
   20 continue
      wnorm=0.0
      rmswt=0.0
c  normalize weights, apply to hypocenter matrix, and place into
c  a-matrix for inversion
      do 30 i=1,nobs
         w(i)=w(i)*wfac
         wnorm=wnorm+w(i)*w(i)
         rmswt=rmswt+(res(i)*w(i))*(res(i)*w(i))
         do 35 j=1,4
            ah(i,j)=dth(i,j)*w(i)
   35    continue
         ah(i,5)=res(i)*w(i)
   30 continue
      rmswt=sqrt(rmswt/wnorm)
      return
  900 write(16,1601) ne
 1601 format(' in WTHYP, event',i5,' HAS LESS THAN 4 READINGS WITH ',
     2 '  NON-ZERO WEIGHT, SHOULD BE DELETED')
      write(16,51)
51    format(1x,4('sta ph   wt   res   ttime delta',2x))
      write(16,53) (stn(isto(j,ne)),phs(intsp(j,ne)+1),w(j),res(j),
     2 (secp(j,ne)-seco(ne)),dlta(j,ne),j=1,kobs(ne))
53    format(4(1x,a4,1x,a1,f6.2,f7.3,2f6.2,1x))
      write(16,1602)
 1602 format(/)
      return
c***** end of subroutine wthyp *****
      end
