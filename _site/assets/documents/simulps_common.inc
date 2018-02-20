c23456789a123456789b123456789c123456789d123456789e123456789f123456789g12
c   common blocks for program SIMULPS12
c
c   10-mar-1994
c
c
c  definitions:
c   nsts = number of stations
c   neqs = number of earthquakes
c   nbls = number of blasts
c   nsht = number of shots
c   nobs = number of observations for a single event.
c   nx  = number of nodes in the x direction.
c   ny = number of velocity nodes in y direction.
c   nz = number of velocity nodes in z direction.
c   nodes2=(nx-2)(ny-2)(nz-2) - number of nodes which are inverted for.
c   nparv=number of velocity parameters, 2*nodes2 if both P and S
c         otherwise=nodes2
c   nparvi=number of velocity parameters allowd to invert for
c   npar = nparv ( station delays not included)
c         nparv + nsts  (station delays,P only)
c         nparv + 2*nsts (station delays,P and S)
c   npari = nparvi (station delays not included)
c           nparvi + nsts (station delays, P only)
c           nparvi + 2*nsts (station delays,P and S)
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  TO CHANGE ARRAY SIZES,
C  PARAMETERS IN THIS SECTION CAN BE CHANGED
c
c   maxev: maximum number of events (earthquakes+blasts+shots)
      parameter(maxev=4000)
c
c   maxobs: maximum number of arrivals for a particular event.  Note
c          that maxobs can generally be set much lower than maxsta
c          since events will not generally be recorded at all
c          stations.
      parameter(maxobs=100)
c
c   maxsta: maximum number of stations in station list
      parameter(maxsta=100)
c
c   mxpari: maximum number of parameters to invert for.  This 
c           includes Vp and Vp/Vs nodes and P and S station 
c           corrections.  mxpari will generally be much lower
c           than maxpar because velocity nodes and station
c           corrections can be fixed. (Station corrections
c           are generally fixed to zero in the region of the
c           velocity inversion.)
      parameter(mxpari=20000)
c
c   maxpar: maximum number of potential parameters that could
c           be included in the inversion.  This includes all 
c           velocity nodes (includes both Vp and Vp/Vs) that
c           define the 3-D velocity model (nparv). Does not 
c           include edge nodes.
c           If station corrections are included in the inversion,
c           maxpar also includes stations (nparv+nsts or 
c           nparv + 2*nsts if P and S-P are used).
c           maxpar is set much larger than mxpari to allow
c           for assigning regional velocity (from prior
c           coarse inversions) values (fixed) in peripheral
c           areas where there is nil resolution.
      parameter(maxpar=25000)
c
c   maxnx: maximum number of nodes in x direction
      parameter(maxnx=24)
c   maxny: maximum number of nodes in y direction
      parameter(maxny=24)     
c   maxnz: maximum number of nodes in z direction
      parameter(maxnz=15)     
c
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c  common block variables:
c
c   control parameters
c
      common/contrl/ nsts,neqs,nobt,wnobt,ires,nitloc,eigtol,rmscut,
     * dvpmx,dvsmx,hitct,nsht,nbls,nevt,i3d,nitmax,snrmct,ihomo,rmstop
     * ,delt1,delt2,ddlt,res1,res2,res3,dres12,dres23,
     * stepl,wtsht,wtsp,kout,kout2,kout3
     * ,ndip,iskip,scale1,scale2,iusep,iuses,invdel,zmin,dxmax,idmp
     * ,vdamp(3),ifixl,xfac,tlim,nitpb(2),nobtp,nobts,nobtex,nobteq
     * ,nwrt,nswrt,rderr,ercof,wnobtp,wnobts
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   event observation and station arrays
      common/events/ mino(maxev),seco(maxev),ltde(maxev),eltm(maxev),
     *  ihr(maxev),lnde(maxev),elnm(maxev),evc(3,maxev),kobs(maxev),
     *  iyrmo(maxev),iday(maxev),rmag(maxev),wrmsr(maxev)
      common/ntemp/ netemp,nbtemp,neb
      common/sepdat/ totrms(maxev),seprms(maxev)
      common/fastrc/ jfl,nco(maxobs,maxev),ndo(maxobs,maxev)
c
      common/observ/ stn(maxsta),ltds(maxsta),sltm(maxsta),lnds(maxsta),
     *  slnm(maxsta),stc(3,maxsta),isto(maxobs,maxev),wt(maxobs,maxev),
     *  secp(maxobs,maxev),intsp(maxobs,maxev),pdl(maxsta),sdl(maxsta),
     *  nfixst(maxsta)
      common/obsiw/ iw(maxobs,maxev)
      common/nrdsta/nrd(maxsta,2)
      common/fpsdat/ az(maxobs),toa(maxobs)
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   velocity arrays
c
c   /vmod3d/
c    vel(nx,ny,2*nz)     ! 1st half is P velocity, 2nd half is S velocity
      parameter(maxnz2=maxnz*2)
c
      common/vmod3d/bld,nx,ny,nz,xn(maxnx),yn(maxny),zn(maxnz),
     *  vel(maxnx,maxny,maxnz2),vpvs(maxnx,maxny,maxnz)
c
c   /locate/
c   ixkms: size of map in x direction.  The distance between 
c            furthest nodes.  This will be in units of km if bld=1.0 
c            or in units of 0.1km if bld=0.1
      parameter(ixkms=1500)
c   iykms=size of map in y dir.
      parameter(iykms=1500)
c   izkms=size of map in z dir.
      parameter(izkms=1500)
c
      common/locate/ xl,yl,zl,ixloc(ixkms),iyloc(iykms),izloc(izkms)
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  inversion arrays
c
      common/hypinv/ res(maxobs),dth(maxobs,4),dthp(maxobs,4),
     *  dlta(maxobs,maxev)
c
c   /modinv/
      parameter(mxdtm=maxobs*mxpari)
      common/modinv/ dtm(mxdtm),dtmp(mxdtm),resp(maxobs),nz2,nx1,
     *  ny1,nz1,nparvi,npari,inf,nfix(maxpar),ndexfx(maxpar),
     *  mdexfx(mxpari),nparv,npar,nodes,nxy,hit(maxpar),
     *  vadj(maxpar),khit(maxpar),nx2,nxy2,nodes2,sumhit(maxnz2)
c
c   /solutn/
      parameter(mgsol=(mxpari*(mxpari+1))/2)
      parameter(mxpri1=mxpari+1)
      common/solutn/ g(mgsol),g1(mgsol),rhs(mxpri1),rhs1(mxpri1),
     *  index(mxpari),jndex(mxpari),mbl,nbtot,ssqr,mbl1,ssqr1,
     *  ssqrw,ssqrw1,ssqrwp,ssqrws,var,varw,var1,varw1,
     *  ndof,ndof1,wndof,wndof1
c
      common/resltn/ covdi(maxpar),drm(maxpar),stderr(maxpar)
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
c  DIMENSIONS FOR RAY TRACING
c   
c    NOTE: This is not designed to be changed by the user.
c          There are other dimension and commons throughout
c          the program besides /rpath/.  /rpath/ is just
c          included here because it also uses maxobs to
c          dimension the array.
c
c    nseg = number of segments in longest ray. 130 is used as
c           programmed by Cliff Thurber.
c
c     nseg=2**(1+nint(3.32193*alog(sep/scale1))). The
c     maximum number of segments allowed (see subroutine setup)
c     is 2**7.
c
c     ncr=1+.5*sep/scale2   - - the maximum number of curves
c                        for the circular rays.
c
c     sep=max dist between event and station. Scale1 is set in
c     input parameters. So is scale2.
c
c    (3*nseg) : strpth,fstpth,trpth1
c    (nseg) : pthsep
c    (3,9) : dipvec
c    (3*nseg,9) : disvec,trpath
c    (9) : trtime
c   /rpath/
c   (3,nseg,nobs) : rp
c   (nobs) : nrp,pl,tt
c
      common/rpath/rp(3,130,maxobs),nrp(maxobs),ttc(maxobs),
     *  xmov,pl(maxobs)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
