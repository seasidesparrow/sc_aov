      program aovsimple
      implicit none

      character*80 :: infile,outfile
      integer (kind=2) :: ifile
      integer (kind=4) :: ntot,i
      integer (kind=4), parameter :: ndatamax=10000000
      real (kind=8), allocatable :: tread(:),xread(:)
      real (kind=8), allocatable :: tdata(:),xdata(:)
      real (kind=8) :: xave

      real (kind=8) :: fmaxout,thetamaxout
      
      write(6,100)
100   format(32hPlease enter an input file name:)
      read*,infile
      ifile=index(infile,' ')
      outfile=infile(1:ifile-1)//'.aov'

      open(unit=1,file=infile,status='old',err=10)
      goto 11
10    stop 'Could not open input file.'
11    open(unit=2,file=outfile,status='unknown',err=12)
      goto 13
12    stop 'Could not open output file.'
13    continue

      allocate(tread(ndatamax),xread(ndatamax))
      ntot=0
      do i=1,ndatamax
       read(1,*,end=99) tread(i),xread(i)
       ntot=ntot+1
      enddo
99    close(1)

      allocate(tdata(ntot),xdata(ntot))
      tdata=tread(1:ntot)
      xdata=xread(1:ntot)
      deallocate(tread,xread)
      xave=sum(xdata)/dfloat(ntot)

      call aovmax(ntot,tdata,xdata,xave,fmaxout,thetamaxout)
      print*,infile,fmaxout,thetamaxout
      end



      subroutine aovmax(ntot,tdata,xdata,xave,fmaxout,thetamaxout)
!------------------------------------------------------------------------------
! AOVMAX: Return the frequency at which theta(AOV) is maximum
!         Call this subroutine with the following variables:
!              ntot = total number of data points
!              tdata = time values of the data points
!              xdata = magnitude, flux, etc
!              xave  = average of xdata
!
!         Return values:
!              fmaxout = frequency where theta(AOV) is maximum
!              thetamaxout = value of theta(AOV) at maximum
!
! Program by M. Templeton, August 15, 2007
! Adapted from A. Schwarzenberg-Czerny 1989 (MNRAS 241,153)
!------------------------------------------------------------------------------
      implicit none

      integer (kind=4) :: ntot,i,j,k
      real (kind=8) :: tdata(ntot),xdata(ntot),xave
      real (kind=8) :: fmaxout,thetamaxout

      real (kind=8) :: freqmin,freqmax,deltafreq,pfold
      integer (kind=4) :: nfreq

      integer (kind=4) :: rbin,ibin
      real (kind=8) :: oversample
      integer (kind=4), allocatable :: nbin(:)
      real (kind=8), allocatable :: tbin(:),xbin(:,:)
      real (kind=8), allocatable :: xbinave(:)
      real (kind=8) :: s0sq,s1sq,s2sq
      real (kind=8), allocatable :: thetaaov(:),frequency(:)
      integer (kind=4) :: im(1)
      
      rbin=7
      oversample=16.d0

      freqmin=1.d-2
      freqmax=5.d-1
      deltafreq=2.d0/(maxval(tdata)-minval(tdata))
      deltafreq=deltafreq/oversample
      nfreq=1+(freqmax-freqmin)/deltafreq
      allocate(frequency(nfreq),thetaaov(nfreq))
      
      do i=1,nfreq
       allocate(tbin(ntot),xbin(ntot,rbin),nbin(rbin))
       tbin=0.d0
       xbin=0.d0
       nbin=0
       pfold=1.d0/(freqmin+(deltafreq*dfloat(i-1)))
       tbin=tdata/pfold
       tbin=mod(tbin,1.d0)
       allocate(xbinave(rbin))
       do j=1,ntot
        if(tbin(j).lt.0.d0) tbin(j)=tbin(j)+1.d0
        ibin=idint(rbin*tbin(j))+1
        nbin(ibin)=nbin(ibin)+1
        xbin(nbin(ibin),ibin)=xbin(nbin(ibin),ibin)+xdata(j)
       enddo
       s1sq=0.d0
       do j=1,rbin
        xbinave(j)=(sum(xbin(1:nbin(j),j)))/dfloat(nbin(j))
        s1sq=s1sq+nbin(j)*(xbinave(j)-xave)**2
        s0sq=0.d0
        s2sq=0.d0
        do k=1,nbin(j)
         s2sq=s2sq+(xbin(k,j)-xbinave(j))**2
         s0sq=s0sq+(xbin(k,j)-xave)**2
        enddo
       enddo
       s1sq=s1sq/dfloat(rbin-1)
       s2sq=s2sq/dfloat(ntot-rbin)
       s0sq=s0sq/dfloat(ntot-1)
       thetaaov(i)=s1sq/s2sq
       frequency(i)=1.d0/pfold
       deallocate(tbin,xbin,nbin,xbinave)
      enddo
      im=maxloc(thetaaov)
      fmaxout=frequency(im(1))
      thetamaxout=thetaaov(im(1))
      return

      end
