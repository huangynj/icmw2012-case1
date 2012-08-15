      program kinematic
C
C AUTHORS: WOJCIECH GRABOWSKI (grabow@ncar.ucar.edu, 303-497-8974)
C      with 2D MPDATA routine written by P. Smolarkiewicz (NCAR)
C
C
c This program can be used in the 2D kinematic tests of the warm
c rain microphysics. The airflow is prescribed using streamfunction
c subroutine that calculates anelastic flow pattern at every
c timestep. Standard stagerred C-grid is assumed, i.e., velocities
c are shifted half the distance from the thermodynamic fields. The
c diagram below shows a sample of the grid for the 4 by 4 setup;
c Xs mark positions of thermodynamic fields, Ws and Us - positions
c of velocity components, and 'o's mark positions of the stream-
c function that is used to calculate nondivergent velocities.
c Note that the grid requires 5 by 4 horizontal velocity points,
c 4 by 5 vertical velocity points, and 5 by 5 streamfunction points.
c
c          o----W----o----W----o----W----o----W----o
c          |         |         |         |         |
c          U    X    U    X    U    X    U    X    U
c          |         |         |         |         |
c          o----W----o----W----o----W----o----W----o
c          |         |         |         |         |
c          U    X    U    X    U    X    U    X    U
c          |         |         |         |         |
c          o----W----o----W----o----W----o----W----o
c          |         |         |         |         |
c          U    X    U    X    U    X    U    X    U
c          |         |         |         |         |
c          o----W----o----W----o----W----o----W----o
c          |         |         |         |         |
c          U    X    U    X    U    X    U    X    U
c          |         |         |         |         |
c          o----W----o----W----o----W----o----W----o
c
c
c
c  Bulk model is used in this example. If more sophisticated
c  thermodynamics is to be used (e.g., detailed microphysics)
c  changes to this code have to be introduced. For example, 
c  in a detailed microphysics case, every class of cloud
c  droplets and raindrops has to be advected independently. Also,
c  terminal velocities have to be added to the airflow velo-
c  cities for every class to account for particle sedimentation. 
c  Once the logic of this test is understood, these changes 
c  should be simple to make.
c
c  NOTE: there is considerable number of comments both in the main
c        program and in subroutines that should help to understand
c        details of this program 
c
c  THIS CODE HAS BEEN WRITTEN WITHOUT EFFICIENCY CONSIDERATIONS;
c  IT CAN BE SIGNIFICANTLY ENHANCED IF REQUIRED.
c
c
c basic parameters of the model
      parameter(nx=75,nz=75,dx=20.,dz=20.)
      parameter(nxp=nx+1,nxm=nx-1,nzp=nz+1,nzm=nz-1)
c
c potential temperature, water vapor mixing ratio, cloud water
c mixing ratio, rain water mixing ratio:
c   IN THE DETAILED CASE: qc and qr should be replaced with
c   qq(nx,nz,ncl) where ncl is number of droplet classes considered 
      dimension theta(nx,nz),qv(nx,nz),qc(nx,nz),qr(nx,nz)
c
c surface precipitation rate (used in analysis of model data)
      dimension sprec(nx)
c Jacobian used in advection (required in runs with topo, here
c redundant)
      dimension gac(nx,nz)
c
c velocity fields needed for advection:
c uxa,uza are air velocities, uzam is used to include sedimentation
c uz for printout
      dimension uxa(nxp,nz),uza(nx,nzp),uzam(nx,nzp),uz(nx,nz)

c velocity for ploting
      dimension uxpl(nx,nz),uzpl(nx,nz)

c
c environmental temperature, potential temperature, pressure and
c water vapor mixing ratio
      common /environ2/ temp_e(nz),theta_e(nz),pres_e(nz),
     1                  qv_e(nz),qc_e(nz)
c base state density profile
      common /environ1/ rho(nz)

ccc
c time step used in the model, number of time steps
      data dt,ntstp /4.,3600/
c
      character*50 lhead
c statement function to define rain terminal velocity:
      vterm(qq,rro)=36.34*sqrt(1./rro)*(rro*1.e-3*qq)**.1346
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    
      call opngks
      call gsclip(0)
      xl=nx*dx*1.e-3
      zl=nz*dz*1.e-3
C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    
C
cc parameters for the MPDATA routine:
      iord=2
      isor=1
      nonos=0
      idiv=0

c INSERT INITIAL FIELDS AND WRITE TO HISTORY TAPE (fort.22):
c -------> prescribe initial fields
        call init(nz,dz)
          do k=1,nz
           do i=1,nx
            theta(i,k)=theta_e(k)
            qv(i,k)   =qv_e(k)
            qc(i,k)   =qc_e(k)
            qr(i,k)   =0.
            gac(i,k)  =rho(k)
           enddo
          enddo 
c -------> call microphysical adjustement routine:

C      write(22) temp_e,theta_e,pres_e,qv_e,qc_e
C      write(22) rho
C
C      write(22) theta
C      write(22) qv
C      write(22) qc
C      write(22) qr
C      write(22) qr  !<-- qr in place of updraft

ccc plot of initial conditions:
C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    
      time=0.
c   theta
      call setusv('LW',3000)
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call perim (3,5,3,5)
      cmx=300.
      cmn=260.
      cnt=.2
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(theta,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc negative values
      ct=-cnt
      call cpcnrc(theta,nx,nx,nz,cmn,ct,cnt,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.50,0.85, 'potential temperature (K)', 0.016,0.,0)
      CALL plchhq(.50,0.16, 'x (km)', 0.016,0.,0)
      CALL plchhq(.14,0.50, 'z (km)', 0.016,90.,0)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.20,0.16, '0.0', 0.016,0.,0)
      CALL plchhq(.80,0.16, '1.5', 0.016,0.,0)
      CALL plchhq(.15,0.20, '0.0', 0.016,0.,0)
      CALL plchhq(.15,0.80, '1.5', 0.016,0.,0)
      call frame

cc   qv
c      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
c      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
c      call periml(3,5,3,5)
c      cmx=20.*1.e-3
c      cmn=2.*1.e-3
c      cnt=0.25*1.e-3
ccc positive values
c      ct=cnt
c      call cpsetc('ILT',' ')
c      call cpcnrc(qv,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
ccc negative values
c      ct=-cnt
c      call cpcnrc(qv,nx,nx,nz,cmn,ct,cnt,-1,-1,-682)
      write(lhead,351) time
  351 format(' QV (kg/kg)         ',f6.0) 
c      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c      CALL plchmq(.50,0.85, lhead(1:50), 0.016,0.,0)
c      CALL plchhq(.50,0.12, 'X (km)', 0.016,0.,0)
c      CALL plchhq(.12,0.50, 'Z (km)', 0.016,90.,0)
c      call frame

c   qc
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call perim (3,5,3,5)
      cmx=10.*1.e-3
      cmn=.1 *1.e-3
      cnt=cmn
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(qc,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc trace:
      call cpcnrc(qc,nx,nx,nz,1.e-5,1.1e-5,1.e-5,-1,-1,682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.50,0.85, 'cloud water (g/kg)', 0.016,0.,0)
      CALL plchhq(.50,0.16, 'x (km)', 0.016,0.,0)
      CALL plchhq(.14,0.50, 'z (km)', 0.016,90.,0)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.20,0.16, '0.0', 0.016,0.,0)
      CALL plchhq(.80,0.16, '1.5', 0.016,0.,0)
      CALL plchhq(.15,0.20, '0.0', 0.016,0.,0)
      CALL plchhq(.15,0.80, '1.5', 0.016,0.,0)
      call frame

cc   qr
c      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
c      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
c      call perim (3,5,3,5)
c      cmx=10.*1.e-3
c      cmn=.002*1.e-3
c      cnt=cmn
ccc positive values
c      ct=cnt
c      call cpsetc('ILT',' ')
c      call cpcnrc(qr,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
ccc trace:
c      call cpcnrc(qr,nx,nx,nz,1.e-6,1.1e-6,1.e-6,-1,-1,682)
      write(lhead,353) time
  353 format(' QR (kg/kg)         ',f6.0) 
c      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c      CALL plchmq(.50,0.85, lhead(1:50), 0.016,0.,0)
c      CALL plchhq(.50,0.12, 'x (km)', 0.016,0.,0)
c      CALL plchhq(.12,0.50, 'z (km)', 0.016,90.,0)
c      call frame
C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    


c -------> get advective Courant numbers:
cc    FLOW DOES NOT CHANGE IN TIME....
        call adv_vel(uxa,uza,nx,nz,dx,dz,dt)

cc plot of the flow field:
C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    
       ax1=-1.e7
       ax2=-1.e7
       an1= 1.e7
       an2= 1.e7
       do i=1,nx
       do k=1,nz
       uxpl(i,k)=.5*(uxa(i,k)+uxa(i+1,k))   /dt*dx  /rho(k)
       uzpl(i,k)=.5*(uza(i,k)+uza(i,k+1))   /dt*dz  /rho(k)
       ax1=amax1(ax1,uxpl(i,k))
       ax2=amax1(ax2,uzpl(i,k))
       an1=amin1(an1,uxpl(i,k))
       an2=amin1(an2,uzpl(i,k))
       enddo
       enddo
        print*,'ux: min,max: ',an1,ax1
        print*,'uz: min,max: ',an2,ax2
c   u
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call perim (3,5,3,5)
cc positive values
      call cpsetc('ILT',' ')
      call cpcnrc(uxpl,nx,nx,nz,0.1,5.,0.1,-1,-1,1)
cc negative values
      call cpcnrc(uxpl,nx,nx,nz,-5.,-0.1,0.1,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.50,0.85,'horizontal flow', 0.016,0.,0)
      CALL plchhq(.50,0.16, 'x (km)', 0.016,0.,0)
      CALL plchhq(.14,0.50, 'z (km)', 0.016,90.,0)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.20,0.16, '0.0', 0.016,0.,0)
      CALL plchhq(.80,0.16, '1.5', 0.016,0.,0)
      CALL plchhq(.15,0.20, '0.0', 0.016,0.,0)
      CALL plchhq(.15,0.80, '1.5', 0.016,0.,0)
      call frame
c   w
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call perim (3,5,3,5)
      call cpsetc('ILT',' ')
      call cpcnrc(uzpl,nx,nx,nz,.2,5.,.2,-1,-1,1)
      call cpcnrc(uzpl,nx,nx,nz,-5.,-.2,.2,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.50,0.85,'vertical flow', 0.016,0.,0)
      CALL plchhq(.50,0.16, 'x (km)', 0.016,0.,0)
      CALL plchhq(.14,0.50, 'z (km)', 0.016,90.,0)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.20,0.16, '0.0', 0.016,0.,0)
      CALL plchhq(.80,0.16, '1.5', 0.016,0.,0)
      CALL plchhq(.15,0.20, '0.0', 0.016,0.,0)
      CALL plchhq(.15,0.80, '1.5', 0.016,0.,0)
      call frame
C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    


c MAIN TIME LOOP STARTS HERE:
cc

                do it=1,ntstp
           time=it*dt

cccccc add relaxation to initial profiles in the lower part of the domain
        relax_depth=200.      ! e-folding relaxation depth
        relax_time=5. * 60.   ! relaxation time scale
        do k=1,nz
        zz=float(k-1)*dz
        tau=relax_time * exp(zz/relax_depth)
        do i=1,nx
        theta(i,k)=theta(i,k) - (theta(i,k)-theta_e(k))*dt/tau
           qv(i,k)=   qv(i,k) - (   qv(i,k)-   qv_e(k))*dt/tau
        enddo
        enddo

c
c -------> call 2D advection for all fields that move
c          with the air
       call mpdata2d(uxa,uza ,theta,gac,nx,nz,iord,isor,nonos,idiv,1)
       call mpdata2d(uxa,uza ,  qv ,gac,nx,nz,iord,isor,nonos,idiv,2)
       call mpdata2d(uxa,uza ,  qc ,gac,nx,nz,iord,isor,nonos,idiv,3)
c
c -------> modify vertical Courant numbers for precipitating field(s):
c          (only rain in this case):
         do k=2,nz
          rho_shift=.5*(rho(k)+rho(k-1))
         do i=1,nx
          qr_shift=amax1(0.,.5*(qr(i,k)+qr(i,k-1)))
          vt=vterm(qr_shift,rho_shift)
           uzam(i,k)=uza(i,k)-rho_shift*vt*dt/dz
         enddo
         enddo
c extrapolate k=1:
          rho_shift=1.5*rho(1)-.5*rho(2)
         do i=1,nx
          qr_shift=amax1(0.,1.5*qr(i,1)-.5*qr(i,2))
          vt=vterm(qr_shift,rho_shift)
           uzam(i,1)=uza(i,1)-rho_shift*vt*dt/dz
         enddo
c
       call mpdata2d(uxa,uzam,  qr ,gac,nx,nz,iord,isor,nonos,idiv,4)
c
c -------> call microphysical adjustement routine:
       call micro(theta,qv,qc,qr,nx,nz,dt)
cc
      if(amod(time,300.).eq.0.) then
cc calculate and output diagnostics:
      qcmx=-10.
      qcmn= 10.
      qrmx=-10.
      qrmn= 10.
      iqc=1
      kqc=1
      iqr=1
      kqr=1
      do i=1,nx
      do k=1,nz
         qcmx=amax1(qcmx,qc(i,k))
         qcmn=amin1(qcmn,qc(i,k))
         qrmx=amax1(qrmx,qr(i,k))
         qrmn=amin1(qrmn,qr(i,k))
           if(qcmx.eq.qc(i,k)) then
             iqc=i
             kqc=k
           endif
           if(qrmx.eq.qr(i,k)) then
             iqr=i
             kqr=k
           endif
      enddo
      enddo
         print*,time,' ******** TIME *******'
         print*,qcmx,qcmn,'      **** qcmx,qcmn: '
         print*,iqc,kqc,'             i,k position of the qcmx'
         print*,qrmx,qrmn,'      **** qrmx,qrmn: '
         print*,iqr,kqr,'             i,k position of the qrmx'
         vt=vterm(qr(iqr,kqr),rho(kqr))
         prec=qr(iqr,kqr)*rho(kqr)*vt
         prec=prec*3600.    ! convert from kg/m**2/s to mm/hr
         print*,prec,'                corresponding precip (mm/hr)'
cc k=1 precipitation:
         k=1
         do i=1,nx
         vt=vterm(qr(i,k),rho(k))
         prec=qr(i,k)*rho(k)*vt
         prec=prec*3600.    ! convert from kg/m**2/s to mm/hr
         sprec(i)=prec
         enddo
         sum=0.
         imx=1
         amx=-1000.
         do i=1,nx
         sum=sum+sprec(i)
         amx=amax1(amx,sprec(i))
           if(amx.eq.sprec(i)) imx=i
         enddo
         sum=sum/float(nx)
          print*,amx,imx,sum,'   SURF PREC: max,imax,aver'
      endif

c -------> write and plot thermodynamic fields:
      print*,' *** done with minute: ',time/60.

Ccc write data every 15 minutes...
C      if(amod(time,900.).eq.0.) then
C      write(22) theta
C      write(22) qv
C      write(22) qc
C      write(22) qr
C        do i=1,nx
C        do k=1,nz
C        uz(i,k)=.5*(uza(i,k)+uza(i,k+1))
C        enddo
C        enddo
C      write(22) uz
C      endif
cc
      if(amod(time,300.).eq.0.) then
C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    
c   theta
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call periml(3,5,3,5)
      cmx=300.
      cmn=260.
      cnt=.1
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(theta,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc negative values
      ct=-cnt
      call cpcnrc(theta,nx,nx,nz,cmn,ct,cnt,-1,-1,-682)
      write(lhead,350) time
  350 format(' Theta (K)         ',f6.0) 
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchmq(.50,0.85, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.50,0.12, 'X (km)', 0.016,0.,0)
      CALL plchhq(.12,0.50, 'Z (km)', 0.016,90.,0)
      call frame

c   qv
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call periml(3,5,3,5)
      cmx=20.*1.e-3
      cmn=2.*1.e-3
      cnt=0.25*1.e-3
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(qv,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc negative values
      ct=-cnt
      call cpcnrc(qv,nx,nx,nz,cmn,ct,cnt,-1,-1,-682)
      write(lhead,351) time
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchmq(.50,0.85, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.50,0.12, 'X (km)', 0.016,0.,0)
      CALL plchhq(.12,0.50, 'Z (km)', 0.016,90.,0)
      call frame

c   qc
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call periml(3,5,3,5)
      cmx=10.*1.e-3
      cmn=.1 *1.e-3
      cnt=cmn
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(qc,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc trace:
      call cpcnrc(qc,nx,nx,nz,1.e-5,1.1e-5,1.e-5,-1,-1,682)
      write(lhead,352) time
  352 format(' QC (kg/kg)         ',f6.0)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchmq(.50,0.85, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.50,0.12, 'X (km)', 0.016,0.,0)
      CALL plchhq(.12,0.50, 'Z (km)', 0.016,90.,0)
      call frame

c   qr
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call periml(3,5,3,5)
      cmx=10.*1.e-3
      cmn=.002*1.e-3
      cnt=cmn
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(qr,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc trace:
      call cpcnrc(qr,nx,nx,nz,1.e-6,1.1e-6,1.e-6,-1,-1,682)
      write(lhead,353) time
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchmq(.50,0.85, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.50,0.12, 'X (km)', 0.016,0.,0)
      CALL plchhq(.12,0.50, 'Z (km)', 0.016,90.,0)
      call frame
C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    
       endif

c MAIN TIME LOOP ENDS HERE:
            enddo

ccc final solution
c   theta
      call setusv('LW',3000)
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call perim (3,5,3,5)
      cmx=300.
      cmn=260.
      cnt=.2
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(theta,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc negative values
      ct=-cnt
      call cpcnrc(theta,nx,nx,nz,cmn,ct,cnt,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.50,0.85, 'potential temperature (K)', 0.016,0.,0)
      CALL plchhq(.50,0.16, 'x (km)', 0.016,0.,0)
      CALL plchhq(.14,0.50, 'z (km)', 0.016,90.,0)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.20,0.16, '0.0', 0.016,0.,0)
      CALL plchhq(.80,0.16, '1.5', 0.016,0.,0)
      CALL plchhq(.15,0.20, '0.0', 0.016,0.,0)
      CALL plchhq(.15,0.80, '1.5', 0.016,0.,0)
      call frame

c   qv
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call periml(3,5,3,5)
      cmx=20.*1.e-3
      cmn=2.*1.e-3
      cnt=0.25*1.e-3
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(qv,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc negative values
      ct=-cnt
      call cpcnrc(qv,nx,nx,nz,cmn,ct,cnt,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.50,0.85, 'water vapor (g/kg)', 0.016,0.,0)
      CALL plchhq(.50,0.16, 'x (km)', 0.016,0.,0)
      CALL plchhq(.14,0.50, 'z (km)', 0.016,90.,0)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.20,0.16, '0.0', 0.016,0.,0)
      CALL plchhq(.80,0.16, '1.5', 0.016,0.,0)
      CALL plchhq(.15,0.20, '0.0', 0.016,0.,0)
      CALL plchhq(.15,0.80, '1.5', 0.016,0.,0)
      call frame

c   qc
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call perim (3,5,3,5)
      cmx=10.*1.e-3
      cmn=.1 *1.e-3
      cnt=cmn
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(qc,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc trace:
      call cpcnrc(qc,nx,nx,nz,1.e-5,1.1e-5,1.e-5,-1,-1,682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.50,0.85, 'cloud water (g/kg)', 0.016,0.,0)
      CALL plchhq(.50,0.16, 'x (km)', 0.016,0.,0)
      CALL plchhq(.14,0.50, 'z (km)', 0.016,90.,0)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.20,0.16, '0.0', 0.016,0.,0)
      CALL plchhq(.80,0.16, '1.5', 0.016,0.,0)
      CALL plchhq(.15,0.20, '0.0', 0.016,0.,0)
      CALL plchhq(.15,0.80, '1.5', 0.016,0.,0)
      call frame

c   qr
      call set(.2,.8,.2,.8, 0.,xl,0.,zl,1)
      call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
      call perim (3,5,3,5)
      cmx=10.*1.e-3
      cmn=.002*1.e-3
      cnt=cmn
cc positive values
      ct=cnt
      call cpsetc('ILT',' ')
      call cpcnrc(qr,nx,nx,nz,ct,cmx,cnt,-1,-1,1)
cc trace:
      call cpcnrc(qr,nx,nx,nz,1.e-6,1.1e-6,1.e-6,-1,-1,682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.50,0.85, 'drizzle water (g/kg)', 0.016,0.,0)
      CALL plchhq(.50,0.16, 'x (km)', 0.016,0.,0)
      CALL plchhq(.14,0.50, 'z (km)', 0.016,90.,0)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.20,0.16, '0.0', 0.016,0.,0)
      CALL plchhq(.80,0.16, '1.5', 0.016,0.,0)
      CALL plchhq(.15,0.20, '0.0', 0.016,0.,0)
      CALL plchhq(.15,0.80, '1.5', 0.016,0.,0)
      call frame

C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    
          call clsgks
C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS  C NCAR GRAPHICS    
        stop
        end

      subroutine init(nz1,dz1)
      parameter(nz=75,dz=20.)
      common /environ1/ rho(nz)
      common /environ2/ temp_e(nz),theta_e(nz),pres_e(nz),
     1                  qv_e(nz),qc_e(nz)
      data g /9.72/

cc note that saturated water vapor pressure and latent heat
cc of vaporozation are given for the reference temperature
cc of 20/10 deg C.
      data rg,cp,rv,t00,e00,hlat 
cc 20 deg C
c     1       /287.,1005.,461.,293.16,2337.,2.45e6/ 
cc 10 deg C
     1       /287.,1005.,461.,283.16,1228.,2.45e6/ 

      PARAMETER(NPIN=nz)
      DIMENSION press(npin),theta_l(npin),q_tot(npin),zin(npin)
      DIMENSION temp(npin),vap(npin),qcin(npin)

cccc  check consistency:
      if(nz.ne.nz1.or.dz.ne.dz1) then
      print*,' *** inconsistent setup in init, stop'
      stop 'init'
      endif
cc

cc initial theta_l and q_tot profiles:
         do k=1,npin
         zin(k)=float(k-1)*dz
           theta_l(k)=289.
           q_tot(k)=7.5e-3
         enddo

cc  we assume that the first level is undersaturated
cc  (no cloud), temp is potential temp here
        press(1)=1015.  ! surface pressure in hPa
        vap(1)=q_tot(1)*1.e3  ! in g/kg
        qcin(1)=0.
        temp(1)=theta_l(1)

cc decompose theta_l and q_tot to get cloud water at higher levels:
cc predictor-corrector approach to get pressure

      a=rg/rv
      b=hlat/(rv*t00)
      c=hlat/cp
      d=hlat/rv
      e=-cp/rg
      f=rg/cp
            do k=2,nz

cc inital guess for pressure:
            km=k-1
            tempkm=temp(km)*(1.e3/press(km))**(-rg/cp) 
     1                            * (1.+.6e-3*vap(km))
             rhokm=press(km)*1.e2/(rg*tempkm)
             press(k)=press(km)*1.e2 - g*dz*rhokm 
             press(k)=press(k)/1.e2

              temp(k)=theta_l(k)
              vap(k)=q_tot(k)*1.e3  ! in g/kg

cc iterate pressure and theta_l/q_tot decomposition
       qcin(k)=0.
               do iter=1,5
      pre=press(k)*1.e2
      thetme=(1.e5/pre)**f
      thi=1./temp(k)
      y=b*thetme*t00*thi
      ees=e00*exp(b-y)
      qvs=a*ees/(pre-ees)
ccc linearized condensation rate is next:
      cf1=thetme*thetme*thi*thi
      cf1=c*cf1*pre/(pre-ees)*d
      delta=(vap(k)*1.e-3-qvs)/(1.+qvs*cf1)
c--->
      delta=amax1(0.,delta) 
         vap(k)=vap(k)-delta*1.e3
         qcin(k)=qcin(k)+delta
         temp(k)=temp(k)+c*thetme*delta

             pressa=(press(km)+press(k))/2.
            tempk=temp(k)*(1.e3/press(k))**(-rg/cp) 
     1                            * (1.+.6e-3*vap(k))
             tempa=(tempkm+tempk)/2.
             rhoa=pressa*1.e2/(rg*tempa)
             press(k)=press(km)*1.e2 - g*dz*rhoa 
             press(k)=press(k)/1.e2
 
                enddo ! iter

       print*,'z,t,vap,qc: ',(k-0.5)*dz,temp(k),vap(k),qcin(k)

            enddo  ! levels


compute environmental profiles from input sounding:
      do 64 k=1,nz
          theta_e(k)=temp(k)
          qv_e(k)=vap(k)*1.e-3
          qc_e(k)=qcin(k)
          presnl=press(k)
          pres_e(k)=presnl*1.e2
          temp_e(k)=theta_e(k) * (1000./presnl)**(-rg/cp)
          te_virt=temp_e(k)*(1.+.6*qv_e(k))
          rho(k)=pres_e(k)/(rg*te_virt)
cccccccc
          print 284,k,theta_e(k),temp_e(k),qv_e(k)*1.e3,qc_e(k)*1.e3,
     1            pres_e(k)*1.e-2,rho(k)
284   format(1x,'k,ts,qs,pre,rho: ',i4,6f10.3)
cccccccc
 64   continue
      return
      end

       subroutine adv_vel(ux,uz,mx,mz,dx,dz,dt)
c  this routine provides flow input to the model
c       ON INPUT:              
c                      mx     - x dimension of scalar arrays
c                      mz     - z dimension of scalar arrays
c                      dx,dz  - grid intervals in x and z direction
c       ON OUTPUT:
c                     ux,uz - rho times velocity components normalized
c                             by dx/dt or dz/dt
c
C
      parameter(nx=75,nz=75)
      parameter(nxp=nx+1,nzp=nz+1)
cc streamfunction and positions of its gridpoints:
      dimension phi(nxp,nzp),xp(nxp),zp(nzp)
cc arrays with velocities 
      dimension ux(mx+1,mz),uz(mx,mz+1)
C
C CHECK DIMENSIONS:
      if(mx.ne.nx.or.mz.ne.nz) then
      print*,' *** dimensions do not match in adv_vel, stop'
      stop
      endif
C
C CALCULATE X AND Z DISTANCES (IN METERS)
      do k=1,nzp
      zp(k)=(k-1)*dz
      enddo
      do i=1,nxp
      xp(i)=(i-1)*dx
      enddo
      XSCALE=xp(nxp)
      ZSCALE=zp(nzp)
      pi=4.*atan(1.)
C
CC INITIAL DATA FOR THE STREAMFUNCTION. AMPL IS WMAX IN M/S,
C
c      AMPL=1.0
      AMPL=0.6
C
C ADOPT THE AMPL FOR STREAMFUNCTION CALCULATION TO BE W_MAX/K_x
      AMPL=AMPL/pi*XSCALE
C
C DEFINE STREAMFUNCTION AS A FUNCTION OF HEIGHT
C ALSO, CENTRALIZE THE UPDRAFT TO OCCUPY ONLY THE INNER XSCALE
C OF THE DOMAIN
C
      ZTOP=zp(nzp)/ZSCALE
      XCEN=.5*xp(nxp)
      X0=(xp(nxp)-XSCALE)/2.
      DO 1 I=1,NXP
      DO 1 K=1,NZP
      PHI(i,k)=-cos(2.*pi*(xp(i)-X0)/XSCALE)*sin(pi*zp(k)/ZSCALE)
      PHI(i,k)=PHI(i,k)*AMPL
 1    CONTINUE
C

cc calculate rho*vel by derivation of streamfunction and normalize
cc rho*ux velocity:
      do i=1,nxp
      do k=1,nz
      ux(i,k)=-(phi(i,k+1)-phi(i,k))/dz  *dt/dx
      enddo
      enddo
cc rho*uz velocity
      do k=1,nzp
      do i=1,nx
      uz(i,k)=(phi(i+1,k)-phi(i,k))/dx  *dt/dz
      enddo
      enddo
ccc
      return
      end

      SUBROUTINE MPDATA2D(U1,U2,X,H,N,M,IORD,ISOR,NONOS,IDIV,IFL)
C THIS SUBROUTINE SOLVES 2-D ADVECTIVE TRANSPORT IN CARTESIAN GEOMETRY
C ON STAGGERRED GRID (X,Y VELOCITIES SHIFTED HALF GRID IN X, Y DIR, RESP)
C*************************************************************************
C ADVECTION ALGORITHM: IORD - NUMBER OF ITERATIONS (IORD=1 OVERWRITES
C CALLS TO OTHER OPTIONS AND GIVES SIMPLE UPSTREAM SCHEME); ISOR=1 2ND ORDER
C COMPUTATIONS WHEREAS ISOR=3 AND IORD=3 3D ORDER SCHEME; IDIV=1 ACTIVATES
C CORRECTION FOR DIVERGENT FLOW; NONOS=1 STRICTLY MONOTONE ADVECTION
C   N O T E: idiv MUST be 0 for a nondivergent flow
C A GOOD POINT TO START WOULD BE:
C      PARAMETER(IORD0=2,ISOR0=1,IDIV0=0,NONO=0)
C IFL IS THE FLAG TO DISTINGUISH BETWEEN FIELDS THAT ARE BEING ADVECTED;
C THAT IS NECESSARY TO PROVIDE LATERAL BOUNDARY CONDITIONS ON INLFOW
C AND LOWER BC FOR SEDIMENTING FIELDS
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C ** NOTE THAT THIS ROUTINE WILL WORK FOR FIELD WITH VARIABLE SIGN **
C ** (AS MOMENTUM) SINCE ABSOLUTE VALUES ARE USE IN THE DEFINITION **
C **             OF ANTIDIFFUSIVE VELOCITIES                       **
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C AUTHOR: Piotr Smolarkiewicz (smolar@ncar.ucar.edu), (303)-497-8972
C  modified for this test by Wojciech Grabowski (grabow@ncar.ucr.edu)
C*************************************************************************
      PARAMETER(N1=76,N2=76)
      PARAMETER(N1M=N1-1,N2M=N2-1)
      DIMENSION U1(N+1,M),U2(N,M+1),X(N,M),H(N,M)
      DIMENSION V1(N1,N2M),V2(N1M,N2),F1(N1,N2M),F2(N1M,N2)
     *         ,CP(N1M,N2M),CN(N1M,N2M)
      REAL MX(N1M,N2M),MN(N1M,N2M)
      DATA EP/1.E-10/
C
ccc for use NOT on a CRAY computer you have to replace CVMGM
ccc with a substitute which gives the same result: below is 
ccc an example of something that will work:
ccc CVMGM(a,b,c)= a if c.lt.0 or b if c.ge.0
      aneg(d)=.5*(1.-sign(1.,d))
      apos(d)=.5*(1.+sign(1.,d))
      cvmgm(a,b,c)=aneg(c)*a + apos(c)*b
ccc
      DONOR(Y1,Y2,A)=CVMGM(Y2,Y1,A)*A
      VDYF(X1,X2,A,R)=(ABS(A)-A**2/R)*(ABS(X2)-ABS(X1))
     1                               /(ABS(X2)+ABS(X1)+EP)
      VCORR(A,B,Y1,Y2,R)=-0.125*A*B*Y1/(Y2*R)
      VCOR31(A,X0,X1,X2,X3,R)= -(A -3.*ABS(A)*A/R+2.*A**3/R**2)/3.
     1                         *(ABS(X0)+ABS(X3)-ABS(X1)-ABS(X2))
     2                         /(ABS(X0)+ABS(X3)+ABS(X1)+ABS(X2)+EP)
      VCOR32(A,B,Y1,Y2,R)=0.25*B/R*(ABS(A)-2.*A**2/R)*Y1/Y2
      VDIV1(A1,A2,A3,R)=0.25*A2*(A3-A1)/R
      VDIV2(A,B1,B2,B3,B4,R)=0.25*A*(B1+B2-B3-B4)/R
      PP(Y)= AMAX1(0.,Y)
      PN(Y)=-AMIN1(0.,Y)
C
cc    test dimensions:
      if(n1m.ne.n.or.n2m.ne.m) then
      print*,' dimensions do not match in advection. stop'
      stop 'mpdata'
      endif
      IF(ISOR.EQ.3) IORD=MAX0(IORD,3)
C
      DO 1 J=1,N2-1
      DO 1 I=1,N1
    1 V1(I,J)=U1(I,J)
      DO 2 J=1,N2
      DO 2 I=1,N1-1
    2 V2(I,J)=U2(I,J)
C
      IF(NONOS.EQ.1) THEN
      DO 400 J=2,N2-2
      DO 400 I=2,N1-2
      MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
  400 MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
      ENDIF
C
                         DO 3 K=1,IORD
C
      
      DO 337 J=1,N2-1
      DO 331 I=2,N1-1
  331 F1(I,J)=DONOR(X(I-1,J),X(I,J),V1(I,J))
      if(k.eq.1) then
      F1( 1,J)=DONOR(X(N1-1,J),X(1,J),V1(1,J))
      F1(N1,J)=DONOR(X(N1-1,J),X(1,J),V1(N1,J))
      else
      F1( 1,J)=0.
      F1(N1,J)=0.
      endif
  337 continue

      DO 338 I=1,N1-1
      DO 332 J=2,N2-1
  332 F2(I,J)=DONOR(X(I,J-1),X(I,J),V2(I,J))
      if(k.eq.1.and.ifl.eq.4) then
      F2(I, 1)=DONOR(0.,X(I,1),V2(I,1))  ! fallout of rain
      F2(I,N2)=0.                        ! no fall-in of rain
      else
      F2(I, 1)=0.
      F2(I,N2)=0.
      endif
  338 continue
C 
      DO 333 J=1,N2-1
      DO 333 I=1,N1-1
  333 X(I,J)=X(I,J)-(F1(I+1,J)-F1(I,J)+F2(I,J+1)-F2(I,J))/H(I,J)
C
      IF(K.EQ.IORD) GO TO 6
      DO 49 J=1,N2-1
      DO 49 I=1,N1
      F1(I,J)=V1(I,J)
   49 V1(I,J)=0.
      DO 50 J=1,N2
      DO 50 I=1,N1-1
      F2(I,J)=V2(I,J)
   50 V2(I,J)=0.
      DO 51 J=2,N2-2
      DO 51 I=2,N1-1
   51 V1(I,J)=VDYF(X(I-1,J),X(I,J),V1(I,J),.5*(H(I-1,J)+H(I,J)))
     *       +VCORR(V1(I,J), F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),
     *   ABS(X(I-1,J+1))+ABS(X(I,J+1))-ABS(X(I-1,J-1))-ABS(X(I,J-1)),
     *   ABS(X(I-1,J+1))+ABS(X(I,J+1))+ABS(X(I-1,J-1))+ABS(X(I,J-1))+EP,
     *                 .5*(H(I-1,J)+H(I,J)))
      IF(IDIV.EQ.1) THEN
      DO 511 J=2,N2-2
      DO 511 I=2,N1-1
  511 V1(I,J)=V1(I,J)
     *    -VDIV1(F1(I-1,J),F1(I,J),F1(I+1,J),.5*(H(I-1,J)+H(I,J)))
     *    -VDIV2(F1(I,J),F2(I-1,J+1),F2(I,J+1),F2(I-1,J),F2(I,J),
     *                 .5*(H(I-1,J)+H(I,J)))
      ENDIF
      DO 52 J=2,N2-1
      DO 52 I=2,N1-2
   52 V2(I,J)=VDYF(X(I,J-1),X(I,J),V2(I,J),.5*(H(I,J-1)+H(I,J)))
     *       +VCORR(V2(I,J), F1(I,J-1)+F1(I,J)+F1(I+1,J)+F1(I+1,J-1),
     *   ABS(X(I+1,J-1))+ABS(X(I+1,J))-ABS(X(I-1,J-1))-ABS(X(I-1,J)),
     *   ABS(X(I+1,J-1))+ABS(X(I+1,J))+ABS(X(I-1,J-1))+ABS(X(I-1,J))+EP,
     *                 .5*(H(I,J-1)+H(I,J)))
      IF(IDIV.EQ.1) THEN
      DO 521 J=2,N2-1
      DO 521 I=2,N1-2
  521 V2(I,J)=V2(I,J)
     *    -VDIV1(F2(I,J-1),F2(I,J),F2(I,J+1),.5*(H(I,J-1)+H(I,J)))
     *    -VDIV2(F2(I,J),F1(I+1,J),F1(I+1,J-1),F1(I,J-1),F1(I,J),
     *                 .5*(H(I,J-1)+H(I,J)))
      ENDIF
      IF(ISOR.EQ.3) THEN
      DO 61 J=2,N2-2
      DO 61 I=3,N1-2
   61 V1(I,J)=V1(I,J)     +VCOR31(F1(I,J),
     1        X(I-2,J),X(I-1,J),X(I,J),X(I+1,J),.5*(H(I-1,J)+H(I,J)))
      DO 62 J=2,N2-2
      DO 62 I=3,N1-2
   62 V1(I,J)=V1(I,J)
     1 +VCOR32(F1(I,J),F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),
     *   ABS(X(I,J+1))-ABS(X(I,J-1))-ABS(X(I-1,J+1))+ABS(X(I-1,J-1)),
     *   ABS(X(I,J+1))+ABS(X(I,J-1))+ABS(X(I-1,J+1))+ABS(X(I-1,J-1))+EP,
     *                   .5*(H(I-1,J)+H(I,J)))
      DO 63 J=3,N2-2
      DO 63 I=2,N1-2
   63 V2(I,J)=V2(I,J)     +VCOR31(F2(I,J),
     1        X(I,J-2),X(I,J-1),X(I,J),X(I,J+1),.5*(H(I,J-1)+H(I,J)))
      DO 64 J=3,N2-2
      DO 64 I=2,N1-2
   64 V2(I,J)=V2(I,J)
     1 +VCOR32(F2(I,J),F1(I,J-1)+F1(I+1,J-1)+F1(I+1,J)+F1(I,J),
     *   ABS(X(I+1,J))-ABS(X(I-1,J))-ABS(X(I+1,J-1))+ABS(X(I-1,J-1)),
     *   ABS(X(I+1,J))+ABS(X(I-1,J))+ABS(X(I+1,J-1))+ABS(X(I-1,J-1))+EP,
     *                   .5*(H(I,J-1)+H(I,J)))
      ENDIF
C
C
      IF(NONOS.EQ.0) GO TO 3
C                 NON-OSSCILATORY OPTION
      DO 401 J=2,N2-2
      DO 401 I=2,N1-2
      MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MX(I,J))
  401 MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MN(I,J))
C
      DO 402 J=2,N2-2 
      DO 402 I=2,N1-1
  402 F1(I,J)=DONOR(X(I-1,J),X(I,J),V1(I,J))
      DO 403 J=2,N2-1
      DO 403 I=2,N1-2
  403 F2(I,J)=DONOR(X(I,J-1),X(I,J),V2(I,J))
      DO 404 J=2,N2-2   
      DO 404 I=2,N1-2
      CP(I,J)=(MX(I,J)-X(I,J))*H(I,J)/
     1(PN(F1(I+1,J))+PP(F1(I,J))+PN(F2(I,J+1))+PP(F2(I,J))+EP)
      CN(I,J)=(X(I,J)-MN(I,J))*H(I,J)/
     1(PP(F1(I+1,J))+PN(F1(I,J))+PP(F2(I,J+1))+PN(F2(I,J))+EP)
  404 CONTINUE
      DO 405 J=3,N2-2 
      DO 405 I=3,N1-2 
      V1(I,J)=PP(V1(I,J))*
     1  ( AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1., X(I-1,J)))
     1   +AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1.,-X(I-1,J))) )
     2       -PN(V1(I,J))*
     2  ( AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1., X(I ,J )))
     2   +AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1.,-X(I ,J ))) )
  405 V2(I,J)=PP(V2(I,J))*
     1  ( AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1., X(I,J-1)))
     1   +AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1.,-X(I,J-1))) )
     1       -PN(V2(I,J))*
     2  ( AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1., X(I ,J )))
     2   +AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1.,-X(I ,J ))) )
C
    3                      CONTINUE
    6 CONTINUE
      RETURN
      END   

      subroutine micro(th,qv,qc,qr,nx,nz,dt)
      dimension th(nx,nz),qv(nx,nz),qc(nx,nz),qr(nx,nz)
      parameter(mz=75)
c environmental temperature, potential temperature, pressure and
c water vapor mixing ratio
      common /environ2/ temp_e(mz),theta_e(mz),pres_e(mz),
     1                  qv_e(mz),qc_e(mz)
c base state density profile
      common /environ1/ rho(mz)

cc note that saturated water vapor pressure and latent heat
cc of vaporozation are given for the reference temperature
cc of 20 deg C.
      data rg,cp,rv,t00,e00,hlat 
     1       /287.,1005.,461.,293.16,2337.,2.45e6/ 

cc coefficients below are for warm rain formulation:
      data rac,qctr,rc/1.e-3, .4e-3, 2.2/

cc coefficients below are for Berry's parameterization:
c      data dconc,ddisp /100.,.32/
      data dconc,ddisp /1000.,.22/
cc
cc  check consistency:
       if(nz.ne.mz) then
       print*,' *** inconsistent setup in micro, stop'
       stop 'micro'
       endif
cc
condensation/evaporation 
      dti=1./dt
      a=rg/rv
      b=hlat/(rv*t00)
      c=hlat/cp
      d=hlat/rv
      e=-cp/rg
      do 10 k=1,nz
      pre=pres_e(k)
      thetme=theta_e(k)/temp_e(k)
      do 10 i=1,nx
      thi=1./th(i,k)
      y=b*thetme*t00*thi
      ees=e00*exp(b-y)
      qvs=a*ees/(pre-ees)
ccc linearized condensation rate is next:
      cf1=thetme*thetme*thi*thi
      cf1=c*cf1*pre/(pre-ees)*d
      delta=(qv(i,k)-qvs)/(1.+qvs*cf1)
c--->
ccc one Newton-Raphson iteration is next:
      thi=1./(th(i,k)+c*thetme*delta)
      y=b*thetme*t00*thi
      ees=e00*exp(b-y)
      qvs=a*ees/(pre-ees)
      fff=qv(i,k)-delta-qvs
      cf1=thetme*thetme*thi*thi
      cf1=c*cf1*pre/(pre-ees)*d
      fffp=-1.-qvs*cf1
      delta=delta-fff/fffp
ccc end of the iteration; if required, it can be repeated
c--->
      delta=amax1(-qc(i,k),delta) 
      qv(i,k)=qv(i,k)-delta
      qc(i,k)=qc(i,k)+delta
      th(i,k)=th(i,k)+c*thetme*delta
   10 continue

compute forces due to rain processes
      do 20 k=1,nz
      pre=pres_e(k)
      thetme=theta_e(k)/temp_e(k)
      do 20 i=1,nx
      thi=1./th(i,k)
      y=b*thetme*t00*thi
      ees=e00*exp(b-y)
      qvs=a*ees/(pre-ees)
      ss=amin1(qv(i,k)/qvs-1., 0.)
c
c Kessler:
c      dcol= rac*amax1(qc(i,k)-qctr, 0.)
c     .         + rc*qc(i,k)*amax1(0.,qr(i,k))**.875
c Berry:
      del2=1.e3*rho(k)*qc(i,k)
      del1=1./rho(k)*1.67e-5*del2*del2 /
     1 (5. + .0366*dconc/(ddisp*(del2+1.E-8)))
      autocon=del1
      dcol=  autocon + rc*qc(i,k)*amax1(0.,qr(i,k))**.875
c
      dcol=amin1(dcol,dti*qc(i,k))
      presmb=pre/100.
      rhqr=rho(k)*1.e-3*amax1(0.,qr(i,k))
      qvs=qv(i,k)/(1.+ss)
      bottom=1.e-3*rho(k)*(5.4e5 + 2.55e6/(presmb*qvs))
      ventc=1.6+124.9*rhqr**.2046
      devp=ss*ventc*rhqr**.525 / bottom 
      devp=amax1(devp, -dti*amax1(0.,qr(i,k)))
cc
      qr(i,k)=qr(i,k) + (devp+dcol)*dt
      qc(i,k)=qc(i,k)-dcol*dt
      qv(i,k)=qv(i,k)-devp*dt
      th(i,k)=th(i,k)+c*devp*thetme*dt
   20 continue

      return
      end

