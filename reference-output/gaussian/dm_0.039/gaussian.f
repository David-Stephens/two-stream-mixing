C     Run the program
      program two_stream

      implicit none

c     number of species
      integer, parameter :: nsp = 2

c     maximum dimension of up and down stream abundance and cell arrays
C     NOTE: This is modified automatically by two-streams-test.py
      integer, parameter :: nmax = 512

C     The dimension of shell quantities
      integer, parameter :: nshellmax = nmax + 1

c     maximum dimension of stitched stream abundance and cell mass arrays
      integer, parameter :: nmax2 = 2*nmax
      integer, parameter :: nshellmax2 = nmax2 + 1

c     number of ghost cells
      integer, parameter :: ng = 2

c     maximum dimension of stitched stream abundance and cell mass
c     arrays with ghost cells
      integer, parameter :: ngmax = 2*(nmax+ng)

C     maximum dimension of stitched stream abundance slopes with ng-1
c     ghost cells
      integer, parameter :: nslopes = 2*(nmax+ng-1)

C     The number of time steps that this simulation will do
C     NOTE: This is modified automatically by two-streams-test.py
      integer, parameter :: steps = 2048

C     For every x steps, we will save an output file for plots
C     NOTE: This is modified automatically by two-streams-test.py
      integer, parameter :: saveevery = 5

C     The dt for a single step. This is constant in all models
C     NOTE: This is modified automatically by two-streams-test.py
      double precision, parameter :: dt = 0.0195312500d0

c     up and down quantities and central quantities, to be read from a
c     file
      double precision, dimension (nsp,nmax) :: Xup, Xdown
      double precision, dimension (nmax) :: gammaadv, dm

C     shell quantities to be read from a file
      double precision, dimension (nshellmax) :: v_shell, rho_shell
      double precision, dimension (nshellmax) :: r_shell

C     characters for output filename
      character*19 :: filename
C     I add a dash here so I can format the string
      character*8, parameter :: testcase = 'gaussian'

C     useful integers throughout
      integer i, j, k, savecount, saved

C     all central quantities are going to be read-in
C     NOTE: This is modified automatically by two-streams-test.py
      open(UNIT=1,FILE=testcase//".central",STATUS="OLD",
     $     FORM='FORMATTED')
 1    format(E18.0,E18.0,E18.0,E18.0)

      do i=1,nmax

C     Xup, Xdown, gammaadv, dm
         read(1,*) Xup(1,i), Xdown(1,i), gammaadv(i), dm(i)

C     Now we can simply write (1-X) for other species
         Xup(2,i) = 1.0d0 - Xup(1,i)
         Xdown(2,i) = 1.0d0 - Xdown(1,i)
      end do

C     close file
      close(1)

C     all shell quantities are going to be read-in
C     NOTE: This is modified automatically by two-streams-test.py
      open(UNIT=2,FILE=testcase//'.shell',STATUS='OLD',FORM='FORMATTED')
 2    format(E18.0,E18.0,E18.0)

      do i=1,nshellmax
C     v, rho, r
         read(2,*) v_shell(i), rho_shell(i), r_shell(i)
      end do

C     close file
      close(2)

C     initialize these ints
      savecount = 0
      saved = 0

C     Now, we loop through all of the steps we need to take
      do i=1,steps
         
C     Count the number of steps since saveevery
         savecount = savecount + 1

C     If we have savecount == saveevery then save file
         if ( (i .EQ. 1) .OR. (savecount .EQ. saveevery) ) then
            
C     create the empty file which we will write to
            write(filename, '(A9, i6.6, ".txt")') testcase//'-', saved
            open(UNIT=saved+100,FILE=filename,STATUS='UNKNOWN',
     $           FORM='FORMATTED')

C     Now we write Xup, (1-Xup), Xdown, (1-Xdown)
            do j=1,nmax
               write(saved+100,'(E18.12,1x,E18.12,1x,E18.12,1x,E18.12)')
     $              Xup(1,j), Xup(2,j), Xdown(1,j), Xdown(2,j)
            end do

C     Close the file
            close(saved+100)

C     Set savecout to zero; increment saved
            savecount = 0
            saved = saved + 1

         end if

C        ok we are good to go and call advect
         call advect(Xup,Xdown,dm,gammaadv,v_shell,rho_shell,r_shell,
     $ nsp,nmax,nmax2,nshellmax,nshellmax2,ng,ngmax,nslopes,dt)

      end do

C     We didnt have a chance to write the last file because of the
C     the way we looped!
      
C     create the empty file which we will write to
      write(filename, '(A9, i6.6, ".txt")') testcase//'-', saved
      open(UNIT=saved+100,FILE=filename,STATUS='UNKNOWN',
     $     FORM='FORMATTED')

C     Now we write Xup, (1-Xup), Xdown, (1-Xdown)
      do j=1,nmax
         write(saved+100,'(E18.12,1x,E18.12,1x,E18.12,1x,E18.12)')
     $        Xup(1,j), Xup(2,j), Xdown(1,j), Xdown(2,j)
      end do

C     Close the file
      close(saved+100)

C     Stop the program!
      stop

      end program two_stream




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC





C     This code was written by David Stephens at the University of 
C     Victoria in 2019. The details of the algorithm is described
C     within the paper:
C     https://ui.adsabs.harvard.edu/abs/2020arXiv200110969S/abstract


      subroutine advect(Xup,Xdown,dm,gammaadv,v_shell,rho_shell,r_shell,
     $     nspadv,nmax,nmax2,nshellmax,nshellmax2,ng,ngmax,nslopes,dt)

C     this subroutine advects the up- and down-streams of species for a
C     single timestep which will modify the Xup and Xdown arrays.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE INPUTS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  Xup, Xdown      double SHAPE(nspadv,nmax). Contains the mass fractions 
C                  of species k along all mass zones, i
C
C  dm              double SHAPE(nmax). A cells mass (this is the
C                  difference between mass coordinates of two shells
C                  divided by two)
C
C  gammaadv        double SHAPE(nmax). The additional horizontal mixing
C                  as defined in equation 17
C
C  v_shell,        double SHAPE(nshellmax). The velocity, density and
C  rho_shell,      radius at a shell. These are the same for the up- and
C  r_shell         down-streams
C
C  nspadv          int. The number of species, k
C
C  nmax            int. The number of mass zones, i
C
C  nmax2           int. The number of mass zones of the stitched arrays
C
C  nshellmax       int. The number of shell interfaces
C
C  nshellmax2      int. The number of stitched shell interfaces
C
C  ng              int. The number of ghost cells to be added to each
C                  end. MUST BE 2 unless you are using a limiter that is
C                  different from the original code
C
C  ngmax           int. The number of cells that are stitched including
C                  the ghost cells
C
C  nslopes         int. The number of slope estimates including ghost
C                  cells
C
C  dt              double. The explicit timestep for runge-kutta3

      implicit none

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INPUT DOUBLE AND INT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      integer nspadv, nmax, nmax2, nshellmax, nshellmax2
      integer ng, ngmax, nslopes
      double precision dt

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INPUT DOUBLE ARRAYS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     up and down stream mass fractions
      double precision, dimension (nspadv,nmax) :: Xup, Xdown

C     cell quantities
      double precision, dimension (nmax) :: gammaadv, dm

C     shell quantities
      double precision, dimension (nshellmax) :: v_shell, rho_shell
      double precision, dimension (nshellmax) :: r_shell

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     CREATE STITCHED ARRAYS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     stitched stream abundance arrays
      double precision, dimension (nspadv,nmax2) :: mX_stitch, X_stitch

C     Temporary rk interpolated steps
      double precision, dimension (nspadv,nmax2) :: rk1, rk2, rk3

c     stitched cell quantities
      double precision, dimension (nmax2) :: beta_stitch
      double precision, dimension (nmax2) :: dm_stitch, gammaadv_stitch

C     stitched shell quantities which are used for fluxes
      double precision, dimension (nshellmax2) :: v_stitch, rho_stitch
      double precision, dimension (nshellmax2) :: r_stitch

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE VARIABLES
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     rk factors
      double precision, parameter :: rkf1 = 2.0d0/9.0d0
      double precision, parameter :: rkf2 = 3.0d0/9.0d0
      double precision, parameter :: rkf3 = 4.0d0/9.0d0

C     Looping integers
      integer i, j, k

C     external subroutines called
      EXTERNAL stitch
      EXTERNAL fluxes

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     BEGIN ADVECTION
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     create the stitched arrays from input arrays
      call stitch(Xup,Xdown,dm,gammaadv,v_shell,rho_shell,r_shell,
     $     X_stitch,dm_stitch,beta_stitch,gammaadv_stitch,v_stitch,
     $     rho_stitch,r_stitch,nspadv,nmax,nmax2,nshellmax,nshellmax2)

C     Now I have the X_stitch, I will create my initial mass of every
C     species in each cell
      do k=1,nspadv
         mX_stitch(k,:) = X_stitch(k,:) * dm_stitch
      end do

C     Now I have the central and shell stitched arrays. These are used
C     in flux calls to get the dmXdt slopes at spatial interpolations
C     which are stored within the runge-kutta arrays to be put together
C     at the end for the total explicit timestep

C     Our first runge-kutta interpolated state (we already have it,
C     it is X_stitch!). Since we don't have explicit time dependence,
C     ignore the time interpolation

C     call fluxes to determine the dmXdt (rk1) of the interpolated state,
C     X_stitch
      call fluxes(rk1,X_stitch,dm_stitch,beta_stitch,gammaadv_stitch,
     $     v_stitch,rho_stitch,r_stitch,nspadv,nmax,nmax2,nshellmax,
     $     nshellmax2,nslopes,ng,ngmax)

C     Our second runge-kutta interpolated state
      do k=1,nspadv
         X_stitch(k,:) = (mX_stitch(k,:) + (1.0d0/2.0d0)*dt*rk1(k,:)) /
     $        dm_stitch
      end do

C     call fluxes to determine the dmXdt (rk2) of the interpolated state,
C     X_stitch
      call fluxes(rk2,X_stitch,dm_stitch,beta_stitch,gammaadv_stitch,
     $     v_stitch,rho_stitch,r_stitch,nspadv,nmax,nmax2,nshellmax,
     $     nshellmax2,nslopes,ng,ngmax)

C     Our last runge-kutta interpolated step
      do k=1,nspadv
         X_stitch(k,:) = (mX_stitch(k,:) + (3.0d0/4.0d0)*dt*rk2(k,:)) /
     $        dm_stitch
      end do

C     call fluxes to determine the dmXdt (rk3) of the interpolated state,
C     X_stitch
      call fluxes(rk3,X_stitch,dm_stitch,beta_stitch,gammaadv_stitch,
     $     v_stitch,rho_stitch,r_stitch,nspadv,nmax,nmax2,nshellmax,
     $     nshellmax2,nslopes,ng,ngmax)

C     stitch all of the rk interpolated states into our new final state
      do k=1,nspadv
         mX_stitch(k,:) = mX_stitch(k,:) + rkf1*dt*rk1(k,:) +
     $        rkf2*dt*rk2(k,:) + rkf3*dt*rk3(k,:)
      end do

C     Now convert this back to mass fractions
      do k=1,nspadv
         do i=1,nmax
            Xup(k,i) = mX_stitch(k,i) / dm_stitch(i)
         end do
C     We go backwards for the downstream
         do i=1,nmax
            Xdown(k,nmax+1-i) = mX_stitch(k,nmax+i) / dm_stitch(nmax+i)
         end do
      end do

      return

C     That's all folks!
      end subroutine advect




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINES
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     The external subroutines that are defined here are used to do all
C     of the hard work.


      subroutine stitch(Xup,Xdown,dm,gammaadv,v_shell,rho_shell,r_shell,
     $     X_stitch,dm_stitch,beta_stitch,gammaadv_stitch,v_stitch,
     $     rho_stitch,r_stitch,nspadv,nmax,nmax2,nshellmax,nshellmax2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE INPUTS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Anything that is NOT defined within "advect" is shown here. If you
C     are unsure of any variable definitions in this subroutine look at
C     subroutine advect

C  X_stitch        double SHAPE(nspadv,nmax2). Contains the mass 
C                  fractions of species k along the stitched mass zones,
C                  i
C
C  dm_stitch,      double SHAPE(nmax2). Stitched arrays containing the
C  beta_stitch,    values of dm, beta, and gamma for every cell
C  gammaadv_stitch
C
C  v_stitch,       double SHAPE(nshellmax2). Stitched arrays containing
C  rho_stitch,     the values of velocity, density and radius at a shell
C  r_stitch
C

      implicit none

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INPUT INT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     The integers read in for array sizes
      integer nspadv, nmax, nmax2, nshellmax, nshellmax2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INPUT DOUBLE ARRAYS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     up and down mass fractions and central quantities
      double precision, dimension (nspadv,nmax) :: Xup, Xdown
      double precision, dimension (nmax) :: gammaadv, dm

C     shell quantities
      double precision, dimension (nshellmax) :: v_shell, rho_shell
      double precision, dimension (nshellmax) :: r_shell

c     stitched stream mass fraction array
      double precision, dimension (nspadv,nmax2) :: X_stitch

c     stitched central quantities (no species)
      double precision, dimension (nmax2) :: beta_stitch
      double precision, dimension (nmax2) :: dm_stitch, gammaadv_stitch

C     stitched shell quantities which are used for fluxes
      double precision, dimension (nshellmax2) :: v_stitch, rho_stitch
      double precision, dimension (nshellmax2) :: r_stitch

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE VARIABLES
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Looping integers
      integer i, j, k

C     pi
      double precision :: pi
      pi = 4.0d0 * ATAN(1.0d0)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     CREATE THE STITCHED ARRAYS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     ASSUMING ONLY WORKING WITH THE CONVECTION ZONE

C     Create the CENTRAL stitched arrays. Start with X
C     stitch flipped Xdown atop Xup into X_stitch
      do k = 1,nspadv
         do i = 1,nmax
            X_stitch(k,i) = Xup(k,i)
         end do
         do i = 1,nmax
C     We go backwards for Xdown
            X_stitch(k,nmax+i) = Xdown(k,nmax+1-i)
         end do
      end do

C     I need to calculate b from the quantities I know
      do i = 1,nmax
         beta_stitch(i) = 2.*pi*(rho_shell(i+1)*v_shell(i+1)*
     $        r_shell(i+1)**2.d0 - rho_shell(i)*v_shell(i)*
     $        r_shell(i)**2.d0)
      end do
      do i = 1,nmax
C     We go backwards for bdown, it is the same as bup
         beta_stitch(nmax+i) = beta_stitch(nmax+1-i)
      end do

C     similarly, stitch dm to follow that pattern
      do i = 1,nmax
         dm_stitch(i) = dm(i)
      end do
      do i = 1,nmax
C     We go backwards on dm now
         dm_stitch(nmax+i) = dm(nmax+1-i)
      end do

C     similarly, stitch gammaadv to follow that pattern
      do i = 1,nmax
         gammaadv_stitch(i) = gammaadv(i)
      end do
      do i = 1,nmax
C     We go backwards on gammaadv now
         gammaadv_stitch(nmax+i) = gammaadv(nmax+1-i)
      end do


C     Create the SHELL stitched arrays. In this case I will NOT
C     repeat the top of the convection zone!

C     Loop through shell quantities and put them into stitch
      do i=1,nshellmax
         v_stitch(i) = v_shell(i)
         rho_stitch(i) = rho_shell(i)
         r_stitch(i) = r_shell(i)
      end do

      do i=1,nmax
C     We go backwards in shell but skip the top of the stream
C     (don't duplicate) therefore j goes from nmax to 1
         j = nshellmax-i
         v_stitch(nshellmax+i) = v_shell(j)
         rho_stitch(nshellmax+i) = rho_shell(j)
         r_stitch(nshellmax+i) = r_shell(j)
      end do

      return

      end subroutine stitch




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




      subroutine setBC(Xng,dmng,X_stitch,dm_stitch,beta_stitch,
     $     gammaadv_stitch,v_stitch,rho_stitch,r_stitch,nspadv,nmax,
     $     nmax2,ng,ngmax,nshellmax,nshellmax2)

C     This subroutine sets up the stitched arrays such that the boundary
C     conditions can be easily set. This includes constructing the
C     "ghost" arrays so that the calculations of slopes are vectorized

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE INPUTS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Anything that is NOT defined within "advect" is shown here. If you
C     are unsure of any variable definitions in this subroutine look at
C     subroutine advect

C  Xng             double SHAPE(nspadv,ngmax). Contains the mass 
C                  fractions of species k along all mass zones which 
C                  include additional ghost cells at the "top" and 
C                  "bottom" of the standard stitched array.
C
C  dmng            double SHAPE(ngmax). Every cells mass which includes
C                  the additional ghost cells at the "top" and "bottom"
C                  of the standard stitched array.
C
C  X_stitch        double SHAPE(nspadv,nmax2). Contains the mass 
C                  fractions of species k along the stitched mass zones,
C                  i
C
C  dm_stitch,      double SHAPE(nmax2). Stitched arrays containing the
C  beta_stitch,    values of dm, beta, and gamma for every cell
C  gammaadv_stitch
C
C  v_stitch,       double SHAPE(nshellmax2). Stitched arrays containing
C  rho_stitch,     the values of velocity, density and radius at a shell
C  r_stitch
C

      implicit none

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INPUT INT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Subroutine read-in integers
      integer nspadv, nmax, nmax2, ng, ngmax, nshellmax, nshellmax2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INPUT DOUBLE ARRAYS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     stitched stream abundance and central quantities
      double precision, dimension (nspadv,nmax2) :: X_stitch
      double precision, dimension (nmax2) :: beta_stitch
      double precision, dimension (nmax2) :: dm_stitch, gammaadv_stitch

c     stitched stream abundance and cell mass arrays with ghost cells
      double precision, dimension (nspadv,ngmax) :: Xng
      double precision, dimension (ngmax) :: dmng

C     stitched shell quantities
      double precision, dimension (nshellmax2) :: v_stitch, rho_stitch
      double precision, dimension (nshellmax2) :: r_stitch

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE VARIABLES
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     looping and other useful integers
      integer i, j, k, uXngi

c     Where is my upper index in my ghost array for the upstream?
      uXngi = ngmax/2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SET BCs IN SHELL AND CENTRAL QUANTITIES
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     WORKING WITH CENTRAL QUANTITIES
C     Put X_stitch into Xng
      do k = 1,nspadv
         do i = ng+1,uXngi
C     j varies between 1 and nmax
            j = i-(ng+1)+1
            Xng(k,i) = X_stitch(k,j)
         end do
         do i = uXngi+1,ngmax-ng
C     j varies between nmax+1 and 2*nmax
            j = i-(uXngi+1)+nmax+1
            Xng(k,i) = X_stitch(k,j)
         end do
      end do

C     put dm_stitch into dmng
      do i = ng+1,ng+nmax2
C     j varies between 1 and nmax
          j = i-(ng+1)+1
          dmng(i) = dm_stitch(j)
      end do
      do i = uXngi+1,ngmax-ng
c     j varies between nmax+1 and nmax+(nmax-1+1)
          j = i-(uXngi+1)+nmax+1
          dmng(i) = dm_stitch(j)
      end do

C     Now to set the boundary conditions using the ghost cells
      dmng(1) = dmng(ngmax-ng-1)
      dmng(2) = dmng(ngmax-ng)

      dmng(ngmax-1) = dmng(ng+1)
      dmng(ngmax) = dmng(ng+2)

      Xng(:,1) = Xng(:,ngmax-ng-1)
      Xng(:,2) = Xng(:,ngmax-ng)

      Xng(:,ngmax-1) = Xng(:,ng+1)
      Xng(:,ngmax) = Xng(:,ng+2)

C     I will ensure no weird things with b and gammaadv. By reading them 
C     in they should have no horizontal transfer at the boundaries but 
C     we should make sure they are explicitly zero here
      beta_stitch(1) = 0d0
      beta_stitch(nmax) = 0d0
      beta_stitch(nmax+1) = 0d0
      beta_stitch(nmax2) = 0d0

      gammaadv_stitch(1) = 0d0
      gammaadv_stitch(nmax) = 0d0
      gammaadv_stitch(nmax+1) = 0d0
      gammaadv_stitch(nmax2) = 0d0

C     WORKING WITH SHELL QUANTITIES
C     Now to include boundary conditions. These enforce FULL mass 
C     transfer from upstream into downstream at the top and at the 
C     bottom, FULL mass transfer from downstream into the upstream

C     flux INTO top upstream cell is equal to flux OUT OF top downstream 
C     cell. We already have a mirror of this except for the last shell
C     coordinate
      v_stitch(nshellmax) = v_stitch(nshellmax-1)
      rho_stitch(nshellmax) = rho_stitch(nshellmax-1)
      r_stitch(nshellmax) = r_stitch(nshellmax-1)

C     flux INTO bottom downstream cell is equal to flux OUT OF bottom
C     upstream cell. We already have a mirror of this except for the
C     first shell coordinate
      v_stitch(1) = v_stitch(2)
      rho_stitch(1) = rho_stitch(2)
      r_stitch(1) = r_stitch(2)

      v_stitch(nshellmax2) = v_stitch(2)
      rho_stitch(nshellmax2) = rho_stitch(2)
      r_stitch(nshellmax2) = r_stitch(2)

      return

      end subroutine setBC




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




      subroutine Xslopes(slopes,Xng,dmng,nspadv,ngmax,nslopes)

C     This subroutine calculates the slopes at the interfaces (shells)
C     of all cells. This requires differences and with the ghost cells
C     the routine can be run in a single loop. The appropriate slope is
C     chosen using minmod (see ...) and upwinding

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE INPUTS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
C     Anything that is NOT defined within "advect" is shown here. If you
C     are unsure of any variable definitions in this subroutine look at
C     subroutine advect

C  slopes          double SHAPE(nspadv,nslopes). Contains the estimate of
C                  the slopes at every shell including within the "ghost"
C                  cells
C
C  Xng             double SHAPE(nspadv,ngmax). Contains the mass 
C                  fractions of species k along all mass zones which 
C                  include additional ghost cells at the "top" and 
C                  "bottom" of the standard stitched array.
C     
C  dmng            double SHAPE(ngmax). Every cells mass which
C                  includes the additional ghost cells at the "top" and 
C                  "bottom" of the standard stitched array.
C     
C  ngmax           int. The number of ghost cells
C
C  nslopes         int. The number of slope estimates including ghost
C                  cells
C     

      implicit none

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INPUT INT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Subroutine read-in integers
      integer nspadv, ngmax, nslopes

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INPUT DOUBLE ARRAYS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     stitched stream abundance and cell mass arrays with ghost cells
      double precision, dimension (nspadv,ngmax) :: Xng
      double precision, dimension (ngmax) :: dmng

C     stitched stream abundance slopes with ng-1 ghost cells
      double precision, dimension (nspadv,nslopes) :: slopes

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE VARIABLES
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     temporary double precision variables
      double precision Lslopes, Rslopes
      double precision Labs, Rabs, product

C     looping integers
      integer i, j, k

C     set all elements of slopes to 0
      slopes = 0d0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     START LOOPING FOR SLOPE ESTIMATES
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do k=1,nspadv
         do i=1,nslopes
C     First get the "Left-sided" and "Right-sided" estimates of
C     interface states
            Lslopes = (Xng(k,i+1) - Xng(k,i)) / dmng(i+1)
            Rslopes = (Xng(k,i+2) - Xng(k,i+1)) / dmng(i+1)

C     If we have a product greater than zero then the slope could will
C     be non zero
            product = Lslopes * Rslopes
            if (product .gt. 0) then

C     depending on the magnitude of slope, we would choose the smaller
C     one
               Labs = ABS(Lslopes)
               Rabs = ABS(Rslopes)
               if (Labs .lt. Rabs) then
                  slopes(k,i) = Lslopes * dmng(i+1)
               else
                  slopes(k,i) = Rslopes * dmng(i+1)
               endif

C     Slope is forced to zero
            else
               slopes(k,i) = 0d0
            endif
         end do
      end do

      return

      end subroutine Xslopes




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




      subroutine fluxes(rk,X_stitch,dm_stitch,beta_stitch,
     $     gammaadv_stitch,v_stitch,rho_stitch,r_stitch,nspadv,nmax,
     $     nmax2,nshellmax,nshellmax2,nslopes,ng,ngmax)

C     This subroutine calculates the fluxes on the cell interfaces
C     (shells) to radially advect material. There are also the
C     horizontal fluxes which contribute to mixing between the up- and
C     down-streams. This will modify a "rk" interpolation to piece
C     together the total explicit step

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE INPUTS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
C     Anything that is NOT defined within "advect" is shown here. If you
C     are unsure of any variable definitions in this subroutine look at
C     subroutine advect

C  rk              double SHAPE(nspadv,nmax2). Contains the "runge-kutta"
C                  interpolated estimate of the slope, dmXdt
C
C  X_stitch        double SHAPE(nspadv,nmax2). Contains the mass 
C                  fractions of species k along the stitched mass zones,
C                  i
C
C  dm_stitch,      double SHAPE(nmax2). Stitched arrays containing the
C  beta_stitch,    values of dm, beta, and gamma for every cell
C  gammaadv_stitch
C
C  v_stitch,       double SHAPE(nshellmax2). Stitched arrays containing
C  rho_stitch,     the values of velocity, density and radius at a shell
C  r_stitch
C
C  nslopes         int. The number of slope estimates including ghost
C                  cells
C     

      implicit none

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     INPUT INT
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
C     integers read in
      integer nspadv, nmax, nmax2, nshellmax, nshellmax2, nslopes, ng
      integer ngmax

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     INPUT DOUBLE ARRAYS
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     stitched arrays read in with species
      double precision, dimension (nspadv,nmax2) :: X_stitch, rk

C     stitched arrays read in without species
      double precision, dimension (nmax2) :: dm_stitch, beta_stitch
      double precision, dimension (nmax2) :: gammaadv_stitch
      double precision, dimension (nshellmax2) :: v_stitch, rho_stitch
      double precision, dimension (nshellmax2) :: r_stitch

c     stitched stream abundance and cell mass arrays with ghost cells
      double precision, dimension (nspadv,ngmax) :: Xng
      double precision, dimension (ngmax) :: dmng

C     stitched stream abundance slopes with ng-1 ghost cells
      double precision, dimension (nspadv,nslopes) :: slopes
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     SUBROUTINE VARIABLES
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     looping and convenience integers
      integer i, j, k, l

C     the flux is a shell quantity but for every species!
      double precision, dimension (nspadv,nshellmax2) :: flux

C     External subroutines called
      EXTERNAL setBC
      EXTERNAL Xslopes

C     Pi
      double precision :: pi
      pi = 4.0d0 * ATAN(1.0d0)

C     This will setup the boundary conditions and ghost values for
C     central and shell quantities!
      call setBC(Xng,dmng,X_stitch,dm_stitch,beta_stitch,
     $     gammaadv_stitch,v_stitch,rho_stitch,r_stitch,nspadv,nmax,
     $     nmax2,ng,ngmax,nshellmax,nshellmax2)

C     This will determine the slopes needed for the interpolated state 
C     for the mass fraction on the shells
      call Xslopes(slopes,Xng,dmng,nspadv,ngmax,nslopes)


C     Now I have everything ready to calculate the fluxes on the 
C     interface and then the rate of change of the mass of species k


C     With "upwinding" the entire array, the "left" interpolated state 
C     is always chosen for the fluxes
      do k=1,nspadv
         do i=1,nshellmax2
            flux(k,i) = 2.0d0 * pi * r_stitch(i)**2 * rho_stitch(i)
     $           * v_stitch(i) * (Xng(k,i+1) + 0.5d0 * slopes(k,i))
         end do
      end do

C     Now I determine the rate of change of the mass of species k in 
C     cell i
      do k=1,nspadv

C     This is working with the UPSTREAM ONLY
         do i=1,nmax

C     j is the index for the downstream rk
            j = nmax2 + 1 - i

C     The radial mass flux
            rk(k,i) = -(flux(k,i+1) - flux(k,i))

C     The required horizontal mass flux (can be positive or negative)
C     IF beta_stitch > 0, mass comes INTO upstream FROM downstream
            if (beta_stitch(i) > 0) then
               rk(k,i) = rk(k,i) + beta_stitch(i) * X_stitch(k,j)
C     IF beta_stitch < 0, mass LEAVES upstream and will go into 
C     downstream
            else
               rk(k,i) = rk(k,i) + beta_stitch(i) * X_stitch(k,i)
            endif

C     The additional gammaadv mass flux INTO upstream FROM downstream
            rk(k,i) = rk(k,i) + gammaadv_stitch(i) * X_stitch(k,j)

C     The additional gammaadv mass flux OUTOF upsream INTO downstream
            rk(k,i) = rk(k,i) - gammaadv_stitch(i) * X_stitch(k,i)
         end do

C     This is working with the DOWNSTREAM ONLY
         do i=1,nmax

C     j is the index for the downstream rk
            j = nmax + i

C     l is the index for the upstream rk
            l = nmax + 1 - i

C     The radial mass flux
            rk(k,j) = -(flux(k,j+1) - flux(k,j))

C     The required horizontal mass flux (can be positive or negative)
C     IF beta_stitch > 0, mass comes INTO upstream FROM downstream
            if (beta_stitch(j) > 0) then
               rk(k,j) = rk(k,j) - beta_stitch(j) * X_stitch(k,j)
C     IF beta_stitch < 0, mass LEAVES upstream and will go into 
C     downstream
            else
               rk(k,j) = rk(k,j) - beta_stitch(j) * X_stitch(k,l)
            endif

C     The additional gammaadv mass flux INTO upstream FROM downstream
            rk(k,j) = rk(k,j) - gammaadv_stitch(j) * X_stitch(k,j)

C     The additional gammaadv mass flux OUTOF upsream INTO downstream
            rk(k,j) = rk(k,j) + gammaadv_stitch(j) * X_stitch(k,l)
         end do
      end do

      return

C     This would end a single rk interpolation
      end subroutine fluxes
