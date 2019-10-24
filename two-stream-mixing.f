!     For a single time step, dt, advect to second-order with runge-kutta3
      subroutine advect(Xup,Xdown,dm,gamma,v_shell,rho_shell,r_shell,
     $     nsp,nmax,nmax2,nshellmax,nshellmax2,ng,ngmax,nslopes,dt)

      implicit none

! *** 

!  Xup, Xdown       YPS arrays that contain separtely the abundances in the 
!                   down and up stream
!  dm               cell mass 
!  b                beta is the enforced horizontal mixing, this is input for 
!                   each time step and would come from the mixing analysis
!  gamma            same for additional horizontal mixing according to 3D hydro
!                   flow assymmetry
! v_shell,rho_shell quantities on cell faces
! r_shell
! ng                number of ghosts
! nslopes           

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     MPPNP INPUT NO-ARRAYS
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      integer nsp, nmax, nmax2, nshellmax, nshellmax2
      integer ng, ngmax, nslopes
      double precision dt

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     MPPNP INPUT ARRAYS
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     up and down quantities and central quantities, taken from mppnp
      double precision, dimension (nsp,nmax) :: Xup, Xdown
      double precision, dimension (nmax) :: gamma, dm

C     shell quantities
      double precision, dimension (nshellmax) :: v_shell, rho_shell
      double precision, dimension (nshellmax) :: r_shell

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     DERIVED ARRAYS
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     stitched stream abundance arrays
      double precision, dimension (nsp,nmax2) :: mX_stitch, X_stitch

C     Temporary rk interpolated steps
      double precision, dimension (nsp,nmax2) :: rk1, rk2, rk3

c     stitched central quantities (no species)
      double precision, dimension (nmax2) :: b_stitch, gamma_stitch
      double precision, dimension (nmax2) :: dm_stitch

C     stitched shell quantities which are used for fluxes
      double precision, dimension (nshellmax2) :: v_stitch, rho_stitch
      double precision, dimension (nshellmax2) :: r_stitch
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     SUBROUTINE VARIABLES
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     rk spatial interpolated states
      double precision, parameter :: rkf1 = 2.0d0/9.0d0
      double precision, parameter :: rkf2 = 3.0d0/9.0d0
      double precision, parameter :: rkf3 = 4.0d0/9.0d0

      integer i, j, k

C     external subroutines called
      EXTERNAL stitch
      EXTERNAL fluxes

C     create the stitched arrays
      call stitch(Xup,Xdown,dm,gamma,v_shell,rho_shell,r_shell,
     $     X_stitch,dm_stitch,b_stitch,gamma_stitch,v_stitch,rho_stitch,
     $     r_stitch,nsp,nmax,nmax2,nshellmax,nshellmax2)

C     Now I have the X_stitch, I will create my initial mass of every species in each cell
      do k=1,nsp
         mX_stitch(k,:) = X_stitch(k,:) * dm_stitch
      end do

C     Now I have the central and shell stitched arrays. These are used in
C     flux calls to get the dmXdt slopes at spatial interpolations which are stored
C     as the runge-kutta arrays

C     Our first runge-kutta interpolated step (we already have it, it is X_stitch!)
C     Since we don't have explicit time dependence, ignore the time interpolation

C     call fluxes to determine the dmXdt (rk1) of the interpolated state, X_stitch
      call fluxes(rk1,X_stitch,dm_stitch,b_stitch,gamma_stitch,
     $     v_stitch,rho_stitch,r_stitch,nsp,nmax,nmax2,nshellmax,
     $     nshellmax2,nslopes,ng,ngmax)

C     Our second runge-kutta interpolated step
C     Since we don't have explicit time dependence, ignore the time interpolation
      do k=1,nsp
         X_stitch(k,:) = (mX_stitch(k,:) + (1.0d0/2.0d0)*dt*rk1(k,:)) / 
     $        dm_stitch
      end do

C     call fluxes to determine the dmXdt (rk2) of the interpolated state, X_stitch
      call fluxes(rk2,X_stitch,dm_stitch,b_stitch,gamma_stitch,
     $     v_stitch,rho_stitch,r_stitch,nsp,nmax,nmax2,nshellmax,
     $     nshellmax2,nslopes,ng,ngmax)

C     Our last runge-kutta interpolated step
C     Since we don't have explicit time dependence, ignore the time interpolation
      do k=1,nsp
         X_stitch(k,:) = (mX_stitch(k,:) + (3.0d0/4.0d0)*dt*rk2(k,:)) /
     $        dm_stitch
      end do

C     call fluxes to determine the dmXdt (rk3) of the interpolated state, X_stitch
      call fluxes(rk3,X_stitch,dm_stitch,b_stitch,gamma_stitch,
     $     v_stitch,rho_stitch,r_stitch,nsp,nmax,nmax2,nshellmax,
     $     nshellmax2,nslopes,ng,ngmax)

C     stitch all of the rk interpolated states into our new final state
      do k=1,nsp
         mX_stitch(k,:) = mX_stitch(k,:) + rkf1*dt*rk1(k,:) +
     $        rkf2*dt*rk2(k,:) + rkf3*dt*rk3(k,:)
      end do

C     Now convert this back to what mppnp uses
      do k=1,nsp
         do i=1,nmax
            Xup(k,i) = mX_stitch(k,i) / dm_stitch(i)
         end do
         do i=1,nmax
C     We go backwards for the downstream
            Xdown(k,nmax+1-i) = mX_stitch(k,nmax+i) / dm_stitch(nmax+i)
         end do
      end do

      return

C     That's all folks!
      end subroutine advect


C     This subroutine will create all of the stitched arrays from the original arrays from mppnp
      subroutine stitch(Xup,Xdown,dm,gamma,v_shell,rho_shell,r_shell,
     $     X_stitch,dm_stitch,b_stitch,gamma_stitch,v_stitch,rho_stitch,
     $     r_stitch,nsp,nmax,nmax2,nshellmax,nshellmax2)

      implicit none

C     The integers read in for array sizes
      integer nsp, nmax, nmax2, nshellmax, nshellmax2

c     up and down quantities and central quantities
      double precision, dimension (nsp,nmax) :: Xup, Xdown
      double precision, dimension (nmax) :: gamma, dm

C     shell quantities
      double precision, dimension (nshellmax) :: v_shell, rho_shell
      double precision, dimension (nshellmax) :: r_shell

c     stitched stream abundance array
      double precision, dimension (nsp,nmax2) :: X_stitch

c     stitched central quantities (no species)
      double precision, dimension (nmax2) :: b_stitch, gamma_stitch
      double precision, dimension (nmax2) :: dm_stitch

C     stitched shell quantities which are used for fluxes
      double precision, dimension (nshellmax2) :: v_stitch, rho_stitch
      double precision, dimension (nshellmax2) :: r_stitch

C     Looping integers
      integer i, j, k

C     Pi
      double precision :: pi

C     Here is pi
      pi = 4.0d0 * ATAN(1.0d0)

C     ASSUMING ONLY WORKING WITH THE CONVECTION ZONE

C     Create the CENTRAL stitched arrays. Start with X
C     stitch flipped Xdown atop Xup into X_stitch
      do k = 1,nsp
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
         b_stitch(i) = 2.*pi*(rho_shell(i+1)*v_shell(i+1)*
     $        r_shell(i+1)**2.d0 - rho_shell(i)*v_shell(i)*
     $        r_shell(i)**2.d0)
      end do
      do i = 1,nmax
C     We go backwards for bdown, it is the same as bup
         b_stitch(nmax+i) = b_stitch(nmax+1-i)
      end do

C     similarly, stitch dm to follow that pattern
      do i = 1,nmax
         dm_stitch(i) = dm(i)
      end do
      do i = 1,nmax
C     We go backwards on dm now
         dm_stitch(nmax+i) = dm(nmax+1-i)
      end do

C     similarly, stitch gamma to follow that pattern
      do i = 1,nmax
         gamma_stitch(i) = gamma(i)
      end do
      do i = 1,nmax
C     We go backwards on gamma now
         gamma_stitch(nmax+i) = gamma(nmax+1-i)
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
C     We go backwards in shell but skip the top of the stream (don't duplicate)
C     therefore j goes from nmax to 1
         j = nshellmax-i
         v_stitch(nshellmax+i) = v_shell(j)
         rho_stitch(nshellmax+i) = rho_shell(j)
         r_stitch(nshellmax+i) = r_shell(j)
      end do

      return

      end subroutine stitch


C     This subroutine will add the ghost cells to the X and dm stitched
C     arrays so that the slope estimates can be easily calculated as well
C     as any other boundary conditions on quantities
      subroutine setBC(Xng,dmng,X_stitch,dm_stitch,b_stitch,
     $     gamma_stitch,v_stitch,rho_stitch,r_stitch,nsp,nmax,
     $     nmax2,ng,nng,ngmax,nshellmax,nshellmax2)

      implicit none

C     Subroutine read-in integers
      integer nsp, nmax, nmax2, ng, nng, ngmax, nshellmax, nshellmax2

c     stitched stream abundance and central quantities
      double precision, dimension (nsp,nmax2) :: X_stitch
      double precision, dimension (nmax2) :: dm_stitch
      double precision, dimension (nmax2) :: b_stitch, gamma_stitch

c     stitched stream abundance and cell mass arrays with ghost cells
      double precision, dimension (nsp,ngmax) :: Xng
      double precision, dimension (ngmax) :: dmng

C     stitched shell quantities
      double precision, dimension (nshellmax2) :: v_stitch, rho_stitch
      double precision, dimension (nshellmax2) :: r_stitch

C     other useful integers
      integer i, j, k, uXngi

c     Where is my upper index in my ghost array for the upstream?
      uXngi = nng/2

C     WORKING WITH CENTRAL QUANTITIES
C     Put X_stitch into Xng
      do k = 1,nsp
         do i = ng+1,uXngi
C     j varies between 1 and nmax
            j = i-(ng+1)+1
            Xng(k,i) = X_stitch(k,j)
         end do
         do i = uXngi+1,nng-ng
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
      do i = uXngi+1,nng-ng
c     j varies between nmax+1 and nmax+(nmax-1+1)
          j = i-(uXngi+1)+nmax+1
          dmng(i) = dm_stitch(j)
      end do

C     Now to set the boundary conditions using the ghost cells
      dmng(1) = dmng(nng-ng-1)
      dmng(2) = dmng(nng-ng)

      dmng(nng-1) = dmng(ng+1)
      dmng(nng) = dmng(ng+2)

      Xng(:,1) = Xng(:,nng-ng-1)
      Xng(:,2) = Xng(:,nng-ng)

      Xng(:,nng-1) = Xng(:,ng+1)
      Xng(:,nng) = Xng(:,ng+2)

C     I will ensure no weird things with b and gamma. By reading them in they
C     should have no horizontal transfer at the boundaries but we should make
C     sure they are explicitly zero here
      b_stitch(1) = 0d0
      b_stitch(nmax) = 0d0
      b_stitch(nmax+1) = 0d0
      b_stitch(nmax2) = 0d0

      gamma_stitch(1) = 0d0
      gamma_stitch(nmax) = 0d0
      gamma_stitch(nmax+1) = 0d0
      gamma_stitch(nmax2) = 0d0

C     WORKING WITH SHELL QUANTITIES
C     Now to include boundary conditions. These enforce FULL mass transfer from upstream
C     into downstream at the top and at the bottom, FULL mass transfer from downstream into
C     the upstream

C     flux INTO top upstream cell is equal to flux OUT OF top downstream cell
C     We already have a mirror of this except for the last shell coordinate
      v_stitch(nshellmax) = v_stitch(nshellmax-1)
      rho_stitch(nshellmax) = rho_stitch(nshellmax-1)
      r_stitch(nshellmax) = r_stitch(nshellmax-1)

C     flux INTO bottom downstream cell is equal to flux OUT OF bottom upstream
C     cell. We already have a mirror of this except for the first shell coordinate
      v_stitch(1) = v_stitch(2)
      rho_stitch(1) = rho_stitch(2)
      r_stitch(1) = r_stitch(2)

      v_stitch(nshellmax2) = v_stitch(2)
      rho_stitch(nshellmax2) = rho_stitch(2)
      r_stitch(nshellmax2) = r_stitch(2)

      return

      end subroutine setBC



      subroutine Xslopes(slopes,Xng,dmng,nsp,ngmax,nslopes)

      implicit none
      integer i, j, k, nsp, ngmax, nslopes

c     stitched stream abundance and cell mass arrays with ghost cells
      double precision, dimension (nsp,ngmax) :: Xng
      double precision, dimension (ngmax) :: dmng

C     stitched stream abundance slopes with ng-1 ghost cells
      double precision, dimension (nsp,nslopes) :: slopes

C     looping values
      logical product_bool, Lbool, Rbool
      double precision Lslopes, Rslopes
      double precision Labs, Rabs, product

      slopes = 0d0

      do k=1,nsp
         do i=1,nslopes
C     First get the "Left-sided" and "Right-sided" estimates of interface states
            Lslopes = (Xng(k,i+1) - Xng(k,i)) / dmng(i+1)
            Rslopes = (Xng(k,i+2) - Xng(k,i+1)) / dmng(i+1)

C     If we have a product greater than zero then the slope could will be non zero
            product = Lslopes * Rslopes
            if (product .gt. 0) then

C     depending on the magnitude of slope, we would choose the smaller one
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


C     This subroutine will calculate the mass fluxes at the shell interfaces
C     so that rk3 can be run on it
      subroutine fluxes(rk,X_stitch,dm_stitch,b_stitch,gamma_stitch,
     $     v_stitch,rho_stitch,r_stitch,nsp,nmax,nmax2,nshellmax,
     $     nshellmax2,nslopes,ng,ngmax)

      implicit none

C     integers read in
      integer nsp, nmax, nmax2, nshellmax, nshellmax2, nslopes, ng
      integer ngmax

C     looping and convenience integers
      integer i, j, k, l, nng

C     stitched arrays read in with species
      double precision, dimension (nsp,nmax2) :: X_stitch, rk

C     stitched arrays read in without species
      double precision, dimension (nmax2) :: dm_stitch, b_stitch
      double precision, dimension (nmax2) :: gamma_stitch
      double precision, dimension (nshellmax2) :: v_stitch, rho_stitch
      double precision, dimension (nshellmax2) :: r_stitch

C     the flux is a shell quantity but for every species!
      double precision, dimension (nsp,nshellmax2) :: flux

c     stitched stream abundance and cell mass arrays with ghost cells
      double precision, dimension (nsp,ngmax) :: Xng
      double precision, dimension (ngmax) :: dmng

C     stitched stream abundance slopes with ng-1 ghost cells
      double precision, dimension (nsp,nslopes) :: slopes

C     Pi
      double precision :: pi

C     External subroutines called
      EXTERNAL setBC
      EXTERNAL Xslopes

      pi = 4.0d0 * ATAN(1.0d0)

c     number of cells in stiched arrays with ghost cells
C     This is left here as a reminder that nng != ngmax for cases
C     where we include more than the convection zone!
      nng = 2*(nmax+ng)

C     This will setup the boundary conditions and ghost values for
C     central and shell quantities!
      call setBC(Xng,dmng,X_stitch,dm_stitch,b_stitch,
     $     gamma_stitch,v_stitch,rho_stitch,r_stitch,nsp,nmax,
     $     nmax2,ng,nng,ngmax,nshellmax,nshellmax2)

C     This will determine the slopes needed for the interpolated state for the
C     mass fraction on the shells
      call Xslopes(slopes,Xng,dmng,nsp,ngmax,nslopes)


C     Now I have everything ready to calculate the fluxes on the interface and then
C     the rate of change of the mass of species k


C     With "upwinding" the entire array, the "left" interpolated state is always
C     chosen for the fluxes
      do k=1,nsp
         do i=1,nshellmax2
            flux(k,i) = 2.0d0 * pi * r_stitch(i)**2 * rho_stitch(i)
     $           * v_stitch(i) * (Xng(k,i+1) + 0.5d0 * slopes(k,i))
         end do
      end do

C     Now I determine the rate of change of the mass of species k in cell i
      do k=1,nsp

C     This is working with the UPSTREAM ONLY
         do i=1,nmax

C     j is the index for the downstream rk
            j = nmax2 + 1 - i

C     The radial mass flux
            rk(k,i) = -(flux(k,i+1) - flux(k,i))

C     The required horizontal mass flux (can be positive or negative)
C     IF b_stitch > 0, mass comes INTO upstream FROM downstream
            if (b_stitch(i) > 0) then
               rk(k,i) = rk(k,i) + b_stitch(i) * X_stitch(k,j)
C     IF b_stitch < 0, mass LEAVES upstream and will go into downstream
            else
               rk(k,i) = rk(k,i) + b_stitch(i) * X_stitch(k,i)
            endif

C     The additional gamma mass flux INTO upstream FROM downstream
            rk(k,i) = rk(k,i) + gamma_stitch(i) * X_stitch(k,j)

C     The additional gamma mass flux OUTOF upsream INTO downstream
            rk(k,i) = rk(k,i) - gamma_stitch(i) * X_stitch(k,i)
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
C     IF b_stitch > 0, mass comes INTO upstream FROM downstream
            if (b_stitch(j) > 0) then
               rk(k,j) = rk(k,j) - b_stitch(j) * X_stitch(k,j)
C     IF b_stitch < 0, mass LEAVES upstream and will go into downstream
            else
               rk(k,j) = rk(k,j) - b_stitch(j) * X_stitch(k,l)
            endif

C     The additional gamma mass flux INTO upstream FROM downstream
            rk(k,j) = rk(k,j) - gamma_stitch(j) * X_stitch(k,j)

C     The additional gamma mass flux OUTOF upsream INTO downstream
            rk(k,j) = rk(k,j) + gamma_stitch(j) * X_stitch(k,l)
         end do
      end do

      return

C     This would end a single rk interpolation
      end subroutine fluxes


