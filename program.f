C     Run the program
      program two_stream

      implicit none

c     number of species
      integer, parameter :: nsp = 2

c     maximum dimension of up and down stream abundance and cell arrays
C     NOTE: This is modified automatically by two-streams-test.py
      integer, parameter :: nmax = 1000

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
      integer, parameter :: steps = 10

C     For every x steps, we will save an output file for plots
C     NOTE: This is modified automatically by two-streams-test.py
      integer, parameter :: saveevery = 2

C     The dt for a single step. This is constant in all models
C     NOTE: This is modified automatically by two-streams-test.py
      double precision, parameter :: dt = 0.05d0

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





