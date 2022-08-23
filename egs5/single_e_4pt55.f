!-----------------------------------------------------------------------
! Step 1: Initialization
!-----------------------------------------------------------------------

      implicit none

!     Main EGS header file      
      include 'include/egs5_h.f'

!     bounds contains the COMMON block BOUNDS with the variables
!       ECUT: An array of region-dependent charged particle cutoff
!       energies in MeV.
!       PCUT: An array of region-dependent photon cutoff energies in
!       MeV.
      include 'include/egs5_bounds.f'
!     media contains the COMMON block MEDIA with the variables
!       NMED number of media being used
!       MEDIA Array containing the density of the media in g/cm3
!       IRAYLM Array of flags for turning on Rayleigh scattering in the
!       various regions.
!       CHARD Array used to define the dimensions in cm (default 1 cm)
!       MED Array containing the medium index for each region
      include 'include/egs5_media.f'
!     usersrc contains the COMMON block USERSC with the variable 
!       EMAXE: The maximum total energy of any electron in the
!       simulation
      include 'include/egs5_usersc.f'
!     edge contains the COMMON block EDGE2 with the variable IEDGFL     
      include 'include/egs5_edge.f'
!     misc contains the COMMON block MISC with the variable IMPACR
      include 'include/egs5_misc.f'
!     randomm contains the COMMON block RLUXDAT with the variables
!       INSEED: The initial seed
!       LUXLEV: Luxury level used by RANLUX
      include 'include/randomm.f'
!     useful contains the COMMON block USEFUL with the variable
!       RM: Rest mass of the electron
      include 'include/egs5_useful.f'
!
      include 'include/egs5_epcont.f'

!     Set the target thickness in cm. This value will be used by the
!     subroutine howfar
      common/geo/tgtdz
      real*8 tgtdz

      common/sim/enp, emax
      integer enp
      real*8 emax

      real*8 ei, ekin, xi, yi, zi, ui, vi, wi, wti
      integer iqi, iri

!     Initialize the array that will contain the names of the media. In
!     this case, the array is of length 1 since only the target material
!     needs to be defined.
!     NOTE: The names of the media need to be exactly 24 characters
!     long.
      integer i, j
      character*24 medarr(1)

      open(UNIT=6, FILE='egs5job.out', STATUS='UNKNOWN')
      open(UNIT=100, FILE='electron.csv', STATUS='UNKNOWN')
    
      print *, 'Preparing for PEGS5 call.'
!-----------------------------------------------------------------------
! Step 2: pegs5-call
!-----------------------------------------------------------------------
!     Initialize some general variables
      call block_set                 

!     Define the media before calling PEGS5
      nmed=1
      medarr(1)='W-RAY                   '

      do j=1,nmed
        do i=1,24
          media(i,j)=medarr(j)(i:i)
        end do
      end do

      chard(1) = 0.0001

      call pegs5

!-----------------------------------------------------------------------
! Step 3: Pre-hatch-call-initialization
!-----------------------------------------------------------------------

!     The number of regions in the geometry. In this case, there will be
!     three regions: | Vacuum | W | Vacuum | 
      nreg = 3

!     Specify what type of media fills each region
      med(1) = 0
      med(2) = 1
      med(3) = 0

!     Terminate electron histories at 0.521 MeV in the W target
      ecut(2) = 0.521
!     Terminate photon histories at 0.001 MeV in the W target
      pcut(2) = 0.001
!     This turns on explicit modeling of K and L-edge fluorescent
!     photons.
      iedgfl(2) = 1
!     This turns on Rayleigh scattering
      iraylr(2) = 1
!     This turns on electron impact ionization
      impacr(2) = 1

!     --------------------------------------------------------
!     Random number seeds.  Must be defined before call hatch
!     or defaults will be used.  inseed (1- 2^31)
!     --------------------------------------------------------
      inseed = 1
!     The "luxury level" sets the level of tests for randomness.
      luxlev = 1

!     Initialize the ranlux random-number generator
      call rluxinit

!-----------------------------------------------------------------------
! Step 4:  Determination-of-incident-particle-parameters
!-----------------------------------------------------------------------

!     Specify that the incident particle is an electron
      iqi = -1
!     Set the coordinates of the incident particle
      xi = 0.0
      yi = 0.0
      zi = 0.0
!     Set the direction cosines
      ui = 0.0
      vi = 0.0
      wi = 1.0
!     Set the incident region to be 2 
      iri = 2
!     Weight factor in importance sampling
      wti = 0
!     Total energy of the incident particle in MeV
      ei=4550.D0
      ekin=ei+iqi*RM

!-----------------------------------------------------------------------
! Step 5:   hatch-call
!-----------------------------------------------------------------------

!     Maximum total energy
      if (iqi.ne.0) then
        emaxe = ei
      else
        emaxe = ei + RM
      end if

      open(UNIT=KMPI,FILE='pgs5job.pegs5dat',STATUS='old')
      open(UNIT=KMPO,FILE='egs5job.dummy',STATUS='unknown')

      call hatch

      close(UNIT=KMPI)
      close(UNIT=KMPO)
!-----------------------------------------------------------------------
! Step 6:  Initialization-for-howfar
!-----------------------------------------------------------------------
      tgtdz = 0.002

!-----------------------------------------------------------------------
! Step 7:  Initialization-for-ausgab
!-----------------------------------------------------------------------

      iausfl(8) = 1

      write(100, *) 'Energy (GeV), Theta z (mrad)'
      do i=1,10000000
!        write(100,*) 'Event: ', i
        enp = 1
        emax = 0
        call shower(iqi, ei, xi, yi, zi, ui, vi, wi, iri, wti)
      end do

      close(UNIT=100)
      close(UNIT=7) 

      stop
      end

!
!
!
      subroutine ausgab(iarg)

      implicit none
!     Main EGS header file
      include 'include/egs5_h.f'
!     The COMMON block STACK contains the following useful variables
!       x, y, z: Position of a particle 
!       u, v, w: Direction cosines of a particle
!       ir: Index of particle's current region
!       np: The particle currently being pointed to
      include 'include/egs5_stack.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_epcont.f'

!     Arguments
      integer iarg

      real*8 kine, thetaz
      
      common/sim/enp, emax
      integer enp
      real*8 emax
      integer i

!      write(100, *) 'index:', np, 'energy:', e(np), 'charge:', iq(np)
!      write(100, *) 'time:', time(np), 'eloss:', deinitial

      if (iarg.eq.7) then
        do i = 1, np
          if (iq(i).eq.-1.and.e(i).gt.emax) then
            enp = i
            emax = e(i)
          end if
        end do
!        write(100, *) 'Brem has occurred.'
!        write(100, *) 'index: ', np, 'energy:', e(np), 'charge:', iq(np)
!        write(100, *) 'i-1:',np-1,'energy:',e(np-1),'charge:',iq(np-1)
!        if (e(np).gt.3000) then
        return
      end if
!      if (np.ne.1) then
!        write(100, *) 'Returning index: ', np
!        return
!      end if

      if (iarg.eq.3.and.ir(np).eq.3.and.iq(np).eq.-1.and.np.eq.enp) then

!          write(100, *) 'Particle index: ', np
!          write(100, *) 'Energy: ', e(np), 'Time: ', time(np)
          
          kine = e(np) - RM
          thetaz = acos(w(np))
          write(100, *) kine, ',', thetaz 
      end if
      return
      end
!
!
!

      subroutine howfar

      implicit none

!     Main EGS header file
      include 'include/egs5_h.f'
!     COMMONs required by EGS5 code
      include 'include/egs5_epcont.f'
!     The COMMON block STACK contains the following useful variables
!       x, y, z: Position of a particle 
!       u, v, w: Direction cosines of a particle
!       ir: Index of particle's current region
!       np: The particle currently being pointed to
      include 'include/egs5_stack.f'

!     The target thickness. 
      common/geo/tgtdz
      real*8 tgtdz
    
!     The distance from the particle to the target boundary
      real*8 d

!     First, check if the particle has exited the target i.e. not in
!     region 2. If it's outside the target, discard the particle.
      if ( ir(np).eq.3 ) then
!     Setting idisc to a non-zero value tells the simulation to
!     discard the particle.       
         
        idisc = 1
        return
!     Now consider the case where the particle is in the target and
!     traversing forward
      else if ( ir(np).eq.2 ) then
        if ( w(np).gt.0.0) then
!     Calculate the distance from the particle to the downstream face of
!     the target (d).  If the  distance is greater than the step size
!     (ustep) then continue stepping. Otherwise, set ustep to d and set
!     the new region (irnew) to 3 to indicate the particle will exit.
          d = (tgtdz - z(np))/w(np)
          if (d.gt.ustep) then
            return
          else
            ustep = d
            irnew = 3
            return
          end if
!     In the case where the particle is in the target and traversing
!     backwards, let the particle continue stepping if the step size is
!     smaller than the distance to the upstream face of the target.
        else if ( w(np).lt.0.0 ) then
          d = -z(np)/w(np)
          if ( d.gt.ustep ) then
            return
          else 
            ustep = d
            irnew = 1
            return
          end if
        else if ( w(np).eq.0.0 ) then
          return
        end if
!     If the particle is upstream of the target, then it's the incident
!     particle.  In this case, as long as the particle is pointed
!     towards the target, set the position to the upstream face of the
!     target and let it propagate.
      else if ( ir(np).eq.1 ) then
        if ( w(np).gt.0.0 ) then
          ustep = 0.0
          irnew = 2
          return
!     The particle must have been reflected so discard it
        else
          idisc = 1
          return
        end if
      end if 
      end
