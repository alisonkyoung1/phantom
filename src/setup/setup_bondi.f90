!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for Bondi flow problem, for both relativistic and non-relativistic solution
!
! :References: Liptai & Price (2019), MNRAS 485, 819-842
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - isol   : *(1 = geodesic flow  |  2 = sonic point flow)*
!   - iswind : *wind option (logical)*
!   - np     : *desired number of particles (stretch-mapping will only give this approx.)*
!   - rmax   : *outer edge*
!   - rmin   : *inner edge*
!
! :Dependencies: bondiexact, centreofmass, dim, externalforces,
!   infile_utils, io, kernel, metric_tools, options, part, physcon,
!   prompting, setup_params, spherical, stretchmap, timestep, units
!
 use dim,            only:gr,maxvxyzu
 use bondiexact,     only:get_bondi_solution,rcrit,isol,iswind
 use externalforces, only:mass1
 implicit none

 public :: setpart

 private

 real :: gamma_eos,rmax,rmin
 integer :: np

 logical, parameter :: set_boundary_particles = .false.

contains

!----------------------------------------------------------------
!+
!  setup for bondi accretion
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,        only:pi
 use metric_tools,   only:imet_schwarzschild,imetric
 use setup_params,   only:rhozero,npart_total
 use io,             only:master,fatal
 use kernel,         only:radkern
 use externalforces, only:accradius1,accradius1_hard
 use spherical,      only:set_sphere
 use options,        only:ieos,iexternalforce
 use timestep,       only:tmax,dtmax
 use centreofmass,   only:reset_centreofmass
 use units,          only:set_units,get_G_code
 use physcon,        only:pc,solarm,gg
 use part,           only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,igas,set_particle_type,iboundary
 use stretchmap,     only:get_mass_r,rho_func
 use infile_utils,   only:get_options
 use kernel,         only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer, parameter :: ntab=10000
 real               :: rhotab(ntab)
 real               :: vol,psep,tff,rhor,vr,ur
 real               :: r,pos(3),cs2,totmass,approx_m,approx_h
 integer            :: i,ierr,nx,nbound
 procedure(rho_func), pointer :: density_func
!
!-- Set code units
!
 call set_units(G=1.d0,c=1.d0)
 print*,' G in code units = ',get_G_code()

!
!--Set general parameters
!
 time           = 0.
 iexternalforce = 1
 hfact = hfact_default
 rmin = 7.
 rmax = 8.
 np   = 10000

 if (gr) then
    if (imetric/=imet_schwarzschild) call fatal('setup_bondi',&
       'You are not using the Schwarzschild metric.')
 endif

!
!-- Read things from setup file
!
 print "(/,1x,63('-'),1(/,1x,a),/,1x,63('-'),/)", 'Bondi Flow.'
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 if (gr) then
    ieos  = 2
    gamma = 5./3.
    polyk = 1.
 else
    gamma = 1.
    ieos  = 1
    cs2   = mass1/(2.*rcrit)
    polyk = cs2
 endif

 gamma_eos       = gamma             ! Note, since non rel bondi is isothermal, solution doesn't depend on gamma
 accradius1      = 0.
 accradius1_hard = 0.

 if (gr) then
    rmin = rmin*mass1
    rmax = rmax*mass1
 endif

 vol      = 4./3.*pi*(rmax**3 - rmin**3)
 nx       = int(np**(1./3.))
 psep     = vol**(1./3.)/real(nx)

 call get_rhotab(rhotab,rmin,rmax,mass1,gamma)

 density_func => rhofunc
 totmass  = get_mass_r(density_func,rmax,rmin)
 approx_m = totmass/np
 approx_h = hfact*(approx_m/rhofunc(rmin))**(1./3.)
 rhozero  = totmass/vol

 tff      = sqrt(3.*pi/(32.*rhozero))
 tmax     = 10.*tff
 dtmax    = tmax/150.

 print*,''
 print*,' Setup for gas: '
 print*,' min,max radius = ',rmin,rmax
 print*,' volume         = ',vol        ,' particle separation = ',psep
 print*,' vol/psep**3    = ',vol/psep**3,' totmass             = ',totmass
 print*,' free fall time = ',tff        ,' tmax                = ',tmax
 print*,''

!--- Add stretched sphere
 npart = 0
 npart_total = 0
 call set_sphere('closepacked',id,master,rmin,rmax,psep,hfact,npart,&
                 xyzh,rhotab=rhotab,nptot=npart_total)
 massoftype(:) = totmass/npart
 print "(a,i0,/)",' npart = ',npart

 nbound = 0
 do i=1,npart

    pos = xyzh(1:3,i)
    r = sqrt(dot_product(pos,pos))
    call get_bondi_solution(rhor,vr,ur,r,mass1,gamma)
    vxyzu(1:3,i) = vr*pos/r
    if (maxvxyzu >= 4) vxyzu(4,i) = ur

    if (set_boundary_particles) then
       if (r + radkern*xyzh(4,i)>rmax .or. r - radkern*xyzh(4,i)<rmin) then
          call set_particle_type(i,iboundary)
          nbound = nbound + 1
       else
          call set_particle_type(i,igas)
       endif
    endif

 enddo

!--- Reset centre of mass to the origin
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 npartoftype(:) = 0
 npartoftype(igas) = int(npart_total-nbound)
 npartoftype(iboundary) = nbound

end subroutine setpart

!----------------------------------------------------------------
!+
!  functional form of density profile, used as argument
!  to set_sphere
!+
!----------------------------------------------------------------
real function rhofunc(r)
 real, intent(in) :: r
 real :: rho,v,u

 call get_bondi_solution(rho,v,u,r,mass1,gamma_eos)
 rhofunc = rho

end function rhofunc

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for bondi setup routine'
 write(iunit,"(/,a)") '# solution type'
 if (gr) then
    call write_inopt(isol,'isol','(1 = geodesic flow  |  2 = sonic point flow)',iunit)
    call write_inopt(iswind,'iswind','wind option (logical)',iunit)
 endif
 call write_inopt(rmin,'rmin','inner edge',iunit)
 call write_inopt(rmax,'rmax','outer edge',iunit)
 call write_inopt(np,'np','desired number of particles (stretch-mapping will only give this approx.)',iunit)
 close(iunit)

end subroutine write_setupfile
!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 use units,        only:select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 if (gr) then
    call read_inopt(isol,  'isol',   db,ierr)
    call read_inopt(iswind,'iswind', db,ierr)
 endif
 call read_inopt(rmin, 'rmin', db,ierr)
 call read_inopt(rmax, 'rmax', db,ierr)
 call read_inopt(np,   'np',   db,ierr)
 call close_db(db)

end subroutine read_setupfile

!----------------------------------------------------------------
!+
!  prompt for parameters
!+
!----------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt

 if (gr) then
    call prompt(' Enter solution type isol (1 = geodesic | 2 = sonic point flow) ',isol,1,2)
    call prompt(' Do you want a wind (y/n)? ',iswind)
 endif
 call prompt(' Enter inner edge: ',rmin,0.)
 call prompt(' Enter outer edge: ',rmax,rmin)
 call prompt(' Enter the desired number of particles: ',np,0)

end subroutine setup_interactive

!----------------------------------------------------------------
!+
!  construct table of density as a function of radius
!+
!----------------------------------------------------------------
subroutine get_rhotab(rhotab,rmin,rmax,mass1,gamma)
 real, intent(out) :: rhotab(:)
 real, intent(in) :: rmin,rmax,mass1,gamma
 integer :: i,n
 real :: dr,r,vr,ur

 n  = size(rhotab)
 dr = (rmax-rmin)/(n-1)

 do i=1,n
    r = rmin + (i-1)*dr
    call get_bondi_solution(rhotab(i),vr,ur,r,mass1,gamma)
 enddo

end subroutine get_rhotab

end module setup
