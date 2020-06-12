      subroutine monte_step(cartstatehandle,ministatehandle,
     1 nhm, nihm, nlig,
     2 ens, phi, ssi, rot, xa, ya, za, morph, dlig,
     3 locrests, has_locrests, ! locrest: location restrains !rex:/transfer/attract -> main branch (bazaar-versionwerw.)
     4 seed, label,
     5 gesa, energies, lablen) !gesa: gesamtenergie, vdw, coulomb, restrainste-E (zB beschränkung der normlmoden, anderes internes), gravity <-> density mask-E  4 oder 5
c
c  variable metric minimizer (Harwell subroutine lib.  as in Jumna with modifications)
c     minimizes a single structure
      
      implicit none
      logical mcdebug
      parameter (mcdebug = .FALSE.)
c     Parameters
      integer cartstatehandle,ministatehandle
      include 'max.fin' ! vars like "maxmover"->3,"maxlig"->100
      integer nlig, seed
      real*8 locrests
      dimension locrests(3,maxlig)
      integer has_locrests
      dimension has_locrests(maxlig)
      real *8 gesa, energies, vbias
      dimension energies(6)
      integer lablen
      character label
      dimension label(lablen)
      integer mcprobs
      dimension mcprobs(maxmover)

      integer nhm
      dimension nhm(maxlig)
      integer nihm
      dimension nihm(maxlig)
      integer ens, scaleens
      dimension ens(maxlig)
      real*8 phi, ssi, rot, dlig, xa, ya, za, morph
      dimension phi(maxlig), ssi(maxlig), rot(maxlig)
      dimension dlig(maxmode+maxindexmode, maxlig)
      dimension xa(maxlig), ya(maxlig), za(maxlig)
      dimension morph(maxlig)

ccc----------------ZHE--------------------
      integer mover, trials
      real*8 cumu_ws,sws
      dimension cumu_ws(maxmover)  ! maxmover=3 before
ccc---------------------------------------
c   JULIAN 
c    !  real(8) energy  ! "!" at beginning not working...
c JULIAN
      

c     Local variables
      real*8 enew, energies0
      real*8 rrot1,rrot2,rrot3,rrot4,sphi,sssi,srot
      dimension energies0(6)
c     integer dseed,i,ii,j,jj,k,kk,itr,nfun
      integer i,ii,j,jj,k,kk,itr,nfun
      integer itra, ieig, iindex, iori, fixre, iscore,imcmax
      integer ju,ju0,jl,jb,nmodes,nimodes, jn, jn0
      integer iab,ijk,iaccept,accepts,inner_i
      real*8 xnull
      real*8 scalecenter,scalemode,ensprob,scalerot,rr
      real*8 rotmat,randrot,newrot,sum
      real*8 xaa,delta,deltamorph,bol,pi,mctemp,accept_rate
      integer ensaa
      real*8 dseed
      dimension xaa(maxdof)
      dimension ensaa(maxlig)
      dimension delta(maxdof), deltamorph(maxlig)
      dimension rr(maxdof),randrot(0:9),rotmat(0:9),newrot(0:9)
      integer nrens
      dimension nrens(maxlig)
      pointer(ptr_nrens,nrens)
      real*8 neomorph
      integer, parameter :: ERROR_UNIT = 0
      pi=3.141592654d0


      
      
      do i=1, maxlig ! für alle ligande mit 0 init.
      ensaa(i) = 0  !ensemble-was?
      enddo
c        call print_struc2(seed,label,gesa,energies,nlig,
c     1  ens,phi,ssi,rot,xa,ya,za,locrests,morph,
c     2  nhm,nihm,dlig,has_locrests,lablen)

      call ministate_f_monte_step(ministatehandle, ! sets vars from ministate refed by ministatehandle (in ministate.cpp)
     1 iscore,imcmax,iori,itra,ieig,iindex,fixre,mctemp,
     2 scalerot,scalecenter,scalemode,ensprob,mcprobs) 

!     cccc-------------------ZHE-------------------- +Julian
!C$$$  nmover = 3 [rigid_body_mover, ensemble_mover, both]
!C$$$
!c       p für wahrscheinlichkeit, dass rigid/ensemble/both, p ist jew. x/GESAMTZAHL 



!      cumu_ws(1) = 70;
!      cumu_ws(2) = 0;   
!      cumu_ws(3) = 0;
!      cumu_ws(4) = 30;
      cumu_ws=mcprobs;
      sws = 0;
      do i=1,maxmover ! =3
         sws = sws+cumu_ws(i)
      enddo
      cumu_ws(1) = cumu_ws(1)/sws
      do i=2,maxmover
         cumu_ws(i) = cumu_ws(i-1)+cumu_ws(i)/sws ! sets upper bounds for  its range
      enddo
!      end ZHE
      
c     always calculate only energies
      iab = 0

c
c  all variables without lig-hm ! -> ligand-hormonic-modes?
c
      jb=3*iori*(nlig-fixre)+3*itra*(nlig-fixre) ! jb and jl are what?
c  all variables including lig-hm
      nmodes = 0
      nimodes = 0
      do 5 i=fixre, nlig
      nmodes = nmodes + nhm(i)
      nimodes = nimodes + nihm(i)
    5 continue
      ju=jb+ieig*nmodes
      jn = ju + iindex*nimodes
      ju0=ju
      jn0 = jn
      do i=1,nlig
      if (morph(i).ge.0) then
      jn = jn + 1
      endif
      enddo

c  only trans or ori
      jl=3*iori*(nlig-fixre)        ! fixre is 0 according to ministate.cpp, iori=1, nlig probably count of ligands (different ligand-positions)

      call ministate_calc_pairlist(ministatehandle,cartstatehandle)
      call cartstate_get_nrens(cartstatehandle, ptr_nrens)

      xnull=0.0d0
      accept_rate=0.0d0
      dseed=seed
!      dseed=12345   !seed konstant gesetzt, normal für jede strukt, müsste normal array sein
      accepts=0
      scalerot=pi*scalerot/180d0 ! scalerot=0.05 according to ministat.cpp -> 5% einer Halb-Kreisbewegung?
      nfun=0
      itr=0
      if (iscore.eq.1) then
        iori = 1
        itra = 1
      endif
c intial energy evaluation
c
      call energy(cartstatehandle,ministatehandle,
     1 iab,iori,itra,ieig,iindex,fixre,
     2 ens,phi,ssi,rot,xa,ya,za,morph,dlig,
     3 locrests, has_locrests, seed,
     4 gesa,energies,delta,deltamorph)
       
c   start Monte Carlo
c      write(*,*) "start MonteCarlo" ! will show up after all printed structures!
      trials=0
      iaccept=1
c      print*, "vmax: ", imcmax
      do 4000 ijk=1,imcmax          ! loop ends at linelabel "4000" -> ganz unten
c         do 4100 inner_i=1, 1000
c      write (ERROR_UNIT,*), ijk, imcmax
c store old Euler angle, position and ligand and receptor coordinates
c
c phi,ssi,rot for first molecule are fixed!

      if(iaccept.eq.1) then

      do i=1,nlig
         ensaa(i)=ens(i) ! sets energy-values for every ligand in ensemble?
      enddo
      if(iori.eq.1) then
      do 118 i=1+fixre,nlig !fixre is 0 or 1, wird evtl bewegt (0), falls true=fix=1: starte erst bei lig2 mit bwegung
      ii=3*(i-fixre-1)
      xaa(ii+1)=phi(i)  !xaa just stores old values for all ligands? Nein, da werden veränderungen drauf gemacht -> passt echt nicht zu iaccept... eher doch store
      xaa(ii+2)=ssi(i)
      xaa(ii+3)=rot(i)
  118 continue
      endif
      if(itra.eq.1) then
      do 122 i=1+fixre,nlig
      ii=jl+3*(i-fixre-1)   ! wie oben, nur um "jl" versetzt? -> jr ist wohl tatsächlich  meist 3*c_lig -> Rotations-"plätze"
      xaa(ii+1)=xa(i)
      xaa(ii+2)=ya(i)
      xaa(ii+3)=za(i)
  122 continue
      endif
      jj = jn0
      do i=1,nlig
       if (morph(i).ge.0) then  ! what is "morph(i)! ? seems to be calculated in energy(...)
        xaa(jj+1) = morph(i)
        jj = jj + 1
       endif
      enddo

c if ligand flex is included store deformation factor in every mode in dlig(j)

      if(ieig.eq.1) then
      jj = 0
      do 130 j=1,nlig
      do 131 i=1,nhm(j) !nhm nr of normal harmonic modes
      xaa(jb+jj+i)=dlig(i,j)    !dlig änderung in hedem F° WERDEN hier gespeichert
  131 continue
      jj = jj + nhm(j)
  130 continue
      endif

      if(iindex.eq.1) then
      jj = 0
      do 140 j=1,nlig
      do 141 i=1,nihm(j)    !iattract harmonic modes, interface-atom-bewegung, alls mit i- ist iattract
      xaa(ju0+jj+i)= dlig(ju0+i,j)
  141 continue
      jj = jj + nihm(j)
  140 continue
      endif

      endif
c       end of storing old values (if last step was accepted)

c old Cartesians are not stored! -> those in cartstate.x?

c generate a total of ju random numbers
      call GGUBS_STEP(dseed,2,rr)
        
      call random_mover(rr(1),cumu_ws,mover,maxmover) !chooses from [0:1] intervall according to rr(1) the approbriate nr for kind-of-move, stores nr in "mover"
c      print*,"mover ",mover

      if ( mover.eq.1 ) trials = trials+1

      if (mover.eq.1 .or. mover.eq.3) then ! not compativle with sidechain atm  .or. ensprob.eq.0) then ! move tranlational/rotational (rigid body), always do this if no ensemble_mover allowed (ensprob=0)
           trials = trials+1
           call rigid_body_mover(nlig,jl,iori,itra,phi,ssi,rot,xa,ya,za,
     1          fixre,scalerot,scalecenter,dseed)
     


      endif
      if (mover.eq.2 .or. mover.eq.3) then    ! HOW DOES THIS WORK? wohl ensemble-austausch, ensembles müssen in reihenfolge ud atomzahl gleich sien
           call GGUBS_STEP(dseed,3,rr)
           do i=1,nlig  
              if (nrens(i).gt.0.and.morph(i).lt.0) then
C                call GGUBS(dseed,3,rr)
                 if (rr(1).lt.ensprob.and.rr(3).lt.float(i)/nlig) then
c	    ens(i) = int(rr(2)*nrens(i))+1
                    call enstrans(cartstatehandle,i-1,ens(i),rr(2),
     2                   ens(i))
                    exit
                 endif
              endif
           enddo
      endif
      if (mover.eq.4) then
      error stop "trying to use unimplemented SidechainSwitcher"
      !call GGUBS(dseed,2,rr)
      !call sidechain_switcher(cartstatehandle,rr(1),rr(2))
      endif
c make a move in HM direction and update x, y(1,i) and y(2,i) and dlig(j)
c     call crand(dseed,ju+1,rr)
      call GGUBS_STEP(dseed,jn+1,rr)
c     dseed = int(10000*rr(ju+1))
      if(ieig.eq.1) then
      kk = 0
      do 1180 k=1,nlig
      do 1200 i=1,nhm(k)
      dlig(i,k)=xaa(i+jb+kk)+scalemode*(rr(i+jb+kk)-0.5d0) ! neue Koords in dlig!!
 1200 continue
      kk = kk + nhm(k)
 1180 continue
      endif
      if(iindex.eq.1) then
      kk = 0
      do 1280 k=1,nlig          
      do 1300 i=1,nihm(k)
      dlig(ju0+i,k)=xaa(i+ju0+kk)+scalemode*(rr(i+ju0+kk)-0.5d0)
 1300 continue
      kk = kk+ nihm(k)
 1280 continue
      endif
c rigid body move, translation and rotation

c      call rigid_body_mover(nlig,jl,iori,itra,phi,ssi,rot,xa,ya,za,
c     1 fixre,scalerot,scalecenter,dseed)

      jj = jn0
      do i=1,nlig
       if (morph(i).ge.0) then  ! morph is legacy code, ignore
        neomorph = morph(i)+scalemode*(0.5d0-rr(ii+1))
        if (neomorph.lt.0) neomorph = 0
        if (neomorph.gt.nrens(i)-1.001) neomorph = nrens(i)-1.001
        morph(i) = neomorph
        jj = jj + 1
       endif
      enddo
      
c   energy calc for this new state
      call energy(cartstatehandle,ministatehandle,      ! sth. here REVERSES MY COORD_CHANGES !!!!
     1 iab,iori,itra,ieig,iindex,fixre,
     2 ens,phi,ssi,rot,xa,ya,za,morph,dlig, !again: dlig are changes to coords
     3 locrests, has_locrests, seed,
     4 enew,energies0,delta,deltamorph)
       
     
     
c  new energy
c      write (ERROR_UNIT,*),'Energy2', enew

c if using wte, bias energy
c      call evaluate_bias( enew, vbias);
                ! make metropolis crit.
      bol=enew-gesa
c      bol = enew + vbias - gesa
      if (mctemp.eq.0) then
      bol=sign(1.0d0,-bol) ! 1 or -1, depending if delta-en ("bol") negative (good) or positive (bad) ->energies are positive if bad...
      else
      bol=exp(-bol/mctemp)
      endif
c      write(*,*)'exp(bol)',enew,gesa,enew-gesa,bol
c     call crand(dseed,2,rr)        ! old (bad) random number generator?
      call GGUBS_STEP(dseed,2,rr)
c     dseed = int(10000*rr(2)) ! 2nd random number WAS used to get new seed
      if(bol.gt.rr(1)) then         ! !!!! LOOP wenn angenommen, always accepted if bol>=1, and that is >=1 with an energy improvement (exp(postive)=...)
      
c      write(ERROR_UNIT,*)'accept the step', bol, rr(1)

c     write(*,*)
c    1 'rrot1,rrot2,rrot3,rrot4,sphi,phi(i),sssi,ssi(i),srot,rot(i)',
c    2 rrot1,rrot2,rrot3,rrot4,sphi,phi(2),sssi,ssi(2),srot,rot(2)
c      gesa=enew+vbias
c      call update_bias( enew );

c energien übernhemen

        gesa=enew
        energies(:)=energies0(:)
      iaccept=1
      if (mcdebug) then 
      write(ERROR_UNIT,*)' accepted, c_acs:',accepts," mover:",mover
      endif 
      if (mover.eq.1.or.mover.eq.3 ) accepts=accepts+1    ! only +1 when mover==1? what about mover==3? ->added
      
c overwrite old xaa variables, see above
      else      !why has "not accepting" many more commands than "accepting" ? -> acc-rate should be 0.3, why not other way around? -> acc-steps are done above (storing in xaa...)
c do not overwrite xaa variables
       if (mcdebug) then 
      write(ERROR_UNIT,*)' step rejected c_acs:',accepts," mover:",mover
      endif
      iaccept=0
      if (mover.eq.4) then
      if (mcdebug) then 
            write(ERROR_UNIT,*) " reversing switch"
      endif 
      endif
      do i=1,nlig
         ens(i)=ensaa(i)
      enddo
      if(iori.eq.1) then
      do 1118 i=1+fixre,nlig
      ii=3*(i-fixre-1)
      phi(i)=xaa(ii+1)
      ssi(i)=xaa(ii+2)
      rot(i)=xaa(ii+3)
 1118 continue
      endif
      if(itra.eq.1) then
      do 1122 i=1+fixre,nlig
      ii=jl+3*(i-fixre-1)
      xa(i)=xaa(ii+1)
      ya(i)=xaa(ii+2)
      za(i)=xaa(ii+3)
 1122 continue
      endif

c if ligand flex is included store deformation factor in every mode in dlig(j)

      if(ieig.eq.1) then ! ob harmonische modes -> eig=eigenvektoren, des elastic network models
      jj = 0
      do 230 j=1,nlig
      do 231 i=1,nhm(j)
      dlig(i,j)=xaa(jb+jj+i)
  231 continue
      jj = jj + nhm(j)
  230 continue
      endif

      if(iindex.eq.1) then ! ob iattract
      jj = 0
      do 240 j=1,nlig
      do 241 i=1,nihm(j)
      dlig(ju0+i,j)=xaa(ju0+jj+i)
  241 continue
      jj = jj + nihm(j)
  240 continue
      endif
      endif

      jj = jn0
      do i=1,nlig
       if (morph(i).ge.0) then
        morph(i) = xaa(jj+1)
        jj = jj + 1
       endif
      enddo

c 4100 continue
c      print*,"accepts ", accepts

      if(mod(trials,25)==0) then
         accept_rate = real(accepts)/trials
         if(accept_rate.gt.0.3) then
            scalecenter=scalecenter*1.1         !! mache größere schritte, wenn accRate hoch
	    if ( scalerot*1.1<pi ) then
c control scalerot range in [0 pi]
               scalerot=scalerot*1.1
	    endif
         endif
         if(accept_rate.lt.0.3) then
            scalecenter=scalecenter*0.9     ! oder runter
            scalerot=scalerot*0.9
         endif
      endif


      if (iscore.eq.2) then     ! wenn --traj, schreibe akzeptierte schritte raus
        call print_struc2(seed,label,gesa,energies,nlig,
     1  ens,phi,ssi,rot,xa,ya,za,locrests,morph,
     2  nhm,nihm,dlig,has_locrests,lablen,cartstatehandle)
      endif
 4000 continue

c      print*,"scalerot scalecenter " , scalerot, scalecenter
c     Clean up
      call ministate_free_pairlist(ministatehandle)
      end
