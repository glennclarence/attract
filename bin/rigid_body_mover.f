      subroutine rigid_body_mover(
     1 nlig,jl,iori, itra, phi, ssi, rot, xa, ya, za,
     2 fixre, scalerot, scalecenter, dseed)

      implicit none

      include 'max.fin'
      integer nlig, seed

      real*8 phi, ssi, rot, xa, ya, za
      dimension phi(maxlig), ssi(maxlig), rot(maxlig)
      dimension xa(maxlig), ya(maxlig), za(maxlig)

c     Local variables
c     integer dseed,i,ii
      integer i,ii
      integer itra, iori, fixre
      integer ju,ju0,jl,jb, jn, jn0
      real*8 xnull
      real*8 scalecenter,scalerot,rr
      real*8 rotmat,randrot,newrot,sum
      real*8 xaa,pi
      real*8 dseed
      dimension xaa(maxdof)
      dimension rr(maxdof),randrot(0:9),rotmat(0:9),newrot(0:9)
      integer, parameter :: ERROR_UNIT = 0
      pi=3.141592654d0

C$$$ scalerot should be in the range [0 pi], this should be controled
      if ( scalerot<0 .or. scalerot>pi ) then
         write(ERROR_UNIT,*)'scalerot range [0,pi]'
         call exit(1)
      endif

c generate a total of ju random numbers
c     write(*,*)'random rr(1),rr(2)..', rr(1),rr(2),rr(3)
c make an Euler rotation

      if(iori.eq.1) then
        do 190 i=1+fixre,nlig
        ii=3*(i-fixre-1)
c       call crand(dseed,5,rr)
        call GGUBS(dseed,5,rr)

        call random_point_on_unit_sphere(dseed,rr)
        call random_rotation_angle(dseed,scalerot,rr(4))
        call axisrot(rr,randrot)
        call euler2rotmat(phi(i),ssi(i),rot(i),rotmat)
        call matmult(rotmat,randrot,newrot)
        phi(i) = atan2(newrot(5),newrot(2))
        ssi(i) = acos(newrot(8))
        rot(i) = atan2(-newrot(7),-newrot(6))
        if (abs(newrot(8)) >= 0.9999) then
c         gimbal lock
          phi(i) = 0.0d0
          if (abs(newrot(0)) >= 0.9999) then
            if (newrot(0) * newrot(8) < 0) then
              rot(i) = pi
            else
              rot(i) = 0.0d0
            endif
            if (newrot(8) < 0) then
              ssi(i) = pi
            else
              ssi(i) = 0.0d0
            endif
            ssi(i) = 0.0d0
            rot(i) = 0.0d0
          else
            if (newrot(8) < 0) then
              ssi(i) = pi
              rot(i) = -acos(-newrot(0))
            else
              ssi(i) = 0.0d0
              rot(i) = acos(newrot(0))
            endif
          endif
          if (newrot(1) < 0) then
            rot(i) = -rot(i)
          endif
        endif

C$$$  phi,ssi,rot in range of [0 2pi]
        if ( phi(i) < 0 ) then
           phi(i) = phi(i) + 2*pi
        endif
        if ( rot(i) < 0 ) then
           rot(i) = rot(i) + 2*pi
        endif
c       write(*,*)'new ii,c,phi,ssi,rot',i,ii,c,phi(i),ssi(i),rot(i)
  190   continue
      endif

c make a translation of the ligand center
      call GGNML(dseed,3*nlig,rr)
      if(itra.eq.1) then
      do 1220 i=1+fixre,nlig
      ii=jl+3*(i-fixre-1)
      xa(i)=xa(i)+scalecenter*(rr(ii+1))
      ya(i)=ya(i)+scalecenter*(rr(ii+2))
      za(i)=za(i)+scalecenter*(rr(ii+3))
c     write(*,*)'trans-step',i,ii,
c    1 0.5d0-rr(ii+1),0.5d0-rr(ii+2),0.5d0-rr(ii+3)
c    1 rr(ii+1),rr(ii+2),rr(ii+3),xaa(ii+1),xaa(ii+2),xaa(ii+3)
 1220 continue
      endif

      END

