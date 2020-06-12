      SUBROUTINE random_point_on_unit_sphere(DSEED,R)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	random_point_on_unit_sphere genrates uniformly
C     distributed direction on unit sphere from gaussian
C     distributed coordinates
C     date: 02-27-2015
C     author: Zhe Zhang
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      REAL*8             R(3), sum
      DOUBLE PRECISION   DSEED
C                              LOCAL

      call GGNML(DSEED,3,R)
      sum=sqrt(R(1)**2+R(2)**2+R(3)**2)
      do while ( sum .eq. 0)
      	 call GGNML(DSEED,NR,R)
	 sum=sqrt(R(1)**2+R(2)**2+R(3)**2)
      end do

      R(1) = 1d0/sum*R(1)
      R(2) = 1d0/sum*R(2)
      R(3) = 1d0/sum*R(3)

      RETURN
      END

