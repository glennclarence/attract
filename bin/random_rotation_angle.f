      SUBROUTINE random_rotation_angle(DSEED,MAG,ANGLE)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     random_rotation_angle genrates random rotation angles in the range
C     of [0 pi] from gamma-distribution-like distribution which is calculated
C     from euler rotation matrix
C     date: 02-27-2015
C     author: Zhe Zhang
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8    EULER,R,AXIS,ANGLE,MAG
      dimension R(0:9),EULER(3),AXIS(3)
      integer   DSEED
      integer            i

C                              LOCAL
      call GGNML(DSEED,3,euler)
      euler(1) = Euler(1)*mag
      euler(2) = euler(2)*mag
      euler(3) = euler(3)*mag
      call euler2rotmat(euler(1),euler(2),euler(3),R)
      call R_to_axis_angle(R,AXIS,ANGLE)

      END
