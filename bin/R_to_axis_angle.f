      subroutine R_to_axis_angle(rotmat,axis,theta)

C$$$    """Convert the rotation matrix into the axis-angle notation.

C$$$    Conversion equations
C$$$    ====================

C$$$    From Wikipedia (Http://en.wikipedia.org/wiki/Rotation_matrix), the conversion is given by::

C$$$        x = Qzy-Qyz
C$$$        y = Qxz-Qzx
C$$$        z = Qyx-Qxy
C$$$        r = hypot(x,hypot(y,z))
C$$$        t = Qxx+Qyy+Qzz
C$$$        theta = atan2(r,t-1)

C$$$    @param matrix:  The 3x3 rotation matrix to update.
C$$$    @type matrix:   3x3 numpy array
C$$$    @return:    The 3D rotation axis and angle.
C$$$    @rtype:     numpy 3D rank-1 array, float
C$$$    """

      real*8 rotmat,axis,theta,sum
      dimension rotmat(0:9),axis(3)

      axis(1) = rotmat(7)-rotmat(5)
      axis(2) = rotmat(2)-rotmat(6)
      axis(3) = rotmat(3)-rotmat(1)

      sum=1d0/sqrt(axis(1)**2+axis(2)**2+axis(3)**2)
      axis(1) = sum*axis(1)
      axis(2) = sum*axis(2)
      axis(3) = sum*axis(3)

      theta = atan2(1d0/sum,rotmat(0)+rotmat(4)+rotmat(8)-1)


      return
      end


