c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotpole  --  rotate multipoles to global frame  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotpole" constructs the set of atomic multipoles in the global
c     frame by applying the correct rotation matrix for each site
c
c
      subroutine rotpole
      use mpole
      implicit none
      integer i
      real*8 angle
      real*8 a(3,3)
      logical planar
c
c
c     rotate the atomic multipoles at each site in turn
c
      do i = 1, npole
         call rotmat (i,a,planar,angle)
         call rotsite (i,a,planar,angle)
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotmat  --  find global frame rotation matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotmat" finds the rotation matrix that rotates the local
c     coordinate system into the global frame at a multipole site
c
c
      subroutine rotmat (i,a,planar,ratio)
      use atoms
      use math
      use mpole
      implicit none
      integer i,ii
      integer ix,iy,iz
      real*8 r,dot
      real*8 eps,ratio
      real*8 xi,yi,zi
      real*8 dx,dy,dz
      real*8 dx1,dy1,dz1
      real*8 dx2,dy2,dz2
      real*8 dx3,dy3,dz3
      real*8 a(3,3)
      logical planar
c
c
c     get coordinates and frame definition for the multipole site
c
      ii = ipole(i)
      xi = x(ii)
      yi = y(ii)
      zi = z(ii)
      iz = zaxis(i)
      ix = xaxis(i)
      iy = abs(yaxis(i))
c
c     set values to mark planarity at Z-Bisect and 3-Fold sites
c
      eps = 25.0d0
      planar = .false.
      ratio = 1.0d0
c
c     use the identity matrix as the default rotation matrix
c
      a(1,1) = 1.0d0
      a(2,1) = 0.0d0
      a(3,1) = 0.0d0
      a(1,3) = 0.0d0
      a(2,3) = 0.0d0
      a(3,3) = 1.0d0
c
c     get Z-Only rotation matrix elements for z-axis only
c
      if (polaxe(i) .eq. 'Z-Only') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = 1.0d0
         dy = 0.0d0
         dz = 0.0d0
         dot = a(1,3)
         if (abs(dot) .gt. 0.866d0) then
            dx = 0.0d0
            dy = 1.0d0
            dot = a(2,3)
         end if
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     get Z-then-X rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. 'Z-then-X') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     get Bisector rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. 'Bisector') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dot*a(1,3)
         dy = dy2 - dot*a(2,3)
         dz = dz2 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     get Z-Bisect rotation matrix elements for z- and x-axes;
c     use alternate x-axis if central atom is close to planar
c
      else if (polaxe(i) .eq. 'Z-Bisect') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(iy) - xi
         dy = y(iy) - yi
         dz = z(iy) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx = dx / r
         dy = dy / r
         dz = dz / r
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
c        ratio = (180.0d0-radian*acos(dot)) / eps
c        if (ratio .lt. 1.0d0) then
c           planar = .true.
c           dx = dy1*dz2 - dz1*dy2
c           dy = dz1*dx2 - dx1*dz2
c           dz = dx1*dy2 - dy1*dx2
c           dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
c           if (dot .lt. 0.0d0) then
c              dx = -dx
c              dy = -dy
c              dz = -dz
c              dot = -dot
c           end if
c        end if
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     get 3-Fold rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. '3-Fold') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = x(iy) - xi
         dy = y(iy) - yi
         dz = z(iy) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx3 = dx / r
         dy3 = dy / r
         dz3 = dz / r
         dx = dx1 + dx2 + dx3
         dy = dy1 + dy2 + dy3
         dz = dz1 + dz2 + dz3
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dot*a(1,3)
         dy = dy2 - dot*a(2,3)
         dz = dz2 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
      end if
c
c     finally, find rotation matrix elements for the y-axis
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotsite  --  rotate multipoles at single site  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotsite" rotates the local frame atomic multipoles at a
c     specified site into the global coordinate frame by applying
c     a rotation matrix
c
c
      subroutine rotsite (isite,a,planar,ratio)
      use atoms
      use mpole
      implicit none
      integer i,j,k,m
      integer isite
      real*8 ratio
      real*8 sigmoid
      real*8 spole(13)
      real*8 a(3,3)
      real*8 mp(3,3)
      real*8 rp(3,3)
      logical planar
c
c
c     store and modify multipoles for planar Z-Bisect sites
c
      if (planar) then
         do i = 1, 13
            spole(i) = pole(i,isite)
         end do
c        pole(2,isite) = 0.0d0
c        pole(7,isite) = 0.0d0
c        pole(11,isite) = 0.0d0
c        ratio = sigmoid(3.5d0,ratio)
c        pole(2,isite) = ratio * pole(2,isite)
c        pole(7,isite) = ratio * pole(7,isite)
c        pole(11,isite) = ratio * pole(11,isite)
c        pole(5,isite) = 0.5d0 * (pole(5,isite)+pole(9,isite))
c        pole(9,isite) = pole(5,isite)
      end if
c
c     monopoles have the same value in any coordinate frame
c
      rpole(1,isite) = pole(1,isite)
c
c     rotate the dipoles to the global coordinate frame
c
      do i = 2, 4
         rpole(i,isite) = 0.0d0
         do j = 2, 4
            rpole(i,isite) = rpole(i,isite) + pole(j,isite)*a(i-1,j-1)
         end do
      end do
c
c     rotate the quadrupoles to the global coordinate frame
c
      k = 5
      do i = 1, 3
         do j = 1, 3
            mp(i,j) = pole(k,isite)
            rp(i,j) = 0.0d0
            k = k + 1
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            if (j .lt. i) then
               rp(i,j) = rp(j,i)
            else
               do k = 1, 3
                  do m = 1, 3
                     rp(i,j) = rp(i,j) + a(i,k)*a(j,m)*mp(k,m)
                  end do
               end do
            end if
         end do
      end do
      k = 5
      do i = 1, 3
         do j = 1, 3
            rpole(k,isite) = rp(i,j)
            k = k + 1
         end do
      end do
c
c     restore multipoles to original at planar Z-Bisect sites
c
      if (planar) then
         do i = 1, 13
            pole(i,isite) = spole(i)
         end do
      end if
      return
      end
