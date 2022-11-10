c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine exrepel3  --  exch repulsion energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erepel3" calculates the exchange repulsion energy and partitions
c     the energy among the atoms
c
c     literature reference:
c
c     TBD
c
c
      subroutine exrepel3
      use limits
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
c      if (use_mlist) then
c         call exrepel3b
c      else
c         call exrepel3a
c      end if
      call exrepel3a
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine exrepel3a  --  exch repulsion analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "exrepel3a" calculates the exchange repulsion energy and also
c     partitions the energy among the atoms using a double loop
c
c
      subroutine exrepel3a
      use action
      use analyz
      use atoms
      use bound
      use couple
      use energi
      use group
      use inter
      use molcul
      use mpole
      use mutant
      use repel
      use reppot
      use shunt
      use xrepel
      use units
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer ind1,ind2,ind3
      real*8 e,eterm
      real*8 fgrp,taper
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 r,r2
      real*8 normi
      real*8 zxri,zxrk
      real*8 vali,valk
      real*8 dmpi,dmpk
      real*8 cis,cks
      real*8 cix,ckx
      real*8 ciy,cky
      real*8 ciz,ckz
      real*8 rcix,rckx
      real*8 rciy,rcky
      real*8 rciz,rckz
      real*8 S,S2
      real*8 bi(3)
      real*8 bj(3)
      real*8 bk(3)
      real*8 dmpik(5)
      real*8, allocatable :: rscale(:)
      logical proceed,usei
      logical muti,mutk,mutik
      logical header,huge
      character*6 mode
c
c
c     zero out the repulsion energy and partitioning terms
c
      ner = 0
      er = 0.0d0
      do i = 1, n
         aer(i) = 0.0d0
      end do
      if (nrep .eq. 0)  return
c
c     determine pseudo orbital coefficients
c
      call solvcoeff
c
c     rotate the coefficient components into the global frame
c
      call rotcoeff
c
c     perform dynamic allocation of some local arrays
c
      allocate (rscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         rscale(i) = 1.0d0
      end do
c
c     calculate the exchange repulsion interaction energy term
c
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         zxri = zpxr(ii)
         dmpi = dmppxr(ii)
         vali = elepxr(ii)
         cis = rcpxr(1,ii)
         cix = rcpxr(2,ii)
         ciy = rcpxr(3,ii)
         ciz = rcpxr(4,ii)
         usei = use(i)
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = r2scale
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = r3scale
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = r4scale
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = r5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               xk = x(k)
               yk = y(k)
               zk = z(k)
               xr = xk - xi
               yr = yk - yi
               zr = zk - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr * xr + yr * yr + zr * zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  zxrk = zpxr(kk)
                  dmpk = dmppxr(kk)
                  valk = elepxr(kk)
                  cks = rcpxr(1,kk)
                  ckx = rcpxr(2,kk)
                  cky = rcpxr(3,kk)
                  ckz = rcpxr(4,kk)
c
c     choose orthogonal 2-body coordinates / solve rotation matrix
c
                  bk(1) = xr / r
                  bk(2) = yr / r
                  bk(3) = zr / r
                  ind1 = maxloc(abs(bk), dim=1)
                  ind2 = mod(ind1,3) + 1
                  ind3 = mod(ind1+1,3) + 1
                  bi(ind1) = -bk(ind2)
                  bi(ind2) = bk(ind1)
                  bi(ind3) = 0.0d0
                  normi = sqrt(bi(1)**2 + bi(2)**2 + bi(3)**2)
                  bi(1) = bi(1) / normi
                  bi(2) = bi(2) / normi
                  bi(3) = bi(3) / normi
                  bj(1) = bk(2)*bi(3) - bk(3)*bi(2)
                  bj(2) = bk(3)*bi(1) - bk(1)*bi(3)
                  bj(3) = bk(1)*bi(2) - bk(2)*bi(1)
c
c     rotate p orbital cofficients to 2-body (prolate spheroid) frame
c
                  rcix = bi(1)*cix + bi(2)*ciy + bi(3)*ciz
                  rciy = bj(1)*cix + bj(2)*ciy + bj(3)*ciz
                  rciz = bk(1)*cix + bk(2)*ciy + bk(3)*ciz
                  rckx = bi(1)*ckx + bi(2)*cky + bi(3)*ckz
                  rcky = bj(1)*ckx + bj(2)*cky + bj(3)*ckz
                  rckz =-bk(1)*ckx - bk(2)*cky - bk(3)*ckz
c
c     compute overlap matrix
c
                  call computeS(dmpi,dmpk,r,dmpik)
                  S = cis*cks*dmpik(1) + cis*rckz*dmpik(2)
     &             + rciz*cks*dmpik(3) + (rcix*rckx+rciy*rcky)*dmpik(4)
     &             + rciz*rckz*dmpik(5)
                  S2 = S * S
                  e = -hartree*(zxri*valk+zxrk*vali)*S2/r*rscale(k)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall exchange repulsion energy component
c
                  if (e .ne. 0.0d0) then
                     ner = ner + 1
                     er = er + e
                     aer(i) = aer(i) + 0.5d0*e
                     aer(k) = aer(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = 1.0d0
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine solvcoeff  --  solve for orbital coefficients  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "solvcoeff" finds the coefficients for the pseudo orbital
c
c
      subroutine solvcoeff
      use mpole
      use xrepel
      implicit none
      integer ii,k,jj
      integer ind1,ind2,ind3
      real*8 cr,cs
      real*8 pcoeff(3)
      real*8 ppole(3)
      real*8 p2p1,p3p1
      logical l1,l2,l3
c
c
c     determine pseudo orbital coefficients
c
      do ii = 1, npole
         do k = 1, 3
            pcoeff(k) = 0.0d0
         end do
         ppole(1) = pole(2,ii)
         ppole(2) = pole(3,ii)
         ppole(3) = pole(4,ii)
         cr = crpxr(ii)
         l1 = (abs(ppole(1)) < 1.0d-10)
         l2 = (abs(ppole(2)) < 1.0d-10)
         l3 = (abs(ppole(3)) < 1.0d-10)
c
c     case for no dipole
c
         if (l1.and.l2.and.l3) then
            cs = 1.0d0
            ind1 = 1
            ind2 = 2
            ind3 = 3
c
c     case for p orbital coefficients set to 0
c
         else if (cr < 1.0d-10) then
            cs = 1.0d0
            ind1 = 1
            ind2 = 2
            ind3 = 3
c
c     determine normalized coefficients
c
         else
            cs = 1.0d0 / sqrt(1.0d0 + cr)
c
c     determine index for largest absolute dipole component
c
            ind1 = maxloc(abs(ppole), dim=1)
            ind2 = mod(ind1,3) + 1
            ind3 = mod(ind1+1,3) + 1
            p2p1 = ppole(ind2) / ppole(ind1)
            p3p1 = ppole(ind3) / ppole(ind1)
            pcoeff(ind1) = cs * sqrt(cr / (1.0d0 + p2p1**2 + p3p1**2))
            if (ppole(ind1) < 0.0d0) then
               pcoeff(ind1) = -pcoeff(ind1)
            end if
            pcoeff(ind2) = pcoeff(ind1) * p2p1
            pcoeff(ind3) = pcoeff(ind1) * p3p1
         end if
         cpxr(1,ii) = cs
         cpxr(ind1+1,ii) = pcoeff(ind1)
         cpxr(ind2+1,ii) = pcoeff(ind2)
         cpxr(ind3+1,ii) = pcoeff(ind3)
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine rotcoeff  --  rotate orbital coefficients  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "rotcoeff" rotates the coefficients for the pseudo orbital
c     to global coordinates
c
c
      subroutine rotcoeff
      use mpole
      use xrepel
      implicit none
      integer isite
      integer i,j,k
      real*8 a(3,3)
      real*8 cpole(4)
      logical planar
c
c
c     rotate pseudo orbital coefficients
c
      do isite = 1, npole
c
c     determine rotation matrix
c
         call rotmat (isite,a,planar)
c
c     copy local frame coefficients
c
         do i = 1, 4
            cpole(i) = cpxr(i,isite)
         end do
c
c     s orbital coefficients have the same value in any coordinate frame
c
         rcpxr(1,isite) = cpole(1)
c
c     rotate the p orbital coefficients to the global coordinate frame
c
         do i = 2, 4
            rcpxr(i,isite) = 0.0d0
            do j = 2, 4
               rcpxr(i,isite) = rcpxr(i,isite) + cpole(j)*a(i-1,j-1)
            end do
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine rotcoeff  --  rotate orbital coefficients  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "rotcoeff" rotates the coefficients for the pseudo orbital
c     to global coordinates
c
c
      subroutine computeS (dmpi,dmpk,r,dmpik)
      implicit none
      real*8 dmpi,dmpk
      real*8 r
      real*8 eps,diff
      real*8 expi,expk
      real*8 dmpidmpk
      real*8 dampi,dampi2
      real*8 dampi3,dampi4
      real*8 dampk,dampk2
      real*8 dampk3,dampk4
      real*8 dampik
      real*8 dampik2,dampik3
      real*8 tau,tau2
      real*8 taup,taum
      real*8 taupm
      real*8 kappa
      real*8 kappap,kappam
      real*8 kappap2,kappam2
      real*8 pre
      real*8 term1,term2
      real*8 ss,spz,pzs
      real*8 pxpx,pzpz
      real*8 dmpik(5)
c
c
c     compute tolerance value for damping exponents
c
      eps = 0.001d0
      diff = abs(dmpi-dmpk)
c
c     treat the case where alpha damping exponents are equal
c
      if (diff .lt. eps) then
         dampi = dmpi * r
         dampi2 = dampi * dampi
         dampi3 = dampi2 * dampi
         dampi4 = dampi3 * dampi
         expi = exp(-dampi)
         ss = (1.0d0 + dampi + dampi2 / 3.0d0) * expi
         spz = 0.5d0 * dampi * (1.0d0 + dampi + dampi2 / 3.0d0) * expi
         pzs = spz
         pxpx = (1.0d0 + dampi + 2.0d0/5.0d0 * dampi2
     &            + dampi3 / 15.0d0) * expi
         pzpz = (-1.0d0 - dampi - dampi2 / 5.0d0
     &            + 2.0d0 / 15.0d0 * dampi3 + dampi4 / 15.0d0) * expi
c
c     treat the case where alpha damping exponents are unequal
c
      else
         dmpidmpk = dmpi + dmpk
         dampik = 0.5d0 * (dmpidmpk) * r
         tau = (dmpi - dmpk) / dmpidmpk
         kappa = 0.5d0 * (tau + 1.0d0 / tau)
         kappap = 1.0d0 + kappa
         kappam = 1.0d0 - kappa
         dampi = dmpi * r
         dampk = dmpk * r
         expi = exp(-dampi)
         expk = exp(-dampk)
         dampik2 = dampik * dampik
         dampik3 = dampik2 * dampik
         tau2 = tau * tau
         taup = sqrt(1.0d0 + tau)
         taum = sqrt(1.0d0 - tau)
         taupm = taup * taum
         kappap2 = kappap * kappap
         kappam2 = kappam * kappam
         dampi2 = dampi * dampi
         dampi3 = dampi2 * dampi
         dampi4 = dampi3 * dampi
         dampk2 = dampk * dampk
         dampk3 = dampk2 * dampk
         dampk4 = dampk3 * dampk
         pre = taupm / (tau * dampik)
         term1 =-kappam * (2.0d0 * kappap + dampi) * expi
         term2 = kappap * (2.0d0 * kappam + dampk) * expk
         ss = pre * (term1 + term2)
         pre = taup / (taum * tau * dampik2)
         term1 =-kappam2 * (6.0d0 * kappap * (1.0d0 + dampi)
     &                      + 2.0d0 * dampi2) * expi
         term2 = kappap * (6.0d0 * kappam2 * (1.0d0 + dampk)
     &                      + 4.0d0 * kappam * dampk2 + dampk3) * expk
         spz = pre * (term1 + term2)

         pre = -taum / (taup * tau * dampik2)
         term1 =-kappap2 * (6.0d0 * kappam * (1.0d0 + dampk)
     &                      + 2.0d0 * dampk2) * expk
         term2 = kappam * (6.0d0 * kappap2 * (1.0d0 + dampi)
     &                      + 4.0d0 * kappap * dampi2 + dampi3) * expi
         pzs = pre * (term1 + term2)
         pre = 1.0d0 / (taupm * tau * dampik3)
         term1 =-kappam2 * (24.0d0 * kappap2 * (1.0d0 + dampi)
     &               + 12.0d0 * kappap * dampi2 + 2.0d0 * dampi3) * expi
         term2 = kappap2 * (24.0d0 * kappam2 * (1.0d0 + dampk)
     &               + 12.0d0 * kappam * dampk2 + 2.0d0 * dampk3) * expk
         pxpx = pre * (term1 + term2)
         term1 =-kappam2 * (48.0d0 * kappap2 * (1.0d0 + dampi
     &               + 0.5d0 * dampi2) + 2.0d0 * (5.0d0 + 6.0d0 * kappa)
     &                      * dampi3 + 2.0d0 * dampi4) * expi
         term2 = kappap2 * (48.0d0 * kappam2 * (1.0d0 + dampk
     &               + 0.5d0 * dampk2) + 2.0d0 * (5.0d0 - 6.0d0 * kappa)
     &                      * dampk3 + 2.0d0 * dampk4) * expk
         pzpz = pre * (term1 + term2)
      end if
      dmpik(1) = ss
      dmpik(2) = spz
      dmpik(3) = pzs
      dmpik(4) = pxpx
      dmpik(5) = pzpz
      return
      end
