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
c     "exrepel3" calculates the exchange repulsion energy and partitions
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
      use atomid
      use atoms
      use bound
      use cell
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mpole
      use mutant
      use repel
      use reppot
      use shunt
      use units
      use usage
      use xrepel
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
      real*8 vali,valk,valik
      real*8 dmpi,dmpk
      real*8 dis,dks
      real*8 dmpip,dmpkp
      real*8 cis,cks
      real*8 cix,ckx
      real*8 ciy,cky
      real*8 ciz,ckz
      real*8 rcix,rckx
      real*8 rciy,rcky
      real*8 rciz,rckz
      real*8 intS,intS2
      real*8 intK,intJ
      real*8 intaAb,intbBa
      real*8 intbAb,intaBa
      real*8 overlapTotal
      real*8 NEaAbTotal,NEaBbTotal,NEbAbTotal,NEaBaTotal
      real*8 coulombTotal,exchangeTotal
      real*8 pre,termS0
      real*8 termS1,termS2
      real*8 bi(3)
      real*8 bj(3)
      real*8 bk(3)
      real*8 coeffi(4),coeffk(4)
      real*8, allocatable :: rscale(:)
      logical proceed,usei
      logical muti,mutk,mutik
      logical header,huge
      logical exact
      character*6 mode
      character*8 example
      real*8 overlapSS
c
c
c
c     for non-exchange integral terms, choose exact or STO
c
      exact = .true.
c     zero out the repulsion energy and partitioning terms
c
      ner = 0
      er = 0.0d0
      do i = 1, n
         aer(i) = 0.0d0
      end do
      if (nrep .eq. 0)  return
c
c     test functions, all evaluated at STO-4G
c
c      call testcoulomb()
c      call testcoulombflip()
c      call testcoulombsto()
c      call testexchange()
c      call teststooverlap()
c      call testoverlap()
c      call testoverlapflip()
c      call teststoNE()
c      call testNE()
c      call testNEflip()
c      call testoverlapsum()
c      call testNEsumbAb()
c      call testcoulombsum()
c      call testexchangesum()
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
         dis = drpxr(ii)
         dmpip = dis * dmpi
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
                  dks = drpxr(kk)
                  dmpkp = dks * dmpk
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
                  rckz = bk(1)*ckx + bk(2)*cky + bk(3)*ckz
                  coeffi = (/ cis, rcix, rciy, rciz /)
                  coeffk = (/ cks, rckx, rcky, rckz /)
                  valik = vali * valk
                  if (xreptyp == "S2R") then
                     intS = overlapTotal (coeffi, coeffk, dmpi, dmpk,
     &                                    dmpip, dmpkp, 0.0d0, r, exact)
                     intS2 = intS * intS
                     e = -hartree*(zxri*valk+zxrk*vali)*intS2/r
     &                                                        *rscale(k)
                  else if (xreptyp == "H2") then
                     intS = overlapTotal (coeffi, coeffk, dmpi, dmpk,
     &                                    dmpip, dmpkp, 0.0d0, r, exact)
                     intS2 = intS**2
                     intJ = coulombTotal (coeffi, coeffk, dmpi, dmpk,
     &                                                  0.0d0, r, exact)
                     intaAb = NEaAbTotal (coeffi, coeffk, dmpi, dmpk,
     &                                                  0.0d0, r, exact)
                     intbBa = NEaBbTotal (coeffi, coeffk, dmpi, dmpk,
     &                                                  0.0d0, r, exact)
                     intbAb = NEbAbTotal (coeffk, dmpk, 0.0d0, r, exact)
                     intaBa = NEaBaTotal (coeffi, dmpi, 0.0d0, r, exact)
                     intK = exchangeTotal (coeffi, coeffk, dmpi, dmpk,
     &                                                         0.0d0, r)
                     pre = hartree / (1.0d0 - intS2)
                     termS0 = -valik * intK
                     termS1 = -intS * (zxri * valk * intaAb
     &                               + zxrk * vali * intbBa)
                     termS2 = intS2 * (zxri * valk * intbAb
     &                               + zxrk * vali * intaBa
     &                               + valik * intJ)
                     e = pre * (termS0 + termS1 + termS2)
                  end if
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
c
c     perform deallocation of some local arrays
c
      deallocate (rscale)
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
      subroutine testcoulomb()
      implicit none
      real*8 a,b,z1,z2
      real*8 coulombSSSS, cSSSS
      real*8 coulombSSPzS, cSSPzS
      real*8 coulombSSSPz, cSSSPz
      real*8 coulombSSPzPz, cSSPzPz
      real*8 coulombSPzSS, cSPzSS
      real*8 coulombSPzSPz, cSPzSPz
      real*8 coulombSPzPzPz, cSPzPzPz
      real*8 coulombPzPzSS, cPzPzSS
      real*8 coulombPzPzSPz, cPzPzSPz
      real*8 coulombPzPzPzPz, cPzPzPzPz
      real*8 coulombSSPxPx, cSSPxPx
      real*8 coulombSPxSPx, cSPxSPx
      real*8 coulombPxPxSS, cPxPxSS
      real*8 coulombPxPxPxPx, cPxPxPxPx
      real*8 coulombSPxPxPz, cSPxPxPz
      real*8 coulombSPzPxPx, cSPzPxPx
      real*8 coulombPxPxSPz, cPxPxSPz
      real*8 coulombPxPxPzPz, cPxPxPzPz
      real*8 coulombPxPzSPx, cPxPzSPx
      real*8 coulombPxPzPxPz, cPxPzPxPz
      real*8 coulombPzPzPxPx, cPzPzPxPx
      real*8 coulombPxPxPyPy, cPxPxPyPy
      real*8 coulombPxPyPxPy, cPxPyPxPy
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      print*, "Coulomb"
c
c     coulombSSSS
c
      cSSSS = coulombSSSS(a,b,z1,z2) - (0.49988666809518023d0)
      print*, "cSSSS a != b", cSSSS
      cSSSS = coulombSSSS(a,a,z1,z2) - (0.49995653530091183d0)
      print*, "cSSSS a = b", cSSSS
      cSSSS = coulombSSSS(b,a,z1,z2) - (0.49988666809518023d0)
      print*, "cSSSS b != a", cSSSS
      print*, ""
c
c     coulombSSSPz
c
      cSSSPz = coulombSSSPz(a,b,z1,z2) - (-0.0828085932448284d0)
      print*, "cSSSPz a != b", cSSSPz
      cSSSPz = coulombSSSPz(a,a,z1,z2) - (-0.07124826952203973d0)
      print*, "cSSSPz a = b", cSSSPz
      cSSSPz = coulombSSSPz(b,a,z1,z2) - (-0.07108221554866272d0)
      print*, "cSSSPz b != a", cSSSPz
      print*, ""
c
c     coulombSSPzPz
c
      cSSPzPz = coulombSSPzPz(a,b,z1,z2) - (0.5388382815880638d0)
      print*, "cSSPzPz a != b", cSSPzPz
      cSSPzPz = coulombSSPzPz(a,a,z1,z2) - (0.5296597904077978d0)
      print*, "cSSPzPz a = b", cSSPzPz
      cSSPzPz = coulombSSPzPz(b,a,z1,z2) - (0.5290770985691252d0)
      print*, "cSSPzPz b != a", cSSPzPz
      print*, ""
c
c     coulombSPzSS
c
      cSPzSS = coulombSPzSS(a,b,z1,z2) - (0.07108221554866272d0)
      print*, "cSPzSS a != b", cSPzSS
      cSPzSS = coulombSPzSS(a,a,z1,z2) - (0.07124826952203973d0)
      print*, "cSPzSS a = b", cSPzSS
      cSPzSS = coulombSPzSS(b,a,z1,z2) - (0.0828085932448284d0)
      print*, "cSPzSS b != a", cSPzSS
      print*, ""
c
c     coulombSPzSPz
c
      cSPzSPz = coulombSPzSPz(a,b,z1,z2) - (-0.022583362072077230d0)
      print*, "cSPzSPz a != b", cSPzSPz
      cSPzSPz = coulombSPzSPz(a,a,z1,z2) - (-0.0198440568868558450d0)
      print*, "cSPzSPz a = b", cSPzSPz
      cSPzSPz = coulombSPzSPz(b,a,z1,z2) - (-0.022583362072077230d0)
      print*, "cSPzSPz b != a", cSPzSPz
      print*, ""
c
c     coulombSPzPzPz
c
      cSPzPzPz = coulombSPzPzPz(a,b,z1,z2) - (0.08387741384927963d0)
      print*, "cSPzPzPz a != b", cSPzPzPz
      cSPzPzPz = coulombSPzPzPz(a,a,z1,z2) - (0.0821342383420232d0)
      print*, "cSPzPzPz a = b", cSPzPzPz
      cSPzPzPz = coulombSPzPzPz(b,a,z1,z2) - (0.09418948144354163d0)
      print*, "cSPzPzPz b != a", cSPzPzPz
      print*, ""
c
c     coulombPzPzSS
c
      cPzPzSS = coulombPzPzSS(a,b,z1,z2) - (0.5290770985691252d0)
      print*, "cPzPzSS a != b", cPzPzSS
      cPzPzSS = coulombPzPzSS(a,a,z1,z2) - (0.5296597904077978d0)
      print*, "cPzPzSS a = b", cPzPzSS
      cPzPzSS = coulombPzPzSS(b,a,z1,z2) - (0.5388382815880638d0)
      print*, "cPzPzSS b != a", cPzPzSS
      print*, ""
c
c     coulombPzPzSPz
c
      cPzPzSPz = coulombPzPzSPz(a,b,z1,z2) - (-0.09418948144354163d0)
      print*, "cPzPzSPz a != b", cPzPzSPz
      cPzPzSPz = coulombPzPzSPz(a,a,z1,z2) - (-0.0821342383420232d0)
      print*, "cPzPzSPz a = b", cPzPzSPz
      cPzPzSPz = coulombPzPzSPz(b,a,z1,z2) - (-0.08387741384927963d0)
      print*, "cPzPzSPz b != a", cPzPzSPz
      print*, ""
c
c     coulombPzPzPzPz
c
      cPzPzPzPz = coulombPzPzPzPz(a,b,z1,z2) - (0.5703689519863779d0)
      print*, "cPzPzPzPz a != b", cPzPzPzPz
      cPzPzPzPz = coulombPzPzPzPz(a,a,z1,z2) - (0.5635236359648949d0)
      print*, "cPzPzPzPz a = b", cPzPzPzPz
      cPzPzPzPz = coulombPzPzPzPz(b,a,z1,z2) - (0.5703689519863779d0)
      print*, "cPzPzPzPz b != a", cPzPzPzPz
      print*, ""
c
c     coulombSSPxPx
c
      cSSPxPx = coulombSSPxPx(a,b,z1,z2) - (0.47912883137076784d0)
      print*, "cSSPxPx a != b", cSSPxPx
      cSSPxPx = coulombSSPxPx(a,a,z1,z2) - (0.48466391688540134d0)
      print*, "cSSPxPx a = b", cSSPxPx
      cSSPxPx = coulombSSPxPx(b,a,z1,z2) - (0.4845990763501186d0)
      print*, "cSSPxPx b != a", cSSPxPx
      print*, ""
c
c     coulombSPxSPx
c
      cSPxSPx = coulombSPxSPx(a,b,z1,z2) - (0.01178206610937095d0)
      print*, "cSPxSPx a != b", cSPxSPx
      cSPxSPx = coulombSPxSPx(a,a,z1,z2) - (0.01015319212807509d0)
      print*, "cSPxSPx a = b", cSPxSPx
      cSPxSPx = coulombSPxSPx(b,a,z1,z2) - (0.01178206610937095d0)
      print*, "cSPxSPx b != a", cSPxSPx
      print*, ""
c
c     coulombPxPxSS
c
      cPxPxSS = coulombPxPxSS(a,b,z1,z2) - (0.4845990763501186d0)
      print*, "cPxPxSS a != b", cPxPxSS
      cPxPxSS = coulombPxPxSS(a,a,z1,z2) - (0.48466391688540134d0)
      print*, "cPxPxSS a = b", cPxPxSS
      cPxPxSS = coulombPxPxSS(b,a,z1,z2) - (0.47912883137076784d0)
      print*, "cPxPxSS b != a", cPxPxSS
      print*, ""
c
c     coulombPxPxPxPx
c
      cPxPxPxPx = coulombPxPxPxPx(a,b,z1,z2) - (0.4692241032238778d0)
      print*, "cPxPxPxPx a != b", cPxPxPxPx
      cPxPxPxPx = coulombPxPxPxPx(a,a,z1,z2) - (0.47343555297830386d0)
      print*, "cPxPxPxPx a = b", cPxPxPxPx
      cPxPxPxPx = coulombPxPxPxPx(b,a,z1,z2) - (0.4692241032238778d0)
      print*, "cPxPxPxPx b != a", cPxPxPxPx
      print*, ""
c
c     coulombSPxPxPz
c
      cSPxPxPz = coulombSPxPxPz(a,b,z1,z2) - (-0.008359792145233679d0)
      print*, "cSPxPxPz a != b", cSPxPxPz
      cSPxPxPz = coulombSPxPxPz(a,a,z1,z2) - (-0.0063400196190138685d0)
      print*, "cSPxPxPz a = b", cSPxPxPz
      cSPxPxPz = coulombSPxPxPz(b,a,z1,z2) - (-0.007228374045474331d0)
      print*, "cSPxPxPz b != a", cSPxPxPz
      print*, ""
c
c     coulombSPzPxPx
c
      cSPzPxPx = coulombSPzPxPx(a,b,z1,z2) - (0.062332604510992d0)
      print*, "cSPzPxPx a != b", cSPzPxPx
      cSPzPxPx = coulombSPzPxPx(a,a,z1,z2) - (0.06473232527751177d0)
      print*, "cSPzPxPx a = b", cSPzPxPx
      cSPzPxPx = coulombSPzPxPx(b,a,z1,z2) - (0.07524608872012231d0)
      print*, "cSPzPxPx b != a", cSPzPxPx
      print*, ""
c
c     coulombPxPxSPz
c
      cPxPxSPz = coulombPxPxSPz(a,b,z1,z2) - (-0.07524608872012231d0)
      print*, "cPxPxSPz a != b", cPxPxSPz
      cPxPxSPz = coulombPxPxSPz(a,a,z1,z2) - (-0.06473232527751177d0)
      print*, "cPxPxSPz a = b", cPxPxSPz
      cPxPxSPz = coulombPxPxSPz(b,a,z1,z2) - (-0.062332604510992d0)
      print*, "cPxPxSPz b != a", cPxPxSPz
      print*, ""
c
c     coulombPxPxPzPz
c
      cPxPxPzPz = coulombPxPxPzPz(a,b,z1,z2) - (0.5164642045478991d0)
      print*, "cPxPxPzPz a != b", cPxPxPzPz
      cPxPxPzPz = coulombPxPxPzPz(a,a,z1,z2) - (0.5089658284234035d0)
      print*, "cPxPxPzPz a = b", cPxPxPzPz
      cPxPxPzPz = coulombPxPxPzPz(b,a,z1,z2) - (0.5013159038041091d0)
      print*, "cPxPxPzPz b != a", cPxPxPzPz
      print*, ""
c
c     coulombPxPzSPx
c
      cPxPzSPx = coulombPxPzSPx(a,b,z1,z2) - (0.0072283740454743310d0)
      print*, "cPxPzSPx a != b", cPxPzSPx
      cPxPzSPx = coulombPxPzSPx(a,a,z1,z2) - (0.00634001961901386850d0)
      print*, "cPxPzSPx a = b", cPxPzSPx
      cPxPzSPx = coulombPxPzSPx(b,a,z1,z2) - (0.0083597921452336790d0)
      print*, "cPxPzSPx b != a", cPxPzSPx
      print*, ""
c
c     coulombPxPzPxPz
c
      cPxPzPxPz = coulombPxPzPxPz(a,b,z1,z2) - (-0.006034801127456413d0)
      print*, "cPxPzPxPz a != b", cPxPzPxPz
      cPxPzPxPz = coulombPxPzPxPz(a,a,z1,z2) - (-0.004848437602982191d0)
      print*, "cPxPzPxPz a = b", cPxPzPxPz
      cPxPzPxPz = coulombPxPzPxPz(b,a,z1,z2) - (-0.006034801127456413d0)
      print*, "cPxPzPxPz b != a", cPxPzPxPz
      print*, ""
c
c     coulombPzPzPxPx
c
      cPzPzPxPx = coulombPzPzPxPx(a,b,z1,z2) - (0.5013159038041091d0)
      print*, "cPzPzPxPx a != b", cPzPzPxPx
      cPzPzPxPx = coulombPzPzPxPx(a,a,z1,z2) - (0.5089658284234035d0)
      print*, "cPzPzPxPx a = b", cPzPzPxPx
      cPzPzPxPx = coulombPzPzPxPx(b,a,z1,z2) - (0.5164642045478991d0)
      print*, "cPzPzPxPx b != a", cPzPzPxPx
      print*, ""
c
c     coulombPxPxPyPy
c
      cPxPxPyPy = coulombPxPxPyPy(a,b,z1,z2) - (0.46572834151377246d0)
      print*, "cPxPxPyPy a != b", cPxPxPyPy
      cPxPxPyPy = coulombPxPxPyPy(a,a,z1,z2) - (0.47076793678829376d0)
      print*, "cPxPxPyPy a = b", cPxPxPyPy
      cPxPxPyPy = coulombPxPxPyPy(b,a,z1,z2) - (0.4657283415137724d0)
      print*, "cPxPxPyPy b != a", cPxPxPyPy
      print*, ""
c
c     coulombPxPyPxPy
c
      cPxPyPxPy = coulombPxPyPxPy(a,b,z1,z2) - (0.0017478808550526867d0)
      print*, "cPxPyPxPy a != b", cPxPyPxPy
      cPxPyPxPy = coulombPxPyPxPy(a,a,z1,z2) - (0.001333808095005061d0)
      print*, "cPxPyPxPy a = b", cPxPyPxPy
      cPxPyPxPy = coulombPxPyPxPy(b,a,z1,z2) - (0.0017478808550526867d0)
      print*, "cPxPyPxPy b != a", cPxPyPxPy
      print*, ""
      return
      end
c
c
      subroutine testcoulombflip()
      implicit none
      real*8 a,b,z1,z2
      real*8 coulombSSSS, cSSSS
      real*8 coulombSSPzS, cSSPzS
      real*8 coulombSSSPz, cSSSPz
      real*8 coulombSSPzPz, cSSPzPz
      real*8 coulombSPzSS, cSPzSS
      real*8 coulombSPzSPz, cSPzSPz
      real*8 coulombSPzPzPz, cSPzPzPz
      real*8 coulombPzPzSS, cPzPzSS
      real*8 coulombPzPzSPz, cPzPzSPz
      real*8 coulombPzPzPzPz, cPzPzPzPz
      real*8 coulombSSPxPx, cSSPxPx
      real*8 coulombSPxSPx, cSPxSPx
      real*8 coulombPxPxSS, cPxPxSS
      real*8 coulombPxPxPxPx, cPxPxPxPx
      real*8 coulombSPxPxPz, cSPxPxPz
      real*8 coulombSPzPxPx, cSPzPxPx
      real*8 coulombPxPxSPz, cPxPxSPz
      real*8 coulombPxPxPzPz, cPxPxPzPz
      real*8 coulombPxPzSPx, cPxPzSPx
      real*8 coulombPxPzPxPz, cPxPzPxPz
      real*8 coulombPzPzPxPx, cPzPzPxPx
      real*8 coulombPxPxPyPy, cPxPxPyPy
      real*8 coulombPxPyPxPy, cPxPyPxPy
c
c
      a = 3.5
      b = 3.0
      z1 = 1.5
      z2 = -0.5
      print*, "Coulomb flip z1 and z2"
c
c     coulombSSSS
c
      cSSSS = coulombSSSS(a,b,z1,z2) - (0.49988666809518023d0)
      print*, "cSSSS a != b", cSSSS
      cSSSS = coulombSSSS(a,a,z1,z2) - (0.49995653530091183d0)
      print*, "cSSSS a = b", cSSSS
      cSSSS = coulombSSSS(b,a,z1,z2) - (0.49988666809518023d0)
      print*, "cSSSS b != a", cSSSS
      print*, ""
c
c     coulombSSSPz
c
      cSSSPz = coulombSSSPz(a,b,z1,z2) - (0.0828085932448284d0)
      print*, "cSSSPz a != b", cSSSPz
      cSSSPz = coulombSSSPz(a,a,z1,z2) - (0.07124826952203973d0)
      print*, "cSSSPz a = b", cSSSPz
      cSSSPz = coulombSSSPz(b,a,z1,z2) - (0.07108221554866272d0)
      print*, "cSSSPz b != a", cSSSPz
      print*, ""
c
c     coulombSSPzPz
c
      cSSPzPz = coulombSSPzPz(a,b,z1,z2) - (0.5388382815880638d0)
      print*, "cSSPzPz a != b", cSSPzPz
      cSSPzPz = coulombSSPzPz(a,a,z1,z2) - (0.5296597904077978d0)
      print*, "cSSPzPz a = b", cSSPzPz
      cSSPzPz = coulombSSPzPz(b,a,z1,z2) - (0.5290770985691252d0)
      print*, "cSSPzPz b != a", cSSPzPz
      print*, ""
c
c     coulombSPzSS
c
      cSPzSS = coulombSPzSS(a,b,z1,z2) - (-0.07108221554866272d0)
      print*, "cSPzSS a != b", cSPzSS
      cSPzSS = coulombSPzSS(a,a,z1,z2) - (-0.07124826952203973d0)
      print*, "cSPzSS a = b", cSPzSS
      cSPzSS = coulombSPzSS(b,a,z1,z2) - (-0.0828085932448284d0)
      print*, "cSPzSS b != a", cSPzSS
      print*, ""
c
c     coulombSPzSPz
c
      cSPzSPz = coulombSPzSPz(a,b,z1,z2) - (-0.02258336207207723d0)
      print*, "cSPzSPz a != b", cSPzSPz
      cSPzSPz = coulombSPzSPz(a,a,z1,z2) - (-0.019844056886855845d0)
      print*, "cSPzSPz a = b", cSPzSPz
      cSPzSPz = coulombSPzSPz(b,a,z1,z2) - (-0.02258336207207723d0)
      print*, "cSPzSPz b != a", cSPzSPz
      print*, ""
c
c     coulombSPzPzPz
c
      cSPzPzPz = coulombSPzPzPz(a,b,z1,z2) - (-0.08387741384927963d0)
      print*, "cSPzPzPz a != b", cSPzPzPz
      cSPzPzPz = coulombSPzPzPz(a,a,z1,z2) - (-0.0821342383420232d0)
      print*, "cSPzPzPz a = b", cSPzPzPz
      cSPzPzPz = coulombSPzPzPz(b,a,z1,z2) - (-0.09418948144354163d0)
      print*, "cSPzPzPz b != a", cSPzPzPz
      print*, ""
c
c     coulombPzPzSS
c
      cPzPzSS = coulombPzPzSS(a,b,z1,z2) - (0.5290770985691252d0)
      print*, "cPzPzSS a != b", cPzPzSS
      cPzPzSS = coulombPzPzSS(a,a,z1,z2) - (0.5296597904077978d0)
      print*, "cPzPzSS a = b", cPzPzSS
      cPzPzSS = coulombPzPzSS(b,a,z1,z2) - (0.5388382815880638d0)
      print*, "cPzPzSS b != a", cPzPzSS
      print*, ""
c
c     coulombPzPzSPz
c
      cPzPzSPz = coulombPzPzSPz(a,b,z1,z2) - (0.09418948144354163d0)
      print*, "cPzPzSPz a != b", cPzPzSPz
      cPzPzSPz = coulombPzPzSPz(a,a,z1,z2) - (0.0821342383420232d0)
      print*, "cPzPzSPz a = b", cPzPzSPz
      cPzPzSPz = coulombPzPzSPz(b,a,z1,z2) - (0.08387741384927963d0)
      print*, "cPzPzSPz b != a", cPzPzSPz
      print*, ""
c
c     coulombPzPzPzPz
c
      cPzPzPzPz = coulombPzPzPzPz(a,b,z1,z2) - (0.5703689519863779d0)
      print*, "cPzPzPzPz a != b", cPzPzPzPz
      cPzPzPzPz = coulombPzPzPzPz(a,a,z1,z2) - (0.5635236359648949d0)
      print*, "cPzPzPzPz a = b", cPzPzPzPz
      cPzPzPzPz = coulombPzPzPzPz(b,a,z1,z2) - (0.5703689519863779d0)
      print*, "cPzPzPzPz b != a", cPzPzPzPz
      print*, ""
c
c     coulombSSPxPx
c
      cSSPxPx = coulombSSPxPx(a,b,z1,z2) - (0.47912883137076784d0)
      print*, "cSSPxPx a != b", cSSPxPx
      cSSPxPx = coulombSSPxPx(a,a,z1,z2) - (0.48466391688540134d0)
      print*, "cSSPxPx a = b", cSSPxPx
      cSSPxPx = coulombSSPxPx(b,a,z1,z2) - (0.4845990763501186d0)
      print*, "cSSPxPx b != a", cSSPxPx
      print*, ""
c
c     coulombSPxSPx
c
      cSPxSPx = coulombSPxSPx(a,b,z1,z2) - (0.01178206610937095d0)
      print*, "cSPxSPx a != b", cSPxSPx
      cSPxSPx = coulombSPxSPx(a,a,z1,z2) - (0.01015319212807509d0)
      print*, "cSPxSPx a = b", cSPxSPx
      cSPxSPx = coulombSPxSPx(b,a,z1,z2) - (0.01178206610937095d0)
      print*, "cSPxSPx b != a", cSPxSPx
      print*, ""
c
c     coulombPxPxSS
c
      cPxPxSS = coulombPxPxSS(a,b,z1,z2) - (0.4845990763501186d0)
      print*, "cPxPxSS a != b", cPxPxSS
      cPxPxSS = coulombPxPxSS(a,a,z1,z2) - (0.48466391688540134d0)
      print*, "cPxPxSS a = b", cPxPxSS
      cPxPxSS = coulombPxPxSS(b,a,z1,z2) - (0.47912883137076784d0)
      print*, "cPxPxSS b != a", cPxPxSS
      print*, ""
c
c     coulombPxPxPxPx
c
      cPxPxPxPx = coulombPxPxPxPx(a,b,z1,z2) - (0.4692241032238778d0)
      print*, "cPxPxPxPx a != b", cPxPxPxPx
      cPxPxPxPx = coulombPxPxPxPx(a,a,z1,z2) - (0.47343555297830386d0)
      print*, "cPxPxPxPx a = b", cPxPxPxPx
      cPxPxPxPx = coulombPxPxPxPx(b,a,z1,z2) - (0.46922410322387775d0)
      print*, "cPxPxPxPx b != a", cPxPxPxPx
      print*, ""
c
c     coulombSPxPxPz
c
      cSPxPxPz = coulombSPxPxPz(a,b,z1,z2) - (0.008359792145233679d0)
      print*, "cSPxPxPz a != b", cSPxPxPz
      cSPxPxPz = coulombSPxPxPz(a,a,z1,z2) - (0.0063400196190138685d0)
      print*, "cSPxPxPz a = b", cSPxPxPz
      cSPxPxPz = coulombSPxPxPz(b,a,z1,z2) - (0.007228374045474331d0)
      print*, "cSPxPxPz b != a", cSPxPxPz
      print*, ""
c
c     coulombSPzPxPx
c
      cSPzPxPx = coulombSPzPxPx(a,b,z1,z2) - (-0.062332604510992d0)
      print*, "cSPzPxPx a != b", cSPzPxPx
      cSPzPxPx = coulombSPzPxPx(a,a,z1,z2) - (-0.06473232527751177d0)
      print*, "cSPzPxPx a = b", cSPzPxPx
      cSPzPxPx = coulombSPzPxPx(b,a,z1,z2) - (-0.07524608872012231d0)
      print*, "cSPzPxPx b != a", cSPzPxPx
      print*, ""
c
c     coulombPxPxSPz
c
      cPxPxSPz = coulombPxPxSPz(a,b,z1,z2) - (0.07524608872012231d0)
      print*, "cPxPxSPz a != b", cPxPxSPz
      cPxPxSPz = coulombPxPxSPz(a,a,z1,z2) - (0.06473232527751177d0)
      print*, "cPxPxSPz a = b", cPxPxSPz
      cPxPxSPz = coulombPxPxSPz(b,a,z1,z2) - (0.062332604510992d0)
      print*, "cPxPxSPz b != a", cPxPxSPz
      print*, ""
c
c     coulombPxPxPzPz
c
      cPxPxPzPz = coulombPxPxPzPz(a,b,z1,z2) - (0.5164642045478991d0)
      print*, "cPxPxPzPz a != b", cPxPxPzPz
      cPxPxPzPz = coulombPxPxPzPz(a,a,z1,z2) - (0.5089658284234035d0)
      print*, "cPxPxPzPz a = b", cPxPxPzPz
      cPxPxPzPz = coulombPxPxPzPz(b,a,z1,z2) - (0.5013159038041091d0)
      print*, "cPxPxPzPz b != a", cPxPxPzPz
      print*, ""
c
c     coulombPxPzSPx
c
      cPxPzSPx = coulombPxPzSPx(a,b,z1,z2) - (-0.007228374045474331d0)
      print*, "cPxPzSPx a != b", cPxPzSPx
      cPxPzSPx = coulombPxPzSPx(a,a,z1,z2) - (-0.0063400196190138685d0)
      print*, "cPxPzSPx a = b", cPxPzSPx
      cPxPzSPx = coulombPxPzSPx(b,a,z1,z2) - (-0.008359792145233679d0)
      print*, "cPxPzSPx b != a", cPxPzSPx
      print*, ""
c
c     coulombPxPzPxPz
c
      cPxPzPxPz = coulombPxPzPxPz(a,b,z1,z2) - (-0.006034801127456413d0)
      print*, "cPxPzPxPz a != b", cPxPzPxPz
      cPxPzPxPz = coulombPxPzPxPz(a,a,z1,z2) - (-0.004848437602982191d0)
      print*, "cPxPzPxPz a = b", cPxPzPxPz
      cPxPzPxPz = coulombPxPzPxPz(b,a,z1,z2) - (-0.006034801127456413d0)
      print*, "cPxPzPxPz b != a", cPxPzPxPz
      print*, ""
c
c     coulombPzPzPxPx
c
      cPzPzPxPx = coulombPzPzPxPx(a,b,z1,z2) - (0.5013159038041091d0)
      print*, "cPzPzPxPx a != b", cPzPzPxPx
      cPzPzPxPx = coulombPzPzPxPx(a,a,z1,z2) - (0.5089658284234035d0)
      print*, "cPzPzPxPx a = b", cPzPzPxPx
      cPzPzPxPx = coulombPzPzPxPx(b,a,z1,z2) - (0.5164642045478991d0)
      print*, "cPzPzPxPx b != a", cPzPzPxPx
      print*, ""
c
c     coulombPxPxPyPy
c
      cPxPxPyPy = coulombPxPxPyPy(a,b,z1,z2) - (0.46572834151377246d0)
      print*, "cPxPxPyPy a != b", cPxPxPyPy
      cPxPxPyPy = coulombPxPxPyPy(a,a,z1,z2) - (0.47076793678829376d0)
      print*, "cPxPxPyPy a = b", cPxPxPyPy
      cPxPxPyPy = coulombPxPxPyPy(b,a,z1,z2) - (0.4657283415137724d0)
      print*, "cPxPxPyPy b != a", cPxPxPyPy
      print*, ""
c
c     coulombPxPyPxPy
c
      cPxPyPxPy = coulombPxPyPxPy(a,b,z1,z2) - (0.0017478808550526867d0)
      print*, "cPxPyPxPy a != b", cPxPyPxPy
      cPxPyPxPy = coulombPxPyPxPy(a,a,z1,z2) - (0.001333808095005061d0)
      print*, "cPxPyPxPy a = b", cPxPyPxPy
      cPxPyPxPy = coulombPxPyPxPy(b,a,z1,z2) - (0.0017478808550526867d0)
      print*, "cPxPyPxPy b != a", cPxPyPxPy
      print*, ""
      return
      end
c
c
      subroutine testcoulombsto()
      implicit none
      real*8 a,b,z1,z2
      real*8 SSSS,sSSSS
      real*8 SSSPz,sSSSPz
      real*8 SSPzPz,sSSPzPz
      real*8 SPzSS,sSPzSS
      real*8 SPzSPz,sSPzSPz
      real*8 SPzPzPz,sSPzPzPz
      real*8 PzPzSS,sPzPzSS
      real*8 PzPzSPz,sPzPzSPz
      real*8 PzPzPzPz,sPzPzPzPz
      real*8 SSPxPx,sSSPxPx
      real*8 SPxSPx,sSPxSPx
      real*8 PxPxSS,sPxPxSS
      real*8 PxPxPxPx,sPxPxPxPx
      real*8 SPxPxPz,sSPxPxPz
      real*8 SPzPxPx,sSPzPxPx
      real*8 PxPxSPz,sPxPxSPz
      real*8 PxPzSPx,sPxPzSPx
      real*8 PxPxPzPz,sPxPxPzPz
      real*8 PxPzPxPz,sPxPzPxPz
      real*8 PzPzPxPx,sPzPzPxPx
      real*8 PxPxPyPy,sPxPxPyPy
      real*8 PxPyPxPy,sPxPyPxPy
      logical exch
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      exch = .false.
      print*, "STO Coulomb"
c
c     stoSSSS
c
      sSSSS = SSSS(a,b,z1,z2,exch) - (0.4998928969427606d0)
      print*, "sSSSS a != b", sSSSS
      sSSSS = SSSS(a,a,z1,z2,exch) - (0.4999610519241599d0)
      print*, "sSSSS a = a", sSSSS
      sSSSS = SSSS(b,a,z1,z2,exch) - (0.4998928969427606d0)
      print*, "sSSSS b != a", sSSSS
      print*, ""
c
c     stoSSSPz
c
      sSSSPz = SSSPz(a,b,z1,z2,exch) - (-0.08280497127320369d0)
      print*, "sSSSPz a != b", sSSSPz
      sSSSPz = SSSPz(a,a,z1,z2,exch) - (-0.07123901220527948d0)
      print*, "sSSSPz a = a", sSSSPz
      sSSSPz = SSSPz(b,a,z1,z2,exch) - (-0.07107200980311534d0)
      print*, "sSSSPz b != a", sSSSPz
      print*, ""
c
c     stoSSPzPz
c
      sSSPzPz = SSPzPz(a,b,z1,z2,exch) - (0.5388052815345308d0)
      print*, "sSSPzPz a != b", sSSPzPz
      sSSPzPz = SSPzPz(a,a,z1,z2,exch) - (0.5296530877663741d0)
      print*, "sSSPzPz a = a", sSSPzPz
      sSSPzPz = SSPzPz(b,a,z1,z2,exch) - (0.5290687427450043d0)
      print*, "sSSPzPz b != a", sSSPzPz
      print*, ""
c
c     stoSPzSS
c
      sSPzSS = SPzSS(a,b,z1,z2,exch) - (0.07107200980311536d0)
      print*, "sSPzSS a != b", sSPzSS
      sSPzSS = SPzSS(a,a,z1,z2,exch) - (0.0712390122052794d0)
      print*, "sSPzSS a = a", sSPzSS
      sSPzSS = SPzSS(b,a,z1,z2,exch) - (0.0828049712732037d0)
      print*, "sSPzSS b != a", sSPzSS
      print*, ""
c
c     stoSPzSPz
c
      sSPzSPz = SPzSPz(a,b,z1,z2,exch) - (-0.022587828521005462d0)
      print*, "sSPzSPz a != b", sSPzSPz
      sSPzSPz = SPzSPz(a,a,z1,z2,exch) - (-0.019849667288250994d0)
      print*, "sSPzSPz a = a", sSPzSPz
      sSPzSPz = SPzSPz(b,a,z1,z2,exch) - (-0.022587828521005465d0)
      print*, "sSPzSPz b != a", sSPzSPz
      print*, ""
c
c     stoSPzPzPz
c
      sSPzPzPz = SPzPzPz(a,b,z1,z2,exch) - (0.08385729862671476d0)
      print*, "sSPzPzPz a != b", sSPzPzPz
      sSPzPzPz = SPzPzPz(a,a,z1,z2,exch) - (0.08210733195729732d0)
      print*, "sSPzPzPz a = a", sSPzPzPz
      sSPzPzPz = SPzPzPz(b,a,z1,z2,exch) - (0.09416347750670924d0)
      print*, "sSPzPzPz b != a", sSPzPzPz
      print*, ""
c
c     stoPzPzSS
c
      sPzPzSS = PzPzSS(a,b,z1,z2,exch) - (0.5290687427450047d0)
      print*, "sPzPzSS a != b", sPzPzSS
      sPzPzSS = PzPzSS(a,a,z1,z2,exch) - (0.5296530877663737d0)
      print*, "sPzPzSS a = a", sPzPzSS
      sPzPzSS = PzPzSS(b,a,z1,z2,exch) - (0.5388052815345293d0)
      print*, "sPzPzSS b != a", sPzPzSS
      print*, ""
c
c     stoPzPzSPz
c
      sPzPzSPz = PzPzSPz(a,b,z1,z2,exch) - (-0.09416347750670932d0)
      print*, "sPzPzSPz a != b", sPzPzSPz
      sPzPzSPz = PzPzSPz(a,a,z1,z2,exch) - (-0.08210733195729736d0)
      print*, "sPzPzSPz a = a", sPzPzSPz
      sPzPzSPz = PzPzSPz(b,a,z1,z2,exch) - (-0.0838572986267148d0)
      print*, "sPzPzSPz b != a", sPzPzSPz
      print*, ""
c
c     stoPzPzPzPz
c
      sPzPzPzPz = PzPzPzPz(a,b,z1,z2,exch) - (0.5703602232449699d0)
      print*, "sPzPzPzPz a != b", sPzPzPzPz
      sPzPzPzPz = PzPzPzPz(a,a,z1,z2,exch) - (0.5634966290593179d0)
      print*, "sPzPzPzPz a = a", sPzPzPzPz
      sPzPzPzPz = PzPzPzPz(b,a,z1,z2,exch) - (0.5703602232449705d0)
      print*, "sPzPzPzPz b != a", sPzPzPzPz
      print*, ""
c
c     stoSSPxPx
c
      sSSPxPx = SSPxPx(a,b,z1,z2,exch) - (0.479121882227756d0)
      print*, "sSSPxPx a != b", sSSPxPx
      sSSPxPx = SSPxPx(a,a,z1,z2,exch) - (0.48465855281958337d0)
      print*, "sSSPxPx a = a", sSSPxPx
      sSSPxPx = SSPxPx(b,a,z1,z2,exch) - (0.4845951729495459d0)
      print*, "sSSPxPx b != a", sSSPxPx
      print*, ""
c
c     stoSPxSPx
c
      sSPxSPx = SPxSPx(a,b,z1,z2,exch) - (0.011779470975834882d0)
      print*, "sSPxSPx a != b", sSPxSPx
      sSPxSPx = SPxSPx(a,a,z1,z2,exch) - (0.010150174909361135d0)
      print*, "sSPxSPx a = a", sSPxSPx
      sSPxSPx = SPxSPx(b,a,z1,z2,exch) - (0.01177947097583488d0)
      print*, "sSPxSPx b != a", sSPxSPx
      print*, ""
c
c     stoPxPxSS
c
      sPxPxSS = PxPxSS(a,b,z1,z2,exch) - (0.48459517294954596d0)
      print*, "sPxPxSS a != b", sPxPxSS
      sPxPxSS = PxPxSS(a,a,z1,z2,exch) - (0.48465855281958353d0)
      print*, "sPxPxSS a = a", sPxPxSS
      sPxPxSS = PxPxSS(b,a,z1,z2,exch) - (0.4791218822277561d0)
      print*, "sPxPxSS b != a", sPxPxSS
      print*, ""
c
c     stoPxPxPxPx
c
      sPxPxPxPx = PxPxPxPx(a,b,z1,z2,exch) - (0.46920846816656386d0)
      print*, "sPxPxPxPx a != b", sPxPxPxPx
      sPxPxPxPx = PxPxPxPx(a,a,z1,z2,exch) - (0.47342184437293383d0)
      print*, "sPxPxPxPx a = a", sPxPxPxPx
      sPxPxPxPx = PxPxPxPx(b,a,z1,z2,exch) - (0.4692084681665642d0)
      print*, "sPxPxPxPx b != a", sPxPxPxPx
      print*, ""
c
c     stoSPxPxPz
c
      sSPxPxPz = SPxPxPz(a,b,z1,z2,exch) - (-0.008357109447580174d0)
      print*, "sSPxPxPz a != b", sSPxPxPz
      sSPxPxPz = SPxPxPz(a,a,z1,z2,exch) - (-0.006339839532298857d0)
      print*, "sSPxPxPz a = a", sSPxPxPz
      sSPxPxPz = SPxPxPz(b,a,z1,z2,exch) - (-0.007228139438835614d0)
      print*, "sSPxPxPz b != a", sSPxPxPz
      print*, ""
c
c     stoSPzPxPx
c
      sSPzPxPx = SPzPxPx(a,b,z1,z2,exch) - (0.06232047248635676d0)
      print*, "sSPzPxPx a != b", sSPzPxPx
      sSPzPxPx = SPzPxPx(a,a,z1,z2,exch) - (0.06471909488278875d0)
      print*, "sSPzPxPx a = a", sSPzPxPx
      sSPzPxPx = SPzPxPx(b,a,z1,z2,exch) - (0.07523723847649672d0)
      print*, "sSPzPxPx b != a", sSPzPxPx
      print*, ""
c
c     stoPxPxSPz
c
      sPxPxSPz = PxPxSPz(a,b,z1,z2,exch) - (-0.07523723847649681d0)
      print*, "sPxPxSPz a != b", sPxPxSPz
      sPxPxSPz = PxPxSPz(a,a,z1,z2,exch) - (-0.0647190948827888d0)
      print*, "sPxPxSPz a = a", sPxPxSPz
      sPxPxSPz = PxPxSPz(b,a,z1,z2,exch) - (-0.0623204724863568d0)
      print*, "sPxPxSPz b != a", sPxPxSPz
      print*, ""
c
c     stoPxPzSPx
c
      sPxPzSPx = PxPzSPx(a,b,z1,z2,exch) - (0.007228139438835614d0)
      print*, "sPxPzSPx a != b", sPxPzSPx
      sPxPzSPx = PxPzSPx(a,a,z1,z2,exch) - (0.00633983953229886d0)
      print*, "sPxPzSPx a = a", sPxPzSPx
      sPxPzSPx = PxPzSPx(b,a,z1,z2,exch) - (0.008357109447580172d0)
      print*, "sPxPzSPx b != a", sPxPzSPx
      print*, ""
c
c     stoPxPxPzPz
c
      sPxPxPzPz = PxPxPzPz(a,b,z1,z2,exch) - (0.5164311229070725d0)
      print*, "sPxPxPzPz a != b", sPxPxPzPz
      sPxPxPzPz = PxPxPzPz(a,a,z1,z2,exch) - (0.5089515640625875d0)
      print*, "sPxPxPzPz a = a", sPxPxPzPz
      sPxPxPzPz = PxPxPzPz(b,a,z1,z2,exch) - (0.5013023902190399d0)
      print*, "sPxPxPzPz b != a", sPxPxPzPz
      print*, ""
c
c     stoPxPzPxPz
c
      sPxPzPxPz = PxPzPxPz(a,b,z1,z2,exch) - (-0.006031364166773978d0)
      print*, "sPxPzPxPz a != b", sPxPzPxPz
      sPxPzPxPz = PxPzPxPz(a,a,z1,z2,exch) - (-0.004843745148161869d0)
      print*, "sPxPzPxPz a = a", sPxPzPxPz
      sPxPzPxPz = PxPzPxPz(b,a,z1,z2,exch) - (-0.006031364166773985d0)
      print*, "sPxPzPxPz b != a", sPxPzPxPz
      print*, ""
c
c     stoPzPzPxPx
c
      sPzPzPxPx = PzPzPxPx(a,b,z1,z2,exch) - (0.5013023902190399d0)
      print*, "sPzPzPxPx a != b", sPzPzPxPx
      sPzPzPxPx = PzPzPxPx(a,a,z1,z2,exch) - (0.5089515640625878d0)
      print*, "sPzPzPxPx a = a", sPzPzPxPx
      sPzPzPxPx = PzPzPxPx(b,a,z1,z2,exch) - (0.5164311229070717d0)
      print*, "sPzPzPxPx b != a", sPzPzPxPx
      print*, ""
c
c     stoPxPxPyPy
c
      sPxPxPyPy = PxPxPyPy(a,b,z1,z2,exch) - (0.46571361092485064d0)
      print*, "sPxPxPyPy a != b", sPxPxPyPy
      sPxPxPyPy = PxPxPyPy(a,a,z1,z2,exch) - (0.470754460919562d0)
      print*, "sPxPxPyPy a = a", sPxPxPyPy
      sPxPxPyPy = PxPxPyPy(b,a,z1,z2,exch) - (0.46571361092485014d0)
      print*, "sPxPxPyPy b != a", sPxPxPyPy
      print*, ""
c
c     stoPxPyPxPy
c
      sPxPyPxPy = PxPyPxPy(a,b,z1,z2,exch) - (0.0017474286208569935d0)
      print*, "sPxPyPxPy a != b", sPxPyPxPy
      sPxPyPxPy = PxPyPxPy(a,a,z1,z2,exch) - (0.0013336917266853177d0)
      print*, "sPxPyPxPy a = a", sPxPyPxPy
      sPxPyPxPy = PxPyPxPy(b,a,z1,z2,exch) - (0.0017474286208569929d0)
      print*, "sPxPyPxPy b != a", sPxPyPxPy
      print*, ""
      return
      end
c
c
      subroutine testexchange()
      implicit none
      real*8 a,b,z1,z2
      real*8 SSSS, eSSSS
      real*8 SSPzS, eSSPzS
      real*8 SSSPz, eSSSPz
      real*8 SSPzPz, eSSPzPz
      real*8 SPzSPz, eSPzSPz
      real*8 SPzPzS, eSPzPzS
      real*8 SPzPzPz, eSPzPzPz
      real*8 PzSPzS, ePzSPzS
      real*8 PzSPzPz, ePzSPzPz
      real*8 PzPzPzPz, ePzPzPzPz
      real*8 SSPxPx, eSSPxPx
      real*8 SPxSPx, eSPxSPx
      real*8 SPxPxS, eSPxPxS
      real*8 PxSPxS, ePxSPxS
      real*8 PxPxPxPx, ePxPxPxPx
      real*8 SPxPxPz, eSPxPxPz
      real*8 SPxPzPx, eSPxPzPx
      real*8 SPzPxPx, eSPzPxPx
      real*8 PxSPxPz, ePxSPxPz
      real*8 PxSPzPx, ePxSPzPx
      real*8 PzSPxPx, ePzSPxPx
      real*8 PxPxPzPz, ePxPxPzPz
      real*8 PxPzPxPz, ePxPzPxPz
      real*8 PxPzPzPx, ePxPzPzPx
      real*8 PzPxPzPx, ePzPxPzPx
      real*8 PxPxPyPy, ePxPxPyPy
      real*8 PxPyPxPy, ePxPyPxPy
      real*8 PxPyPyPx, ePxPyPyPx
      logical exch
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      exch = .true.
      print*, "STO Exchange"
c
c     exchangeSSSS
c
      eSSSS = SSSS(a,b,z1,z2,exch) - (0.0012601594158343323d0)
      print*, "eSSSS a != b", eSSSS
      eSSSS = SSSS(a,a,z1,z2,exch) - (0.000560943327754073d0)
      print*, "eSSSS a = a", eSSSS
      eSSSS = SSSS(b,a,z1,z2,exch) - (0.001260159415834332d0)
      print*, "eSSSS b != a", eSSSS
      print*, ""
c
c     exchangeSSPzS
c
      eSSPzS = SSPzS(a,b,z1,z2,exch) - (0.003681916344317933d0)
      print*, "eSSPzS a != b", eSSPzS
      eSSPzS = SSPzS(a,a,z1,z2,exch)  - (0.002022526132079656d0)
      print*, "eSSPzS a = a", eSSPzS
      eSSPzS = SSPzS(b,a,z1,z2,exch) - (0.004483574272439649d0)
      print*, "eSSPzS b != a", eSSPzS
      print*, ""
c
c     exchangeSSSPz
c
      eSSSPz = SSSPz(a,b,z1,z2,exch) - (-0.004483574272439654d0)
      print*, "eSSSPz a != b", eSSSPz
      eSSSPz = SSSPz(a,a,z1,z2,exch)  - (-0.0020225261320796556d0)
      print*, "eSSSPz a = a", eSSSPz
      eSSSPz = SSSPz(b,a,z1,z2,exch) - (-0.003681916344317934d0)
      print*, "eSSSPz b != a", eSSSPz
      print*, ""
c
c     exchangeSSPzPz
c
      eSSPzPz = SSPzPz(a,b,z1,z2,exch) - (-0.008628848714794682d0)
      print*, "eSSPzPz a != b", eSSPzPz
      eSSPzPz = SSPzPz(a,a,z1,z2,exch)  - (-0.0047951167424194075d0)
      print*, "eSSPzPz a = a", eSSPzPz
      eSSPzPz = SSPzPz(b,a,z1,z2,exch) - (-0.008628848714794682d0)
      print*, "eSSPzPz b != a", eSSPzPz
      print*, ""
c
c     exchangeSPzSPz
c
      eSPzSPz = SPzSPz(a,b,z1,z2,exch) - (0.01692508955217397d0)
      print*, "eSPzSPz a != b", eSPzSPz
      eSPzSPz = SPzSPz(a,a,z1,z2,exch)  - (0.007979784495560714d0)
      print*, "eSPzSPz a = a", eSPzSPz
      eSPzSPz = SPzSPz(b,a,z1,z2,exch) - (0.012244871046188573d0)
      print*, "eSPzSPz b != a", eSPzSPz
      print*, ""
c
c     exchangeSPzPzS
c
      eSPzPzS = SPzPzS(a,b,z1,z2,exch) - (-0.011907945489374907d0)
      print*, "eSPzPzS a != b", eSPzPzS
      eSPzPzS = SPzPzS(a,a,z1,z2,exch)  - (-0.006621271704245987d0)
      print*, "eSPzPzS a = a", eSPzPzS
      eSPzPzS = SPzPzS(b,a,z1,z2,exch) - (-0.011907945489374907d0)
      print*, "eSPzPzS b != a", eSPzPzS
      print*, ""
c
c     exchangeSPzPzPz
c
      eSPzPzPz = SPzPzPz(a,b,z1,z2,exch) - (0.029737406917215128d0)
      print*, "eSPzPzPz a != b", eSPzPzPz
      eSPzPzPz = SPzPzPz(a,a,z1,z2,exch)  - (0.01720058618401024d0)
      print*, "eSPzPzPz a = a", eSPzPzPz
      eSPzPzPz = SPzPzPz(b,a,z1,z2,exch) - (0.026228929037370376d0)
      print*, "eSPzPzPz b != a", eSPzPzPz
      print*, ""
c
c     exchangePzSPzS
c
      ePzSPzS = PzSPzS(a,b,z1,z2,exch) - (0.012244871046188575d0)
      print*, "ePzSPzS a != b", ePzSPzS
      ePzSPzS = PzSPzS(a,a,z1,z2,exch)  - (0.007979784495560718d0)
      print*, "ePzSPzS a = a", ePzSPzS
      ePzSPzS = PzSPzS(b,a,z1,z2,exch) - (0.01692508955217398d0)
      print*, "ePzSPzS b != a", ePzSPzS
      print*, ""
c
c     exchangePzSPzPz
c
      ePzSPzPz = PzSPzPz(a,b,z1,z2,exch) - (-0.026228929037370365d0)
      print*, "ePzSPzPz a != b", ePzSPzPz
      ePzSPzPz = PzSPzPz(a,a,z1,z2,exch)  - (-0.017200586184010234d0)
      print*, "ePzSPzPz a = a", ePzSPzPz
      ePzSPzPz = PzSPzPz(b,a,z1,z2,exch) - (-0.029737406917215124d0)
      print*, "ePzSPzPz b != a", ePzSPzPz
      print*, ""
c
c     exchangePzPzPzPz
c
      ePzPzPzPz = PzPzPzPz(a,b,z1,z2,exch) - (0.06253163254020173d0)
      print*, "ePzPzPzPz a != b", ePzPzPzPz
      ePzPzPzPz = PzPzPzPz(a,a,z1,z2,exch)  - (0.04237787847359695d0)
      print*, "ePzPzPzPz a = a", ePzPzPzPz
      ePzPzPzPz = PzPzPzPz(b,a,z1,z2,exch) - (0.06253163254020176d0)
      print*, "ePzPzPzPz b != a", ePzPzPzPz
      print*, ""
c
c     exchangeSSPxPx
c
      eSSPxPx = SSPxPx(a,b,z1,z2,exch) - (0.0020032104297817988d0)
      print*, "eSSPxPx a != b", eSSPxPx
      eSSPxPx = SSPxPx(a,a,z1,z2,exch)  - (0.0009761184659617603d0)
      print*, "eSSPxPx a = a", eSSPxPx
      eSSPxPx = SSPxPx(b,a,z1,z2,exch) - (0.002003210429781801d0)
      print*, "eSSPxPx b != a", eSSPxPx
      print*, ""
c
c     exchangeSPxSPx
c
      eSPxSPx = SPxSPx(a,b,z1,z2,exch) - (0.0002774061712578313d0)
      print*, "eSPxSPx a != b", eSPxSPx
      eSPxSPx = SPxSPx(a,a,z1,z2,exch)  - (0.0001586806871235956d0)
      print*, "eSPxSPx a = a", eSPxSPx
      eSPxSPx = SPxSPx(b,a,z1,z2,exch) - (0.0003846727099939042d0)
      print*, "eSPxSPx b != a", eSPxSPx
      print*, ""
c
c     exchangeSPxPxS
c
      eSPxPxS = SPxPxS(a,b,z1,z2,exch) - (0.00032647074114874094d0)
      print*, "eSPxPxS a != b", eSPxPxS
      eSPxPxS = SPxPxS(a,a,z1,z2,exch)  - (0.00015845001211751495d0)
      print*, "eSPxPxS a = a", eSPxPxS
      eSPxPxS = SPxPxS(b,a,z1,z2,exch) - (0.000326470741148741d0)
      print*, "eSPxPxS b != a", eSPxPxS
      print*, ""
c
c     exchangePxSPxS
c
      ePxSPxS = PxSPxS(a,b,z1,z2,exch) - (0.00038467270999390424d0)
      print*, "ePxSPxS a != b", ePxSPxS
      ePxSPxS = PxSPxS(a,a,z1,z2,exch)  - (0.00015868068712359556d0)
      print*, "ePxSPxS a = a", ePxSPxS
      ePxSPxS = PxSPxS(b,a,z1,z2,exch) - (0.000277406171257831d0)
      print*, "ePxSPxS b != a", ePxSPxS
      print*, ""
c
c     exchangePxPxPxPx
c
      ePxPxPxPx = PxPxPxPx(a,b,z1,z2,exch) - (0.003765299410180266d0)
      print*, "ePxPxPxPx a != b", ePxPxPxPx
      ePxPxPxPx = PxPxPxPx(a,a,z1,z2,exch)  - (0.0019991582352288694d0)
      print*, "ePxPxPxPx a = a", ePxPxPxPx
      ePxPxPxPx = PxPxPxPx(b,a,z1,z2,exch) - (0.003765299410180263d0)
      print*, "ePxPxPxPx b != a", ePxPxPxPx
      print*, ""
c
c     exchangeSPxPxPz
c
      eSPxPxPz = SPxPxPz(a,b,z1,z2,exch) - (-0.0011314383048050068d0)
      print*, "eSPxPxPz a != b", eSPxPxPz
      eSPxPxPz = SPxPxPz(a,a,z1,z2,exch)  - (-0.0005642253147881832d0)
      print*, "eSPxPxPz a = a", eSPxPxPz
      eSPxPxPz = SPxPxPz(b,a,z1,z2,exch) - (-0.0009760376315431994d0)
      print*, "eSPxPxPz b != a", eSPxPxPz
      print*, ""
c
c     exchangeSPxPzPx
c
      eSPxPzPx = SPxPzPx(a,b,z1,z2,exch) - (0.0008241859222528157d0)
      print*, "eSPxPzPx a != b", eSPxPzPx
      eSPxPzPx = SPxPzPx(a,a,z1,z2,exch)  - (0.0005571195789497128d0)
      print*, "eSPxPzPx a = a", eSPxPzPx
      eSPxPzPx = SPxPzPx(b,a,z1,z2,exch) - (0.0013275177608087609d0)
      print*, "eSPxPzPx b != a", eSPxPzPx
      print*, ""
c
c     exchangeSPzPxPx
c
      eSPzPxPx = SPzPxPx(a,b,z1,z2,exch) - (-0.0070491060326243445d0)
      print*, "eSPzPxPx a != b", eSPzPxPx
      eSPzPxPx = SPzPxPx(a,a,z1,z2,exch)  - (-0.0035322353520320722d0)
      print*, "eSPzPxPx a = a", eSPzPxPx
      eSPzPxPx = SPzPxPx(b,a,z1,z2,exch) - (-0.005983763782569305d0)
      print*, "eSPzPxPx b != a", eSPzPxPx
      print*, ""
c
c     exchangePxSPxPz
c
      ePxSPxPz = PxSPxPz(a,b,z1,z2,exch) - (-0.0013275177608087607d0)
      print*, "ePxSPxPz a != b", ePxSPxPz
      ePxSPxPz = PxSPxPz(a,a,z1,z2,exch)  - (-0.000557119578949713d0)
      print*, "ePxSPxPz a = a", ePxSPxPz
      ePxSPxPz = PxSPxPz(b,a,z1,z2,exch) - (-0.0008241859222528157d0)
      print*, "ePxSPxPz b != a", ePxSPxPz
      print*, ""
c
c     exchangePxSPzPx
c
      ePxSPzPx = PxSPzPx(a,b,z1,z2,exch) - (0.0009760376315431994d0)
      print*, "ePxSPzPx a != b", ePxSPzPx
      ePxSPzPx = PxSPzPx(a,a,z1,z2,exch)  - (0.0005642253147881825d0)
      print*, "ePxSPzPx a = a", ePxSPzPx
      ePxSPzPx = PxSPzPx(b,a,z1,z2,exch) - (0.001131438304805008d0)
      print*, "ePxSPzPx b != a", ePxSPzPx
      print*, ""
c
c     exchangePzSPxPx
c
      ePzSPxPx = PzSPxPx(a,b,z1,z2,exch) - (0.005983763782569308d0)
      print*, "ePzSPxPx a != b", ePzSPxPx
      ePzSPxPx = PzSPxPx(a,a,z1,z2,exch)  - (0.0035322353520320722d0)
      print*, "ePzSPxPx a = a", ePzSPxPx
      ePzSPxPx = PzSPxPx(b,a,z1,z2,exch) - (0.007049106032624344d0)
      print*, "ePzSPxPx b != a", ePzSPxPx
      print*, ""
c
c     exchangePxPxPzPz
c
      ePxPxPzPz = PxPxPzPz(a,b,z1,z2,exch) - (-0.013694893664449407d0)
      print*, "ePxPxPzPz a != b", ePxPxPzPz
      ePxPxPzPz = PxPxPzPz(a,a,z1,z2,exch)  - (-0.00829283834156019d0)
      print*, "ePxPxPzPz a = a", ePxPxPzPz
      ePxPxPzPz = PxPxPzPz(b,a,z1,z2,exch) - (-0.0136948936644494d0)
      print*, "ePxPxPzPz b != a", ePxPxPzPz
      print*, ""
c
c     exchangePxPzPxPz
c
      ePxPzPxPz = PxPzPxPz(a,b,z1,z2,exch) - (0.005073439663329332d0)
      print*, "ePxPzPxPz a != b", ePxPzPxPz
      ePxPzPxPz = PxPzPxPz(a,a,z1,z2,exch)  - (0.00226940825085089d0)
      print*, "ePxPzPxPz a = a", ePxPzPxPz
      ePxPzPxPz = PxPzPxPz(b,a,z1,z2,exch) - (0.0029683923453331248d0)
      print*, "ePxPzPxPz b != a", ePxPzPxPz
      print*, ""
c
c     exchangePxPzPzPx
c
      ePxPzPzPx = PxPzPzPx(a,b,z1,z2,exch) - (-0.0028786346120677412d0)
      print*, "ePxPzPzPx a != b", ePxPzPzPx
      ePxPzPzPx = PxPzPzPx(a,a,z1,z2,exch)  - (-0.0016982516334179468d0)
      print*, "ePxPzPzPx a = a", ePxPzPzPx
      ePxPzPzPx = PxPzPzPx(b,a,z1,z2,exch) - (-0.0028786346120677417d0)
      print*, "ePxPzPzPx b != a", ePxPzPzPx
      print*, ""
c
c     exchangePzPxPzPx
c
      ePzPxPzPx = PzPxPzPx(a,b,z1,z2,exch) - (0.002968392345333127d0)
      print*, "ePzPxPzPx a != b", ePzPxPzPx
      ePzPxPzPx = PzPxPzPx(a,a,z1,z2,exch)  - (0.002269408250850889d0)
      print*, "ePzPxPzPx a = a", ePzPxPzPx
      ePzPxPzPx = PzPxPzPx(b,a,z1,z2,exch) - (0.005073439663329334d0)
      print*, "ePzPxPzPx b != a", ePzPxPzPx
      print*, ""
c
c     exchangePxPxPyPy
c
      ePxPxPyPy = PxPxPyPy(a,b,z1,z2,exch) - (0.0033861961455178813d0)
      print*, "ePxPxPyPy a != b", ePxPxPyPy
      ePxPxPyPy = PxPxPyPy(a,a,z1,z2,exch)  - (0.0017998204940336629d0)
      print*, "ePxPxPyPy a = a", ePxPxPyPy
      ePxPxPyPy = PxPxPyPy(b,a,z1,z2,exch) - (0.0033861961455178805d0)
      print*, "ePxPxPyPy b != a", ePxPxPyPy
      print*, ""
c
c     exchangePxPyPxPy
c
      ePxPyPxPy = PxPyPxPy(a,b,z1,z2,exch) - (0.00018955163233119358d0)
      print*, "ePxPyPxPy a != b", ePxPyPxPy
      ePxPyPxPy = PxPyPxPy(a,a,z1,z2,exch)  - (0.00009966887059760348d0)
      print*, "ePxPyPxPy a = a", ePxPyPxPy
      ePxPyPxPy = PxPyPxPy(b,a,z1,z2,exch) - (0.00018955163233119366d0)
      print*, "ePxPyPxPy b != a", ePxPyPxPy
      print*, ""
c
c     exchangePxPyPyPx
c
      ePxPyPyPx = PxPyPyPx(a,b,z1,z2,exch) - (0.00018955163233119358d0)
      print*, "ePxPyPyPx a != b", ePxPyPyPx
      ePxPyPyPx = PxPyPyPx(a,a,z1,z2,exch)  - (0.00009966887059760348d0)
      print*, "ePxPyPyPx a = a", ePxPyPyPx
      ePxPyPyPx = PxPyPyPx(b,a,z1,z2,exch) - (0.00018955163233119366d0)
      print*, "ePxPyPyPx b != a", ePxPyPyPx
      print*, ""
      return
      end
c
c
      subroutine teststooverlap()
      implicit none
      real*8 a,b,z1,z2
      real*8 stoOSS, sSS
      real*8 stoOSPz, sSPz
      real*8 stoOPzS, sPzS
      real*8 stoOPzPz, sPzPz
      real*8 stoOPxPx, sPxPx
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      print*, "STO Overlap"
c
c     stoOSS
c
      sSS = stoOSS(a,b,z1,z2) - (0.03259146113402994d0)
      print*, "sSS a != b", sSS
      sSS = stoOSS(a,a,z1,z2) - (0.02136399226329799d0)
      print*, "sSS a = a", sSS
      sSS = stoOSS(b,a,z1,z2) - (0.03259146113402994d0)
      print*, "sSS b != a", sSS
      print*, ""
c
c     stoOSPz
c
      sSPz = stoOSPz(a,b,z1,z2) - (-0.11491683162998949d0)
      print*, "sSPz a != b", sSPz
      sSPz = stoOSPz(a,a,z1,z2) - (-0.07850320352440002d0)
      print*, "sSPz a = a", sSPz
      sSPz = stoOSPz(b,a,z1,z2) - (-0.09952502840177815d0)
      print*, "sSPz b != a", sSPz
      print*, ""
c
c     stoOPzS
c
      sPzS = stoOPzS(a,b,z1,z2) - (0.09952502840177814d0)
      print*, "sPzS a != b", sPzS
      sPzS = stoOPzS(a,a,z1,z2) - (0.07850320352440004d0)
      print*, "sPzS a = a", sPzS
      sPzS = stoOPzS(b,a,z1,z2) - (0.1149168316299895d0)
      print*, "sPzS b != a", sPzS
      print*, ""
c
c     stoOPzPz
c
      sPzPz = stoOPzPz(a,b,z1,z2) - (-0.2079653261681678d0)
      print*, "sPzPz a != b", sPzPz
      sPzPz = stoOPzPz(a,a,z1,z2) - (-0.17110584149438524d0)
      print*, "sPzPz a = a", sPzPz
      sPzPz = stoOPzPz(b,a,z1,z2) - (-0.20796532616816776d0)
      print*, "sPzPz b != a", sPzPz
      print*, ""
c
c     stoOPxPx
c
      sPxPx = stoOPxPx(a,b,z1,z2) - (0.06483001236088108d0)
      print*, "sPxPx a != b", sPxPx
      sPxPx = stoOPxPx(a,a,z1,z2) - (0.04629809464396494d0)
      print*, "sPxPx a = a", sPxPx
      sPxPx = stoOPxPx(b,a,z1,z2) - (0.06483001236088104d0)
      print*, "sPxPx b != a", sPxPx
      print*, ""
      return
      end
c
c
      subroutine testoverlap()
      implicit none
      real*8 a,b,z1,z2
      real*8 overlapSS, oSS
      real*8 overlapSPz, oSPz
      real*8 overlapPzS, oPzS
      real*8 overlapPzPz, oPzPz
      real*8 overlapPxPx, oPxPx
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      print*, "Overlap"
c
c     overlapSS
c
      oSS = overlapSS(a,b,z1,z2) - (0.03316251188725393d0)
      print*, "oSS a != b", oSS
      oSS = overlapSS(a,a,z1,z2) - (0.02218912782849323d0)
      print*, "oSS a = a", oSS
      oSS = overlapSS(b,a,z1,z2) - (0.03316251188725393d0)
      print*, "oSS b != a", oSS
      print*, ""
c
c     overlapSPz
c
      oSPz = overlapSPz(a,b,z1,z2) - (-0.11449159982005137d0)
      print*, "oSPz a != b", oSPz
      oSPz = overlapSPz(a,a,z1,z2) - (-0.07766194739972629d0)
      print*, "oSPz a = a", oSPz
      oSPz = overlapSPz(b,a,z1,z2) - (-0.09856405008738535d0)
      print*, "oSPz b != a", oSPz
      print*, ""
c
c     overlapPzS
c
      oPzS = overlapPzS(a,b,z1,z2) - (0.09856405008738535d0)
      print*, "oPzS a != b", oPzS
      oPzS = overlapPzS(a,a,z1,z2) - (0.07766194739972629d0)
      print*, "oPzS a = a", oPzS
      oPzS = overlapPzS(b,a,z1,z2) - (0.11449159982005137d0)
      print*, "oPzS b != a", oPzS
      print*, ""
c
c     overlapPzPz
c
      oPzPz = overlapPzPz(a,b,z1,z2) - (-0.20875548850469805d0)
      print*, "oPzPz a != b", oPzPz
      oPzPz = overlapPzPz(a,a,z1,z2) - (-0.17143380952424905d0)
      print*, "oPzPz a = a", oPzPz
      oPzPz = overlapPzPz(b,a,z1,z2) - (-0.20875548850469805d0)
      print*, "oPzPz b != a", oPzPz
      print*, ""
c
c     overlapPxPx
c
      oPxPx = overlapPxPx(a,b,z1,z2) - (0.06463289175879093d0)
      print*, "oPxPx a != b", oPxPx
      oPxPx = overlapPxPx(a,a,z1,z2) - (0.04601964319498459d0)
      print*, "oPxPx a = a", oPxPx
      oPxPx = overlapPxPx(b,a,z1,z2) - (0.06463289175879093d0)
      print*, "oPxPx b != a", oPxPx
      print*, ""
      return
      end
c
c
      subroutine testoverlapflip()
      implicit none
      real*8 a,b,z1,z2
      real*8 overlapSS, oSS
      real*8 overlapSPz, oSPz
      real*8 overlapPzS, oPzS
      real*8 overlapPzPz, oPzPz
      real*8 overlapPxPx, oPxPx
c
c
      a = 3.5
      b = 3.0
      z1 = 1.5
      z2 = -0.5
      print*, "Overlap flip z1 and z2"
c
c     overlapSS
c
      oSS = overlapSS(a,b,z1,z2) - (0.03316251188725393d0)
      print*, "oSS a != b", oSS
      oSS = overlapSS(a,a,z1,z2) - (0.02218912782849323d0)
      print*, "oSS a = a", oSS
      oSS = overlapSS(b,a,z1,z2) - (0.03316251188725393d0)
      print*, "oSS b != a", oSS
      print*, ""
c
c     overlapSPz
c
      oSPz = overlapSPz(a,b,z1,z2) - (0.11449159982005137d0)
      print*, "oSPz a != b", oSPz
      oSPz = overlapSPz(a,a,z1,z2) - (0.07766194739972629d0)
      print*, "oSPz a = a", oSPz
      oSPz = overlapSPz(b,a,z1,z2) - (0.09856405008738535d0)
      print*, "oSPz b != a", oSPz
      print*, ""
c
c     overlapPzS
c
      oPzS = overlapPzS(a,b,z1,z2) - (-0.09856405008738535d0)
      print*, "oPzS a != b", oPzS
      oPzS = overlapPzS(a,a,z1,z2) - (-0.07766194739972629d0)
      print*, "oPzS a = a", oPzS
      oPzS = overlapPzS(b,a,z1,z2) - (-0.11449159982005137d0)
      print*, "oPzS b != a", oPzS
      print*, ""
c
c     overlapPzPz
c
      oPzPz = overlapPzPz(a,b,z1,z2) - (-0.20875548850469805d0)
      print*, "oPzPz a != b", oPzPz
      oPzPz = overlapPzPz(a,a,z1,z2) - (-0.17143380952424905d0)
      print*, "oPzPz a = a", oPzPz
      oPzPz = overlapPzPz(b,a,z1,z2) - (-0.20875548850469805d0)
      print*, "oPzPz b != a", oPzPz
      print*, ""
c
c     overlapPxPx
c
      oPxPx = overlapPxPx(a,b,z1,z2) - (0.06463289175879093d0)
      print*, "oPxPx a != b", oPxPx
      oPxPx = overlapPxPx(a,a,z1,z2) - (0.04601964319498459d0)
      print*, "oPxPx a = a", oPxPx
      oPxPx = overlapPxPx(b,a,z1,z2) - (0.06463289175879093d0)
      print*, "oPxPx b != a", oPxPx
      print*, ""
      return
      end
c
c
      subroutine teststoNE()
      implicit none
      real*8 a,b,z1,z2
      real*8 stoNESaS, NESaS
      real*8 stoNESaPz, NESaPz
      real*8 stoNEPzaS, NEPzaS
      real*8 stoNEPzaPz, NEPzaPz
      real*8 stoNEPxaPx, NEPxaPx
      real*8 stoNESbS, NESbS
      real*8 stoNESbPz, NESbPz
      real*8 stoNEPzbS, NEPzbS
      real*8 stoNEPzbPz, NEPzbPz
      real*8 stoNEPxbPx, NEPxbPx
      real*8 stoNEaSS, NEaSS
      real*8 stoNEaSPz, NEaSPz
      real*8 stoNEaPzPz, NEaPzPz
      real*8 stoNEaPxPx, NEaPxPx
      real*8 stoNESSb, NESSb
      real*8 stoNESPzb, NESPzb
      real*8 stoNEPzPzb, NEPzPzb
      real*8 stoNEPxPxb, NEPxPxb
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      print*, "STO NE"
c
c     NESaS
c
      NESaS = stoNESaS(a,b,z1,z2) - (0.04309604354339577d0)
      print*, "NESaS a != b", NESaS
      NESaS = stoNESaS(a,a,z1,z2) - (0.024414302419061208d0)
      print*, "NESaS a = a", NESaS
      NESaS = stoNESaS(b,a,z1,z2) - (0.031508944525227636d0)
      print*, "NESaS b != a", NESaS
      print*, ""
c
c     NESaPz
c
      NESaPz = stoNESaPz(a,b,z1,z2) - (-0.1905029299553451d0)
      print*, "NESaPz a != b", NESaPz
      NESaPz = stoNESaPz(a,a,z1,z2) - (-0.12087413033035013d0)
      print*, "NESaPz a = a", NESaPz
      NESaPz = stoNESaPz(b,a,z1,z2) - (-0.13536642810481384d0)
      print*, "NESaPz b != a", NESaPz
      print*, ""
c
c     NEPzaS
c
      NEPzaS = stoNEPzaS(a,b,z1,z2) - (0.081221901783265d0)
      print*, "NEPzaS a != b", NEPzaS
      NEPzaS = stoNEPzaS(a,a,z1,z2) - (0.06042771221408536d0)
      print*, "NEPzaS a = a", NEPzaS
      NEPzaS = stoNEPzaS(b,a,z1,z2) - (0.0807931008266918d0)
      print*, "NEPzaS b != a", NEPzaS
      print*, ""
c
c     NEPzaPz
c
      NEPzaPz = stoNEPzaPz(a,b,z1,z2) - (-0.2134639646520742d0)
      print*, "NEPzaPz a != b", NEPzaPz
      NEPzaPz = stoNEPzaPz(a,a,z1,z2) - (-0.16899389944110002d0)
      print*, "NEPzaPz a = a", NEPzaPz
      NEPzaPz = stoNEPzaPz(b,a,z1,z2) - (-0.19366613504699962d0)
      print*, "NEPzaPz b != a", NEPzaPz
      print*, ""
c
c     NEPxaPx
c
      NEPxaPx = stoNEPxaPx(a,b,z1,z2) - (0.057328228552561d0)
      print*, "NEPxaPx a != b", NEPxaPx
      NEPxaPx = stoNEPxaPx(a,a,z1,z2) - (0.039106906279815196d0)
      print*, "NEPxaPx a = a", NEPxaPx
      NEPxaPx = stoNEPxaPx(b,a,z1,z2) - (0.04949516342206589d0)
      print*, "NEPxaPx b != a", NEPxaPx
      print*, ""
c
c     NESbS
c
      NESbS = stoNESbS(a,b,z1,z2) - (0.03150894452522763d0)
      print*, "NESbS a != b", NESbS
      NESbS = stoNESbS(a,a,z1,z2) - (0.02441430241906121d0)
      print*, "NESbS a = a", NESbS
      NESbS = stoNESbS(b,a,z1,z2) - (0.043096043543395766d0)
      print*, "NESbS b != a", NESbS
      print*, ""
c
c     NESbPz
c
      NESbPz = stoNESbPz(a,b,z1,z2) - (-0.08079310082669183d0)
      print*, "NESbPz a != b", NESbPz
      NESbPz = stoNESbPz(a,a,z1,z2) - (-0.06042771221408537d0)
      print*, "NESbPz a = a", NESbPz
      NESbPz = stoNESbPz(b,a,z1,z2) - (-0.08122190178326498d0)
      print*, "NESbPz b != a", NESbPz
      print*, ""
c
c     NEPzbS
c
      NEPzbS = stoNEPzbS(a,b,z1,z2) - (0.13536642810481386d0)
      print*, "NEPzbS a != b", NEPzbS
      NEPzbS = stoNEPzbS(a,a,z1,z2) - (0.12087413033035013d0)
      print*, "NEPzbS a = a", NEPzbS
      NEPzbS = stoNEPzbS(b,a,z1,z2) - (0.1905029299553451d0)
      print*, "NEPzbS b != a", NEPzbS
      print*, ""
c
c     NEPzbPz
c
      NEPzbPz = stoNEPzbPz(a,b,z1,z2) - (-0.19366613504699964d0)
      print*, "NEPzbPz a != b", NEPzbPz
      NEPzbPz = stoNEPzbPz(a,a,z1,z2) - (-0.1689938994411d0)
      print*, "NEPzbPz a = a", NEPzbPz
      NEPzbPz = stoNEPzbPz(b,a,z1,z2) - (-0.21346396465207423d0)
      print*, "NEPzbPz b != a", NEPzbPz
      print*, ""
c
c     NEPxbPx
c
      NEPxbPx = stoNEPxbPx(a,b,z1,z2) - (0.049495163422065884d0)
      print*, "NEPxbPx a != b", NEPxbPx
      NEPxbPx = stoNEPxbPx(a,a,z1,z2) - (0.039106906279815196d0)
      print*, "NEPxbPx a = a", NEPxbPx
      NEPxbPx = stoNEPxbPx(b,a,z1,z2) - (0.05732822855256101d0)
      print*, "NEPxbPx b != a", NEPxbPx
      print*, ""
c
c     NEaSS
c
      NEaSS = stoNEaSS(a,z1,z2) - (0.4999985446647136d0)
      print*, "NEaSS a", NEaSS
      NEaSS = stoNEaSS(b,z1,z2) - (0.49998411943117493d0)
      print*, "NEaSS b", NEaSS
      print*, ""
c
c     NEaSPz
c
      NEaSPz = stoNEaSPz(a,z1,z2) - (-0.07138902598070776d0)
      print*, "NEaSPz a", NEaSPz
      NEaSPz = stoNEaSPz(b,z1,z2) - (-0.08317633408219458d0)
      print*, "NEaSPz b", NEaSPz
      print*, ""
c
c     NEaPzPz
c
      NEaPzPz = stoNEaPzPz(a,z1,z2) - (0.5303903982086646d0)
      print*, "NEaPzPz a", NEaPzPz
      NEaPzPz = stoNEaPzPz(b,z1,z2) - (0.5404919821032068d0)
      print*, "NEaPzPz b", NEaPzPz
      print*, ""
c
c     NEaPxPx
c
      NEaPxPx = stoNEaPxPx(a,z1,z2) - (0.48468907742221023d0)
      print*, "NEaPxPx a", NEaPxPx
      NEaPxPx = stoNEaPxPx(b,z1,z2) - (0.47917468722489526d0)
      print*, "NEaPxPx b", NEaPxPx
      print*, ""
c
c     NESSb
c
      NESSb = stoNESSb(a,z1,z2) - (0.4999985446647136d0)
      print*, "NESSb a", NESSb
      NESSb = stoNESSb(b,z1,z2) - (0.499984119431175d0)
      print*, "NESSb b", NESSb
      print*, ""
c
c     NESPzb
c
      NESPzb = stoNESPzb(a,z1,z2) - (0.07138902598070769d0)
      print*, "NESPzb a", NESPzb
      NESPzb = stoNESPzb(b,z1,z2) - (0.08317633408219449d0)
      print*, "NESPzb b", NESPzb
      print*, ""
c
c     NEPzPzb
c
      NEPzPzb = stoNEPzPzb(a,z1,z2) - (0.5303903982086646d0)
      print*, "NEPzPzb a", NEPzPzb
      NEPzPzb = stoNEPzPzb(b,z1,z2) - (0.5404919821032068d0)
      print*, "NEPzPzb b", NEPzPzb
      print*, ""
c
c     NEPxPxb
c
      NEPxPxb = stoNEPxPxb(a,z1,z2) - (0.48468907742221023d0)
      print*, "NEPxPxb a", NEPxPxb
      NEPxPxb = stoNEPxPxb(b,z1,z2) - (0.47917468722489526d0)
      print*, "NEPxPxb b", NEPxPxb
      print*, ""
      return
      end
c
c
      subroutine testNE()
      implicit none
      real*8 a,b,z1,z2
      real*8 NESaS, nNESaS
      real*8 NESaPz, nNESaPz
      real*8 NEPzaS, nNEPzaS
      real*8 NEPzaPz, nNEPzaPz
      real*8 NEPxaPx, nNEPxaPx
      real*8 NESbS, nNESbS
      real*8 NESbPz, nNESbPz
      real*8 NEPzbS, nNEPzbS
      real*8 NEPzbPz, nNEPzbPz
      real*8 NEPxbPx, nNEPxbPx
      real*8 NEaSS, nNEaSS
      real*8 NEaSPz, nNEaSPz
      real*8 NEaPzPz, nNEaPzPz
      real*8 NEaPxPx, nNEaPxPx
      real*8 NESSb, nNESSb
      real*8 NESPzb, nNESPzb
      real*8 NEPzPzb, nNEPzPzb
      real*8 NEPxPxb, nNEPxPxb
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      print*, "Nuclear-Electron"
c
c     nNESaS
c
      nNESaS = NESaS(a,b,z1,z2) - (0.043232760593498634d0)
      print*, "nNESaS a != b", nNESaS
      nNESaS = NESaS(a,a,z1,z2) - (0.025532695035526454d0)
      print*, "nNESaS a = a", nNESaS
      nNESaS = NESaS(b,a,z1,z2) - (0.03247519342015251d0)
      print*, "nNESaS b != a", nNESaS
      print*, ""
c
c     nNESaPz
c
      nNESaPz = NESaPz(a,b,z1,z2) - (-0.1906499804413834d0)
      print*, "nNESaPz a != b", nNESaPz
      nNESaPz = NESaPz(a,a,z1,z2) - (-0.11915257683245681d0)
      print*, "nNESaPz a = a", nNESaPz
      nNESaPz = NESaPz(b,a,z1,z2) - (-0.1337546158060454d0)
      print*, "nNESaPz b != a", nNESaPz
      print*, ""
c
c     nNEPzaS
c
      nNEPzaS = NEPzaS(a,b,z1,z2) - (0.08020434697287639d0)
      print*, "nNEPzaS a != b", nNEPzaS
      nNEPzaS = NEPzaS(a,a,z1,z2) - (0.059576288416228404d0)
      print*, "nNEPzaS a = a", nNEPzaS
      nNEPzaS = NEPzaS(b,a,z1,z2) - (0.08020434697287639d0)
      print*, "nNEPzaS b != a", nNEPzaS
      print*, ""
c
c     nNEPzaPz
c
      nNEPzaPz = NEPzaPz(a,b,z1,z2) - (-0.2150998208606128d0)
      print*, "nNEPzaPz a != b", nNEPzaPz
      nNEPzaPz = NEPzaPz(a,a,z1,z2) - (-0.16968603575693622d0)
      print*, "nNEPzaPz a = a", nNEPzaPz
      nNEPzaPz = NEPzaPz(b,a,z1,z2) - (-0.19441184620070445d0)
      print*, "nNEPzaPz b != a", nNEPzaPz
      print*, ""
c
c     nNEPxaPx
c
      nNEPxaPx = NEPxaPx(a,b,z1,z2) - (0.057245799910025694d0)
      print*, "nNEPxaPx a != b", nNEPxaPx
      nNEPxaPx = NEPxaPx(a,a,z1,z2) - (0.038830973699863144d0)
      print*, "nNEPxaPx a = a", nNEPxaPx
      nNEPxaPx = NEPxaPx(b,a,z1,z2) - (0.049282025043692675d0)
      print*, "nNEPxaPx b != a", nNEPxaPx
      print*, ""
c
c     nNESbS
c
      nNESbS = NESbS(a,b,z1,z2) - (0.03247519342015251d0)
      print*, "nNESbS a != b", nNESbS
      nNESbS = NESbS(a,a,z1,z2) - (0.025532695035526454d0)
      print*, "nNESbS a = a", nNESbS
      nNESbS = NESbS(b,a,z1,z2) - (0.043232760593498634d0)
      print*, "nNESbS b != a", nNESbS
      print*, ""
c
c     nNESbPz
c
      nNESbPz = NESbPz(a,b,z1,z2) - (-0.08020434697287639d0)
      print*, "nNESbPz a != b", nNESbPz
      nNESbPz = NESbPz(a,a,z1,z2) - (-0.059576288416228404d0)
      print*, "nNESbPz a = a", nNESbPz
      nNESbPz = NESbPz(b,a,z1,z2) - (-0.08020434697287639d0)
      print*, "nNESbPz b != a", nNESbPz
      print*, ""
c
c     nNEPzbS
c
      nNEPzbS = NEPzbS(a,b,z1,z2) - (0.1337546158060454d0)
      print*, "nNEPzbS a != b", nNEPzbS
      nNEPzbS = NEPzbS(a,a,z1,z2) - (0.11915257683245681d0)
      print*, "nNEPzbS a = a", nNEPzbS
      nNEPzbS = NEPzbS(b,a,z1,z2) - (0.1906499804413834d0)
      print*, "nNEPzbS b != a", nNEPzbS
      print*, ""
c
c     nNEPzbPz
c
      nNEPzbPz = NEPzbPz(a,b,z1,z2) - (-0.19441184620070445d0)
      print*, "nNEPzbPz a != b", nNEPzbPz
      nNEPzbPz = NEPzbPz(a,a,z1,z2) - (-0.16968603575693622d0)
      print*, "nNEPzbPz a = a", nNEPzbPz
      nNEPzbPz = NEPzbPz(b,a,z1,z2) - (-0.2150998208606128d0)
      print*, "nNEPzbPz b != a", nNEPzbPz
      print*, ""
c
c     nNEPxbPx
c
      nNEPxbPx = NEPxbPx(a,b,z1,z2) - (0.049282025043692675d0)
      print*, "nNEPxbPx a != b", nNEPxbPx
      nNEPxbPx = NEPxbPx(a,a,z1,z2) - (0.038830973699863144d0)
      print*, "nNEPxbPx a = a", nNEPxbPx
      nNEPxbPx = NEPxbPx(b,a,z1,z2) - (0.057245799910025694d0)
      print*, "nNEPxbPx b != a", nNEPxbPx
      print*, ""
c
c     nNEaSS
c
      nNEaSS = NEaSS(a,z1,z2) - (0.4999966738851236d0)
      print*, "nNEaSS a", nNEaSS
      nNEaSS = NEaSS(b,z1,z2) - (0.49997849525676336d0)
      print*, "nNEaSS b", nNEaSS
      print*, ""
c
c     nNEaSPz
c
      nNEaSPz = NEaSPz(a,z1,z2) - (-0.07140148735029206d0)
      print*, "nNEaSPz a", nNEaSPz
      nNEaSPz = NEaSPz(b,z1,z2) - (-0.08317921600680402d0)
      print*, "nNEaSPz b", nNEaSPz
      print*, ""
c
c     nNEaPzPz
c
      nNEaPzPz = NEaPzPz(a,z1,z2) - (0.5303892212590073d0)
      print*, "nNEaPzPz a", nNEaPzPz
      nNEaPzPz = NEaPzPz(b,z1,z2) - (0.5405450919033362d0)
      print*, "nNEaPzPz b", nNEaPzPz
      print*, ""
c
c     nNEaPxPx
c
      nNEaPxPx = NEaPxPx(a,z1,z2) - (0.4846963551672039d0)
      print*, "nNEaPxPx a", nNEaPxPx
      nNEaPxPx = NEaPxPx(b,z1,z2) - (0.4791836912550624d0)
      print*, "nNEaPxPx b", nNEaPxPx
      print*, ""
c
c     nNESSb
c
      nNESSb = NESSb(a,z1,z2) - (0.4999966738851236d0)
      print*, "nNESSb a", nNESSb
      nNESSb = NESSb(b,z1,z2) - (0.49997849525676336d0)
      print*, "nNESSb b", nNESSb
      print*, ""
c
c     nNESPzb
c
      nNESPzb = NESPzb(a,z1,z2) - (0.07140148735029206d0)
      print*, "nNESPzb a", nNESPzb
      nNESPzb = NESPzb(b,z1,z2) - (0.08317921600680402d0)
      print*, "nNESPzb b", nNESPzb
      print*, ""
c
c     nNEPzPzb
c
      nNEPzPzb = NEPzPzb(a,z1,z2) - (0.5303892212590073d0)
      print*, "nNEPzPzb a", nNEPzPzb
      nNEPzPzb = NEPzPzb(b,z1,z2) - (0.5405450919033362d0)
      print*, "nNEPzPzb b", nNEPzPzb
      print*, ""
c
c     nNEPxPxb
c
      nNEPxPxb = NEPxPxb(a,z1,z2) - (0.4846963551672039d0)
      print*, "nNEPxPxb a", nNEPxPxb
      nNEPxPxb = NEPxPxb(b,z1,z2) - (0.4791836912550624d0)
      print*, "nNEPxPxb b", nNEPxPxb
      print*, ""
      return
      end
c
c
      subroutine testNEflip()
      implicit none
      real*8 a,b,z1,z2
      real*8 NESaS, nNESaS
      real*8 NESaPz, nNESaPz
      real*8 NEPzaS, nNEPzaS
      real*8 NEPzaPz, nNEPzaPz
      real*8 NEPxaPx, nNEPxaPx
      real*8 NESbS, nNESbS
      real*8 NESbPz, nNESbPz
      real*8 NEPzbS, nNEPzbS
      real*8 NEPzbPz, nNEPzbPz
      real*8 NEPxbPx, nNEPxbPx
      real*8 NEaSS, nNEaSS
      real*8 NEaSPz, nNEaSPz
      real*8 NEaPzPz, nNEaPzPz
      real*8 NEaPxPx, nNEaPxPx
      real*8 NESSb, nNESSb
      real*8 NESPzb, nNESPzb
      real*8 NEPzPzb, nNEPzPzb
      real*8 NEPxPxb, nNEPxPxb
c
c
      a = 3.5
      b = 3.0
      z1 = 1.5
      z2 = -0.5
      print*, "Nuclear-Electron flip z1 and z2"
c
c     nNESaS
c
      nNESaS = NESaS(a,b,z1,z2) - (0.043232760593498634d0)
      print*, "nNESaS a != b", nNESaS
      nNESaS = NESaS(a,a,z1,z2) - (0.025532695035526454d0)
      print*, "nNESaS a = a", nNESaS
      nNESaS = NESaS(b,a,z1,z2) - (0.03247519342015251d0)
      print*, "nNESaS b != a", nNESaS
      print*, ""
c
c     nNESaPz
c
      nNESaPz = NESaPz(a,b,z1,z2) - (0.1906499804413834d0)
      print*, "nNESaPz a != b", nNESaPz
      nNESaPz = NESaPz(a,a,z1,z2) - (0.11915257683245681d0)
      print*, "nNESaPz a = a", nNESaPz
      nNESaPz = NESaPz(b,a,z1,z2) - (0.1337546158060454d0)
      print*, "nNESaPz b != a", nNESaPz
      print*, ""
c
c     nNEPzaS
c
      nNEPzaS = NEPzaS(a,b,z1,z2) - (-0.08020434697287639d0)
      print*, "nNEPzaS a != b", nNEPzaS
      nNEPzaS = NEPzaS(a,a,z1,z2) - (-0.059576288416228404d0)
      print*, "nNEPzaS a = a", nNEPzaS
      nNEPzaS = NEPzaS(b,a,z1,z2) - (-0.08020434697287639d0)
      print*, "nNEPzaS b != a", nNEPzaS
      print*, ""
c
c     nNEPzaPz
c
      nNEPzaPz = NEPzaPz(a,b,z1,z2) - (-0.2150998208606128d0)
      print*, "nNEPzaPz a != b", nNEPzaPz
      nNEPzaPz = NEPzaPz(a,a,z1,z2) - (-0.16968603575693622d0)
      print*, "nNEPzaPz a = a", nNEPzaPz
      nNEPzaPz = NEPzaPz(b,a,z1,z2) - (-0.19441184620070445d0)
      print*, "nNEPzaPz b != a", nNEPzaPz
      print*, ""
c
c     nNEPxaPx
c
      nNEPxaPx = NEPxaPx(a,b,z1,z2) - (0.057245799910025694d0)
      print*, "nNEPxaPx a != b", nNEPxaPx
      nNEPxaPx = NEPxaPx(a,a,z1,z2) - (0.038830973699863144d0)
      print*, "nNEPxaPx a = a", nNEPxaPx
      nNEPxaPx = NEPxaPx(b,a,z1,z2) - (0.049282025043692675d0)
      print*, "nNEPxaPx b != a", nNEPxaPx
      print*, ""
c
c     nNESbS
c
      nNESbS = NESbS(a,b,z1,z2) - (0.03247519342015251d0)
      print*, "nNESbS a != b", nNESbS
      nNESbS = NESbS(a,a,z1,z2) - (0.025532695035526454d0)
      print*, "nNESbS a = a", nNESbS
      nNESbS = NESbS(b,a,z1,z2) - (0.043232760593498634d0)
      print*, "nNESbS b != a", nNESbS
      print*, ""
c
c     nNESbPz
c
      nNESbPz = NESbPz(a,b,z1,z2) - (0.08020434697287639d0)
      print*, "nNESbPz a != b", nNESbPz
      nNESbPz = NESbPz(a,a,z1,z2) - (0.059576288416228404d0)
      print*, "nNESbPz a = a", nNESbPz
      nNESbPz = NESbPz(b,a,z1,z2) - (0.08020434697287639d0)
      print*, "nNESbPz b != a", nNESbPz
      print*, ""
c
c     nNEPzbS
c
      nNEPzbS = NEPzbS(a,b,z1,z2) - (-0.1337546158060454d0)
      print*, "nNEPzbS a != b", nNEPzbS
      nNEPzbS = NEPzbS(a,a,z1,z2) - (-0.11915257683245681d0)
      print*, "nNEPzbS a = a", nNEPzbS
      nNEPzbS = NEPzbS(b,a,z1,z2) - (-0.1906499804413834d0)
      print*, "nNEPzbS b != a", nNEPzbS
      print*, ""
c
c     nNEPzbPz
c
      nNEPzbPz = NEPzbPz(a,b,z1,z2) - (-0.19441184620070445d0)
      print*, "nNEPzbPz a != b", nNEPzbPz
      nNEPzbPz = NEPzbPz(a,a,z1,z2) - (-0.16968603575693622d0)
      print*, "nNEPzbPz a = a", nNEPzbPz
      nNEPzbPz = NEPzbPz(b,a,z1,z2) - (-0.2150998208606128d0)
      print*, "nNEPzbPz b != a", nNEPzbPz
      print*, ""
c
c     nNEPxbPx
c
      nNEPxbPx = NEPxbPx(a,b,z1,z2) - (0.049282025043692675d0)
      print*, "nNEPxbPx a != b", nNEPxbPx
      nNEPxbPx = NEPxbPx(a,a,z1,z2) - (0.038830973699863144d0)
      print*, "nNEPxbPx a = a", nNEPxbPx
      nNEPxbPx = NEPxbPx(b,a,z1,z2) - (0.057245799910025694d0)
      print*, "nNEPxbPx b != a", nNEPxbPx
      print*, ""
c
c     nNEaSS
c
      nNEaSS = NEaSS(a,z1,z2) - (0.4999966738851236d0)
      print*, "nNEaSS a", nNEaSS
      nNEaSS = NEaSS(b,z1,z2) - (0.49997849525676336d0)
      print*, "nNEaSS b", nNEaSS
      print*, ""
c
c     nNEaSPz
c
      nNEaSPz = NEaSPz(a,z1,z2) - (0.07140148735029206d0)
      print*, "nNEaSPz a", nNEaSPz
      nNEaSPz = NEaSPz(b,z1,z2) - (0.08317921600680402d0)
      print*, "nNEaSPz b", nNEaSPz
      print*, ""
c
c     nNEaPzPz
c
      nNEaPzPz = NEaPzPz(a,z1,z2) - (0.5303892212590073d0)
      print*, "nNEaPzPz a", nNEaPzPz
      nNEaPzPz = NEaPzPz(b,z1,z2) - (0.5405450919033362d0)
      print*, "nNEaPzPz b", nNEaPzPz
      print*, ""
c
c     nNEaPxPx
c
      nNEaPxPx = NEaPxPx(a,z1,z2) - (0.4846963551672039d0)
      print*, "nNEaPxPx a", nNEaPxPx
      nNEaPxPx = NEaPxPx(b,z1,z2) - (0.4791836912550624d0)
      print*, "nNEaPxPx b", nNEaPxPx
      print*, ""
c
c     nNESSb
c
      nNESSb = NESSb(a,z1,z2) - (0.4999966738851236d0)
      print*, "nNESSb a", nNESSb
      nNESSb = NESSb(b,z1,z2) - (0.49997849525676336d0)
      print*, "nNESSb b", nNESSb
      print*, ""
c
c     nNESPzb
c
      nNESPzb = NESPzb(a,z1,z2) - (-0.07140148735029206d0)
      print*, "nNESPzb a", nNESPzb
      nNESPzb = NESPzb(b,z1,z2) - (-0.08317921600680402d0)
      print*, "nNESPzb b", nNESPzb
      print*, ""
c
c     nNEPzPzb
c
      nNEPzPzb = NEPzPzb(a,z1,z2) - (0.5303892212590073d0)
      print*, "nNEPzPzb a", nNEPzPzb
      nNEPzPzb = NEPzPzb(b,z1,z2) - (0.5405450919033362d0)
      print*, "nNEPzPzb b", nNEPzPzb
      print*, ""
c
c     nNEPxPxb
c
      nNEPxPxb = NEPxPxb(a,z1,z2) - (0.4846963551672039d0)
      print*, "nNEPxPxb a", nNEPxPxb
      nNEPxPxb = NEPxPxb(b,z1,z2) - (0.4791836912550624d0)
      print*, "nNEPxPxb b", nNEPxPxb
      print*, ""
      return
      end
c
c
      subroutine testoverlapsum()
      implicit none
      real*8 a,b,z1,z2
      real*8 overlap,overlapTotal
      real*8 coeff1(4),coeff2(4)
      logical exact
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      coeff1(1) = 0.5345224838248488d0
      coeff1(2) = 0.25210259d0
      coeff1(3) = 0.31686303d0
      coeff1(4) = 0.74184083d0
      coeff2(1) = 0.7071067811865475d0
      coeff2(2) = 0.22280876d0
      coeff2(3) = 0.14315011d0
      coeff2(4) = 0.65564038d0
      z1 = -0.5
      z2 = 1.5
      print*, "Overlap Sum"
      exact = .true.
      overlap = overlapTotal(coeff1, coeff2, a, b, a, b, z1, z2, exact)
      print*, "Exact: ", overlap - (-0.0708595303952499d0)
      exact = .false.
      overlap = overlapTotal(coeff1, coeff2, a, b, a, b, z1, z2, exact)
      print*, "STO: ", overlap - (-0.07031596688785882d0)
      exact = .true.
      overlap = overlapTotal(coeff2, coeff1, b, a, b, a, z1, z2, exact)
      print*, "Exact: ", overlap - (-0.09401698262812486d0)
      exact = .false.
      overlap = overlapTotal(coeff2, coeff1, b, a, b, a, z1, z2, exact)
      print*, "STO: ", overlap - (-0.09418355317621996d0)
      print*, "Flip z1 and z2"
      z1 = 1.5
      z2 = -0.5
      print*, "Overlap Sum"
      exact = .true.
      overlap = overlapTotal(coeff1, coeff2, a, b, a, b, z1, z2, exact)
      print*, "Exact: ", overlap - (-0.09401698262812484d0)
      exact = .false.
      overlap = overlapTotal(coeff1, coeff2, a, b, a, b, z1, z2, exact)
      print*, "STO: ", overlap - (-0.09418355317621999d0)
      exact = .true.
      overlap = overlapTotal(coeff2, coeff1, b, a, b, a, z1, z2, exact)
      print*, "Exact: ", overlap - (-0.0708595303952499d0)
      exact = .false.
      overlap = overlapTotal(coeff2, coeff1, b, a, b, a, z1, z2, exact)
      print*, "STO: ", overlap - (-0.07031596688785881d0)
      return
      end
c
c
      subroutine testNEsumbAb()
      implicit none
      real*8 a,b,z1,z2
      real*8 NEbAb,NEbAbTotal
      real*8 NEaBa,NEaBaTotal
      real*8 NEaAb,NEaAbTotal
      real*8 NEaBb,NEaBbTotal
      real*8 coeff1(4),coeff2(4)
      logical exact
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      coeff1(1) = 0.5345224838248488d0
      coeff1(2) = 0.25210259d0
      coeff1(3) = 0.31686303d0
      coeff1(4) = 0.74184083d0
      coeff2(1) = 0.7071067811865475d0
      coeff2(2) = 0.22280876d0
      coeff2(3) = 0.14315011d0
      coeff2(4) = 0.65564038d0
      print*, "NE bAb Sum"
      z1 = -0.5
      z2 = 1.5
      exact = .true.
      NEbAb = NEbAbTotal(coeff2, b, z1, z2, exact)
      print*, "Exact: ", NEbAb - (0.4388331120603912d0)
      exact = .false.
      NEbAb = NEbAbTotal(coeff2, b, z1, z2, exact)
      print*, "STO: ", NEbAb - (0.4388151348014041d0)
      print*, "Flip z1 and z2"
      z1 = 1.5
      z2 = -0.5
      exact = .true.
      NEbAb = NEbAbTotal(coeff2, b, z1, z2, exact)
      print*, "Exact: ", NEbAb - (0.5930832316796388d0)
      exact = .false.
      NEbAb = NEbAbTotal(coeff2, b, z1, z2, exact)
      print*, "STO: ", NEbAb - (0.5930599100902159d0)
      print*, "NE aBa Sum"
      z1 = -0.5
      z2 = 1.5
      exact = .true.
      NEaBa = NEaBaTotal(coeff1, a, z1, z2, exact)
      print*, "Exact: ", NEaBa - (0.5708396789695409d0)
      exact = .false.
      NEaBa = NEaBaTotal(coeff1, a, z1, z2, exact)
      print*, "STO: ", NEaBa - (0.570829785313584d0)
      print*, "Flip z1 and z2"
      z1 = 1.5
      z2 = -0.5
      exact = .true.
      NEaBa = NEaBaTotal(coeff1, a, z1, z2, exact)
      print*, "Exact: ", NEaBa - (0.45758817961760334d0)
      exact = .false.
      NEaBa = NEaBaTotal(coeff1, a, z1, z2, exact)
      print*, "STO: ", NEaBa - (0.45759805121923114d0)
      print*, "NE aAb Sum"
      z1 = -0.5
      z2 = 1.5
      exact = .true.
      NEaAb = NEaAbTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "Exact: ", NEaAb - (-0.10720993185653073d0)
      exact = .false.
      NEaAb = NEaAbTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "STO: ", NEaAb - (-0.10587228428335016d0)
      exact = .true.
      NEaAb = NEaAbTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "Exact: ", NEaAb - (-0.11933457036134326d0)
      exact = .false.
      NEaAb = NEaAbTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "STO: ", NEaAb - (-0.1199546005679584d0)
      print*, "Flip z1 and z2"
      z1 = 1.5
      z2 = -0.5
      exact = .true.
      NEaAb = NEaAbTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "Exact: ", NEaAb - (-0.057725729213826155d0)
      exact = .false.
      NEaAb = NEaAbTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "STO: ", NEaAb - (-0.05755868914645285d0)
      exact = .true.
      NEaAb = NEaAbTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "Exact: ", NEaAb - (-0.035225698080198105d0)
      exact = .false.
      NEaAb = NEaAbTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "STO: ", NEaAb - (-0.0345674044748482d0)
      print*, "NE aBb Sum"
      z1 = -0.5
      z2 = 1.5
      exact = .true.
      NEaBb = NEaBbTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "Exact: ", NEaBb - (-0.0352256980801981d0)
      exact = .false.
      NEaBb = NEaBbTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "STO: ", NEaBb - (-0.03456740447484821d0)
      exact = .true.
      NEaBb = NEaBbTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "Exact: ", NEaBb - (-0.057725729213826155d0)
      exact = .false.
      NEaBb = NEaBbTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "STO: ", NEaBb - (-0.05755868914645283d0)
      print*, "Flip z1 and z2"
      z1 = 1.5
      z2 = -0.5
      exact = .true.
      NEaBb = NEaBbTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "Exact: ", NEaBb - (-0.11933457036134326d0)
      exact = .false.
      NEaBb = NEaBbTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "STO: ", NEaBb - (-0.11995460056795842d0)
      exact = .true.
      NEaBb = NEaBbTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "Exact: ", NEaBb - (-0.10720993185653074d0)
      exact = .false.
      NEaBb = NEaBbTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "STO: ", NEaBb - (-0.1058722842833502d0)
      return
      end
c
c
      subroutine testcoulombsum()
      implicit none
      real*8 a,b,z1,z2
      real*8 coulomb,coulombTotal
      real*8 coeff1(4),coeff2(4)
      logical exact
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      coeff1(1) = 0.5345224838248488d0
      coeff1(2) = 0.25210259d0
      coeff1(3) = 0.31686303d0
      coeff1(4) = 0.74184083d0
      coeff2(1) = 0.7071067811865475d0
      coeff2(2) = 0.22280876d0
      coeff2(3) = 0.14315011d0
      coeff2(4) = 0.65564038d0
      z1 = -0.5
      z2 = 1.5
      print*, "Coulomb Sum"
      exact = .true.
      coulomb = coulombTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "Exact: ", coulomb - (0.49174301620966243d0)
      exact = .false.
      coulomb = coulombTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "STO: ", coulomb - (0.4917333207032612d0)
      exact = .true.
      coulomb = coulombTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "Exact: ", coulomb - (0.5334252923825443d0)
      exact = .false.
      coulomb = coulombTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "STO: ", coulomb - (0.5334069433534098d0)
      print*, "Flip z1 and z2"
      z1 = 1.5
      z2 = -0.5
      print*, "Coulomb Sum"
      exact = .true.
      coulomb = coulombTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "Exact: ", coulomb - (0.5334252923825444d0)
      exact = .false.
      coulomb = coulombTotal(coeff1, coeff2, a, b, z1, z2, exact)
      print*, "STO: ", coulomb - (0.5334069433534099d0)
      exact = .true.
      coulomb = coulombTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "Exact: ", coulomb - (0.49174301620966254d0)
      exact = .false.
      coulomb = coulombTotal(coeff2, coeff1, b, a, z1, z2, exact)
      print*, "STO: ", coulomb - (0.4917333207032608d0)
      return
      end
c
c
      subroutine testexchangesum()
      implicit none
      real*8 a,b,z1,z2
      real*8 exchange,exchangeTotal
      real*8 coeff1(4),coeff2(4)
c
c
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      coeff1(1) = 0.5345224838248488d0
      coeff1(2) = 0.25210259d0
      coeff1(3) = 0.31686303d0
      coeff1(4) = 0.74184083d0
      coeff2(1) = 0.7071067811865475d0
      coeff2(2) = 0.22280876d0
      coeff2(3) = 0.14315011d0
      coeff2(4) = 0.65564038d0
      z1 = -0.5
      z2 = 1.5
      print*, "Exchange Sum"
      exchange = exchangeTotal(coeff1, coeff2, a, b, z1, z2)
      print*, "STO: ", exchange - (0.008989681212052883d0)
      exchange = exchangeTotal(coeff2, coeff1, b, a, z1, z2)
      print*, "STO: ", exchange - (0.014931157864616185d0)
      print*, "Flip z1 and z2"
      z1 = 1.5
      z2 = -0.5
      print*, "Exchange Sum"
      exchange = exchangeTotal(coeff1, coeff2, a, b, z1, z2)
      print*, "STO: ", exchange - (0.01493115786461617d0)
      exchange = exchangeTotal(coeff2, coeff1, b, a, z1, z2)
      print*, "STO: ", exchange - (0.008989681212052885d0)
      return
      end