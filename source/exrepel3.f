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
      real*8 e,eterm
      real*8 fgrp,taper
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 r,r2
      real*8 zxri,zxrk
      real*8 vali,valk,valik
      real*8 dmpi,dmpk
      real*8 intS,intS2
      real*8 intK,intJ
      real*8 intaAb,intbBa
      real*8 intbAb,intaBa
      real*8 termS0,termS1
      real*8 termS2
      real*8, allocatable :: rscale(:)
      logical proceed,usei
      logical muti,mutk,mutik
      logical header,huge
      character*6 mode
c
c
c     zero out the repulsion energy and partitioning terms
c
      call testcoulomb()
      call testcoulombsto()
      call testexchange()
      ner = 0
      er = 0.0d0
      do i = 1, n
         aer(i) = 0.0d0
      end do
      if (nrep .eq. 0)  return
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
                  valik = vali * valk
                  call computeInt (dmpi,dmpk,r,r2,xi,yi,zi,xk,yk,zk,
     &                       intS,intK,intaAb,intbBa,intbAb,intaBa,intJ)
                  intS2 = intS * intS
                  termS0 = -valik * intK
                  termS1 = -intS * (zxri * valk * intaAb
     &                            + zxrk * vali * intbBa)
                  termS2 = intS2 * (zxri * valk * intbAb
     &                            + zxrk * vali * intaBa
     &                            + valik * intJ)
                  e = hartree * (termS0 + termS1 + termS2)
     &                          / (1.0d0 - intS2)
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine computeInt  --  exchange repulsion integrals  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "computeIntegrals" computes the six integrals needed for exchange
c     repulsion calculation
c
c
      subroutine computeInt (dmpi,dmpk,r,r2,xi,yi,zi,xk,yk,zk,
     &                       intS,intK,intaAb,intbBa,intbAb,intaBa,intJ)
      use xrepel
      implicit none
      real*8 dmpi,dmpk
      real*8 r,r2
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 intS
      real*8 intK,intJ
      real*8 intaAb,intbBa
      real*8 intbAb,intaBa
      real*8 rho,rho2,rho3
      real*8 expi,expk
      real*8 expi2,expk2
      real*8 expik,expik2
      real*8 eps,diff
      real*8 dmpik
      real*8 dmpi2,dmpk2
      real*8 rhoi,rhok
      real*8 tau,tau2
      real*8 kappa
      real*8 kappap,kappam
      real*8 pre1,pre2
      real*8 term1,term2
      real*8 SSSS
c
c
      if (xreptyp .eq. 'DEFAULT') then
c
c     compute tolerance value for damping exponents
c
         eps = 0.001d0
         diff = abs(dmpi-dmpk)
c
c     treat the case where alpha damping exponents are equal
c
         if (diff .lt. eps) then
            rho = dmpi * r
            rho2 = rho * rho
            rho3 = rho2 * rho
            expi = exp(-rho)
            expi2 = expi * expi
            intS = (1.0d0 + rho + 1.0d0 / 3.0d0 * rho2) * expi
            intaAb = dmpi * (1.0d0 + rho) * expi
            intbBa = intaAb
            intbAb = (1.0d0 - (1.0d0 + rho) * expi2) / r
            intaBa = intbAb
            intJ = (1.0d0 - (1.0d0 + 11.0d0/8.0d0 * rho
     &           + 3.0d0/4.0d0 * rho2 + 1.0d0/6.0d0 * rho3) * expi2) / r
            intK = SSSS(dmpi,dmpi,0.0d0,r,.true.)
c
c     treat the case where alpha damping exponents are not equal
c
         else
            dmpik = (dmpi + dmpk) / 2.0d0
            dmpi2 = dmpi * dmpi
            dmpk2 = dmpk * dmpk
            rhoi = dmpi * r
            rhok = dmpk * r
            rho = dmpik * r
            expi = exp(-rhoi)
            expk = exp(-rhok)
            expi2 = expi * expi
            expk2 = expk * expk
            tau = (dmpi - dmpk) / (dmpi + dmpk)
            tau2 = tau * tau
            kappa = (dmpi2 + dmpk2) / (dmpi2 - dmpk2)
            kappap = 1.0d0 + kappa
            kappam = 1.0d0 - kappa
            pre1 = sqrt(1.0d0 - tau2) / (tau * rho)
            term1 = -kappam * (2.0d0 * kappap + rhoi) * expi
            term2 = kappap * (2.0d0 * kappam + rhok) * expk
            intS = pre1 * (term1 + term2)
            pre2 = dmpik * (1.0d0 + tau) * pre1
            term1 = -kappam * expi
            term2 = (kappam + rhok) * expk
            intaAb = pre2 * (term1 + term2)
            pre2 = -dmpik * (1.0d0 - tau) * pre1
            term1 = -kappap * expk
            term2 = (kappap + rhoi) * expi
            intbBa = pre2 * (term1 + term2)
            intbAb = (1.0d0 - (1.0d0 + rhok) * expk2) / r
            intaBa = (1.0d0 - (1.0d0 + rhoi) * expi2) / r
            term1 = -kappam * kappam * (2.0d0 + kappa + rhoi) * expi2
            term2 = -kappap * kappap * (2.0d0 - kappa + rhok) * expk2
            intJ = (1.0d0 + (term1 + term2) / 4.0d0) / r
            intK = SSSS(dmpi,dmpk,0.0d0,r,.true.)
         end if
      else if (xreptyp .eq. 'UNITED ') then
         dmpik = sqrt(dmpi * dmpk)
         rho = dmpik * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         expik = exp(-rho)
         expik2 = expik * expik
         intS = (1.0d0 + rho + 1.0d0 / 3.0d0 * rho2) * expik
         intaAb = dmpik * (1.0d0 + rho) * expik
         intbBa = intaAb
         intbAb = (1.0d0 - (1.0d0 + rho) * expik2) / r
         intaBa = intbAb
         intJ = (1.0d0 - (1.0d0 + 11.0d0/8.0d0 * rho
     &          + 3.0d0/4.0d0 * rho2 + 1.0d0/6.0d0 * rho3) * expik2) / r
         intK = SSSS(dmpik,dmpik,0.0d0,r,.true.)
      end if
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
      a = 3.5
      b = 3.0
      z1 = -0.5
      z2 = 1.5
      print*, ""
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
