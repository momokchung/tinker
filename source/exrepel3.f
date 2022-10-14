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
                  e = hartree * (termS0 + termS1 + termS2) /
     &                          (1.0d0 - intS2)
                  print*, e
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
                  print*, xreptyp
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
            intJ = (1.0d0 - (1.0d0 + 11.0d0/8.0d0 * rho + 
     &           3.0d0/4.0d0 * rho2 + 1.0d0/6.0d0 * rho3) * expi2) / r
            call computeK (dmpi,dmpi,r2,xi,yi,zi,xk,yk,zk,intK)
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
            call computeK (dmpi,dmpk,r2,xi,yi,zi,xk,yk,zk,intK)
         end if
      else if (xreptyp .eq. 'UNITED') then
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
         intJ = (1.0d0 - (1.0d0 + 11.0d0/8.0d0 * rho + 
     &           3.0d0/4.0d0 * rho2 + 1.0d0/6.0d0 * rho3) * expik2) / r
         call computeK (dmpik,dmpik,r2,xi,yi,zi,xk,yk,zk,intK)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine computeK  --  two electron exchange integral  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "computeK" computes the two electron exchange integral using
c     STO-3G type basis
c
c
      subroutine computeK (dmpi,dmpj,r2,xi,yi,zi,xj,yj,zj,intK)
      use math
      implicit none
      integer i,im,j
      integer index
      real*8 dmpi,dmpj
      real*8 dmpi2,dmpj2
      real*8 r,r2
      real*8 xi,yi,zi
      real*8 xj,yj,zj
      real*8 xp,yp,zp
      real*8 xq,yq,zq
      real*8 xpq,ypq,zpq
      real*8 rpq2
      real*8 term1,term2
      real*8 term3
      real*8 ai,aj
      real*8 aiaj,akal
      real*8 aiajakal
      real*8 aij,akl
      real*8 ci,cj
      real*8 pre, preEx
      real*8 intK, boysF
      real*8 exptermi(3)
      real*8 exptermj(3)
      real*8 coeff(3)
      real*8 cArray(9)
      real*8 aijArray(9)
      real*8 aiajArray(9)
      real*8 xpArray(9)
      real*8 ypArray(9)
      real*8 zpArray(9)
c
c
c     compute exponent scaling
c
      dmpi2 = dmpi * dmpi
      dmpj2 = dmpj * dmpj
c
c     set STO-3G exponent and scaling (Levine 15.4: Basis Functions)
c
      coeff(1) = 0.444615d0
      coeff(2) = 0.535336d0
      coeff(3) = 0.154340d0
      exptermi(1) = 0.109814d0 * dmpi2
      exptermi(2) = 0.40575d0 * dmpi2
      exptermi(3) = 2.22746d0 * dmpi2
      exptermj(1) = 0.109814d0 * dmpj2
      exptermj(2) = 0.40575d0 * dmpj2
      exptermj(3) = 2.22746d0 * dmpj2
c
c     intermediate pairwise calculations
c
      do i = 1, 3
         im = i - 1
         do j = 1, 3
            index = 3 * im + j
            ai = exptermi(i)
            aj = exptermj(j)
            aiaj = ai + aj
            ci = coeff(i) * ai**(0.75d0)
            cj = coeff(j) * aj**(0.75d0)
            xp = (xi * ai + xj * aj) / aiaj
            yp = (yi * ai + yj * aj) / aiaj
            zp = (zi * ai + zj * aj) / aiaj
            cArray(index) = ci * cj
            aijArray(index) = ai * aj
            aiajArray(index) = aiaj
            xpArray(index) = xp
            ypArray(index) = yp
            zpArray(index) = zp
         end do
      end do
c
c     compute exchange integral
c
      intK = 0.0d0
      do i = 1, 9
         do j = i, 9
            pre = cArray(i) * cArray(j)
            aij = aijArray(i)
            akl = aijArray(j)
            aiaj = aiajArray(i)
            akal = aiajArray(j)
            aiajakal = aiaj + akal
            xp = xpArray(i)
            yp = ypArray(i)
            zp = zpArray(i)
            xq = xpArray(j)
            yq = ypArray(j)
            zq = zpArray(j)
            xpq = xp - xq
            ypq = yp - yq
            zpq = zp - zq
            rpq2 = xpq**2 + ypq**2 + zpq**2
            term1 = 1.0d0 / (aiaj * akal * sqrt(aiajakal))
            term2 = exp(-(aij / aiaj + akl / akal) * r2)
            term3 = boysF(aiaj * akal / aiajakal * rpq2)
            preEx = pre * term1 * term2 * term3
            if (i /= j) then
               preEx = 2.0d0 * preEx
            end if
            intK = intK + preEx
         end do
      end do
      intK = intK * (2.0d0 / pi)**3.0d0 * 2.0d0 * pi**2.5d0
      return
      end
c
c
c     #################################################
c     ##                                             ##
c     ##  function boysF  --  compute boys function  ##
c     ##                                             ##
c     #################################################
c
c
      function boysF (t)
      use math
      use xrepel
      implicit none
      integer j
      real*8 t,t2,d,x
      real*8 boysF
c
c
c     compute boys integral using Chebyshev expansion
c
      d = 0.0005d0
      x = t/d
      if (t >= 0.0d0 .and. t < 34.0d0) then
         j = int(x) * 4 + 1
         boysF = boysCoeff(j) + x * (boysCoeff(j+1) + 
     &         x * (boysCoeff(j+2) + x * boysCoeff(j+3)))
      else
         t2 = t * 2.0d0
         boysF = sqrt(pi / 2.0d0 / t2)
      end if
      return
      end
