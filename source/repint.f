c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     #################################################
c     ##                                             ##
c     ##  function boysF  --  compute boys function  ##
c     ##                                             ##
c     #################################################
c
c
c     "boysF" computes the boys integral using Chebyshev expansion
c
c     literature reference:
c
c     L. L. Shipman and R. E. Christoffersen, "High Speed Evaluation
c     of F0(x), Computer Physics Communications, 2, 201-206 (1971)
c
c     P. M. W. Gill, B. G. Johnson, and J. A. Pople, "Two-Electron
c     Repulsion Integrals Over Gaussian s Functions", International
c     Journal of Quantum Chemistry, 40, 745-752 (1991)
c
c     S. Obara and A. Saika, "Efficient Recursive Computation of
c     Molecular Integrals Over Cartesian Gaussian Functions",
c     The Journal of Chemical Physics, 84, 3963-3974 (1986)
c
c
      function boysF (m, t)
      use math
      use xrepel
      implicit none
      integer m,i,j,mj
      real*8 t,t2,d,x
      real*8 expt
      real*8 boysF
c
c
c     compute boys integral using Chebyshev expansion
c
      d = 0.0005d0
      x = t/d
      if (t >= 0.0d0 .and. t < 34.0d0) then
         j = int(x) * 4 + 1
         mj = m * ncoeff + j
         boysF = boysCoeff(mj) + x * (boysCoeff(mj+1)
     &         + x * (boysCoeff(mj+2) + x * boysCoeff(mj+3)))
      else
         t2 = t * 2.0d0
         boysF = sqrt(pi / 2.0d0 / t2)
         expt = exp(-t)
         do i = 0, m-1
            boysF = (boysF * (2.0d0 * i + 1.0d0) - expt) / t2
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine contract2  --  general 2-electron contraction  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "contract2" performs contraction over 2-electron integrals
c
c
      subroutine contract2 (ms, m, n, ssss, alphai, alphaj,
     &                                    alphak, alphal, result)
      use math
      use xrepel
      implicit none
      integer m,ms
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ni2,nj2,nk2,nl2
      real*8 ai,aj,ak,al
      real*8 ci,cj,ck,cl
      real*8 corri,corrj,corrk,corrl,corr
      real*8 ssss(m,nsto,nsto,nsto,nsto)
      real*8 result(m)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      ni2 = ni / 2.0d0
      nj2 = nj / 2.0d0
      nk2 = nk / 2.0d0
      nl2 = nl / 2.0d0
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      do i = 1, nsto
         ai = stoexp(nip,i) * alphai2
         ci = stocoeff(nip,i)
         corri = ci * (2.0d0 * ai / pi)**(0.75d0) * (4.0d0 * ai)**ni2
         do j = 1, nsto
            aj = stoexp(njp,j) * alphaj2
            cj = stocoeff(njp,j)
            corrj = cj * (2.0d0 * aj / pi)**(0.75d0) * (4.0d0 * aj)**nj2
            do k = 1, nsto
               ak = stoexp(nkp,k) * alphak2
               ck = stocoeff(nkp,k)
               corrk = ck * (2.0d0 * ak / pi)**(0.75d0)
     &                                               * (4.0d0 * ak)**nk2
               do l = 1, nsto
                  al = stoexp(nlp,l) * alphal2
                  cl = stocoeff(nlp,l)
                  corrl = cl * (2.0d0 * al / pi)**(0.75d0)
     &                                               * (4.0d0 * al)**nl2
                  corr = corri*corrj*corrk*corrl
                  do ii = ms, m
                     result(ii) = result(ii) + corr*ssss(ii,i,j,k,l)
                  end do
               end do
            end do
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine integralSSSS  --  primitive SSSS integral  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "integralSSSS" computes 2-electron integral over primitive basis
c
c
      subroutine integralSSSS (ms, m, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      use math
      use xrepel
      implicit none
      integer m,ms
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zij,zkl
      real*8 rij2,rkl2
      real*8 zp,zq,zpq,rpq2
      real*8 pre, exp_term
      real*8 fii, boysF
      real*8 ssss(m,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      zij = zi - zj
      zkl = zk - zl
      rij2 = zij**2
      rkl2 = zkl**2
      do i = 1, nsto
         ai = stoexp(nip,i) * alphai2
         do j = 1, nsto
            aj = stoexp(njp,j) * alphaj2
            aiaj = ai + aj
            zp = (ai * zi + aj * zj) / aiaj
            do k = 1, nsto
               ak = stoexp(nkp,k) * alphak2
               do l = 1, nsto
                  al = stoexp(nlp,l) * alphal2
                  akal = ak + al
                  aiajakal = aiaj + akal
                  zq = (ak * zk + al * zl) / akal
                  zpq = zp - zq
                  rpq2 = zpq**2
                  pre = 2.0d0*pi**2.5d0 / (aiaj * akal * sqrt(aiajakal))
                  exp_term = exp(-(ai*aj/aiaj*rij2 + ak*al/akal*rkl2))
                  pre = pre * exp_term
                  do ii = ms, m
                     fii = boysF(ii-1, aiaj * akal / aiajakal * rpq2)
                     ssss(ii,i,j,k,l) = pre * fii
                  end do
               end do
            end do
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine VRR1z  --  first vertical recursion relation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "VRR1z" is the first vertical recursion relation from OS recursion
c     in the "z" direction
c
c
      subroutine VRR1z (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                         zi, zj, zk, zl, eri1n, eri2n, eri1, eri2)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = (zp - zz)
     &       * eri1(ii,i,j,k,l) + (zw - zp) * eri1(ii+1,i,j,k,l)
                     end do
                  end do
               end do
            end do
         end do
      else if ((vn .eq. 3) .or. (vn .eq. 4)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = (zq - zz)
     &       * eri1(ii,i,j,k,l) + (zw - zq) * eri1(ii+1,i,j,k,l)
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR2Az  --  "2Az" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR2Az" is the second "A" vertical recursion relation from OS
c     recursion in the "z" direction
c
c
      subroutine VRR2Az (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &            zi, zj, zk, zl, eri1n, eri2n, eri3n, eri1, eri2, eri3)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n,eri3n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
      real*8 eri3(eri3n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri3(ii,i,j,k,l) = (zp - zz) * eri2(ii,i,j,k,l)
     &                 + (zw - zp) * eri2(ii+1,i,j,k,l)
     &                 + 1.0d0 / (2.0d0 * aiajakal) * eri1(ii+1,i,j,k,l)
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR2Ax  --  "2Ax" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR2Ax" is the second "A" vertical recursion relation from OS
c     recursion in the "x" direction
c
c
      subroutine VRR2Ax (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                         zi, zj, zk, zl, eri1n, eri2n, eri1, eri2)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = 1.0d0 / (2.0d0 * aiajakal)
     &                                   * eri1(ii+1,i,j,k,l)
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR2Bz  --  "2Bz" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR2Bz" is the second "B" vertical recursion relation from OS
c     recursion in the "z" direction
c
c
      subroutine VRR2Bz (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &            zi, zj, zk, zl, eri1n, eri2n, eri3n, eri1, eri2, eri3)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n,eri3n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
      real*8 eri3(eri3n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri3(ii,i,j,k,l) = (zp - zz) * eri2(ii,i,j,k,l)
     &                 + (zw - zp) * eri2(ii+1,i,j,k,l)
     &                 + 1.0d0 / (2.0d0 * aiaj) * (eri1(ii,i,j,k,l)
     &                 - akal / aiajakal * eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      else if ((vn .eq. 3) .or. (vn .eq. 4)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri3(ii,i,j,k,l) = (zq - zz) * eri2(ii,i,j,k,l)
     &                 + (zw - zq) * eri2(ii+1,i,j,k,l)
     &                 + 1.0d0 / (2.0d0 * akal) * (eri1(ii,i,j,k,l)
     &                 - aiaj / aiajakal * eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR2Bx  --  "2Bx" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR2Bx" is the second "B" vertical recursion relation from OS
c     recursion in the "x" direction
c
c
      subroutine VRR2Bx (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                         zi, zj, zk, zl, eri1n, eri2n, eri1, eri2)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = 1.0d0 / (2.0d0 * aiaj)
     &                           * (eri1(ii,i,j,k,l)
     &                           - akal / aiajakal * eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      else if ((vn .eq. 3) .or. (vn .eq. 4)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = 1.0d0 / (2.0d0 * akal)
     &                           * (eri1(ii,i,j,k,l)
     &                           - aiaj / aiajakal * eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR3Az  --  "3Az" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR3Az" is the third "A" vertical recursion relation from OS
c     recursion in the "z" direction
c
c
      subroutine VRR3Az (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                   zi, zj, zk, zl, eri1n, eri2n, eri3n, eri4n,
     &                                           eri1, eri2, eri3, eri4)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n,eri3n,eri4n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
      real*8 eri3(eri3n,nsto,nsto,nsto,nsto)
      real*8 eri4(eri4n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri4(ii,i,j,k,l) = (zp - zz) * eri3(ii,i,j,k,l)
     &                 + (zw - zp) * eri3(ii+1,i,j,k,l)
     &                 + 1.0d0 / (2.0d0 * aiajakal)
     &                 * (eri2(ii+1,i,j,k,l) + eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      else if ((vn .eq. 3) .or. (vn .eq. 4)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri4(ii,i,j,k,l) = (zq - zz) * eri3(i,i,j,k,l)
     &                 + (zw - zq) * eri3(ii+1,i,j,k,l)
     &                 + 1.0d0 / (2.0d0 * aiajakal)
     &                 * (eri2(ii+1,i,j,k,l) + eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR3Bz  --  "3Bz" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR3Bz" is the third "B" vertical recursion relation from OS
c     recursion in the "z" direction
c
c
      subroutine VRR3Bz (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                   zi, zj, zk, zl, eri1n, eri2n, eri1, eri2)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = (zp - zz) * eri1(ii,i,j,k,l)
     &                                  + (zw - zp) * eri1(ii+1,i,j,k,l)
                     end do
                  end do
               end do
            end do
         end do
      else if ((vn .eq. 3) .or. (vn .eq. 4)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = (zq - zz) * eri1(ii,i,j,k,l)
     &                                  + (zw - zq) * eri1(ii+1,i,j,k,l)
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR3Bx  --  "3Bx" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR3Bx" is the third "B" vertical recursion relation from OS
c     recursion in the "x" direction
c
c
      subroutine VRR3Bx (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                   zi, zj, zk, zl, eri1n, eri2n, eri1, eri2)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      do i = 1, nsto
         ai = stoexp(nip,i) * alphai2
         do j = 1, nsto
            aj = stoexp(njp,j) * alphaj2
            aiaj = ai + aj
            zp = (ai * zi + aj * zj) / aiaj
            do k = 1, nsto
               ak = stoexp(nkp,k) * alphak2
               do l = 1, nsto
                  al = stoexp(nlp,l) * alphal2
                  akal = ak + al
                  aiajakal = aiaj + akal
                  zq = (ak * zk + al * zl) / akal
                  zw = (aiaj * zp + akal * zq) / aiajakal
                  do ii = ms, m
                     eri2(ii,i,j,k,l) = 1.0d0 / (2.0d0 * aiajakal)
     &                               * eri1(ii+1,i,j,k,l)
                  end do
               end do
            end do
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR4Az  --  "4Az" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR4Az" is the fourth "A" vertical recursion relation from OS
c     recursion in the "z" direction
c
c
      subroutine VRR4Az (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                zi, zj, zk, zl, eri1n, eri2n, eri3n, eri4n, eri5n,
     &                                     eri1, eri2, eri3, eri4, eri5)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n,eri3n,eri4n,eri5n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
      real*8 eri3(eri3n,nsto,nsto,nsto,nsto)
      real*8 eri4(eri4n,nsto,nsto,nsto,nsto)
      real*8 eri5(eri5n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri5(ii,i,j,k,l) = (zp - zi) * eri4(ii,i,j,k,l)
     &                  + (zw - zp) * eri4(ii+1,i,j,k,l)
     &                  + 1.0d0 / (2.0d0 * aiaj) * (eri3(ii,i,j,k,l)
     &                  - akal / aiajakal * eri3(ii+1,i,j,k,l))
     &                  + 1.0d0 / (2.0d0 * aiajakal)
     &                  * (eri2(ii+1,i,j,k,l) + eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR4Ax  --  "4Ax" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR4Ax" is the fourth "A" vertical recursion relation from OS
c     recursion in the "x" direction
c
c
      subroutine VRR4Ax (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                       zi, zj, zk, zl, eri1n, eri2n, eri3n, eri4n,
     &                                           eri1, eri2, eri3, eri4)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n,eri3n,eri4n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
      real*8 eri3(eri3n,nsto,nsto,nsto,nsto)
      real*8 eri4(eri4n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri4(ii,i,j,k,l) = 1.0d0 / (2.0d0 * aiaj) 
     &                  * (eri3(ii,i,j,k,l) - akal / aiajakal
     &                  * eri3(ii+1,i,j,k,l))
     &                  + 1.0d0 / (2.0d0 * aiajakal)
     &                  * (eri2(ii+1,i,j,k,l) + eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR4Bx  --  "4Bx" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR4Bx" is the fourth "B" vertical recursion relation from OS
c     recursion in the "x" direction
c
c
      subroutine VRR4Bx (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, eri1n, eri2n, eri1, eri2)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = 1.0d0 / (2.0d0 * aiaj)
     &                  * (eri1(ii,i,j,k,l)
     &                  - akal / aiajakal * eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      else if ((vn .eq. 3) .or. (vn .eq. 4)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = 1.0d0 / (2.0d0 * akal)
     &                  * (eri1(ii,i,j,k,l)
     &                  - aiaj / aiajakal * eri1(ii+1,i,j,k,l))
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine VRR4Cx  --  "4Cx" vertical recursion relation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "VRR4Cx" is the fourth "C" vertical recursion relation from OS
c     recursion in the "x" direction
c
c
      subroutine VRR4Cx (ms, m, n, vn, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, eri1n, eri2n, eri1, eri2)
      use math
      use xrepel
      implicit none
      integer m,ms,vn
      integer eri1n,eri2n
      integer ni,nj,nk,nl
      integer nip,njp,nkp,nlp
      integer i,j,k,l,ii
      integer n(*)
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 alphai2,alphaj2,alphak2,alphal2
      real*8 ai,aj,ak,al
      real*8 aiaj,akal,aiajakal
      real*8 zz,zp,zq,zw
      real*8 eri1(eri1n,nsto,nsto,nsto,nsto)
      real*8 eri2(eri2n,nsto,nsto,nsto,nsto)
c
c
      alphai2 = alphai * alphai
      alphaj2 = alphaj * alphaj
      alphak2 = alphak * alphak
      alphal2 = alphal * alphal
      ni = n(1)
      nj = n(2)
      nk = n(3)
      nl = n(4)
      if (vn .eq. 1) then
         zz = zi
      else if (vn .eq. 2) then
         zz = zj
      else if (vn .eq. 3) then
         zz = zk
      else if (vn .eq. 4) then
         zz = zl
      end if
      nip = ni + 1
      njp = nj + 1
      nkp = nk + 1
      nlp = nl + 1
      if ((vn .eq. 1) .or. (vn .eq. 2)) then
         do i = 1, nsto
            ai = stoexp(nip,i) * alphai2
            do j = 1, nsto
               aj = stoexp(njp,j) * alphaj2
               aiaj = ai + aj
               zp = (ai * zi + aj * zj) / aiaj
               do k = 1, nsto
                  ak = stoexp(nkp,k) * alphak2
                  do l = 1, nsto
                     al = stoexp(nlp,l) * alphal2
                     akal = ak + al
                     aiajakal = aiaj + akal
                     zq = (ak * zk + al * zl) / akal
                     zw = (aiaj * zp + akal * zq) / aiajakal
                     do ii = ms, m
                        eri2(ii,i,j,k,l) = 1.0d0 / (2.0d0 * aiajakal)
     &                  * eri1(ii+1,i,j,k,l)
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
      end
c
c
c     #################################################
c     ##                                             ##
c     ##  function gSSSS  --  general SSSS integral  ##
c     ##                                             ##
c     #################################################
c
c
c     "gSSSS" computes 2-electron SSSS integral
c
c
      function gSSSS (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSSSS
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 0
      n(3) = 0
      n(4) = 0
      call integralSSSS(1, m, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call contract2(1,m,n,ssss,alphai,alphaj,alphak,alphal,result)
      gSSSS = result(1)
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  function gSSSPz  --  general SSSPz integral  ##
c     ##                                               ##
c     ###################################################
c
c
c     "gSSSPz" computes 2-electron SSSPz integral
c
c
      function gSSSPz (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSSSPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(2,nsto,nsto,nsto,nsto)
      real*8 sssp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 0
      n(3) = 0
      n(4) = 1
      call integralSSSS(1, m+1, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssss, sssp)
      call contract2(1,m,n,sssp,alphai,alphaj,alphak,alphal,result)
      gSSSPz = result(1)
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  function gSSPzS  --  general SSPzS integral  ##
c     ##                                               ##
c     ###################################################
c
c
c     "gSSPzS" computes 2-electron SSPzS integral
c
c
      function gSSPzS (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSSPzS
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(2,nsto,nsto,nsto,nsto)
      real*8 ssps(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 0
      n(3) = 1
      n(4) = 0
      call integralSSSS(1, m+1, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssss, ssps)
      call contract2(1,m,n,ssps,alphai,alphaj,alphak,alphal,result)
      gSSPzS = result(1)
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  function gSPzSS  --  general SPzSS integral  ##
c     ##                                               ##
c     ###################################################
c
c
c     "gSPzSS" computes 2-electron SPzSS integral
c
c
      function gSPzSS (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSPzSS
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(2,nsto,nsto,nsto,nsto)
      real*8 spss(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 1
      n(3) = 0
      n(4) = 0
      call integralSSSS(1, m+1, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssss, spss)
      call contract2(1,m,n,spss,alphai,alphaj,alphak,alphal,result)
      gSPzSS = result(1)
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  function gSSPzPz  --  general SPzPz integral  ##
c     ##                                                ##
c     ####################################################
c
c
c     "gSSPzPz" computes 2-electron SSPzPz integral
c
c
      function gSSPzPz (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSSPzPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 ssps(2,nsto,nsto,nsto,nsto)
      real*8 sspp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 0
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+1, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, ssps)
      call VRR2Bz(1, m, n, 4, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 3, 2, 1, ssss, ssps, sspp)
      call contract2(1,m,n,sspp,alphai,alphaj,alphak,alphal,result)
      gSSPzPz = result(1)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function gSSPxPx  --  general SSPxPx integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "gSSPxPx" computes 2-electron SSPxPx integral
c
c
      function gSSPxPx (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSSPxPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(2,nsto,nsto,nsto,nsto)
      real*8 sspp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 0
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+1, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Bx(1, m, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssss, sspp)
      call contract2(1,m,n,sspp,alphai,alphaj,alphak,alphal,result)
      gSSPxPx = result(1)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function gSPzSPz  --  general SPzSPz integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "gSPzSPz" computes 2-electron SPzSPz integral
c
c
      function gSPzSPz (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSPzSPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 sssp(2,nsto,nsto,nsto,nsto)
      real*8 spsp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 1
      n(3) = 0
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, sssp)
      call VRR2Az(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 3, 2, 1, ssss, sssp, spsp)
      call contract2(1,m,n,spsp,alphai,alphaj,alphak,alphal,result)
      gSPzSPz = result(1)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function gSPxSPx  --  general SPxSPx integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "gSPxSPx" computes 2-electron SPxSPx integral
c
c
      function gSPxSPx (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSPxSPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(2,nsto,nsto,nsto,nsto)
      real*8 spsp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 1
      n(3) = 0
      n(4) = 1
      call integralSSSS(2, m+1, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Ax(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssss, spsp)
      call contract2(1,m,n,spsp,alphai,alphaj,alphak,alphal,result)
      gSPxSPx = result(1)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function gSPzPzS  --  general SPzPzS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "gSPzPzS" computes 2-electron SPzPzS integral
c
c
      function gSPzPzS (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSPzPzS
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 ssps(2,nsto,nsto,nsto,nsto)
      real*8 spps(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 1
      n(3) = 1
      n(4) = 0
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+1, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, ssps)
      call VRR2Az(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 3, 2, 1, ssss, ssps, spps)
      call contract2(1,m,n,spps,alphai,alphaj,alphak,alphal,result)
      gSPzPzS = result(1)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function gSPxPxS  --  general SPxPxS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "gSPxPxS" computes 2-electron SPxPxS integral
c
c
      function gSPxPxS (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSPxPxS
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(2,nsto,nsto,nsto,nsto)
      real*8 spps(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 1
      n(3) = 1
      n(4) = 0
      call integralSSSS(2, m+1, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Ax(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssss, spps)
      call contract2(1,m,n,spps,alphai,alphaj,alphak,alphal,result)
      gSPxPxS = result(1)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function gPzSPzS  --  general PzSPzS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "gPzSPzS" computes 2-electron PzSPzS integral
c
c
      function gPzSPzS (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gPzSPzS
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 ssps(2,nsto,nsto,nsto,nsto)
      real*8 psps(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 0
      n(3) = 1
      n(4) = 0
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+1, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, ssps)
      call VRR2Az(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 3, 2, 1, ssss, ssps, psps)
      call contract2(1,m,n,psps,alphai,alphaj,alphak,alphal,result)
      gPzSPzS = result(1)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function gPxSPxS  --  general PxSPxS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "gPxSPxS" computes 2-electron PxSPxS integral
c
c
      function gPxSPxS (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gPxSPxS
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(2,nsto,nsto,nsto,nsto)
      real*8 psps(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 0
      n(3) = 1
      n(4) = 0
      call integralSSSS(2, m+1, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Ax(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssss, psps)
      call contract2(1,m,n,psps,alphai,alphaj,alphak,alphal,result)
      gPxSPxS = result(1)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function gPzPzSS  --  general PzPzSS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "gPzPzSS" computes 2-electron PzPzSS integral
c
c
      function gPzPzSS (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gPzPzSS
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 spss(2,nsto,nsto,nsto,nsto)
      real*8 ppss(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 0
      n(4) = 0
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, spss)
      call VRR2Bz(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 3, 2, 1, ssss, spss, ppss)
      call contract2(1,m,n,ppss,alphai,alphaj,alphak,alphal,result)
      gPzPzSS = result(1)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function gPxPxSS  --  general PxPxSS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "gPxPxSS" computes 2-electron PxPxSS integral
c
c
      function gPxPxSS (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gPxPxSS
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(2,nsto,nsto,nsto,nsto)
      real*8 ppss(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 0
      n(4) = 0
      call integralSSSS(1, m+1, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Bx(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssss, ppss)
      call contract2(1,m,n,ppss,alphai,alphaj,alphak,alphal,result)
      gPxPxSS = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gSPzPzPz  --  general SPzPzPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gSPzPzPz" computes 2-electron SPzPzPz integral
c
c
      function gSPzPzPz (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSPzPzPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(4,nsto,nsto,nsto,nsto)
      real*8 ssps(3,nsto,nsto,nsto,nsto)
      real*8 sssp(2,nsto,nsto,nsto,nsto)
      real*8 sspp(2,nsto,nsto,nsto,nsto)
      real*8 sppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+3, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+2, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 3, ssss, ssps)
      call VRR1z(2, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 2, ssss, sssp)
      call VRR2Bz(1, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 4, 3, 2, ssss, ssps, sspp)
      call VRR3Az(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &               zi, zj, zk, zl, 3, 2, 2, 1, ssps, sssp, sspp, sppp)
      call contract2(1,m,n,sppp,alphai,alphaj,alphak,alphal,result)
      gSPzPzPz = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gSPzPxPx  --  general SPzPxPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gSPzPxPx" computes 2-electron SPzPxPx integral
c
c
      function gSPzPxPx (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSPzPxPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 sspp(2,nsto,nsto,nsto,nsto)
      real*8 sppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Bx(1, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, sspp)
      call VRR3Bz(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, sspp, sppp)
      call contract2(1,m,n,sppp,alphai,alphaj,alphak,alphal,result)
      gSPzPxPx = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gSPxPzPx  --  general SPxPzPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gSPxPzPx" computes 2-electron SPxPzPx integral
c
c
      function gSPxPzPx (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSPxPzPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 ssps(2,nsto,nsto,nsto,nsto)
      real*8 sppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(2, m+1, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, ssps)
      call VRR3Bx(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssps, sppp)
      call contract2(1,m,n,sppp,alphai,alphaj,alphak,alphal,result)
      gSPxPzPx = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gSPxPxPz  --  general SPxPxPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gSPxPxPz" computes 2-electron SPxPxPz integral
c
c
      function gSPxPxPz (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      real*8 gSPxPxPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 sssp(2,nsto,nsto,nsto,nsto)
      real*8 sppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 0
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(2, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, sssp)
      call VRR3Bx(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, sssp, sppp)
      call contract2(1,m,n,sppp,alphai,alphaj,alphak,alphal,result)
      gSPxPxPz = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gPzSPzPz  --  general PzSPzPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gPzSPzPz" computes 2-electron PzSPzPz integral
c
c
      function gPzSPzPz (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPzSPzPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(4,nsto,nsto,nsto,nsto)
      real*8 ssps(3,nsto,nsto,nsto,nsto)
      real*8 sssp(2,nsto,nsto,nsto,nsto)
      real*8 sspp(2,nsto,nsto,nsto,nsto)
      real*8 pspp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 0
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+3, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+2, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 3, ssss, ssps)
      call VRR1z(2, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 2, ssss, sssp)
      call VRR2Bz(1, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 4, 3, 2, ssss, ssps, sspp)
      call VRR3Az(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &               zi, zj, zk, zl, 3, 2, 2, 1, ssps, sssp, sspp, pspp)
      call contract2(1,m,n,pspp,alphai,alphaj,alphak,alphal,result)
      gPzSPzPz = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gPzSPxPx  --  general PzSPxPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gPzSPxPx" computes 2-electron PzSPxPx integral
c
c
      function gPzSPxPx (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPzSPxPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 sspp(2,nsto,nsto,nsto,nsto)
      real*8 pspp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 0
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Bx(1, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, sspp)
      call VRR3Bz(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, sspp, pspp)
      call contract2(1,m,n,pspp,alphai,alphaj,alphak,alphal,result)
      gPzSPxPx = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gPxSPzPx  --  general PxSPzPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gPxSPzPx" computes 2-electron PzSPzPx integral
c
c
      function gPxSPzPx (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxSPzPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 ssps(2,nsto,nsto,nsto,nsto)
      real*8 pspp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 0
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(2, m+1, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, ssps)
      call VRR3Bx(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ssps, pspp)
      call contract2(1,m,n,pspp,alphai,alphaj,alphak,alphal,result)
      gPxSPzPx = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gPxSPxPz  --  general PxSPxPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gPxSPxPz" computes 2-electron PzSPzPz integral
c
c
      function gPxSPxPz (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxSPxPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 sssp(2,nsto,nsto,nsto,nsto)
      real*8 pspp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 0
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(2, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, sssp)
      call VRR3Bx(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, sssp, pspp)
      call contract2(1,m,n,pspp,alphai,alphaj,alphak,alphal,result)
      gPxSPxPz = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gPzPzSPz  --  general PzPzSPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gPzPzSPz" computes 2-electron PzPzSPz integral
c
c
      function gPzPzSPz (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPzPzSPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(4,nsto,nsto,nsto,nsto)
      real*8 psss(3,nsto,nsto,nsto,nsto)
      real*8 spss(2,nsto,nsto,nsto,nsto)
      real*8 ppss(2,nsto,nsto,nsto,nsto)
      real*8 ppsp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 0
      n(4) = 1
      call integralSSSS(1, m+3, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+2, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 3, ssss, psss)
      call VRR1z(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 2, ssss, spss)
      call VRR2Bz(1, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 4, 3, 2, ssss, psss, ppss)
      call VRR3Az(1, m, n, 4, alphai, alphaj, alphak, alphal,
     &               zi, zj, zk, zl, 3, 2, 2, 1, psss, spss, ppss, ppsp)
      call contract2(1,m,n,ppsp,alphai,alphaj,alphak,alphal,result)
      gPzPzSPz = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gPxPzSPx  --  general PxPzSPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gPxPzSPx" computes 2-electron PxPzSPx integral
c
c
      function gPxPzSPx (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxPzSPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 spss(2,nsto,nsto,nsto,nsto)
      real*8 ppsp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 0
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, spss)
      call VRR3Bx(1, m, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, spss, ppsp)
      call contract2(1,m,n,ppsp,alphai,alphaj,alphak,alphal,result)
      gPxPzSPx = result(1)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function gPxPxSPz  --  general PxPxSPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gPxPxSPz" computes 2-electron PxPxSPz integral
c
c
      function gPxPxSPz (alphai, alphaj, alphak, alphal, zi, zj, zk, zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxPxSPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 ppss(2,nsto,nsto,nsto,nsto)
      real*8 ppsp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 0
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Bx(1, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, ppss)
      call VRR3Bz(1, m, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ppss, ppsp)
      call contract2(1,m,n,ppsp,alphai,alphaj,alphak,alphal,result)
      gPxPxSPz = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPzPzPzPz  --  general PzPzPzPz integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPzPzPzPz" computes 2-electron PzPzPzPz integral
c
c
      function gPzPzPzPz (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPzPzPzPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(5,nsto,nsto,nsto,nsto)
      real*8 ssps(4,nsto,nsto,nsto,nsto)
      real*8 sssp(3,nsto,nsto,nsto,nsto)
      real*8 sspp(3,nsto,nsto,nsto,nsto)
      real*8 spps(2,nsto,nsto,nsto,nsto)
      real*8 spsp(2,nsto,nsto,nsto,nsto)
      real*8 sppp(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+4, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+3, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 5, 4, ssss, ssps)
      call VRR1z(2, m+2, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 5, 3, ssss, sssp)
      call VRR2Bz(1, m+2, n, 4, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 5, 4, 3, ssss, ssps, sspp)
      call VRR2Az(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 5, 4, 2, ssss, ssps, spps)
      call VRR2Az(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 5, 3, 2, ssss, sssp, spsp)
      call VRR3Az(1, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &               zi, zj, zk, zl, 4, 3, 3, 2, ssps, sssp, sspp, sppp)
      call VRR4Az(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &      zi, zj, zk, zl, 2, 2, 3, 2, 1, spps, spsp, sspp, sppp, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPzPzPzPz = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPzPzPxPx  --  general PzPzPxPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPzPzPxPx" computes 2-electron PxPxPzPz integral
c
c
      function gPzPzPxPx (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPzPzPxPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(4,nsto,nsto,nsto,nsto)
      real*8 psss(3,nsto,nsto,nsto,nsto)
      real*8 ppss(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+3, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+2, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 3, ssss, psss)
      call VRR2Bz(1, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 4, 3, 2, ssss, psss, ppss)
      call VRR4Bx(1, m, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, ppss, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPzPzPxPx = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPzPxPzPx  --  general PzPxPzPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPzPxPzPx" computes 2-electron PxPxPzPz integral
c
c
      function gPzPxPzPx (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPzPxPzPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(4,nsto,nsto,nsto,nsto)
      real*8 ssps(3,nsto,nsto,nsto,nsto)
      real*8 psps(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(2, m+3, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(2, m+2, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 3, ssss, ssps)
      call VRR2Az(2, m+1, n, 1, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 4, 3, 2, ssss, ssps, psps)
      call VRR4Cx(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, psps, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPzPxPzPx = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPxPzPzPx  --  general PxPzPzPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPxPzPzPx" computes 2-electron PxPxPzPz integral
c
c
      function gPxPzPzPx (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxPzPzPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(4,nsto,nsto,nsto,nsto)
      real*8 ssps(3,nsto,nsto,nsto,nsto)
      real*8 spps(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(2, m+3, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(2, m+2, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 3, ssss, ssps)
      call VRR2Az(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 4, 3, 2, ssss, ssps, spps)
      call VRR4Cx(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, spps, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPxPzPzPx = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPxPzPxPz  --  general PxPzPxPz integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPxPzPxPz" computes 2-electron PxPxPxPz integral
c
c
      function gPxPzPxPz (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxPzPxPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(4,nsto,nsto,nsto,nsto)
      real*8 sssp(3,nsto,nsto,nsto,nsto)
      real*8 spsp(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(2, m+3, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(2, m+2, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 3, ssss, sssp)
      call VRR2Az(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 4, 3, 2, ssss, sssp, spsp)
      call VRR4Cx(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, spsp, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPxPzPxPz = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPxPxPzPz  --  general PxPxPzPz integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPxPxPzPz" computes 2-electron PxPxPzPz integral
c
c
      function gPxPxPzPz (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxPxPzPz
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(4,nsto,nsto,nsto,nsto)
      real*8 ssps(3,nsto,nsto,nsto,nsto)
      real*8 sspp(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+3, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR1z(1, m+2, n, 3, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 4, 3, ssss, ssps)
      call VRR2Bz(1, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                        zi, zj, zk, zl, 4, 3, 2, ssss, ssps, sspp)
      call VRR4Bx(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, sspp, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPxPxPzPz = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPxPxPxPx  --  general PxPxPxPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPxPxPxPx" computes 2-electron PxPxPxPx integral
c
c
      function gPxPxPxPx (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxPxPxPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 sspp(2,nsto,nsto,nsto,nsto)
      real*8 spps(2,nsto,nsto,nsto,nsto)
      real*8 spsp(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Bx(1, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, sspp)
      call VRR2Ax(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, spps)
      call VRR2Ax(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, spsp)
      call VRR4Ax(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &      zi, zj, zk, zl, 2, 2, 2, 1, spps, spsp, sspp, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPxPxPxPx = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPxPxPyPy  --  general PxPxPyPy integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPxPxPyPy" computes 2-electron PxPxPyPy integral
c
c
      function gPxPxPyPy (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxPxPyPy
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 sspp(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(1, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Bx(1, m+1, n, 4, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, sspp)
      call VRR4Bx(1, m, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, sspp, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPxPxPyPy = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPxPyPxPy  --  general PxPyPxPy integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPxPyPxPy" computes 2-electron PxPyPxPy integral
c
c
      function gPxPyPxPy (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxPyPxPy
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 spsp(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(2, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Ax(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, spsp)
      call VRR4Cx(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, spsp, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPxPyPxPy = result(1)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function gPxPyPyPx  --  general PxPyPyPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "gPxPyPyPx" computes 2-electron PxPyPyPx integral
c
c
      function gPxPyPyPx (alphai, alphaj, alphak, alphal, zi,zj,zk,zl)
      use xrepel
      implicit none
      integer m
      integer n(4)
      integer ii,i,j,k,l
      real*8 gPxPyPyPx
      real*8 alphai,alphaj,alphak,alphal
      real*8 zi,zj,zk,zl
      real*8 result(1)
      real*8 ssss(3,nsto,nsto,nsto,nsto)
      real*8 spps(2,nsto,nsto,nsto,nsto)
      real*8 pppp(1,nsto,nsto,nsto,nsto)
c
c
      result(1) = 0.0d0
      m = 1
      n(1) = 1
      n(2) = 1
      n(3) = 1
      n(4) = 1
      call integralSSSS(2, m+2, n, alphai, alphaj, alphak, alphal,
     &                                             zi, zj, zk, zl, ssss)
      call VRR2Ax(2, m+1, n, 2, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 3, 2, ssss, spps)
      call VRR4Cx(1, m, n, 1, alphai, alphaj, alphak, alphal,
     &                                 zi, zj, zk, zl, 2, 1, spps, pppp)
      call contract2(1,m,n,pppp,alphai,alphaj,alphak,alphal,result)
      gPxPyPyPx = result(1)
      return
      end
c
c
c     ##################################################
c     ##                                              ##
c     ##  function SSSS  --  computes SSSS  integral  ##
c     ##                                              ##
c     ##################################################
c
c
c     "SSSS" computes 2-electron integral using STO basis
c
c
      function SSSS (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SSSS,gSSSS
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SSSS = gSSSS (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SSSS = gSSSS (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  function SSSPz  --  computes SSSPz integral  ##
c     ##                                               ##
c     ###################################################
c
c
c     "SSSPz" computes 2-electron integral using STO basis
c
c
      function SSSPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SSSPz,gSSSPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SSSPz = gSSSPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SSSPz = gSSSPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  function SSPzS  --  computes SSPzS integral  ##
c     ##                                               ##
c     ###################################################
c
c
c     "SSPzS" computes 2-electron integral using STO basis
c
c
      function SSPzS (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SSPzS,gSSPzS
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SSPzS = gSSPzS (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SSPzS = gSSPzS (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  function SPzSS  --  computes SPzSS integral  ##
c     ##                                               ##
c     ###################################################
c
c
c     "SPzSS" computes 2-electron integral using STO basis
c
c
      function SPzSS (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SPzSS,gSPzSS
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SPzSS = gSPzSS (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SPzSS = gSPzSS (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function SSPzPz  --  computes SSPzPz integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "SSPzPz" computes 2-electron integral using STO basis
c
c
      function SSPzPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SSPzPz,gSSPzPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SSPzPz = gSSPzPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SSPzPz = gSSPzPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function SSPxPx  --  computes SSPxPx integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "SSPxPx" computes 2-electron integral using STO basis
c
c
      function SSPxPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SSPxPx,gSSPxPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SSPxPx = gSSPxPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SSPxPx = gSSPxPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function SPzSPz  --  computes SPzSPz integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "SPzSPz" computes 2-electron integral using STO basis
c
c
      function SPzSPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SPzSPz,gSPzSPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SPzSPz = gSPzSPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SPzSPz = gSPzSPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function SPxSPx  --  computes SPxSPx integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "SPxSPx" computes 2-electron integral using STO basis
c
c
      function SPxSPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SPxSPx,gSPxSPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SPxSPx = gSPxSPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SPxSPx = gSPxSPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function SPzPzS  --  computes SPzPzS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "SPzPzS" computes 2-electron integral using STO basis
c
c
      function SPzPzS (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SPzPzS,gSPzPzS
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SPzPzS = gSPzPzS (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SPzPzS = gSPzPzS (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function SPxPxS  --  computes SPxPxS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "SPxPxS" computes 2-electron integral using STO basis
c
c
      function SPxPxS (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SPxPxS,gSPxPxS
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SPxPxS = gSPxPxS (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SPxPxS = gSPxPxS (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function PzSPzS  --  computes PzSPzS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "PzSPzS" computes 2-electron integral using STO basis
c
c
      function PzSPzS (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PzSPzS,gPzSPzS
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PzSPzS = gPzSPzS (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PzSPzS = gPzSPzS (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function PxSPxS  --  computes PxSPxS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "PxSPxS" computes 2-electron integral using STO basis
c
c
      function PxSPxS (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxSPxS,gPxSPxS
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxSPxS = gPxSPxS (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxSPxS = gPxSPxS (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function PzPzSS  --  computes PzPzSS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "PzPzSS" computes 2-electron integral using STO basis
c
c
      function PzPzSS (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PzPzSS,gPzPzSS
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PzPzSS = gPzPzSS (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PzPzSS = gPzPzSS (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  function PxPxSS  --  computes PxPxSS integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "PxPxSS" computes 2-electron integral using STO basis
c
c
      function PxPxSS (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPxSS,gPxPxSS
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPxSS = gPxPxSS (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPxSS = gPxPxSS (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function SPzPzPz  --  computes SPzPzPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "SPzPzPz" computes 2-electron integral using STO basis
c
c
      function SPzPzPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SPzPzPz,gSPzPzPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SPzPzPz = gSPzPzPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SPzPzPz = gSPzPzPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function SPzPxPx  --  computes SPzPxPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "SPzPxPx" computes 2-electron integral using STO basis
c
c
      function SPzPxPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SPzPxPx,gSPzPxPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SPzPxPx = gSPzPxPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SPzPxPx = gSPzPxPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function SPxPzPx  --  computes SPxPzPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "SPxPzPx" computes 2-electron integral using STO basis
c
c
      function SPxPzPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SPxPzPx,gSPxPzPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SPxPzPx = gSPxPzPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SPxPzPx = gSPxPzPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function SPxPxPz  --  computes SPxPxPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "SPxPxPz" computes 2-electron integral using STO basis
c
c
      function SPxPxPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 SPxPxPz,gSPxPxPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         SPxPxPz = gSPxPxPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         SPxPxPz = gSPxPxPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function PzSPzPz  --  computes PzSPzPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "PzSPzPz" computes 2-electron integral using STO basis
c
c
      function PzSPzPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PzSPzPz,gPzSPzPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PzSPzPz = gPzSPzPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PzSPzPz = gPzSPzPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function PzSPxPx  --  computes PzSPxPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "PzSPxPx" computes 2-electron integral using STO basis
c
c
      function PzSPxPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PzSPxPx,gPzSPxPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PzSPxPx = gPzSPxPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PzSPxPx = gPzSPxPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function PxSPzPx  --  computes PxSPzPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "PxSPzPx" computes 2-electron integral using STO basis
c
c
      function PxSPzPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxSPzPx,gPxSPzPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxSPzPx = gPxSPzPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxSPzPx = gPxSPzPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function PxSPxPz  --  computes PxSPxPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "PxSPxPz" computes 2-electron integral using STO basis
c
c
      function PxSPxPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxSPxPz,gPxSPxPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxSPxPz = gPxSPxPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxSPxPz = gPxSPxPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function PzPzSPz  --  computes PzPzSPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "PzPzSPz" computes 2-electron integral using STO basis
c
c
      function PzPzSPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PzPzSPz,gPzPzSPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PzPzSPz = gPzPzSPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PzPzSPz = gPzPzSPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function PxPzSPx  --  computes PxPzSPx integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "PxPzSPx" computes 2-electron integral using STO basis
c
c
      function PxPzSPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPzSPx,gPxPzSPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPzSPx = gPxPzSPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPzSPx = gPxPzSPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function PxPxSPz  --  computes PxPxSPz integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "PxPxSPz" computes 2-electron integral using STO basis
c
c
      function PxPxSPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPxSPz,gPxPxSPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPxSPz = gPxPxSPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPxSPz = gPxPxSPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PzPzPzPz  --  computes PzPzPzPz integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PzPzPzPz" computes 2-electron integral using STO basis
c
c
      function PzPzPzPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PzPzPzPz,gPzPzPzPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PzPzPzPz = gPzPzPzPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PzPzPzPz = gPzPzPzPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PzPzPxPx  --  computes PzPzPxPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PzPzPxPx" computes 2-electron integral using STO basis
c
c
      function PzPzPxPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PzPzPxPx,gPzPzPxPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PzPzPxPx = gPzPzPxPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PzPzPxPx = gPzPzPxPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PzPxPzPx  --  computes PzPxPzPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PzPxPzPx" computes 2-electron integral using STO basis
c
c
      function PzPxPzPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PzPxPzPx,gPzPxPzPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PzPxPzPx = gPzPxPzPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PzPxPzPx = gPzPxPzPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PxPzPzPx  --  computes PxPzPzPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PxPzPzPx" computes 2-electron integral using STO basis
c
c
      function PxPzPzPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPzPzPx,gPxPzPzPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPzPzPx = gPxPzPzPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPzPzPx = gPxPzPzPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PxPzPxPz  --  computes PxPzPxPz integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PxPzPxPz" computes 2-electron integral using STO basis
c
c
      function PxPzPxPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPzPxPz,gPxPzPxPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPzPxPz = gPxPzPxPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPzPxPz = gPxPzPxPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PxPxPzPz  --  computes PxPxPzPz integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PxPxPzPz" computes 2-electron integral using STO basis
c
c
      function PxPxPzPz (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPxPzPz,gPxPxPzPz
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPxPzPz = gPxPxPzPz (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPxPzPz = gPxPxPzPz (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PxPxPxPx  --  computes PxPxPxPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PxPxPxPx" computes 2-electron integral using STO basis
c
c
      function PxPxPxPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPxPxPx,gPxPxPxPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPxPxPx = gPxPxPxPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPxPxPx = gPxPxPxPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PxPxPyPy  --  computes PxPxPyPy integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PxPxPyPy" computes 2-electron integral using STO basis
c
c
      function PxPxPyPy (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPxPyPy,gPxPxPyPy
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPxPyPy = gPxPxPyPy (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPxPyPy = gPxPxPyPy (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PxPyPxPy  --  computes PxPyPxPy integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PxPyPxPy" computes 2-electron integral using STO basis
c
c
      function PxPyPxPy (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPyPxPy,gPxPyPxPy
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPyPxPy = gPxPyPxPy (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPyPxPy = gPxPyPxPy (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function PxPyPyPx  --  computes PxPyPyPx integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "PxPyPyPx" computes 2-electron integral using STO basis
c
c
      function PxPyPyPx (alphai, alphaj, zi, zj, exch)
      use xrepel
      implicit none
      real*8 PxPyPyPx,gPxPyPyPx
      real*8 alphai,alphaj,zi,zj
      logical exch
c
c
      if (exch) then
         PxPyPyPx = gPxPyPyPx (alphai,alphaj,alphai,alphaj,zi,zj,zi,zj)
      else
         PxPyPyPx = gPxPyPyPx (alphai,alphai,alphaj,alphaj,zi,zi,zj,zj)
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function coulombSSSS  --  SSSS coulomb integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "coulombSSSS" computes 2-electron coulomb integral
c
c
      function coulombSSSS (a, b, z1, z2)
      implicit none
      real*8 coulombSSSS
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 term1,term2
c
c
      diff = abs(a - b)
      eps = 0.005d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         coulombSSSS = 1.0d0 / r * (1.0d0 - (1.0d0 + 11.0d0/8.0d0*rho
     &             + 3.0d0/4.0d0*rho2 + rho3/6.0d0) * exp(-2.0d0 * rho))
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         term1 = -(1.0d0 - kappa)**2 / 4.0d0 * (2.0d0 + kappa + rhoA)
     &                                                 * exp(-2. * rhoA)
         term2 = -(1.0d0 + kappa)**2 / 4.0d0 * (2.0d0 - kappa + rhoB)
     &                                                 * exp(-2. * rhoB)
         coulombSSSS = 1.0d0 / r * (1.0d0 + term1 + term2)
      end if
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function coulombSSPzS  --  SSPzS coulomb integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "coulombSSPzS" computes 2-electron coulomb integral
c
c
      function coulombSSPzS (a, b, z1, z2)
      implicit none
      real*8 coulombSSPzS,coul
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoB2,rhoB3
      real*8 kappa2,kappam,kappap
      real*8 kappam3,kappap2
c
c
      diff = abs(a - b)
      eps = 0.005d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         coul = a / rho2 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &                   + 59.0d0/48.0d0 * rho3 + 11.0d0/24.0d0 * rho4
     &                   + rho5 / 12.0d0) * exp(-2.0d0 * rho))
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho2 = rho**2
         rhoA2 = rhoA**2
         rhoB2 = rhoB**2
         rhoB3 = rhoB2 * rhoB
         kappa2 = kappa**2
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappam3 = kappam**3
         kappap2 = kappap**2
         pre = alpha / ((1.0d0 - tau) * rho2)
         term1 = -kappam3 * (1.0d0/16.0d0 * (5.0d0 + 3.0d0 * kappa)
     &                    * (1.0d0 + 2.0d0 * rhoA) + rhoA2 / 4.0d0)
     &                    * exp(-2.0d0 * rhoA)
         term2 = -kappap2 * (1.0d0/16.0d0 * (11.0d0 - 10.0d0 * kappa
     &                     + 3.0d0 * kappa2) * (1.0d0 + 2.0d0 * rhoB)
     &                     + 0.5d0 * (2.0d0 - kappa) * rhoB2 + rhoB3
     &                     / 4.0d0) * exp(-2.0d0 * rhoB)
         coul = pre * (1.0d0 + term1 + term2)
         end if
      coulombSSPzS = -coul
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  function coulombSSPzPz  --  SSPzPz coulomb integral  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "coulombSSPzPz" computes 2-electron coulomb integral
c
c
      function coulombSSPzPz (a, b, z1, z2)
      implicit none
      real*8 coulombSSPzPz,coul1,coul2
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3
      real*8 rhoB2,rhoB3,rhoB4,rhoB5
      real*8 kappa2,kappa3
      real*8 kappap,kappap2
      real*8 kappam,kappam3,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.005d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         exp2 = exp(-2.0d0 * rho)
         coul1 = 1.0d0 / r * (1.0d0 - (1.0d0 + 25.0d0/16.0d0 * rho
     &        + 9.0d0/8.0d0 * rho2 + 23.0d0/48.0d0 * rho3 + rho4 / 8.0d0
     &        + rho5 / 60.0d0) * exp2)
         coul2 = a / rho3 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &        + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4 + 31.0d0/120.0d0
     &        * rho5 + 13.0d0/180.0d0 * rho6 + rho7 / 90.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap2 = kappap**2
         kappam3 = kappam**3
         kappam4 = kappam3 * kappam
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 1.0d0 / r
         term1 = -kappam3 * (1.0d0/16.0d0 * (1.0d0 - 5.0d0 * kappa
     &           - 4.0d0 * kappa2) - 1.0d0/8.0d0 * kappa * rhoA) * exp2A
         term2 = -kappap2 * (1.0d0/16.0d0 * (15.0d0 - 22.0d0 * kappa
     &           + 15.0d0 * kappa2 - 4.0d0 * kappa3) + 3.0d0/8.0d0
     &           * (3.0d0 - 3.0d0 * kappa + kappa2) * rhoB + 1.0d0/4.0d0
     &           * (2.0d0 - kappa) * rhoB2 + rhoB3/12.0d0) * exp2B
         coul1 = pre * (1.0d0 + term1 + term2)
         pre = alpha / ((1.0d0 - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/16.0d0 * (3.0d0 + 2.0d0 * kappa)
     &         * (1.0d0 + 2.0d0 * rhoA) + 1.0d0/24.0d0 * (7.0d0 + 4.0d0
     &         * kappa) * rhoA2 + rhoA3 / 12.0d0) * exp2A
         term2 = -kappap2 * (1.0d0/16.0d0 * (13.0d0 - 16.0d0 * kappa
     &       + 9.0d0 * kappa2 - 2.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &       + 1.0d0/24.0d0 * (37.0d0 - 42.0d0 * kappa + 21.0d0 * kappa2
     &       - 4.0d0 * kappa3) * rhoB2 + 1.0d0/12.0d0 * (11.0d0 - 10.0d0
     &       * kappa + 3.0d0 * kappa2) * rhoB3 + 1.0d0/6.0d0 * (2.0d0
     &       - kappa) * rhoB4 + rhoB5 / 18.0d0) * exp2B
         coul2 = pre * (1.0d0 + term1 + term2)
         end if
      coulombSSPzPz = coul1 + 3.0d0 * coul2
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  function coulombSSPxPx  --  SSPxPx coulomb integral  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "coulombSSPxPx" computes 2-electron coulomb integral
c
c
      function coulombSSPxPx (a, b, z1, z2)
      implicit none
      real*8 coulombSSPxPx,coul1,coul2
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3
      real*8 rhoB2,rhoB3,rhoB4,rhoB5
      real*8 kappa2,kappa3
      real*8 kappap,kappap2
      real*8 kappam,kappam3,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         exp2 = exp(-2.0d0 * rho)
         coul1 = 1.0d0 / r * (1.0d0 - (1.0d0 + 25.0d0/16.0d0 * rho
     &      + 9.0d0/8.0d0 * rho2 + 23.0d0/48.0d0 * rho3 + rho4 / 8.0d0
     &      + rho5 / 60.0d0) * exp2)
         coul2 = a / rho3 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &      + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4 + 31.0d0/120.0d0
     &      * rho5 + 13.0d0/180.0d0 * rho6 + rho7 / 90.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap2 = kappap**2
         kappam3 = kappam**3
         kappam4 = kappam3 * kappam
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 1.0d0 / r
         term1 = -kappam3 * (1.0d0/16.0d0 * (1.0d0 - 5.0d0 * kappa
     &      - 4.0d0 * kappa2) - 1.0d0/8.0d0 * kappa * rhoA) * exp2A
         term2 = -kappap2 * (1.0d0/16.0d0 * (15.0d0 - 22.0d0 * kappa
     &      + 15.0d0 * kappa2 - 4.0d0 * kappa3) + 3.0d0/8.0d0
     &      * (3.0d0- 3.0d0 * kappa + kappa2) * rhoB + 1.0d0/4.0d0
     &      * (2.0d0 - kappa) * rhoB2 + rhoB3 / 12.0d0) * exp2B
         coul1 = pre * (1.0d0 + term1 + term2)
         pre = alpha / ((1.0d0 - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/16.0d0 * (3.0d0 + 2.0d0 * kappa)
     &      * (1.0d0 + 2.0d0 * rhoA) + 1.0d0/24.0d0 * (7.0d0 + 4.0d0
     &      * kappa) * rhoA2 + rhoA3 / 12.0d0) * exp2A
         term2 = -kappap2 * (1.0d0/16.0d0 * (13.0d0 - 16.0d0 * kappa
     &      + 9.0d0 * kappa2 - 2.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/24.0d0 * (37.0d0 - 42.0d0 * kappa + 21.0d0 * kappa2
     &      - 4.0d0 * kappa3) * rhoB2 + 1.0d0/12.0d0 * (11.0d0 - 10.0d0
     &      * kappa + 3.0d0 * kappa2) * rhoB3 + 1.0d0/6.0d0
     &      * (2.0d0 - kappa) * rhoB4 + rhoB5 / 18.0d0) * exp2B
         coul2 = pre * (1.0d0 + term1 + term2)
         end if
      coulombSSPxPx = coul1 - 1.5d0 * coul2
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  function coulombSPzSPz  --  SPzSPz coulomb integral  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "coulombSPzSPz" computes 2-electron coulomb integral
c
c
      function coulombSPzSPz (a, b, z1, z2)
      implicit none
      real*8 coulombSPzSPz,coul
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4
      real*8 rhoB2,rhoB3,rhoB4
      real*8 kappa2,kappa3
      real*8 kappap,kappap3
      real*8 kappam,kappam3
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.005d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         exp2 = exp(-2.0d0 * rho)
         coul = 2.0d0 * a / rho3 * (1.0d0 - (1.0d0 + 2.0d0 * rho
     &        + 2.0d0 * rho2 + 263.0d0/192.0d0 * rho3 + 71.0d0/96.0d0
     &        * rho4 + 77.0d0/240.0d0 * rho5 + rho6 / 10.0d0
     &        + rho7 / 60.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)

         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappam3 = kappam**3
         kappap3 = kappap**3
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 2.0d0 * alpha / ((1.0d0 + tau) * (1. 0d0- tau) * rho3)
         term1 = -kappam3 * (1.0d0/16.0d0 * (8.0d0 + 9.0d0 * kappa
     &        + 3.0d0 * kappa2) * (1.0d0 + 2.0d0 * rhoA + 2.0d0 * rhoA2)
     &        + 3.0d0/16.0d0 * (3.0d0 + 2.0d0 * kappa) * rhoA3
     &        + rhoA4 / 8.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/16.0d0 * (8.0d0 - 9.0d0 * kappa
     &        + 3.0d0 * kappa2) * (1.0d0 + 2.0d0 * rhoB + 2.0d0 * rhoB2)
     &        + 3.0d0/16.0d0 * (3.0d0 - 2.0d0 * kappa) * rhoB3
     &        + rhoB4 / 8.0d0) * exp2B
         coul = pre * (1.0d0 + term1 + term2)
         end if
      coulombSPzSPz = -coul
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  function coulombSPxSPx  --  SPxSPx coulomb integral  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "coulombSPxSPx" computes 2-electron coulomb integral
c
c
      function coulombSPxSPx (a, b, z1, z2)
      implicit none
      real*8 coulombSPxSPx,coul
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3
      real*8 rhoB2,rhoB3
      real*8 kappa2
      real*8 kappap,kappap3
      real*8 kappam,kappam3
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         exp2 = exp(-2.0d0 * rho)
         coul = a / rho3 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &      + 121.0d0/96.0d0 * rho3 + 25.0d0/48.0d0 * rho4
     &      + 2.0d0/15.0d0 * rho5 + rho6 / 60.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         kappa2 = kappa * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam3 = kappam**3
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = alpha / ((1.0d0 + tau) * (1. - tau) * rho3)
         term1 = -kappam3 * (1.0d0/16.0d0 * (8.0d0 + 9.0d0 * kappa
     &      + 3.0d0 * kappa2) * (1.0d0 + 2.0d0 * rhoA) + 1.0d0/8.0d0
     &      * (5.0d0 + 3.0d0 * kappa) * rhoA2 + rhoA3 / 8.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/16.0d0 * (8.0d0 - 9.0d0 * kappa
     &      + 3.0d0 * kappa2) * (1.0d0 + 2.0d0 * rhoB) + 1.0d0/8.0d0
     &      * (5.0d0 - 3.0d0 * kappa) * rhoB2 + rhoB3 / 8.0d0) * exp2B
         coul = pre * (1.0d0 + term1 + term2)
         end if
      coulombSPxSPx = coul
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function coulombSPzPzPz  --  SPzPzPz coulomb integral  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "coulombSPzPzPz" computes 2-electron coulomb integral
c
c
      function coulombSPzPzPz (a, b, z1, z2)
      implicit none
      real*8 coulombSPzPzPz,coul1,coul2
      real*8 a,b,z1,z2
      real*8 c
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5
      real*8 rhoB2,rhoB3,rhoB4,rhoB5,rhoB6
      real*8 kappa2,kappa3
      real*8 kappap,kappap3
      real*8 kappam,kappam3,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         rho8 = rho7 * rho
         rho9 = rho8 * rho
         exp2 = exp(-2.0d0 * rho)
         coul1 = a / rho2 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &    + 83.0d0/64.0d0 * rho3 + 19.0d0/32.0d0 * rho4 + 47.0d0/240.0d0
     &    * rho5 + 2.0d0/45.0d0 * rho6 + rho7 / 180.0d0) * exp2)
         coul2 = 3.0d0 * a / rho4 * (1.0d0 - (1.0d0 + 2.0d0 * rho
     &    + 2.0d0 * rho2 + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4
     &    + 583.0d0/2160.0d0 * rho5 + 103.0d0/1080.0d0 * rho6
     &    + 11.0d0/360.0d0 * rho7 + 13.0d0/1620.0d0 * rho8
     &    + rho9 / 810.0d0) * exp2)
      else
         c = b
         b = a
         a = c
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho2 = rho**2
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam3 = kappam**3
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = alpha / ((1.0d0 - tau) * rho2)
         term1 = -kappam3 * (1.0d0/16.0d0 * (13.0d0 + 24.0d0 * kappa
     &      + 18.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/8.0d0 * (11.0d0 + 15.0d0 * kappa + 6.0d0 * kappa2)
     &      * rhoA2 + 1.0d0/24.0d0 * (13.0d0 + 9.0d0 * kappa) * rhoA3
     &      + rhoA4 / 12.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/16.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/8.0d0 * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2)
     &      * rhoB2 + 1.0d0/8.0d0 * kappa * rhoB3) * exp2B
         coul1 = pre * (1.0d0 + term1 + term2)
         c = b
         b = a
         a = c
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho4 = rho**4
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         rhoB6 = rhoB5 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam4 = kappam**4
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 3.0d0 * alpha / ((1.0d0 + tau) * (1. - tau)**2 * rho4)
         term1 = -kappam4 * (1.0d0/32.0d0 * (11.0d0 + 14.0d0 * kappa
     &      + 5.0d0 * kappa2) * (1.0d0 + 2.0d0 * rhoA) + 1.0d0/72.0d0
     &      * (47.0d0 + 58.0d0 * kappa + 20.0d0 * kappa2) * rhoA2
     &      + 1.0d0/36.0d0 * (14.0d0 + 16.0d0 * kappa + 5.0d0 * kappa2)
     &      * rhoA3 + 1.0d0/72.0d0 * (11.0d0 + 8.0d0 * kappa) * rhoA4
     &      + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (21.0d0 - 33.0d0 * kappa
     &      + 21.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/72.0d0 * (92.0d0 - 141.0d0 * kappa + 87.0d0 * kappa2
     &      - 20.0d0 * kappa3) * rhoB2 + 1.0d0/36.0d0 * (29.0d0 - 42.0d0
     &      * kappa + 24.0d0 * kappa2 - 5.0d0 * kappa3) * rhoB3
     &      + 1.0d0/24.0d0 * (9.0d0 - 11.0d0 * kappa + 4.0d0 * kappa2)
     &      * rhoB4 + 1.0d0/108.0d0 * (13.0d0 - 9.0d0 * kappa) * rhoB5
     &      + rhoB6 / 54.0d0) * exp2B
         coul2 = pre * (1.0d0 + term1 + term2)
         end if
      coulombSPzPzPz = coul1 + 3.0d0 * coul2
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function coulombSPzPxPx  --  SPzPxPx coulomb integral  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "coulombSPzPxPx" computes 2-electron coulomb integral
c
c
      function coulombSPzPxPx (a, b, z1, z2)
      implicit none
      real*8 coulombSPzPxPx,coul1,coul2
      real*8 a,b,z1,z2
      real*8 c
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5
      real*8 rhoB2,rhoB3,rhoB4,rhoB5,rhoB6
      real*8 kappa2,kappa3
      real*8 kappap,kappap3
      real*8 kappam,kappam3,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         rho8 = rho7 * rho
         rho9 = rho8 * rho
         exp2 = exp(-2.0d0 * rho)
         coul1 = a / rho2 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &     + 83.0d0/64.0d0 * rho3 + 19.0d0/32.0d0 * rho4
     &     + 47.0d0/240.0d0 * rho5 + 2.0d0/45.0d0 * rho6
     &     + rho7 / 180.0d0) * exp2)
         coul2 = 3.0d0 * a / rho4 * (1.0d0 - (1.0d0 + 2.0d0 * rho
     &     + 2.0d0 * rho2 + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4
     &     + 583.0d0/2160.0d0 * rho5 + 103.0d0/1080.0d0 * rho6
     &     + 11.0d0/360.0d0 * rho7 + 13.0d0/1620.0d0 * rho8
     &     + rho9 / 810.0d0) * exp2)
      else
         c = b
         b = a
         a = c
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho2 = rho**2
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam3 = kappam**3
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = alpha / ((1.0d0 - tau) * rho2)
         term1 = -kappam3 * (1.0d0/16.0d0 * (13.0d0 + 24.0d0 * kappa
     &      + 18.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/8.0d0 * (11.0d0 + 15.0d0 * kappa + 6.0d0 * kappa2)
     &      * rhoA2 + 1.0d0/24.0d0 * (13.0d0 + 9.0d0 * kappa) * rhoA3
     &      + rhoA4 / 12.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/16.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/8.0d0 * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2)
     &      * rhoB2 + 1.0d0/8.0d0 * kappa * rhoB3) * exp2B
         coul1 = pre * (1.0d0 + term1 + term2)
         c = b
         b = a
         a = c
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho4 = rho**4
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         rhoB6 = rhoB5 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam4 = kappam**4
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 3.0d0 * alpha / ((1.0d0 + tau) * (1.0d0 - tau)**2 * rho4)
         term1 = -kappam4 * (1.0d0/32.0d0 * (11.0d0 + 14.0d0 * kappa
     &      + 5.0d0 * kappa2) * (1.0d0 + 2.0d0 * rhoA) + 1.0d0/72.0d0
     &      * (47.0d0 + 58.0d0 * kappa + 20.0d0 * kappa2) * rhoA2
     &      + 1.0d0/36.0d0 * (14.0d0 + 16.0d0 * kappa + 5.0d0 * kappa2)
     &      * rhoA3 + 1.0d0/72.0d0 * (11.0d0 + 8.0d0 * kappa) * rhoA4
     &      + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (21.0d0 - 33.0d0 * kappa
     &      + 21.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/72.0d0 * (92.0d0 - 141.0d0 * kappa + 87.0d0 * kappa2
     &      - 20.0d0 * kappa3) * rhoB2 + 1.0d0/36.0d0 * (29.0d0 - 42.0d0
     &      * kappa + 24.0d0 * kappa2 - 5.0d0 * kappa3) * rhoB3
     &      + 1.0d0/24.0d0 * (9.0d0 - 11.0d0 * kappa + 4.0d0 * kappa2)
     &      * rhoB4 + 1.0d0/108.0d0 * (13.0d0 - 9.0d0 * kappa) * rhoB5
     &      + rhoB6 / 54.0d0) * exp2B
         coul2 = pre * (1.0d0 + term1 + term2)
      end if
      coulombSPzPxPx = coul1 - 1.50d0 * coul2
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function coulombSPxPxPz  --  SPxPxPz coulomb integral  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "coulombSPxPxPz" computes 2-electron coulomb integral
c
c
      function coulombSPxPxPz (a, b, z1, z2)
      implicit none
      real*8 coulombSPxPxPz,coul
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7,rho8
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4
      real*8 rhoB2,rhoB3,rhoB4,rhoB5
      real*8 kappa2,kappa3
      real*8 kappap,kappap3
      real*8 kappam,kappam3,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         rho8 = rho7 * rho
         exp2 = exp(-2.0d0 * rho)
         coul = sqrt(3.0d0) * a / rho4 * (1.0d0 - (1.0d0 + 2.0d0 * rho
     &      + 2.0d0 * rho2 + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4
     &      + 377.0d0/1440.0d0 * rho5 + 19.0d0/240.0d0 * rho6
     &      + rho7 / 60.0d0 + rho8 / 540.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho4 = rho**4
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam4 = kappam**4
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = sqrt(3.0d0) * alpha / ((1.0d0 + tau)
     &                                * (1.0d0 - tau)**2 * rho4)
         term1 = -kappam4 * (1.0d0/32.0d0 * (11.0d0 + 14.0d0 * kappa
     &      + 5.0d0 * kappa2) * (1.0d0 + 2.0d0 * rhoA) + 1.0d0/24.0d0
     &      * (14.0d0 + 16.0d0 * kappa + 5.0d0 * kappa2) * rhoA2
     &      + 1.0d0/12.0d0 * (3.0d0 + 2.0d0 * kappa) * rhoA3
     &      + rhoA4 / 24.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (21.0d0 - 33.0d0 * kappa
     &      + 21.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/24.0d0 * (29.0d0 - 42.0d0 * kappa + 24.0d0 * kappa2
     &      - 5.0d0 * kappa3) * rhoB2 + 1.0d0/12.0d0 * (8.0d0 - 9.0d0
     &      * kappa + 3.0d0 * kappa2) * rhoB3 + 1.0d0/24.0d0
     &      * (5.0d0 - 3.0d0 * kappa) * rhoB4 + rhoB5 / 36.0d0) * exp2B
         coul = pre * (1.0d0 + term1 + term2)
      end if
      coulombSPxPxPz = -1.5d0 * sqrt(3.0d0) * coul
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function coulombPzPzPzPz  --  PzPzPzPz coulomb integral  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "coulombPzPzPzPz" computes 2-electron coulomb integral
c
c
      function coulombPzPzPzPz (a, b, z1, z2)
      implicit none
      real*8 coulombPzPzPzPz,coul1,coul2,coul3,coul4
      real*8 a,b,z1,z2
      real*8 c
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5,rhoA6,rhoA7
      real*8 rhoB2,rhoB3,rhoB4,rhoB5,rhoB6,rhoB7
      real*8 kappa2,kappa3,kappa4
      real*8 kappap,kappap3,kappap4
      real*8 kappam,kappam3,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         rho8 = rho7 * rho
         rho9 = rho8 * rho
         rho10 = rho9 * rho
         rho11 = rho10 * rho
         exp2 = exp(-2.0d0 * rho)
         coul1 = 1.0d0 / r * (1.0d0 - (1.0d0 + 419.0d0/256.0d0 * rho
     &      + 163.0d0/128.0d0 * rho2 + 119.0d0/192.0d0 * rho3
     &      + 5.0d0/24.0d0 * rho4 + rho5 / 20.0d0 + rho6 / 120.0d0
     &      + rho7 / 1260.0d0) * exp2)
         coul2 = a / rho3 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &      + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4 + 191.0d0/720.0d0
     &      * rho5 + 31.0d0/360.0d0 * rho6 + 19.0d0/840.0d0 * rho7
     &      + 17.0d0/3780.0d0 * rho8 + rho9 / 1890.0d0) * exp2)
         coul3 = coul2
         coul4 = 6.0d0 * a / rho5 * (1.0d0 - (1.0d0 + 2.0d0 * rho
     &      + 2.0d0 * rho2 + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4
     &      + 511.0d0/1920.0d0 * rho5 + 253.0d0/2880.0d0 * rho6
     &      + 2237.0d0/90720.0d0 * rho7 + 71.0d0/11340.0d0 * rho8
     &      + rho9 / 630.0d0 + 13.0d0/34020.0d0 * rho10
     &      + rho11 / 17010.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rho5 = rho**5
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoA6 = rhoA5 * rhoA
         rhoA7 = rhoA6 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         rhoB6 = rhoB5 * rhoB
         rhoB7 = rhoB6 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappa4 = kappa3 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappap4 = kappap3 * kappap
         kappam3 = kappam**3
         kappam4 = kappam3 * kappam
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 1.0d0 / r
         term1 = -kappam3 * (1.0d0/16.0d0 * (8.0d0 - kappa - 27.0d0
     &      * kappa2 - 30.0d0 * kappa3 - 10.0d0 * kappa4) + 1.0d0/32.0d0
     &      * (11.0d0 - 19.0d0 * kappa - 44.0d0 * kappa2 - 20.0d0
     &      * kappa3) * rhoA + 1.0d0/16.0d0 * (1.0d0 - 5.0d0 * kappa
     &      - 4.0d0 * kappa2) * rhoA2 - 1.0d0/24.0d0 * kappa * rhoA3)
     &      * exp2A
         term2 = -kappap3 * (1.0d0/16.0d0 * (8.0d0 + kappa - 27.0d0
     &      * kappa2 + 30.0d0 * kappa3 - 10.0d0 * kappa4) + 1.0d0/32.0d0
     &      * (11.0d0 + 19.0d0 * kappa - 44.0d0 * kappa2 + 20.0d0
     &      * kappa3) * rhoB + 1.0d0/16.0d0 * (1.0d0 + 5.0d0 * kappa
     &      - 4.0d0 * kappa2) * rhoB2 + 1.0d0/24.0d0 * kappa * rhoB3)
     &      * exp2B
         coul1 = pre * (1.0d0 + term1 + term2)
         pre = alpha / ((1.0d0 - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/32.0d0 * (21.0d0 + 44.0d0 * kappa
     &      + 35.0d0 * kappa2 + 10.0d0 * kappa3) * (1.0d0 + 2.0d0
     &      * rhoA) + 1.0d0/24.0d0 * (29.0d0 + 56.0d0 * kappa + 40.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA2 + 1.0d0/12.0d0 * (8.0d0
     &      + 12.0d0 * kappa + 5.0d0 * kappa2) * rhoA3 + 1.0d0/24.0d0
     &      * (5.0d0 + 4.0d0 * kappa) * rhoA4 + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (11.0d0 + 7.0d0 * kappa
     &      - 39.0d0 * kappa2 + 35.0d0 * kappa3 - 10.0d0 * kappa4)
     &      * (1.0d0 + 2.0d0 * rhoB) + 1.0d0/24.0d0 * (14.0d0 + 13.0d0
     &      * kappa - 51.0d0 * kappa2 + 40.0d0 * kappa3 - 10.0d0 *
     &      kappa4) * rhoB2 + 1.0d0/12.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * rhoB3 + 1.0d0/24.0d0
     &      * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2) * rhoB4
     &      + 1.0d0/36.0d0 * kappa * rhoB5) * exp2B
         coul2 = pre * (1.0d0 + term1 + term2)
         pre = 6.0d0 * alpha / ((1.0d0 + tau)**2 * (1.0d0 - tau)**2
     &               * rho5)
         term1 = -kappam4 * (1.0d0/32.0d0 * (16.0d0 + 29.0d0 * kappa
     &      + 20.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/144.0d0 * (139.0d0 + 246.0d0 * kappa + 165.0d0
     &      * kappa2 + 40.0d0 * kappa3) * rhoA2 + 1.0d0/72.0d0 * (43.0d0
     &      + 72.0d0 * kappa + 45.0d0 * kappa2 + 10.0d0 * kappa3)
     &      * rhoA3 + 1.0d0/216.0d0 * (57.0d0 + 88.0d0 * kappa + 50.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA4 + 1.0d0/432.0d0
     &      * (39.0d0 + 52.0d0 * kappa + 20.0d0 * kappa2) * rhoA5
     &      + 1.0d0/216.0d0 * (5.0d0 + 4.0d0 * kappa) * rhoA6
     &      + rhoA7 / 324.0d0) * exp2A
         term2 = -kappap4 * (1.0d0/32.0d0 * (16.0d0 - 29.0d0 * kappa
     &      + 20.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/144.0d0 * (139.0d0 - 246.0d0 * kappa + 165.0d0
     &      * kappa2 - 40.0d0 * kappa3) * rhoB2 + 1.0d0/72.0d0 * (43.0d0
     &      - 72.0d0 * kappa + 45.0d0 * kappa2 - 10.0d0 * kappa3)
     &      * rhoB3 + 1.0d0/216.0d0 * (57.0d0 - 88.0d0 * kappa + 50.0d0
     &      * kappa2 - 10.0d0 * kappa3) * rhoB4 + 1.0d0/432.0d0
     &      * (39.0d0 - 52.0d0 * kappa + 20.0d0 * kappa2) * rhoB5
     &      + 1.0d0/216.0d0 * (5.0d0 - 4.0d0 * kappa) * rhoB6
     &      + rhoB7 / 324.0d0) * exp2B
         coul4 = pre * (1.0d0 + term1 + term2)
         c = b
         b = a
         a = c
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappa4 = kappa3 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam4 = kappam**4
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = alpha / ((1. - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/32.0d0 * (21.0d0 + 44.0d0 * kappa
     &      + 35.0d0 * kappa2 + 10.0d0 * kappa3) * (1.0d0 + 2.0d0
     &      * rhoA) + 1.0d0/24.0d0 * (29.0d0 + 56.0d0 * kappa + 40.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA2 + 1.0d0/12.0d0 * (8.0d0
     &      + 12.0d0 * kappa + 5.0d0 * kappa2) * rhoA3 + 1.0d0/24.0d0
     &      * (5.0d0 + 4.0d0 * kappa) * rhoA4 + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (11.0d0 + 7.0d0 * kappa
     &      - 39.0d0 * kappa2 + 35.0d0 * kappa3 - 10.0d0 * kappa4)
     &      * (1.0d0 + 2.0d0 * rhoB) + 1.0d0/24.0d0 * (14.0d0 + 13.0d0
     &      * kappa - 51.0d0 * kappa2 + 40.0d0 * kappa3 - 10.0d0 *
     &      kappa4) * rhoB2 + 1.0d0/12.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * rhoB3 + 1.0d0/24.0d0
     &      * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2) * rhoB4
     &      + 1.0d0/36.0d0 * kappa * rhoB5) * exp2B
         coul3 = pre * (1.0d0 + term1 + term2)
         c = b
         b = a
         a = c
         end if
      coulombPzPzPzPz = coul1 + 3.0d0 * coul2 + 3.0d0 * coul3
     &                        + 9.0d0 * coul4
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function coulombPxPzPxPz  --  PxPzPxPz coulomb integral  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "coulombPxPzPxPz" computes 2-electron coulomb integral
c
c
      function coulombPxPzPxPz (a, b, z1, z2)
      implicit none
      real*8 coulombPxPzPxPz,coul
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5,rhoA6
      real*8 rhoB2,rhoB3,rhoB4,rhoB5,rhoB6
      real*8 kappa2,kappa3,kappa4
      real*8 kappap,kappap4
      real*8 kappam,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         rho8 = rho7 * rho
         rho9 = rho8 * rho
         rho10 = rho9 * rho
         exp2 = exp(-2.0d0 * rho)
         coul = 4.0d0 * a / rho5 * (1.0d0 - (1.0d0 + 2.0d0 * rho
     &      + 2.0d0 * rho2 + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4
     &      + 1027.0d0/3840.0d0 * rho5 + 521.0d0/5760.0d0 * rho6
     &      + 67.0d0/2520.0d0 * rho7 + 67.0d0/10080.0d0 * rho8
     &      + 19.0d0/15120.0d0 * rho9 + rho10 / 7560.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho5 = rho**5
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoA6 = rhoA5 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         rhoB6 = rhoB5 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappa4 = kappa3 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap4 = kappap**4
         kappam4 = kappam**4
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 4.0d0*alpha/((1.0d0 + tau)**2 * (1.0d0 - tau)**2 * rho5)
         term1 = -kappam4 * (1.0d0/32.0d0 * (16.0d0 + 29.0d0 * kappa
     &      + 20.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/96.0d0 * (91.0d0 + 159.0d0 * kappa + 105.0d0
     &      * kappa2 + 25.0d0 * kappa3) * rhoA2 + 1.0d0/48.0d0
     &      * (27.0d0 + 43.0d0 * kappa + 25.0d0 * kappa2 + 5.0d0
     &      * kappa3) * rhoA3 + 1.0d0/48.0d0 * (11.0d0 + 14.0d0 * kappa
     &      + 5.0d0 * kappa2) * rhoA4 + 1.0d0/288.0d0 * (17.0d0 + 12.0d0
     &      * kappa) * rhoA5 + rhoA6 / 144.0d0) * exp2A
         term2 = -kappap4 * (1.0d0/32.0d0 * (16.0d0 - 29.0d0 * kappa
     &      + 20.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/96.0d0 * (91.0d0 - 159.0d0 * kappa + 105.0d0
     &      * kappa2 - 25.0d0 * kappa3) * rhoB2 + 1.0d0/48.0d0
     &      * (27.0d0 - 43.0d0 * kappa + 25.0d0 * kappa2 - 5.0d0
     &      * kappa3) * rhoB3 + 1.0d0/48.0d0 * (11.0d0 - 14.0d0 * kappa
     &      + 5.0d0 * kappa2) * rhoB4 + 1.0d0/288.0d0 * (17.0d0 - 12.0d0
     &      * kappa) * rhoB5 + rhoB6 / 144.0d0) * exp2B
         coul = pre * (1.0d0 + term1 + term2)
      end if
      coulombPxPzPxPz = -27.0d0/4.0d0 * coul
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function coulombPxPxPzPz  --  PxPxPzPz coulomb integral  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "coulombPxPxPzPz" computes 2-electron coulomb integral
c
c
      function coulombPxPxPzPz (a, b, z1, z2)
      implicit none
      real*8 coulombPxPxPzPz,coul1,coul2,coul3,coul4
      real*8 a,b,z1,z2
      real*8 c
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5,rhoA6,rhoA7
      real*8 rhoB2,rhoB3,rhoB4,rhoB5,rhoB6,rhoB7
      real*8 kappa2,kappa3,kappa4
      real*8 kappap,kappap3,kappap4
      real*8 kappam,kappam3,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         c = b
         b = a
         a = c
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         rho8 = rho7 * rho
         rho9 = rho8 * rho
         rho10 = rho9 * rho
         rho11 = rho10 * rho
         exp2 = exp(-2.0d0 * rho)
         coul1 = 1.0d0 / r * (1.0d0 - (1.0d0 + 419.0d0/256.0d0 * rho
     &      + 163.0d0/128.0d0 * rho2 + 119.0d0/192.0d0 * rho3
     &      + 5.0d0/24.0d0 * rho4 + rho5 / 20.0d0 + rho6 / 120.0d0
     &      + rho7 / 1260.0d0) * exp2)
         coul2 = a / rho3 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &      + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4 + 191.0d0/720.0d0
     &      * rho5 + 31.0d0/360.0d0 * rho6 + 19.0d0/840.0d0 * rho7
     &      + 17.0d0/3780.0d0 * rho8 + rho9 / 1890.0d0) * exp2)
         coul3 = coul2
         coul4 = 6.0d0 * a / rho5 * (1.0d0 - (1.0d0 + 2.0d0 * rho
     &      + 2.0d0 * rho2 + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4
     &      + 511.0d0/1920.0d0 * rho5 + 253.0d0/2880.0d0 * rho6
     &      + 2237.0d0/90720.0d0 * rho7 + 71.0d0/11340.0d0 * rho8
     &      + rho9 / 630.0d0 + 13.0d0/34020.0d0 * rho10
     &      + rho11 / 17010.0d0) * exp2)
         c = b
         b = a
         a = c
      else
         c = b
         b = a
         a = c
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rho5 = rho**5
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoA6 = rhoA5 * rhoA
         rhoA7 = rhoA6 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         rhoB6 = rhoB5 * rhoB
         rhoB7 = rhoB6 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappa4 = kappa3 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappap4 = kappap3 * kappap
         kappam3 = kappam**3
         kappam4 = kappam3 * kappam
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 1.0d0 / r
         term1 = -kappam3 * (1.0d0/16.0d0 * (8.0d0 - kappa - 27.0d0
     &      * kappa2 - 30.0d0 * kappa3 - 10.0d0 * kappa4) + 1.0d0/32.0d0
     &      * (11.0d0 - 19.0d0 * kappa - 44.0d0 * kappa2 - 20.0d0
     &      * kappa3) * rhoA + 1.0d0/16.0d0 * (1.0d0 - 5.0d0 * kappa
     &      - 4.0d0 * kappa2) * rhoA2 - 1.0d0/24.0d0 * kappa * rhoA3)
     &      * exp2A
         term2 = -kappap3 * (1.0d0/16.0d0 * (8.0d0 + kappa - 27.0d0
     &      * kappa2 + 30.0d0 * kappa3 - 10.0d0 * kappa4) + 1.0d0/32.0d0
     &      * (11.0d0 + 19.0d0 * kappa - 44.0d0 * kappa2 + 20.0d0
     &      * kappa3) * rhoB + 1.0d0/16.0d0 * (1.0d0 + 5.0d0 * kappa
     &      - 4.0d0 * kappa2) * rhoB2 + 1.0d0/24.0d0 * kappa * rhoB3)
     &      * exp2B
         coul1 = pre * (1.0d0 + term1 + term2)
         pre = alpha / ((1.0d0 - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/32.0d0 * (21.0d0 + 44.0d0 * kappa
     &      + 35.0d0 * kappa2 + 10.0d0 * kappa3) * (1.0d0 + 2.0d0
     &      * rhoA) + 1.0d0/24.0d0 * (29.0d0 + 56.0d0 * kappa + 40.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA2 + 1.0d0/12.0d0
     &      * (8.0d0 + 12.0d0 * kappa + 5.0d0 * kappa2) * rhoA3
     &      + 1.0d0/24.0d0 * (5.0d0 + 4.0d0 * kappa) * rhoA4
     &      + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (11.0d0 + 7.0d0 * kappa
     &      - 39.0d0 * kappa2 + 35.0d0 * kappa3 - 10.0d0 * kappa4)
     &      * (1.0d0 + 2.0d0 * rhoB) + 1.0d0/24.0d0 * (14.0d0 + 13.0d0
     &      * kappa - 51.0d0 * kappa2 + 40.0d0 * kappa3 - 10.0d0
     &      * kappa4) * rhoB2 + 1.0d0/12.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * rhoB3 + 1.0d0/24.0d0
     &      * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2) * rhoB4
     &      + 1.0d0/36.0d0 * kappa * rhoB5) * exp2B
         coul2 = pre * (1.0d0 + term1 + term2)
         pre = 6. * alpha / ((1.0d0 + tau)**2 * (1.0d0 - tau)**2 * rho5)
         term1 = -kappam4 * (1.0d0/32.0d0 * (16.0d0 + 29.0d0 * kappa
     &      + 20.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/144.0d0 * (139.0d0 + 246.0d0 * kappa + 165.0d0
     &      * kappa2 + 40.0d0 * kappa3) * rhoA2 + 1.0d0/72.0d0
     &      * (43.0d0 + 72.0d0 * kappa + 45.0d0 * kappa2 + 10.0d0
     &      * kappa3) * rhoA3 + 1.0d0/216.0d0 * (57.0d0 + 88.0d0 * kappa
     &      + 50.0d0 * kappa2 + 10.0d0 * kappa3) * rhoA4 + 1.0d0/432.0d0
     &      * (39.0d0 + 52.0d0 * kappa + 20.0d0 * kappa2) * rhoA5
     &      + 1.0d0/216.0d0 * (5.0d0 + 4.0d0 * kappa) * rhoA6
     &      + rhoA7 / 324.0d0) * exp2A
         term2 = -kappap4 * (1.0d0/32.0d0 * (16.0d0 - 29.0d0 * kappa
     &      + 20.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/144.0d0 * (139.0d0 - 246.0d0 * kappa + 165.0d0
     &      * kappa2 - 40.0d0 * kappa3) * rhoB2 + 1.0d0/72.0d0
     &      * (43.0d0 - 72.0d0 * kappa + 45.0d0 * kappa2 - 10.0d0
     &      * kappa3) * rhoB3 + 1.0d0/216.0d0 * (57.0d0 - 88.0d0 * kappa
     &      + 50.0d0 * kappa2 - 10.0d0 * kappa3) * rhoB4 + 1.0d0/432.0d0
     &      * (39.0d0 - 52.0d0 * kappa + 20.0d0 * kappa2) * rhoB5
     &      + 1.0d0/216.0d0 * (5.0d0 - 4.0d0 * kappa) * rhoB6
     &      + rhoB7 / 324.0d0) * exp2B
         coul4 = pre * (1.0d0 + term1 + term2)
         c = b
         b = a
         a = c
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappa4 = kappa3 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam4 = kappam**4
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = alpha / ((1.0d0 - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/32.0d0 * (21.0d0 + 44.0d0 * kappa
     &      + 35.0d0 * kappa2 + 10.0d0 * kappa3) * (1.0d0 + 2.0d0
     &      * rhoA) + 1.0d0/24.0d0 * (29.0d0 + 56.0d0 * kappa + 40.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA2 + 1.0d0/12.0d0
     &      * (8.0d0 + 12.0d0 * kappa + 5.0d0 * kappa2) * rhoA3
     &      + 1.0d0/24.0d0 * (5.0d0 + 4.0d0 * kappa) * rhoA4
     &      + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (11.0d0 + 7.0d0 * kappa
     &      - 39.0d0 * kappa2 + 35.0d0 * kappa3 - 10.0d0 * kappa4)
     &      * (1.0d0 + 2.0d0 * rhoB) + 1.0d0/24.0d0 * (14.0d0 + 13.0d0
     &      * kappa - 51.0d0 * kappa2 + 40.0d0 * kappa3 - 10.0d0
     &      * kappa4) * rhoB2 + 1.0d0/12.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * rhoB3 + 1.0d0/24.0d0
     &      * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2) * rhoB4
     &      + 1.0d0/36.0d0 * kappa * rhoB5) * exp2B
         coul3 = pre * (1.0d0 + term1 + term2)
      end if
      coulombPxPxPzPz = coul1 - 1.5 d0* coul2
     &                + 3.0d0 * coul3 - 4.5d0 * coul4
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function coulombPxPxPxPx  --  PxPxPxPx coulomb integral  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "coulombPxPxPxPx" computes 2-electron coulomb integral
c
c
      function coulombPxPxPxPx (a, b, z1, z2)
      implicit none
      real*8 coulombPxPxPxPx,coul1,coul2,coul3,coul4,coul5
      real*8 a,b,z1,z2
      real*8 c
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5,rhoA6,rhoA7
      real*8 rhoB2,rhoB3,rhoB4,rhoB5,rhoB6,rhoB7
      real*8 kappa2,kappa3,kappa4
      real*8 kappap,kappap3,kappap4
      real*8 kappam,kappam3,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         rho8 = rho7 * rho
         rho9 = rho8 * rho
         rho10 = rho9 * rho
         rho11 = rho10 * rho
         exp2 = exp(-2.0d0 * rho)
         coul1 = 1.0d0/ r * (1.0d0 - (1.0d0 + 419.0d0/256.0d0 * rho
     &      + 163.0d0/128.0d0 * rho2 + 119.0d0/192.0d0 * rho3
     &      + 5.0d0/24.0d0 * rho4 + rho5 / 20.0d0 + rho6 / 120.0d0
     &      + rho7 / 1260.0d0) * exp2)
         coul2 = a / rho3 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &      + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4 + 191.0d0/720.0d0
     &      * rho5 + 31.0d0/360.0d0 * rho6 + 19.0d0/840.0d0 * rho7
     &      + 17.0d0/3780.0d0 * rho8 + rho9 / 1890.0d0) * exp2)
         coul3 = coul2
         coul4 = 6.0d0 * a / rho5 * (1.0d0 - (1.0d0 + 2.0d0 * rho
     &      + 2.0d0 * rho2 + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4
     &      + 511.0d0/1920.0d0 * rho5 + 253.0d0/2880.0d0 * rho6
     &      + 2237.0d0/90720.0d0 * rho7 + 71.0d0/11340.0d0 * rho8
     &      + rho9 / 630.0d0 + 13.0d0/34020.0d0 * rho10
     &      + rho11 / 17010.0d0) * exp2)
         coul5 = a / rho5 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &      + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4 + 253.0d0/960.0d0
     &      * rho5 + 119.0d0/1440.0d0 * rho6 + 11.0d0/560.0d0 * rho7
     &      + rho8 / 315.0d0 + rho9 / 3780.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rho5 = rho**5
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoA6 = rhoA5 * rhoA
         rhoA7 = rhoA6 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         rhoB6 = rhoB5 * rhoB
         rhoB7 = rhoB6 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappa4 = kappa3 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappap4 = kappap3 * kappap
         kappam3 = kappam**3
         kappam4 = kappam3 * kappam
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 1.0d0 / r
         term1 = -kappam3 * (1.0d0/16.0d0 * (8.0d0 - kappa - 27.0d0
     &      * kappa2 - 30.0d0 * kappa3 - 10.0d0 * kappa4) + 1.0d0/32.0d0
     &      * (11.0d0 - 19.0d0 * kappa - 44.0d0 * kappa2 - 20.0d0
     &      * kappa3) * rhoA + 1.0d0/16.0d0 * (1.0d0 - 5.0d0 * kappa
     &      - 4.0d0 * kappa2) * rhoA2 - 1.0d0/24.0d0 * kappa * rhoA3)
     &      * exp2A
         term2 = -kappap3 * (1.0d0/16.0d0 * (8.0d0 + kappa - 27.0d0
     &      * kappa2 + 30.0d0 * kappa3 - 10.0d0 * kappa4) + 1.0d0/32.0d0
     &      * (11.0d0 + 19.0d0 * kappa - 44.0d0 * kappa2 + 20.0d0
     &      * kappa3) * rhoB + 1.0d0/16.0d0 * (1.0d0 + 5.0d0 * kappa
     &      - 4.0d0 * kappa2) * rhoB2 + 1.0d0/24.0d0 * kappa * rhoB3)
     &      * exp2B
         coul1 = pre * (1.0d0 + term1 + term2)
         pre = alpha / ((1.0d0 - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/32.0d0 * (21.0d0 + 44.0d0 * kappa
     &      + 35.0d0 * kappa2 + 10.0d0 * kappa3) * (1.0d0 + 2.0d0
     &      * rhoA) + 1.0d0/24.0d0 * (29.0d0 + 56.0d0 * kappa + 40.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA2 + 1.0d0/12.0d0 * (8.0d0
     *      + 12.0d0 * kappa + 5.0d0 * kappa2) * rhoA3 + 1.0d0/24.0d0
     &      * (5.0d0 + 4.0d0 * kappa) * rhoA4 + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (11.0d0 + 7.0d0 * kappa
     &      - 39.0d0 * kappa2 + 35.0d0 * kappa3 - 10.0d0 * kappa4)
     &      * (1.0d0 + 2.0d0 * rhoB) + 1.0d0/24.0d0 * (14.0d0 + 13.0d0
     &      * kappa - 51.0d0 * kappa2 + 40.0d0 * kappa3 - 10.0d0
     *      * kappa4) * rhoB2 + 1.0d0/12.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * rhoB3 + 1.0d0/24.0d0
     &      * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2) * rhoB4
     &      + 1.0d0/36.0d0 * kappa * rhoB5) * exp2B
         coul2 = pre * (1.0d0 + term1 + term2)
         pre = 6.0d0*alpha/((1.0d0 + tau)**2 * (1.0d0 - tau)**2 * rho5)
         term1 = -kappam4 * (1.0d0/32.0d0 * (16.0d0 + 29.0d0 * kappa
     &      + 20.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/144.0d0 * (139.0d0 + 246.0d0 * kappa + 165.0d0
     &      * kappa2 + 40.0d0 * kappa3) * rhoA2 + 1.0d0/72.0d0 * (43.0d0
     &      + 72.0d0 * kappa + 45.0d0 * kappa2 + 10.0d0 * kappa3)
     &      * rhoA3 + 1.0d0/216.0d0 * (57.0d0 + 88.0d0 * kappa + 50.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA4 + 1.0d0/432.0d0
     &      * (39.0d0 + 52.0d0 * kappa + 20.0d0 * kappa2) * rhoA5
     &      + 1.0d0/216.0d0 * (5.0d0 + 4.0d0 * kappa) * rhoA6
     &      + rhoA7 / 324.0d0) * exp2A
         term2 = -kappap4 * (1.0d0/32.0d0 * (16.0d0 - 29.0d0 * kappa
     &      + 20.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/144.0d0 * (139.0d0 - 246.0d0 * kappa + 165.0d0
     &      * kappa2 - 40.0d0 * kappa3) * rhoB2 + 1.0d0/72.0d0 * (43.0d0
     &      - 72.0d0 * kappa + 45.0d0 * kappa2 - 10.0d0 * kappa3)
     &      * rhoB3 + 1.0d0/216.0d0 * (57.0d0 - 88.0d0 * kappa + 50.0d0
     &      * kappa2 - 10.0d0 * kappa3) * rhoB4 + 1.0d0/432.0d0
     &      * (39.0d0 - 52.0d0 * kappa + 20.0d0 * kappa2) * rhoB5
     &      + 1.0d0/216.0d0 * (5.0d0 - 4.0d0 * kappa) * rhoB6
     &      + rhoB7 / 324.0d0) * exp2B
         coul4 = pre * (1.0d0 + term1 + term2)
         pre = alpha / ((1.0d0 + tau)**2 * (1.0d0 - tau)**2 * rho5)
         term1 = -kappam4 * (1.0d0/32.0d0 * (16.0d0 + 29.0d0 * kappa
     &      + 20.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/48.0d0 * (43.0d0 + 72.0d0 * kappa + 45.0d0 * kappa2
     &      + 10.0d0 * kappa3) * rhoA2 + 1.0d0/24.0d0 * (11.0d0 + 14.0d0
     &      * kappa + 5.0d0 * kappa2) * rhoA3 + 1.0d0/24.0d0 * (3.0d0
     &      + 2.0d0 * kappa) * rhoA4 + rhoA5 / 72.0d0) * exp2A
         term2 = -kappap4 * (1.0d0/32.0d0 * (16.0d0 - 29.0d0 * kappa
     &      + 20.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/48.0d0 * (43.0d0 - 72.0d0 * kappa + 45.0d0 * kappa2
     &      - 10.0d0 * kappa3) * rhoB2 + 1.0d0/24.0d0 * (11.0d0 - 14.0d0
     &      * kappa + 5.0d0 * kappa2) * rhoB3 + 1.0d0/24.0d0 * (3.0d0
     &      - 2.0d0 * kappa) * rhoB4 + rhoB5 / 72.0d0) * exp2B
         coul5 = pre * (1.0d0 + term1 + term2)
         c = b
         b = a
         a = c
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappa4 = kappa3 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam4 = kappam**4
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = alpha / ((1.0d0 - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/32.0d0 * (21.0d0 + 44.0d0 * kappa
     &      + 35.0d0 * kappa2 + 10.0d0 * kappa3) * (1.0d0 + 2.0d0
     &      * rhoA) + 1.0d0/24.0d0 * (29.0d0 + 56.0d0 * kappa + 40.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA2 + 1.0d0/12.0d0 * (8.0d0
     *      + 12.0d0 * kappa + 5.0d0 * kappa2) * rhoA3 + 1.0d0/24.0d0
     &      * (5.0d0 + 4.0d0 * kappa) * rhoA4 + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (11.0d0 + 7.0d0 * kappa
     &      - 39.0d0 * kappa2 + 35.0d0 * kappa3 - 10.0d0 * kappa4)
     &      * (1.0d0 + 2.0d0 * rhoB) + 1.0d0/24.0d0 * (14.0d0 + 13.0d0
     &      * kappa - 51.0d0 * kappa2 + 40.0d0 * kappa3 - 10.0d0
     *      * kappa4) * rhoB2 + 1.0d0/12.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * rhoB3 + 1.0d0/24.0d0
     &      * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2) * rhoB4
     &      + 1.0d0/36.0d0 * kappa * rhoB5) * exp2B
         coul3 = pre * (1.0d0 + term1 + term2)
         c = b
         b = a
         a = c
      end if
      coulombPxPxPxPx = coul1 - 1.5d0 * (coul2 + coul3)
     &                + 9.0d0/4.0d0 * coul4 + 27.0d0/4.0d0 * coul5
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function coulombPxPxPyPy  --  PxPxPyPy coulomb integral  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "coulombPxPxPyPy" computes 2-electron coulomb integral
c
c
      function coulombPxPxPyPy (a, b, z1, z2)
      implicit none
      real*8 coulombPxPxPyPy,coul1,coul2,coul3,coul4,coul5
      real*8 a,b,z1,z2
      real*8 c
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5,rhoA6,rhoA7
      real*8 rhoB2,rhoB3,rhoB4,rhoB5,rhoB6,rhoB7
      real*8 kappa2,kappa3,kappa4
      real*8 kappap,kappap3,kappap4
      real*8 kappam,kappam3,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         rho8 = rho7 * rho
         rho9 = rho8 * rho
         rho10 = rho9 * rho
         rho11 = rho10 * rho
         exp2 = exp(-2.0d0 * rho)
         coul1 = 1.0d0/ r * (1.0d0 - (1.0d0 + 419.0d0/256.0d0 * rho
     &      + 163.0d0/128.0d0 * rho2 + 119.0d0/192.0d0 * rho3
     &      + 5.0d0/24.0d0 * rho4 + rho5 / 20.0d0 + rho6 / 120.0d0
     &      + rho7 / 1260.0d0) * exp2)
         coul2 = a / rho3 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &      + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4 + 191.0d0/720.0d0
     &      * rho5 + 31.0d0/360.0d0 * rho6 + 19.0d0/840.0d0 * rho7
     &      + 17.0d0/3780.0d0 * rho8 + rho9 / 1890.0d0) * exp2)
         coul3 = coul2
         coul4 = 6.0d0 * a / rho5 * (1.0d0 - (1.0d0 + 2.0d0 * rho
     &      + 2.0d0 * rho2 + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4
     &      + 511.0d0/1920.0d0 * rho5 + 253.0d0/2880.0d0 * rho6
     &      + 2237.0d0/90720.0d0 * rho7 + 71.0d0/11340.0d0 * rho8
     &      + rho9 / 630.0d0 + 13.0d0/34020.0d0 * rho10
     &      + rho11 / 17010.0d0) * exp2)
         coul5 = a / rho5 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &      + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4 + 253.0d0/960.0d0
     &      * rho5 + 119.0d0/1440.0d0 * rho6 + 11.0d0/560.0d0 * rho7
     &      + rho8 / 315.0d0 + rho9 / 3780.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rho5 = rho**5
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoA6 = rhoA5 * rhoA
         rhoA7 = rhoA6 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         rhoB6 = rhoB5 * rhoB
         rhoB7 = rhoB6 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappa4 = kappa3 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappap4 = kappap3 * kappap
         kappam3 = kappam**3
         kappam4 = kappam3 * kappam
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = 1.0d0 / r
         term1 = -kappam3 * (1.0d0/16.0d0 * (8.0d0 - kappa - 27.0d0
     &      * kappa2 - 30.0d0 * kappa3 - 10.0d0 * kappa4) + 1.0d0/32.0d0
     &      * (11.0d0 - 19.0d0 * kappa - 44.0d0 * kappa2 - 20.0d0
     &      * kappa3) * rhoA + 1.0d0/16.0d0 * (1.0d0 - 5.0d0 * kappa
     &      - 4.0d0 * kappa2) * rhoA2 - 1.0d0/24.0d0 * kappa * rhoA3)
     &      * exp2A
         term2 = -kappap3 * (1.0d0/16.0d0 * (8. 0d0+ kappa - 27.0d0
     &      * kappa2 + 30.0d0 * kappa3 - 10.0d0 * kappa4) + 1.0d0/32.0d0
     &      * (11.0d0 + 19.0d0 * kappa - 44.0d0 * kappa2 + 20.0d0
     &      * kappa3) * rhoB + 1.0d0/16.0d0 * (1.0d0 + 5.0d0 * kappa
     &      - 4.0d0 * kappa2) * rhoB2 + 1.0d0/24.0d0 * kappa * rhoB3)
     &      * exp2B
         coul1 = pre * (1.0d0 + term1 + term2)
         pre = alpha / ((1.0d0 - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/32.0d0 * (21.0d0 + 44.0d0 * kappa
     &      + 35.0d0 * kappa2 + 10.0d0 * kappa3) * (1.0d0 + 2.0d0
     &      * rhoA) + 1.0d0/24.0d0 * (29.0d0 + 56.0d0 * kappa + 40.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA2 + 1.0d0/12.0d0 * (8.0d0
     &      + 12.0d0 * kappa + 5.0d0 * kappa2) * rhoA3 + 1.0d0/24.0d0
     &      * (5.0d0 + 4.0d0 * kappa) * rhoA4 + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (11.0d0 + 7.0d0 * kappa
     &      - 39.0d0 * kappa2 + 35.0d0 * kappa3 - 10.0d0 * kappa4)
     &      * (1.0d0 + 2.0d0 * rhoB) + 1.0d0/24.0d0 * (14.0d0 + 13.0d0
     &      * kappa - 51.0d0 * kappa2 + 40.0d0 * kappa3 - 10.0d0
     &      * kappa4) * rhoB2 + 1.0d0/12.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * rhoB3 + 1.0d0/24.0d0
     &      * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2) * rhoB4
     &      + 1.0d0/36.0d0 * kappa * rhoB5) * exp2B
         coul2 = pre * (1.0d0 + term1 + term2)
         pre = 6.0d0*alpha/((1.0d0 + tau)**2 * (1.0d0 - tau)**2 * rho5)
         term1 = -kappam4 * (1.0d0/32.0d0 * (16.0d0 + 29.0d0 * kappa
     &      + 20.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/144.0d0 * (139.0d0 + 246.0d0 * kappa + 165.0d0
     &      * kappa2 + 40.0d0 * kappa3) * rhoA2 + 1.0d0/72.0d0 * (43.0d0
     &      + 72.0d0 * kappa + 45.0d0 * kappa2 + 10.0d0 * kappa3)
     &      * rhoA3 + 1.0d0/216.0d0 * (57.0d0 + 88.0d0 * kappa + 50.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA4 + 1.0d0/432.0d0
     &      * (39.0d0 + 52.0d0 * kappa + 20.0d0 * kappa2) * rhoA5
     &      + 1.0d0/216.0d0 * (5.0d0 + 4.0d0 * kappa) * rhoA6
     &      + rhoA7 / 324.0d0) * exp2A
         term2 = -kappap4 * (1.0d0/32.0d0 * (16.0d0 - 29.0d0 * kappa
     &      + 20.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/144.0d0 * (139.0d0 - 246.0d0 * kappa + 165.0d0
     &      * kappa2 - 40.0d0 * kappa3) * rhoB2 + 1.0d0/72.0d0 * (43.0d0
     &      - 72.0d0 * kappa + 45.0d0 * kappa2 - 10.0d0 * kappa3)
     &      * rhoB3 + 1.0d0/216.0d0 * (57.0d0 - 88.0d0 * kappa + 50.0d0
     &      * kappa2 - 10.0d0 * kappa3) * rhoB4 + 1.0d0/432.0d0
     &      * (39.0d0 - 52.0d0 * kappa + 20.0d0 * kappa2) * rhoB5
     &      + 1.0d0/216.0d0 * (5.0d0 - 4.0d0 * kappa) * rhoB6
     &      + rhoB7 / 324.0d0) * exp2B
         coul4 = pre * (1.0d0 + term1 + term2)
         pre = alpha / ((1.0d0 + tau)**2 * (1.0d0 - tau)**2 * rho5)
         term1 = -kappam4 * (1.0d0/32.0d0 * (16.0d0 + 29.0d0 * kappa
     &      + 20.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/48.0d0 * (43.0d0 + 72.0d0 * kappa + 45.0d0 * kappa2
     &      + 10.0d0 * kappa3) * rhoA2 + 1.0d0/24.0d0 * (11.0d0 + 14.0d0
     &      * kappa + 5.0d0 * kappa2) * rhoA3 + 1.0d0/24.0d0 * (3.0d0
     &      + 2.0d0 * kappa) * rhoA4 + rhoA5 / 72.0d0) * exp2A
         term2 = -kappap4 * (1.0d0/32.0d0 * (16.0d0 - 29.0d0 * kappa
     &      + 20.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/48.0d0 * (43.0d0 - 72.0d0 * kappa + 45.0d0 * kappa2
     &      - 10.0d0 * kappa3) * rhoB2 + 1.0d0/24.0d0 * (11.0d0 - 14.0d0
     &      * kappa + 5.0d0 * kappa2) * rhoB3 + 1.0d0/24.0d0 * (3.0d0
     &      - 2.0d0 * kappa) * rhoB4 + rhoB5 / 72.0d0) * exp2B
         coul5 = pre * (1.0d0 + term1 + term2)
         c = b
         b = a
         a = c
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappa4 = kappa3 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap3 = kappap**3
         kappam4 = kappam**4
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = alpha / ((1.0d0 - tau)**2 * rho3)
         term1 = -kappam4 * (1.0d0/32.0d0 * (21.0d0 + 44.0d0 * kappa
     &      + 35.0d0 * kappa2 + 10.0d0 * kappa3) * (1.0d0 + 2.0d0
     &      * rhoA) + 1.0d0/24.0d0 * (29.0d0 + 56.0d0 * kappa + 40.0d0
     &      * kappa2 + 10.0d0 * kappa3) * rhoA2 + 1.0d0/12.0d0 * (8.0d0
     &      + 12.0d0 * kappa + 5.0d0 * kappa2) * rhoA3 + 1.0d0/24.0d0
     &      * (5.0d0 + 4.0d0 * kappa) * rhoA4 + rhoA5 / 36.0d0) * exp2A
         term2 = -kappap3 * (1.0d0/32.0d0 * (11.0d0 + 7.0d0 * kappa
     &      - 39.0d0 * kappa2 + 35.0d0 * kappa3 - 10.0d0 * kappa4)
     &      * (1.0d0 + 2.0d0 * rhoB) + 1.0d0/24.0d0 * (14.0d0 + 13.0d0
     &      * kappa - 51.0d0 * kappa2 + 40.0d0 * kappa3 - 10.0d0
     &      * kappa4) * rhoB2 + 1.0d0/12.0d0 * (3.0d0 + 6.0d0 * kappa
     &      - 12.0d0 * kappa2 + 5.0d0 * kappa3) * rhoB3 + 1.0d0/24.0d0
     &      * (1.0d0 + 5.0d0 * kappa - 4.0d0 * kappa2) * rhoB4
     &      + 1.0d0/36.0d0 * kappa * rhoB5) * exp2B
         coul3 = pre * (1.0d0 + term1 + term2)
         c = b
         b = a
         a = c
      end if
      coulombPxPxPyPy = coul1 - 1.5d0 * (coul2 + coul3)
     &                + 9.0d0/4.0d0 * coul4 - 27.0d0/4.0d0 * coul5
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function coulombPxPyPxPy  --  PxPyPxPy coulomb integral  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "coulombPxPyPxPy" computes 2-electron coulomb integral
c
c
      function coulombPxPyPxPy (a, b, z1, z2)
      implicit none
      real*8 coulombPxPyPxPy,coul
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9
      real*8 exp2
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5
      real*8 rhoB2,rhoB3,rhoB4,rhoB5
      real*8 kappa2,kappa3
      real*8 kappap,kappap4
      real*8 kappam,kappam4
      real*8 exp2A,exp2B
c
c
      diff = abs(a - b)
      eps = 0.01d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         rho5 = rho4 * rho
         rho6 = rho5 * rho
         rho7 = rho6 * rho
         rho8 = rho7 * rho
         rho9 = rho8 * rho
         exp2 = exp(-2.0d0 * rho)
         coul = a / rho5 * (1.0d0 - (1.0d0 + 2.0d0 * rho + 2.0d0 * rho2
     &      + 4.0d0/3.0d0 * rho3 + 2.0d0/3.0d0 * rho4 + 253.0d0/960.0d0
     &      * rho5 + 119.0d0/1440.0d0 * rho6 + 11.0d0/560.0d0 * rho7
     &      + rho8 / 315.0d0 + rho9 / 3780.0d0) * exp2)
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho5 = rho**5
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoA5 = rhoA4 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         rhoB5 = rhoB4 * rhoB
         kappa2 = kappa * kappa
         kappa3 = kappa2 * kappa
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappap4 = kappap**4
         kappam4 = kappam**4
         exp2A = exp(-2.0d0 * rhoA)
         exp2B = exp(-2.0d0 * rhoB)
         pre = alpha / ((1.0d0 + tau)**2 * (1.0d0 - tau)**2 * rho5)
         term1 = -kappam4 * (1.0d0/32.0d0 * (16.0d0 + 29.0d0 * kappa
     &      + 20.0d0 * kappa2 + 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoA)
     &      + 1.0d0/48.0d0 * (43.0d0 + 72.0d0 * kappa + 45.0d0 * kappa2
     &      + 10.0d0 * kappa3) * rhoA2 + 1.0d0/24.0d0 * (11.0d0 + 14.0d0
     &      * kappa + 5.0d0 * kappa2) * rhoA3 + 1.0d0/24.0d0 * (3.0d0
     &      + 2.0d0 * kappa) * rhoA4 + rhoA5 / 72.0d0) * exp2A
         term2 = -kappap4 * (1.0d0/32.0d0 * (16.0d0 - 29.0d0 * kappa
     &      + 20.0d0 * kappa2 - 5.0d0 * kappa3) * (1.0d0 + 2.0d0 * rhoB)
     &      + 1.0d0/48.0d0 * (43.0d0 - 72.0d0 * kappa + 45.0d0 * kappa2
     &      - 10.0d0 * kappa3) * rhoB2 + 1.0d0/24.0d0 * (11.0d0 - 14.0d0
     &      * kappa + 5.0d0 * kappa2) * rhoB3 + 1.0d0/24.0d0 * (3.0d0
     &      - 2.0d0 * kappa) * rhoB4 + rhoB5 / 72.0d0) * exp2B
         coul = pre * (1. + term1 + term2)
      end if
      coulombPxPyPxPy = 27.0d0/4.0d0 * coul
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function coulombSSSPz  --  SSSPz coulomb integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "coulombSSSPz" computes 2-electron coulomb integral
c
c
      function coulombSSSPz (a, b, z1, z2)
      implicit none
      real*8 coulombSSSPz,coulombSSPzS
      real*8 a,b,z1,z2
c
c
      coulombSSSPz = coulombSSPzS(a,b,z1,z2)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  function coulombSPzSS  --  SPzSS coulomb integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "coulombSPzSS" computes 2-electron coulomb integral
c
c
      function coulombSPzSS (a, b, z1, z2)
      implicit none
      real*8 coulombSPzSS,coulombSSSPz
      real*8 a,b,z1,z2
c
c
      coulombSPzSS = -coulombSSSPz(b, a, z2, z1)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  function coulombPzPzSS  --  PzPzSS coulomb integral  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "coulombPzPzSS" computes 2-electron coulomb integral
c
c
      function coulombPzPzSS (a, b, z1, z2)
      implicit none
      real*8 coulombPzPzSS,coulombSSPzPz
      real*8 a,b,z1,z2
c
c
      coulombPzPzSS = coulombSSPzPz(b, a, z2, z1)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  function coulombPxPxSS  --  PxPxSS coulomb integral  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "coulombPxPxSS" computes 2-electron coulomb integral
c
c
      function coulombPxPxSS (a, b, z1, z2)
      implicit none
      real*8 coulombPxPxSS,coulombSSPxPx
      real*8 a,b,z1,z2
c
c
      coulombPxPxSS = coulombSSPxPx(b, a, z2, z1)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function coulombPzPzSPz  --  PzPzSPz coulomb integral  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "coulombPzPzSPz" computes 2-electron coulomb integral
c
c
      function coulombPzPzSPz (a, b, z1, z2)
      implicit none
      real*8 coulombPzPzSPz,coulombSPzPzPz
      real*8 a,b,z1,z2
c
c
      coulombPzPzSPz = -coulombSPzPzPz(b, a, z2, z1)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function coulombPxPzSPx  --  PxPzSPx coulomb integral  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "coulombPxPzSPx" computes 2-electron coulomb integral
c
c
      function coulombPxPzSPx (a, b, z1, z2)
      implicit none
      real*8 coulombPxPzSPx,coulombSPxPxPz
      real*8 a,b,z1,z2
c
c
      coulombPxPzSPx = -coulombSPxPxPz(b, a, z2, z1)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function coulombPxPxSPz  --  PxPxSPz coulomb integral  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "coulombPxPxSPz" computes 2-electron coulomb integral
c
c
      function coulombPxPxSPz (a, b, z1, z2)
      implicit none
      real*8 coulombPxPxSPz,coulombSPzPxPx
      real*8 a,b,z1,z2
c
c
      coulombPxPxSPz = -coulombSPzPxPx(b, a, z2, z1)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function coulombPzPzPxPx  --  PzPzPxPx coulomb integral  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "coulombPzPzPxPx" computes 2-electron coulomb integral
c
c
      function coulombPzPzPxPx (a, b, z1, z2)
      implicit none
      real*8 coulombPzPzPxPx,coulombPxPxPzPz
      real*8 a,b,z1,z2
c
c
      coulombPzPzPxPx = coulombPxPxPzPz(b, a, z2, z1)
      return
      end
