c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module xrepel  --  exch repulsion for current structure  ##        
c     ##                                                           ##
c     ###############################################################
c
c
c     nsto       number of primitives for STO-nG
c     n2sto      nsto * nsto
c     ncoeff     number of Chebyshev coefficients for Boys Function
c     zpxr       nuclear charge parameter value at each site
c     dmppxr     exchange repulsion alpha damping value at each site
c     elepxr     electron charge parameter value at each site
c     stocoeff   coefficients for STO-nG basis
c     stoexp     exponents for STO-nG basis
c     boysCoeff  Boys function Chebyshev coefficients
c     xreptyp    exchange repulsion type
c     stong      STO-nG
c
c
      module xrepel
      implicit none
      integer nsto
      integer n2sto
      integer ncoeff
      real*8, allocatable :: zpxr(:)
      real*8, allocatable :: dmppxr(:)
      real*8, allocatable :: elepxr(:)
      real*8, allocatable :: stocoeff(:,:)
      real*8, allocatable :: stoexp(:,:)
      real*8, allocatable :: boysCoeff(:)
      character*7 xreptyp
      character*6 stong
      save
      end
