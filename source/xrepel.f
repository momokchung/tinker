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
c     zpxr       nuclear charge parameter value at each site
c     dmppxr     exchange repulsion alpha damping value at each site
c     elepxr     electron charge parameter value at each site
c     crpxr      ratio of p/s orbital cofficients at each site
c     cpxr       local coefficient for pseudo wavefunction at each site
c     rcpxr      global coefficient for pseudo wavefunction at each site
c     xreptyp    exchange repulsion type
c
c
      module xrepel
      implicit none
      real*8, allocatable :: zpxr(:)
      real*8, allocatable :: dmppxr(:)
      real*8, allocatable :: elepxr(:)
      real*8, allocatable :: crpxr(:)
      real*8, allocatable :: cpxr(:,:)
      real*8, allocatable :: rcpxr(:,:)
      character*7 xreptyp
      save
      end
