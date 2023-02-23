c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module kxrepl  --  exch repulsion forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c

c     pxrz       nuclear charge parameter value for each atom class
c     pxrdmp     exch repulsion alpha parameter for each atom class
c     pxrele     electron charge parameter value for each atom class
c     pxrcr      ratio of p/s orbital cofficients for each atom class
c     pxrdr      ratio of p/s orbital exponents for each atom class
c
c
      module kxrepl
      implicit none
      real*8, allocatable :: pxrz(:)
      real*8, allocatable :: pxrdmp(:)
      real*8, allocatable :: pxrele(:)
      real*8, allocatable :: pxrcr(:)
      real*8, allocatable :: pxrdr(:)
      save
      end
