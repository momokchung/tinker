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
c
c
      module kxrepl
      implicit none
      real*8, allocatable :: pxrz(:)
      real*8, allocatable :: pxrdmp(:)
      real*8, allocatable :: pxrele(:)
      save
      end
