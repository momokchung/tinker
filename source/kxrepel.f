c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kxrepel  --  exch repulsion term assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kxrepel" assigns the nuclear charge parameter and exponential
c     parameter for exchange repulsion interaction and processes any
c     new or changed values for these parameters
c
c
      subroutine kxrepel
      use atomid
      use atoms
      use inform
      use iounit
      use kxrepl
      use keys
      use mpole
      use potent
      use xrepel
      use repel
      use reppot
      use sizes
      implicit none
      integer i,k,ii
      integer iboys
      integer freeunit
      integer ia,ic,next
      real*8 zpr,ele,apr,cpr,dpr
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     exit if using 'Pauli' repulsion
c
      if (use_repuls) return
c
c     process keywords containing exch repulsion parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'XREPULSION ') then
            k = 0
            zpr = 0.0d0
            ele = 0.0d0
            apr = 0.0d0
            cpr = 0.0d0
            dpr = 1.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  zpr,ele,apr,cpr,dpr
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Exchange Repulsion',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',13x,'Core',11x,'Charge',
     &                       9x,'Damp',11x,'P/S Coeff',6x,'P/S Damp'/)
               end if
               if (k .le. maxclass) then
                  pxrz(k) = zpr
                  pxrele(k) = -abs(ele)
                  pxrdmp(k) = apr
                  pxrcr(k) = cpr
                  pxrdr(k) = dpr
                  if (.not. silent) then
                     write (iout,30)  k,zpr,ele,apr,cpr,dpr
   30                format (6x,i6,7x,5f15.4)
                  end if
               else
                  write (iout,40)
   40             format (/,' KXREPEL  --  Too many Exchange Repulsion',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(zpxr))  deallocate (zpxr)
      if (allocated(dmppxr))  deallocate (dmppxr)
      if (allocated(elepxr))  deallocate (elepxr)
      if (allocated(crpxr))  deallocate (crpxr)
      if (allocated(drpxr))  deallocate (drpxr)
      if (allocated(cpxr))  deallocate (cpxr)
      if (allocated(rcpxr))  deallocate (rcpxr)
      allocate (zpxr(n))
      allocate (dmppxr(n))
      allocate (elepxr(n))
      allocate (crpxr(n))
      allocate (drpxr(n))
      allocate (cpxr(4,n))
      allocate (rcpxr(4,n))
c
c     assign the core charge and alpha parameters 
c
      nrep = n
      do i = 1, n
         zpxr(i) = 0.0d0
         dmppxr(i) = 0.0d0
         elepxr(i) = 0.0d0
         crpxr(i) = 0.0d0
         drpxr(i) = 1.0d0
         ic = class(i)
         if (ic .ne. 0) then
            zpxr(i) = pxrz(ic)
            elepxr(i) = pxrele(ic)
            dmppxr(i) = pxrdmp(ic)
            crpxr(i) = pxrcr(ic)
            drpxr(i) = pxrdr(ic)
         end if
      end do
c
c     process keywords containing atom specific exchange repulsion
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'XREPULSION ') then
            ia = 0
            zpr = 0.0d0
            ele = 0.0d0
            apr = 0.0d0
            cpr = 0.0d0
            dpr = 1.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,zpr,ele,apr,cpr,dpr
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Exchange Repulsion Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',16x,'Core',11x,'Charge',9x,
     &                       'Damp',11x,'P/S Coeff',6x,'P/S Damp'/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,zpr,ele,apr,cpr,dpr
   60             format (6x,i6,7x,5f15.4)
               end if
               zpxr(ia) = zpr
               elepxr(ia) = -abs(ele)
               dmppxr(ia) = apr
               crpxr(ia) = cpr
               drpxr(ia) = dpr
            end if
   70       continue
         end if
      end do
c
c     condense repulsion sites to the list of multipole sites
c
      use_repuls = .true.
      if (use_repuls) then
         nrep = 0
         do ii = 1, npole
            i = ipole(ii)
            if (zpxr(i) .ne. 0)  nrep = nrep + 1
            zpxr(ii) = zpxr(i)
            dmppxr(ii) = dmppxr(i)
            elepxr(ii) = elepxr(i)
            crpxr(ii) = crpxr(i)
            drpxr(ii) = drpxr(i)
         end do
      end if
c
c     turn off the exchange repulsion potential if not used
c
      if (nrep .eq. 0)  use_repuls = .false.
c
c     change repulsion type
c
      if (use_repuls) reptyp = 'EXCHANGE'
c
c     load Boys coefficients and STO-nG coefficients and exponents
c
c      if (use_repuls) then
c         ncoeff = 271996
c         nsto = 3
c         if (stong .eq. 'STO-3G') then
c            nsto = 3
c         else if (stong .eq. 'STO-4G') then
c            nsto = 4
c         else if (stong .eq. 'STO-5G') then
c            nsto = 5
c         else if (stong .eq. 'STO-6G') then
c            nsto = 6
c         end if
c         n2sto = nsto * nsto
cc
cc     perform dynamic allocation of some global arrays
cc
c         if (allocated(boysCoeff))  deallocate (boysCoeff)
c         allocate (boysCoeff(ncoeff * 7))
c         if (allocated(stocoeff))  deallocate (stocoeff)
c         allocate (stocoeff(2,nsto))
c         if (allocated(stoexp))  deallocate (stoexp)
c         allocate (stoexp(2,nsto))
cc
cc     read in boys coefficients
cc
c         iboys = freeunit ()
c         open (unit=iboys,file="/home/kchung25/tinker/tmp/f0006.bin",
c     &      form='unformatted',access='stream',status='old',
c     &      action='read')
c         read (iboys) boysCoeff
c         close (iboys)
cc
cc     read in STO-nG coefficients and exponents
cc
c         if (nsto .eq. 3) then
c            stocoeff(1,1) = 0.444635d0
c            stocoeff(1,2) = 0.535328d0
c            stocoeff(1,3) = 0.154329d0
c            stocoeff(2,1) = 0.391957d0
c            stocoeff(2,2) = 0.607684d0
c            stocoeff(2,3) = 0.155916d0
c            stoexp(1,1) = 0.109818d0
c            stoexp(1,2) = 0.405771d0
c            stoexp(1,3) = 2.227660d0
c            stoexp(2,1) = 0.0751386d0
c            stoexp(2,2) = 0.231031d0
c            stoexp(2,3) = 0.994203d0
c         else if (nsto .eq. 4) then
c            stocoeff(1,1) = 0.291626d0
c            stocoeff(1,2) = 0.532846d0
c            stocoeff(1,3) = 0.260141d0
c            stocoeff(1,4) = 0.0567523d0
c            stocoeff(2,1) = 0.246313d0
c            stocoeff(2,2) = 0.583575d0
c            stocoeff(2,3) = 0.286379d0
c            stocoeff(2,4) = 0.0436843d0
c            stoexp(1,1) = 0.0880187d0
c            stoexp(1,2) = 0.265204d0
c            stoexp(1,3) = 0.954620d0
c            stoexp(1,4) = 5.21686d0
c            stoexp(2,1) = 0.0628104d0
c            stoexp(2,2) = 0.163541d0
c            stoexp(2,3) = 0.502989d0
c            stoexp(2,4) = 2.32350d0
c         else if (nsto .eq. 5) then
c            stocoeff(1,1) = 0.193572d0
c            stocoeff(1,2) = 0.482570d0
c            stocoeff(1,3) = 0.331816d0
c            stocoeff(1,4) = 0.113541d0
c            stocoeff(1,5) = 0.0221406d0
c            stocoeff(2,1) = 0.156828d0
c            stocoeff(2,2) = 0.510240d0
c            stocoeff(2,3) = 0.373598d0
c            stocoeff(2,4) = 0.107558d0
c            stocoeff(2,5) = 0.0125561d0
c            stoexp(1,1) = 0.0744527d0
c            stoexp(1,2) = 0.197572d0
c            stoexp(1,3) = 0.578648d0
c            stoexp(1,4) = 2.07173d0
c            stoexp(1,5) = 11.3056d0
c            stoexp(2,1) = 0.0544949d0
c            stoexp(2,2) = 0.127920d0
c            stoexp(2,3) = 0.329060d0
c            stoexp(2,4) = 1.03250d0
c            stoexp(2,5) = 5.03629d0
c         else if (nsto .eq. 6) then
c            stocoeff(1,1) = 0.130334d0
c            stocoeff(1,2) = 0.416492d0
c            stocoeff(1,3) = 0.370563d0
c            stocoeff(1,4) = 0.168538d0
c            stocoeff(1,5) = 0.0493615d0
c            stocoeff(1,6) = 0.00916360d0
c            stocoeff(2,1) = 0.101708d0
c            stocoeff(2,2) = 0.425860d0
c            stocoeff(2,3) = 0.418036d0
c            stocoeff(2,4) = 0.173897d0
c            stocoeff(2,5) = 0.0376794d0
c            stocoeff(2,6) = 0.00375970d0
c            stoexp(1,1) = 0.0651095d0
c            stoexp(1,2) = 0.158088d0
c            stoexp(1,3) = 0.407099d0
c            stoexp(1,4) = 1.18506d0
c            stoexp(1,5) = 4.23592d0
c            stoexp(1,6) = 23.1030d0
c            stoexp(2,1) = 0.0485690d0
c            stoexp(2,2) = 0.105960d0
c            stoexp(2,3) = 0.243977d0
c            stoexp(2,4) = 0.634142d0
c            stoexp(2,5) = 2.04036d0
c            stoexp(2,6) = 10.3087d0
c         end if
c      end if
      return
      end
