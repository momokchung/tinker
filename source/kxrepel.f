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
      integer freeunit
      integer ia,ic,next
      real*8 zpr,apr,cpr
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
            apr = 0.0d0
            cpr = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  zpr,apr,cpr
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Exchange Repulsion',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',15x,'Core',11x,'Damp',
     &                       8x,'P/S Ratio'/)
               end if
               if (k .le. maxclass) then
                  pxrz(k) = abs(zpr)
                  pxrdmp(k) = abs(apr)
                  pxrcr(k) = abs(cpr)
                  if (.not. silent) then
                     write (iout,30)  k,zpr,apr,cpr
   30                format (6x,i6,7x,2f15.4,f15.3)
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
      if (allocated(cpxr))  deallocate (cpxr)
      if (allocated(rcpxr))  deallocate (rcpxr)
      allocate (zpxr(n))
      allocate (dmppxr(n))
      allocate (elepxr(n))
      allocate (crpxr(n))
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
         ic = class(i)
         if (ic .ne. 0) then
            zpxr(i) = pxrz(ic)
            dmppxr(i) = pxrdmp(ic)
            elepxr(i) = pole(1,i) - zpxr(i)
            crpxr(i) = pxrcr(ic)
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
            apr = 0.0d0
            cpr = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,zpr,apr,cpr
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Exchange Repulsion Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',17x,'Core',12x,'Damp',
     &                       8x,'P/S Ratio'/)
               end if
               zpxr(ia) = abs(zpr)
               dmppxr(ia) = abs(apr)
               elepxr(ia) = pole(1,ia) - zpxr(ia)
               crpxr(ia) = abs(cpr)
               if (.not. silent) then
                  write (iout,60)  ia,zpr,apr,cpr
   60             format (6x,i6,7x,2f15.4,f15.3)
               end if
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
      return
      end
