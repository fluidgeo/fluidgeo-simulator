!================================================================================
!         programa de elementos finitos  
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
! 
!         + implementacoes de Abimael Loula
!
!         + novo projeto: modular e fortran 90, por
!         Eduardo Garcia,        bidu@lncc.br 
!         Tuane Lopes,           tuane@lncc.br
!
!         LNCC/MCT
!         Petropolis, 07.2013
!=================================================================================
!
      module mAlgMatricial

        integer              :: neq_BM, nalhs_BM, ned_BM
        integer              :: neq_B, nalhs_B, ned_B
        real*8,  allocatable :: alhs_BM(:), brhs_BM(:), dlhs_BM(:)
        real*8,  allocatable :: alhs_B(:), brhs_B(:), dlhs_B(:)
        integer, allocatable :: id_BM(:,:), idiag_BM(:), lm_BM(:,:,:)
        integer, allocatable :: id_B(:,:), idiag_B(:), lm_B(:,:,:)       
     

!funcoes e subrotinas
        public :: backns, factns, back, factor
        public :: diag, load, addnsl, addlhs, addrhs
        public :: btod, kdbc, pivots, ftod, btdb, colht

        public :: coldot, rowdot
        public :: matadd

      contains
!
!**** new **********************************************************************
!
      subroutine solverDiretoSkyLine(alhs, brhs, idiag, nalhs, neq,label, opt)
!
         implicit none
!
         integer, intent(in)   :: nalhs, neq, idiag(*), opt
         real*8, intent(inout) :: alhs(*), brhs(*)
         character(len=3) :: label
         real*8 :: t1,t2
	
	 if (opt .eq. 1) goto 111
	 if (opt .eq. 0) goto 222
111      call timing(t1)
         call factor(alhs,idiag,nalhs,neq)
         call timing(t2)
#ifdef mostrarTempos
         print*, "solver skyline: ", label, ", tempo de factorization=", t2-t1
#endif

222      call timing(t1)
         call back  (alhs,brhs, idiag,neq)
         call timing(t2)
#ifdef mostrarTempos
         print*, "solver skyline: ", label, ", tempo de backsubstitution=", t2-t1
#endif

      end subroutine solverDiretoSkyLine

!
!**** new **********************************************************************
!
       subroutine solverGaussSkyline(opt,alhs, brhs, idiag, id, nalhs, neq, ndof, numnp)

         implicit none

         integer*4, intent(in)   :: nalhs, neq, opt
         integer*4, intent(in)   :: ndof, numnp
         integer*4, intent(in)   :: idiag(neq), id(ndof, numnp)
         real*8, intent(inout) :: alhs(nalhs), brhs(neq)
         
         if (opt .eq. 1) call factor(alhs,idiag,nalhs,neq)
         call back  (alhs,brhs,idiag,neq)

      end subroutine solverGaussSkyline      
      
!
!**** new **********************************************************************
!
      subroutine diag(idiag,neq,n)
!
      implicit none
!
!.... program to compute diagonal addresses of left-hand-side matrix
!
      integer :: neq, n
      integer :: idiag(neq)
!
      integer :: i

      n = 1
      idiag(1) = 1 
      if (neq.eq.1) return
!
      do 100 i=2,neq
      idiag(i) = idiag(i) + idiag(i-1) + 1
  100 continue
      n = idiag(neq)
!
      return
      end subroutine
!
!**** new **********************************************************************
!
      subroutine load(id,f,brhs,ndof,numnp,nlvect)
!
!.... program to accumulate nodal forces and transfer into
!        right-hand-side vector
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: id(ndof,*)
      real*8  :: f(ndof,numnp,*),brhs(*)
      integer :: ndof, numnp, nlvect
!
      integer :: nlv
      integer :: i, j, k

      do i=1,ndof
         do j=1,numnp
            k = id(i,j)
            if (k.gt.0) then
               do nlv=1,nlvect
                      brhs(k) = brhs(k) + f(i,j,nlv)
                  !write(1,*) brhs(k)
               end do
            endif
         end do
      end do

      return
      end subroutine
      
!
!**** new **********************************************************************
!
      subroutine dirichletConditions(id,d,f,ndof,numnp,nlvect)
!
!.... program to compute displacement boundary conditions
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: ndof, numnp, nlvect
      integer*4:: id(ndof,*)
      real*8  :: d(ndof,*),f(ndof,numnp,*)
!
      integer*4:: i, j, k, lv
      real*8  :: val
!
!
      do 300 i=1,ndof
!
            do 200 j=1,numnp
!
            k = id(i,j)
            if (k.gt.0) go to 200
            val = 0.d0
                  do 100 lv=1,nlvect
                  val = val + f(i,j,lv)
100               continue
!
            d(i,j) = val
!
  200       continue
!
 300  continue
      return
      end subroutine

!**** new*********************************************************
      subroutine addnsl(alhs,clhs,eleffm,idiag,lm,nee,ldiag)
!
!         program to add element left-hand-side matrix to          
!                global left-hand-side matrix                      
!                                                                  
!        ldiag = .true.,  add diagonal element matrix              
!                                                                  
!        ldiag = .false, then                                     
!        add full nonsymmetric element matrix                   
!                                                                  
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: alhs(*),clhs(*),eleffm(nee,*)
      integer :: idiag(*),lm(*), nee
      logical ldiag
!
      integer :: i,j,k,l,m
!
      if (ldiag) then
!
         do 100 j=1,nee
             k = iabs(lm(j))
            if (k.gt.0) then
               l = idiag(k)
               alhs(l) = alhs(l) + eleffm(j,j)
            endif
  100    continue
!
      else
!
         do 400 j=1,nee
             k = iabs(lm(j))
            if (k.gt.0) then
               do 200 i=1,nee
                   m = iabs(lm(i))
                  if (m.gt.0) then
                     if (k.gt.m) then
                        l = idiag(k) - k + m
                        alhs(l) = alhs(l) + eleffm(i,j)
                     else
                        l = idiag(m) - m + k
                        clhs(l) = clhs(l) + eleffm(i,j)
                     endif
                     if (k.eq.m) then
                        l = idiag(k)
                        alhs(l) = alhs(l) + eleffm(i,j)
                        clhs(l) = alhs(l)
                     endif
                  endif
  200          continue
            endif
  400    continue
!
      endif
!
      return
      end subroutine

!**** new*********************************************************
      subroutine backns(a,c,b,idiag,neq)
!
!.... program to perform forward reduction and back substitution
!
      implicit none
!
!.... deactivate above card(s) for single-precision operation
!
      real*8  :: a(*),c(*),b(*)
      integer :: idiag(*)
      integer :: neq
!
      integer :: i, j, jj, jcolht, jjlast, jjnext, istart, jtemp
      real*8  :: ajj, bj
!
!.... forward reduction
!
      jj = 0
!
      do 100 j=1,neq
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1) then
           b(j) = b(j) - coldot(c(jjlast+1),b(j-jcolht+1),jcolht-1)
      endif
  100 continue
!
!.... diagonal scaling
!
      do 200 j=1,neq
      ajj = a(idiag(j))
!
!.... warning: diagonal scaling is not performed if ajj equals zero
!
      if (ajj.ne.0.0d0) b(j) = b(j)/ajj
  200 continue
!
!.... back substitution
!
      if (neq.eq.1) return
      jjnext = idiag(neq)
!
      do 400 j=neq,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
         bj = b(j)
         istart = j - jcolht + 1
         jtemp  = jjnext - istart + 1
!
         do 300 i=istart,j-1
         b(i) = b(i) - a(jtemp+i)*bj
  300    continue
!
      endif
!
  400 continue
!
      return
      end subroutine
!**** new *********************************************************
      subroutine factns(a,c,idiag,neq)
!
!.... program to perform crout factorization: a = l * d * u
!
!        a(i):  coefficient matrix stored in compacted column form;
!               after factorization contains d and u
!
!        c(i):  non-symmetric lower triangular coefficient matrix stored in
!                compacted row form; after factorization contains l
! 
!
      implicit none
!
!.... deactivate above card(s) for single-precision operation
!
      real*8  :: a(*),c(*)
      integer :: idiag(*)
      integer :: neq
!
      integer :: i, j, ii, jj, ij, iilast, jjlast
      integer :: istart,  icolht, jcolht, jtemp, jm1
      integer :: length
!
      jj = 0
!
      do 300 j=1,neq
!
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
!
      if (jcolht.gt.2) then
!
!....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j)
!
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
!
         do 100 i=istart,jm1
!
         iilast = ii
         ii     = idiag(i)
         icolht = ii - iilast
         length = min0(icolht-1,i - istart + 1)
         if (length.gt.0)  then
            a(ij) = a(ij) - coldot(a(ij-length),c(ii-length),length)
            c(ij) = c(ij) - coldot(c(ij-length),a(ii-length),length)
         endif
         ij = ij + 1
  100    continue
!
      endif
!
      if (jcolht.ge.2) then
!
!....... for column j and i.le.j-1, replace a(i,j) with u(i,j);
!           replace a(j,j) with d(j,j).
!
         jtemp = j - jj
!
         do 200 ij=jjlast+1,jj-1
!
         ii = idiag(jtemp + ij)
!
!....... warning: the following calculations are skipped 
!                 if a(ii) equals zero
!
         if (a(ii).ne.0.0d0) then
             c(ij) = c(ij)/a(ii)
             a(jj) = a(jj) - c(ij)*a(ij)
             a(ij) = a(ij)/a(ii)
         endif
  200    continue
!
      endif
!
  300 continue
!
      return
      end subroutine

!**** NEW ****************************************************************** 
      subroutine back(a,b,idiag,neq)

! 
!.... program to perform forward reduction and back substitution 
! 
      implicit none
! 
!.... remove above card for single-precision operation 
! 
      real*8  :: a(*)
      integer :: neq
      integer :: idiag(neq)
      real*8  :: b(neq)
!
      integer :: i,j,jcolht,istart,jtemp,jj,jjnext,jjlast
      real*8  :: ajj, bj
      integer :: iniA, iniB, fimA, fimB

!
!.... forward reduction 
! 
      jj = 0
! 

      do j=1,neq
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1) then
           iniA=jjlast+1
           iniB=j-jcolht+1
           fimA=iniA+jcolht-1-1
           fimB=iniB+jcolht-1-1
            b(j) = b(j) - coldot(a(jjlast+1),b(j-jcolht+1),jcolht-1)
      end if
      enddo
!
!.... diagonal scaling 
! 
      do j=1,neq
      ajj = a(idiag(j))
      if (ajj.ne.0.0d0) then
             b(j) = b(j)/ajj
      end if
      end do

! 
!.... back substitution 
! 
      if (neq.eq.1) return
      jjnext = idiag(neq)
! 
      do j=neq,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
        bj = b(j)
        istart = j - jcolht + 1
        jtemp  = jjnext - istart + 1
        do i=istart,j-1
           b(i) = b(i) - a(jtemp+i)*bj
        enddo
      endif
! 
      end do
! 
      return
      end subroutine

   
!**** new *************************************************************** 
      subroutine factor(a,idiag,nalhs,neq)
! 
!.... program to perform crout factorization: a = u(transpose) * d * u 
! 
!        a(i):  coefficient matrix stored in compacted column form; 
!               after factorization contains d and u 
! 
      implicit none
! 
!.... remove above card for single-precision operation 
! 
      integer, intent(in) :: neq, nalhs
      real*8, intent(inout)  :: a(nalhs)

      integer, intent(in) :: idiag(neq)
!
      integer :: i, j, jlast, icolht, jcolht, istart, jm1, jtemp
      integer :: ii, ij, jj, jlngth, length, iilast, jjnext
      integer :: jjlast
      real*8  :: ajj, bj, temp 
      integer :: iniA, iniB, fimA, fimB
! 
      jj = 0
      i=0; j=0; jlast=0; icolht=0; jcolht=0; istart=0; jm1=0; jtemp=0;
      ii=0;ij=0;jj=0;jlngth=0;length=0;iilast=0;jjnext=0;jjlast=0;
      ajj=0.0d0;bj=0.0d0;temp=0.0d0;

      do 300 j=1,neq
!
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
! 
      if (jcolht.gt.2) then
! 
!....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j) 
! 
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
! 
         do 100 i=istart,jm1
! 
         iilast = ii
         ii     = idiag(i)
         icolht = ii - iilast
         jlngth = i - istart + 1
         length = min0(icolht-1,jlngth)
         if (length.gt.0) then
           iniA=ii-length
           iniB=ij-length
           fimA=iniA+length-1
           fimB=iniB+length-1
           a(ij) = a(ij) - coldot(a(ii-length),a(ij-length),length)

         end if
         ij = ij + 1
  100    continue
!
      endif
! 
      if (jcolht.ge.2) then
! 
!....... for column j and i.le.j-1, replace a(i,j) with u(i,j); 
!           replace a(j,j) with d(j,j). 
! 
         jtemp = j - jj
! 
         do 200 ij=jjlast+1,jj-1
! 
         ii = idiag(jtemp + ij)
         if (a(ii).ne.0.0d0) then
            temp  = a(ij)
            a(ij) = temp/a(ii)
            a(jj) = a(jj) - temp*a(ij)
         endif
  200    continue
! 
      endif
! 
  300 continue

      return
      end subroutine


!**** new **********************************************************************
      subroutine btod(id,d,brhs,ndof,numnp)
!
!.... program to perform transfer from r.h.s. to displacement array
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: id(ndof,*)
      real*8  :: d(ndof,*),brhs(*)
      integer :: ndof, numnp
!
      integer :: i, j, k
!
      do 200 i=1,ndof
!
         do 100 j=1,numnp
         k = id(i,j)
         if (k.gt.0) then 
             d(i,j) = brhs(k)
         end if
  100    continue
!
  200    continue
 !
      return
      end subroutine


!**** new **********************************************************************
      subroutine kdbc(eleffm,elresf,dl,nee)
!
!.... program to adjust load vector for prescribed displacement
!     boundary condition
      use mGlobaisEscalares, only: zero
!
       implicit real*8 (a-h,o-z) 
!
!.... remove above card for single-precision operation
!
      integer :: nee
      real*8  :: eleffm(nee,*),elresf(*),dl(*)
!
      integer :: i,j
      real*8  :: val
!
      do 200 j=1,nee
!
      val=dl(j)
      if(val.eq.zero) go to 200
!
      do 100 i=1,nee
      elresf(i)=elresf(i)-eleffm(i,j)*val
100   continue
!
200   continue
!
      return
      end subroutine

!**** new **********************************************************************
      subroutine pivots(a,idiag,neq,nsq,iecho,*)
!
!.... program to determine the number of zero and negative terms in
!        array d of factorization a = u(transpose) * d * u
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                
      real*8  :: a(*)                                                       
      integer :: idiag(*)
      integer :: neq, nsq, iecho
!
      integer :: iz, in, n, i
!
      iz = 0
      in = 0
!
      do 100 n=1,neq
      i = idiag(n)
      if (a(i).eq.0.) iz = iz + 1
      if (a(i).lt.0.) in = in + 1
  100 continue
!
      write(iecho,1000) nsq,iz,in
!
      return 1
!
 1000 format(' ',&
     ' zero and/or negative pivots encountered                ', ///5x,&
     ' time sequence number   . . . . . . . . . . . (nsq  ) = ',i10//5x,&
     ' number of zeroes . . . . . . . . . . . . . . . . . . = ',i10//5x,&
     ' number of negatives  . . . . . . . . . . . . . . . . = ',i10//5x)
!
      end subroutine

!**** new **********************************************************************
      subroutine ftod(id,d,f,ndof,numnp,nlvect)
!
!.... program to compute displacement boundary conditions
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: id(ndof,*)
      real*8  :: d(ndof,*),f(ndof,numnp,*)
      integer :: ndof, numnp, nlvect
!
      integer :: i, j, k, lv
      real*8  :: val
!
!
      do 300 i=1,ndof
!
            do 200 j=1,numnp
!
            k = id(i,j)
            if (k.gt.0) go to 200
            val = 0.0d0
                  do 100 lv=1,nlvect
                  val = val + f(i,j,lv)
100               continue
!
            d(i,j) = val
!
  200       continue
!
 300  continue
      return
      end subroutine

!**** new **********************************************************************
      subroutine btdb(elstif,b,db,nee,nrowb,nstr)
!
!.... program to multiply b(transpose) * db taking account of symmetry
!        and accumulate into element stiffness matrix
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: elstif(nee,*),b(nrowb,*),db(nrowb,*)
      integer :: nee,nrowb,nstr
!
      integer :: i,j
!
      do 200 j=1,nee
!
      do 100 i=1,j
      elstif(i,j) = elstif(i,j) + coldot(b(1,i),db(1,j),nstr)   
  100 continue
!
  200 continue
!
      return
      end subroutine

!**** new **********************************************************************
      subroutine colht(idiag,lm,ned,nen,numel,neq)
!
!.... program to compute column heights in global left-hand-side matrix
!
      implicit none
      integer :: idiag(*),lm(ned,nen,*)
      integer :: ned, nen, numel, neq
!
      integer :: i, j, k
      integer :: m, min, num
!
      do 500 k=1,numel
      min = neq
!
      do 200 j=1,nen
!
      do 100 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) min = min0(min,num)
  100 continue
!
  200 continue
!
      do 400 j=1,nen
!
      do 300 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) then
         m = num - min
         if (m.gt.idiag(num)) then
              idiag(num) = m
         endif
      endif
!
  300 continue

  400 continue
!
  500 continue
!
      return
      end subroutine
! **** new *********************************************************************
       subroutine addlhs(alhs,eleffm,idiag,lm,nee,ldiag,lsym)
! 
! .... program to add element left-hand-side matrix to
!         global left-hand-side matrix
! 
!         ldiag = .true.,  add diagonal element matrix
! 
!         ldiag = .false., then
!           lsym = .true., add upper triangle of full element matrix
!           lsym = .false., add full element matrix
! 
       implicit none
! 
! .... deactivate above card(s) for single-precision operation
! 
       real*8  :: alhs(*),eleffm(nee,*)
       integer :: idiag(*),lm(*)
       integer :: nee
       logical :: ldiag,lsym
!
       integer :: i, j, l, k, m

       if (ldiag) then
! 
          do 100 j=1,nee
          k = lm(j)
          if (k.gt.0) then
             l = idiag(k)
             alhs(l) = alhs(l) + eleffm(j,j)
          endif
  100     continue
! 
       else
! 
          do 400 j=1,nee
          k = lm(j)
          if (k.gt.0) then
! 
             do 200 i=1,j
             m = lm(i)
             if (m.gt.0) then
                if (k.ge.m) then
                   l = idiag(k) - k + m
                else
                   l = idiag(m) - m + k
                endif
                alhs(l) = alhs(l) + eleffm(i,j)
             endif
  200       continue
! 
             if (.not. lsym) then
                do 300 i = j,nee
                m = lm(i)
               if (m .gt. 0) then
                 if (k .ge. m) then
                   l = idiag(k) - k + m
                else
                   l = idiag(m) - m + k
                 endif
               endif
  300         continue
             endif
          endif
  400    continue
! 
       endif
! 
       return
       end subroutine

!      
!:::: new ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
       subroutine addlhsGEO(alhs,clhs,eleffm,idiag,lm,nee,ldiag,lsym) 
! 
! .... program to add element left-hand-side matrix to 
!         global left-hand-side matrix 
! 
!         ldiag = .true.,  add diagonal element matrix 
! 
!         ldiag = .false., then 
!           lsym = .true., add upper triangle of full element matrix 
!           lsym = .false., add full element matrix 
! 
       implicit real*8 (a-h,o-z) 
! 
! .... deactivate above card(s) for single-precision operation 
! 
       integer :: nee
       real*8  :: alhs(*),clhs(*),eleffm(nee,*)
       integer :: idiag(*),lm(*) 
       logical :: ldiag,lsym 
!
       integer :: i,j,k,l,m
! 
       if (ldiag) then 
! 
          do 100 j=1,nee 
          k = lm(j) 
          if (k.gt.0) then 
             l = idiag(k) 
             alhs(l) = alhs(l) + eleffm(j,j) 
          endif 
  100     continue 
! 
       else 
! 
          do 400 j=1,nee 
          k = lm(j) 
          if (k.gt.0) then 
! 
             do 200 i=1,j 
             m = lm(i) 
             if (m.gt.0) then 
                if (k.ge.m) then 
                   l = idiag(k) - k + m 
                else 
                   l = idiag(m) - m + k 
                endif 
                alhs(l) = alhs(l) + eleffm(i,j) 
             endif 
  200       continue 
! 
             if (.not. lsym) then 
                do 300 i = j,nee 
                m = lm(i) 
               if (m .gt. 0) then 
                 if (k .ge. m) then 
                   l = idiag(k) - k + m 
                else 
                   l = idiag(m) - m + k 
                 endif 
                 clhs(l) = clhs(l) + eleffm(i,j) 
               endif 
  300         continue 
             endif 
          endif 
  400    continue 
! 
       endif 
! 
       return 
       end subroutine
!:
!**** new **********************************************************************
      subroutine addrhs (brhs,elresf,lm,nee)
!
!.... program to add element residual-force vector to
!        global right-hand-side vector
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: brhs(*),elresf(*)
      integer :: lm(*)
      integer :: nee
!
      integer :: k, j
!
      do 100 j=1,nee
      k = lm(j)
      if (k.gt.0) brhs(k) = brhs(k) + elresf(j)
  100 continue
!
      return
      end subroutine

!**** new **********************************************************************
      function coldot(a,b,n)
!
!.... program to compute the dot product of vectors stored column-wise
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: a(*),b(*)
      integer :: n
!
      real*8  :: coldot
      integer :: i
!

      coldot = 0.0d0
!
       coldot=dot_product(a(1:n),b(1:n))

!   100 continue

!
      return
      end function
!**** new **********************************************************************
      subroutine matadd(a,b,c,ma,mb,mc,m,n,iopt)
!
!.... program to add rectangular matrices
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: a(ma,*),b(mb,*),c(mc,*)
      integer :: ma,mb,mc,m,n,iopt
!
      integer :: i,j
!
      go to (1000,2000,3000),iopt
!
!.... iopt = 1, add entire matrices
!
 1000 do 1200 j=1,n
!
      do 1100 i=1,m 
      c(i,j) = a(i,j) + b(i,j)
 1100 continue
!
 1200 continue
      return
!
!.... iopt = 2, add lower triangular and diagonal elements
!
 2000 do 2200 j=1,n
!
      do 2100 i=j,m 
      c(i,j) = a(i,j) + b(i,j)
 2100 continue
!
 2200 continue
      return
!
!.... iopt = 3, add upper triangular and diagonal elements
!
 3000 do 3200 j=1,n
!
      do 3100 i=1,j 
      c(i,j) = a(i,j) + b(i,j)
 3100 continue
!
 3200 continue
      return
!
      end subroutine

!**** new **********************************************************************
      function rowdot(a,b,ma,mb,n)
!
!.... program to compute the dot product of vectors stored row-wise
!
      implicit none
!
!.... remove above card for single precision operation
!
      real*8  :: a(ma,*),b(mb,*)
      integer :: ma, mb, n
!
      real*8  :: rowdot
      integer :: i
!
      rowdot = 0.0d00
!
      do i=1,n
      rowdot = rowdot + a(1,i)*b(1,i)
      enddo
!
      return
      end function

  end module

