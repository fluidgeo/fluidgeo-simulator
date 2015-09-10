!=================================================================================
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

      module mMalha
      
         real*8, allocatable :: x_BM(:,:), xc_BM(:,:)
         real*8, allocatable :: x_B(:,:), xc_B(:,:)

         integer, allocatable :: conecNodaisElem_BM(:,:), conecLadaisElem_BM(:,:)
         integer, allocatable :: conecNodaisElem_B(:,:), conecLadaisElem_B(:,:)

         integer, allocatable :: listaDosElemsPorNo_BM(:,:), listaDosElemsPorFace_BM(:,:)
         integer, allocatable :: listaDosElemsPorNo_B(:,:), listaDosElemsPorFace_B(:,:)

         integer :: nsd_BM, numel_BM, numnp_BM, numLados_BM, nen_BM, numLadosELem_BM
         integer :: nsd_B, numel_B, numnp_B, numLados_B, nen_B, numLadosELem_B

         integer :: numnpD

!          integer :: numnpReserv, numelReserv
!          integer :: numnpReserv, numelReserv
         
         integer :: nelx_BM, nely_BM, nelz_BM
         integer :: nelx_B, nely_B, nelz_B
!
         public :: criarListaVizinhos, nbound, formlm
         public :: genel, genel1, geneld, geneli, genfl, genfl1
         public :: gensh, gensh1, gensh2, gensh3, igen

!
       contains
!
!=======================================================================
!     
      subroutine criarListaVizinhos(nen,numnp,numel,conecElem,listaDosElemsPorNo)
      implicit none
!     
!     Objetivo: cria a matriz listaDosElemsPorNo: .... indices dos elementos que possuem o no
!
      integer nen,numnp,numel
      integer, dimension(nen,numel)  :: conecElem
      integer, dimension(nen,numnp)  :: listaDosElemsPorNo
      integer :: no,nel,l

      do nel=1, numel
         do l=1, nen
            no=conecElem(l,nel)
            listaDosElemsPorNo(l,no) = nel
         end do
      end do
    
      end subroutine
!
!=======================================================================
!     
      subroutine nconec0(nen,numnp,numel,conecElem,listaDosElemsPorNo)
      implicit none
!     
!     Objetivo: cria a matriz listaDosElemsPorNo: .... indices dos elementos que possuem o no
!
      integer nen,numnp,numel
      integer, dimension(nen,numel)  :: conecElem
      integer, dimension(nen,numnp)  :: listaDosElemsPorNo
      integer i,j,no,nel,ni,ns

      character(len=22) :: nomeArq
      integer :: ilistaDosElemsPorNo
      logical :: existe


      if(nen==6) then
         nomeArq="listaDosElemsPorFc.dat"
         ilistaDosElemsPorNo=445
      else
         nomeArq="listaDosElemsPorNo.dat"
         ilistaDosElemsPorNo=444
      end if

      inquire (file=nomeArq,exist=existe)
      open(UNIT=ilistaDosElemsPorNo,FILE=nomeArq)
!      open(UNIT=ilistaDosElemsPorNo,FILE=nomeArq,FORM='UNFORMATTED')

      if(existe.eqv..false.) then
!     
!      nmax=max(nelx,nely)
!      nmin=min(nelx,nely)
!      nmax=nelx*nely
!
      listaDosElemsPorNo=0

!
!$OMP PARALLEL DO PRIVATE(J, NO, NEL)
      do i=1,numnp
!     
!     centra no noh
!
!         ni=i*nmax
!         if(ni.lt.0) ni=1
!     
         ni=1
!
!         ns=i+nmax
!         if(ns.gt.numel) ns=numel
!
         ns=numel
!
         do nel=ni,ns
!     
            do j=1,nen
!     
               no=conecElem(j,nel)
               if(no.eq.i) listaDosElemsPorNo(j,i) = nel
!     
            end do
!     
         end do
!     
      end do
!$OMP END PARALLEL DO

       print*, "escrevendo arquivo: ", nomeArq
       do i=1, numnp
       write(ilistaDosElemsPorNo,*), listaDosElemsPorNo(1:nen,i)
       end do
!      write(ilistaDosElemsPorNo,*), listaDosElemsPorNo ! teste com arq.formatado
!      write(ilistaDosElemsPorNo), listaDosElemsPorNo ! teste com arq.nao-formatado


      else
          print*, "lendo arquivo: ", nomeArq
          do i=1, numnp
            read(ilistaDosElemsPorNo,*) (listaDosElemsPorNo(j,i),j=1,nen)
          end do
          
!       read(ilistaDosElemsPorNo) listaDosElemsPorNo
      end if

      close(ilistaDosElemsPorNo)
    
      end subroutine
!
!=======================================================================
!     
      subroutine nbound(nen,numel,nit,ni,conecElem,listaDosElemsPorNo)
!     
!     identifica os volumes (e nos) interiores e do contorno
!     
      implicit none
!     
      integer nen,numel
      integer, dimension(*)       :: ni
      integer, dimension(nen,*)   :: conecElem
      integer, dimension(nen,*)   :: listaDosElemsPorNo
      integer, dimension(13)      :: nit
!     
      integer, dimension(:,:), allocatable :: inaux
!
      integer i,no,nel,nflag,nicont,nbcont,naux
      integer n,n1,n2,n3,n4,ncont
      integer nIn,nU,nD,nL,nR,nLU,nLD,nRu,nRD
!     
!     contadores: nicont = interior
!                 nbcont = boundary
!     
!     identificacao do treze tipos de elementos de fronteira
!     I  =1;  L  =2;  D  =3;  R  =4;    U=5;
!     LD =6;  RD =7;  RU =8;  LU =9;
!     LDI=10; RDI=11; RUI=12; LUI=13
!
      allocate(inaux(2,numel))
!
      nIn=0
      nU=0
      nD=0
      nL=0
      nR=0
      nLU=0
      nLD=0
      nRu=0
      nRD=0
!
      do i=1,13
         nit(i)=0
      end do
!
      nicont=0
      nbcont=numel+1
!
      do nel=1,numel
         nicont=nicont+1
         nflag=1
         do i=1,nen
            n=0
            no=conecElem(i,nel)
            n1=0
            n2=0
            n3=0
            n4=0
            ncont=0
            do n1=1,nen
               naux=listaDosElemsPorNo(n1,no)
!                print*, "n1",  n1, no, "naux=", naux
               if(naux.eq.0)then
                  n=n+1
                  nflag=0
                  if(n.eq.3)then
                     go to 113
                  else
                     if(n1.eq.4) go to 113
                     do n2=n1+1,nen
                        naux=listaDosElemsPorNo(n2,no)
!                         print*, "n2",  n2, no, "naux=", naux
                        if(naux.eq.0)then
                           n=n+1
                           if(n.eq.3)then
                              go to 113
                           else
                              if(n2.eq.4) go to 113
                              do n4=n2+1,nen
                                 naux=listaDosElemsPorNo(n4,no)
!                                  print*, "n4",  n4, no, "naux=",naux
                                 if(naux.eq.0)then
                                    n=n+1
                                    n3=n4
                                    go to 113
                                 else
                                    n3=0
                                    go to 113
                                 end if
                              end do 
                           end if
                        end if
                     end do
                  end if
               end if
            end do
         end do
!     
!
 113     continue
!
         if(nflag.eq.1) then
            inaux(1,nicont)=nel
            inaux(2,nicont)=1
            nit(1)=nit(1)+1
            nIn=nIn+1
            ni(nIn)=nel
         else
!            write(*,*)'%%%',nel,n1,n2,n3,n,i
            if(n1.eq.1)then
               if(n2.eq.2)then
                  nit(5)=nit(5)+1
                  inaux(1,nicont)=nel
                  inaux(2,nicont)=5
               else
                  no=conecElem(3,nel)
                  if(listaDosElemsPorNo(2,no).eq.0)then
!                     write(*,*)'RU=',nel
                     nit(8)=nit(8)+1
                     inaux(1,nicont)=nel
                     inaux(2,nicont)=8                  
                  else
                     nit(4)=nit(4)+1
                     inaux(1,nicont)=nel
                     inaux(2,nicont)=4
                  end if
               end if
            else
               if(n1.eq.3)then
                  no=conecElem(2,nel)
                  if(listaDosElemsPorNo(1,no).eq.0)then
!                     write(*,*)'RD=',nel
                     nit(7)=nit(7)+1
                     inaux(1,nicont)=nel
                     inaux(2,nicont)=7
                  else
                     nit(3)=nit(3)+1
                     inaux(1,nicont)=nel
                     inaux(2,nicont)=3
                  end if
               else
                  if(n3.eq.4)then
!                    write(*,*)'LD=',nel
                     nit(6)=nit(6)+1
                     inaux(1,nicont)=nel
                     inaux(2,nicont)=6
                  else
                     no=conecElem(4,nel)
                     if(listaDosElemsPorNo(1,no).eq.0)then
!                        write(*,*)'LU=',nel
                        nit(9)=nit(9)+1
                        inaux(1,nicont)=nel
                        inaux(2,nicont)=9
                     else
                        nit(2)=nit(2)+1
                        inaux(1,nicont)=nel
                        inaux(2,nicont)=2
                     end if
               end if           
            end if
         end if
      end if
!
      end do                    !nel->numel
!
!
! ordenacao dos valores em ni
!
!     ajuste dos nites > marcam a posicao final
      nit(2)=nit(1)+nit(2)
      nit(3)=nit(2)+nit(3)
      nit(4)=nit(3)+nit(4)
      nit(5)=nit(4)+nit(5)
      nit(6)=nit(5)+nit(6)
      nit(7)=nit(6)+nit(7)
      nit(8)=nit(7)+nit(8)
      nit(9)=nit(8)+nit(9)
!
!      write(*,*)nit(1),nit(2),nit(3),nit(4),nit(5),
!     & nit(6),nit(7),nit(8),nit(9)
      do nel=1,numel
!
!     interiores foram feitos anteriormente
!
!     L ###################################        
         if(inaux(2,nel).eq.2)then
            nL=nL+1
            ni(nit(1)+nL)=nel
            go to 123
         end if
!
!     D ###################################        
         if(inaux(2,nel).eq.3)then
            nD=nD+1
            ni(nit(2)+nD)=nel
            go to 123
         end if
!     R ###################################        
         if(inaux(2,nel).eq.4)then
            nR=nR+1
            ni(nit(3)+nR)=nel
            go to 123
         end if
!     U ###################################        
         if(inaux(2,nel).eq.5)then
            nU=nU+1
            ni(nit(4)+nU)=nel
            go to 123
         end if
!     LD###################################        
         if(inaux(2,nel).eq.6)then
            ni(nit(5)+1)=nel
            go to 123
         end if
!     RD###################################        
         if(inaux(2,nel).eq.7)then
            ni(nit(6)+1)=nel
            go to 123
         end if
!     RU###################################        
         if(inaux(2,nel).eq.8)then
            ni(nit(7)+1)=nel
            go to 123
         end if
!     LU###################################        
         if(inaux(2,nel).eq.9)then
            ni(nit(8)+1)=nel
            go to 123
         end if
!       ###################################
 123     continue
      end do                    !numel
!
      deallocate(inaux)
!
      end subroutine
!      
!**** new **********************************************************************
      subroutine formlm (id,conecElem,lm,ndof,ned,nen,numel)
!
!.... program to form lm array
!
      integer :: ndof,ned,nen,numel
      integer :: id(ndof,*),conecElem(nen,*),lm(ned,nen,*)
!
      integer :: i,j,k,node      
!
      do 300 k=1,numel
!
      do 200 j=1,nen
      node=conecElem(j,k)
!
      do 100 i=1,ndof

      lm(i,j,k) = id(i,node)

  100 continue
!
  200 continue
!
  300 continue

!
      return
      end subroutine

!******************************************************************************
      subroutine genel(conecElem,mat,nen, iin)   
!                                                                       
!.... program to read and generate element node and material numbers    
!                                                                       
!         conecElem(nen,numel) = element node numbers                         
!         mat(numel)     = element material numbers                     
!         nen            = number of element nodes (le.27)              
!         n              = element number                               
!         ng             = generation parameter                         
!         nel(i)         = number of elements in direction i            
!         incel(i)       = element number increment for direction i     
!         inc(i)         = node number increment for direction i        
!   
      integer :: nen, iin                                                                    
      integer :: conecElem(nen,*),mat(*),itemp(27) 
!
      integer :: n,m,ng,i
      common /genelc/ n,nel(3),incel(3),inc(3)   
      
      character(len=21) :: testeStr
!                                                                       
  100 continue                                                          
      read(iin,1000) n,m,(itemp(i),i=1,nen),ng     
                
      if (n.eq.0) return                                                
!       call imove(conecElem(1,n),itemp,nen)                                    
      conecElem(1:nen,n)=itemp(1:nen)                                    
      mat(n)=m                                                          
      if (ng.ne.0) then
!                                                                       
!....... generate data                                                     
!                                                                       
         read(iin,1000) (nel(i),incel(i),inc(i),i=1,3)     
!           write(*,1000) (nel(i),incel(i),inc(i),i=1,3)                  
         call genel1(conecElem,mat,nen)                                          
      endif
      go to 100       
      
      read(iin, 1100) testeStr
!                                                                       
 1000 format(16i10,10x,14i10)                                             
 1100 format(20a4)

!                                                                       
      end subroutine          

!:::: new ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine imove(ia,ib,n) 
! 
!.... program to move an integer array 
! 
      integer :: ia(1),ib(*) 
      integer :: n
!
      integer :: i
! 
      do 100 i=1,n 
      ia(i)=ib(i) 
  100 continue 
 
      return 
      end subroutine
!:::: new :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
      subroutine move(a,b,n)
! 
!.... program to move a floating-point array 
! 
      implicit real*8(a-h,o-z) 
! 
!.... remove above card for single-precision operation 
! 
      real*8 :: a(*),b(*) 
      integer :: n
!
      integer i
! 
      do 100 i=1,n 
      a(i) = b(i) 
  100 continue 
! 
      return 
      end subroutine
                                                     
!******************************************************************************
      subroutine genel1(conecElem,mat,nen)                                    
!                                                                       
!.... program to generate element node and material numbers             
!                                        
      integer :: nen                               
      integer:: conecElem(nen,*),mat(*)                                       
!
      integer :: i,j,k,ii,jj,kk,ie,je,ke,le
      common /genelc/ n,nel(3),incel(3),inc(3)                          
!                                                                       
!.... set defaults                                                      
!                                                                       
      call geneld
!                                                                       
!.... generation algorithm                                              
!                                                                       
      ie = n                                                            
      je = n                                                            
      ke = n                                                            
!                                                                      
      ii = nel(1)                                                       
      jj = nel(2)                                                       
      kk = nel(3)                                                       
!                                                                       
      do 300 k=1,kk                                                     
!
      do 200 j=1,jj                                                     
!
      do 100 i=1,ii                                                     
!                                                                       
      if (i.ne.ii) then
         le = ie                                                           
         ie = le + incel(1)                                                
         call geneli(conecElem(1,ie),conecElem(1,le),inc(1),nen)                       
         mat(ie) = mat(le)                                                 
      endif
  100 continue                                                          
!                                                                       
      if (j.ne.jj) then
         le = je                                                           
         je = le + incel(2)                                                
         call geneli(conecElem(1,je),conecElem(1,le),inc(2),nen)                       
         mat(je) = mat(le)                                                 
         ie = je                                                           
      endif
  200 continue                                                          
!                                                                       
      if (k.ne.kk) then
         le = ke                                                           
         ke = le + incel(3)                                                
         call geneli(conecElem(1,ke),conecElem(1,le),inc(3),nen)                       
         mat(ke) = mat(le)                                                 
         ie = ke
         je=ke                                                           
      endif
  300 continue                                                          
!                                                                      
      return                                                            
      end subroutine                                                               
!******************************************************************************
      subroutine geneld                                                 
!                                                                       
!.... program to set defaults for element node       
!        and material number generation                              
!                                                                       
      common /genelc/ n,nel(3),incel(3),inc(3)
!                                                                       
      if (nel(1).eq.0) nel(1) = 1                                       
      if (nel(2).eq.0) nel(2) = 1                                       
      if (nel(3).eq.0) nel(3) = 1                                       
!                                                                       
      if (incel(1).eq.0) incel(1) = 1                                   
      if (incel(2).eq.0) incel(2) = nel(1)                              
      if (incel(3).eq.0) incel(3) = nel(1)*nel(2)                       
!                                                                       
      if (inc(1).eq.0) inc(1) = 1                                       
      if (inc(2).eq.0) inc(2) = (1+nel(1))*inc(1)                       
      if (inc(3).eq.0) inc(3) = (1+nel(2))*inc(2)                       
!                                                                       
      return                                                            
      end subroutine                                                               
!******************************************************************************
      subroutine geneli(conecElem2,conecElem1,inc,nen)                              
!                                                                       
!.... program to increment element node numbers                         
!    
      integer :: inc, nen                                                                   
      integer:: conecElem1(*),conecElem2(*)                                         
!
      integer :: i
!                                                                       
      do 100 i=1,nen                                                    
      if (conecElem1(i).eq.0) then
         conecElem2(i) = 0
      else
         conecElem2(i) = conecElem1(i) + inc                         
      endif
  100 continue                                                          
!                                                                       
      return                                                            
      end  subroutine                                                              
!******************************************************************************
      subroutine genfl(a,nra,iin)   
      use mGlobaisEscalares  
!                                                                       
!.... program to read and generate floating-point nodal data            
!                                                                       
!         a       = input array                                         
!         nra     = number of rows in a (le.6)                          
!         n       = node number                                         
!         numgp   = number of generation points                         
!         ninc(i) = number of increments for direction i                
!         inc(i)  = increment for direction i                           
!                                                                       
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                        
      integer :: nra, iin                                               
      real*8  :: a(nra,*)
!
      real*8  :: temp(6,20)
      integer :: n,numgp,ninc(3),inc(3)   
      integer :: i, j, m, mgen
!                                                                       
  100 continue                                                          
      read(iin,1000) n,numgp,(temp(i,1),i=1,nra)      
!       write(*,1000) n,numgp,(temp(i,1),i=1,nra)      
 
      if (n.eq.0) return                                                
      call move(a(1,n),temp,nra)           
 !     a(1:nra,n) = temp(1:nra,1)

      if (numgp.ne.0) then
         do 200 j=2,numgp
                                                                       
         read(iin,1000) m,mgen,(temp(i,j),i=1,nra)
!           write(*,1000) m,mgen,(temp(i,j),i=1,nra)        

        if (mgen.ne.0) call move(temp(1,j),a(1,m),nra) 
!          if (mgen.ne.0) temp(1:nra,j)=a(1:nra,m) 
!                                                                       
  200    continue  

         read(iin,2000) (ninc(i),inc(i),i=1,3)
!          write(*,2000) (ninc(i),inc(i),i=1,3)
 
         call genfl1(a,nra, temp, n, numgp, ninc, inc)                                                
      endif
      go to 100                                                         
!                                                                       
! 1000 format(2i10,6f10.0)                                                
 1000 format(2i10,6f10.0)                                                
 2000 format(16i10)                                                      
!                                                                       
      end  subroutine                                                              
!******************************************************************************
      subroutine genfl1(a, nra, temp, n, numgp, ninc, inc)    
      use mGlobaisEscalares 
!                                                                       
!.... program to generate floating-point nodal data 
!        via isoparametric interpolation         
!                                                                       
!         iopt = 1, generation along a line                             
!              = 2, generation over a surface                           
!              = 3, generation within a volume                            
!                                                                       
      implicit none                                          
!                                                                       
!.... remove above card for single-precision operation                  
!
      integer :: nra
      real*8  :: a(nra,*)
      real*8  :: temp(6,20)
      integer :: n,numgp,ninc(3),inc(3)                 
!
      real*8  :: sh(20,1), dr, ds, dt, r, s, t
      integer :: iopt, ni, nj, nk, ii, jj, kk, i, j, k
      real*8  :: resultado(6,1)
!
      iopt = 3                                                          
      if (ninc(3).eq.0) iopt = 2                                        
      if (ninc(2).eq.0) iopt = 1                                        
!                                                                       
      dr = 0.0
      ds = 0.0
      dt = 0.0
!                                                                       
      if (ninc(1).ne.0) dr = two/ninc(1)                                
      if (ninc(2).ne.0) ds = two/ninc(2)                                
      if (ninc(3).ne.0) dt = two/ninc(3)                                
!                                                                       
      ii = ninc(1)+1                                                    
      jj = ninc(2)+1                                                    
      kk = ninc(3)+1                                                    
!                                                                       
      ni = n                                                            
      nj = n                                                            
      nk = n                                                            
!                                                                       
      t = -one                                                          
      do 300 k=1,kk                                                     
!
      s = -one                                                          
      do 200 j=1,jj                                                     
!
      r = -one                                                          
      do 100 i=1,ii                                                     
!                                                                       
      call gensh(r,s,t,sh,numgp,iopt)                                   
 
!       resultado=matmul(temp,sh)
!       a(1:nra,ni)=resultado(1:nra,1)
      call multab(temp,sh,a(1,ni),6,20,nra,numgp,nra,1,1)       

      ni = ni + inc(1)                                                      
      r = r + dr                                                            
  100 continue                                                          
!                                                                       
      nj = nj + inc(2)                                                      
      ni = nj                                                             
      s = s + ds                                                            
  200 continue                                                          
!                                                                       
      nk = nk + inc(3)                                                      
      ni = nk                                                             
      t = t + dt                                                            
  300 continue                                                          
!                                                                       
      return                                                            
      end  subroutine          
!**** new **********************************************************************
      subroutine multab(a,b,c,ma,mb,mc,l,m,n,iopt)
      use mAlgMatricial, only: coldot, rowdot
!
!.... program to multiply two matrices
!
!        l = range of dot-product index
!        m = number of active rows in c
!        n = number of active columns in c
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                                                                       
      real*8  :: a(ma,*),b(mb,*),c(mc,*)
      integer :: ma,mb,mc,l,m,n,iopt
!
      integer :: i,j
!
      go to (1000,2000,3000,4000),iopt
!
!.... iopt = 1, c(i,j) = a(i,k)*b(k,j) , (c = a * b)
!
 1000 do 1200 i=1,m
!
      do 1100 j=1,n
      c(i,j) = rcdot(a(i,1),b(1,j),ma,l)
 1100 continue
!
 1200 continue
      return
!                                            t
!.... iopt = 2, c(i,j) = a(k,i)*b(k,j) (c = a  * b)
!
 2000 do 2200 i=1,m
!
      do 2100 j=1,n
      c(i,j) = coldot(a(1,i),b(1,j),l)
 2100 continue
!
 2200 continue
      return
!                                                t
!.... iopt = 3, c(i,j) = a(i,k)*b(j,k) (c = a * b )
!
 3000 do 3200 i=1,m
!
      do 3100 j=1,n
      c(i,j) = rowdot(a(i,1),b(j,1),ma,mb,l)
 3100 continue
!
 3200 continue
      return
!                                            t    t
!.... iopt = 4, c(i,j) = a(k,i)*b(j,k) (c = a  * b )
!
 4000 do 4200 i=1,m
!
      do 4100 j=1,n
      c(i,j) = rcdot(b(j,1),a(1,i),mb,l)
 4100 continue
!
 4200 continue

!
      return
      end subroutine

!**** new **********************************************************************
      function rcdot(a,b,ma,n)
!
!.... program to compute the dot product of a vector stored row-wise
!        with a vector stored column-wise
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                                                                       
      real*8  :: a(ma,*),b(*)
      integer :: ma, n
!
      real*8  :: rcdot
      integer :: i
!
      rcdot = 0.0
!
      do 100 i=1,n
      rcdot = rcdot + a(1,i)*b(i)
  100 continue
!
      return
      end function
                                                    
!****************************************************************************** 
    
      subroutine gensh(r,s,t,sh,numgp,iopt)                             
!                                                                       
!.... program to call shape function routines         
!        for isoparametric generation         
!                                                                       
      implicit none                                      
!                                                                       
!.... modify above card for single-precision operation               
!               
      real*8  :: r, s, t, sh(*)                                                        
      integer :: numgp, iopt
!                                                                       
      go to (100,200,300),iopt                                                
!                                                                       
  100 call gensh1(r,sh,numgp)                                           
      return                                                            
!                                                                       
  200 call gensh2(r,s,sh,numgp)                                         
      return                                                            
!                                                                       
  300 call gensh3(r,s,t,sh,numgp)                                       
      return                                                            
!                                                                       
      end subroutine                                                               
!******************************************************************************
      subroutine gensh1(r,sh,n)  
      use mGlobaisEscalares                                       
!                                                                       
!.... program to compute 1d shape functions           
!        for isoparametric generation                     
!                                                                       
      implicit none                                          
!                                                                       
!.... modify above card(s) for single-precision operation               
!                                                                       
      real*8  :: r, sh(*)                                                   
      integer :: n
!                                                                       
      sh(2) = pt5*r                                                       
      sh(1) = pt5 - sh(2)                                                   
      sh(2) = pt5 + sh(2)                                                   
      if (n.eq.3) then
         sh(3) = one - r*r                                                     
         sh(1) = sh(1) - pt5*sh(3)                                             
         sh(2) = sh(2) - pt5*sh(3)                                             
      endif
!                                                                       
      return                                                            
      end subroutine                                                               
!******************************************************************************
      subroutine gensh2(r,s,sh,n)     
      use mGlobaisEscalares                                  
!
!.... program to compute 2d shape functions 
!        for isoparametric generation    
!                                                                       
      implicit none                                         
!
!.... modify above card for single-precision operation               
!                                                                       
      real*8  :: r, s, sh(*)                                                   
      integer :: n    
!
      real*8  :: r1, r2, r3, s1, s2, s3
!
      r2 = pt5*r                                                          
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s                                                          
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      sh(1) = r1*s1                                                       
      sh(2) = r2*s1                                                       
      sh(3) = r2*s2                                                       
      sh(4) = r1*s2                                                       
      if (n.eq.4) return                                                
!                                                                       
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      sh(5) = r3*s1                                                       
      sh(6) = s3*r2                                                       
      sh(7) = r3*s2                                                       
      sh(8) = s3*r1                                                       
      sh(1) = sh(1) - pt5*(sh(5) + sh(8))
      sh(2) = sh(2) - pt5*(sh(6) + sh(5))
      sh(3) = sh(3) - pt5*(sh(7) + sh(6))
      sh(4) = sh(4) - pt5*(sh(8) + sh(7))
!                                                                       
      return
      end subroutine                                                               
!******************************************************************************
      subroutine gensh3(r,s,t,sh,n) 
      use mGlobaisEscalares                                    
!                                                                       
!.... program to compute 3d shape functions            
!        for isoparametric generation   
!                                                                       
      implicit none                                         
!                                                                       
!.... modify above card for single-precision operation               
!                                                                       
      real*8  :: r, s, t, sh(*)                                                   
      integer :: n    
!
      real*8  :: r1, r2, r3, rs1, rs2, rs3, rs4
      real*8  :: s1, s2, s3, t1, t2, t3                                              
!                                                                       
      r2 = pt5*r
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      t2 = pt5*t
      t1 = pt5 - t2                                                         
      t2 = pt5 + t2                                                         
!                                                                       
      rs1 = r1*s1                                                         
      rs2 = r2*s1                                                         
      rs3 = r2*s2                                                         
      rs4 = r1*s2                                                         
      sh(1) = rs1*t1                                                      
      sh(2) = rs2*t1                                                      
      sh(3) = rs3*t1                                                      
      sh(4) = rs4*t1                                                      
      sh(5) = rs1*t2                                                      
      sh(6) = rs2*t2                                                      
      sh(7) = rs3*t2                                                      
      sh(8) = rs4*t2                                                      
      if (n.eq.8) return                                                 
!                                                                       
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      t3 = one - t*t                                                        
      sh(17) = t3*rs1                                                     
      sh(18) = t3*rs2                                                     
      sh(19) = t3*rs3                                                     
      sh(20) = t3*rs4                                                     
      rs1 = r3*s1                                                         
      rs2 = s3*r2                                                         
      rs3 = r3*s2                                                         
      rs4 = s3*r1                                                         
      sh( 9) = rs1*t1                                                     
      sh(10) = rs2*t1                                                     
      sh(11) = rs3*t1                                                     
      sh(12) = rs4*t1                                                     
      sh(13) = rs1*t2                                                     
      sh(14) = rs2*t2                                                     
      sh(15) = rs3*t2                                                     
      sh(16) = rs4*t2                                                     
!                                                                       
      sh(1) = sh(1) - pt5*(sh( 9) + sh(12) + sh(17))
      sh(2) = sh(2) - pt5*(sh( 9) + sh(10) + sh(18))
      sh(3) = sh(3) - pt5*(sh(10) + sh(11) + sh(19))
      sh(4) = sh(4) - pt5*(sh(11) + sh(12) + sh(20))
      sh(5) = sh(5) - pt5*(sh(13) + sh(16) + sh(17))
      sh(6) = sh(6) - pt5*(sh(13) + sh(14) + sh(18))
      sh(7) = sh(7) - pt5*(sh(14) + sh(15) + sh(19))
      sh(8) = sh(8) - pt5*(sh(15) + sh(16) + sh(20))
!                                                                       
      return                                                            
      end subroutine                                                              
!**** new **********************************************************************
      subroutine igen(ia,m,iin)
      use mGlobaisEscalares
!
!.... program to read and generate integer nodal data
!
!        ia = input array
!         m = number of rows in ia
!         n = node number
!        ne = end node in generation sequence
!        ng = generation increment
!    
      integer :: m, ia(m,*), iin
!
      integer :: ib(m)
      integer :: n, ne, ng
      integer :: i
!
  100 continue
      read(iin,1000) n,ne,ng,(ib(i),i=1,m)

      if (n.eq.0) return

      if (ng.eq.0) then
         ne = n
         ng = 1
      else
         ne = ne - mod(ne-n,ng)
      endif
!
      do 200 i=n,ne,ng
!     call imove(ia(1,i),ib,m)
      ia(:,i)=ib
  200 continue
!
      go to 100
!
 1000 format(16i10)
      end subroutine
!     

!**** new **********************************************************************
      subroutine local(conectElem,x,xl,nen,nrowx,nrowxl)
!
!.... program to localize a global array
!
!        note: it is assumed nrowxl.le.nrowx
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: conectElem(*)
      integer :: nrowx, nrowxl, nen
      double precision :: x(nrowx,*),xl(nrowxl,*)
!
      integer :: i, j, node
!
      do 200 j=1,nen
      node = conectElem(j)
!
      do 100 i=1,nrowxl
      xl(i,j)= x(i,node)
  100 continue
!
  200 continue
!
      return
      end subroutine

         
! ! !******************************************************************************
! !       subroutine genelFaces(conecElem,nen,nelx, nely, nelz, iin)   
! !       use mGlobaisEscalares 
! !   
! !       implicit none
! ! !                                                                       
! ! !.... program to read and generate element faces and material numbers    
! ! !                                                                       
! ! !         conecElem(nen,numel) = element node numbers                         
! ! !         mat(numel)     = element material numbers                     
! ! !         nen            = number of element nodes (le.27)              
! ! !         n              = element number                               
! ! !         ng             = generation parameter                         
! ! !                                                                       
! !       integer :: conecElem(nen,*)
! !       integer :: nen, nelx, nely, nelz, iin
! ! !
! !       integer :: ng, n, m, nel, i
! !       integer :: condicao, condicao2                      
! ! !                                                                       
! !       read(iin,1000) n,m,(conecElem(i,1),i=1,nen),ng     
! ! !      write(*,1000) n,m,(conecElem(i,1),i=1,nen),ng
! ! 
! !       condicao=0
! !       condicao2=0
! ! 
! !       do nel=2, numel
! ! 
! !         if(condicao==0.and.condicao2==0) then
! !           do i=1, nen
! !             conecElem(i,nel)=conecElem(i,nel-1)+1
! !           end do
! !         else
! !           if(condicao==1.and.condicao2==0) then
! !             do i=1, nen
! !             if(i<=4) then
! !               conecElem(i,nel)=conecElem(i,nel-1)+(nelx+1)+1
! !             else
! !               conecElem(i,nel)=conecElem(i,nel-1)+1
! !             endif
! !             end do
! !           else
! !           if(condicao==1.and.condicao2==1) then
! !             do i=1, nen
! !             if(i<=4) then
! !                 conecElem(i,nel)=conecElem(i,nel-1)+(nelx+1)+nelx+(nelx*nely)+1
! !             else
! !                 conecElem(i,nel)=conecElem(i,nel-1)+(nelx*(nely+1))+(nely*(nelx+1))+1
! !             endif
! !             enddo 
! !           end if 
! !           end if
! !         end if
! !         
! !         if(mod(nel, nelx)==0) then
! !             condicao=1
! !         else 
! !             condicao=0
! !         end if
! ! 
! !         if(mod(nel, nelx*nely)==0) then
! !             condicao2=1
! !         else 
! !             condicao2=0
! !         end if
! ! 
! !       end do
! ! !                                                                       
! !   1000 format(16i10,10x,14i10)                                             
! ! !                                                                       
! !       end subroutine          
! ! 
! ! 
      end module
