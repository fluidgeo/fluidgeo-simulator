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
      subroutine ztest(a,n,lzero)
!
!.... program to determine if an array contains only zero entries
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer, intent(in)    :: n 
      real*8, intent(in)     :: a(n)
      logical, intent(inout) :: lzero
!
      integer :: i
!
      lzero = .true.
!
      do 100 i=1,n
      if (a(i).ne.0.0d0) then
         lzero = .false.
         return
      endif
  100 continue
!
      end     
!
!**** new **********************************************************************
!
      subroutine timing(tempoDeParede)
      implicit none     
!
!.... program to determine elapsed cpu time
!
!
      real*8, intent(inout) :: tempoDeParede
!
      character(LEN=8)  :: date
      character(LEN=10) :: minhaHora
      character(LEN=5)  :: zone
      integer,dimension(8) :: values
      integer ::  horas, minutos, segundos, milesimosSeg
!
      call date_and_time(date,minhaHora,zone,values);
!
         horas=values(5);      minutos=values(6); 
      segundos=values(7); milesimosSeg=values(8);    
!
      tempoDeParede = (60*horas+minutos)*60.0+segundos+milesimosSeg/1000.00
!
      return
      end
!
!=======================================================================
!
!       subroutine prtvB3D(nen,nsd,numel,conecNodaisElem,t0,velocLadal,ndofV, numLados, conecLadaisElem, numLadosElem, iunit)
! !
!       use mMalha, only: x
!       use mMalha, only: local
! !
!       implicit none
! !
! !     imprime campos vetoriais para o gnuplot ou para o matlab
! !
!       integer                   :: nen,numel,nsd, ndofV ,numLados,numLadosElem
!       real(8), dimension(ndofV,numLados) :: velocLadal
!       real(8)                   :: t0
!       integer                   :: conecNodaisElem(nen,numel),conecLadaisElem(numLadosElem,numel)
! !
!       integer  :: nel
!       real*8   :: xl(nsd,nen)
!       real(8)  :: xg,yg,zg,vcx,vcy,vcz
!       integer :: lado1,lado2,lado3,lado4,lado5,lado6
! !
!       integer :: iunit
! !
!       write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
! !
!       write(iunit,*)
! !
!       do nel=1,numel
! !
!         call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
!         xg = sum(xl(1,1:nen))/nen
!         yg = sum(xl(2,1:nen))/nen
!         zg = sum(xl(3,1:nen))/nen
! 
!         lado1 = conecLadaisElem(1,nel);  lado2 = conecLadaisElem(2,nel);
!         lado3 = conecLadaisElem(3,nel);  lado4 = conecLadaisElem(4,nel);
!         lado5 = conecLadaisElem(5,nel);  lado6 = conecLadaisElem(6,nel);
!         vcx = (velocLadal(1,lado2)+velocLadal(1,lado4))/2.0
!         vcy = (velocLadal(1,lado1)+velocLadal(1,lado3))/2.0
!         vcz = (velocLadal(1,lado5)+velocLadal(1,lado6))/2.0
!         write(iunit,"(6(f25.15,2x))") xg,yg,zg,vcx,vcy,vcz
! 
!       end do
! !
!       write(iunit,*)
! !
!       end subroutine
! !
!=======================================================================
!
      subroutine clear(a,m) 
! 
!.... program to clear a floating-point array 
! 
      implicit real*8 (a-h,o-z) 
! 
!.... remove above card for single-precision operation 
! 
      INTEGER :: M
      REAL*8 :: a(*) 
      integer :: i
! 
      do 100 i=1,m 
      a(i) = 0.0d0 
  100 continue 
! 
      return 
      end subroutine
!
!======================================================================
! 
         subroutine dividirTrabalho   (posI_, posF_, nth_, id_, inicio_, fim_)
         implicit none 
! 
         integer, intent(in)   :: posI_, posF_, nth_, id_
         integer, intent(out)  :: inicio_, fim_ 
! 
         integer :: tamanhoBloco, resto, elementos
         LOGICAL :: distribuirResto = .true.

         tamanhoBloco = int((posF_-posI_+1)/nth_)
         resto = mod((posF_-posI_+1),nth_)
         inicio_ = posI_ + id_*tamanhoBloco 
         fim_    = inicio_ + tamanhoBloco - 1

         if(resto > 0) then 
            if (distribuirResto) then
                if(id_<resto) then
                  if(id_ /= 0) inicio_ = inicio_+ id_
                  fim_    = fim_ + id_ + 1 
                else 
                  inicio_= inicio_ + resto
                  fim_   = fim_    + resto 
                end if
            else
               if(id_ == nth_ - 1) fim_ = fim_ + resto      
            end if 
         end if
         elementos = fim_ - inicio_ + 1
  ! write(*,'(a,4(a,i0))')rotulo,', id=', id_,  ' em dividirTrabalho, ' , inicio_,' ',  fim_, ' ', elementos
  !  write(*,'(4(a,i0))')' id=', id_,  ' em dividirTrabalho, ' , inicio_,' ',  fim_, ' ', elementos
 
         end subroutine dividirTrabalho 
!
!**** new **********************************************************************
!
      subroutine gerarLabel(label,tempo)
!
      use mGlobaisEscalares, only: tt
!
      implicit none


      character(LEN=21), intent(out) :: label
      real*8, intent(in) :: tempo
!
      character(LEN=21) :: labelAux, num     
      integer :: i
            
      write(num,'(f20.4)') tempo
!     labelAux="t="//ADJUSTL(num)
      labelAux="t="//num
      do i = 1, 21
         if(labelAux(i:i) .ne. ' ') then
            label(i:i) = labelAux(i:i)
         else
            label(i:i) = '0'
         end if
      end do
     
      end subroutine
