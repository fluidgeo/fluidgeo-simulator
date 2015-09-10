
  module mCSR
        integer*4,  allocatable :: AiPotencial(:), ApPotencial(:)
        integer*4,  allocatable :: LMstencilEqPotencial(:,:)

        integer*4,  allocatable :: AiFluxo(:),     ApFluxo(:)
        integer*4,  allocatable :: LMstencilEqFluxo(:,:)


        integer*4 :: posPonteiro, contPonteiro, posColunas, posCoef
        integer*4 :: lda,  nonzerosEst
        integer*4 :: nonzeros

      contains
!
!=======================================================================
!    
      subroutine criarPonteirosMatEsparsa_CSR(nsd, ndof, neq, numCoefPorLinha,  &
              conectsElem, listaDosElems, id, numConexoes, nen, numConexoesPorElem,  &
              nonzeros, simetria,tipo)

      use mGlobaisEscalares, only: optSolver

      implicit none 
!
      integer*4,  intent(in) :: nsd, numConexoes, nen, numConexoesPorElem,ndof
      integer*4,  intent(in) :: neq, nonzeros
      logical, intent(in) :: simetria
      integer*4,  intent(out) :: numCoefPorLinha
      integer*4 :: conectsElem(nen,*), listaDosElems(numConexoesPorElem,*)
      integer*4 :: id(ndof,*)
      character :: tipo

      write(*,*) " em subroutine criarPonteirosMatEsparsa_CSR(nsd, ndof, neq, num"
!
      if(tipo=='p') then 
         if(.not.allocated(LMstencilEqPotencial))allocate(LMstencilEqPotencial(neq,numCoefPorLinha)); 
         LMstencilEqPotencial=0

         call montarLmStencilNodal_CSR    (LMstencilEqPotencial,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, nen, numConexoesPorElem, neq, simetria)

         allocate(ApPotencial(neq+1));    ApPotencial=0  
         call montarPonteiroAp_CSR(ApPotencial, LMstencilEqPotencial, numCoefPorLinha, neq, nonzeros)
!

         allocate(AiPotencial(nonzeros)); AiPotencial=0
         call montarPonteiroAi_CSR(AiPotencial, LMstencilEqPotencial, numCoefPorLinha, neq, nonzeros)

!        call escreverEstruturaEsparsa( AlhsP, brhsP, ApPotencial, AiPotencial, neq, nonzeros)

       else

         if(.not.allocated(LMstencilEqFluxo))allocate(LMstencilEqFluxo(neq,numCoefPorLinha)); 
         LMstencilEqFluxo=0

         call montarLmStencilNodal_CSR    (LMstencilEqFluxo,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, nen, numConexoesPorElem, neq, simetria)

         allocate(ApFluxo(neq+1));    ApFluxo=0  
         call montarPonteiroAp_CSR(ApFluxo, LMstencilEqFluxo, numCoefPorLinha, neq, nonzeros)
!

         allocate(AiFluxo(nonzeros)); AiFluxo=0
         call montarPonteiroAi_CSR(AiFluxo, LMstencilEqFluxo, numCoefPorLinha, neq, nonzeros)

     !    call escreverEstruturaEsparsa( ApFluxo, AiFluxo, neq, nonzeros)
       endif

! 
      contains
      !
!**** new *************************************************************
!
      subroutine montarLmStencilNodal_CSR(LMstencilEq, listaDosElems, id, &
        conectsElem, numCoefPorLinha, ndof, numConexoes, nen, numConexoesPorElem, neq, simetria) 
!
      implicit none
!
      integer*4,  intent(in)    :: neq
      integer*4,  intent(in)    :: numCoefPorLinha, numConexoes,numConexoesPorElem,ndof,nen
      integer*4,  intent(inout) :: LMstencilEq(neq,numCoefPorLinha)
      integer*4,  intent(in)    :: listaDosElems(numConexoesPorElem,*), conectsElem(nen,*)
      integer*4,  intent(in)    :: id(ndof,*)
      logical    :: simetria
      !logical, intent(in)    :: simetria
!
      integer*4 :: nc, cont, nel, i, j, k, l, m,  numEq, dir
      logical :: propriaEq, difZero
      integer*4 :: numCoefAux, numViz, menorQueNumEq(ndof)
      integer*4,  allocatable :: LMstencilEqAux(:)
      integer*4 :: shift
!
    !  simetria=.false.
      write(*,*) " ++ subroutine montarLmStencilNodal_CSR(...., numCoefPorLinha = ", numCoefPorLinha
      shift=10 !8 !+10
      numCoefAux=numCoefPorLinha+shift ! + 12
      allocate(LMstencilEqAux(numCoefAux))

              numViz=nen
               numEq=0
         LMstencilEq=0
      LMstencilEqAux=0

      do nc=1, numConexoes

         do dir=1, ndof 
   
         if(id(dir,nc).ne.0) then
            numEq=numEq+1
            LMstencilEqAux=0
            cont=1
            propriaEq=.false.

            do k=1, numViz
               if(listaDosElems(k,nc).ne.0) then
                  nel=listaDosElems(k,nc)

                  do i=1, nen!4!numConexoesPorElem

                      menorQueNumEq=0
                      difZero=.false.
                      do j=1, ndof
                         if(id(j,conectsElem(i,nel)).ne.0) difZero=.true.
                         if(id(j,conectsElem(i,nel)).ne.numEq) menorQueNumEq(j)=1
                      enddo

                      if(difZero.eqv..true.) then !condicional verdadeira se algum dos id's for diferentes de zero

                        if(sum(menorQueNumEq)==ndof) then !condicional verdadeiro se os id's de todas as direcoes nao for o numero da equacao que estah sendo avaliada
                           
                           !inclui no LMStencil as equacoes que contribuem para a linha da matriz exceto o id que contem a propria equacao
                           do m=1, ndof

                              if(simetria.eqv..true.) then
                                 if(id(m,conectsElem(i,nel))>=numEq) then
                                    LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                    cont=cont+1
                                 endif
                              else
                                 if(id(m,conectsElem(i,nel)).ne.0) then
                                    LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                    cont=cont+1
                                 endif
                              endif
                           enddo

                        else
                          !inclui no LMStencil as equacoes onde o id contem a propria equacao 
                            if(propriaEq.eqv..false.) then
                               do m=1, ndof
                                  propriaEq=.true.
                                  if(simetria.eqv..true.) then
                                     if(id(m,conectsElem(i,nel))>=numEq) then
                                        LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                        cont=cont+1
                                     endif
                                  else
                                     if(id(m,conectsElem(i,nel)).ne.0) then
                                        LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                        cont=cont+1
                                     endif
                                  endif
                                  
                               enddo
                            endif

                        end if
                     end if
                  end do
!
               end if
            end do !k

         call ordenarLMstencil(LMstencilEqAux(:),numCoefAux)

         do i=1, numCoefAux-1
             if(LMstencilEqAux(i)==LMstencilEqAux(i+1)) LMstencilEqAux(i)=0
         end do

         call ordenarLMstencil(LMstencilEqAux(:),numCoefAux)

         LMstencilEq(numEq,1:numCoefPorLinha)=LMstencilEqAux(shift+1:numCoefAux)

        ! write(*,'(a,i0,a,81(1x,i0))'), "LmStencil ", numEq, " ->",LMstencilEq(numEq,1:numCoefPorLinha)
         end if

       end do !dir
      end do !gl

      deallocate(LMstencilEqAux)
      
      end subroutine montarLmStencilNodal_CSR      !
      !
      !----------------------------------------------------------------------
      !
      subroutine ordenarLMstencil(LMstencilEq,numCoefPorLinha)
      
      implicit none
!
      integer*4,  intent(in)    :: numCoefPorLinha      
      integer*4,  intent(inout) :: LMstencilEq(numCoefPorLinha)
!
      integer*8 :: menorEq, n, nn , tmp

      do n = 1, numCoefPorLinha
        menorEq=n
        do nn = n+1, numCoefPorLinha
           if(LMstencilEq(nn)<LMstencilEq(menorEq)) menorEq = nn
        end do

        if(n == menorEq) cycle
        tmp                  = LMstencilEq(n)
        LMstencilEq(n)       = LMstencilEq(menorEq)
        LMstencilEq(menorEq) = tmp
                
      enddo

      end subroutine ordenarLMstencil
      !
      !**** new *************************************************************
      !
      subroutine montarPonteiroAp_CSR (Ap, LMstencilEq, numCoefPorLinha, neq, nonzeros)

      implicit none

      integer*4 :: Ap(*)
      integer*4 :: nonzeros, neq
      integer*4  :: numCoefPorLinha
      integer*4 :: LMstencilEq(neq,numCoefPorLinha)
!
      integer*8 :: l,j

      ! Montando Ap

      call montarListaPonteiros(Ap, LMstencilEq,neq,numCoefPorLinha)

      ! Contando os valores nao nulos
      nonzeros=0
      do l=2, neq+1
            nonzeros=nonzeros+(Ap(l)-Ap(l-1))
      end do

       end subroutine montarPonteiroAp_CSR
      !
      !**** new *************************************************************
      !
      subroutine montarPonteiroAi_CSR (Ai, LMstencilEq, numCoefPorLinha, neq, nonzeros)

      implicit none

      integer*4,  intent(out) :: Ai(*)
      integer*4 :: nonzeros, neq
      integer*4 :: numCoefPorLinha
      integer*4 :: LMstencilEq(neq,numCoefPorLinha)
!
      integer*8 :: p

      call montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)

      end subroutine montarPonteiroAi_CSR
      !
      !**** new *************************************************************
      !
      subroutine montarListaPonteiros(Ap, LMstencilEq, neq, numCoefPorLinha)

      implicit none
      
      integer*4,  intent(out) :: Ap(*)
      integer*4,  intent(in) :: neq
      integer*4,  intent(in) :: numCoefPorLinha
      integer*4,  intent(in) :: LMstencilEq(neq,numCoefPorLinha)
      integer*4 :: LMstencilEqTemp(0:numCoefPorLinha)
!
      integer*8 :: n, l

      posPonteiro=0
      contPonteiro=0
      do l=1, neq
         posPonteiro=posPonteiro+1
         Ap(posPonteiro)=contPonteiro+1
         LMstencilEqTemp=0
         LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(l,:)

         do n = 1, numCoefPorLinha
            if(LMstencilEqTemp(n) == LMstencilEqTemp(n-1)) cycle
            contPonteiro= contPonteiro+1
         enddo    

         if(l==neq) then
            posPonteiro=posPonteiro+1
            Ap(posPonteiro)=contPonteiro+1
         end if
      end do

      end subroutine montarListaPonteiros
      !
      !**** new *************************************************************
      !
      subroutine montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)
!
      implicit none
      integer*4,  intent(in)    :: nonzeros, neq
      integer*4,  intent(inout) :: Ai(nonzeros)
      integer*4,  intent(in)    :: numCoefPorLinha
      integer*4,  intent(in)    :: LMstencilEq(neq,numCoefPorLinha)
!
      integer*4 :: LMstencilEqTemp(0:numCoefPorLinha)
      integer*4 :: i, n, posColunas
!
      posColunas=0
      do i=1, neq
         LMstencilEqTemp=0
         LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(i,:)
         do n = 1, numCoefPorLinha
             if(LMstencilEqTemp(n).ne.LMstencilEqTemp(n-1).and.LMstencilEqTemp(n).ne.0 ) then
                 posColunas=posColunas+1
                 Ai(posColunas)= LMstencilEqTemp(n)
             end if
         enddo
      end do
!
      end subroutine montarListaIndices
!
      end subroutine criarPonteirosMatEsparsa_CSR
!
      subroutine escreverSistemaCSR(alhs, brhs, Ap, Ai,  nonzeros, neq, nomeArq)
      implicit none 
      real*8,  intent(in)   :: Alhs(:), brhs(:)
      integer*4,  intent(in)  :: Ap(:), Ai(:)
      integer*8,  intent(in)  :: neq, nonzeros
      character (len=*)  :: nomeArq
!
      integer*8 :: i, j, k
      integer*8 :: luSist = 1836 
      character(len=40), parameter :: formatoEscritaA='(2(i0,1x),e23.16)'
      character(len=40), parameter :: formatoEscritaB='(1(i0,1x),e23.16)'
      real*8 :: t1, t2, t3, t4
!
      write(*,*) " em subroutine escreverMatrizCSR(alhs, brhs, Ap, Ai, neq, nonzeros)"
      call timing(t1)
      open(file=nomeArq, unit=luSist) 
!
      write(luSist,'(a)')'%% matriz A de coeficientes reais simetrica positiva definida ' 
      write(luSist,'(a)')'% produzida pelo metodo classico de galerkin para o metodo de elementos finitos '
      write(luSist,'(a)')'% armazenamento esparso CSR '
      write(luSist,'(a,3(i0,a))' )  '%  matriz ',neq, 'X',neq, ' com ', nonzeros, ' coefs diferentes de zero' 
      write(luSist,'(  3(i0,1x))' )  neq, neq, nonzeros 
      k = 1
      do i = 1, neq
          !write(luA, *) Ap(i) ,  Ap(i+1) 
           do j = Ap(i) ,  Ap(i+1) - 1  
                write(luSist, formatoEscritaA   ) i, Ai(j), alhs(k)
                !write(luA+luAux, '( 3(i10,2x), e20.10, 2x,i5)'     ) i, j,  Ai(j), alhs(k), k
                k = k + 1
           end do
      end do
!      write(luSist,'( (a,i0,a))') '% lado direito com ',neq,' elementos' 
!      write(luSist,'(  3(i10))' )  neq 
      i = 1
      write(luSist, * )  brhs(i), "   BRHS "
      do i = 2, neq
          write(luSist, * )  brhs(i)
         ! write(luSist, formatoEscritaB ) i, brhs(i)
      end do
      close(luSist)
      call timing(t2)
      write(*,*) " tempo de escrita =", t2 - t1

     end subroutine escreverSistemaCSR

!     
!     
     subroutine lerSistemaAlgMTXemCSR(alhs, brhs, Ap, Ai, neq, nonzeros,  nomeArq) 

     implicit none 
     real*8   ,intent(inout) :: alhs(:) ! matriz do sistema 
     real*8   ,intent(inout) :: brhs(:) ! lado direito do sistema
     integer*4,intent(inout) :: Ap(:) ! vetor apontadores para  inicio de uma linha
     integer*4,intent(inout) :: Ai(:) ! vetor com as colunas dos coefcientes de alhs
     integer*8,intent(inout) :: nonzeros
     integer*8,intent(inout) :: neq 
     character(len=*), intent(in) :: nomeArq
!
     integer   :: luSist = 1000, luSistLido = 2000 
     integer*8 :: i, j, k, ieq, ieqAnterior
     character(len=100) :: label 
     character(len=40), parameter :: formatoEscritaA='(2(i0,1x),e15.8)'
     character(len=40), parameter :: formatoEscritaB='(1(i0,1x),e15.8)'
!
     write(*,*) " em subroutine lerSistemaMTXemCSR(alhs, brhs, Ap, Ai, neq, nonzeros,",nomeArq,")"
     write(*,*) " FUNCIONANDO PARA MATRIZES SIMETRICAS COM ELEMENTOS FORNECIDOS POR LINHAS"
     open(file=nomeArq, unit=luSist) 

     read(luSist,'(a)' ) label; write(*,* ) label 
     read(luSist,'(a)' ) label; write(*,* ) label 
     read(luSist,'(a)' ) label; write(*,* ) label 
     read(luSist,'(a)' ) label; write(*,* ) label 
     read(luSist,* )  neq, neq, nonzeros
     write(*,'(a,3(i10))' ) "+++", neq, neq, nonzeros

! o arquivo que será lido estah escrito com os coeficientes de uma linha 
!      considerando matrix simetrica
! as linhas estao em ordem crescente de 1 ateh neq (num de equacoes) 

     Ap(:) = 1
     k = 0; ieqAnterior = 0
     do i = 1, nonzeros
          k = k + 1 ! contador de nao zeros lidos
          read(luSist, * ) ieq,  Ai(k), alhs(k)
          if(ieq > ieqAnterior) Ap(ieq) = k  
          ieqAnterior = ieq
     end do
     k = k + 1
     Ap(ieq+1) = k
   
     write(*,'( (a,i0,a))') 'lado direito com ',neq,' elementos' 
     do i = 1, neq
         read(luSist, * ) brhs(i)
     end do

      write(luSistLido,'(  3(i0,1x))' )  neq, neq, nonzeros 
      k = 1
      do i = 1, neq
           do j = Ap(i) ,  Ap(i+1) - 1  
                write(luSistLido, formatoEscritaA   ) i, Ai(j), alhs(k)
                k = k + 1
           end do
      end do
      do i = 1, neq
          write(luSistLido, formatoEscritaB ) i, brhs(i)
      end do

     close(luSist); 
     end subroutine lerSistemaAlgMTXemCSR

!     
     end module mCSR

