! !-------------------------------------------------------------------------------------     
! Descricao: dados fisicos do problema 
!-------------------------------------------------------------------------------------
      module mBlocoMacro
!
      implicit none 
! 
       integer :: iin, iecho, iinPassosPressaoBlocoMacro
!        integer :: iechoPressao, iechoProducao


       REAL*8, ALLOCATABLE :: solucao_BM(:,:), solucaoTmpAnt_BM(:,:), solucaoNaoLinearAnt_BM(:,:)      

       REAL*8, ALLOCATABLE :: f_BM(:,:), flux_BM(:,:)       !?
       REAL*8, ALLOCATABLE :: th_BM(:), rho_BM(:)           !?

       INTEGER :: IHARD_BM, NITERMAX_BM, NUSTEP_BM, NDIMC_BM
       REAL*8  :: TEMPO_BM, DTEMPO_BM
       REAL*8  :: flMassico_BM

       integer :: NESD_BM, NEE_BM, NEESQ_BM, NENP_BM, NUMITER_BM
       logical :: CHCKC_BM

       INTEGER :: NITER_BM, metodoLinear_BM !?
 
       integer :: ndof_BM, nlvect_BM
           
       INTEGER :: SUM_NUMITER_BM;   
       
       INTEGER, allocatable :: passosTempoParaGravarSaida(:)
       INTEGER :: numPassos
       
             
      contains
           
!***NEW  (Patricia - 23/08/13) **************************************************************
      SUBROUTINE INPUTadim(F,NDOF,NUMNP,J,NLVECT,IPRTIN,ID)
!
!.... PROGRAM TO READ, GENERATE AND WRITE NODAL INPUT DATA
!
!        F(NDOF,NUMNP,NLVECT) = PRESCRIBED FORCES/KINEMATIC DATA (J=0)
!                             = NODAL BODY FORCES(J=1)

      use mLeituraEscrita,   only:  printf, printd
      use mMalha,            only: genfl
      use mParametros,       only: condContAdim


      IMPLICIT NONE
      
!      IMPLICIT REAL*8 (A-H,O-Z)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION
!
      LOGICAL :: LZERO
      INTEGER :: NDOF, NUMNP, NLVECT, IPRTIN
      REAL*8  :: F(NDOF,NUMNP,*)
      INTEGER :: ID(NDOF,*)
!
      integer :: i, j, nlv, luEscrita

      CALL CLEAR(F,NLVECT*NUMNP*NDOF)
!
      luEscrita = iecho ! variavel de module

      DO 100 NLV=1,NLVECT
         CALL GENFL(F(1,1,NLV),NDOF, iin)
         DO I=1, NUMNP
            if ( ID(1,I).EQ.0 ) then
               F(1,I,NLVECT) = condContAdim;
            endif 
         ENDDO           
         CALL ZTEST(F(1,1,NLV),NDOF*NUMNP,LZERO)

      IF (IPRTIN.EQ.0) THEN
         IF (LZERO) THEN
            IF (J.EQ.0) WRITE(IECHO,1000) NLV
            IF (J.EQ.1) WRITE(IECHO,2000)
         ELSE
            IF (J.EQ.0) CALL PRINTF(F,NDOF,NUMNP,NLV)!,luEscrita)
            IF (J.EQ.1) &
     &      CALL PRINTD(' N O D A L  B O D Y  F O R C E S  ', &
     &                  F,NDOF,NUMNP,luEscrita)
         ENDIF
      ENDIF
!
  100 CONTINUE
!
      RETURN
 1000 FORMAT('1'//,' THERE ARE NO NONZERO PRESCRIBED FORCES AND ', &
     &    'KINEMATIC BOUNDARY CONDITIONS FOR LOAD VECTOR NUMBER ',I5)
 2000 FORMAT('1'//,' THERE ARE NO NONZERO NODAL BODY FORCES')
      END SUBROUTINE
      
!*** NEW ************************************************************************
!
      SUBROUTINE PLANMX(RHO,TH,C,NUMAT,ulEscrita) 
!
!     PROGRAM TO READ WRITE AND STORE PROPERTIES
!     FOR PLANE STRESS MIXED ELEMENTS
!
      IMPLICIT REAL*8    (A-H,O-Z)
!
!     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION
!
      real*8 :: RHO(*),TH(*),C(7,*)
      integer :: M, N, numat,ulEscrita
!
      DO 100 M=1,NUMAT
      IF (MOD(N,50).EQ.1) WRITE(ulEscrita,1000) NUMAT
!
      READ (iin,2000) PHI,BETA,RK 
      C(1,M)=PHI 
      C(2,M)=BETA
      C(3,M)=RK  
      WRITE(ulEscrita,3000) PHI,BETA,RK
      WRITE(*,3000) PHI,BETA,RK  
!
 100  CONTINUE
!
      RETURN
!
1000  FORMAT(8X,'MATERIAL SET DATA',//,2X,'NUMBER OF THE MATERIAL SET=',I1,//)
!
2000  FORMAT(8F10.0)
!
3000  FORMAT(3X,'PHI',3X,'BETA',3X,'RK' /3X,F10.3,3X,F10.3,3X,F10.3) 
      END SUBROUTINE
!
!**** NEW **********************************************************************
!
!
      SUBROUTINE GHARD(FORCE,X,NSD,NUMNP,IHARD)
!
!     PROGAM TO GENERATE NODAL BODY FORCE DATA FROM HARD-WIRED
!     FUNCTION F(X,Y)Á/C
!       IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT NONE
      
      REAL*8  :: FORCE(NSD,*), X(NSD,*)
      INTEGER :: NSD,NUMNP,IHARD
      INTEGER :: IPROB, I, J

      IPROB= IHARD-1
      DO  I=1,NUMNP
         DO J=1,NSD
            FORCE(J,I)=F( X(1,I), J, NSD, IPROB, IHARD)
         ENDDO
      ENDDO

      END SUBROUTINE
!
!*******************************************************************


      FUNCTION F(XV, J, NSD, IPRCB, IHARD)
!
!     HARD WIRED FUNCTION
      USE mGlobaisEscalares
!
      IMPLICIT REAL*8  (A-H,O-Z) 
!
      REAL*8 :: XV(NSD,*)
      INTEGER :: J, NSD, IPRCB, IHARD

      REAL*8 :: X, Y
!
      X=XV(1,1)
      Y=XV(2,1)
!
      GO TO (100,200,300,400) IPRCB
!
!   PROBLEM 1 - BODY FORCE ENGELMAN ET AL.
!
 100  CONTINUE
      IF(J.EQ.1) THEN
      F = G(X,Y) + Y - PT5
       ELSE
       F = - G(Y,X) + X - PT5
      END IF
      GO TO 1001
!
!   PROBLEM 2 - BODY FORCE ENGELMAN CUBIC PRESSURE
!
 200  CONTINUE
      PT125 = PT25 * PT5
      X2 = X * X
      X3 = X2* X
      Y2 = Y * Y
      Y3 = Y2* Y
      IF(J.EQ.1) THEN
        F = G(X,Y) + THREE * X2 * ( Y3 - PT125 )
       ELSE
       F = - G(Y,X) + THREE * Y2 * ( X3 - PT125 )
       END IF
       GO TO 1001
!
!   PROBLEM 3 - BODY FORCE ENGELMAN ET AL. WITH 100* PRESSURE CONTRIBUTION
!
 300    CONTINUE
        IF(J.EQ.1) THEN
          F = G(X,Y) + 100.D0 * (Y - PT5)
        ELSE
           F = - G(Y,X) + 100.D0 * (X - PT5)
        END IF
        GO TO 1001
!
! ERROR MESSAGE ON SCREEN AND STOP PROGRAM
!
 400    CONTINUE
        IHARD = IPRCB + 1
        print 2000, IHARD
        STOP
!
 1001   RETURN
 2000   FORMAT(10X,'WRONG IHARD =',I3)

        END function

!*******************************************************************************

         FUNCTION G(X,Y)
!
!  ENGELMAN FPART OF THE BODY FORCE
! 
      USE mGlobaisEscalares

      IMPLICIT REAL*8  (A-H,O-Z) 
!
!.... REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION
!
      real*8 :: X,Y
!
      X2  = X*X
      XM1 = X - ONE
!
      TEMP = SIX * X2 * XM1 * XM1
      TEMP1= TWO * Y * (Y-ONE) * (SIX*X2-SIX*X+ONE)
      TEMP = TEMP + TEMP1
!     G    = TEMP * 256D0 * (TWO * Y - ONE) 
      G    = TEMP * (TWO * Y - ONE) 
!
      RETURN
      END FUNCTION
!

!
!********************************************************************
!
      SUBROUTINE CHECKCONV_BM(D,DD,NDOF,NED,NUMNP,NUMITER,CHKCONV,NUSTEP)
!
!     Esta rotina faz o controle de convergencia das iteracoes
!     de newton. retorna .TRUE. if iterations converged
!
      IMPLICIT NONE
!       IMPLICIT REAL*8  (A-H,O-Z) 
      LOGICAL CHKCONV
      INTEGER :: NDOF,NED,NUMNP,NUMITER,NUSTEP
!
!     REMOVE ABOVE CARD FOR SI7NGLE PRECISION OPERATION
!
      REAL*8 :: D(NDOF,NUMNP),DD(NDOF,NUMNP),CTOL(2),DTOL
!
      INTEGER ::I,J, II
!   
      DTOL=1.0D-04
!
      CHKCONV=.FALSE.
!
      DO I=1,2
         CTOL(I)=0.0D+00
      ENDDO
!
      DO J=1,NUMNP
         CTOL(1)=CTOL(1)+(D(1,J)-DD(1,J))**2.0
         CTOL(2)=CTOL(2)+(0.5D0*(D(1,J)+DD(1,J)))**2.0
      ENDDO    
!
      CTOL(1)=DSQRT(CTOL(1))/(1.0+DSQRT(CTOL(2)))
      
      write(10,*) "iter",numiter 
      do II=1,numnp
         write(10,*) "ccw",ii,CTOL(1)
      end do 
      
      IF(CTOL(1).GT.DTOL) RETURN
!  
      CHKCONV=.TRUE.
!
      RETURN
!
      END SUBROUTINE

!*** NEW Patricia (setembro/2013) ********************************************************************** 
!    Para facilitar a comparação entre o modelo 1D e o 2D,
!    no modelo 2D, vamos imprimir o numero do no/2, ou seja, para o 2D, 
!    a numeraco de 1 a 201 corresponde a 1 a 401 de 2 em 2.
!    Modificada por Aline para o caso 2D

      SUBROUTINE PRINTSOL_BM(solucao,X,NUMNP,TEMPO,idx)
!
!     Esta rotina imprime a solucao do Bloco Macro

      use mGlobaisEscalares, only: dimModelo
      use mMalha,            only: nelx_BM, nely_BM
      use mLeituraEscrita,   only: iechoPressao
      
      IMPLICIT NONE
      
      REAL*8   TEMPO, solucao(1,*),X(2,*), L, pressaoPa  
      INTEGER  NUMNP, I, J,N, idx
      if (dimModelo=='1D') then
        DO I=1,NUMNP 
          N = I 
          WRITE(iechoPressao,210) N, TEMPO, X(1,N), X(2,N), solucao(1,N)     
        ENDDO      
      else ! 2D para criar plot no gnuplot
        write(iechoPressao,100) idx
        DO I=1,nelx_BM + 1
          DO J=1,nely_BM + 1
            N = I+(J-1)*(nelx_BM + 1)
            WRITE(iechoPressao,210) N, TEMPO, X(1,N), X(2,N), solucao(1,N)
!             WRITE(iechoPressao,220) X(1,N), X(2,N), solucao(1,N)

          ENDDO
          write(iechoPressao,*) ' '
        ENDDO        
          write(iechoPressao,*) ' '
          write(iechoPressao,*) ' ' 
       endif   
       
 100  FORMAT('#TEMPO',I5)
 200  FORMAT(4X,I5,10x,4(1PE15.8,2X))
 ! 4 espaços, inteiro max 5 posicoes, 10 espacos, 3 floats 8.2 com espaco de 2 entre eles
 210  FORMAT(4X,I5,10x,4(1PE15.8,2X))
 220  FORMAT(4(1PE15.8,2X))

      END subroutine
      
      
!*** NEW Patricia (setembro/2013) ********************************************************************** 
!    Para facilitar a comparação entre o modelo 1D e o 2D,
!    no modelo 2D, vamos imprimir o numero do no/2, ou seja, para o 2D, 
!    a numeraco de 1 a 201 corresponde a 1 a 401 de 2 em 2.
!    Modificada por Aline para o caso 2D
!    Modificado por Diego para a utilizacao do Python como visualizacao.

      SUBROUTINE PRINTSOL_BM2(solucao,X,NUMNP,TEMPO,idx)
!
!     Esta rotina imprime a solucao do Bloco Macro

      use mGlobaisEscalares, only: dimModelo
      use mMalha,            only: nelx_BM, nely_BM
      use mLeituraEscrita,   only: iechoPressao
      
      IMPLICIT NONE
      
      REAL*8   TEMPO, solucao(1,*),X(2,*), L, pressaoPa  
      INTEGER  NUMNP, I, J,N, idx, idr, idr2, idr3
      character*6  solP, idxStr, solP_BC
      
      idr = idx*17
      idr2 = idx*11
      idr3 = idx*13
      write(idxStr,'(i0)') idx
      solP = 'solP.'//idxStr
      OPEN(UNIT=idr, FILE= solP)
      
      OPEN(UNIT=idr2, FILE= 'solP_BC.'//idxStr)
      OPEN(UNIT=idr3, FILE= 'solP_C.'//idxStr)
      
      if (dimModelo=='1D') then
        DO I=1,NUMNP 
          N = I 
          WRITE(idr,210) N, TEMPO, X(1,N), X(2,N), solucao(1,N)     
        ENDDO      
      else ! 2D para criar plot no gnuplot
        DO I=1,nelx_BM + 1
          DO J=1,nely_BM + 1
            N = I+(J-1)*(nelx_BM + 1)
	      WRITE(idr,210) N, TEMPO, X(1,N), X(2,N), solucao(1,N)
          ENDDO
          write(iechoPressao,*) ' '
        ENDDO        
          write(iechoPressao,*) ' '
          write(iechoPressao,*) ' ' 
       endif
       
       DO I=nely_BM/2 + 1,nely_BM + 1, nely_BM/2
          DO J=1,nelx_BM + 1
            N = J+(I-1)*(nely_BM + 1)
            if (I .eq. (nely_BM + 1)) then
		WRITE(idr2,210) N, TEMPO, X(1,N), X(2,N), solucao(1,N)
	    else
	        WRITE(idr3,210) N, TEMPO, X(1,N), X(2,N), solucao(1,N)
	    endif
          ENDDO
       ENDDO   
          
 200  FORMAT(4X,I5,10x,4(1PE15.8,2X))
 ! 4 espaços, inteiro max 5 posicoes, 10 espacos, 3 floats 8.2 com espaco de 2 entre eles
 210  FORMAT(4X,I5,10x,4(1PE15.8,2X))
 220  FORMAT(4(1PE15.8,2X))
 223  FORMAT(4X,I5,10x,1(1PE15.8,2X))
      close(idr)
      close((11*idx))
      close(idr2)
      close(idr3)

      END subroutine

!*** NEW Diego (Abril/2016) ********************************************************************** 
!   Sub-rotina para impressao do resultado do campo de pressao 2D no
!   formato matricial.

      SUBROUTINE fieldP_BM(solucao,X,NUMNP,TEMPO,idx)
!
!     Esta rotina imprime a solucao do Bloco Macro

      use mGlobaisEscalares, only: dimModelo
      use mMalha,            only: nelx_BM, nely_BM
      use mLeituraEscrita,   only: iechoPressao, genUnit
      
      IMPLICIT NONE
      
      REAL*8   TEMPO, solucao(1,*),X(2,*), L, pressaoPa, teste
      INTEGER  NUMNP, I, J,N, idx, idr, idr2, idr3
      character*30  solP, idxStr, solP_BC
      logical :: fileCheck
      
      !idr = idx*17

      !call random_number(teste)
      !teste = teste*1d5
      
      !fileCheck = .true.
      !do while (fileCheck .eqv. .true.)
        !call random_number(teste)
        !teste = teste*1e5
        !idr = nint(teste)
        !inquire(unit=idr, opened=fileCheck)
      !enddo

      !write(*,*) nint(teste); stop
      idr = 1111
      inquire(unit=idr, opened=fileCheck)
      if (fileCheck) then
          write(*,*) "Unidade de escrita ja aberta: idr pressao"; stop
      endif
      !call genUnit(idr)!; write(*,*) idr; stop
      write(idxStr,'(i0)') idx
      solP = 'fieldP.'//idxStr
      OPEN(UNIT=idr, FILE=solP)
      
!      if (dimModelo=='2D') then
        DO I=1,nelx_BM + 1
          DO J=1,nely_BM + 1
            N = I+(J-1)*(nelx_BM + 1)
	      WRITE(idr,210) N, TEMPO, X(1,N), X(2,N), solucao(1,N)
          ENDDO
        ENDDO        
!       endif
       
          
 200  FORMAT(4X,I5,10x,4(1PE15.8,2X))
 ! 4 espaços, inteiro max 5 posicoes, 10 espacos, 3 floats 8.2 com espaco de 2 entre eles
 210  FORMAT(4X,I5,10x,4(1PE15.8,2X))
 220  FORMAT(4(1PE15.8,2X))
 223  FORMAT(4X,I5,10x,1(1PE15.8,2X))
      close(idr)
      !close((11*idx))
      !close(idr2)
      !close(idr3)

      END subroutine

!*** Diego (Abril/2015) ********************************************************************** 
!    Esta rotina tem como finalidade o calculo e armazenamento do fluxo massico
!    e a velocidade em todos os nos.
!    TODO: 1) Imprimir os resultados de fluxo massico para todos os nos, que dependem de P/Z(P).
!    2) Corrigir o calculo da velocidade em y (nao imprime).

      SUBROUTINE fieldV_BM(flux,X,NUMNP,TEMPO,idx)
!
!     Esta rotina imprime a solucao do Bloco Macro

      use mGlobaisEscalares, only: dimModelo
      use mGlobaisArranjos,  only: k_bm, knp_bm
      use mMalha,            only: nelx_BM, nely_BM
      use mParametros,       only: constK_BM, constMu, p_Ref, widthBlocoMacro, tamBlocoMacro
      use mLeituraEscrita,   only: genUnit
      
      IMPLICIT NONE
      
      REAL*8   TEMPO, flux(2,*),X(2,*), L, pressaoPa, vx, vy, v, Lx, Ly
      INTEGER  NUMNP, I, J,N, idx, idrx, idry, idrv, idrvx, idrvy, I2, J2, N2
      CHARACTER*30  idxStr, gradPx, gradPy, solVelocity, solVelocity_x, solVelocity_y
      REAL*8  :: Keff
      logical :: fileCheck
            
!       Abrindo os arquivos para saída

      write(idxStr,'(i0)') idx
      
      !idrx = 11*idx
      idrx = 2221
      inquire(unit=idrx, opened=fileCheck)
      if (fileCheck) then
          write(*,*) "Unidade de escrita ja aberta: idrx"; stop
      endif
      !call genUnit(idrx)!; write(*,*) idrx
      gradPx = 'gradPx.'//idxStr
      OPEN(UNIT=idrx, FILE= gradPx)
      
      !idry = 12*idx
      idry = 2222
      inquire(unit=idry, opened=fileCheck)
      if (fileCheck) then
          write(*,*) "Unidade de escrita ja aberta: idry"; stop
      endif
      !call genUnit(idry)!; write(*,*) idry
      gradPy = 'gradPy.'//idxStr
      OPEN(UNIT=idry, FILE= gradPy)
      
      if (dimModelo=='1D') then
!	idrv = 13*idx
      idrv = 2223
      inquire(unit=idrv, opened=fileCheck)
      if (fileCheck) then
          write(*,*) "Unidade de escrita ja aberta: idrv"; stop
      endif
!      call genUnit(idrv)!; write(*,*) idrv
	solVelocity = 'solVelocity.'//idxStr
	OPEN(UNIT=idrv, FILE= solVelocity)
      else
!	idrvx = 14*idx
      idrvx = 2224
      inquire(unit=idrvx, opened=fileCheck)
      if (fileCheck) then
          write(*,*) "Unidade de escrita ja aberta: idrvx"; stop
      endif
!      call genUnit(idrvx)!; write(*,*) idrvx
	solVelocity_x = 'solVelocity_x.'//idxStr
	OPEN(UNIT=idrvx, FILE= solVelocity_x)
      
!	idrvy = 17*idx
      idrvy = 2225
      inquire(unit=idrvy, opened=fileCheck)
      if (fileCheck) then
          write(*,*) "Unidade de escrita ja aberta: idrvy"; stop
      endif
!      call genUnit(idrvy)!; write(*,*) idrvy
	solVelocity_y = 'solVelocity_y.'//idxStr
	OPEN(UNIT=idrvy, FILE= solVelocity_y)
      endif
!       ***********************************
      
      Lx    = widthBlocoMacro
      Ly    = tamBlocoMacro
      if (dimModelo=='1D') then
        DO I=1,NUMNP 
          N = I 
	    WRITE(idrx,210) N, TEMPO, X(1,N), X(2,N), flux(1,N)
! 	    v = (constK_BM/constMu)*(p_Ref/Lx)*flux(1,N)
	    v = (constK_BM/constMu)*flux(1,N)
	    WRITE(idrv,210) N, TEMPO, X(1,N), X(2,N), v
        ENDDO      
      else
        DO I=1,nelx_BM + 1
          DO J=1,nely_BM + 1
            N = I+(J-1)*(nelx_BM + 1)
            Keff = Knp_bm(N)
            !write(*,*) "N2 = ", N2
            !write(*,*) "Keff = ", Keff(N2)
            if (Keff .le. 1.0d-50) Keff = 0.0
	      WRITE(idrx,210) N, TEMPO, X(1,N), X(2,N), flux(1,N)
	      WRITE(idry,210) N, TEMPO, X(1,N), X(2,N), flux(2,N)
	      vx = Keff*flux(1,N)
	      vy = Keff*flux(2,N)
	      WRITE(idrvx,210) N, TEMPO, X(1,N), X(2,N), vx
	      WRITE(idrvy,210) N, TEMPO, X(1,N), X(2,N), vy
          ENDDO
          !stop
          !write(1574,*) "N2 = ", N2; 
          !write(1577,*) "Keff(N2) = ", Keff; 
        ENDDO
        !stop        
       endif   
      !stop
  
 200  FORMAT(4X,I5,10x,4(1PE15.8,2X))
! 4 espaços, inteiro max 5i posicoes, 10 espacos, 3 floats 8.2 com espaco de 2 entre eles
 210  FORMAT(4X,I5,10x,4(1PE15.8,2X))
 220  FORMAT(4(1PE15.8,2X))
 
      close(idrx)
      close(idry)
      close(idrv)
      close(idrvx)
      close(idrvy)

      END subroutine

!
!*** NEW Patricia (set/2013) ********************************************************************** 
!    Modificada por Aline para o caso 2D

      SUBROUTINE calcularFluxoMassico(FLUX, solucao, X,  NUSTEP, TEMPO, NUMNP, DT, fluxoAnt, volAcumulado)
      
      use mMalha,            only: NSD_BM, nelx_BM, nely_BM
      use mGlobaisEscalares, only: dimModelo
      use mCoeficientes,     only: calcularZ_P
      use mParametros,       only: p_Ref, tamBlocoMacro, widthBlocoMacro, constK_BM, T, R_, constMu, M_m
      use mParametros,       only: gasTotalKg, gasRecuperavelKg, areaContatoBlocoMacroFratura, gasProduzidoKg
      use mLeituraEscrita,   only: iechoProducao
      
      
      IMPLICIT NONE
      
      REAL*8   FLUX(2*NDOF_BM,*), solucao(NDOF_BM,*), X(NSD_BM,*), TEMPO, DT
      INTEGER  NUSTEP, NUMNP, I, N
      REAL*8  :: volAcumulado, fluxoAnt

      REAL *8  Lx, Ly, Z(nely_BM), pStar(nely_BM), flMassico(NSD_BM,nely_BM), volInst, &
     &         gradP(NSD_BM,nely_BM), volAcumuladoKg
      
!     fluxoAnt é iniciado com 0
      Lx    = widthBlocoMacro
      Ly    = tamBlocoMacro      ! BM_y
      
      DO I = 1,nely_BM                ! Laco nos nós verticais
        N = (I-1)*(nelx_BM+1) + 1
        pStar(I) = solucao(1, N)
        gradP(1,I) = FLUX(1,N)*p_Ref/Lx
        gradP(2,I) = FLUX(2,N)*p_Ref/Ly
        call calcularZ_P(pStar(I)*p_Ref, Z(I))
 
      ENDDO

        DO I = 1,nely_BM
          flMassico(1,I) = - (pStar(I)*p_ref*constK_BM*M_m*gradP(1,I))/(Z(I)*R_*T*constMu)  ! Kg / (m^2 s)
          flMassico(2,I) = - (pStar(I)*p_ref*constK_BM*M_m*gradP(2,I))/(Z(I)*R_*T*constMu)  ! Kg / (m^2 s)
        ENDDO  

      volInst        = (fluxoAnt + flMassico(1,1))/2.0D0 * DT      ! regra do trapezio
      volAcumulado   = volAcumulado + volInst             ! kg/m^2
      fluxoAnt       = flMassico(1,1)
      volAcumuladoKg = volAcumulado * areaContatoBlocoMacroFratura
         
      gasProduzidoKg = volAcumuladoKg;
   
      write(iechoProducao,'(i8,7e15.7)') NUSTEP, TEMPO, DT, flMassico(1,1), volAcumulado,      &
     &                                   volAcumuladoKg, volAcumuladoKg/gasRecuperavelKg, &
     &                                   volAcumuladoKg/gasTotalKg
      
     END SUBROUTINE
!
!*** NEW Patricia (set/2013) ********************************************************************** 
!    Modificada por Aline para o caso 2D
!    Ratificado por Diego em Maio/2015.

      SUBROUTINE calcularFluxoMassico2D(FLUX, solucao, X,  NUSTEP, TEMPO, NUMNP, DT, fluxoAnt, volAcumulado)
      
      use mMalha,            only: NSD_BM, nelx_BM, nely_BM
      use mGlobaisEscalares, only: dimModelo
      use mGlobaisArranjos,  only: k_bm, Knp_bm
      use mCoeficientes,     only: calcularZ_P
      use mParametros,       only: p_Ref, tamBlocoMacro, widthBlocoMacro, constK_BM, T, R_, constMu, M_m
      use mParametros,       only: gasTotalKg, gasRecuperavelKg, areaContatoBlocoMacroFratura, gasProduzidoKg
      use mLeituraEscrita,   only: iechoProducao
      
      
      IMPLICIT NONE
      
      REAL*8   FLUX(2*NDOF_BM,*), solucao(NDOF_BM,*), X(NSD_BM,*), TEMPO, DT
      INTEGER  NUSTEP, NUMNP, I, N, idx, idxr, I2, N2
      REAL*8  :: volAcumulado, fluxoAnt
      CHARACTER*30  idxStr, prod

      REAL*8  Lx, Ly, Z(nely_BM), pStar(nely_BM), flMassico(NSD_BM,nely_BM), volInst, &
     &         gradP(NSD_BM,nely_BM), volAcumuladoKg, flMassico_soma(NSD_BM,1)
      REAL*8  :: Keff(nely_BM), flInterface            
!     fluxoAnt é iniciado com 0
      Lx    = widthBlocoMacro
      Ly    = tamBlocoMacro      ! BM_y
      flInterface = 0.0d0
      
        DO I = 1,nely_BM                ! Laco nos nós verticais
	  N = (I-1)*(nelx_BM+1) + 1
	  pStar(I) = solucao(1, N)
	  gradP(1,I) = FLUX(1,N)!*p_Ref/Lx
      gradP(2,I) = FLUX(2,N)!*p_Ref/Ly
        Keff(I)  = Knp_bm(N)
        !write(*,*) "Keff = ", Keff(I2)
        if (Keff(I) .le. 1.0d-50) Keff(I) = 0.0
	  call calcularZ_P(pStar(I), Z(I))
	  !flMassico(1,I) = - (pStar(I)*Keff(I2)*M_m*gradP(1,I))/(Z(I)*R_*T)  ! Kg / (m^2 s)
	  !flMassico(2,I) = - (pStar(I)*Keff(I2)*M_m*gradP(2,I))/(Z(I)*R_*T)  ! Kg / (m^2 s)
	  flMassico(1,I) = - (pStar(I)*Keff(I)*M_m*gradP(1,I))/(Z(I)*R_*T)  ! Kg / (m^2 s)
	  flMassico(2,I) = - (pStar(I)*Keff(I)*M_m*gradP(2,I))/(Z(I)*R_*T)  ! Kg / (m^2 s)
            if (flMassico(1,I) .le. 1.0d-30) flMassico(1,I) = 0.0
            if (flMassico(2,I) .le. 1.0d-30) flMassico(2,I) = 0.0
            flInterface = flInterface + &
      &     flMassico(1,I)*(areaContatoBlocoMacroFratura/nely_bm)
        ENDDO
      !stop 
      volInst        = (fluxoAnt + flInterface)/2.0D0 * DT      ! regra do trapezio
      ! volAcumulado é iniciado com 0
      volAcumulado   = volAcumulado + volInst             ! kg
      fluxoAnt       = flInterface
      volAcumuladoKg = volAcumulado
         
      gasProduzidoKg = volAcumuladoKg;
   
      write(iechoProducao,'(i8,7e15.7)') NUSTEP, TEMPO, DT, flMassico(1,1), volAcumulado,      &
     &                                   volAcumuladoKg, volAcumuladoKg/gasRecuperavelKg, &
     &                                   volAcumuladoKg/gasTotalKg
     
!      close(idxr)
      
      END SUBROUTINE

!
!*** Diego (maio/2015) ********************************************************************** 
!    Essa rotina calcula o fluxo mássico em todos os nós.

      SUBROUTINE fieldJ_BM(FLUX, solucao, solucaoant, X, TEMPO, DT, NUMNP, idx)
      
      use mMalha,            only: NSD_BM, nelx_BM, nely_BM, NUMNP_BM
      use mGlobaisEscalares, only: dimModelo
      use mGlobaisArranjos,  only: k_bm, knp_bm 
      use mCoeficientes,     only: calcularZ_P
      use mParametros,       only: p_Ref, tamBlocoMacro, widthBlocoMacro, constK_BM, T, R_, constMu, M_m, phi_BM
      use mParametros,       only: gasTotalKg, gasRecuperavelKg, areaContatoBlocoMacroFratura, gasProduzidoKg
      use mLeituraEscrita,   only: genUnit
      
      
      IMPLICIT NONE
      
      REAL*8   solucao(NDOF_BM,*), solucaoant(NDOF_BM,*), TEMPO, flux(2,*), X(2,*), DT
      INTEGER  NUSTEP, NUMNP, I, N, J, idxf, idyf, idx, idxr, idxe, N2
      INTEGER  idxk

      REAL*8  Lx, Ly, Z(NUMNP_BM), pStar(NUMNP_BM), flMassico(NSD_BM,NUMNP_BM), residJ, &
     &         gradP(NSD_BM,NUMNP_BM),pStar_ant(NUMNP_BM),residt,Zant(NUMNP_BM),resid(NUMNP_BM),&
     &	       erro(NUMNP_BM)

     CHARACTER*30  idxStr, nodeFlux_x, nodeFlux_y, residueFlux_x, k_n
     logical :: fileCheck

      Lx    = widthBlocoMacro
      Ly    = tamBlocoMacro      ! BM_y
      
!       Abrindo os arquivos para saída
      write(idxStr,'(i0)') idx
      
!      idxf = 19*idx
      idxf = 3331
      inquire(unit=idxf, opened=fileCheck)
      if (fileCheck) then
          write(*,*) "Unidade de escrita ja aberta: idxf"; stop
      endif
!      call genUnit(idxf); !write(*,*) idr; stop
      nodeFlux_x = 'nodeFlux_x.'//idxStr
      OPEN(UNIT=idxf, FILE= nodeFlux_x)
      
!      idxr = 23*idx
!      residueFlux_x = 'residueFlux_x.'//idxStr
!      OPEN(UNIT=idxr, FILE= residueFlux_x)

!      idyf = 7*idx
      idyf = 3332
      inquire(unit=idyf, opened=fileCheck)
      if (fileCheck) then
          write(*,*) "Unidade de escrita ja aberta: idyf"; stop
      endif
!      call genUnit(idyf); !write(*,*) idr; stop
      nodeFlux_y = 'nodeFlux_y.'//idxStr
      OPEN(UNIT=idyf, FILE= nodeFlux_y)

!      idxk = 11*idx
      idxk = 3333
      inquire(unit=idxk, opened=fileCheck)
      if (fileCheck) then
          write(*,*) "Unidade de escrita ja aberta: idxk"; stop
      endif
!      call genUnit(idxk); !write(*,*) idr; stop
      k_n = 'condHydraulic.'//idxStr
      OPEN(UNIT=idxk, FILE= k_n)
!       *********************************************
      
      DO I=1,nelx_BM + 1
          DO J=1,nely_BM + 1
            N = I+(J-1)*(nelx_BM + 1)
            pStar(N) = solucao(1, N)
            pStar_ant(N) = solucaoant(1, N)
            gradP(1,N) = FLUX(1,N)!*p_Ref/Lx
            gradP(2,N) = FLUX(2,N)!*p_Ref/Ly
            call calcularZ_P(pStar(N), Z(N))
            call calcularZ_P(pStar_ant(N), Zant(N))
            flMassico(1,N) = - (pStar(N)*Knp_BM(N)*M_m*gradP(1,N))/(Z(N)*R_*T)  ! Kg / (m^2 s)
            flMassico(2,N) = - (pStar(N)*Knp_BM(N)*M_m*gradP(2,N))/(Z(N)*R_*T)  ! Kg / (m^2 s)
            WRITE(idxf,222) N, TEMPO, X(1,N), X(2,N), flMassico(1,N)  
            WRITE(idyf,222) N, TEMPO, X(1,N), X(2,N), flMassico(2,N)  
            ! Condutividade hidráulica
            WRITE(idxk,222) N, TEMPO, X(1,N), X(2,N), (pStar(N)*Knp_BM(N)*M_m)/(Z(N)*R_*T)  
          ENDDO
      ENDDO
      
      ! Para a conservação de massa espacial local (colocar como saída da subrotina)
!      DO I=1,nely_BM + 1
!          DO J=2,nelx_BM
!            N = J+(I-1)*(nely_BM + 1)
!            residJ = (flMassico(1,N+1) - flMassico(1,N-1))/((X(1,N+1)-X(1,N-1))) ! Ver o espaçamento da malha depois
!            residt = ((phi_BM*M_m)/(R_*T*DT))*((pStar(N)/Z(N))-(pStar_ant(N)/Zant(N)))
!            resid(N) = residt + residJ
!          ENDDO
!      ENDDO
      
!      DO I=2,nelx_BM
!          DO J=1,nely_BM + 1
!            N = I+(J-1)*(nelx_BM + 1)
!            WRITE(idxr,222) N, TEMPO, X(1,N), X(2,N), resid(N)
!          ENDDO
!      ENDDO

      222  FORMAT(4X,I5,10x,4(1PE15.8,2X))
      close(idxf)
      close(idyf)
      close(idxk)
      
      END SUBROUTINE

      
!-------------------------------------------------------------------------------------     
      SUBROUTINE montarSistema_BM (NED, NDOF,trU, &
          trUAnt, DTEMPO)
!-------------------------------------------------------------------------------------     
      
!
!.... PROGRAM TO CALCULATE STIFNESS MATRIX AND FORCE ARRAY FOR THE
!        STOKE'S DISPLACEMENT  ELEMENT AND
!        ASSEMBLE INTO THE GLOBAL LEFT-HAND-SIDE MATRIX
!        AND RIGHT-HAND SIDE VECTOR

      use mGlobaisEscalares, only : dimModelo, ntype, numat_BM, npint_BM, nicode_BM, iprtin, nrowsh_BM, lambda, mu
      use mGlobaisArranjos,  only : mat_BM, c_BM, grav_BM, bf_BM, phi_n, coupling_mode, phi_n0, k_bm, celast
      USE mLeituraEscrita,   only : printd, prntel
      use mMalha,            only : numel_BM, numnp_BM, nsd_BM, nen_BM, genfl, genel, local
      use mMalha,            only : conecNodaisElem_BM, x_BM
      use mFuncoesDeForma,   only : oneshl, oneshg,  shlq, SHGQ
      use mAlgMatricial,     only : colht, kdbc, addrhs, addnsl, addlhs, idiag_BM, lm_BM, alhs_BM
      use mAlgMatricial,     only : brhs_BM, dlhs_BM
      use mCoeficientes,     only : calcularZ_P
      use mParametros,       only : p_Ref, tamBlocoMacro,widthBlocoMacro, constMu, phi_BM, constK_BM, R_, T, K_re, K_abs
      use mParametros,       only : fraVol_BM, M_m, alpha_r, Kbulk, k_s, beta_r, Sg, Sw !, Se, Swr, Sgr, VL, PL
      use mParametros,       only : VL, PL, phi_n0_Num, K_absRnd
      
!
      IMPLICIT NONE
      
      INTEGER :: NDOF, NED
      REAL*8  :: trU(2*ndof,numnp_bm), trUAnt(2*ndof,numnp_bm),DIVU,DIVU_ANT,rho_sc, rho_r
      REAL*8  :: DTEMPO
            
      REAL*8  :: UU, UUP, GRXUU, GRYUU, tamElem, Se, Swr, Sgr 
      REAL*8  :: Z_UU, Z_UUP, R_UU, R_UUP, K_tmp, Keff_tmp
      REAL*8  :: M, ALPHA, CVRHO, GF1, GF2
      REAL*8  :: H, H2, EE, XNI, C1, T0, RK, RKCON, DIX, DIY, DIN, DJN, DJX, DJY
      
      integer :: iopt

      real*8  :: xl(nesd_BM,nen_BM), dl(ned,nen_BM),dpl(ned,nen_BM)
      real*8  :: ul(2*ned,nen_BM),dul(2*ned,nen_BM)
      real*8, dimension(nrowsh_BM,nenp_BM,npint_BM) :: SHLP, SHGP
      real*8, dimension(nrowsh_BM,nen_BM, npint_BM) :: SHL, SHG
      real*8  :: det(npint_BM), detp(npint_BM), W(npint_BM), WP(npint_BM)
      REAL*8  :: flNoAnt, flNoPost
      
      logical :: lsym

      REAL*8 :: fonteMassaDeBlocoParaBlocoMacro
      
      integer :: i,j,k,l,nel
      
! 
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION 
! 
      LOGICAL DIAG,QUAD,ZERODL
      real*8 :: ELEFFM(NEE_BM,NEE_BM),ELRESF(NEE_BM) ! bidu 20ago 2015

      shl = 0.0
      if(dimModelo=='1D') then
          call oneshl(shl,w,npint_BM,nen_BM)
      end if
      if(dimModelo=='2D')  then
          CALL SHLQ(SHL,W,NPINT_BM,NEN_BM)
        !       CALL SHLQ(SHLP,WP,NPINT,NENP)
      end if 	
!
!      CONSISTENT MATRIX
!   
      DIAG = .FALSE.
      
      Sw        = 0.10d0
      Sg        = 1.00d0 - Sw
      Swr       = 0.05               
      Sgr       = 0.05
      Se        = (Sw-Swr)/(1-Swr-Sgr);
      K_re      = (1-Se)*(1-Se)*(1-Se**2);
      K_abs     = constK_BM*(3.0d0-phi_n0_Num)/(2.0d0*phi_n0_Num);

      lambda = (celast(1)*celast(2))/((1.0+celast(2))*(1.0-2.0*celast(2)))
      mu = (celast(1))/(2.0*(1.0+celast(2)))
      rho_r = celast(3)*10.0**3.0d0
      rho_sc = 0.67
      Kbulk = (lambda + 2.0/3.0*mu)
      alpha_r = 1.0d0 - Kbulk/k_s
      !write(*,*) rho_r; write(*,*) VL; write(*,*) PL; stop
      write(*,*) "Kbulk =", Kbulk;! stop 1 
!      write(*,*) "k_s =", k_s
!      write(*,*) "E =", celast(1), "nu =", celast(2)
      write(*,*) "alpha_r =", alpha_r; stop

      DO 500 NEL=1,NUMEL_BM
!
!      CLEAR STIFFNESS MATRIX AND FORCE ARRAY
!
      CALL CLEAR(ELEFFM,NEESQ_BM)
      CALL CLEAR(ELRESF,NEE_BM)           
!
!      LOCALIZE COORDINATES AND DIRICHLET B.C.
!  
      CALL LOCAL(conecNodaisElem_BM(1,NEL),X_BM,XL,NEN_BM,NSD_BM, NESD_BM)
      CALL LOCAL(conecNodaisElem_BM(1,NEL),solucao_BM,DL,NEN_BM,NDOF,NED)
      CALL LOCAL(conecNodaisElem_BM(1,NEL),solucaoTmpAnt_BM,DPL,NEN_BM,NDOF,NED)
      call local(conecNodaisElem_BM(1,nel),trU,ul,nen_bm,2*ndof,2*ndof)
      call local(conecNodaisElem_BM(1,nel),trUAnt,dul,nen_bm,2*ndof,2*ndof)
      !stop 
      M = MAT_BM(NEL)
! 
         shg  = 0.0
         if(dimModelo=='1D') then 
            call oneshg(xl,det,shl,shg,nen_BM,npint_BM,nsd_BM,1,nel,1) !BIDU
         end if
         if(dimModelo=='2D') then
            QUAD = .TRUE.
            CALL SHGQ(XL,DET,SHL,SHG,NPINT_BM,NEL,QUAD,NEN_BM)        
         end if
!
!....... FORM STIFFNESS MATRIX
!
       tamElem  = X_BM(2,NEL+1)-X_BM(2,NEL);
       beta_r = ((alpha_r-phi_n0(nel))/k_s + (alpha_r**2.0)/Kbulk)	! Computing total compressibility (Diego, set/2015)
!      K_absRnd(NEL) = constK_BM*(3.0d0-phi_n(nel))/(2.0d0*phi_n(nel));
      DO 400 L=1,NPINT_BM
         C1=DET(L)*W(L) 
         UU=0.D00
         UUP=0.D00
         GRXUU=0.D00
         GRYUU=0.D00
         DIVU=0.d0
         DIVU_ANT=0.d0
         DO J=1,NEN_BM
            UU   =UU   +SHG(3,J,L)*DL(1,J)  
            UUP  =UUP  +SHG(3,J,L)*DPL(1,J)
            GRXUU=GRXUU+SHG(1,J,L)*DL(1,J)
            GRYUU=GRYUU+SHG(2,J,L)*DL(1,J)
            DIVU=DIVU+SHG(1,J,L)*UL(1,J)+SHG(2,J,L)*UL(2,J)
            DIVU_ANT=DIVU_ANT+SHG(1,J,L)*DUL(1,J)+SHG(2,J,L)*DUL(2,J)
         ENDDO
         !stop
         CALL calcularZ_P(UU, Z_UU)         
         CALL calcularZ_P(UUP,Z_UUP)
         R_UU      = Sg*phi_n(nel)*M_m/(Z_UU *  R_ * T)
         R_UUP     = Sg*phi_n(nel)*M_m/(Z_UUP *  R_ * T) 
         k_bm(nel) = K_re*K_absRnd(nel)/constMu*(2.0d0*phi_n(nel))/(3.0d0-phi_n(nel));   
         Keff_Tmp  = K_bm(nel)*UU*M_m/(R_*T*Z_UU) !* UU/(R_*T*Z_UU) * 1.0 !* fc! + D * ( Gamma*rhoL/H + p*rhoL * GammaDivH_Linha );

         if (coupling_mode .eq. "oneway") then
         DO J=1,NEN_BM
            DJN=SHG(3,J,L)*C1
            DJX=SHG(1,J,L)*C1
            DJY=SHG(2,J,L)*C1

            ELRESF(NED*J)=ELRESF(NED*J)+DJN*R_UUP*UUP  &
!                &        +djn*rho_sc*rho_r*VL*PL/((PL+UU)**2.0d0)*UUP  &
	&		+ beta_r*djn*M_m*(UU/Z_UU)*UUP/(R_*T)	! Diego, 1-way (set/2015) 
         ENDDO        
!
!.... ELEMENT STIFFNESS
!
         DO J=1,NEN_BM
            DJX=SHG(1,J,L)*C1
            DJY=SHG(2,J,L)*C1
            DJN=SHG(3,J,L)*C1
            DO I=1,NEN_BM
               DIX=SHG(1,I,L)
               DIY=SHG(2,I,L) 
               DIN=SHG(3,I,L)
       
              ELEFFM(NED*J,NED*I) = ELEFFM(NED*J,NED*I)                                  &
       &                            + R_UU*DJN*DIN                                        &
       &                            + Keff_tmp*DTEMPO & 
!       &                            + constK_BM*DTEMPO*UU*M_m/(R_*T*Z_UU) & 
       &                            * (DIX*DJX+DIY*DJY) &
!       &                            * (DIX*DJX+DIY*DJY)/constMu &
!                &        +djn*rho_sc*rho_r*VL*PL/((PL+UU)**2.0d0)*din  &
       &			    + beta_r*djn*M_m*(UU/Z_UU)*(DIN)/(R_*T)	! Diego, 1-way (set/2015)
       

       
!             write(*,*) "======================="
!             write(*,*) "phi", phi_F
!             write(*,*) "K",   constK_F
!             write(*,*) "mu",  constMu
!             write(*,*) "DT",  DTEMPO_F
!             write(*,*) "UU",  UU
!             write(*,*) "UUP", UUP
!             write(*,*) "Z_UU",  Z_UU
!             write(*,*) "Z_UUP", Z_UUP
!             write(*,*) "======================="

            ENDDO
         ENDDO
         
         else if (coupling_mode .eq. "twoway") then
         DO J=1,NEN_BM
            DJN=SHG(3,J,L)*C1
            DJX=SHG(1,J,L)*C1
            DJY=SHG(2,J,L)*C1

        !write(*,*)  beta_r,djn,M_m,UU,Z_UU,UUP,R_*T    ! Diego, 2-way (dec/2015)
           
            ELRESF(ndof*J)=ELRESF(ndof*J)+DJN*R_UUP*UUP  &
                  &+ beta_r*djn*M_m*(UU/Z_UU)*UUP/(R_*T) &
!                &        +djn*rho_sc*rho_r*VL*PL/((PL+UU)**2.0d0)*UUP  &
       &+ alpha_r*djn*M_m*(UU/Z_UU)*((DIVU_ANT-(alpha_r/Kbulk)*UUP)-(DIVU-(alpha_r/Kbulk)*UU))/(R_*T)
! 	&		+ ((alpha_r**2.0)/Kbulk)*djn*M_m*(UU/Z_UU)*(UU-UUP)/(R_*T)  
            !write(19,*) (trU(nel)-trUAnt(nel))
         ENDDO        
!
!.... ELEMENT STIFFNESS
!
         DO J=1,NEN_BM
            DJX=SHG(1,J,L)*C1
            DJY=SHG(2,J,L)*C1
            DJN=SHG(3,J,L)*C1
            DO I=1,NEN_BM
               DIX=SHG(1,I,L)
               DIY=SHG(2,I,L) 
               DIN=SHG(3,I,L)
       
              ELEFFM(ndof*J,ndof*I) = ELEFFM(ndof*J,ndof*I)                                  &
       &                            + R_UU*DJN*DIN                                        &
!       &                            + constK_BM*DTEMPO*UU*M_m/(R_*T*Z_UU) & 
       &                            + Keff_tmp*DTEMPO & 
       &                            * (DIX*DJX+DIY*DJY) &
!       &                            * (DIX*DJX+DIY*DJY)/constMu &
!                &        +djn*rho_sc*rho_r*VL*PL/((PL+UU)**2.0d0)*din  &
       &                            +beta_r*djn*M_m*(UU/Z_UU)*(DIN)/(R_*T)! &  ! Diego, 2-way (dec/2015)
!       &                            + alpha_r*djn*M_m*(UU/Z_UU)*trU(nel)/(R_*T)
       
       !write(*,*)  "R_UU,DJN,DIN, constK_BM,DTEMPO,UU,M_m,R_,T,Z_UU,DIX,DJX,DIY,DJY,constMu "
       !write(*,*)  R_UU,DJN,DIN, constK_BM,DTEMPO,UU,M_m,R_,T,Z_UU,DIX,DJX,DIY,DJY,constMu 
       !write(*,*)  "beta_r,djn,M_m,UU,Z_UU,DIN,R_*T"; 
       !write(*,*)  beta_r,djn,M_m,UU,Z_UU,DIN,R_*T; !stop


            ENDDO
         ENDDO
         else  
         DO J=1,NEN_BM
            DJN=SHG(3,J,L)*C1
            DJX=SHG(1,J,L)*C1
            DJY=SHG(2,J,L)*C1

            ELRESF(NED*J)=ELRESF(NED*J)+DJN*R_UUP*UUP 
!                &        +djn*rho_sc*rho_r*VL*PL/((PL+UU)**2.0d0)*UUP 
         ENDDO        
!
!.... ELEMENT STIFFNESS
!
         DO J=1,NEN_BM
            DJX=SHG(1,J,L)*C1
            DJY=SHG(2,J,L)*C1
            DJN=SHG(3,J,L)*C1
            DO I=1,NEN_BM
               DIX=SHG(1,I,L)
               DIY=SHG(2,I,L) 
               DIN=SHG(3,I,L)
       
              ELEFFM(NED*J,NED*I) = ELEFFM(NED*J,NED*I)                                  &
       &                            + R_UU*DJN*DIN                                        &
!       &                            + constK_BM*DTEMPO*UU*M_m/(R_*T*Z_UU) & 
       &                            + Keff_tmp*DTEMPO & 
!       &                            * (DIX*DJX+DIY*DJY)/constMu
!                &        +djn*rho_sc*rho_r*VL*PL/((PL+UU)**2.0d0)*din  &
       &                            * (DIX*DJX+DIY*DJY)
       

       
!             write(*,*) "======================="
!             write(*,*) "phi", phi_F
!             write(*,*) "K",   constK_F
!             write(*,*) "mu",  constMu
!             write(*,*) "DT",  DTEMPO_F
!             write(*,*) "UU",  UU
!             write(*,*) "UUP", UUP
!             write(*,*) "Z_UU",  Z_UU
!             write(*,*) "Z_UUP", Z_UUP
!             write(*,*) "======================="

            ENDDO
         ENDDO
         endif
!
  400 CONTINUE
         !stop

!      COMPUTATION OF DIRICHLET B.C. CONTRIBUTION
!   
      CALL ZTEST(DL,NEE_BM,ZERODL)
!
!       IF(.NOT.ZERODL) CALL KDBC(ELEFFM,ELRESF,DL,NEE,LM(1,1,NEL),NEL)
      IF(.NOT.ZERODL) CALL KDBC2(ELEFFM,ELRESF,DL,NEE_BM,LM_BM(1,1,NEL),NEL)

!
!.... ASSEMBLE ELEMENT STIFNESS MATRIX AND FORCE ARRAY INTO GLOBAL
!        LEFT-HAND-SIDE MATRIX AND RIGHT-HAND SIDE VECTOR

       lsym=.true.
!        CALL ADDLHS(ALHSP,ELEFFM,idiagP,LMP(1,1,NEL),NEE,DIAG,lsym)
       CALL ADDNSL(ALHS_BM,DLHS_BM,ELEFFM,IDIAG_BM,LM_BM(1,1,NEL),NEE_BM,DIAG)
!
       CALL ADDRHS(BRHS_BM,ELRESF,LM_BM(1,1,NEL),NEE_BM)
!                
 500   CONTINUE

      RETURN
      END SUBROUTINE

!-------------------------------------------------------------------------------------     
      SUBROUTINE montarSistema_T (NED, NDOF,trU, &
          trUAnt, DTEMPO)
!-------------------------------------------------------------------------------------     

!    Montagem do sistema do problema de consolidacao de Terzaghi.
!    (Diego, jan/2016)

!
!.... PROGRAM TO CALCULATE STIFNESS MATRIX AND FORCE ARRAY FOR THE
!        STOKE'S DISPLACEMENT  ELEMENT AND
!        ASSEMBLE INTO THE GLOBAL LEFT-HAND-SIDE MATRIX
!        AND RIGHT-HAND SIDE VECTOR

      use mGlobaisEscalares, only : dimModelo, ntype, numat_BM, npint_BM, nicode_BM, iprtin, nrowsh_BM, lambda, mu
      use mGlobaisArranjos,  only : mat_BM, c_BM, grav_BM, bf_BM, phi_n, coupling_mode, phi_n0, k_bm, celast
      USE mLeituraEscrita,   only : printd, prntel
      use mMalha,            only : numel_BM, numnp_BM, nsd_BM, nen_BM, genfl, genel, local
      use mMalha,            only : conecNodaisElem_BM, x_BM
      use mFuncoesDeForma,   only : oneshl, oneshg,  shlq, SHGQ
      use mAlgMatricial,     only : colht, kdbc, addrhs, addnsl, addlhs, idiag_BM, lm_BM, alhs_BM
      use mAlgMatricial,     only : brhs_BM, dlhs_BM
      use mCoeficientes,     only : calcularZ_P
      use mParametros,       only : p_Ref, tamBlocoMacro,widthBlocoMacro, constMu, phi_BM, constK_BM, R_, T, K_re, K_abs
      use mParametros,       only : fraVol_BM, M_m, alpha_r, Kbulk, k_s, beta_r, Sg, Sw !, Se, Swr, Sgr, VL, PL
      use mParametros,       only : VL, PL
      
!
      IMPLICIT NONE
      
      INTEGER :: NDOF, NED
      REAL*8  :: trU(2*ndof,numnp_bm), trUAnt(2*ndof,numnp_bm),DIVU,DIVU_ANT,rho_sc, rho_r
      REAL*8  :: DTEMPO
            
      REAL*8  :: UU, UUP, GRXUU, GRYUU, tamElem, Se, Swr, Sgr 
      !REAL*8  :: Z_UU, Z_UUP, R_UU, R_UUP, K_tmp, Keff_tmp
      !REAL*8  :: M, ALPHA, CVRHO, GF1, GF2
      REAL*8  :: H, H2, EE, XNI, C1, T0, RK, RKCON, DIX, DIY, DIN, DJN, DJX, DJY
      
      integer :: iopt

      real*8  :: xl(nesd_BM,nen_BM), dl(ned,nen_BM),dpl(ned,nen_BM)
      real*8  :: ul(2*ned,nen_BM),dul(2*ned,nen_BM)
      real*8, dimension(nrowsh_BM,nenp_BM,npint_BM) :: SHLP, SHGP
      real*8, dimension(nrowsh_BM,nen_BM, npint_BM) :: SHL, SHG
      real*8  :: det(npint_BM), detp(npint_BM), W(npint_BM), WP(npint_BM)
      REAL*8  :: flNoAnt, flNoPost
      
      logical :: lsym

      REAL*8 :: fonteMassaDeBlocoParaBlocoMacro
      
      integer :: i,j,k,l,nel
      
! 
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION 
! 
      LOGICAL DIAG,QUAD,ZERODL
      real*8 :: ELEFFM(NEE_BM,NEE_BM),ELRESF(NEE_BM) ! bidu 20ago 2015
      real*8  :: biotM, Ku, Kfluid, Cf, pi, const_a

      pi = 3.141592653589793
      const_a = 1.0d0
  
      shl = 0.0
      if(dimModelo=='1D') then
          call oneshl(shl,w,npint_BM,nen_BM)
      end if
      if(dimModelo=='2D')  then
          CALL SHLQ(SHL,W,NPINT_BM,NEN_BM)
        !       CALL SHLQ(SHLP,WP,NPINT,NENP)
      end if 	
!
!      CONSISTENT MATRIX
!   
      DIAG = .FALSE.
      
      Sw        = 0.30d0
      Sg        = 1.00d0 - Sw
      Swr       = 0.05               
      Sgr       = 0.05
      Se        = (Sw-Swr)/(1-Swr-Sgr);
      K_re      = (1-Se)*(1-Se)*(1-Se**2);
      !K_abs     = constK_BM; 
      K_abs     = constMu*constK_BM/(K_re)*(3.0d0-sum(phi_n0)/numel_bm)/(2.0d0*sum(phi_n0)/numel_bm);
      !k_tmp     = K_re*K_abs/constMu*(2.0d0*phi_n(nel))/(3.0d0-phi_n(nel));
      !k_tmp     = K_re*K_abs/constMu*(2.0d0*sum(phi_n)/numel_bm)/(3.0d0-sum(phi_n)/numel_bm);   
      !Keff_Tmp  =  K_tmp !* UU/(R_*T*Z_UU) * 1.0 !* fc! + D * ( Gamma*rhoL/H + p*rhoL * GammaDivH_Linha );

      lambda = (celast(1)*celast(2))/((1.0+celast(2))*(1.0-2.0*celast(2)))
      mu = (celast(1))/(2.0*(1.0+celast(2)))
      rho_r = celast(3)*10.0**3.0d0
      rho_sc = 0.67
      Kbulk = (lambda + 2.0/3.0*mu)
      alpha_r = 1.0d0 - Kbulk/k_s
      !write(*,*) rho_r; write(*,*) VL; write(*,*) PL; stop
      !write(*,*) "Kbulk =", Kbulk
      !write(*,*) "k_s =", k_s
      !write(*,*) "alpha_r =", alpha_r; stop

      DO 500 NEL=1,NUMEL_BM
!
!      CLEAR STIFFNESS MATRIX AND FORCE ARRAY
!
      CALL CLEAR(ELEFFM,NEESQ_BM)
      CALL CLEAR(ELRESF,NEE_BM)           
!
!      LOCALIZE COORDINATES AND DIRICHLET B.C.
!  

      CALL LOCAL(conecNodaisElem_BM(1,NEL),X_BM,XL,NEN_BM,NSD_BM, NESD_BM)
      CALL LOCAL(conecNodaisElem_BM(1,NEL),solucao_BM,DL,NEN_BM,NDOF,NED)
      CALL LOCAL(conecNodaisElem_BM(1,NEL),solucaoTmpAnt_BM,DPL,NEN_BM,NDOF,NED)
      !stop 
! 
         shg  = 0.0
         if(dimModelo=='1D') then 
            call oneshg(xl,det,shl,shg,nen_BM,npint_BM,nsd_BM,1,nel,1) !BIDU
         end if
         if(dimModelo=='2D') then
            QUAD = .TRUE.
            CALL SHGQ(XL,DET,SHL,SHG,NPINT_BM,NEL,QUAD,NEN_BM)        
         end if
!
!....... FORM STIFFNESS MATRIX
!
       Kfluid = 1.0d00
       biotM = 1.0d0/(((alpha_r-0.25)/k_s) + 0.25/Kfluid)
       Ku = Kbulk + (alpha_r**2.0d0)*biotM
       Cf = constK_BM*biotM*((Kbulk+4.0d0/3.0d0*mu)/(Ku+4.0d0/3.0d0*mu))
       !write(*,*) "Cf =", Cf; write(*,*) "Ku = ", Ku; write(*,*) "alpha = ", alpha_r; 
       !write(*,*) "Kbulk =", Kbulk; write(*,*) "M =", biotM; stop

      DO 400 L=1,NPINT_BM
         C1=DET(L)*W(L) 
         UU=0.D00
         UUP=0.D00
         GRXUU=0.D00
         GRYUU=0.D00
         DO J=1,NEN_BM
            UU   =UU   +SHG(3,J,L)*DL(1,J)  
            UUP  =UUP  +SHG(3,J,L)*DPL(1,J)
            GRXUU=GRXUU+SHG(1,J,L)*DL(1,J)
            GRYUU=GRYUU+SHG(2,J,L)*DL(1,J)
         ENDDO

         DO J=1,NEN_BM
            DJN=SHG(3,J,L)*C1
            DJX=SHG(1,J,L)*C1
            DJY=SHG(2,J,L)*C1

            ELRESF(NED*J)=ELRESF(NED*J)+DJN*UUP !+ djn*(dsin(pi*XL(1,J))*dsin(pi*XL(2,J)))
!                &        +djn*rho_sc*rho_r*VL*PL/((PL+UU)**2.0d0)*UUP 
         ENDDO       
!
!.... ELEMENT STIFFNESS
!
         DO J=1,NEN_BM
            DJX=SHG(1,J,L)*C1
            DJY=SHG(2,J,L)*C1
            DJN=SHG(3,J,L)*C1
            DO I=1,NEN_BM
               DIX=SHG(1,I,L)
               DIY=SHG(2,I,L) 
               DIN=SHG(3,I,L)
       
              ELEFFM(NED*J,NED*I) = ELEFFM(NED*J,NED*I)                                  &
       &                            + DJN*DIN                                        &
!       &                            + DTEMPO*const_a*(DIX*DJX+DIY*DJY)
       &                            + DTEMPO*Cf*(DIX*DJX+DIY*DJY)
       

       
!             write(*,*) "======================="
!             write(*,*) "phi", phi_F
!             write(*,*) "K",   constK_F
!             write(*,*) "mu",  constMu
!             write(*,*) "DT",  DTEMPO_F
!             write(*,*) "UU",  UU
!             write(*,*) "UUP", UUP
!             write(*,*) "Z_UU",  Z_UU
!             write(*,*) "Z_UUP", Z_UUP
!             write(*,*) "======================="

            ENDDO
         ENDDO
!
  400 CONTINUE
         !stop

!      COMPUTATION OF DIRICHLET B.C. CONTRIBUTION
!   
      CALL ZTEST(DL,NEE_BM,ZERODL)
!
!       IF(.NOT.ZERODL) CALL KDBC(ELEFFM,ELRESF,DL,NEE,LM(1,1,NEL),NEL)
      IF(.NOT.ZERODL) CALL KDBC2(ELEFFM,ELRESF,DL,NEE_BM,LM_BM(1,1,NEL),NEL)

!
!.... ASSEMBLE ELEMENT STIFNESS MATRIX AND FORCE ARRAY INTO GLOBAL
!        LEFT-HAND-SIDE MATRIX AND RIGHT-HAND SIDE VECTOR

       lsym=.true.
!        CALL ADDLHS(ALHSP,ELEFFM,idiagP,LMP(1,1,NEL),NEE,DIAG,lsym)
       CALL ADDNSL(ALHS_BM,DLHS_BM,ELEFFM,IDIAG_BM,LM_BM(1,1,NEL),NEE_BM,DIAG)
!
       CALL ADDRHS(BRHS_BM,ELRESF,LM_BM(1,1,NEL),NEE_BM)
!                
 500   CONTINUE

      RETURN
      END SUBROUTINE
!************************************************************** 
!*** NEW ******************************************************
!*** SUBROUTINA QUE COMPUTA : grad* p*
!************************************************************** 
      subroutine calcularGradiente_BM(NED, NDOF)
!
      use mMalha,            only: conecNodaisElem_BM, x_BM
      use mGlobaisEscalares, only: dimModelo, ntype, numat_BM, npint_BM, nicode_BM, iprtin, nrowsh_BM
      use mMalha,            only: numel_BM, numnp_BM, nsd_BM, nen_BM, genfl, genel, local
      use mFuncoesDeForma,   only: oneshl, oneshg, shlq, SHGQ

      IMPLICIT NONE

      integer :: NED, NDOF

      LOGICAL QUAD

      REAL*8 :: GRADPX, GRADPY
     
      INTEGER   :: I, J, L, NEL, NELG
  
      real*8 ::  HNM(NUMNP_BM)

      real*8 :: xl(nesd_BM,nen_BM), dl(ned,nen_BM)
      real*8 :: det(npint_BM), detp(npint_BM), W(npint_BM), WP(npint_BM)
      real*8, dimension(nrowsh_BM,nenp_BM,npint_BM) :: SHLP, SHGP
      real*8, dimension(nrowsh_BM,nen_BM, npint_BM) :: SHL, SHG

!
!      GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
!
      shl = 0.0d0
      if(dimModelo=='1D') call oneshl(shl,w,npint_BM,nen_BM)
      if(dimModelo=='2D') CALL SHLQ(SHL,W,NPINT_BM,NEN_BM)
      shlp = shl
      
      CALL CLEAR(flux_BM,2*NUMNP_BM)       

      HNM=0.D0
!
!.... BEGIN LOOP OVER THE MESH ELEMENTS
!     
      DO NEL=1,NUMEL_BM
!
!....    LOCALIZE UNKNOWNS AND COORDINATES

         CALL LOCAL(conecNodaisElem_BM(1,NEL),x_BM,XL,NEN_BM,NSD_BM,NESD_BM)              
         CALL LOCAL(conecNodaisElem_BM(1,NEL),solucao_BM,DL,NEN_BM,NDOF,NED)

!....    EVALUATE GLOBAL SHAPE FUNCTION 

         shgp = 0.0
         shg  = 0.0
         if(dimModelo=='1D') &
            call oneshg(xl,det,shl,shg,nen_BM,npint_BM,nsd_BM,1,nel,1) !BIDU
         if(dimModelo=='2D') &
             CALL SHGQ(XL,DET,SHL,SHG,NPINT_BM,NEL,QUAD,NEN_BM)        
         SHGP=SHG     
!
!....    COMPUTE FIRST DERIVATIVES OF NODAL PRESSURE
!
!....    BEGIN LOOP OVER ELEMENT NODES
!    
         DO L=1,NEN_BM
!
!           COMPUTE FIRST DERIVATIVES OF NODAL DISPLACEMENT
!           Para 1D nao tem derivada em Y
            GRADPX=0.D0
            GRADPY=0.D0
            
            ! SHG(1,1),SHG(1,2), SHG(1,3) e SHG(1,4) avaliadas no ponto de integracao L
            DO J=1,NEN_BM
               GRADPX = GRADPX + SHGP(1,J,L)*DL(1,J)
               GRADPY = GRADPY + SHGP(2,J,L)*DL(1,J)
            ENDDO
            ! fornece GRADPX e GRADPY avaliados no ponto de integracao L
            ! por enquanto nao tem grady porque a condicao de contorno nao tem gradiente em y 

!           STORAGING STRESSES ON GLOBAL ARRAY
            NELG=conecNodaisElem_BM(L,NEL)
            
            ! esse menos vem de fluxo = -gradP ??? 
            flux_BM(1,NELG)= flux_BM(1,NELG) -  1.D0*GRADPX
            flux_BM(2,NELG)= flux_BM(2,NELG) -  1.D0*GRADPY
!
            HNM(NELG)=HNM(NELG)+1.D0 ! usado para cálculo da média
            !write(*,*) "NELG = ", NELG, "HNM = ", HNM(NELG)!; stop

         ENDDO !    END OF LOOP OVER ELEMENT NODES                   
!
      ENDDO ! END OF LOOP OVER MESH ELEMENTS
!
!       COMPUTING AVERAGE STRESSES AT NODAL POINTS
!
      DO J=1,NUMNP_BM
         DO I=1,2
             flux_BM(I,J)=flux_BM(I,J)/HNM(J) 
         END DO
      END DO   
 
11000 FORMAT(I2)
12000 FORMAT(2X,14(1PE15.8,2X)) 
!
      end subroutine
!*** NEW ******************************************************
!*** SUBROUTINA QUE COMPUTA : K nos elementos (Diego, mar/2016)
!************************************************************** 
      subroutine kElem_BM(NED, NDOF)
!
      use mGlobaisArranjos,  only: K_BM, Knp_BM    
      use mMalha,            only: conecNodaisElem_BM, x_BM
      use mGlobaisEscalares, only: dimModelo, ntype, numat_BM, npint_BM, nicode_BM, iprtin, nrowsh_BM
      use mMalha,            only: numel_BM, numnp_BM, nsd_BM, nen_BM, genfl, genel, local
      use mFuncoesDeForma,   only: oneshl, oneshg, shlq, SHGQ

      IMPLICIT NONE

      integer :: NED, NDOF

      LOGICAL QUAD

!      REAL*8 :: Kmean(NUMNP_BM)
     
      INTEGER   :: I, J, L, NEL, NELG
  
      real*8 ::  HNM(NUMNP_BM)

      real*8 :: xl(nesd_BM,nen_BM), dl(ned,nen_BM)
      real*8 :: det(npint_BM), detp(npint_BM), W(npint_BM), WP(npint_BM)
      real*8, dimension(nrowsh_BM,nenp_BM,npint_BM) :: SHLP, SHGP
      real*8, dimension(nrowsh_BM,nen_BM, npint_BM) :: SHL, SHG

      Knp_BM = 0.0d0

      HNM=0.D0
!
!.... BEGIN LOOP OVER THE MESH ELEMENTS
!     
      DO NEL=1,NUMEL_BM
!
!....    LOCALIZE UNKNOWNS AND COORDINATES

         CALL LOCAL(conecNodaisElem_BM(1,NEL),x_BM,XL,NEN_BM,NSD_BM,NESD_BM)              
!
!....    BEGIN LOOP OVER ELEMENT NODES
!    
         DO L=1,NEN_BM
!
!           STORAGING STRESSES ON GLOBAL ARRAY
            NELG=conecNodaisElem_BM(L,NEL)
            
            ! esse menos vem de fluxo = -gradP ??? 
            Knp_BM(NELG)= Knp_BM(NELG) +  1.D0*K_BM(NEL)
!
            HNM(NELG)=HNM(NELG)+1.D0 ! usado para cálculo da média
            !write(*,*) "NELG = ", NELG, "HNM = ", HNM(NELG)!; stop

         ENDDO !    END OF LOOP OVER ELEMENT NODES                   
!
      ENDDO ! END OF LOOP OVER MESH ELEMENTS
!
!       COMPUTING AVERAGE STRESSES AT NODAL POINTS
!
!      open(UNIT=2,FILE="k1.dat")
      DO I=1,NUMNP_BM
         Knp_BM(I)=Knp_BM(I)/HNM(I) 
!         write(2,*) Knp_BM(I)
      END DO   
!      close(2)
 
11000 FORMAT(I2)
12000 FORMAT(2X,14(1PE15.8,2X)) 
!
      end subroutine
!
!**** NEW **********************************************************************
      SUBROUTINE KDBC2(ELEFFM,ELRESF,DL,NEE,LM,NEL)
!
!.... PROGRAM TO ADJUST LOAD VECTOR FOR PRESCRIBED DISPLACEMENT
!     BOUNDARY CONDITION
      use mGlobaisEscalares, only: zero
!
      IMPLICIT NONE !REAL*8 (A-H,O-Z)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION
!
      real*8 :: ELEFFM(NEE,*),ELRESF(*),DL(*), VAL
      integer :: LM(*), nel, nee
!
      integer :: i, j, L

!
!    THIS VERSION OF KDBC IS ONLY VALID FOR THERMOELASTIC CONS.
!
!
      DO 200 J=1,NEE
         L=LM(J) 
         VAL=DL(J)
!
         IF(L.GT.0) GO TO 200
         IF(VAL.EQ.ZERO) GO TO 200
!
         DO 100 I=1,NEE
            ELRESF(I)=ELRESF(I)-ELEFFM(I,J)*VAL
100   CONTINUE
!
200   CONTINUE
!
      RETURN
      END SUBROUTINE

!**** new **********************************************************************
      subroutine preprocessadorBlocoMacro()
!
       use mGlobaisArranjos,  only: title_BM
       use mGlobaisEscalares, only: dimModelo, exec, iprtin, ndofD
       use mAlgMatricial,     only: id_BM, idiag_BM, neq_BM, nalhs_BM
       use mMalha,            only: nen_BM, nsd_BM, numel_BM, numnp_BM, x_BM, numLados_BM
       use mMalha,            only: nelx_BM, nely_BM, nelz_BM
       use mMalha,            only: conecNodaisElem_BM, listaDosElemsPorNo_BM
       use mLeituraEscrita,   only: echo, abrirArquivo
       use mLeituraEscrita,   only: leiturageracaocoordenadas, leituracodigoscondcontorno, leituraValoresCondContorno
       use mParametros,       only: inicializarParametrosBlocoMacro
       use mcoeficientes,     only: calcularQuantidadeGasBlocoMacro
       
      implicit none

      character(len=21) :: label
      integer :: n, i
      CHARACTER*200 ::  linhaAux
      integer :: naoUsado
      logical :: exist
      
      
!.... input phase
!

      call echo(iin)
           
!
      ! 1D-2D
      read(iin,1000)  title_BM
      write(*,*) iin, title_BM
      if (title_BM(1).eq.'*end') return
!
      read(iin,'(5i10, 10a)') exec, iprtin, nsd_BM, nen_BM, ndof_BM, label  
      write(*,'(5i10, a)')    exec, iprtin, nsd_BM, nen_BM, ndof_BM, label  
      dimModelo=trim(adjustl(label))
      write(*,*) '...',dimModelo,'...'
      read(iin,'(7i10)') numnp_BM, numel_BM, numLados_BM, nelx_BM, nely_BM, nelz_BM
      read(iin,'(3i10)') nlvect_BM, naoUsado, naoUsado
!
      write(iecho,1000) title_BM 
      write(iecho,3000) exec, iprtin, nsd_BM, nen_BM, dimModelo
      write(iecho,4000) numnp_BM, numLados_BM, numnp_BM, &
                        ndof_BM, naoUsado, nlvect_BM, naoUsado
      WRITE(IECHO,5000) nelx_BM,nely_BM,nelz_BM

      nesd_BM = 2            ! +++ quem é este e pq não é lido de arquivo de entrada?

!       numBlocos = numel_BM   ! +++ aqui mesmo devo inicializar ncelr?

!---------------------------------------------------------------------

     call inicializarParametrosBlocoMacro()
     
     call calcularQuantidadeGasBlocoMacro()
          
     call alocarMemoriaBlocoMacro()
!
!.... input coordinate data
!
     call leituraGeracaoCoordenadas(x_BM,nsd_BM,numnp_BM, iin, iecho, iprtin)
!
!.... input boundary condition data and establish equation numbers
!
     call leituraCodigosCondContorno(id_BM,ndof_BM,numnp_BM,n,iin,iecho,iprtin)
     neq_BM=n        
     allocate(idiag_BM(neq_BM));  idiag_BM=0

     print*, "ndof_BM=,", ndof_BM, "neq_BM=", neq_BM
!
!.... input nodal force and prescribed kinematic boundary-value data
!
!      write(*,*) 'F(1)', f_F(1,1), 'F(2)', f_F(1,2)
     IF (NLVECT_BM.GT.0) CALL INPUTadim(f_BM,NDOF_BM,NUMNP_BM,0,NLVECT_BM, IPRTIN, id_BM)
!
!.... input element data
!
      print*, "topologiaMalhaSistEquacoes"
      ! esta rotina substitui HEAT1 (não é identica, mas faz seu papel
      call topologiaMalhaSistEquacoesBlocoMacro(NALHS_BM, NEQ_BM) 

 1000 format(20a4)
 3000 format(/&
     ' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,&
     ' execution code  . . . . . . . . . . . . . . (exec ) = ',i10//5x,&
     '    eq. 0, data check                                   ',   /5x,&
     '    eq. 1, execution                                    ',  //5x,&
     ' input data print code . . . . . . . . . . . (iprtin) = ',i10//5x,&
     '    eq. 0, print nodal and element input data           ',   /5x,&
     '    eq. 1, do not print nodal and element input data    ',   /5x, &
     ' number of space dimensions  . . . . . . . . (nsd   ) = ',i10/5x, &
     ' number of element nodes     . . . . . . . . (nen   ) = ',i10/5x, &
     ' dimensao do modelo        . . . . . . . (dimModelo ) = ',a)
 4000 format(5x,&
     ' number of nodal points  . . . . . . . . . . (numnpP) = ',i10//5x,&
     ' number of Lados         . . . . . . . . . (numLados) = ',i10//5x,&
     ' number of nodal points  . . . . . . . . . . (numnpD) = ',i10//5x,&
     ' number of nodal degrees-of-freedom  . . . . (ndofc ) = ',i10//5x,&
     ' number of nodal degrees-of-freedom  . . . . (ndofV ) = ',i10//5x,&
     ' number of nodal degrees-of-freedom  . . . . (ndofD ) = ',i10//5x,&
     ' number of load vectors  . . . . . . . . . . (nlvectP) = ',i10//5x,&
     ' number of load vectors  . . . . . . . . . . (nlvectV) = ',i10//5x,&
     ' number of load vectors  . . . . . . . . . . (nlvectD) = ',i10//5x)
 5000 FORMAT(5X,  &
     &' MESH DATA FOR RESERVOIR AND OVERBUDEN DOMAINS:        '//5X,     &
     &'    ELEMENTS IN X-DIRECTION GLOBAL DOMAIN. . (NELX ) = ',I10//5X, &
     &'    ELEMENTS IN Y-DIRECTION GLOBAL DOMAIN. . (NELY ) = ',I10//5X, &
     &'    ELEMENTS IN Z-DIRECTION GLOBAL DOMAIN. . (NELZ ) = ',I10//5X)

!
      inquire(file="passosPressaoBlocoMacro.dat", exist=exist)
      if (exist) then
          iinPassosPressaoBlocoMacro = abrirArquivo("passosPressaoBlocoMacro.dat", 19, 'old' ) 
          read(iinPassosPressaoBlocoMacro, '(i10)'), numPassos          
          write(*,*) "Gravar pressão para", numPassos, "passos"
          allocate(passosTempoParaGravarSaida(numPassos))
          DO I = 1, numPassos 
              read(iinPassosPressaoBlocoMacro,'(i20)') passosTempoParaGravarSaida(I)
              write(*,*) passosTempoParaGravarSaida(I)
          ENDDO  
          close(iinPassosPressaoBlocoMacro)          
      else
           numPassos = 0
      endif     
           
      end subroutine preprocessadorBlocoMacro
           
!
!**** new **********************************************************************
!
     subroutine alocarMemoriaBlocoMacro()

     use mGlobaisArranjos,  only: mat_BM, grav_BM, bf_BM, phi_n, phi_n0
     use mGlobaisArranjos,  only: Knp_BM
     use mGlobaisEscalares, only: ndofD, nlvectD
     use mMalha,            only: nsd_BM, numel_BM, numnp_BM, numLados_BM, nen_BM, numLadosElem_BM, x_BM, xc_BM
     use mMalha,            only: listaDosElemsPorNo_BM, conecNodaisElem_BM
     use mAlgMatricial,     only: id_BM, lm_BM
     use mParametros,       only: phi_BM, K_absRnd

     implicit none

!malha
     allocate(x_BM(nsd_BM,numnp_BM));                      x_BM = 0.0d0
     allocate(xc_BM(nsd_BM,numel_BM));                     xc_BM= 0.0d0
     allocate(conecNodaisElem_BM(nen_BM,numel_BM));        conecNodaisElem_BM=0
     allocate(listaDosElemsPorNo_BM(nen_BM,numnp_BM));     listaDosElemsPorNo_BM=0
!       allocate(listaDosElemsPorFace(numLadosElem,numLados)); listaDosElemsPorFace=0

!material
     allocate(mat_BM(numel_BM)); mat_BM=0.d0

!gravidade
     allocate(grav_BM(3)); grav_BM=0.d0

!!!!!!! SHALE GAS

     allocate(flux_BM (2*ndof_BM, numnp_BM)); flux_BM =0.d0
     allocate(solucao_BM (ndof_BM, numnp_BM)); solucao_BM =0.d0
     allocate(phi_n(numel_BM));        phi_n= 0.0d0
     allocate(Knp_BM(numnp_BM));        Knp_BM= 0.0d0
     allocate(K_absRnd(numel_BM));        K_absRnd= 0.0d0
     allocate(phi_n0(numel_BM));        !phi_n= phi_BM
     allocate(solucaoNaoLinearAnt_BM(ndof_BM, numnp_BM)); solucaoNaoLinearAnt_BM=0.d0
     allocate(solucaoTmpAnt_BM(ndof_BM, numnp_BM)); solucaoTmpAnt_BM=0.d0
     allocate(BF_BM(nesd_BM, numnp_BM));  BF_BM=0.d0

     if (nlvect_BM.ne.0)  then
        allocate(f_BM(ndof_BM,numnp_BM))
        f_BM = 0.0d0
     endif

     allocate(id_BM(ndof_BM,numnp_BM));      id_BM = 0
     allocate(lm_BM(ndof_BM,nen_BM,numel_BM));   lm_BM =0

     end subroutine alocarMemoriaBlocoMacro
!
!**** new **********************************************************************
!
      subroutine topologiaMalhaSistEquacoesBlocoMacro(NALHS_BM, NEQ_BM)
      use mGlobaisEscalares
      use mGlobaisArranjos
!
      use mAlgMatricial,   only: lm_BM, id_BM, idiag_BM, ned_BM
      use mAlgMatricial,   only: ALHS_BM, BRHS_BM, DLHS_BM
!
      use mMalha,          only: numel_BM,numnp_BM,nsd_BM,nen_BM
      use mMalha,          only: x_BM, xc_BM, conecNodaisElem_BM
      use mMalha,          only: listaDosElemsPorNo_BM
!
      use mMalha,          only: genel, formlm, genfl, criarListaVizinhos, local
      use mLeituraEscrita, only: prntel, printD
      use mAlgMatricial,   only: colht, diag

!  
     implicit none
!
!.... program to set storage and call tasks for  
!        an one dimensional conveccion-diffusion problem   
!             in non  dimensional form     
!             displacement formulation
!
!
     integer, intent(in) ::NALHS_BM, NEQ_BM
     integer :: i, m, n
     real*8, allocatable :: xl(:,:) 
!
!.... set element parameters
!
     allocate(npar_BM(numParElem))
     read(iin,'(i10,14i10)') (npar_BM(i),i=1,numParElem)
     ntype   = npar_BM( 1)
     numat_BM = npar_BM( 2)
     nen_BM   = npar_BM( 3)
     nicode_BM  = npar_BM( 4)
     if (nsd_BM==2) nrowsh_BM = 3
     if (nsd_BM==3) nrowsh_BM = 4
     ned_BM   = 1
     NNP = 0

     NENP_BM   = nen_BM
     IHARD_BM  = 0
     NDIMC_BM  = 7
     NEE_BM    = NEN_BM*NED_BM
     NEESQ_BM  = NEE_BM*NEE_BM

     tmGeo=0.0
     tmVel=0.0

     if (nicode_BM.eq.0) nicode_BM=nen_BM
     npint_BM   = nicode_BM 
!
!....... set memory pointers
! 
     allocate(c_BM(7,numat_BM)); c_BM=0.d0
!
!
!.... input element data ('input___')

     write(iecho,1000) ntype,numel_BM,numat_BM,nen_BM,npint_BM
!
!      read material properties
!
     do 400 n=1,numat_BM
     if (mod(n,50).eq.1) write(iecho,4000) numat_BM
     read (iin,  5000) m,(c_BM(i,m),i=1,3)
     write(iecho,6000) m,(c_BM(i,m),i=1,3)
 400 continue
!
!     constant body forces
!
     read  (iin,  7000) (grav_BM(i),i=1,3)
     write (iecho,8000) (grav_BM(i),i=1,3)
!
!      READ MATERIAL PROPERTIES
! 
     allocate(rho_BM(numat_BM))
     allocate(th_BM (numat_BM))
!
     CALL PLANMX(rho_BM,th_BM,C_BM,NUMAT_BM,iecho)
!
!      GENERATION OF NODAL BODY FORCES
!
      READ(iin,7500) IHARD_BM
      WRITE(IECHO,8500) IHARD_BM
      IF(IHARD_BM.GT.0) THEN
           IF(IHARD_BM.EQ.1) THEN
      CALL CLEAR(BF_BM,NESD_BM*NUMNP_BM)
      CALL GENFL(BF_BM,NESD_BM,iin)
      ELSE 
      CALL GHARD(BF_BM,X_BM,NESD_BM,NUMNP_BM,IHARD_BM)
      END IF
!        CALL PRINTD(' N O D A L  B O D Y  F O R C E S       ', &
!      &              BF,NESD,NUMNP,IECHO,NUSTEP,TEMPO)

      END IF
!
!     DATA FOR THE TRANSIENT PROBLEM
!
      READ (iin,11000) DTEMPO_BM,NITER_BM
      WRITE(IECHO,12000) DTEMPO_BM,NITER_BM

      WRITE(*,*) 'DTEMPO,NITER XXX', DTEMPO_BM,NITER_BM
!      call adimTempo(DTEMPO) colocar ese numero em um modulo

!    generation of conectivities

     call genel(conecNodaisElem_BM,mat_BM,nen_BM,iin)
!      call criarListaVizinhos(nen_F,numnp_F,numel_F,conecNodaisElem_F,listaDosElemsPorNo_F)

!#ifdef debug
!      if(iprtin.eq.0) then
!         call prntel(mat_BM,conecNodaisElem_BM,nen_BM,numel_BM, 1, iecho)
!         call prntel(mat_BM,conecNodaisElem_BM,nen_BM,numel_BM, 1)
!      end if
!#endif
!
!     generation of lm array
!
     call formlm(id_BM,conecNodaisElem_BM,lm_BM,ned_BM,ned_BM,nen_BM,numel_BM)
!
!     modification of idiag array
!
     call colht(idiag_BM,lm_BM,ned_BM,nen_BM,numel_BM,neq_BM)
!
     call diag(idiag_BM,neq_BM,nalhs_BM)

     allocate(ALHS_BM(NALHS_BM))
     allocate(DLHS_BM(NALHS_BM))
     allocate(BRHS_BM(NEQ_BM))

!     calculo de xc (centros dos elementos)
     allocate(xl(nsd_BM,nen_BM)); xl=0.0
     do i=1,numel_BM
        call local(conecNodaisElem_BM(1,i),x_BM,xl,nen_BM,nsd_BM,nsd_BM)
        xc_BM(1,i) = sum(xl(1,1:nen_BM))/nen_BM
        xc_BM(2,i) = sum(xl(2,1:nen_BM))/nen_BM
        if(nsd_BM==3) xc_BM(3,i) = sum(xl(3,1:nen_BM))/nen_BM
     end do
!
   
     return
!
1000 format(//,&
    ' two/three-n o d e    e l e m e n t s ',//,5x,&
    ' element type number . . . . . . . . . . . (ntype ) = ',i10,//5x,&
    ' number of elements  . . . . . . . . . . . (numel ) = ',i10,//5x,&
    ' number of element material sets . . . . . (numat ) = ',i10,//5x,&
    ' number of element nodes . . . . . . . . . (nen   ) = ',i10,//5x,&
    ' number of integration points. . . . . . . (npint  ) = ',i10)
4000  format(///,&
    ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
    ' number of material sets . . . . . . . . . . (numat ) = ',i10///,2x,'set',4x,'Kx ',&
      10x,'Ky',10x,'Kz')
5000  format(i10,5x,5f10.0)
6000  format(2x,i3,1x,5(1x,1pe11.4))
7000 format(8f10.0)
7500  FORMAT(I5)
8000 format(///,&
    ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
    ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
    ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
    ' exemplo 3............................ = ',      1pe15.8,//)
8500  FORMAT(5X, &
    &' NODAL BODY FORCE OPTION . . . . . . . . . . (IHARD ) = ',I5//5X, &
    &'    EQ.0, NO NODAL BODY FORCES . . . . . . .            ',   /5X, &
    &'    EQ.1, NODAL GENERATION                              ',   /5X, &
    &'    GE.2, HARD WIRED PROBLEM NUMBER . . .  .            ',   /5X)
9000 format('1',a///&
    ' e q u a t i o n    s y s t e m    d a t a              ',  //5x,&
    ' number of equations . . . . . . . . . . . . (neq    ) = ',i8//5x,&
    ' number of terms in left-hand-side matrix  . (nalhs  ) = ',i12//5x,&
    ' mean half bandwidth . . . . . . . . . . . . (meanbw ) = ',i8//5x,&
    ' memoria necessaria para a matriz do sistema  (Mbytes)  = ',e10.2)
11000 FORMAT(F10.6,I10)
12000 FORMAT(5X,  &
    &' STEP TIME . . . . . . . . . . . . . . . . .(DTEMPO)=',1PE15.8 //5X, &
    &'NUMBER OF STEPS. . . . . . . . . . . . . . . (NITER ) = ',I5)
!
     end  subroutine topologiaMalhaSistEquacoesBlocoMacro
     
! **************************************************************************
      end module mBlocoMacro
