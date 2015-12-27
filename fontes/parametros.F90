!     ************************************************************
!     SETEMBRO/2014
!     Patricia de Araujo Pereira
!     LNCC/MCT
!  
!     ************************************************************
!


    module mParametros
       
       implicit none


!-------------------------------------------------------------------------------------
!      Parametros para fratura e bloco
!-------------------------------------------------------------------------------------
       REAL*8 :: nAv          ! numero de avogrado 
       REAL*8 :: R_           ! constante universal dos gases
       REAL*8 :: aV, bV;      ! atraction and repulsion forces (?)                 
       REAL*8 :: Tc           ! temperatura do gás no ponto crítico
       REAL*8 :: Pc           ! pressao do gás metano no ponto crítico
       REAL*8 :: T            ! temperatura do gás 
       REAL*8 :: Tr           ! temperatura de referencia
       REAL*8 :: constMu      ! +++ viscodidade; deveria variar com a pressão
       REAL*8 :: M_m          ! massa molar do metano kg/mol 

       REAL*8  :: p_Reservatorio
       REAL*8  :: p_Ref
       REAL*8  :: p_Poco
       REAL*8  :: condContAdim

       REAL*8  :: gasTotalKg
       REAL*8  :: gasRecuperavelKg
       REAL*8  :: gasProduzidoKg
       
       
!-------------------------------------------------------------------------------------
!      Parametros para bloco
!-------------------------------------------------------------------------------------
       REAL*8 :: rhoL;        ! massa específica da fase líquida
       REAL*8 :: Hch;         ! constante de Henry
       REAL*8 :: v2_inf;      ! partial molar volume of CH4 in the liquid phase 
       REAL*8 :: PL_sat;      ! saturation pressure of the water 
       REAL*8 :: D_CH4_H2O;   ! coeficiente de difusao do metano na água
       REAL*8 :: gammaMax;    ! numero de sitios disponiveis para adsorção de moléculas de metano
       REAL*8 :: A_sf;        ! área superficial   |∂Y^W_fs|/|Y^W| 
       REAL*8 :: P_L;         ! reference value of the Langmuir pressure which provides 
                              ! one-half of the total adsorption capacity.                                                        
             
       REAL*8 :: rhoK   ! massa específica do querogênio, usado par calculo do TOC
       REAL*8 :: rhoI   ! massa específica da fase inorgânica, usado par calculo do TOC 
       
       REAL*8 :: Sw, phi_M, phi_N, TOC, phi_K;
       REAL*8 :: kEff, K_, D;
       REAL*8 :: K_abs, K_re;

       REAL*8 :: tamBloco; 
       REAL*8 :: alturaBloco; ! usada para calcular o volume do bloco e a quantidade de gás total
       REAL*8  :: gasTotalKg_B
       REAL*8  :: gasRecuperavelKg_B

       ! parametros especificos para cada uma das equaçoes de estado 
       ! estamos usando van der Walls
!        REAL*8 :: epsilon_
!        REAL*8 :: sigma
!        REAL*8 :: gamma_
!        REAL*8 :: omega
!        REAL*8 :: w
!        REAL*8 :: alfa      

!-------------------------------------------------------------------------------------
!      Parametros para o bloco macroscopico
!-------------------------------------------------------------------------------------

       REAL*8  :: phi_BM
!        REAL*8, allocatable  :: phi_n(:) 
       REAL*8  :: fraVol_BM 
       REAL*8  :: constK_BM
       REAL*8  :: k_s	! Bulk modulus da matriz shale, Pa (Diego, set/2015)
       REAL*8  :: beta_r	! Compressibilidade total da rocha (Diego, set/2015)
       REAL*8  :: Kbulk		! Bulk modulus do bloco macro (Diego, set/2015)
       REAL*8  :: alpha_r	! Coeficiente de Biot do bloco macro (Diego, set/2015)
    
       REAL*8  :: tamBlocoMacro      ! altura do bloco macro    (F_y = BM_y) 
       REAL*8  :: widthBlocoMacro    ! espessura do bloco macro (BM_x)
       REAL*8  :: dimZ_FB            ! dimensão z da fratura hidraulica e do bloco macro (F_z, BM_z)
       REAL*8  :: areaContatoFratPoco ! area de contato entre a fratura e o poço, usada para calcular a produção de gas
       REAL*8  :: areaContatoBlocoMacroFratura ! area de contato entre o bloco macro e a fratura hidráulica, usada para calcular a produção de gás
       REAL*8  :: gasTotalKg_BM
       REAL*8  :: gasRecuperavelKg_BM

    contains


!     ************************************************************
    subroutine inicializarParametros()

       implicit none
 
       REAL*8  kB    ! constante de Boltzmann
       REAL*8  Swr, Sgr, Se
       
       kB        = 1.381d-23;        
       nAv       = 6.0221415d23
       R_        = nAv*kB    ! 8.31d0   (Pa m3)/(mol ºK)
       aV        = 0.225;    !  Pa.m^6/mol^2 -- wiki for methane
       bV        = 4.28d-5;  ! m^3/mol -- wiki for methane
       Tc        = 8*aV/27/bV/R_;    ! DUNG       ! antigo Tc        = 190.6d0
       Pc        = aV/27/(bV*bV);   ! DUNG       ! antigo Pc        = 4.599d6
!        T         = 320.d0  ! DUNG                ! antigo T         = 323.15d0  
       Tr        = T/Tc;   
       constMu   = 1.2d-5; ! Pa.s       
       M_m       = 16.01D-3   ! kg/mol 

!        p_Reservatorio = 6.32d7;           ! Pa -> Barnett
!        p_Reservatorio = 8.72d7;		  ! Pa -> Eagle Ford
       p_Ref          = p_Reservatorio;   ! Pa 
!        p_Poco         = 5.d5;             ! Pa      
       condContAdim   = p_Poco!/p_Ref;     
       
       dimZ_FB       = 1.d0 !m
                    
     end subroutine inicializarParametros
!     ************************************************************
    
    
    
!     ************************************************************
    subroutine inicializarParametrosBloco()

       implicit none
 
       REAL*8  :: Swr, Sgr, Se              
       REAL*8  :: y, Tw, D_l, D_l0
       INTEGER :: I;
       
       rhoL      = 5.556d4;   ! = 1.d6/18       mol/m^3
       Hch       = 5.d9       ! Pa, DUNG        antigo Hch        = 4.185d9;
       PL_sat    = 7.d5       ! Pa, DUNG        antigo PL_sat     = 1.23d4 ! Pa
       v2_inf    = 3.4501d-5  ! DUNG            antigo v2_inf     = 4.d-5  ! m^3/mol 
       D_CH4_H2O = 1.d-9      ! m^2/s;   ! eu estava usando 1.8d-9, mas fiz isso para ficar igual ao codigo Dung
       gammaMax  = 10.d18/nAv 

       rhoK      = 1.2; !  DUNG                 antigo rhoK      = 1.3; 
       rhoI      = 2.7; !  DUNG                 ! antigo rhoI      = 2.36;
       ! a unidade não importa pq faremos rhoK/rhoI

! !        !___________ PARA TESTE
! !        Swr=0.2
! !        Sgr=0.2;
! !        Sw =0.0
! !        K_abs = 1.d-19; ! 1 microD = 1.e-18 m^2 %Daniel J. Soeder, Scty of Ptrleum engineers, 1988
! !        phi_M = 0.1
! !        phi_N = 0.05
! !        
! !        DO I=1,11
! !         Se   = (Sw-Swr)/(1-Swr-Sgr);
! !         K_re = (1-Se)*(1-Se)*(1-Se**2);
! !         K_ = K_re*K_abs/constMu*2*phi_M/(3-phi_M);  ! self-consitent
! !         write(*,*) Sw, K_re, K_
! !         Sw= Sw + 0.1  
! !        END DO
! !        stop       
! !        !___________
          
       ! 1. se trocar o Sw, tem que trocar o y usado para calcular D=D(Sw), de acordo com a tabela 
       ! 2. Se Sw = 0.0, entao Swr = 0.0, pq Sw >= Swr
       
        Sw        = 0.3d0
        Swr       = 0.2;
        
!          Sw        = 0.0d0      
!          Swr       = 0.0;
               
        Sgr       = 0.2;
        Se        = (Sw-Swr)/(1-Swr-Sgr);
       
        phi_M = 0.1
        phi_N = 0.05
        TOC   = 0.05
       
       
       phi_K = (1-phi_M)/(1+(1-phi_N)*rhoK*(1-TOC)/rhoI/TOC);  
       
       ! Permeability
       K_abs = 1.d-19; ! 1 microD = 1.e-18 m^2 %Daniel J. Soeder, Scty of Ptrleum engineers, 1988
!        K_abs = 1.d-18; ! 1 microD = 1.e-18 m^2 %Daniel J. Soeder, Scty of Ptrleum engineers, 1988
       
       K_re = (1-Se)*(1-Se)*(1-Se**2);
       
       write(*,*) "Sw, K relativo", Sw, K_re;
       K_ = K_re*K_abs/constMu*2*phi_M/(3-phi_M);  ! self-consitent
       
       write(*,*) "Sw, K absoluto", Sw, K_;
                    
       ! Diffusion  
       D=D_CH4_H2O*2*phi_M/(3-phi_M);              ! self-consitent ! old way
       
       ! Diffusion D=D(S_w)
       !%y  = fsolve(@(y) (Sw^2)^y+(1-Sw)^y-1,0.5 );
       ! |-------------------|
       ! |  Sw   |       y   |
       ! |-------|-----------|
       ! | 0.0   |   0.5     |
       ! | 0.05  |   0.5871  |     
       ! | 0.1   |   0.6048  |
       ! | 0.15  |   0.6186  |
       ! | 0.2   |   0.6308  |
       ! | 0.25  |   0.6420  |
       ! | 0.3   |   0.6527  |
       ! | 0.4   |   0.6734  |
       ! | 0.5   |   0.6942  |
       ! | 0.6   |   0.7161  | 
       ! | 0.7   |   0.7401  |
       ! | 0.8   |   0.7684  |
       ! | 0.9   |   0.8062  |
       ! | 1.0   |   0.5     |
       ! |-------------------|
       y    = 0.6527;   ! +++ y varia em função de Sw!!!!
!         y    = 0.5;      ! +++ y varia em função de Sw!!!!
       Tw   = Sw**(2*y+1); ! Tortuosity
       D_l0 = 2.4e-9;
       D_l  = D_l0*Sw*Tw;
       D    = D_l*2*phi_M/(3-phi_M);
            
       write(*,*) "Coeficiente de difusão nos blocos D(Sw)", D
       
       ! para resolver o problema via formulação concentração 
!         p_Ref         =  4.36579470d3      ! mol/m^3    
!         condCont     =  1.89361775d2;       ! mol/m^3    

       tamBloco    = 10.d0;  
       alturaBloco = tamBlocoMacro;
       
       
     end subroutine inicializarParametrosBloco
!     ************************************************************


!     ************************************************************
      SUBROUTINE inicializarParametrosBlocoMacro()
      
      use mGlobaisArranjos,  only: phi_n, phi_n0, phi_range
      use mGlobaisEscalares, only: random_porosity
      
      IMPLICIT NONE
      
      CHARACTER*200 :: linhaAux
      INTEGER       :: I
                
      phi_BM      =  0.2d0        ! adim - porosidade das fraturas naturais
      fraVol_BM   =  1.0d0       ! adim - fracao de volume das fraturas naturais
!       constK_F   =  1.D-15      ! m^2
!       constK_BM   =  1.0D-18    ! m^2
!       constK_BM   =  1.0D-18    ! m^2
!       k_s = 25.0d9	! Pa, Tobiloluwa (Diego, set/2015)
!       Kbulk = (6894.75729)*2.6d6	! Shale gas revolution (Diego, set/2015)
!       Kbulk = 5.28d9	! Pa, Skalle (Diego, set/2015)
      alpha_r = 1.0d0 - Kbulk/k_s
!       alpha_r = 1.0d0
!       alpha_r = 0.001

!       tamBlocoMacro          = 50.d0  ! metros - altura bloco macro (BM_y)
!       widthBlocoMacro        = 10.D0  ! metros - espessura bloco macro (BM_x)
!       areaContatoFratPoco = widthBlocoMacro * dimZ_FB;
      areaContatoFratPoco = 0;
      areaContatoBlocoMacroFratura = tamBlocoMacro * dimZ_FB
      
!       if (random_porosity .eqv. .true.) then
!       call random_number(phi_n0)
!       phi_n0 = (phi_n0*(0.25d0-0.15d0)) + 0.15d0
!       phi_n = phi_n0
!       else
!       phi_n0 = 0.25d0
!       phi_n = phi_n0
!       endif
      
      call random_number(phi_n0)
      phi_n0 = (phi_n0*(phi_range(2)-phi_range(1))) + phi_range(1)
      phi_n = phi_n0
          
      END SUBROUTINE inicializarParametrosBlocoMacro
      


!     ************************************************************
    subroutine calcularPhiI(phiM, phiN, TOC, phiI)
!
!        use mParametros;
        implicit none                    

       INTEGER    I;
       REAL*8     phiM, phiN, phiI, TOC;

       phiI = 1 - phiM - (1-phiM)/(1 + ( ((1-phiN)*rhoK*(1-TOC))/(rhoI*TOC) ))
              
!        write(*,220) TOC, phiI
         
 220  format(F9.5, F9.5)
 

    end subroutine calcularPhiI
!     ************************************************************


   !   
   end module mParametros


