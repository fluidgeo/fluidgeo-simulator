
!!     ************************************************************
!     SETEMBRO/2014
!     Patricia de Araujo Pereira
!     LNCC/MCT
!     Modificado por Diego Volpatto. Nov (2015)
!  
!     ************************************************************
!

    module mCoeficientes
       
       implicit none

!        REAL*8              :: z0,  z1,  z2,  z3,  z4,  z5,  z6 
       REAL*8              :: z0,  z1,  z2,  z3,  z4,  z5,  z6, z7, z8, z9
       REAL*8              :: gm0, gm1, gm2, gm3
       REAL*8              :: g12, g11, g10, g9, g8, g7, g6, g5, g4, g3,  g2,  g1,  g0

       REAL*8              :: gmLinha3, gmLinha2, gmLinha1, gmLinha0
       integer             :: iechoQuantidadeGas
 
       contains


!     ***********************************************************
       subroutine inicializarCoeficientes_FormP()

          use mGlobaisArranjos,  only: reservoir_case
          implicit none      
          
	  if (reservoir_case .eq. "Barnett") then
		z9 = -1.87107151262658931063530292883052060553757276234474d-69
		z8 =  4.61904974863756544871671744183262036246196976705905d-61
		z7 = -4.74966407523363820528916577631702818378758980133890d-53
		z6 =  2.57059856765925816157611341245935475535833353616401d-45
		z5 = -7.37364199704143266168241625354514632553541863947761d-38
		z4 =  8.60397932783036217170980456336014633101296514147358d-31
		z3 =  9.20635847128269864669945643703931564721528582181228d-25
		z2 =  3.68047931147880489270481476821021104296660013226183d-16
		z1 = -1.40612027122484796410775253604902301773904582660180d-08
		z0 =  1.00002521342627681555370600108290091156959533691406d+00
          else if (reservoir_case .eq. "Marcellus") then
                z9 = -1.41198588449330946930195173913772931514821145181013d-70
		z8 =  4.08952559506348956786000715311376733535803178795308d-62
		z7 = -4.51117795072239156925468754041945827666398072194733d-54
		z6 =  2.04063492760344119267079403311173657857483418100168d-46
		z5 =  7.97876080525611670962674089732074998323977842827591d-40
		z4 = -4.31083994696858647558194747846842879124848660025549d-31
		z3 =  1.18124319921963000907932160166645689043561760623650d-23
		z2 =  2.78838592506377415762739646902775903815341340674649d-16
		z1 = -1.19418631483071175177721137835272491312110787475831d-08
		z0 =  9.99877401853861536018541755765909329056739807128906d-01
          else if (reservoir_case .eq. "Default") then
		z9 = -3.68431902256159053689437065251333093037795407893241d-69
		z8 =  9.14946624686535986182556642360366660211415974666699d-61
		z7 = -9.44338011856949795344994323397196930574391246777314d-53
		z6 =  5.13687001799859286834788551460446931005592377556025d-45
		z5 = -1.50312026743859811014673427238691130569623659438624d-37
		z4 =  1.95593931572846383927528426545236377126803207845334d-30
		z3 = -2.61495319029415052586608495843491372166791807335120d-24
		z2 =  4.04285221837828608518775198752794922853351047191753d-16
		z1 = -1.71143833232033277489744777085864391175107357412344d-08
		z0 =  1.00002845080220614804034084954764693975448608398438d+00
          endif

           gm3 =  4.69685620646729666887173708310525515307394743146910d-24
           gm2 =  1.11154469983911601384394559546778573976707712479828d-16
           gm1 = -1.56446443490768436927189653567679683554558778268984d-08
           gm0 =  1.00001782339149025702340622956398874521255493164062d+00

           gmLinha3 = 5.80012011155713811161051814852137957791755881683127d-31
           gmLinha2 = 5.29699871739355908979459025069218808404007231271524d-24
           gmLinha1 = 2.60672926117014932258083850664859388895914962380046d-16
           gmLinha0 = -1.5688610161285264504953311059769205382252721392433d-08

              g12 = -7.63045236378476089810643306380680483955051598244986d-87
              g11 =  1.11067443561618653831513234064030962976641511221974d-78
              g10 = -5.71545284175409686986785258958874441356429583564553d-71
              g9  =  5.00709244753032033749739285163102243983589921091405d-64
              g8  =  8.07757243336225834020277915547222824898422890564172d-56
              g7  = -4.56359653437325984822649836578985781932518055767335d-48
              g6  =  1.26990222746843495546527654976406747647066782595377d-40
              g5  = -2.20335133623880615127261102846821081317358359421996d-33
              g4  =  2.53486439585973867288620212745978302523810041234362d-26
              g3  = -1.96776286454299110177072541113156410960386028289043d-19
              g2  =  1.03708649703153532170629685189980462149163253915418d-12
              g1  = -3.81578248014776313365220214435247214623814215883613d-06
              g0  =  1.08535647409796247586655226768925786018371582031250d+01

       end subroutine inicializarCoeficientes_FormP
       
!     ************************************************************
      SUBROUTINE calcularZ_P(p, Z)
      
         IMPLICIT NONE     
         REAL*8   :: p, Z
         
!          Z = 1; 
!          return;

!          Z     =  z5*(p**5) + z4*(P**4) + z3 *p*p*p + z2 *p*p + z1 *p + z0
         Z = z9*(p**9) + z8*(p**8) + z7*(p**7) + z6*(p**6) + z5*(p**5) + & 
         & z4*(p**4) + z3*(p**3) + z2*(p**2)+ z1*p + z0
           
      END SUBROUTINE calcularZ_P
!     ************************************************************


!     ************************************************************
      SUBROUTINE calcularZ_C(c, Z)
      
         use mParametros        
         IMPLICIT NONE
      
         REAL*8   :: c, Z
         
         Z  = 1./(1-bV*c) - aV*c/R_/T;
           
      END SUBROUTINE calcularZ_C
!     ************************************************************


!     ************************************************************
       SUBROUTINE calcularZLinha(c,ZLinha)
       
          use mParametros,   only: R_, T, aV, bV       
          REAL*8   :: c, ZLinha
          
          ZLinha     =  bV/((1-bV*c)*(1-bV*c)) - aV/(R_*T)
                   
!         write(*,*) 'ZLinha',  ZLinha

      END SUBROUTINE calcularZLinha
!     ************************************************************


!     ************************************************************
      SUBROUTINE calcularQuantidadeGasBloco()
      
         use mParametros,   only: p_Ref, p_Poco, M_m   
         use mParametros,   only: gasTotalKg_B, gasRecuperavelKg_B
         use mParametros,   only: tamBloco, alturaBloco, dimZ_FB

         IMPLICIT NONE
         
         REAL*8  :: gasRetiravel, volBloco         
         REAL*8  :: p_Inicial, gasLivreI, gasDissolvidoI, gasAdsorvidoI, gasInicial
         REAL*8  :: p_Final,   gasLivreF, gasDissolvidoF, gasAdsorvidoF, gasFinal
!        ainda não temos bloco
!        não fazemos nada aqui por enquanto  
         gasTotalKg_B = 0
         gasRecuperavelKg_B = 0 
         return 
         volBloco = tamBloco*alturaBloco*dimZ_FB        

         p_Inicial = p_Ref
         p_Final   = p_Poco
         
         call calcularQuantidadeGas(p_Inicial, gasLivreI, gasDissolvidoI, gasAdsorvidoI)
         call calcularQuantidadeGas(p_Final,   gasLivreF, gasDissolvidoF, gasAdsorvidoF)
      
         gasInicial = (gasLivreI + gasDissolvidoI + gasAdsorvidoI)*M_m
         gasFinal   = (gasLivreF + gasDissolvidoF + gasAdsorvidoF)*M_m

         gasRetiravel        = (gasInicial - gasFinal)
         
         gasTotalKg_B        = gasInicial  *volBloco 
         gasRecuperavelKg_B  = gasRetiravel*volBloco
         
         write(iechoQuantidadeGas,*)    "-----------------------  B L O C O  -------------------------"
         write(iechoQuantidadeGas,300)  "     pressao(Pa)      gasLivre         gasDissolvido    gasAdsorvido  &
       &                                      gasTotal(kg/m^3)  gasTotal(kg)"

         write(iechoQuantidadeGas,200)  p_Inicial, gasLivreI*M_m, gasDissolvidoI*M_m, gasAdsorvidoI*M_m, &
       &                                gasInicial, gasInicial*volBloco
         write(iechoQuantidadeGas,200)  p_Final  , gasLivreF*M_m, gasDissolvidoF*M_m, gasAdsorvidoF*M_m, &
       &                                gasFinal,  gasFinal*volBloco

!         write(iechoQuantidadeGas,210)  "A quantidade total de gas no bloco, com pressao", p_inicial, "pa", gasInicial, "kg", gasInicial*volBloco, "kg"         
!         write(iechoQuantidadeGas,210)  "A quantidade total de gas no bloco, com pressao", p_inicial, "(Pa) e'", gasInicial, "(kg/m^3) ou", gasInicial*volBloco, "(kg)"
         write(iechoQuantidadeGas,210)  "A quantidade de gas recuperavel, reduzindo a pressao para", & 
       &                                p_Final, "(Pa) e'", gasRetiravel, "(kg/m^3) ou", gasRetiravel*volBloco, "(kg)"
         write(iechoQuantidadeGas,220)  "O volume do bloco e'", tamBloco, " X", alturaBloco, " X", dimZ_FB, " =", volBloco, "(m^3)"
         write(iechoQuantidadeGas,*)    ""
         
 300  FORMAT(A)                          
 200  FORMAT(6(1PE15.5,2X))         
 210  FORMAT(3(A,1PE15.5),A)               
 220  FORMAT(4(A,1PE10.3),A)               
 
      END SUBROUTINE calcularQuantidadeGasBloco
!     ************************************************************
      
      
!     ************************************************************
      SUBROUTINE calcularFluxoGasLivre(p1Star, p2Star, gradP, J_gasLivre)
      
         use mParametros,   only: K_, R_, T, phi_M, Sw, phi_K, phi_N, rhoL, p_Ref, tamBloco
         IMPLICIT NONE                 
         
         REAL*8  ::  p1Star, p2Star, gradP, J_gasLivre
         REAL*8  ::  p, Z
         
         p = (p2Star+p1Star)/2
         !p = p1
         
!          call calcularZ_P(p*p_Ref,Z)
         call calcularZ_P(p,Z)
!          J_gasLivre = p/(Z*R_*T) *p_Ref * K_ * gradP * p_Ref/tamBloco         
         J_gasLivre = p/(Z*R_*T) * K_ * gradP/tamBloco         

      END SUBROUTINE calcularFluxoGasLivre
      
!     ************************************************************
      SUBROUTINE calcularQuantidadeGas(p, gasLivre, gasDissolvido, gasAdsorvido)
      
         use mParametros,   only: Hch, v2_inf, PL_sat, R_, T, phi_M, Sw, phi_K, phi_N, rhoL
         IMPLICIT NONE                 
      
         REAL*8   :: p, gasLivre, gasDissolvido, gasAdsorvido
         REAL*8   :: Z, Gamma, G, H, cB

         call calcularZ_P(p,Z)
         Gamma =  gm3*p*p*p + gm2*p*p + gm1*p + gm0
         G     =  g12*(p**12) + g11*(p**11) + g10*(p**10) + g9*(p**9) + g8*(p**8) +       &
       &          g7*(p**7) + g6*(p**6) +         &
       &          g5*(p**5) + g4*(p**4) + g3*(p**3) + g2*(p**2) + g1*p    + g0      
         H     =  Hch*exp(v2_inf*(p-PL_sat)/(R_*T));
         cB    =  p/(Z*R_*T)             
         
         gasLivre      =  cB * (1-Sw) * phi_M
         
         gasDissolvido = gamma*rhoL/H * p * Sw * phi_M
         
         gasAdsorvido  = G * cB * phi_N * phi_K

!          write(*,*) 'p             ', p
!          write(*,*) 'gasLivre      ',gasLivre
!          write(*,*) 'gasDissolvido ',gasDissolvido
!          write(*,*) 'gasAdsorvido  ',gasAdsorvido

           
      END SUBROUTINE calcularQuantidadeGas
!     ************************************************************

!
!**** new **********************************************************************
!
      
      SUBROUTINE calcularQuantidadeGasBlocoMacro()
      
         use mParametros,   only: R_, T, M_m, p_Poco, phi_BM, M_m, p_Ref;   
         use mParametros,   only: tamBlocoMacro, widthBlocoMacro, dimZ_FB, gasTotalKg_BM, gasRecuperavelKg_BM    
         use mMalha,        only: numel_bm
         use mGlobaisArranjos, only: phi_n0
         
         implicit none;

         REAL*8 :: Z, gas_pInicial, gas_pFinal, gasRecuperavel, volBlocoMacro;
         REAL*8 :: pInicial, pFinal;
                 
         volBlocoMacro = tamBlocoMacro*widthBlocoMacro*dimZ_FB        
         
         pInicial = p_Ref;
         pFinal   = p_Poco;
         
         call calcularZ_P(pInicial, Z);         
         gas_pInicial   = pInicial /(Z*R_*T) * (sum(phi_n0)/numel_bm) * M_m       ! kg/m^3
         
         call calcularZ_P(pFinal, Z);
         gas_pFinal     = pFinal/(Z*R_*T) * (sum(phi_n0)/numel_bm) * M_m          ! kg/m^3
          
         gasRecuperavel      = gas_pInicial - gas_pFinal          ! kg/m^3  
                  
         gasTotalKg_BM       = gas_pInicial   * volBlocoMacro     ! kg
         gasRecuperavelKg_BM = gasRecuperavel * volBlocoMacro     ! kg
         
         write(iechoQuantidadeGas,*)   "---------------------- BLOCO MACRO ------------------------"
         
         write(iechoQuantidadeGas,300)  "    pressao(Pa)    gasTotal(kg/m^3)  gasTotal(kg)"
         write(iechoQuantidadeGas,200)  pInicial,  gas_pInicial, gas_pInicial*volBlocoMacro
         write(iechoQuantidadeGas,200)  pFinal  ,  gas_pFinal,   gas_pFinal  *volBlocoMacro
         
         write(iechoQuantidadeGas,210)  "A quantidade total de gás no bloco macro, com pressao", pInicial, & 
       &                                "(Pa) é", gas_pInicial, "(kg/m^3) ou", gasTotalKg_BM, "(kg)"
         write(iechoQuantidadeGas,210)  "A quantidade de gás recuperável, reduzindo a pressão para", pFinal, &
       &                                "(Pa) é", gasRecuperavel, "(kg/m^3) ou", gasRecuperavelKg_BM, "(kg)"
         write(iechoQuantidadeGas,220)  "O volume do bloco macro é", widthBlocoMacro, " X", tamBlocoMacro, & 
       &                                " X", dimZ_FB, " =", volBlocoMacro, "(m^3)"
         write(iechoQuantidadeGas,*)    ""

 300  FORMAT(A)                          
 200  FORMAT(3(1PE15.5,2X))         
 210  FORMAT(3(A,1PE15.5),A)               
 220  FORMAT(4(A,1PE10.3),A)               
 
      end subroutine calcularQuantidadeGasBlocoMacro



!     ************************************************************
      SUBROUTINE calcularQuantidadeGasReservatorio()

         use mParametros,   only: gasTotalKg, gasRecuperavelKg, p_Ref, p_Poco
         use mParametros,   only: gasTotalKg_B, gasRecuperavelKg_B, gasTotalKg_BM, gasRecuperavelKg_BM
         
      
          call calcularQuantidadeGasBlocoMacro();     
          call calcularQuantidadeGasBloco();     
          
          gasTotalKg        =  2*gasTotalKg_B       + gasTotalKg_BM
          gasRecuperavelKg  =  2*gasRecuperavelKg_B + gasRecuperavelKg_BM                    
          
         write(iechoQuantidadeGas,*)     "------------------   R E S E R V A T O R I O  -------------------------"
         write(iechoQuantidadeGas,210)  "A quantidade total de gas em uma fratura e dois blocos, com pressao", &
      &                                 p_Ref, "(Pa) é", gasTotalKg, "(kg)"
         write(iechoQuantidadeGas,210)  "A quantidade de gás recuperavel, reduzindo a pressao para", &
      &                                 p_Poco, "(Pa) é", gasRecuperavelKg, "(kg)"
         write(iechoQuantidadeGas,*)    ""

 210  FORMAT(2(A,1PE15.5),A)                        

      END SUBROUTINE calcularQuantidadeGasReservatorio
!     ************************************************************
      
end module mCoeficientes

