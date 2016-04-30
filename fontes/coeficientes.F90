
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

          use mGlobaisArranjos,  only: reservoir_case, reservoir_depth

          implicit none      
          
	  if (reservoir_case .eq. "Barnett") then
          if (reservoir_depth .eq. "deep") then
		z9 = -1.46813275864625275340248499807375363330232413388946d-69
		z8 =  3.89286860815022995767815577802742829208759485046757d-61
		z7 = -4.20216264893913085875639833749768734092830317190706d-53
		z6 =  2.34652438552974914994651591014486670859497542832351d-45
		z5 = -6.83328474629057939848682573913470030424436097710665d-38
		z4 =  7.82407936202825200588374484099196994797136786767157d-31
		z3 =  1.57469723499054256160833559235538845010094679741240d-24
		z2 =  3.65108201531390095727445040254060942955032995071996d-16
		z1 = -1.40552212649040271292603956573320778122848651037202d-08
		z0 =  1.00002150782930088190880724141607061028480529785156d+00
    else if (reservoir_depth .eq. "shallow") then
		z9 =  1.05051859508609798315934360560149325516963157643906d-68
		z8 = -1.32693006912533947290680956442010463612415806251938d-60
		z7 =  6.18175691488432948566783690907570356642321568248978d-53
		z6 = -1.09594687947361750662850961963980876318205028866178d-45
		z5 = -4.30114796343324441277816616037572494910330940396453d-40
		z4 = -2.74846951182484554271033057507003686292117172573491d-32
		z3 =  7.25051970710602153749825409482109115638560860091129d-24
		z2 =  3.43512429508569771754599155163368566803903399577494d-16
		z1 = -1.40173709345126285792148873501211658876286492159124d-08
		z0 =  1.00000068863746860436947372363647446036338806152344d+00
           endif
        endif
    if (reservoir_case .eq. "Marcellus") then
          if (reservoir_depth .eq. "shallow") then
		z9 = -7.98895291698732445674490481624196333419077971226648d-70
		z8 =  2.06734839141633471293566515214478334980052870185833d-61
		z7 = -2.20815538212100272945411713233799356777385338483419d-53
		z6 =  1.21865612971548572856744516373821946426709358479981d-45
		z5 = -3.38418744517029589308353965148872685037737752255913d-38
		z4 =  2.78403647000864843211164514481374410926427028326407d-31
		z3 =  3.36622544823788784878876183247780168761056853814598d-24
		z2 =  3.32509948163453924952944528677456822984425684587920d-16
		z1 = -1.20940215101529927208345326409243136733806522897794d-08
		z0 =  1.00000262779538195978545900288736447691917419433594d+00
          else if (reservoir_depth .eq. "deep") then
		z9 = -2.42926517319393920651606665298459161434867679132808d-70
		z8 =  6.96584174374744927048653426883823314263335677413527d-62
		z7 = -7.90995787006734393072263251251053413138923617742665d-54
		z6 =  4.21756086831535086775256556143298098377555103824524d-46
		z5 = -7.40016953399249021667823040632784926032215181166411d-39
		z4 = -2.46929568577390678055712987160325876321312377285599d-31
		z3 =  9.42231835327529024372629441226386885009463920848079d-24
		z2 =  2.95290097079410864772724180594317997860333752145959d-16
		z1 = -1.19919829040115086206268166821135856547897446944262d-08
		z0 =  9.99921144125617722409060661448165774345397949218750d-01
          endif
        endif
        if (reservoir_case .eq. "Default") then
		z9 = 2.44607753144010275912270835506986340579952760811386d-68
		z8 = -3.17995473655437825551062457959541722421961094200843d-60
		z7 = 1.55412845527681961006050437355374362038786050116595d-52
		z6 = -3.14360996436033073913200620508303555086909725473188d-45
		z5 = 1.10534570471081858395271201170316210212011387045732d-38
		z4 =  8.59628983449214775531714125802186192211141107401304d-32
		z3 = 9.80554282268042934512911051817996381727225584150280d-24
		z2 =  3.61252647056522632950428200479880896105125550558240d-16
		z1 = -1.70506264451730330194726358537804511428248588345014d-08
		z0 =  1.00000378745711571149001883895834907889366149902344d+00          
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
      SUBROUTINE calcularMu_P(p, mu)
      
         IMPLICIT NONE     
         REAL*8   :: p, mu

         mu = 1.2d-5;
         

!          mu     =  mu2*(p**2) + mu1*(p) + mu0
           
      END SUBROUTINE calcularMu_P

     ! a pressao de entrada, p, eh p* X p_ref
!     ************************************************************
!      SUBROUTINE calcularK_P(p, phi, K_eff)
      
!         use mParametros,   only: K_abs, K_re, Hch, v2_inf, PL_sat, R_, T, rhoL, D, D_kn
!         use mParametros,   only: M_m, pi, phi_M
            
!         REAL*8   :: p, K_eff, K_eff_Old, Z, Gamma, G, H, HLinha, GammaLinha, GammaDivH_Linha, FLinha, D_Knudsen, ZLinha
!         REAL*8   :: r, fc1, fc2, K_eff1, K_eff2, K_, mu, phi
!          REAL*8   :: fc6, fc7, fc8, fc9, K_eff6, K_eff7, K_eff8, K_eff9;
!         REAL*8   :: fc
              
!         call calcularMu_P(p, mu);
!         K_     = K_re*K_abs/mu*2*phi/(3-phi);  ! self-consitent                                         
     
!         call calcularZ_P(p, Z);
!         Gamma =  gm5*(p**5) + gm4*(p**4) + gm3*(p**3) + gm2*(p**2) + gm1*p  + gm0
!         H     =  Hch*exp(v2_inf*(p-PL_sat)/(R_*T));
!
!         HLinha     =  H*v2_inf/(R_*T);
!         GammaLinha =  gmLinha10*(p**10) + gmLinha9*(p**9) + gmLinha8*(p**8) + &
!      &		       gmLinha7*(p**7) + gmLinha6*(p**6) + gmLinha5*(p**5) + &
!      &                gmLinha4*(p**4) + gmLinha3*(p**3) + gmLinha2*(p**2) + gmLinha1*(p) + gmLinha0 
         
!         GammaDivH_Linha =  (GammaLinha*H - HLinha*Gamma)/(H*H)

!         r  = 1.d-8;
!         call calcular_fc(p, r, fc);

!         K_eff    =  K_ * p/(R_*T*Z) * fc! + D * ( Gamma*rhoL/H + p*rhoL * GammaDivH_Linha );
        
!100  FORMAT(5(2X,1PE15.8))               
!150  FORMAT(6(2X,1PE15.8))               
!200  FORMAT(8(2X,1PE15.8))               
        
!      END SUBROUTINE calcularK_P

!     ************************************************************
!      SUBROUTINE calcular_fc(p, r, fc)
!
!         use mParametros,   only: R_, T, correcaoPerm 
!         use mParametros,   only: M_m, pi
            
!         REAL*8   :: p, r, fc, Z, mu
!         REAL*8   :: lambda, A, B, alpha, alpha0, Kn;

!          lambda = constMu/p * sqrt((pi*R_*T)/(2*M_m));

!         call calcularZ_P(p,Z);
!         call calcularMu_P(p,mu);
!         Z=1;
!         lambda = mu*Z/p * sqrt((pi*R_*T)/(2*M_m));
!         Kn = lambda/r;
         
!         A = 0.17;
!         B = 0.4348;
!         alpha0 = 1.358;
         
!         alpha = alpha0/(1+ A/Kn**B);

!         fc = (1+alpha*Kn)*(1+(4*Kn)/(1+Kn));

!         if ( correcaoPerm .EQV. .FALSE. ) fc = 1.d0


!      END SUBROUTINE calcular_fc
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
      
         use mParametros,   only: R_, T, M_m, p_Poco, phi_BM, M_m, p_Ref, Sg;   
         use mParametros,   only: tamBlocoMacro, widthBlocoMacro, dimZ_FB, gasTotalKg_BM, gasRecuperavelKg_BM    
         use mMalha,        only: numel_bm
         use mGlobaisArranjos, only: phi_n0
         
         implicit none;

         REAL*8 :: Z, gas_pInicial, gas_pFinal, gasRecuperavel, volBlocoMacro;
         REAL*8 :: pInicial, pFinal;
         integer :: nel
                 
         volBlocoMacro = tamBlocoMacro*widthBlocoMacro*dimZ_FB        
         
         pInicial = p_Ref;
         pFinal   = p_Poco;
         
         call calcularZ_P(pInicial, Z);
         !gas_pInicial = 0.0d0
         !do nel=1,numel_bm         
         !gas_pInicial   = gas_pInicial + pInicial /(Z*R_*T) * phi_n0(nel) * M_m * Sg      ! kg/m^3
         !enddo
         gas_pInicial   = pInicial /(Z*R_*T) * sum(phi_n0)/numel_bm * M_m * Sg      ! kg/m^3
         !write(*,*) "gas inicial", gas_pInicial, Sg; stop

         call calcularZ_P(pFinal, Z);
         gas_pFinal     = pFinal/(Z*R_*T) * (sum(phi_n0)/numel_bm) * M_m * Sg         ! kg/m^3
          
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

