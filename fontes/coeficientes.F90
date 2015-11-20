
!!     ************************************************************
!     SETEMBRO/2014
!     Patricia de Araujo Pereira
!     LNCC/MCT
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

          implicit none      
          
!            z5 = -4.06322962185335869981427622055802231171594155539599d-39
!            z4 =  1.14189959008041363237063527655648895160774149357430d-30
!            z3 =  3.91944878312003202672584073771593219565327535366406d-24
!            z2 =  2.83388110051296987909531716821951945116758484291032d-17
!            z1 = -1.57095513997258713424683964840489203140805329894647d-08
!            z0 =  1.00001023872218075538853554462548345327377319335938d+00

! 	  Com 8000 ft (Diego, Nov/2015)
	  z9 = -2.90573129727953312175829509245978417825986028996227d-70
          z8 =  8.13499075103263527060986550907300972154425807711961d-62
          z7 = -8.30225038048283519096658596614353319650038800235658d-54
          z6 =  2.83261422992340905091425760414792047206865882681902d-46
          z5 =  1.13472518847971988355465960830442672897432875130282d-38
          z4 = -1.26730975062993999865671634668364234118677853399583d-30
          z3 =  3.45227116546012348292606028195057766129119900949936d-23
          z2 =  1.79079325291113158679825354307108140333365178939115d-16
          z1 = -1.65273196375611158744394684014938912142156368645374d-08
          z0 =  9.99638754126298656466076408833032473921775817871094d-01
          
! 	  Com 11000 ft (Diego, Nov/2015)
! 	  z9 =  6.39483237420864910930183612678877089113028670654775d-71
!           z8 = -3.20253937685551764585985428715584100372702813815494d-62
!           z7 =  6.85714376284594195601576499610992440423411077162451d-54
!           z6 = -8.14915057589326184806939890111108464167770383535834d-46
!           z5 =  5.80217012974160089911722871660584617598017851222426d-38
!           z4 = -2.44438237174924652922747971215093361239997923877945d-30
!           z3 =  5.14742096171312987967394159435763691603295384960230d-23
!           z2 =  5.30570274192737297934501708967198304505487473732225d-17
!           z1 = -1.61419004953716593758516654284404512331718706263928d-08
!           z0 =  9.99390437848697854583690514118643477559089660644531d-01
          
!           z9 =  2.78608152547984387275994827989103503800895961244004d-64
!           z8 = -6.09323649801521485478153842790223057288978465738106d-57
!           z7 =  5.51799587603258674629858885157506796462970653070561d-50
!           z6 = -2.67454928998484373579014365893370223537765295772653d-43
!           z5 =  7.34636534777857941167560256235747313585751756848711d-37
!           z4 = -1.09664363705645203342903776415549662167531218147291d-30
!           z3 =  1.15415577070169971115279526227165936146371448235381d-23
!           z2 =  3.57037602392517225719877310943287917820229552309103d-16
!           z1 = -1.70434786570412165040379803750830722819387119670864d-08
!           z0 =  9.99999954649377187365644203964620828628540039062500d-01


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
!           g12 = -1.38160052222750001347197687667775071690919318320758d-80
!           g11 =  9.04333358064339883989882967669974811981731392092592d-73
!           g10 = -2.63167651977365443016057214137174114935020565537365d-65
!           g9  =  4.48977528416433969038538211715924750138053902045959d-58
!           g8  = -4.98613133019607965681575774156510643670153571610224d-51
!           g7  =  3.78441540696730459073274089202485282874411203744973d-44
!           g6  = -2.00473515565129994751061576269462897416988939942025d-37
!           g5  =  7.42508351366224362116338418165276716857543157944909d-31
!           g4  = -1.88564052183513636782716397785672528146411938512465d-24
!           g3  =  3.08078180695835491695697759838335903580542330030385d-18
!           g2  = -2.48015979772586369599702129678306107259458390679185d-12
!           g1  = -1.74086184794445392775620075342946080354522564448416d-06
!           g0  =  1.03564690440976665541938928072340786457061767578125d+01

       end subroutine inicializarCoeficientes_FormP
!     ************************************************************


     
!     ************************************************************
      SUBROUTINE calcularR_P(p,R)
      
         use mParametros,   only: Hch, v2_inf, PL_sat, R_, T, phi_M, Sw, phi_K, phi_N, rhoL
      
         REAL*8   :: p, R, Z, Gamma, G, H

         call calcularZ_P(p,Z)
         Gamma =  gm3*p*p*p + gm2*p*p + gm1*p + gm0
!          Gamma =  gm5*(p**5) + gm4*(p**4) + gm3*(p**3) + gm2*(p**2) + gm1*(p**1)  + gm0
         
         G     =  g12*(p**12) + g11*(p**11) + g10*(p**10) + g9*(p**9) + g8*(p**8) +       &
       &          g7*(p**7) + g6*(p**6) +         &
       &          g5*(p**5) + g4*(p**4) + g3*(p**3) + g2*(p**2) + g1*p    + g0      
         H     =  Hch*exp(v2_inf*(p-PL_sat)/(R_*T));

!          write(*,*) 'Z', Z
!          write(*,*) 'gamma', Gamma
!          write(*,*) 'G', G
!          write(*,*) 'H', H

         R     =  phi_M*(1-Sw)/(R_*T*Z) + phi_K*phi_N*G/(R_*T*Z)  + phi_M*Sw*rhoL*Gamma/H;     

!           write(*,*) 'R(P)', R
           
      END SUBROUTINE calcularR_P
!     ************************************************************
! 


     ! a pressao de entrada, p, é p* X p_ref
!     ************************************************************
      SUBROUTINE calcularK_P(p, K_eff)
      
         use mParametros,   only: K_, Hch, v2_inf, PL_sat, R_, T, rhoL, D
            
         REAL*8   :: p, K_eff, K_eff_Old, Z, Gamma, G, H, HLinha, GammaLinha, GammaDivH_Linha

!          Z     =  z5*(p**5)  + z4*(p**4)  + z3 *(p**3) + z2 *(p**2) + z1*p + z0
         Z = z9*(p**9) + z8*(p**8) + z7*(p**7) + z6*(p**6) + z5*(p**5) + & 
         & z4*(p**4) + z3*(p**3) + z2*(p**2)+ z1*p + z0
         
         Gamma =  gm3*(p**3) + gm2*(p**2) + gm1*p + gm0
         H     =  Hch*exp(v2_inf*(p-PL_sat)/(R_*T));

         HLinha     =  H*v2_inf/(R_*T);
         GammaLinha =  gmLinha3*(p**3) + gmLinha2*(p**2) + gmLinha1*(p) + gmLinha0 
         
         GammaDivH_Linha =  (GammaLinha*H - HLinha*Gamma)/(H*H)
         
!          K_eff_Old =  K_ * p/(R_*T*Z) + D * ( Gamma*rhoL/H )

         K_eff     =  K_ * p/(R_*T*Z) + D * ( Gamma*rhoL/H + p*rhoL * GammaDivH_Linha ) 

         if ( K_eff_Old-K_eff .GT. 1.d-12 ) write(1231,*) "+X+X p, ", p, K_eff_Old-K_eff

           
      END SUBROUTINE calcularK_P
!     ************************************************************

!      
! !     ************************************************************
!       SUBROUTINE calcularK_P(p, K_eff)
!       
!          use mParametros,   only: K_, Hch, v2_inf, PL_sat, R_, T, rhoL, D
!             
!          REAL*8   :: p, K_eff, Z, Gamma, G, H
! 
!          Z     =  z5*(p**5)  + z4*(p**4)  + z3 *(p**3) + z2 *(p**2) + z1*p + z0
!          Gamma =  gm3*(p**3) + gm2*(p**2) + gm1*p + gm0
!          H     =  Hch*exp(v2_inf*(p-PL_sat)/(R_*T));
!          K_eff =  K_ * p/(R_*T*Z) + D * Gamma*rhoL/H
! 
! !          write(*,*) "p    ", p
! !          write(*,*) "Z    ", Z
! !          write(*,*) "Gamma", Gamma
! !          write(*,*) "Henry", H
! 
! !          write(*,*) "K_", K_
! !          write(*,*) "D", D
! !          write(*,*) "K_eff  ", K_eff
! !          K_eff = 4.d-11/(R_*T);
!            
!       END SUBROUTINE calcularK_P
! !     ************************************************************


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


!     ************************************************************
      SUBROUTINE calcularFluxoGasDissolvido(p1Star, p2Star, gradP, J_gasDissolvido)
      
         use mParametros,   only: Hch, v2_inf, PL_sat, rhoL, p_Ref, tamBloco, D, R_, T
         IMPLICIT NONE                 
         
         REAL*8  ::  p1Star, p2Star, gradP, J_gasDissolvido
         REAL*8  ::  p, Gamma, GammaLinha,  H, HLinha, GammaDivH_Linha
         
         p = (p2Star+p1Star)/2
!          p = p*p_Ref
         
         Gamma      =  gm3*p*p*p + gm2*p*p + gm1*p + gm0
         GammaLinha =  gmLinha3*(p**3) + gmLinha2*(p**2) + gmLinha1*(p) + gmLinha0 
         H          =  Hch*exp(v2_inf*(p-PL_sat)/(R_*T));
         HLinha     =  Hch*exp(v2_inf*(p-PL_sat)/(R_*T))*v2_inf/(R_*T);
         
         GammaDivH_Linha =  (GammaLinha*H - HLinha*Gamma)/(H*H)
         
!          J_Old           = D * ( gamma*rhoL/H )                            * gradP * p_Ref / tamBloco         
!          J_gasDissolvido = D * ( gamma*rhoL/H + p*rhoL * GammaDivH_Linha ) * gradP * p_Ref / tamBloco         
         J_gasDissolvido = D * ( gamma*rhoL/H + p*rhoL * GammaDivH_Linha ) * gradP / tamBloco         
         
         
!          J_gasDissolvido = Gamma*rhoL/H * D * gradP * p_Ref/tamBloco       
         
 200  FORMAT(5(E25.18, 2X))
 

      END SUBROUTINE calcularFluxoGasDissolvido
!     ************************************************************


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
         
         implicit none;

         REAL*8 :: Z, gas_pInicial, gas_pFinal, gasRecuperavel, volBlocoMacro;
         REAL*8 :: pInicial, pFinal;
                 
         volBlocoMacro = tamBlocoMacro*widthBlocoMacro*dimZ_FB        
         
         pInicial = p_Ref;
         pFinal   = p_Poco;
         
         call calcularZ_P(pInicial, Z);         
         gas_pInicial   = pInicial /(Z*R_*T) * phi_BM * M_m       ! kg/m^3
         
         call calcularZ_P(pFinal, Z);
         gas_pFinal     = pFinal/(Z*R_*T) * phi_BM * M_m          ! kg/m^3
          
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

