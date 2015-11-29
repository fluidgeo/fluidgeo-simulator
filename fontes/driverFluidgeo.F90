!  
!         programa de elementos finitos em fortran 90 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!
!         Eduardo Garcia, Tuane Lopes e Patrícia Pereira
!         bidu@lncc.br, tuane@lncc.br, 
!         Modificacoes para incluir fraturas naturais: Aline Rocha
!         Modificacoes e adicoes para incluir acoplamento hidro-geomecanico: Diego Volpatto
!         LNCC/MCT+
!     ************************************************************
!     *                                                          *
!     *                                                          *
!     *         A LINEAR STATIC FINITE ELEMENT PROGRAM FOR       *
!     *                                                          *
!     *                 GALERKIN METHOD                     *
!     *                                                          *
!     *                                                          *
!     ************************************************************
!
	program reservoirSimulator
	
	      use mGlobaisEscalares, only : exec
	      use mLeituraEscrita,   only : fecharArquivos, abrirArquivosDS
	      use mLeituraEscrita,   only : fecharArquivosDS, abrirArquivo
	      use mLeituraEscrita,   only : iechoPressao, iechoProducao
	      use mInputReader,      only : readInputParametersDS 
! 	      use mBloco,            only : iinBloco  =>iin, iechoBloco  =>iecho
! 	      use mBloco,            only : preprocessadorBloco !, processadorBloco 
	      use mBlocoMacro,       only : iinBlocoMacro=>iin, iechoBlocoMacro=>iecho
	      use mBlocoMacro,       only : preprocessadorBlocoMacro
! 	      use mBlocoMacro,       only : preprocessamentoBlocoMacro
	      use mCoeficientes,     only : inicializarCoeficientes_FormP
	      use mParametros,       only : inicializarParametros
              use mCoeficientes,     only : iechoQuantidadeGas
	!
	      implicit none
	      
	      logical :: comDS = .true.
	      
!
	!.... initialization phase
	
	!.... arquivos de entrada de dados
        !.... blocoMacro: blocos e fraturas naturais homogeneizados 
        !.... bloco: agregados orgânicos e microporos homogeneizados 
! 	      iinBlocoMacro      = abrirArquivo("inputBlocoMacro.dat", 16, 'old' ) 
! 	      iinBloco           = abrirArquivo("inputBloco.dat",   15, 'old' ) 
	      
	      
	!.... arquivos de saida 
! 	      iechoBlocoMacro    = abrirArquivo("echoBlocoMacro.dat",       18, 'replace' ) 
! 	      iechoBloco         = abrirArquivo("echoBloco.dat",         17, 'replace' )       
! 	      iechoQuantidadeGas = abrirArquivo("echoQuantidadeGas.dat", 5010, 'replace' ) 
! 	      iechoPressao       = abrirArquivo("echoPressao.dat",       123, 'replace' ) 
! 	      iechoProducao      = abrirArquivo("echoProducao.dat",      5001, 'replace' ) 
	      
	      
	!...  Este arquivo contém os passos de tempo para os quais o campo de pressão será gravado no arquivo de saída fort.123      
	!...  Se este arquivo não existir, todos os passos de tempo serão gravados. Como confiro se o arquivo foi aberto???
! 	      iinPassosPressaoBlocoMacro = abrirArquivo("passosPressaoBlocoMacro.dat", 19, 'old' ) 
	      
!               write(*,*) iechoBloco
!               write(*,*) iechoBlocoMacro

	      comDS = .TRUE.
	      if(comDS) then
		call abrirArquivosDS(); print*, " call abrirArquivosDS ()"
	      else
		call readInputFiles(); print*, " call readInputFiles ()"
	      endif

        !...  call abrirArquivosInformacoesMalha ()
	      print*, ""
	      print*, "PREPROCESSAMENTO" 
          
        !     parametros comuns a bloco e fratura 
	      write(*,*) "call readInputParametersDS()"
	      call readInputParametersDS()
	      write(*,*) "call readSetupParametersDS()"
	      call readSetupParametersDS()
	      write(*,*) "call inicializarParametros()"
              call inicializarParametros();
         !     coeficientes para os ajustes feitos para as curvas Z, gamma e G
              write(*,*) "call inicializarCoeficientes_FormP()"
              call inicializarCoeficientes_FormP();  

		if(comDS) then
		  call preprocessamentoDS(); print*, " call preprocessadorDS ()"
		else
		  call preprocessadorBlocoMacro(); print*, " call preprocessador ()"
		endif
	!
	!.... solution phase
	!
	      if(exec==1) then
		   print*, ""
 		   write(*,*) "Iniciando o PROCESSAMENTO.: Acoplamento hidro-geomecanico";  
 		   call processadorFluidgeo()
 	      endif
	!
	      call fecharArquivos()
	      if(comDS) then
		call fecharArquivosDS();! print*, " call fecharArquivosDS ()"
	      else
		call fecharArquivos();! print*, " call fecharArquivos ()"
	      endif
	 
	end program reservoirSimulator

	
! ********** new ********************************************************
!
! 	Essa subrotina tem como finalidade condensar as entradas padroes em uma so
! 	sub-rotina, visando a modularizacao.
! 	Diego (ago/2015)

	subroutine readInputFiles()
		
		use mGlobaisEscalares, only : exec
		use mLeituraEscrita,   only : abrirArquivosInformacoesMalha, fecharArquivos, abrirArquivo
! 		use mBloco,            only : iinBloco  =>iin, iechoBloco  =>iecho
		use mBlocoMacro,       only : iinBlocoMacro=>iin, iechoBlocoMacro=>iecho
		use mLeituraEscrita,   only : iechoPressao, iechoProducao
		use mCoeficientes,     only : iechoQuantidadeGas
		
		implicit none
		
		!.... initialization phase
		
		!.... arquivos de entrada de dados
		!.... blocoMacro: blocos e fraturas naturais homogeneizados 
		!.... bloco: agregados orgânicos e microporos homogeneizados 
		iinBlocoMacro      = abrirArquivo("inputBlocoMacro.dat", 16, 'old' ) 
	! 	      iinBloco           = abrirArquivo("inputBloco.dat",   15, 'old' ) 
		
		
		!.... arquivos de saida 
		iechoBlocoMacro    = abrirArquivo("echoBlocoMacro.dat",       18, 'replace' ) 
! 		iechoBloco         = abrirArquivo("echoBloco.dat",         17, 'replace' )       
		iechoQuantidadeGas = abrirArquivo("echoQuantidadeGas.dat", 5010, 'replace' ) 
		iechoPressao       = abrirArquivo("echoPressao.dat",       123, 'replace' ) 
		iechoProducao      = abrirArquivo("echoProducao.dat",      5001, 'replace' ) 
		
		
		!...  Este arquivo contém os passos de tempo para os quais o campo de pressão será gravado no arquivo de saída fort.123      
		!...  Se este arquivo não existir, todos os passos de tempo serão gravados. Como confiro se o arquivo foi aberto???
	! 	      iinPassosPressaoBlocoMacro = abrirArquivo("passosPressaoBlocoMacro.dat", 19, 'old' ) 
		
! 		write(*,*) iechoBloco
		write(*,*) iechoBlocoMacro
		
	end subroutine

! 
! 	
! 	
    subroutine preprocessamentoDS! (optSolver)
        use mInputReader,      only: readInputFileDS, readInputParametersDS
        use mInputReader,      only: leituraGeracaoCoordenadasDS
        use mInputReader,      only: leituraCodigosCondContornoDS
        use mInputReader,      only: leituraValoresCondContornoDS
        use mInputReader,      only: leituraGeracaoConectividadesDS
        use mGlobaisArranjos,  only: etime, title_BM, mat_BM
        use mGlobaisEscalares, only: exec, iprtin, npint_BM, ndofD, nlvectD
        use mBlocoMacro,       only: NDOF_BM, NLVECT_BM,alocarMemoriaBlocoMacro
        use mBlocoMacro,       only: passosTempoParaGravarSaida, numPassos, f_BM
        use mBlocoMacro,       only: iinPassosPressaoBlocoMacro
        use mParametros,       only: inicializarParametrosBlocoMacro
        use mMalha,            only: NUMNP_BM, NUMEL_BM, nen_bm, nsd_bm
        use mMalha,            only: x_BM, conecNodaisElem_BM
        use mMalha,            only: nelx_BM, nely_BM, nelz_BM
        use mLeituraEscrita,   only: iin, iecho, icoords, echo
        use mLeituraEscrita,   only: prntel, abrirArquivo
        use mCoeficientes,     only: calcularQuantidadeGasBlocoMacro
        use mAlgMatricial,     only: idiag_BM, id_BM, NEQ_BM,NALHS_BM
        use mElasticidade,     only: fDeslocamento, neqD, nalhsD, idDeslocamento, idiagD
        use mElasticidade,     only: alocarMemoriaElasticidade
!         use mPotencial,        only: numPassos,iinPassosPressaoHidro
!         use mPotencial,        only: passosTempoParaGravarSaida
!         use mFluxo,            only: fFluxo,     neqF, nalhsF, idFluxo,     idiagF
! #ifdef withCSR
!         use mMalha,            only: numConexoesPorElem
!         use mGlobaisEscalares, only: numCoefPorLinhaFluxo, numCoefPorLinhaPotencial
! #endif
        !
        implicit none

        !character(LEN=20), intent(in)  :: optSolver
        character(len=50) keyword_name
        !
        logical :: simetria, exist
        integer :: i, n

!.... input phase

!        call echo
        call readInputFileDS()
        call readSetupPhaseDS()

        etime = 0.0

!.... initialization phase
!
!
!....    set memory pointers for static data arrays,&
!        and call associated input routines
!
!       call alocarMemoria()	! Antigo, pro código pc com leitura do DS (Diego, ago/2015)
      call alocarMemoriaBlocoMacro()
      
      call alocarMemoriaElasticidade()
      
      call inicializarParametrosBlocoMacro()
     
      call calcularQuantidadeGasBlocoMacro()
          
!
!.... input coordinate data
      write(*,*) "call leituraGeracaoCoordenadasDS"
      call leituraGeracaoCoordenadasDS(x_BM,nsd_bm,NUMNP_BM, icoords, iprtin)

!    generation of conectivities
      write(*,*) "call leituraGeracaoConectividadesDS"
      keyword_name = "conectividades_nodais_bm"
      call leituraGeracaoConectividadesDS(keyword_name, conecNodaisElem_BM,mat_BM,nen_bm)

      if (iprtin.eq.0) call prntel(mat_BM,conecNodaisElem_BM,nen_bm,NUMEL_BM,1_4)
!
!.... input boundary condition data and establish equation numbers

        keyword_name = "codigos_cond_contorno_blocoMacro"
        write(*,*) "call leituraCodigosCondContornoDS(id_BM"
!       call leituraCodigosCondContorno(idPotencial,ndofP,numnp,neqP,iin, iecho, iprtin)
        call leituraCodigosCondContornoDS(keyword_name, id_BM, NDOF_BM, NUMNP_BM, n, iecho, iprtin)
        
        keyword_name = "codigos_cond_contorno_elasticidade"
        write(*,*) "call leituraCodigosCondContornoDS(idDeslocamento)"
!       call leituraCodigosCondContorno(idDeslocamento,ndofP,numnp,neqP,iin, iecho, iprtin)
        call leituraCodigosCondContornoDS(keyword_name, idDeslocamento, ndofD, NUMNP_BM, neqD, iecho, iprtin)
        
!         Eh necessario ou n? (Diego)
	neq_BM=n        
        allocate(idiag_BM(neq_BM));  idiag_BM=0

      if (NLVECT_BM.gt.0) then
        keyword_name = "valores_cond_contorno_blocoMacro"
        write(*,*) "call leituraValoresCondContornoDS(f_bm,ndof_bm"
!       call leituraValoresCondContorno(fPotencial,ndofP,numnp,0,nlvectP,iprtin)
        call leituraValoresCondContornoDS(keyword_name, f_BM, ndof_bm, NUMNP_BM, 1_4, NLVECT_BM, iprtin)
      end if
      
      if (nlvectD.gt.0) then
        keyword_name = "valores_cond_contorno_elasticidade"
        write(*,*) "call leituraValoresCondContornoDS(fElasticidade,ndofP,numnp,0,nlvectP,iprtin)"
!       call leituraValoresCondContorno(fElasticidade,ndofP,numnp,0,nlvectP,iprtin)
        call leituraValoresCondContornoDS(keyword_name, fDeslocamento, ndofD, NUMNP_BM, 1_4, nlvectD, iprtin)
      end if

!.... input element data
     write (*,*) "call leituraParamNumericosPropFisicaDS()"
     call leituraParamNumericosPropFisicaDS(NALHS_BM, NEQ_BM)
     
! 
! #ifdef withCSR
!       numConexoesPorElem=nen
!       if(nsd==2)numCoefPorLinhaPotencial=9
!       if(nsd==2)numCoefPorLinhaFluxo    =18
! 
!       if(nsd==3)numCoefPorLinhaPotencial=27
!       if(nsd==3)numCoefPorLinhaFluxo    =81
! #endif
 1000 format(20a4)
 2000 format(4i10,i10,11i10)
 3000 format(5x,&
     ' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,&
     ' execution code  . . . . . . . . . . . . . . (exec  ) = ',i10//5x,&
     '    eq. 0, data check                                   ',   //5x,&
     '    eq. 1, execution                                    ',  //5x,&
     ' input data print code . . . . . . . . . . . (iprtin ) = ',i10//5x,&
     '    eq. 0, print nodal and element input data           ',   //5x,&
     '    eq. 1, do not print nodal and element input data    ',   //5x, &
     ' number of space dimensions  . . . . . . . . (nsd    ) = ',i10)
 4000 format(5x,&
     ' number of nodal points  . . . . . . . . . . (numnp  ) = ',i10//5x,&
     ' number of elements      . . . . . . . . . . (numel  ) = ',i10//5x,&
     ' number of potencial load vectors   . . . . . (nlvectP) = ',i10//5x,&
     ' number of fluxos   load vectors   . . . . . (nlvectF) = ',i10//5x)
     
     
     !      TODO: ADAPTAR ISSO AQUI AO FORMATO NOVO (Diego, jul/2015)
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
     
    end subroutine preprocessamentoDS
    
    !> Efetua a leitura completa dos dados da etapa de setup.
    subroutine readSetupPhaseDS
        use mInputReader,      only:readStringKeywordValue, readIntegerKeywordValue, readIntegerKeywordValue
        use mGlobaisArranjos,  only: title_BM
        use mGlobaisEscalares, only: exec, iprtin, npint_BM, dimModelo, ndofD, nlvectD
        use mMalha,            only: nsd_bm, numnp_bm, numel_bm, nen_bm, numLados_BM
        use mMalha,            only: nelx_BM, nely_BM, nelz_BM
!         use mBloco,            only: numBlocos
        use mBlocoMacro,       only: nesd_BM, NDOF_BM, NLVECT_BM

        implicit none
        character(len=50) keyword_name
        character*4  default_title_value(20)

        !Reads the title
        keyword_name = "title_bm"
        default_title_value = "unknown title"
        call readStringKeywordValue(keyword_name, title_bm, default_title_value)
        !print*, "title_pc: ", title
        !Reads exec
        keyword_name = "exec_bm"
        call readIntegerKeywordValue(keyword_name, exec, 0_4)
        !print*, "exec_pc: ", exec
        !Reads iprtin
        keyword_name = "iprtin_bm"
        call readIntegerKeywordValue(keyword_name, iprtin, 0_4)
        !print*, "iprtin_pc: ", iprtin
        !Reads nsd
        keyword_name = "nsd_bm"
        call readIntegerKeywordValue(keyword_name, nsd_bm, 0_4)
        if (nsd_BM .eq. 1) then
		dimModelo = '1D'
	else if (nsd_BM .eq. 2) then
		dimModelo = '2D'
	endif
        !print*, "nsd_pc: ", nsd
        !Reads numnp
        keyword_name = "numnp_bm"
        call readIntegerKeywordValue(keyword_name, numnp_bm, 0_4)
        !Reads ndof
        keyword_name = "ndof_bm"
        call readIntegerKeywordValue(keyword_name, NDOF_BM, 0_4)
        !print*, "numnp_pc: ", numnp
        !Reads ndof
        keyword_name = "ndofD"
        call readIntegerKeywordValue(keyword_name, NDOFD, 0_4)
        !print*, "numnp_pc: ", numnp
        !Reads numel
        keyword_name = "numel_bm"
        call readIntegerKeywordValue(keyword_name, numel_bm, 0_4)
        !print*, "numel_pc: ", numel
        !Reads nen_pc
        keyword_name = "nen_bm"
        call readIntegerKeywordValue(keyword_name, nen_bm, 0_4)
        !print*, "nen_pc: ", nen
        !Reads npint_pc
        keyword_name = "numLados_BM"
        call readIntegerKeywordValue(keyword_name, numLados_BM, 0_4)
        !Reads npint_pc
        keyword_name = "npint_bm"
        call readIntegerKeywordValue(keyword_name, npint_BM, 0_4)
        !print*, "npint_pc: ", npint
        !Reads nlvectP
        keyword_name = "nlvect_bm"
        call readIntegerKeywordValue(keyword_name, NLVECT_BM, 0_4)
        !Reads nlvectE
        keyword_name = "nlvectD"
        call readIntegerKeywordValue(keyword_name, NLVECTD, 0_4)
        !Reads nelx_BM
        keyword_name = "nelx_BM"
        call readIntegerKeywordValue(keyword_name, nelx_BM, 0_4)
        !Reads nely_BM
        keyword_name = "nely_BM"
        call readIntegerKeywordValue(keyword_name, nely_BM, 0_4)
        !Reads nelz_BM
        keyword_name = "nelz_BM"
        call readIntegerKeywordValue(keyword_name, nelz_BM, 0_4)
        !print*, "nlvectP_pc: ", nlvectP
        
        nesd_BM = 2            ! +++ quem é este e pq não é lido de arquivo de entrada?

! 	numBlocos = numel_BM   ! +++ aqui mesmo devo inicializar ncelr?
	
	
        
        return
    end subroutine readSetupPhaseDS !**********************************************************************************

    !> Efetua a leitura completa dos dados da etapa de setup dos parametros. Diego (Nov/2015)
    subroutine readSetupParametersDS
        use mInputReader,      only: readStringKeywordValue, readrealkeywordvalue, readlogicalkeywordvalue
        use mInputReader,      only: findKeyword, file_lines
        use mParametros,       only: p_reservatorio, p_poco, tamBlocoMacro, widthBlocoMacro
        use mParametros,       only: constK_BM, k_s, Kbulk
        use mGlobaisArranjos,  only: reservoir_case, phi_range, coupling_mode
        use mGlobaisEscalares, only: random_porosity
        

        implicit none
        character(len=50) keyword_name
        character*4  default_input_value(20)
        integer*4 :: i, keyword_line, nLinhaArqInput

!       Reads the reservoir basis
        keyword_name = "reservoir_case"
        default_input_value = "unknown"
        call readStringKeywordValue(keyword_name, reservoir_case, default_input_value)
!         write(*,*) reservoir_case

!       Reads the reservoir basis
        keyword_name = "coupling_mode"
        default_input_value = "oneway"
        call readStringKeywordValue(keyword_name, coupling_mode, default_input_value)
!         write(*,*) reservoir_case
        
!       Reads reservoir pressure
        keyword_name = "p_reservatorio"
	call readRealKeywordValue(keyword_name, p_reservatorio, p_reservatorio)
! 	write(*,*) p_reservatorio
	
!       Reads well pressure
        keyword_name = "p_poco"
	call readRealKeywordValue(keyword_name, p_poco, p_poco)
! 	write(*,*) p_poco

!       Reads reservoir height, the same of fracture
        keyword_name = "tamBlocoMacro"
	call readRealKeywordValue(keyword_name, tamBlocoMacro, tamBlocoMacro)
! 	write(*,*) tamBlocoMacro
	
!       Reads reservoir width/2
        keyword_name = "widthBlocoMacro"
	call readRealKeywordValue(keyword_name, widthBlocoMacro, widthBlocoMacro)
! 	write(*,*) widthBlocoMacro
	
!       Reads reservoir permeability
        keyword_name = "k_bm"
	call readRealKeywordValue(keyword_name, constK_BM, constK_BM)
! 	write(*,*) constK_BM
	
!       Reads k matrix
        keyword_name = "k_s"
	call readRealKeywordValue(keyword_name, k_s, k_s)
! 	write(*,*) k_s

!       Reads K bulk
        keyword_name = "Kbulk"
	call readRealKeywordValue(keyword_name, Kbulk, Kbulk)
! 	write(*,*) Kbulk

	keyword_name = "phi_range"
        keyword_line = findKeyword(keyword_name)
        nLinhaArqInput = keyword_line
        read (file_lines(nLinhaArqInput:),7000) (phi_range(i),i=1,2)
        nLinhaArqInput = nLinhaArqInput + 1
!         write(*,*) phi_range
	
	7000 format(8f10.0)
        return
    end subroutine readSetupParametersDS !**********************************************************************************

    
!**** new **********************************************************************
!	
      subroutine processadorFluidgeo()

      use mBlocoMacro,       only : NITER_BM, NDOF_BM, NLVECT_BM, flux_BM      
      use mBlocoMacro,       only : DTEMPO_BM, solucao_BM, f_BM, SUM_NUMITER_BM, calcularFluxoMassico2d   
      use mBlocoMacro,       only : passosTempoParaGravarSaida, numPassos, calcularFluxoMassico2
      use mBlocoMacro,       only : PRINTSOL_BM, PRINTFLUXV_BM, PRINTSOL_BM2, calcularFluxoMassico, calcularGradiente_BM   
      use mMalha,            only : NUMNP_BM, NUMNP_B, NUMEL_BM, nen_bm, nsd_bm, nelx_BM, nely_BM
      use mMalha,            only : X_BM, conecNodaisElem_BM
      use mAlgMatricial,     only : NED_BM, ID_BM
      use mAlgMatricial,     only : FTOD 
      use mParametros,       only : p_Ref, p_Poco, gasTotalKg, gasRecuperavelKg, gasProduzidoKg, phi_BM
      use mCoeficientes,     only : calcularQuantidadeGasReservatorio
      use mGlobaisEscalares, only : ndofD, nlvectD
      use mGlobaisArranjos,  only : phi_n, phi_n0, coupling_mode
      use mLeituraEscrita,   only : PRINTDISP, PRINTSTRESS, PRINTPORO
!
      use mElasticidade,        only: montarSistEqAlgElasticidade, calcStressPoro!, CHECKEQ!, calcStressPoro2
      use mElasticidade,        only: deslocamento, fDeslocamento, neqD, nalhsD
      use mElasticidade,        only: alhsD=>alhs, brhsD=>brhs
      use mElasticidade,        only: idDeslocamento, idiagD, lmD
!
      
      implicit none

      INTEGER  :: NITERMAX;            !   NUMERO MAXIMO DE ITERACOES DE PICARD ITERATIONS PARA O LOOP DA NAO LINEARIDADE ENTRE A FRATURA E O BLOCO 
      INTEGER  :: NUMITER;             !   contador de iteracoes de picard PARA O LOOP DA NAO LINEARIDADE ENTRE A FRATURA E O BLOCO
      INTEGER  :: NUSTEP, N;           !   contador de iteracoes (passos) de tempo
      INTEGER  :: NITER;
      REAL*8   :: TEMPO, DTEMPO;
      REAL*8   :: solucaoNaoLinearAnt (NDOF_BM, NUMNP_BM)
      REAL*8   :: solucaoTmpAnt       (NDOF_BM, NUMNP_BM) 
      REAL*8   :: RELAX_FB
      REAL*8   :: tempoSimulacaoTotal, DTEMPO_BASE
      LOGICAL  :: flagTempoCaract,flagConvergencia, ERRO
      
      INTEGER  :: numNosFratura, numBlocos
      
      REAL*8, allocatable :: fonteDeBlocoParaFratura (:)
      REAL*8, allocatable :: fluxoMassicoDeBlocoParaBlocoMacro (:)
      REAL*8, allocatable :: condContornoDeBlocoMacroParaBloco (:)
      
      REAL*8  :: maiorFonte, F, dia, mes, ano, ctolAnt
      REAL*8  :: volGasAcumulado, fluxoAnt
      INTEGER :: I, J, SUM_NUMITER, SUM_NUMITER_ANT_B, SUM_NUMITER_ANT_BM
      INTEGER :: idx   ! indice para os PassosTempoParaGravarSaida
      
      logical :: firstD=.true.
      character(LEN=20)  :: optSolver='GaussSkyline'
      real*8 :: t1, t2, t3, t4
      real*8 :: stressD(nen_bm,numel_bm), gradpElem(nsd_BM,NUMEL_BM)
      logical :: simetria, tflag, onlydisp
      character(LEN=12) :: label, etapa
      integer :: optType = 1
      integer :: optType_BM = 1
      
      onlydisp = .false.
      
      call printporo(phi_n0,numel_bm, 0)
      
      maiorFonte = -10d-9;
      numBlocos = NUMNP_BM;
      allocate(fluxoMassicoDeBlocoParaBlocoMacro(numBlocos));
      fluxoMassicoDeBlocoParaBlocoMacro = 0.d0;
      allocate(condContornoDeBlocoMacroParaBloco(numBlocos));
      
      volGasAcumulado = 0.d0
      fluxoAnt        = 0.d0      
      
      dia = 24.d0*60.d0*60.d0;
      mes = 30.d0*dia;
      ano = 12*mes;
      
      idx = 1;
       
      call calcularQuantidadeGasReservatorio()
             
      NITERMAX = 500    
      NUMITER  = 0       
      TEMPO    = 0.D0  
      SUM_NUMITER = 0;
      SUM_NUMITER_BM = 0;
      SUM_NUMITER_ANT_BM = 0;
      
      CALL     colocarCondInicial(solucao_BM, NDOF_BM, NUMNP_BM,  1*p_Ref);  ! solucao_BM = 1.0d0

      IF (NLVECT_BM.GT.0) CALL FTOD(id_BM,solucao_BM,f_BM,NDOF_BM,NUMNP_BM,NLVECT_BM)
      
      solucaoNaoLinearAnt = solucao_BM
      solucaoTmpAnt       = solucao_BM  

!C....    loop das iteracoes no tempo

       NITER  = NITER_BM
       DTEMPO = DTEMPO_BM
       DO NUSTEP=1,NITER
       
          SUM_NUMITER_ANT_BM = SUM_NUMITER_BM
          if ( (numPassos .EQ. 0) .OR. (passosTempoParaGravarSaida(idx) .EQ. NUSTEP) ) then
                 tflag = .true.
          else
		 tflag = .false.
          end if
             
          TEMPO=TEMPO+DTEMPO;                    
          write(*,*) 'PASSO', NUSTEP, 'TEMPO', TEMPO, 'DT', DTEMPO  
          
          flagConvergencia = .FALSE. 
          NUMITER          = 0
              condContornoDeBlocoMacroParaBloco=0  
              CALL processadorBlocoMacro    (solucaoTmpAnt, NUSTEP, TEMPO, fluxoMassicoDeBlocoParaBlocoMacro, &
       &                                     condContornoDeBlocoMacroParaBloco, DTEMPO)
          IF ( NUMITER.GT.NITERMAX ) THEN
               write(*,*) 'LOOP PRINCIPAL ABORT: NUMBER OF ITERATIONS,', NUMITER, ', EXCEEDED no PASSO', NUSTEP;  
               STOP;
          ELSE 
          
!          Se chegar aqui é pq convergiu, então : guardar o fluxo, resolver novamente o problema nos blocos, atualizar variaveis e imprimir           
!          +++++++ calcular e armzenar o fluxo aqui 
             call calcularGradiente_BM(NED_BM, NDOF_BM)

             CALL calcularFluxoMassico2d(flux_BM, solucao_BM, X_BM, NUSTEP, TEMPO, NUMNP_BM, DTEMPO, &
      &                                 fluxoAnt, volGasAcumulado)          
      
             SUM_NUMITER_ANT_BM = SUM_NUMITER_BM  
             
             SUM_NUMITER = SUM_NUMITER + NUMITER

!           if (coupling_mode .eq. "oneway") then
!           !
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULO DA ELASTICIDADE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	  !
! 	  deslocamento = 0.0
!           if(firstD) then
!                 etapa = 'full'
! 		call montarSistEqAlgElasticidade(optSolver, optType, deslocamento, fDeslocamento, alhsD, brhsD,  &
! 					idDeslocamento, lmD, idiagD, numnp_bm, ndofD, nlvectD, nalhsD, neqD,etapa)
! 		firstD=.false.
! 		optType = 0
! 	  else
! ! 		CLEAR necessario para nao carregar a solucao anterior dos deslocamentos nas solucoes
! ! 		seguintes. Diego (Ago/2015)
! ! 		deslocamento = 0.0
!                 etapa = 'back'
! 		call montarSistEqAlgElasticidade(optSolver, optType, deslocamento, fDeslocamento, alhsD, brhsD,  &
! 					idDeslocamento, lmD, idiagD, numnp_bm, ndofD, nlvectD, nalhsD, neqD,etapa)
! 	  endif
! 
! 	  call timing(t2)
! 	  call timing(t3)
!           label='elasticidade linear'
!           call solverD() 
!           call timing(t4)          
! !
! !	Fim do calculo da ELASTICIDADE
! !
          stressD = 0.0d0
          
          call calcStressPoro(stressD, deslocamento,solucao_BM, flux_BM, x_bm, conecNodaisElem_bm, &
          &                   numnp_BM, numel_bm, nen_bm, nsd_bm, ndofD, tflag, idx)
!           endif   
             if ( (numPassos .EQ. 0) .OR. (passosTempoParaGravarSaida(idx) .EQ. NUSTEP) &
             & .and. coupling_mode .eq. "oneway" ) then
! 		 Rotinas de impressao
                 CALL PRINTSOL_BM2(solucao_BM,    x_BM,NUMNP_BM,TEMPO,idx)
                 CALL PRINTFLUXV_BM(flux_BM,    x_BM,NUMNP_BM,TEMPO,idx)
                 CALL calcularFluxoMassico2(FLUX_BM, solucao_BM, solucaoTmpAnt, X_BM, TEMPO, DTEMPO, NUMNP_BM, idx)
                 call PRINTDISP(deslocamento,X_BM,NUMNP_BM,nelx_BM,nely_BM,idx)
                 call PRINTSTRESS(stressD,X_BM,NUMEL_BM,nelx_BM,nely_BM, idx)
                 call printporo(phi_n,numel_bm, idx)
                 idx = idx + 1
             else
! 		 Rotinas de impressao
                 CALL PRINTSOL_BM2(solucao_BM,    x_BM,NUMNP_BM,TEMPO,idx)
                 CALL PRINTFLUXV_BM(flux_BM,    x_BM,NUMNP_BM,TEMPO,idx)
                 CALL calcularFluxoMassico2(FLUX_BM, solucao_BM, solucaoTmpAnt, X_BM, TEMPO, DTEMPO, NUMNP_BM, idx)
                 idx = idx + 1
             end if
             
             solucaoTmpAnt = solucao_BM                     
             write(*,*) 'FIM DO TIME STEP =', NUSTEP
             
          END IF    

      ENDDO 

      write(*,*) " "
      write(*,150) "---> Pressão no reservatório:       ", p_Ref, "(Pa)"
      write(*,150) "---> Pressão no poço:               ", p_Poco, "(Pa)"
      write(*,150) "---> Quantidade de gás total:       ", gasTotalKg, "(kg)"
      write(*,150) "---> Quantidade de gás recuperável: ", gasRecuperavelKg, "(kg)"
      write(*,150) "---> Quantidade de gás produzido:   ", gasProduzidoKg, "(kg)"
      write(*,150) "---> Valor da pressão no último nó: ", solucao_BM(1,NUMNP_BM), "(Pa)"
      write(*,150) "---> Tempo de simulação             ", TEMPO/ano, "(anos)"
      write(*,100) "---> Fator de recuperação (sobre gás total):       ", gasProduzidoKg/gasTotalKg
      write(*,100) "---> Fator de recuperação (sobre gás recuperável): ", gasProduzidoKg/gasRecuperavelKg      
      
 100  FORMAT(A,1PE15.5)               
 150  FORMAT(A,1PE15.5,A)               
      
      RETURN
      
      contains
      
      subroutine solverD() 
      use mMalha,            only: NUMEL_BM, nen_BM
      use mElasticidade, only:rows, solver_id, precond_id
      use mAlgMatricial, only: solverGaussSkyline, btod, factor, back

      REAL*8,  allocatable :: initialGuessP(:)

      integer * 4 ::  num_iterations, localMyID = 0, printSol = 2
      real*8 ::  final_res_norm, elapsedT, tol

      
      write(*,*) 'solver ', optSolver ,', ',  label
      write(*,*) 'etapa ', etapa

!       etapa = 'factor'
      if (etapa .eq. 'full' .or. etapa .eq. 'factor') then
	call factor(alhsD, idiagD, nalhsD, neqD)
      endif
!       etapa = 'back'
      if (etapa .eq. 'full' .or. etapa .eq. 'back') then
	call back(alhsD, brhsD, idiagD, neqD)
      endif
      
      call btod(idDeslocamento,deslocamento, brhsD,ndofD,numnp_bm)

      brhsD = 0.0 ! ladoB
!       alhsD = 0.0
     
      end subroutine solverD
!       end subroutine processadorFluidgeo
!
!**** new **********************************************************************
!
      subroutine processadorBlocoMacro(solucaoTmpAnt, NUSTEP, TEMPO, fluxoMassicoDeBlocoParaBlocoMacro, &
      &                                condContorno, DT)

      use mBlocoMacro,       only : solucao_BM, f_BM, NDOF_BM, NLVECT_BM, flux_BM, SUM_NUMITER_BM
      use mBlocoMacro,       only : solucaoNaoLinearAnt_BM, solucaoTmpAnt_BM
      use mAlgMatricial,     only : id_BM, ALHS_BM, DLHS_BM, BRHS_BM, NED_BM, IDIAG_BM, NEQ_BM
      use mMalha,            only : NUMNP_BM, X_BM, NUMEL_BM
      use mAlgMatricial,     only : FTOD, FACTNS, BACKNS, BTOD, LOAD
      use mBlocoMacro,       only : montarSistema_BM

      implicit none

      REAL*8      :: solucaoTmpAnt(NDOF_BM, NUMNP_BM) 
      INTEGER     :: NUSTEP;
      REAL*8      :: TEMPO, relax_BM, DT;
      LOGICAL     :: ERRO      
      
      REAL*8      :: fluxoMassicoDeBlocoParaBlocoMacro(NUMNP_BM)
      REAL*8      :: condContorno(NUMEL_BM)      
      INTEGER     :: NITERMAX_BM, NUMITER_BM, I, J
      LOGICAL     :: flagConvergencia, flagTempoCaract, tflag2
      
      tflag2 = .false.
      
      NITERMAX_BM = 5000  !     Numero maximo de iteracoes para SET NEWTON/PICARD 
      
      flagConvergencia=.FALSE.             
      flagTempoCaract=.FALSE.
      NUMITER_BM=0                    !CONTADOR DE ITERACOES DE PICARD/NEWTON
      
      ! Isso aqui eh desnecessario.
      solucaoNaoLinearAnt_BM = solucao_BM;
      solucaoTmpAnt_BM       = solucaoTmpAnt;

      !=========================================== loop das iteracoes de Picard/Newton 
      DO WHILE ( (flagConvergencia .eqv. .FALSE.) .AND. (NUMITER_BM.LE.NITERMAX_BM) )      
      
         NUMITER_BM=NUMITER_BM+1                  
         write(*,*) "---------ITERACAO Bloco Macro", NUMITER_BM
           
!           CLEAR LEFT AND RIGHT HAND SIDE FOR THE TRANSIENT PROBLEM
         ALHS_BM=0.d0
         DLHS_BM=0.d0
         BRHS_BM=0.d0
         
!           ACCOUNT THE NODAL FORCES IN THE R.H.S.
         IF (NLVECT_BM.GT.0) CALL LOAD(id_BM,                 f_BM, BRHS_BM,NDOF_BM,NUMNP_BM,NLVECT_BM)

         IF (NLVECT_BM.GT.0) CALL FTOD(id_BM,solucao_BM,      f_BM,        NDOF_BM,NUMNP_BM,NLVECT_BM)        
         !! +++ acho que aqui solucao_F ja deveria ter o valor da cc correta, nao sendo necessario fazer FTOD. Checar, but it does not matter

         
         CALL montarSistema_BM(NED_BM,  NDOF_BM, fluxoMassicoDeBlocoParaBlocoMacro, DT)  ! ASSEMBLE STIFNESS MATRIX
         
         CALL factns          (ALHS_BM, DLHS_BM,    idiag_BM,neq_BM)            ! FACTORIZATION OF THE STIFFNES MATRIX

         call backns          (ALHS_BM, DLHS_BM,    BRHS_BM, idiag_BM,neq_BM)     ! BACK SUBSTITUTION

         CALL BTOD            (id_BM,   solucao_BM, BRHS_BM,NDOF_BM,NUMNP_BM)
         
         
         DO J=1, NDOF_BM
            DO I=1, NUMNP_BM
               IF ( solucao_BM(J,I) .NE. solucao_BM(J,I) ) then
                  ERRO = .TRUE.
                  write(*,*) 'ABORT (Fratura) - NAN: pressao no noh', I, solucao_BM(J,I)
                  STOP;
               END IF    
            END DO   
         END DO 

         CALL CHECKCONV   (solucao_BM, solucaoNaoLinearAnt_BM,NDOF_BM,NED_BM,NUMNP_BM,NUMITER_BM, &
     &                     flagConvergencia, NUSTEP)
          relax_BM = 0.5
          solucaoNaoLinearAnt_BM = solucao_BM*relax_BM + solucaoNaoLinearAnt_BM*(1-relax_BM)
          
          if (coupling_mode .eq. "oneway") then
          !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULO DA ELASTICIDADE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !
	  deslocamento = 0.0
          if(firstD) then
                etapa = 'full'
		call montarSistEqAlgElasticidade(optSolver, optType, deslocamento, fDeslocamento, alhsD, brhsD,  &
					idDeslocamento, lmD, idiagD, numnp_bm, ndofD, nlvectD, nalhsD, neqD,etapa)
		firstD=.false.
		optType = 0
	  else
! 		CLEAR necessario para nao carregar a solucao anterior dos deslocamentos nas solucoes
! 		seguintes. Diego (Ago/2015)
! 		deslocamento = 0.0
                etapa = 'back'
		call montarSistEqAlgElasticidade(optSolver, optType, deslocamento, fDeslocamento, alhsD, brhsD,  &
					idDeslocamento, lmD, idiagD, numnp_bm, ndofD, nlvectD, nalhsD, neqD,etapa)
	  endif

	  call timing(t2)
	  call timing(t3)
          label='elasticidade linear'
          call solverD() 
          call timing(t4)          
          stressD = 0.0d0
          
          call calcStressPoro(stressD, deslocamento,solucao_BM, flux_BM, x_bm, conecNodaisElem_bm, &
          &                   numnp_BM, numel_bm, nen_bm, nsd_bm, ndofD, tflag2, idx)
!
!	Fim do calculo da ELASTICIDADE
!
	  endif
      END DO      
      !===========================================  fim do loop das iteracoes de Picard/NEWTON 
         
      IF (NUMITER_BM.GT.NITERMAX_BM) THEN 
         write(*,*) 'Bloco Macro processadorBlocoMacro() - NUMBER OF ITERATIONS EXCEEDED'
         STOP
      ELSE 
         write(*,*) NUSTEP, "---FIM  ITERACAO Bloco Macro ", NUMITER_BM, "iteracoes";
         
         DO I=1, NUMNP_BM
             condContorno(I) = solucao_BM(1,I)
         END DO        
         SUM_NUMITER_BM = SUM_NUMITER_BM + NUMITER_BM;         
      END IF   

      end subroutine processadorBlocoMacro
      
      end subroutine processadorFluidgeo
!**** new **********************************************************************
      subroutine colocarCondInicial(solucao, NDOF, NUMNP, condInicial)
      
         IMPLICIT   NONE

         INTEGER :: NDOF, NUMNP;
         REAL*8  :: solucao(NDOF, NUMNP);
         REAL*8  :: condInicial

         INTEGER I,J;

         DO I=1,NDOF
            DO J=1,NUMNP
               solucao(I,J) = condInicial;
            ENDDO
         ENDDO 

      end subroutine colocarCondInicial
!*******************************************************************************


!
!********************************************************************
!
      SUBROUTINE CHECKCONV(D,DD,NDOF,NED,NUMNP,NUMITER,CHKCONV,NUSTEP)
!
!     Esta rotina faz o controle de convergencia das iteracoes
!     de newton. retorna .TRUE. if iterations converged
!
      IMPLICIT REAL*8  (A-H,O-Z) 
      LOGICAL CHKCONV
      INTEGER :: NDOF,NED,NUMNP,NUMITER,NUSTEP
!
!     REMOVE ABOVE CARD FOR SI7NGLE PRECISION OPERATION
!
      REAL*8 :: D(NDOF,NUMNP),DD(NDOF,NUMNP),CTOL(2)
!
      INTEGER ::I,J
!   
      DTOL=1.0D-4
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

      IF(CTOL(1).GT.DTOL) RETURN
!  
      CHKCONV=.TRUE.
!
      RETURN
!
      END SUBROUTINE

!**** new **********************************************************************
    !> Faz a leitura de constant body forces.
    !> @param keyword_name  Keyword especifica para constant body forces.
!     Adaptado por Diego (Ago/2015).
    subroutine leituraParamNumericosPropFisicaDS(NALHS_BM, NEQ_BM)
        use mLeituraEscrita,   only: iecho
!         
        use mBlocoMacro,       only: DTEMPO_BM, NITER_BM, NEE_BM,NEESQ_BM,NENP_BM
!         
        use mInputReader,       only: findKeyword, file_lines
        use mInputReader,       only: readrealkeywordvalue,readintegerkeywordvalue
!         
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
	use mLeituraEscrita, only: prntel, printD, iin
	use mAlgMatricial,   only: colht, diag


        !
        !.... program to read, generate and write element data
        !
        implicit none
        !
        !.... remove above card for single precision operation
        !
        !
        integer*4 :: m, n, i, keyword_line, nLinhaArqInput
        character(len=80) :: formatoLeitura
        character(len=50) keyword_name
        integer, intent(in) ::NALHS_BM, NEQ_BM
	real*8, allocatable :: xl(:,:) 

        keyword_name = "nummat_bm"
        keyword_line = findKeyword(keyword_name)
        nLinhaArqInput = keyword_line
        if (keyword_line.eq.-1) return

        nrowsh_bm = 3
        if (nsd_BM==3) nrowsh_bm=nrowsh_bm+1
        ned_BM = 1
        NEE_BM    = NEN_BM*NED_BM
	NEESQ_BM  = NEE_BM*NEE_BM
	NENP_BM   = nen_BM
! 

        allocate(npar_BM(numParElem))
        allocate(celast(3))
        !
        formatoLeitura='(16I10)'
        read(file_lines(nLinhaArqInput:), formatoLeitura) (npar_bm(i),i=1,numParElem)
        nLinhaArqInput = nLinhaArqInput + 1

        nicode_BM = npint_bm
        if (nicode_bm.eq.0) nicode_bm=nen_bm
        !
        numat_bm  = npar_bm( 1)
        allocate(c_bm(6,numat_bm))
        !
        write(iecho,1000) numel_bm,numat_bm,nen_bm,npint_bm
        !
        !      read material properties
        !
        keyword_name = "prop_fisica_meio_bm"
        keyword_line = findKeyword(keyword_name)
        nLinhaArqInput = keyword_line
        do 400 n=1,numat_bm
        if (mod(n,50).eq.1) write(iecho,4000) numat_bm
        read(file_lines(nLinhaArqInput:),5000) m,(c_bm(i,m),i=1,3)
        nLinhaArqInput = nLinhaArqInput + 1
        write(iecho,6000) m,(c_bm(i,m),i=1,3)
        400 continue
        !
        !     constant body forces
        !
        keyword_name = "grav_bm"
        keyword_line = findKeyword(keyword_name)
        nLinhaArqInput = keyword_line
        read (file_lines(nLinhaArqInput:),7000) (grav_bm(i),i=1,3)
        nLinhaArqInput = nLinhaArqInput + 1
        write (iecho,8000) (grav_bm(i),i=1,3)
        !
        !     transient data
        !
        keyword_name = "transient_data_dt_bm"
	call readRealKeywordValue(keyword_name, DTEMPO_BM, DTEMPO_BM)
	
	keyword_name = "transient_data_nstep_bm"
	call readIntegerKeywordValue(keyword_name, NITER_BM, NITER_BM)
        
        write(iecho,12000) DTEMPO_BM, NITER_BM
        
        !
        !     material constants (to linear elasticity problem)
        !
        keyword_name = "constantes_elasticidade"
        keyword_line = findKeyword(keyword_name)
        nLinhaArqInput = keyword_line
        read (file_lines(nLinhaArqInput:),7000) (celast(i),i=1,3)
        nLinhaArqInput = nLinhaArqInput + 1
        write (iecho,8500) (celast(i),i=1,3)
	
#ifdef debug
!      if(iprtin.eq.0) then
!         call prntel(mat_BM,conecNodaisElem_BM,nen_BM,numel_BM, 1, iecho)
        call prntel(mat_BM,conecNodaisElem_BM,nen_BM,numel_BM, 1)
!      end if
#endif
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

        1000 format(//,&
        ' two/three-n o d e    e l e m e n t s ',//,5x,&
        ' number of elements  . . . . . . . . . . . (numel ) = ',i8,//5x,&
        ' number of element material sets . . . . . (numat ) = ',i5,//5x,&
        ' number of element nodes . . . . . . . . . (nen   ) = ',i5,//5x,&
        ' number of integration points. . . . . . . (npint  ) = ',i5)
        2000  format(16i5)
        4000  format(///,&
        ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
        ' number of material sets . . . . . . . . . . (numat ) = ', i5///,&
        2x,'set',4x,'Kx ', 10x,'Ky',10x,'Kz')
        5000  format(i10,5x,5f10.0)
        6000  format(2x,i3,1x,5(1x,1pe14.4))
        7000 format(8f10.0)
        8000 format(///,&
        ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
        ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
        ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
        ' exemplo 3............................ = ',      1pe15.8,//)
        8500 format(///,&
        ' m a t e r i a l   c o n s t a n t s                   ',//5x,&
        ' Young . . . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
        ' Poisson . . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
        ' Rho.................................. = ',      1pe15.8,//)
        9000  format(i5)
        12000 FORMAT(5X,  &
	&'STEP TIME . . . . . . . . . . . . . . . . .(DTEMPO)=',1PE15.8 //5X, &
	&'NUMBER OF STEPS. . . . . . . . . . . . . . . (NITER ) = ',I5)
        !
    end subroutine leituraParamNumericosPropFisicaDS
