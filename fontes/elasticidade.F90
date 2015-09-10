!=================================================================================
!
!         programa de elementos finitos em fortran 90, 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
!
!         novo projeto: 
!         Eduardo Garcia,        bidu@lncc.br [4]
!         Tuane Lopes,           tuane@lncc.br [5]
!
!         LNCC/MCT
!         Petropolis, 02.2014
!=================================================================================
!     
!   Adaptado por Diego (jul/2015) para o caso da Elasticidade Linear estacionária.
! 
      module mElasticidade

        real*8,    allocatable :: deslocamento(:,:),  fDeslocamento(:,:,:)
        integer*4, allocatable :: idDeslocamento(:,:)
        integer*4, pointer     :: idiagD(:), lmD(:,:,:)
        real*8,    pointer     :: alhs(:), brhs(:)
        integer*4              :: neqD, nalhsD
        integer*8              :: A_HYPRE, parcsr_A, b_HYPRE, par_b, u_HYPRE, par_u, solver
        integer*8              :: HYPRE_SYSTEM_ID
        integer*4              :: solver_id, precond_id
        integer*8              :: Clower, Cupper
        integer*4, pointer     :: rows(:)

        contains 
!
!**** new **********************************************************************
!	
!       Alterada por Diego para comportar o caso da elasticidade e as iteracoes no tempo.
!       (Ago/2015)
      subroutine montarSistEqAlgElasticidade(optsolver, opttype, u, f, alhs, brhs, &
                         id, lm, idiag, numnp, ndof, nlvect, nalhs, neq)
!
      use mAlgMatricial, only: load, dirichletConditions
      use mMalha,        only: x_bm, conecNodaisElem_BM, numel_bm, nen_bm, nsd_bm
!       use mHYPREsistema, only: finalizarMontagemSistemaHYPRE
     
!
      implicit none

      character(len=10), intent(in) :: optSolver
      integer*4, intent(in) :: ndof, neq, numnp, nlvect, optType
      real*8  :: u(ndof, numnp)
      real*8  :: f(ndof, numnp, nlvect)
      integer*4:: id(ndof, numnp) 
      !real*8,  pointer, intent(inout) :: alhs(:), brhs(:)
      !integer*4, pointer, intent(inout) :: lm(:,:,:), idiag(:)
      real*8,  pointer :: alhs(:), brhs(:)
      integer*4, pointer :: lm(:,:,:), idiag(:)
      integer*4            :: nalhs
      real*8 :: t1, t2, t3
      
      if (optType .eq. 1) goto 100
!       A montagem da estrutura se torna dispensavel em iteracoes maiores que 1. Diego (Ago/2015)
      if (optType .eq. 0) then
! 	elasticidade = 0.0
! 	fElasticidade = 0.0
	goto 200
      end if
!
      100 print*, " ..... montando sistema de equacoes da elasticidade linear"

      call montarEstrutDadosSistEqAlq(optSolver, u, f, alhs, brhs, id, lm, idiag, &
                                      numnp, ndof, nlvect, nalhs, neq)
      200 write(*,*) " +++ apos call montarEstrutDadosSistEqAlq("
!       write(222,*) "+++ optType = ", optType, "nlvect = ", nlvect
      if (nlvect.gt.0) call load                (id,f,brhs,ndof,numnp,nlvect)
    
      if (nlvect.gt.0) call dirichletConditions (id,u,f,ndof,numnp,nlvect)
      write(*,*) " +++ apos call dirichletConditions("
!       stop
! 
      300 call timing(t1)
!       stop
      call calcCoefSistAlgElasticidade (u, x_bm, alhs, brhs, idiag, lm, conecNodaisElem_bm, &
                                     numnp, numel_bm, nen_bm, nsd_bm, ndof, nalhs, neq )
      write(*,*) " +++ apos call calcCoefSistAlgElasticidade("
      call timing(t2)
      write(*,*) " calculo dos coeficientes, tempo = ", t2 - t1

!       call finalizarMontagemSistemaHYPRE  (A_HYPRE, parcsr_A, b_HYPRE, par_b, u_HYPRE, par_u)
!
      end subroutine montarSistEqAlgElasticidade
!
!**** new **********************************************************************
!
      subroutine montarEstrutDadosSistEqAlq(optSolver, u, f, alhs, brhs, id, lm, idiag, &
                                      numnp, ndof, nlvect, nalhs, neq)
!
      use mGlobaisEscalares, only: ndofD!, ndofF
      use mMalha,            only: x_bm, conecNodaisElem_BM, numel_bm, nen_bm, nsd_bm
      use mMalha,            only: formlm
      use mLeituraEscrita,   only: iecho
#ifdef withGaussSkyline
      use mAlgMatricial,     only: colht, diag
#endif
#ifdef withCSR
!       use mGlobaisEscalares, only: numCoefPorLinhaElasticidade, numCoefPorLinhaFluxo
!       use mMalha,            only: criarListaVizinhos, listaDosElemsPorNoCSR
!       use mMalha,            only: numConexoesPorElem
!       use mPardisoSistema,   only: criarPonteirosMatEsparsa_CSR
#endif
#ifdef withHYPRE
!       use mHYPREsistema, only: criarSistemaAlgHYPRE
!       use mHYPREsistema, only: mpi_comm
#endif
!
      implicit none
      integer*4                :: numnp, ndof, nlvect, neq
      integer*4, intent(inout) :: nalhs
      real*8                   :: u(ndof, numnp), f(ndof, numnp, nlvect)
      real*8,  pointer         :: alhs(:), brhs(:)
      integer*4, pointer       :: lm(:,:,:), idiag(:)
      integer*4                :: id(ndof, numnp)
      character(len=10)        :: optSolver
      integer*4                :: meanbwP 
      logical                  :: simetria
      real*8 :: t1, t2, t3
      integer*4 :: i, localSize
!
      write(*,'(a)',advance='NO') " em montarEstrutDadosSistEqAlq, "
      allocate(lm(ndof,nen_bm,numel_bm))
      call formlm(id,conecNodaisElem_BM,lm,ndof,ndof,nen_bm,numel_bm)

      allocate(idiag(neq)); idiag = 0
#ifdef withGaussSkyline
      call colht(idiag,lm,ndof,nen_bm,numel_bm,neq)
      call diag(idiag,neq,nalhs) ! only to cholesky+skyline solver 
      if(.not.associated(alhs)) allocate(alhs(nalhs)); alhs=0.0
#endif 

#ifdef withCSR
!       call timing(t1)
!       write(*,*) "withCSR, call criarPonteirosMatEsparsa_CSR" 
!       call criarListaVizinhos(nen, numnp, numel, conecNodaisElem, listaDosElemsPorNoCSR)
!       simetria=.true.
!       call criarPonteirosMatEsparsa_CSR(nsd, ndof, neq, numCoefPorLinhaElasticidade, &
!                                         conecNodaisElem, listaDosElemsPorNoCSR, id,&
!                                         numnp, nen, numConexoesPorElem, nalhs, simetria, 'p')
!       call timing(t2)
!       write(*,*) "criar lista de vizinhos e ponteiros matriz esparsa,  tempo = ", t2 - t1
!       if(.not.associated(alhs)) allocate(alhs(nalhs)); alhs=0.0
#endif 
! 
#ifdef withHYPRE
!       Clower = 1 - 1
!       Cupper = neq-1
!       localSize=CUpper-Clower+1
!       if(.not.associated(rows)) allocate(rows(localSize))
!       do i = 1, localSize
!        rows(i) = i - 1
!       end do
!       call criarSistemaAlgHYPRE(A_HYPRE, parcsr_A, b_HYPRE, par_b, u_HYPRE, par_u, &
!                                  solver, Clower, Cupper, mpi_comm)
#endif 

      if(.not.associated(brhs)) allocate(brhs(neq));
       brhs=0.0
! 
!       meanbwP = nalhsP/neqP
!       write(*,    6000) 'Calculo do elasticidade nodal', neqP, nalhsP, meanbwP, (8.0*nalhsP)/1000.0/1000.0
!       write(iecho,6000) 'Calculo do elasticidade nodal', neqP, nalhsP, meanbwP, (8.0*nalhsP)/1000.0/1000.0

 6000 format(a///&
     ' e q u a t i o n    s y s t e m    d a t a              ',  //5x,&
     ' number of equations . . . . . . . . . . . . (neq    ) = ',i8//5x,&
     ' number of terms in left-hand-side matrix  . (nalhs  ) = ',i12//5x,&
     ' mean half bandwidth . . . . . . . . . . . . (meanbw ) = ',i8//5x,&
     ' memoria necessaria para a matriz do sistema (Mbytes)  = ',e10.2)

      end subroutine montarEstrutDadosSistEqAlq
!
!**** new **********************************************************************
!
      subroutine calcCoefSistAlgElasticidade(d, x, alhs, brhs, idiag, lm, conecNodaisElem, &
                                     numnp, numel, nen, nsd, ndof, nalhs, neq )
!

      use mGlobaisEscalares, only: nrowsh_bm, npint_bm	! OK
      use mGlobaisArranjos,  only: mat_bm, c_bm, grav_bm, celast	! OK
      use mAlgMatricial,     only: kdbc, addrhs, addlhs	! OK
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shlq3d, shg3d, shgq, shlq	! OK
      use mMalha,            only: local	! OK
      use mBlocoMacro,       only: solucao_BM, ndof_bm
      use mParametros,       only: p_Ref, tamBlocoMacro, widthBlocoMacro
      !use mMalha,         only:  conecNodaisElem
!       use mPardisoSistema, only: addlhsCSR, ApElasticidade, AiElasticidade
!       use mHYPRESistema,   only: addnslHYPRE, addrhsHYPRE
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!       
      integer*4, intent(in)  :: numnp, numel, nen, nsd, ndof, nalhs, neq 
!
      real*8,  intent(inout) :: alhs(nalhsD),   brhs(neqD)
      real*8,  intent(in)    :: d(ndof, numnp), x(nsd, numnp)
      integer*4, intent(in)    :: idiag(neqD),    lm(ndof,nen,numel) 
      integer*4, intent(in)    :: conecNodaisElem(nen,numel)
!
      real*8 :: xl(nsd,nen), dl(ndof,nen), gpl(1,nen)
      real*8 :: shg(nrowsh_bm,nen,npint_bm), shl(nrowsh_bm,nen,npint_bm)
      real*8 :: det(npint_bm), w(npint_bm)
!
      integer*4 :: nee
      real*8   :: elresf(ndof*nen), eleffm(ndof*nen,ndof*nen)
!
      integer*4:: nel, m, l, i, j, k, ni, nj
      integer*4, parameter :: um = 1
      real*8  :: pi, Kx, Ky, Kz
      real*8  :: temp1, gf1, gf2, gf3
      real*8  :: pss
      real*8  :: djx, djy, djz, djn, dix, diy, diz, gpx, gpy
      logical :: diag,zerodl,quad,lsym
!       Variaveis para o estado plano de tensoes (TODO: Automatizar leitura)
      real*8  :: D1, CMATRIX11, CMATRIX12, CMATRIX33, YOUNG, POISSON, RHOMAT
!       *******************************************
      real*8  :: Lx, Ly

!
      Lx    = widthBlocoMacro
      Ly    = tamBlocoMacro 
	
      nee = nen*ndof
      
      diag = .false.
      pi=4.d00*datan(1.d00)
      YOUNG = celast(1)
      POISSON = celast(2)
      RHOMAT = celast(3)
!
      w=0.0
      shl=0.0

      if(nen==2) call oneshl(shl,w,npint_bm,nen) 
      if(nen==3) call   shlt(shl,w,npint_bm,nen)
      if(nen==4) call   shlq(shl,w,npint_bm,nen)
      if(nen==8) call shlq3d(shl,w,npint_bm,nen)
	
      do 500 nel=1,numel
!
!      clear stiffness matrix and force array
!
      eleffm=0.0
      elresf=0.0
!
!      LOCALIZE COORDINATes and Dirichlet b.c.
!
      call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
      call local(conecNodaisElem(1,nel),d,dl,nen,ndof,ndof)
      call local(conecNodaisElem(1,nel),solucao_BM,gpl,nen,ndof_bm,ndof_bm)
!
      m = mat_bm(nel)
      quad = .true.
      if (nen.eq.4.and.conecNodaisElem(3,nel).eq.conecNodaisElem(4,nel)) quad = .false.
!
      if(nen==2) call oneshg(xl,det,shl,shg,nen,npint_bm,nsd,um,nel,um)
      if(nen==4) call shgq  (xl,det,shl,shg,npint_bm,nel,quad,nen)
      if(nen==8) call shg3d (xl,det,shl,shg,npint_bm,nel,nen)
      
!
!.... SETUP ELASTICITY TENSOR COEFFICINTS: PLANE STRAIN 
!
         D1=YOUNG/((1.0d0+POISSON)*(1.0d0-2.0d0*POISSON))
         CMATRIX11 = D1*(1.0d0-POISSON)
         CMATRIX12=D1*POISSON
         CMATRIX33=D1*(1.0d0-2.0d0*POISSON)*0.5d0
!
!....... form stiffness matrix
!
!... length of the element
!
      Kx=c_bm(1,m)
      Ky=c_bm(2,m)
      if(nsd==3) Kz=c_bm(3,m)
!
!.... loop on integration points
!
      do 400 l=1,npint_bm

      temp1 = w(l)*det(l)
! 
!.... vetor de carga - RHS - f  =  gf0  
!
                 gf1 = grav_bm(1)
                 gf2 = grav_bm(2)
      if(nsd==3) gf3 = grav_bm(3)
      pss = 0.0
!       
      gpx = 0.0d0
      gpy = 0.0d0
      do 301 k=1,nen
	gpx = gpx + shg(1,k,l)*gpl(1,k)
	gpy = gpy + shg(2,k,l)*gpl(1,k)
      301 continue
!
      do 300 j=1,nen
           nj=ndof*j
                 djx=shg(1,j,l)*temp1
                 djy=shg(2,j,l)*temp1
      if(nsd==3) djz=shg(3,j,l)*temp1
      djn=shg(nrowsh_bm,j,l)*temp1

!     
!.... source terms      
!
      ELRESF(nj-1)= ELRESF(nj-1) + RHOMAT*GRAV_BM(1)*djn + gpx*djx*(p_Ref/Lx) !djx*gpl(1,J)*p_Ref/Lx
!
      ELRESF(nj)  = ELRESF(nj) + RHOMAT*GRAV_BM(2)*djn + gpy*djy*(p_Ref/Ly) !djy*gpl(1,J)*p_Ref/Ly
!
      do 300 i=1,nen
      ni = ndof*i
                 dix=shg(1,i,l)
                 diy=shg(2,i,l) 
      if(nsd==3) diz=shg(3,i,l)       
      
      ELEFFM(ni-1,nj-1)=ELEFFM(ni-1,nj-1)+DIX*(CMATRIX11*DJX)+DIY*(CMATRIX33*DJY)
!
      ELEFFM(ni-1,nj)=ELEFFM(ni-1,nj)+DIX*(CMATRIX12*DJY)+DIY*(CMATRIX33*DJX)
!
      ELEFFM(ni,nj-1)=ELEFFM(ni,nj-1)+DIY*(CMATRIX12*DJX)+DIX*(CMATRIX33*DJY)
!
      ELEFFM(ni,nj)=ELEFFM(ni,nj)+DIY*(CMATRIX11*DJY)+DIX*(CMATRIX33*DJX)
!
  300 continue
  400 continue  
!
!      computation of Dirichlet b.c. contribution
!
       call ztest(dl,nee,zerodl)
!
      if(.not.zerodl) then
          call kdbc(eleffm,elresf,dl,nee)
      endif
!
!.... assemble element stifness matrix and force array into global
!        left-hand-side matrix and right-hand side vector
      lsym=.true.
#ifdef withCSR
! !       call addlhsCSR(alhs,eleffm,      lm(1,1,nel),ApElasticidade,AiElasticidade,nee) 
!       call addlhsCSR(alhs,eleffm,      lm(1,1,nel),ApElasticidade,AiElasticidade,nee)
#endif
#ifdef withGaussSkyline
!       call addlhs(alhs,eleffm,idiag,lm(1,1,nel),nee,diag,lsym) 
      call addlhs(alhs,eleffm,idiag,lm(1,1,nel),nee,diag,lsym) 
#endif
#ifdef withHYPRE
! !       call addnslHYPRE(A_HYPRE, eleffm, lm(1,1,nel), nee, lsym) 
!       call addnslHYPRE(A_HYPRE, eleffm, lm(1,1,nel), nee, lsym)
!       !call addrhsHYPRE(B_HYPRE, brhs, elresf,lm(1,1,nel),nee)
#endif
!       call addrhs(brhs,elresf,lm(1,1,nel),nee)
      call addrhs(brhs,elresf,lm(1,1,nel),nee)

  500 continue

      return
      end subroutine calcCoefSistAlgElasticidade
      
      
      subroutine alocarMemoriaElasticidade()
      
!      Alocacao de memoria da Elasticidade Linear. Diego (Ago/2015)

     use mGlobaisEscalares, only: ndofD, nlvectD
     use mMalha,            only: numnp_BM

     implicit none
     
     allocate(  deslocamento(ndofD,NUMNP_BM));   deslocamento=0.0
     allocate(idDeslocamento(ndofD,NUMNP_BM));  idDeslocamento=0
     
     if (nlvectD.ne.0) then 
        allocate(fDeslocamento(nlvectD,ndofD,numnp_BM))
        fDeslocamento = 0.0
     endif
     
     end subroutine alocarMemoriaElasticidade

     end module
!