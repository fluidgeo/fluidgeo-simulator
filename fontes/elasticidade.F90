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

        real*8,    allocatable :: deslocamento(:,:), deslocamentoAnt(:,:), fDeslocamento(:,:,:)
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
                    id, lm, idiag, numnp, ndof, nlvect, nalhs, neq, &
                     label, flagT)
!
      use mAlgMatricial, only: load, dirichletConditions
      use mMalha,        only: x_bm, conecNodaisElem_BM, numel_bm, nen_bm, nsd_bm
!       use mHYPREsistema, only: finalizarMontagemSistemaHYPRE
     
!
      implicit none

      character(len=10), intent(in) :: optSolver
      character(len=12), intent(in) :: label
      logical, intent(in)  :: flagT
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
      if (optType .eq. 0) goto 200
!
      100 print*, " ..... montando sistema de equacoes da elasticidade linear"

      call montarEstrutDadosSistEqAlq(optSolver, u, f, alhs, brhs, id, lm, idiag, &
                                      numnp, ndof, nlvect, nalhs, neq, label)
      200 write(*,*) " +++ apos call montarEstrutDadosSistEqAlq("
      if (nlvect.gt.0) call load                (id,f,brhs,ndof,numnp,nlvect)
      if (nlvect.gt.0) call dirichletConditions (id,u,f,ndof,numnp,nlvect)
      !write(1,*) f; stop
      !stop 
      write(*,*) " +++ apos call dirichletConditions("

! 
      call timing(t1)
      if (label .eq. 'factor' .or. label .eq. 'full') alhs=0.0;
      if (flagT .eqv. .true.) then
      call calcCoefSistAlgTerzaghi (u, x_bm, alhs, brhs, idiag, lm, conecNodaisElem_bm, &
                                     numnp, numel_bm, nen_bm, nsd_bm, ndof, nalhs, neq, label)
      else
      call calcCoefSistAlgElasticidade (u, x_bm, alhs, brhs, idiag, lm, conecNodaisElem_bm, &
                                     numnp, numel_bm, nen_bm, nsd_bm, ndof, nalhs, neq, label)
      endif
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
                                      numnp, ndof, nlvect, nalhs, neq, label)
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
      character(len=12), intent(in) :: label
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
      if(.not.associated(alhs)) allocate(alhs(nalhs)); 
!       if (label .eq. 'factor' .or. label .eq. 'full') alhs=0.0;
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
                                     numnp, numel, nen, nsd, ndof, nalhs, neq, label )
!

      use mGlobaisEscalares, only: nrowsh_bm, npint_bm	! OK
      use mGlobaisArranjos,  only: mat_bm, c_bm, grav_bm, celast	! OK
      use mAlgMatricial,     only: kdbc, addrhs, addlhs	! OK
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shlq3d, shg3d, shgq, shlq	! OK
      use mMalha,            only: local	! OK
      use mBlocoMacro,       only: solucao_BM, ndof_bm
      use mParametros,       only: p_Ref, tamBlocoMacro, widthBlocoMacro, alpha_r
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
      character(LEN=12), intent(in) :: label
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
      real*8  :: djx, djy, djz, djn, dix, diy, diz, p, gpx, gpy
      logical :: diag,zerodl,quad,lsym
!       Variaveis para o estado plano de tensoes (TODO: Automatizar leitura)
      real*8  :: D1, CMATRIX11, CMATRIX12, CMATRIX33, YOUNG, POISSON, RHOMAT
!       *******************************************
      real*8  :: Lx, Ly
	
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
      p = 0.0d0
      gpx = 0.0d0
      gpy = 0.0d0
      do 301 k=1,nen
      p = p + shg(3,k,l)*gpl(1,k)
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
      ELRESF(nj-1)= ELRESF(nj-1) + alpha_r*p*djx!+ RHOMAT*GRAV_BM(1)*djn
      ELRESF(nj)  = ELRESF(nj) + alpha_r*p*djy!+ RHOMAT*GRAV_BM(2)*djn 
!       
      do 300 i=1,nen
      if (label .eq. 'full' .or. label .eq. 'factor') then
      ni = ndof*i
                 dix=shg(1,i,l)
                 diy=shg(2,i,l) 
      if(nsd==3) diz=shg(3,i,l)       
      
      ELEFFM(ni-1,nj-1)=ELEFFM(ni-1,nj-1)+(DIX)*(CMATRIX11*DJX)+(DIY)*(CMATRIX33*DJY)
      ELEFFM(ni-1,nj)=ELEFFM(ni-1,nj)+(DIX)*(CMATRIX12*(DJY))+(DIY)*(CMATRIX33*(DJX))
      ELEFFM(ni,nj-1)=ELEFFM(ni,nj-1)+(DIY)*(CMATRIX12*(DJX))+(DIX)*(CMATRIX33*(DJY))
      ELEFFM(ni,nj)=ELEFFM(ni,nj)+(DIY)*(CMATRIX11*DJY)+(DIX)*(CMATRIX33*DJX)
      
      else 
	continue
      endif
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
      if (label .eq. 'full' .or. label .eq. 'factor') call addlhs(alhs,eleffm,idiag,lm(1,1,nel),nee,diag,lsym) 
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
    
!
!**** new **********************************************************************
!
      subroutine calcCoefSistAlgTerzaghi(d, x, alhs, brhs, idiag, lm, conecNodaisElem, &
                                     numnp, numel, nen, nsd, ndof, nalhs, neq, label )
!

      use mGlobaisEscalares, only: nrowsh_bm, npint_bm	! OK
      use mGlobaisArranjos,  only: mat_bm, c_bm, grav_bm, celast	! OK
      use mAlgMatricial,     only: kdbc, addrhs, addlhs	! OK
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shlq3d, shg3d, shgq, shlq	! OK
      use mMalha,            only: local	! OK
      use mBlocoMacro,       only: solucao_BM, ndof_bm
      use mParametros,       only: p_Ref, tamBlocoMacro, widthBlocoMacro
      use mParametros,       only: k_s, constMu
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
      character(LEN=12), intent(in) :: label
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
      real*8  :: djx, djy, djz, djn, dix, diy, diz, p, gpx, gpy
      logical :: diag,zerodl,quad,lsym
!       Variaveis para o estado plano de tensoes (TODO: Automatizar leitura)
      real*8  :: D1, CMATRIX11, CMATRIX12, CMATRIX33, YOUNG, POISSON, RHOMAT
!       *******************************************
      real*8  :: Lx, Ly
      real*8  :: alpha_r, Kbulk, lambda, mu, omega_T
	
      nee = nen*ndof
      
      diag = .false.
      pi=4.d00*datan(1.d00)
      YOUNG = celast(1)
      POISSON = celast(2)
      RHOMAT = celast(3)
!
      lambda = (celast(1)*celast(2))/((1.0+celast(2))*(1.0-2.0*celast(2)))
      mu = (celast(1))/(2.0*(1.0+celast(2)))
      Kbulk = (lambda + 2.0/3.0*mu)
      alpha_r = 1.0d0 - Kbulk/k_s
      omega_T = 1.0d0

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
      p = 0.0d0
      gpx = 0.0d0
      gpy = 0.0d0
      do 301 k=1,nen
      p = p + shg(3,k,l)*gpl(1,k)
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
      ELRESF(nj-1)= ELRESF(nj-1) &
               &+ ((alpha_r*p-omega_T)/(Kbulk+4.0/3.0*mu))*djn!+ RHOMAT*GRAV_BM(1)*djn
      ELRESF(nj)  = ELRESF(nj) + ((alpha_r*p-omega_T)/(Kbulk+4.0/3.0*mu))*djn!+ RHOMAT*GRAV_BM(2)*djn 
!       
      do 300 i=1,nen
      if (label .eq. 'full' .or. label .eq. 'factor') then
      ni = ndof*i
                 dix=shg(1,i,l)
                 diy=shg(2,i,l) 
      if(nsd==3) diz=shg(3,i,l)       
      
      ELEFFM(ni-1,nj-1)=ELEFFM(ni-1,nj-1)+(DIX)*(DJX)+(DIY)*(DJY)
      ELEFFM(ni-1,nj)=ELEFFM(ni-1,nj)+(DIX)*(DJY)+(DIY)*(DJX)
      ELEFFM(ni,nj-1)=ELEFFM(ni,nj-1)+(DIY)*(DJX)+(DIX)*(DJY)
      ELEFFM(ni,nj)=ELEFFM(ni,nj)+(DIY)*(DJY)+(DIX)*(DJX)
      
      else 
	continue
      endif
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
      if (label .eq. 'full' .or. label .eq. 'factor') call addlhs(alhs,eleffm,idiag,lm(1,1,nel),nee,diag,lsym) 
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
      end subroutine calcCoefSistAlgTerzaghi
    
      

!**** new **********************************************************************
!
!     Essa sub-rotina calcula os sigmas por elemento considerando o estado
!     plano de deformacoes. Funciona apenas para o caso bidimensional.
!     Diego (Set/2015)
!
      subroutine calcStressPoro(stress, d, p, gp, x,conecNodaisElem, &
                            numnp, numel, nen, nsd, ndof, tflag)
!

      use mGlobaisEscalares, only: nrowsh_bm, npint_bm	! OK
      use mGlobaisArranjos,  only: mat_bm, c_bm, grav_bm, celast, phi_n, phi_n0	! OK
      use mGlobaisArranjos,  only: trEps
      use mAlgMatricial,     only: kdbc, addrhs, addlhs	! OK
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shlq3d, shg3d, shgq, shlq	! OK
      use mMalha,            only: local	! OK
      use mMalha,            only: nelx_BM, nely_BM, nen_bm
      use mBlocoMacro,       only: solucao_BM, ndof_bm
      use mParametros,       only: p_Ref, tamBlocoMacro, widthBlocoMacro, alpha_r
      use mParametros,       only: k_s, phi_BM
      !use mMalha,         only:  conecNodaisElem
!       use mPardisoSistema, only: addlhsCSR, ApElasticidade, AiElasticidade
!       use mHYPRESistema,   only: addnslHYPRE, addrhsHYPRE
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!       
      integer*4, intent(in)  :: numnp, numel, nen, nsd, ndof
      logical,   intent(in)  :: tflag 
!
      real*8,  intent(inout) :: stress(nen,numel)!, phi(nel)
      real*8,  intent(in)    :: d(ndof, numnp), x(nsd, numnp),p(ndof, numnp)!, p_ant(ndof, numnp)
!       integer*4, intent(in)    :: idiag(neqD),    lm(ndof,nen,numel) 
      integer*4, intent(in)    :: conecNodaisElem(nen,numel)
!
      real*8 :: xl(nsd,nen), dl(ndof,nen), pl(1,nen), gp(nsd,numnp), gpl(ndof,nen), gp_mean(nsd,numel)
      real*8 :: shg(nrowsh_bm,nen,npint_bm), shl(nrowsh_bm,nen,npint_bm)
      real*8 :: det(npint_bm), w(npint_bm)
!
      integer*4 :: nee, idsx, idx, idsr
      real*8   :: strain(nen), p_mean(numel), p_mean_ant, residbc(nelx_BM)
      real*8   :: denStress_x, denStress_y, valueError_x, valueError_y, denGrad_x, denGrad_y
      REAL*8 :: residS_x, residS_y, sig_down,sig_up,sig_left,sig_right, valueError(nsd,numel)
      CHARACTER*30  idsStr, sigma, rsigma
!
      integer*4:: nel, m, l, i, j, k, ni, nj, n, N_ydown,N_yup,N_xleft,N_xright
      integer*4, parameter :: um = 1
!       real*8  :: pi, Kx, Ky, Kz
      real*8  :: temp1, Area!, gf1, gf2, gf3
!       real*8  :: pss
      real*8  :: djx, djy, djz, djn, gpx, gpy, p_mean_tmp, sigma_bc, divu
      logical :: diag,zerodl,quad,lsym
      real*8  :: D1, CMATRIX11, CMATRIX12, CMATRIX33, YOUNG, POISSON, RHOMAT
!       *******************************************
      real*8  :: Lx, Ly
!
      nee = nen*ndof

      diag = .false.
      YOUNG = celast(1)
      POISSON = celast(2)
      RHOMAT = celast(3)
      w=0.0
      shl=0.0

      call   shlq(shl,w,npint_bm,nen)
	
      do 500 nel=1,numel
!
!      LOCALIZE COORDINATes and Dirichlet b.c.
!
      call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
      call local(conecNodaisElem(1,nel),d,dl,nen,ndof,ndof)
      call local(conecNodaisElem(1,nel),p,pl,nen,ndof_bm,ndof_bm)
      call local(conecNodaisElem(1,nel),gp,gpl,nen,2*ndof_bm,2*ndof_bm)
!
      quad = .true.
      if (nen.eq.4.and.conecNodaisElem(3,nel).eq.conecNodaisElem(4,nel)) quad = .false.
!
      call shgq  (xl,det,shl,shg,npint_bm,nel,quad,nen)
!
!.... SETUP ELASTICITY TENSOR COEFFICIENTS: PLANE STRAIN 
!
      D1=YOUNG/((1.0d0+POISSON)*(1.0d0-2.0d0*POISSON))
      CMATRIX11=D1*(1.0d0-POISSON)
      CMATRIX12=D1*POISSON
      CMATRIX33=D1*(1.0d0-2.0d0*POISSON)*0.5d0
!       
!.... Clear initial strain and grad p
!
      strain = 0.0d0
!       
!.... Area element
!
      Area = 0.0d0
!
!.... loop on integration points
!
      do 400 l=1,npint_bm

      temp1 = w(l)*det(l)
      Area = Area + temp1
      DIVU=0.d0
      DO J=1,NEN
         DIVU=DIVU+SHG(1,J,L)*DL(1,J)+SHG(2,J,L)*DL(2,J)
      ENDDO

      do 300 j=1,nen
!            nj=ndof*j
                 djx=shg(1,j,l)*temp1
                 djy=shg(2,j,l)*temp1
      if(nsd==3) djz=shg(3,j,l)*temp1
      djn=shg(nrowsh_bm,j,l)*temp1
!     
!.... compute strains     
!
      STRAIN(1) = STRAIN(1)+DJX*DL(1,J)
      STRAIN(2) = STRAIN(2)+DJY*DL(2,J)
      STRAIN(3) = STRAIN(3)+DJX*DL(2,J)+DJY*DL(1,J)
!
  300 continue
  400 continue

!
!.... compute mean strain over element
!
      strain = strain/area
! 
!.... compute strain's trace
! 
      trEps(nel) = (strain(1) + strain(2))
!
!.... compute mean pressure over element
! 
      p_mean(nel) = sum(pl)/nen!*p_Ref
      gp_mean(1,nel) = sum(gpl(1,:))/nen!*p_Ref
      gp_mean(2,nel) = sum(gpl(2,:))/nen!*p_Ref
!     
!.... compute mean porosity over element
!
      phi_n(nel) = phi_n0(nel) + alpha_r*divu + &
      &            ((alpha_r-phi_n0(nel))/k_s)*(p_mean(nel) - p_Ref)
      !write(*,*) nel, divu
!
!.... compute element volumetric stresses 
!
      STRESS(1,NEL)=CMATRIX11*STRAIN(1)+CMATRIX12*STRAIN(2)
      STRESS(2,NEL)=CMATRIX12*STRAIN(1)+CMATRIX11*STRAIN(2)
      STRESS(3,NEL)=CMATRIX33*STRAIN(3)
      
  500 continue  

      return
      end subroutine calcStressPoro
      

!**** new **********************************************************************
!
!     Essa sub-rotina imprime os sigmas por elemento considerando o estado
!     plano de tensoes. Funciona apenas para o caso bidimensional.
!     Diego (Dez/2015)
!
      subroutine printStressPoro(stress, d, p, gp, x,conecNodaisElem, &
                            numnp, numel, nen, nsd, ndof, tflag, idx)
!

      use mGlobaisEscalares, only: nrowsh_bm, npint_bm	! OK
      use mGlobaisArranjos,  only: mat_bm, c_bm, grav_bm, celast, phi_n, phi_n0	! OK
      use mAlgMatricial,     only: kdbc, addrhs, addlhs	! OK
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shlq3d, shg3d, shgq, shlq	! OK
      use mMalha,            only: local	! OK
      use mMalha,            only: nelx_BM, nely_BM, nen_bm
      use mBlocoMacro,       only: solucao_BM, ndof_bm
      use mParametros,       only: p_Ref, tamBlocoMacro, widthBlocoMacro, alpha_r
      use mParametros,       only: k_s, phi_BM
      !use mMalha,         only:  conecNodaisElem
!       use mPardisoSistema, only: addlhsCSR, ApElasticidade, AiElasticidade
!       use mHYPRESistema,   only: addnslHYPRE, addrhsHYPRE
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!       
      integer*4, intent(in)  :: numnp, numel, nen, nsd, ndof
      logical,   intent(in)  :: tflag 
!
      real*8,  intent(inout) :: stress(nen,numel)!, phi(nel)
      real*8,  intent(in)    :: d(ndof, numnp), x(nsd, numnp),p(ndof, numnp)!, p_ant(ndof, numnp)
!       integer*4, intent(in)    :: idiag(neqD),    lm(ndof,nen,numel) 
      integer*4, intent(in)    :: conecNodaisElem(nen,numel)
!
      real*8 :: xl(nsd,nen), dl(ndof,nen), pl(1,nen), gp(nsd,numnp), gpl(ndof,nen), gp_mean(nsd,numel)
      real*8 :: shg(nrowsh_bm,nen,npint_bm), shl(nrowsh_bm,nen,npint_bm)
      real*8 :: det(npint_bm), w(npint_bm)
!
      integer*4 :: nee, idsx, idx, idsr
      real*8   :: strain(nen), p_mean(numel), p_mean_ant, trEps, residbc(nelx_BM)
      real*8   :: denStress_x, denStress_y, valueError_x, valueError_y, denGrad_x, denGrad_y
      REAL*8 :: residS_x, residS_y, sig_down,sig_up,sig_left,sig_right, valueError(nsd,numel)
      CHARACTER*30  idsStr, sigma, rsigma
!
      integer*4:: nel, m, l, i, j, k, ni, nj, n, N_ydown,N_yup,N_xleft,N_xright
      integer*4, parameter :: um = 1
!       real*8  :: pi, Kx, Ky, Kz
      real*8  :: temp1, Area!, gf1, gf2, gf3
!       real*8  :: pss
      real*8  :: djx, djy, djz, djn, gpx, gpy, p_mean_tmp, sigma_bc, divu(numel)!, gpl
      logical :: diag,zerodl,quad,lsym
      real*8  :: D1, CMATRIX11, CMATRIX12, CMATRIX33, YOUNG, POISSON, RHOMAT
!       *******************************************
      real*8  :: Lx, Ly
!	
      w=0.0
      shl=0.0
      call   shlq(shl,w,npint_bm,nen)
      divu = 0.0

  do 500 nel=1,numel
!
!      LOCALIZE COORDINATes and Dirichlet b.c.
!
      call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
      call local(conecNodaisElem(1,nel),p,pl,nen,ndof_bm,ndof_bm)
      call local(conecNodaisElem(1,nel),gp,gpl,nen,2*ndof_bm,2*ndof_bm)

      quad = .true.
      if (nen.eq.4.and.conecNodaisElem(3,nel).eq.conecNodaisElem(4,nel)) quad = .false.
!
      call shgq  (xl,det,shl,shg,npint_bm,nel,quad,nen)
!
!.... loop on integration points
!
      do 400 l=1,npint_bm

      !DIVU=0.d0
      DO J=1,NEN
         DIVU(nel)=DIVU(nel)+SHG(1,J,L)*DL(1,J)+SHG(2,J,L)*DL(2,J)
      ENDDO
      !write(*,*) nel, divu(nel)
  400 continue
      !write(*,*) nel, divu(nel)
!
!.... compute mean pressure over element
! 
      p_mean(nel) = sum(pl)/nen!*p_Ref
      gp_mean(1,nel) = sum(gpl(1,:))/nen!*p_Ref
      gp_mean(2,nel) = sum(gpl(2,:))/nen!*p_Ref
      
  500 continue  

  if (tflag .eqv. .true.) then
        !       Abrindo os arquivos para saída
	write(idsStr,'(i0)') idx

	idsx = 11*idx
	sigma = 'sigy.'//idsStr
	OPEN(UNIT=idsx, FILE= sigma)
	
	idsr = 12*idx
	rsigma = 'eqdiv.'//idsStr
	OPEN(UNIT=idsr, FILE= rsigma)
	
	OPEN(UNIT=(13*idx), FILE= 'eqdivL.'//idsStr)
	OPEN(UNIT=(17*idx), FILE= 'eqdivC.'//idsStr)
	OPEN(UNIT=(19*idx), FILE= 'eqdivR.'//idsStr)
	OPEN(UNIT=(201), FILE= 'sigmasL.'//idsStr)
	OPEN(UNIT=(202), FILE= 'sigmasR.'//idsStr)
	OPEN(UNIT=(203), FILE= 'sigmasC.'//idsStr)
	OPEN(UNIT=(301), FILE= 'sigmasLx.'//idsStr)
	OPEN(UNIT=(302), FILE= 'sigmasRx.'//idsStr)
	OPEN(UNIT=(303), FILE= 'sigmasCx.'//idsStr)
	OPEN(UNIT=(997), FILE= 'divu.'//idsStr)
	residbc = 0.0
	DO J=1,nelx_BM
	!             Element localization
		N = J+(nely_BM-1)*(nely_BM)
	        write(idsx, 223) N, (STRESS(2,N) - alpha_r*p_mean(N))
	ENDDO
!       Computation check of element's momentum balance
      DO I=1,nely_BM - 2
          DO J=2,nelx_BM - 1
!             Element localization
            N = J+(I)*(nely_BM)
            N_xleft = N - 1 
            N_xright = N + 1
            N_ydown = J+(I-2)*(nely_BM + 1)
            N_yup = J+(I)*(nely_BM + 1)

!             Stress N,S,E,W terms
            sig_left = (STRESS(1,N_xleft)+(STRESS(3,N_xleft)))
            sig_right = (STRESS(1,N_xright)+(STRESS(3,N_xright)))
            sig_down = (STRESS(2,N_ydown)+(STRESS(3,N_ydown)))
            sig_up = (STRESS(2,N_yup)+(STRESS(3,N_yup)))
            
!             Stress finite difference computation at element
            residS_x = (sig_right - sig_left)/((X(1,N_xright)-X(1,N_xleft)))
            residS_y = (sig_up - sig_down)/((X(2,N_yup)-X(2,N_ydown)))
            
	    denStress_x = dabs(STRESS(1,N))
	    denStress_y = dabs(STRESS(2,N))
	    denGrad_x = dabs(alpha_r*gp_mean(1,N))
	    denGrad_y = dabs(alpha_r*gp_mean(2,N))
	    valueError(1,N) = dabs(residS_x - alpha_r*gp_mean(1,N))/denStress_x
	    valueError(2,N) = dabs(residS_y - alpha_r*gp_mean(2,N))/denStress_y
            write(idsr,224) N, valueError(1,N), valueError(2,N) 
          ENDDO
      ENDDO
      DO I=1,nely_BM
          DO J=1,nelx_BM
            N = J+(I-1)*(nely_BM)
            !write(*,*) "Aqui 1"
            write(997,226) N, X(1,N), X(2,N), DIVU(N)
            !write(*,*) N, "div u = ", divu(n)
          ENDDO
      ENDDO
      
      DO I=2,nelx_BM, (nelx_BM - 3)/2
          DO J=2,nely_BM - 1
            N = I+(J-1)*(nelx_BM)
            if (I .eq. 2) then
		WRITE((13*idx),225) N, X(1,N), X(2,N), valueError(1,N), valueError(2,N)
	    else if (I .eq. (nelx_BM-2)) then
	        WRITE((19*idx),225) N, X(1,N), X(2,N), valueError(1,N), valueError(2,N)
	    else if (I .eq. ((nelx_BM)/2)) then
	        WRITE((17*idx),225) N, X(1,N), X(2,N), valueError(1,N), valueError(2,N)
            endif
          ENDDO
      ENDDO  
      DO I=1,nelx_BM, (nelx_BM - 1)/2
          DO J=1,nely_BM
            N = I+(J-1)*(nelx_BM)
            if (I .eq. 1) then
		write(201,226) N, STRESS(1,N), STRESS(2,N), STRESS(3,N)
	    else if (I .eq. (nelx_BM-1)) then
	        write(202,226) N, STRESS(1,N), STRESS(2,N), STRESS(3,N)
	    else if (I .eq. ((nelx_BM)/2)) then
	        write(203,226) N, STRESS(1,N), STRESS(2,N), STRESS(3,N)
            endif
          ENDDO
      ENDDO
      DO I=1,nely_BM, (nely_BM - 1)/2
          DO J=1,nelx_BM
            N = J+(I-1)*(nely_BM)
            if (I .eq. 1) then
		write(301,226) N, STRESS(1,N), STRESS(2,N), STRESS(3,N)
	    else if (I .eq. (nely_BM-1)) then
	        write(302,226) N, STRESS(1,N), STRESS(2,N), STRESS(3,N)
	    else if (I .eq. ((nely_BM)/2)) then
	        write(303,226) N, STRESS(1,N), STRESS(2,N), STRESS(3,N)
            endif
          ENDDO
      ENDDO
  endif
  
  close(idsx)
  close(idsr)
  close(995)
  close(201)
  close(202)
  close(203)
  close((13*idx))
  close((17*idx))
  close((19*idx))
  close((1*idx))
 
  223  FORMAT(4X,I5,10x,1(1PE15.8,2X))
  224  FORMAT(4X,I5,10x,2(1PE15.8,2X))
  225  FORMAT(4X,I5,10x,4(1PE15.8,2X))
  226  FORMAT(4X,I5,10x,3(1PE15.8,2X))
      return
      end subroutine printStressPoro
      
!       
!     *********************************************************  
!       
      subroutine alocarMemoriaElasticidade()
      
!      Alocacao de memoria da Elasticidade Linear. Diego (Ago/2015)

     use mGlobaisEscalares, only: ndofD, nlvectD
     use mMalha,            only: numnp_BM

     implicit none
     
     allocate(  deslocamento(ndofD,NUMNP_BM));   deslocamento=0.0
     allocate(  deslocamentoAnt(ndofD,NUMNP_BM));   deslocamentoAnt=0.0
     allocate(idDeslocamento(ndofD,NUMNP_BM));  idDeslocamento=0
     
     if (nlvectD.ne.0) then 
        allocate(fDeslocamento(nlvectD,ndofD,numnp_BM))
        fDeslocamento = 0.0
     endif
     
     end subroutine alocarMemoriaElasticidade

     end module
!
