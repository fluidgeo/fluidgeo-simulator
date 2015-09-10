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
      module mLeituraEscrita
!
      implicit none
!
      integer :: contadorUnidadeLogica = 10 
	  
      integer :: iin,iecho,icoords,iconects,iconectsL
      integer :: iechoB, iechoF, iechoQuantidadeGas
      integer :: ignuplot,iparaviewS,iparaviewP,iparaviewV
      integer :: isatTransiente
      integer :: nprint
      integer :: qtdImpSat

      integer :: ipres  = 20
      integer :: isat   = 21
      integer :: iphi   = 22
      integer :: ivel   = 23
      integer :: imass  = 24
      integer :: ifiles = 25
      integer :: imasl  = 27
      integer :: IMASG  = 28
      integer :: ISIGT  = 29
      integer :: ifdata = 150
      integer :: ifrand = 151
      integer :: ipress = 155
      integer :: iperm  = 156
      integer :: iporo  = 157
      integer :: iveloc = 158
      integer :: isaturacao = 159
      INTEGER :: IFEDX = 160
      INTEGER :: IFNOD = 203
      INTEGER :: IFMSH = 204
      INTEGER :: IRBAL = 205

      integer :: iflag_pres,iflag_sat,iflag_phi,iflag_vel
      integer :: iflag_mass, iflag_masl, IFLAG_MASG, IFLAG_SIGT
      integer :: iflag_tipoPrint

!
      character(len=128) :: ifpres_out
      character(len=128) :: ifsat_out
      character(len=128) :: ifphi_out
      character(len=128) :: ifperm_out
      character(len=128) :: ifvel_out
      character(len=128) :: ifmass_out
      character(len=128) :: ifmasl_out
      character(len=128) :: IFMASG_OUT
      character(len=128) :: IFSIGT_OUT


      real(8) :: tprt_pres, dtprt_pres
      real(8) :: tprt_sat , dtprt_sat
      real(8) :: tprt_phi , dtprt_phi
      real(8) :: tprt_vel , dtprt_vel
      real(8) :: tprt_mass, dtprt_mass
      real(8) :: tprt_masl, dtprt_masl
      real(8) :: TPRT_MASG, DTPRT_MASG
      real(8) :: TPRT_SIGT, DTPRT_SIGT
      integer :: nppres, npsat, npphi
      integer :: npvel, npmass, npmasl, NPMASG, NPSIGT

!
      contains
!
!=================================================================================
!
     function abrirArquivo(nomeArquivo, unidLogica, estado) result (novaUnidLogica)
      character(len=*)           :: nomeArquivo
	  	  integer         , optional :: unidLogica 
      character(len=*), optional :: estado
	 	  	  integer   ::  novaUnidLogica

    if (present(unidLogica)) then
        if(unidLogica==contadorUnidadeLogica) unidLogica = unidLogica+1
    else
        contadorUnidadeLogica = contadorUnidadeLogica + 10
	unidLogica = contadorUnidadeLogica 
    end if
	
    if (.not.present(estado)) then
	   estado='old'
    endif
	
     if(unidLogica==contadorUnidadeLogica) contadorUnidadeLogica = contadorUnidadeLogica + 10
     write(*,*) 'Abertura do arquivo ',  nomeArquivo,', ',  estado, ' em ', unidLogica
     novaUnidLogica=unidLogica
     write(*,*) "ANtes"
     open(unit=unidLogica, file= nomeArquivo, status=estado)
     write(*,*) "Após"
     
     end function 
   
    subroutine abrirArquivosInformacoesMalha()
!  
!        iin    = input unit number
!        iecho  = output unit of input data
!        iouter  = output unit of error norms
!
      character(len=2) :: nRank
      character(len=20) :: nomeIecho
      character(len=20) :: nomeIechoB, nomeIechoF
!
      iin           = 15
      iecho         = 116
      iechoB        = 117
      iechoF        = 118
      icoords       = 18
      iconects      = 19
      iconectsL     = 25
      isatTransiente= 43
!
!      nomeIechoB='echoB.dat'
!      open(unit=iechoB , file= nomeIechoB)
!      nomeIechoF='echoF.dat'
!      open(unit=iechoF , file= nomeIechoF)
      
!
#ifdef debug
      open(unit=icoords    , file= 'coordenadas.dat')
      open(unit=iconects   , file= 'conectsNodais.dat')
      open(unit=iconectsL  , file= 'conectsLadais.dat')
#endif

   end subroutine 

!
!**** new **********************************************************************
    subroutine fecharArquivos()
!
      close(iin   )
      close(iechoB )
      close(iechoF )
      close(iechoQuantidadeGas )
      close(icoords )
      close(ignuplot)
      close(iconectsL)
      close(iconects )
      close(isat)
      close(ipres)
      close(ivel)
      close(iphi)
      close(iperm)
      close(imass)
      close(imasl)
      close(IMASG)
 
    end subroutine fecharArquivos
!
!*****************************************************
!
! !      subroutine imprimirCondicoesIniciais(pressaoElem, velocLadal, phi, perm, satElem)
! !        use mGlobaisEscalares, only : ndof_F, ndofV, tTransporte, numdx
! !        use mMalha,            only : x, conecNodaisElem, conecLadaisElem, nen, nsd, numel, numLados, numLadosElem, numnp
! ! 
! ! 
! !       implicit none
! ! 
! !       real*8, intent(in) :: pressaoElem(ndofP, numel), velocLadal(ndofV,numLados)
! !       real*8, intent(in) :: phi(numel), perm(numel), satElem(numel)
! !       integer :: i, k
! ! !
! ! !.... imprime a condicao inicial: pressao
! ! !
! !        if(iflag_pres==1) then
! !           if(iflag_tipoPrint==0) then
! !              call prt(nsd,numel,tTransporte,pressaoElem,ipres)
! !           else
! !              call escreverArqParaview(ipres, pressaoElem, ndofP, numel, nen, conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
! !           endif
! !        endif
! ! !  
! ! !.... imprime a condicao inicial: velocidade
! ! !
! !       if(iflag_vel==1) then
! !          if(iflag_tipoPrint==0) then
! !             call prtvB(nsd,numel,tTransporte,velocLadal,ndofV, conecLadaisElem, numLadosElem,ivel) 
! !          else
! !             write(ivel,*) "impressao da velocidade ladal nao implementada para o paraview"
! !          endif
! !       end if
! ! !     
! ! !.... imprime a condicao inicial: saturacao
! ! !
! !       if(iflag_sat==1) then    
! !          if(iflag_tipoPrint==0) then 
! !             call prt (nsd,numel,tTransporte,satElem,isat)
! !          else
! !             call escreverArqParaview(isat, satElem, ndofV, numel, nen, conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
! !          endif
! !       endif
! ! 
! ! !
! ! !.... imprime a condicao inicial: porosidade
! ! !
! !        if(iflag_phi==1) then
! !           if(iflag_tipoPrint==0) then
! !              call prt(nsd,numel,tTransporte,phi,iphi)
! !           else
! !              call escreverArqParaview(iphi, phi, ndofV, numel, nen, conecNodaisElem, 1, 'porosidade', len('porosidade')) 
! !           endif
! !        endif
! ! 
! ! !
! ! !.... imprime a condicao inicial: permeabilidade
! ! !
! !        if(iflag_phi==1) then
! !           if(iflag_tipoPrint==0) then
! !              call prt(nsd,numel,tTransporte,perm,iperm)
! !           else
! !              call escreverArqParaview(iperm, perm, ndofV, numel, nen, conecNodaisElem, 1, 'permeabilidade', len('permeabilidade')) 
! !           endif
! !        endif
! ! 
! ! 
! ! !
! !      IF (NUMDX.EQ.0) RETURN
! ! ! 
! ! ! .... PRINT DISPLACEMENTS DATA FOR OPEN-DX FILE
! ! ! 
! !      DO 10 I=1,NUMNP 
! ! 	 WRITE(IFNOD,1900) (X(K,I),K=1,NSD)  
! ! 10   CONTINUE 
! ! ! 
! ! ! .... PRINT CONECTIVITIES DATA FOR OPEN-DX FILE
! ! ! 
! !      DO 20 I=1,NUMEL 
! !         WRITE(IFMSH,2000) conecNodaisElem(1,I)-1, conecNodaisElem(2,I)-1, &
! !     &                     conecNodaisElem(4,I)-1, conecNodaisElem(3,I)-1 
! ! 20   CONTINUE 
! ! !                     
! ! ! .... OPEN-DX INFORMATION DATA FILES
! ! ! 
! ! ! .... NODAL POINTS OF FINITE ELEMENT DISCRETIZATION
! !       IF (NUMDX.GT.0) CALL PRINT_DXINFO('OPEN__FEDX_FILE',NUMNP,NUMEL)
! ! !
! !       RETURN
! ! !
! ! 1900  FORMAT(6(1PE11.4,2X)) 
! ! 2000  FORMAT(27I6) 
! ! ! 
! ! 
! ! 
! !        end subroutine
!
!*****************************************************
!
! !      subroutine imprimirSolucaoNoTempo(pressaoElem, velocLadal, tempo)
! !        use mGlobaisEscalares, only : ndofP, ndofV
! !        use mMalha,            only : nsd, numel, numLados, numLadosElem, conecLadaisElem
! ! 
! !       implicit none
! ! 
! !       real*8, intent(in) :: pressaoElem(ndofP,numel), velocLadal(ndofV,numLados)
! !       real*8, intent(in) :: tempo
! ! !
! !       character(21) :: label
! !       character(21) :: num
! !    
! ! ! 
! ! !.... imprime a velocidade no centro
! ! !
! !       if(iflag_vel==1) then
! !          if(iflag_tipoPrint==0) then
! !             call prtvB(nsd,numel,tempo,velocLadal,ndofV, conecLadaisElem, numLadosElem,ivel) 
! !          endif
! !       end if
! ! 
! ! !
! ! !.... imprime a pressao no tempo 
! ! !
! !       if(iflag_pres==1) then
! !          if(iflag_tipoPrint==0) then
! !             call prt(nsd,numel,tempo,pressaoElem,ipres)
! !          else
! !             write(num,'(f12.5)') tempo
! !             label="t="//ADJUSTL(num)
! !             call escreverArqParaviewIntermed(ipres, pressaoElem, ndofP, numel, trim(label), len(trim(label)))
! !          endif
! !       endif
! ! 
! !       end subroutine 
! ! 

!**** new **********************************************************************
      subroutine echo(iin)

      implicit none
!
!.... program to echo input data
!
      character*4 ia(20)
      integer, intent(in)  :: iin
      integer :: iech, i
!      common /iounit/ iin,iout,iecho,iouter,ignuplot
!
!     cabeçalho
      write(iecho,500)

      read(iin,1000) iech
      if (iech.eq.0) return
!
      write(iecho,2000) iech
      backspace iin
!
      do 100 i=1,100000
      read(iin,3000,end=200) ia
      if (mod(i,50).eq.1) write(iecho,4000)
      write(iecho,5000) ia

  100 continue
!
  200 continue
      rewind iin
      read(iin,1000) iech
!
      return
!
 500  format('programa de elementos finitos em fortran 90 baseado em:',// &
      'The Finite Element Method, Hughes, T. J. R., (2003)'//)
 1000 format(16i10)
 2000 format('1',' i n p u t   d a t a   f i l e               ',  //5x,&
     ' echo print code . . . . . . . . . . . . . . (iecho ) = ',i10//5x,&
     '    eq. 0, no echo of input data                        ',   /5x,&
     '    eq. 1, echo input data                              ',   ///)
 3000 format(20a4)
 4000 format(' ',8('123456789*'),//)
 5000 format(' ',20a4)

      end subroutine echo

!******************************************************************************
      subroutine genel1(conecElem,mat,nen, n,nel,incel,inc) 
!                                                                       
!.... program to generate element node and material numbers             
!                                                                       
      integer*4 ::  nen, conecElem(nen,*),mat(*)                                       
      integer*4:: n,nel(3),incel(3),inc(3)                          

      integer*4:: i,j,k,ii,jj,kk,ie,le,je,ke
!                                                                       
!.... set defaults                                                      
!                                                                       
      call geneld                                                 
!                                                                       
!.... generation algorithm                                              
!                                                                       
      ie = n                                                            
      je = n                                                            
      ke = n                                                            
!                                                                      
      ii = nel(1)                                                       
      jj = nel(2)                                                       
      kk = nel(3)                                                       
!                                                                       
      do 300 k=1,kk                                                     
!
      do 200 j=1,jj                                                     
!
      do 100 i=1,ii                                                     
!                                                                       
      if (i.ne.ii) then
         le = ie                                                           
         ie = le + incel(1)                                                
         call geneli(conecElem(1,ie),conecElem(1,le),inc(1),nen)                       
         mat(ie) = mat(le)                                                 
      endif
  100 continue                                                          
!                                                                       
      if (j.ne.jj) then
         le = je                                                           
         je = le + incel(2)                                                
         call geneli(conecElem(1,je),conecElem(1,le),inc(2),nen)                       
         mat(je) = mat(le)                                                 
         ie = je                                                           
      endif
  200 continue                                                          
!                                                                       
      if (k.ne.kk) then
         le = ke                                                           
         ke = le + incel(3)                                                
         call geneli(conecElem(1,ke),conecElem(1,le),inc(3),nen)                       
         mat(ke) = mat(le)                                                 
         ie = ke
         je=ke                                                           
      endif
  300 continue                                                          
!                                                                      
      return                                                            
      contains 
!******************************************************************************
      subroutine geneld                                                 
!                                                                       
!.... program to set defaults for element node       
!        and material number generation                              
!                                                                       
      if (nel(1).eq.0) nel(1) = 1                                       
      if (nel(2).eq.0) nel(2) = 1                                       
      if (nel(3).eq.0) nel(3) = 1                                       
!                                                                       
      if (incel(1).eq.0) incel(1) = 1                                   
      if (incel(2).eq.0) incel(2) = nel(1)                              
      if (incel(3).eq.0) incel(3) = nel(1)*nel(2)                       
!                                                                       
      if (inc(1).eq.0) inc(1) = 1                                       
      if (inc(2).eq.0) inc(2) = (1+nel(1))*inc(1)                       
      if (inc(3).eq.0) inc(3) = (1+nel(2))*inc(2)                       
!                                                                       
      return                                                            
      end subroutine
!
!******************************************************************************
!
      subroutine geneli(conecElem2,conecElem1,inc,nen)                              
!                                                                       
!.... program to increment element node numbers                         
!                                                                       
      integer*4 ::  conecElem1(*),conecElem2(*)
      integer*4:: inc, nen

      integer*4:: i
!
      do 100 i=1,nen                                                    
      if (conecElem1(i).eq.0) then
         conecElem2(i) = 0
      else
         conecElem2(i) = conecElem1(i) + inc                         
      endif
  100 continue                                                          
!                                                                       
      return                                                            
      end  subroutine
      end subroutine genel1
      
!**** new **********************************************************************
      subroutine leituraGeracaoCoordenadas(x, nsd, numnp, iin, icoords, iprtin)
      use mMalha, only: genfl
!
!.... program to read, generate and write coordinate data
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8, intent(inout) :: x(nsd,*)
      integer, intent(in)   :: nsd, numnp, iin, icoords, iprtin
!
      integer :: i, n
!      
      call genfl(x,nsd,iin)
!
      if (iprtin.eq.1) return
!
#ifdef debug
         write(icoords,*) "# Coordenadas ", nsd
         do n=1,numnp
            write(icoords,2000) n,(x(i,n),i=1,nsd) 
         end do
#endif
!
      return
!
!  1000 format('1',' n o d a l   c o o r d i n a t e   d a t a '///5x,&
!      ' node no.',3(13x,' x',i1,' ',:)//)
 2000 format(6x,i12,10x,3(1pe15.8,2x))
      end subroutine
!**** new **********************************************************************
      subroutine leituraCodigosCondContorno(id, ndof, numnp, neq, iin, iecho, iprtin)
!
!.... program to read, generate and write boundary condition data
!        and establish equation numbers
!
      use mMalha, only: igen

      integer, intent(in) :: ndof, numnp, iin, iecho, iprtin
      integer:: neq
      integer, intent(inout) :: id(ndof,numnp)
!
      integer :: nn, n, i
      logical pflag
!
      id = 0
      call igen(id,ndof, iin)
!
      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
         pflag = .false.
!
         do 100 i=1,ndof
         if (id(i,n).ne.0) pflag = .true.
  100    continue
!
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
            write(iecho,2000) n,(id(i,n),i=1,ndof)
         endif
  200    continue
      endif
      
!
!.... establish equation numbers
!
      neq = 0
!
      do 400 n=1,numnp
!
      do 300 i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         id(i,n) = 1 - id(i,n)
      endif
!
  300 continue
!
  400 continue
!
      return
!
 1000 format('1',' n o d a l   b o u n d a r y   c o n d i t i o n & 
     &         c o  d e s'/// &
      5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
!
      end subroutine
!
!**** new **********************************************************************
!
      subroutine leituraValoresCondContorno(f,ndof,numnp,j, nlvect,iprtin, luEscrita)
! 	Alterado por Diego para ter compatibilidade com as rotinas do DS.
! 	(Ago/2015)
!
!.... program to read, generate and write nodal input data
!
!        f(ndof,numnp,nlvect) = prescribed forces/kinematic data (j=0)
!                             = nodal body forces(j=1)
!
      use mMalha, only : genfl
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: ndof, numnp, j, nlvect, iprtin, luEscrita
      real*8 f(ndof,numnp,nlvect)

      logical lzero
      integer nlv
      character(len=35) :: rotulo
!
!     call clear(f,nlvect*numnp*ndof)
      f(1:nlvect,1:numnp,1:ndof)=0.0

      do 100 nlv=1,nlvect
      call genfl(f(1,1,nlv),ndof,iin)
!       
      call ztest(f(1,1,nlv),ndof*numnp,lzero)
! 
      if (iprtin.eq.0) then
!
         if (lzero) then
!             if (j.eq.0) write(luEscrita,1000) nlv
            if (j.eq.0) write(iecho,1000) nlv
!             if (j.eq.1) write(luEscrita,2000)
            if (j.eq.1) write(iecho,2000)
         else
            if (j.eq.0) call printf(f,ndof,numnp,nlv)!, luEscrita)
!
            if (j.eq.1) then
               rotulo=" n o d a l  b o d y  f o r c e s"
               call printd (rotulo, f,ndof,numnp,luEscrita)
            end if
!
         endif
      endif
!
  100 continue
!
      return
 1000 format('1'//,' there are no nonzero prescribed forces and ',&
         'kinematic boundary conditions for load vector number ',i10)
 2000 format('1'//,' there are no nonzero nodal body forces')
      end subroutine

!**** new **********************************************************************
      subroutine printf(f,ndof,numnp,nlv)!, luEscrita)
!       Alterado por Diego para ter compatibilidade com as rotinas do DS. (Ago/2015)
!
!.... program to print prescribed force and boundary condition data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer ndof, numnp, nlv!,luEscrita
      real*8 :: f(ndof,numnp,*)
!
      logical lzero
      integer :: nn, n, i
!
      nn = 0
!
      do 100 n=1,numnp
      call ztest(f(1,n,nlv),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
!          if (mod(nn,50).eq.1) write(luEscrita,1000) nlv,(i,i=1,ndof)
         if (mod(nn,50).eq.1) write(iecho,1000) nlv,(i,i=1,ndof)
!          write(luEscrita,2000) n,(f(i,n,nlv),i=1,ndof)
         write(iecho,2000) n,(f(i,n,nlv),i=1,ndof)
      endif
  100 continue
!
      return
!
 1000 format('1',&
     ' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',&
     '  b o u n d a r y   c o n d i t i o n s'//5x,&
     ' load vector number = ',i10///5x,&
     ' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i10,10x,6(1pe15.8,2x))
      end subroutine
      
!**** new **********************************************************************
      subroutine printd(name,dva,ndof,numnp,luEscrita)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: ndof,numnp, luEscrita
      character (LEN=*) ::  name
      real*8 dva(ndof,*)
!
      logical lzero
      integer nn, n, i
!
      nn = 0
!
      do 100 n=1,numnp
      call ztest(dva(1,n),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) &
           write(luEscrita,1000) name,(i,i=1,ndof)
         write(luEscrita,2000) n,(dva(i,n),i=1,ndof)
      endif
  100 continue
!
      return
!
 1000 format('1',11a4//1x,'node',6(11x,'dof',i1)/)
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end subroutine
!
!**** new **********************************************************************
!
      subroutine printp(a,idiag,neq,nsq,luEscrita)
      use mGlobaisEscalares
!
!.... program to print array d after Crout factorization 
!        a = u(transpose) * d * u
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: neq, nsq,luEscrita
      real*8 :: a(*)
      integer :: idiag(*)
!
      integer :: n, i
!
      do 100 n=1,neq
      if (mod(n,50).eq.1) write(luEscrita,1000) nsq
      write(luEscrita,1000)
      i = idiag(n)
      write(luEscrita,2000) n,a(i)
  100 continue
!
      return 
!
 1000 format('1',' array d of factorization',/&
     ' a = u(transpose) * d * u ',                                //5x,&
     ' time sequence number   . . . . . . . . . . . . (nsq) = ',i10//5x)
 2000 format(1x,i10,4x,1pe20.8)
      end subroutine
! !
! !**** new **********************************************************************
! !
!       subroutine prntel(mat,conectElem,nen,numel,tipo, unidLogicaEscrita)
!       implicit none
! !
! !.... program to print data for element with "nen" nodes
! !
! !        note: presently the label formats are limited to
! !              elements with one to nine nodes
! !
!       integer :: nen, numel, unidLogicaEscrita
!       integer :: mat(*),conectElem(nen,*)
!       integer :: tipo
! !
!       integer n, i
! !
!       if(tipo==1) then
!       write(unidLogicaEscrita,*) "# Conectividades nodais"
!       do n=1,numel
!         write(unidLogicaEscrita,2000) n,mat(n),(conectElem(i,n),i=1,nen)
!       end do
!       end if
! 
!       if(tipo==2) then
!       write(unidLogicaEscrita,*) "# Conectividades ladais"
!       do  n=1,numel
!         write(unidLogicaEscrita,3000) n,mat(n),(conectElem(i,n),i=1,nen)
!       end do
!       end if
! !
!       return
! !
!  2000 format(1x,i10,9(2x,i10))
!  3000 format(1x,i10,7(2x,i10))
!       end subroutine


!
!**** new **********************************************************************
!
      subroutine prntel(mat,conectElem,nen,numel,tipo)
      implicit none
!
!.... program to print data for element with "nen" nodes
!
!        note: presently the label formats are limited to
!              elements with one to nine nodes
!
      integer*4:: nen, numel
      integer*4:: mat(*),conectElem(nen,*)
      integer*4:: tipo
!
      integer*4 n, i
!
      if(tipo==1) then
	write(iconects,*) "# Conectividades nodais"
	do n=1,numel
	  write(iconects,2000) n,mat(n),(conectElem(i,n),i=1,nen)
	end do
      end if

      if(tipo==2) then
	write(iconectsL,*) "# Conectividades ladais"
	do  n=1,numel
	  write(iconectsL,3000) n,mat(n),(conectElem(i,n),i=1,nen)
	end do
      end if
!
      return
!
 1000 format(///10x,&
     ' e l e m e n t   d a t a',//1x,&
     ' element   material ',9('node ',i1,1x),/1x,&
     '  number    number',//)
 2000 format(1x,i10,9(2x,i10))
 3000 format(1x,i10,7(2x,i10))
      end subroutine


!**** new **********************************************************************
      subroutine printResultado(dva, ndof, numnp, inicio, fim, ulEscrita)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: ndof, numnp, inicio, fim, ulEscrita
      real*8  :: dva(ndof,numnp)
!
      integer :: n, i
!
      write(ulEscrita,*) "# Solucao"
      do 100 n=inicio,fim
         write(ulEscrita,2000) n,(dva(i,n),i=1,ndof)
         !write(*,*) n,(dva(i,n),i=1,ndof)
  100 continue
!
      return
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end subroutine

!**** new **********************************************************************
      subroutine prtgnup(name,x,dva,nsd,ndof,numnp,ulEscrita)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: nsd, ndof, numnp, ulEscrita
      character*4 name(11)
      real*8 :: x(nsd,*),dva(ndof,*)
!
      integer :: n, j, i
!
      write(ulEscrita,*) name
      do 100 n=1,numnp
         write(ulEscrita,2000) (x(j,n),j=1,nsd), (dva(i,n),i=1,ndof)
  100 continue
!
      return
!
 2000 format(6(1pe13.6,2x))
      end subroutine
!
!**** new **********************************************************************
!
!     subroutine escreverArqParaview(arquivo, campo, dim1, dim2, nen, conectElem, tipo, rotulo, tamRot)
!     use mMalha, only: x, nsd, numel,numnp
! 
!     implicit none
!     integer, intent(in) :: arquivo,dim1, dim2
!     double precision, intent(in) :: campo(dim1, dim2)
!     integer :: nen
!     integer :: conectElem(nen,numel)
!     integer :: tipo  !tipo=1 para elemento, e tipo=2 para no
!     integer :: tamRot
! 
!     character(len=tamRot) :: rotulo
! 
!   
!     write(arquivo,'(a)')'# vtk DataFile Version 3.0'
!     write(arquivo,'(a)')'vtk output'
!     write(arquivo,'(a)')'ASCII'
!     write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
!     write(arquivo,'(a,i10,a)')'POINTS', numnp,' float '
! 
!     call escreverPontosNodais  (arquivo,x, numnp, nsd)
! ! 
!     write(arquivo,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
!     call escreverConectividades(arquivo,conectElem, numel, nen, nsd)
! ! 
!     write(arquivo,'(a,i10)')'CELL_TYPES ', numel
!     call escreverTiposElementos(arquivo,numel,nsd)
! ! 
! 
!     if(tipo==1) write(arquivo,'(a,i10)')'CELL_DATA ', numel
! 
!     if(tipo==2) write(arquivo,'(a,i10)')'POINT_DATA',  dim1*dim2
! 
!     write(arquivo,'(3a)')'SCALARS ', trim(rotulo), ' float '
!     write(arquivo,'(a)')'LOOKUP_TABLE default'
! 
! !     call escreverEscalaresNodais(arquivo,campo, dim1, dim2,rotulo,tamRot)
!     call escreverEscalaresPorElemento(arquivo,campo, dim2,tamRot)
! 
! 
! 
!     end subroutine escreverArqParaview

!**** new **********************************************************************
      subroutine escreverPontosNodais  (ulEscrita,coords, numnp, nsd)
      implicit none
      integer, intent(in) :: ulEscrita,numnp, nsd
      real*8,  intent(in) :: coords(nsd,numnp)
!
      real*8  :: coordZ = 0.0 
      integer :: d, i
!
        do i=1,numnp
            write(ulEscrita,'(3(1x, 1pe15.8))') coords(1:nsd,i)  !BD
        end do

      end subroutine escreverPontosNodais


!**** new **********************************************************************
      subroutine escreverConectividades(ulEscrita,conectElem, numel, nen, nsd)
      implicit none
      integer, intent(in)  :: ulEscrita,numel, nen, nsd
      integer, intent(in)  :: conectElem(nen,numel)
!
      integer n, i
!
      if(nsd==2) then
      do  n=1,numel
        write(ulEscrita,'(i10,9(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
      end do
      end if

      if(nsd==3) then
      do  n=1,numel
        write(ulEscrita,'(i10,18(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
      end do
      end if

 end subroutine escreverConectividades

!**** new **********************************************************************
      subroutine escreverTiposElementos(ulEscrita,numel, nsd)
      implicit none
      integer, intent(in)   :: ulEscrita, numel, nsd
!
      integer :: i
!
      if(nsd==2) then
      do  i =1,numel
        write(ulEscrita,'(a)') '9'!trim(adjustl(tipo))
      end do
      end if 
!
      if(nsd==3) then
      do  i =1,numel
        write(ulEscrita,'(a)') '12'!trim(adjustl(tipo))
      end do
      end if 

      end subroutine escreverTiposElementos

!**** new **********************************************************************
      subroutine escreverEscalaresNodais(ulEscrita,v, tam1, tam2, rotulo, tamRot)
      implicit none
      integer, intent(in)  :: ulEscrita,tam1,tam2
      real*8, intent(in)   :: v(tam1,tam2)
      integer :: tamRot
      character(len=tamRot) :: rotulo
!
      character(len=tamRot+5) ::  rotuloN
      integer :: i,j
      character(len=5):: eixo
      real*8 :: limite,zero

      limite=1.e-15
      zero=0.0d0
      do i=1,tam1

        if(i>1) then
           write(eixo,'(i0)') i
           rotuloN=trim(rotulo)//trim(eixo)
           write(ulEscrita,'(3a)')'SCALARS ', trim(rotuloN), ' float '
           write(ulEscrita,'(a)')'LOOKUP_TABLE default'
        endif

        do j=1, tam2
             if(v(i,j).lt.limite) then
                write(ulEscrita,*) zero
             else 
                write(ulEscrita,*) v(i,j)
             end if
        end do
       
      end do

      end subroutine escreverEscalaresNodais

!**** new **********************************************************************
      subroutine escreverEscalaresPorElemento(ulEscrita,v, tam, tamRot)
      implicit none
      integer, intent(in)  :: ulEscrita,tam
      real*8, intent(in)   :: v(tam)
      integer :: tamRot
!
      integer :: j
      real*8 :: limite,zero

      limite=1.e-15
      zero=0.0d0

      do j=1, tam
           if(v(j).lt.limite) then
                write(ulEscrita,*) zero
             else
                write(ulEscrita,*) v(j)
             end if
      end do

      end subroutine 

!**** new **********************************************************************

!     subroutine escreverArqParaviewIntermed(arquivo, campo, dim1, dim2, rotulo, tamRot)
!     use mMalha, only: x, nsd, numel, numnp
! 
!     implicit none
!     integer, intent(in) :: arquivo,dim1, dim2
!     double precision, intent(in) :: campo(dim1, dim2)
! 
!     integer :: tamRot
!     character(len=tamRot) :: rotulo
! 
!     write(arquivo,'(3a)')'SCALARS ', trim(rotulo), ' float '
!     write(arquivo,'(a)')'LOOKUP_TABLE default'
! 
!      call escreverEscalaresNodais(arquivo,campo, dim1, dim2, rotulo, tamRot)
! 
!      end subroutine escreverArqParaviewIntermed
!
!----------------------------------------------------------------------
!
      subroutine paraview_geraCase(steps)

      implicit none

      integer :: steps
!
      integer :: numInicial, incremento
      real*8  :: incTempo
      integer :: i

      numInicial=0
      incremento=1
      incTempo =0.0

      open(unit=124,file="./out/transiente.case",status="unknown")

      write(124, "('FORMAT',/,'type:',2x,'ensight')")
      write(124, *)
      write(124, "('GEOMETRY',/,'model:',2x,'solucao.geo')")
      write(124, *)
      write(124, "('VARIABLE',/,'scalar per element:', 2x, 'Saturacao', 2x, 'solucao.***' )")
      write(124, *)
      write(124, "('TIME',/,'time set: 1')")
      write(124, "('number of steps:', i10)"), steps
      write(124, "('filename start number:', i10)"), numInicial
      write(124, "('filename increment:', i10)"), incremento
      write(124, "('time values:')")

      do i=1, steps+1
           write(124, *), incTempo
           incTempo=incTempo+1.0
      end do

      end subroutine
!
!----------------------------------------------------------------------
!
    subroutine paraview_geometria(numel,numnp,nsd,x,conecNodaisElem)

    implicit none

      integer :: numel,numnp,nsd
      real(8), dimension(nsd,*) :: x   
      integer :: conecNodaisElem(8,numel)
!
      integer :: i
      open(unit=125,file="./out/solucao.geo",status="unknown")

      write(125,'(a)')'Title1'
      write(125,'(a)')'Title2'
      write(125,'(a)')'node id given'
      write(125,'(a)')'element id given'
      write(125,'(a)')'coordinates'
      write(125,'(i8)')  numnp

      do i = 1, numnp
      WRITE (125,'(I8,3E12.5)') I,x(1,i),x(2,i),x(3,i)
      enddo

      WRITE (125,'(A,/,A,/,A,/,I8)')                     &
                                'part 1'           ,    &
                                'malha'            ,    &
                                'hexa8'            ,    &
                                 numel

      WRITE (125,'(9I8)')  (I,conecNodaisElem(1,i),conecNodaisElem(2,i),conecNodaisElem(3,i), &
                                        conecNodaisElem(4,i),conecNodaisElem(5,i),conecNodaisElem(6,i), & 
                                        conecNodaisElem(7,i),conecNodaisElem(8,i),i=1, numel ) 

     end subroutine paraview_geometria

!
!*****************************************************
!
!       subroutine imprimirCaseParaview(x, conecNodaisElem, pressaoElem, satElem, phi, perm)
! 
!       use mMalha,               only: numnp,nsd,numel,nen
! 
!       implicit none
! 
!       real*8 :: x(nsd, numnp)
!       integer :: conecNodaisElem(nen,*)
!       real*8 ::  pressaoElem(*), satElem(*), phi(*), perm(*)
! !
!       if(iflag_sat==1) then    
!          if(iflag_tipoPrint==1) then 
!             open(unit=ipress    , file= './out/solucao.P0001')
!             open(unit=iperm     , file= './out/solucao.K0001')
!             open(unit=iporo     , file= './out/solucao.PHI0001')
! !             open(unit=iveloc      , file= 'solucao.V0001')
!             open(unit=isaturacao, file= './out/solucao.S0001')
! 
!             call paraview_geraCase(qtdImpSat)
!             call paraview_geometria(numel,numnp,nsd,x, conecNodaisElem)
!             call paraview_escalarPorElemento(numel, pressaoElem,ipress)
!             call paraview_escalarPorElemento(numel, perm,     iperm)
!             call paraview_escalarPorElemento(numel, phi,        iporo)
!             call paraview_escalarPorElemento(numel, satElem,    isaturacao)
!          endif
!       endif
! 
!       end subroutine
! 
!
!----------------------------------------------------------------------
!
      subroutine paraview_vetorPorElemento(numel,campo,iarq)

! ainda precisa implementar
      implicit none
!
      integer :: numel
      real(8), dimension(*) :: campo
      integer :: iarq

      integer :: i

!
      write(iarq,"('Ensight Scalar passo     1')")
      write(iarq,"('part 1')")
      write(iarq,"('hexa8')")
!
      write(iarq,"(6e12.5)") (campo(i),i=1,numel)
!      
      close(iarq)
!
      end subroutine    
!
!----------------------------------------------------------------------
!
      subroutine paraview_escalarPorElemento(numel,campo,iarq)
      implicit none
!
      integer :: numel
      real(8), dimension(*) :: campo
      integer :: iarq

      integer :: i

      print*, "gerando", iarq
!
      write(iarq,"('Ensight Scalar passo     1')")
      write(iarq,"('part 1')")
      write(iarq,"('hexa8')")
!
      write(iarq,"(6e12.5)") (campo(i),i=1,numel)
!      
      close(iarq)
!
      end subroutine    
!
!
!----------------------------------------------------------------------
!
      subroutine paraview_escalarPorElementoTransiente(numel,campo,passo,iarq)
      implicit none
!
      integer :: numel,passo,iarq
      real(8), dimension(*) :: campo   
      character(len=128) :: name,sol
      character(len=8)   :: c
      integer :: i
      real(4) :: x     
!      
      x=0.001

      sol="solucao"

      write(c,"(f7.3)") x*passo
      c=adjustl(c)
      name='./out/'//trim(sol)//c 
!      
      open(unit=iarq,file=name,status="unknown")
!
      write(iarq,"('Ensight Scalar passo ',i5)") passo
      write(iarq,"('part 1')")
      write(iarq,"('hexa8')")
!      
!     imprime as coordenadas
!
       write(iarq,"(6(e12.5))") (real(campo(i)),i=1,numel)
!      
      close(iarq)
!
      passo=passo+1
!
      end subroutine          
!
!=======================================================================
!     
!       subroutine prt(nsd,numel,t0,u,iunit)
! !      
!       use mMalha, only: xc
! !      
!       implicit none
! !     
! !     imprime campos escalares para o gnuplot ou para o matlab
! !
!       integer                   :: nsd,numel
!       real(8), dimension(*)     :: u
!       real(8)                   :: t0
!       integer :: iunit
! !     
!       integer :: nel
! !     
!       write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!       write(iunit,*)
! !     
!       do nel=1,numel
!          write(iunit,"(5(f25.15,2x))") xc(1:nsd,nel),u(nel)
!       end do
! !     
!       write(iunit,*)
! !     
!       end subroutine
!     
!=======================================================================
!
!       subroutine prtvB(nsd,numel,t0,velocLadal,ndofV,  conecLadaisElem, numLadosElem, iunit)
! !
!       use mMalha, only: xc
! !
!       implicit none
! !
! !     imprime campos vetoriais para o gnuplot ou para o matlab
! !
!       integer                   :: numel,nsd, ndofV ,numLadosElem
!       real(8), dimension(ndofV,*) :: velocLadal
!       real(8)                   :: t0
!       integer                   :: conecLadaisElem(numLadosElem,numel)
! !
!       integer :: nel
!       real(8) :: vc(nsd)
!       real*8 :: mediaCentro
! !
!       integer :: iunit
! 
!       vc=0.0
! !
!       write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
! !
!       write(iunit,*)
! !
!       do nel=1,numel
! !
!         vc(1) = (velocLadal(1,conecLadaisElem(2,nel))+velocLadal(1,conecLadaisElem(4,nel)))/2.0
!         vc(2) = (velocLadal(1,conecLadaisElem(1,nel))+velocLadal(1,conecLadaisElem(3,nel)))/2.0
!         if(nsd==3) vc(3) = (velocLadal(1,conecLadaisElem(5,nel))+velocLadal(1,conecLadaisElem(6,nel)))/2.0
!         mediaCentro=sum(vc)/nsd
!         write(iunit,"(6(f25.15,2x))")xc(1:nsd,nel),mediaCentro
! 
!       end do
! !
!       write(iunit,*)
! !
!       end subroutine

!     
!=======================================================================
!     
!       subroutine prtv(nsd,numel,ndofV,numLados,t0,u,iunit)
! !      
!       use mMalha, only: xc
! !
!       implicit none
! !     
! !     imprime campos vetoriais para o gnuplot ou para o matlab
! !
!       integer                   :: numel,nsd, ndofV, numLados
!       real(8), dimension(ndofV,numLados) :: u
!       real(8)                   :: t0
! !
!       integer :: nel
! !     
!       integer :: iunit
! !     
!       write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
! !
!       write(iunit,*)
! !     
!       do nel=1,numel
! 
!          write(iunit,"(5(f25.15,2x))") xc(1:nsd,nel),u(1,nel),u(2,nel)
!       end do
! !     
!       write(iunit,*)
! !     
!       end subroutine

!
!**** NEW **** FOR DATA EXPLORER OUT PUT ************************************* 
!
!       SUBROUTINE PRINT_DXINFO(TASK,NUMNP,NUMEL)
! 
!       use mGlobaisEscalares, only: NUMDX, NNP, NVEL
! !
! !..... PROGRAM TO SET-UP AND WRITE DATA ON OPEN-DX FORMAT
! !
! !
!       IMPLICIT REAL*8 (A-H,O-Z)
! !
!       CHARACTER*15 TASK
! !
!       INTEGER JJ, NUMNP, NUMEL, NINDX, NNSTEP, NTINDX
! !
!       IF (NUMDX.EQ.0) RETURN
! 
!         IF(TASK=='OPEN__FEDX_FILE') THEN
! ! 
! !..... PRINT NODAL AND MESH INFORMATION FOR DATA EXPLORER FILE
! ! 
! 	WRITE(IFEDX,1000) '## OpenDX format File' 
! 	WRITE(IFEDX,1000) '## OutPut Data  at Nodal Points in the' 
! 	WRITE(IFEDX,1000) '## sense of Finite Element Method'
!         WRITE(IFEDX,1000) '##=========================================='
!         WRITE(IFEDX,1000) '##=========================================='
! ! 
! !..... Nodes of finite element mesh 
! !
! 	WRITE(IFEDX,1000) '# ' 
! 	WRITE(IFEDX,1000) '## Nodes locations'
! 	WRITE(IFEDX,1000) '# '  
! 	WRITE(IFEDX,1500)  &
!      & 'object 1 class array type float rank 1 shape 2 items ', &
!      & NUMNP,' data file "fnodes.stoc"' 
! !
! !..... Conectivity of finite element mesh
! !
! 	WRITE(IFEDX,1000) '# ' 
! 	WRITE(IFEDX,1000) '## Connectivity' 
! 	WRITE(IFEDX,1000) '# ' 
! 	WRITE(IFEDX,1500) &
!      & 'object 2 class array type int rank 1 shape 4 items ', &
!      & NUMEL,' data file "femesh.stoc"' 
! 	WRITE(IFEDX,1000) 'attribute "element type" string "quads"' 
! 	WRITE(IFEDX,1000) 'attribute "ref" string "positions"' 
! 	WRITE(IFEDX,1000) '#  '
! ! 
!        ENDIF
! !
!        IF(TASK=='WRITE_FEDX_DATA') THEN
! !
! !..... PRINT NODAL DATA FOR OPEN-DX FILE
! ! 
! !        NINDX = NUSTEP/10
! !
!         NINDX = NNP/NUMDX
! !
! !        NINDX = IDINT(TPRT_PHI/DTPRT_PHI)
! !
!         NNSTEP=24*NINDX+3
! ! 
! !...... DISPLACEMENTS: VECTOR FIELD
! ! 
! 	WRITE(IFEDX,1000)'# Vector field : DISPLACEMENTS' 
! 	WRITE(IFEDX,2000) NNSTEP, NUMNP, NINDX 
! ! 
! !...... DISPLACEMENT SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   DISPLACEMENT series'
! 	WRITE(IFEDX,3000) NNSTEP+1, NNSTEP
! ! 
! !.....  PRESSURE: SCALAR FIELD
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : PRESSURE' 
! 	WRITE(IFEDX,4000) NNSTEP+2, NUMEL, NINDX 
! ! 
! !...... PRESSURE SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar PRESSURE series'
! 	WRITE(IFEDX,3000) NNSTEP+3, NNSTEP+2
! ! 
! !.....  GEOMECHANIC POROSITY: PORE SCALAR FIELD
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : PORE' 
! 	WRITE(IFEDX,4010) NNSTEP+4, NUMEL, NINDX 
! ! 
! !...... GEOMECHANIC POROSITY: PORE SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar PORE series'
! 	WRITE(IFEDX,3000) NNSTEP+5, NNSTEP+4
! ! 
! !.....  CREEP SCALAR FIELD
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : ' 
! 	WRITE(IFEDX,4020) NNSTEP+6, NUMEL, NINDX 
! ! 
! !...... CREEP SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar CREEP series'
! 	WRITE(IFEDX,3000) NNSTEP+7, NNSTEP+6
! ! 
! !.....  SATURATION: SCALAR FIELD
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : SATURATION' 
! 	WRITE(IFEDX,4030) NNSTEP+8, NUMEL, NINDX 
! ! 
! !...... SATURATION SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar SATURATION series'
! 	WRITE(IFEDX,3000) NNSTEP+9, NNSTEP+8
! ! 
! !.....  VELOCITY VECTOR FIELD
! ! 
! 	WRITE(IFEDX,1000)'# Vector field: VELOCITY' 
! 	WRITE(IFEDX,4040) NNSTEP+10, NUMEL, NINDX 
! ! 
! !...... VELOCITY SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Vector Field VELOCITY series'
! 	WRITE(IFEDX,3000) NNSTEP+11, NNSTEP+10
! ! 
! !.....  PERMEABILITY: SCALAR FIELD
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : PERMEABILITY' 
! 	WRITE(IFEDX,4050) NNSTEP+12, NUMEL, NINDX 
! ! 
! !...... PERMEABILITY SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar PERMEABILITY series'
! 	WRITE(IFEDX,3000) NNSTEP+13, NNSTEP+12
! ! 
! !.....  STRESS ON X DIRECTION
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : Stress_X' 
! 	WRITE(IFEDX,4060) NNSTEP+14, NUMEL, NINDX 
! ! 
! !...... STRESS_X SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar Stress_X series'
! 	WRITE(IFEDX,3000) NNSTEP+15, NNSTEP+14
! ! 
! !.....  STRESS ON Y DIRECTION
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : Stress_Y' 
! 	WRITE(IFEDX,4070) NNSTEP+16, NUMEL, NINDX 
! ! 
! !...... STRESS_Y SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar Stress_Y series'
! 	WRITE(IFEDX,3000) NNSTEP+17, NNSTEP+16
! ! 
! !.....  SHEAR STRESS XY
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : Stress_XY'
! 	WRITE(IFEDX,4080) NNSTEP+18, NUMEL, NINDX 
! ! 
! !...... SHEAR STRESS_XY SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar Stress_XY series'
! 	WRITE(IFEDX,3000) NNSTEP+19, NNSTEP+18
! ! 
! !.....  STRESS Z
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : Stress_Z'
! 	WRITE(IFEDX,4090) NNSTEP+20, NUMEL, NINDX 
! ! 
! !...... SHEAR STRESS_XY SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar Stress_Z series'
! 	WRITE(IFEDX,3000) NNSTEP+21, NNSTEP+20
! ! 
! !.....  YOUNG MODULUS 
! ! 
! 	WRITE(IFEDX,1000)'# Scalar field : Young Modulus'
! 	WRITE(IFEDX,4100) NNSTEP+22, NUMEL, NINDX 
! ! 
! !...... YOUNG MODULUS SERIES INFORMATION 
! ! 
! 	WRITE(IFEDX,1000)'# Next object is a member of the: '
! 	WRITE(IFEDX,1000)'#   Scalar Young modulus series'
! 	WRITE(IFEDX,3000) NNSTEP+23, NNSTEP+22
! !
!       ENDIF
! !
!       IF(TASK=='CLOSE_FEDX_FILE') THEN
! ! C
! ! C..... PRINT SERIES LINKS INFORMATION FOR DATA EXPLORER FILE
! ! C
! ! C
! ! c        NTINDX = NP
! ! C        NTINDX = NVEL
!         NTINDX=(NVEL/NUMDX) 
! 
! !        NTINDX = IDINT(TPRT_PHI/DTPRT_PHI)+1
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the DISPLACEMENT series object'
!       WRITE(IFEDX,1000)'object "displacement" class series'
!       DO 301 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-20,JJ-1
!  301  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the PRESSURE series object'
!       WRITE(IFEDX,1000)'object "pressure" class series'
!       DO 302 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-18,JJ-1
!  302  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the PORE series object'
!       WRITE(IFEDX,1000)'object "pore" class series'
! !
!       DO 303 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-16,JJ-1
!  303  CONTINUE
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the CREEP series object'
!       WRITE(IFEDX,1000)'object "creep" class series'
!       DO 304 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-14,JJ-1
!  304  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the SATURATION series object'
!       WRITE(IFEDX,1000)'object "saturation" class series'
!       DO 305 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-12,JJ-1
!  305  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the VELOCITY series object'
!       WRITE(IFEDX,1000)'object "velocity" class series'
!       DO 306 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-10,JJ-1
!  306  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the PERMEABILITY series object'
!       WRITE(IFEDX,1000)'object "permeability" class series'
!       DO 307 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-8,JJ-1
!  307  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the STRESS X series object'
!       WRITE(IFEDX,1000)'object "stress_x" class series'
!       DO 308 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-6,JJ-1
!  308  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the STRESS Y series object'
!       WRITE(IFEDX,1000)'object "stress_y" class series'
!       DO 309 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-4,JJ-1
!  309  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the SHEAR STRESS series object'
!       WRITE(IFEDX,1000)'object "stress_xy" class series'
!       DO 310 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ-2,JJ-1
!  310  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the STRESS Z series object'
!       WRITE(IFEDX,1000)'object "stress_z" class series'
!       DO 311 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ,JJ-1
!  311  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Here we create the YOUNG MODULUS serie object'
!       WRITE(IFEDX,1000)'object "young" class series'
!       DO 312 JJ=1,NTINDX
!         WRITE(IFEDX,7000) JJ-1,24*JJ+2,JJ-1
!  312  CONTINUE
! !
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'# Structure of VARIAVEL OF DATA FILE'
!       WRITE(IFEDX,1000)'object "campos" class group'
!       WRITE(IFEDX,1000)'member "displacement" value "displacement"'
!       WRITE(IFEDX,1000)'member "pressure" value "pressure"'
!       WRITE(IFEDX,1000)'member "pore" value "pore"'
!       WRITE(IFEDX,1000)'member "creep" value "creep"'
!       WRITE(IFEDX,1000)'member "saturation" value "saturation"'
!       WRITE(IFEDX,1000)'member "velocity" value "velocity"' 
!       WRITE(IFEDX,1000)'member "permeability" value "permeability"' 
!       WRITE(IFEDX,1000)'member "stress_x" value "stress_x"' 
!       WRITE(IFEDX,1000)'member "stress_y" value "stress_y"' 
!       WRITE(IFEDX,1000)'member "stress_xy" value "stress_xy"' 
!       WRITE(IFEDX,1000)'member "stress_z" value "stress_z"' 
!       WRITE(IFEDX,1000)'member "young" value "young"' 
!       WRITE(IFEDX,1000)'#  ' 
!       WRITE(IFEDX,1000)'end'
!       WRITE(IFEDX,1000)'#  ' 
! !
!       ENDIF
! !
!       IF(TASK=='OTHERS__ANOTHER') THEN
! !
!       ENDIF
! !
! !.... FORMATOS DE SAIDA  OPEN-DX
! !
!  1000 FORMAT(A) 
!  1500 FORMAT(A,I7,A)
!  1800 FORMAT(27I6)
!  2000 FORMAT('object ',I5,                                 &
!      & ' class array type float rank 1 shape 2 items', I8, &
!      &' data file "disp',I3.3,'.stoc"'/                    &
!      &'attribute "dep" string "positions"'/'#  ') 
! !
!  3000 FORMAT('object ',I5,' class field'/    &
!      &'component "positions" value 1'/       &
!      &'component "connections" value 2'/     &
!      &'component "data" value ',I5/'#  ')
! !
!  4000 FORMAT('object ',I5,                     &
!      &' class array type float rank 0 items ', &
!      & I8,' data file "prsr',I3.3,'.stoc"'/    &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
!  4010 FORMAT('object ',I5,                       &
!      &' class array type float rank 0 items ',   &
!      & I8,' data file "pore',I3.3,'.stoc"'/      &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
!  4020 FORMAT('object ',I5,                       &
!      &' class array type float rank 0 items ',   &
!      & I8,' data file "crep',I3.3,'.stoc"'/      &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
!  4030 FORMAT('object ',I5,                       &
!      &' class array type float rank 0 items ',   &
!      & I8,' data file "satr',I3.3,'.stoc"'/      &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
!  4040 FORMAT('object ',I5,                                 &
!      & ' class array type float rank 1 shape 2 items', I8, &
!      &' data file "velt',I3.3,'.stoc"'/                    &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
!  4050 FORMAT('object ',I5,                       &
!      &' class array type float rank 0 items ',   &
!      & I8,' data file "perm',I3.3,'.stoc"'/      &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
!  4060 FORMAT('object ',I5,                       &
!      &' class array type float rank 0 items ',   &
!      & I8,' data file "sigx',I3.3,'.stoc"'/      &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
!  4070 FORMAT('object ',I5,                       &
!      &' class array type float rank 0 items ',   &
!      & I8,' data file "sigy',I3.3,'.stoc"'/      &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
! 4080  FORMAT('object ',I5,                         &
!      & ' class array type float rank 0 items', I8, &
!      &' data file "sigt',I3.3,'.stoc"'/            &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
! 4090  FORMAT('object ',I5,                         &
!      & ' class array type float rank 0 items', I8, &
!      &' data file "sigz',I3.3,'.stoc"'/            &
!      &'attribute "dep" string "connections"'/'#  ') 
! !
! 4100  FORMAT('object ',I5,                          &
!      & ' class array type float rank 0 items', I8,  &
!      &' data file "yung',I3.3,'.stoc"'/             &
!      &'attribute "dep" string "connections"'/'#  ') 
! 7000  FORMAT('member ',I5,' value ',I5,' position ',I5)
! !
!       END SUBROUTINE
!
!*** NEW ***************************************************************
!

!*** Diego (Jul/2015) ********************************************************************** 
!    Esta rotina tem como finalidade a escrita dos valores dos deslocamentos calculados
!    para o caso da elasticidade linear.

      SUBROUTINE PRINTDISP(U,X,NUMNP,idx)
      
      IMPLICIT NONE
      
      REAL*8   u(2,*),X(2,*)
      INTEGER  NUMNP, N, idx, idrx
      CHARACTER*30  idxStr, disp
            
!       Abrindo os arquivos para saída

      write(idxStr,'(i1)') idx
      
      idrx = 11*idx
      disp = 'disp.'//idxStr
      OPEN(UNIT=idrx, FILE= disp)
      
      do N=1,NUMNP
	WRITE(idrx,210) N, X(1,N), X(2,N), u(1,N), u(2,N)
      enddo
      
 ! 4 espaços, inteiro max 5 posicoes, 10 espacos, 4 floats 8.2 com espaco de 2 entre eles
 210  FORMAT(4X,I5,10x,5(1PE15.8,2X))

      close(idrx)

      END subroutine    
!
!=================================================================================
!
    subroutine abrirArquivosDS()

!        iin    = input unit number
!        iecho  = output unit of input data
!        iouter  = output unit of error norms
!
    character(len=20) :: nomeIn, nomeEcho
!
      iin        = 15
      iecho      = 16 
      icoords    = 18
      iconects   = 19
!
      ignuplot = 30
!       ignuplotFluxo   = 31
!       iparaview       = 32
!
      nomeIn='input.dat'
      nomeEcho='echo.dat'
!
      open(unit=iin,    file=nomeIn, status='old', err=100) 
      open(unit=iecho , file=nomeEcho)
!
      open(unit=icoords   ,file= 'coordenadas.dat')
      open(unit=iconects  ,file= 'conectsNodais.dat')
!
      open(unit=ignuplot ,file= 'resultadoP.dat')
!       open(unit=ignuplotFluxo     ,file= 'resultadoF.dat')
!       open(unit=iparaview ,file= 'resultado.vtk')

      return 

      100 continue
      write(*,*) ' arquivo: ', nomeIn, ' NAO encontrado'
      stop 1
       
   end subroutine abrirArquivosDS
!
!=================================================================================
!
    subroutine fecharArquivosDS()

      close(iin      )
      close(iecho    )
      close(icoords  )
      close(iconects )
      close(ignuplot)
!       close(ignuplotFluxo   )
!       close(iparaview)
    end subroutine fecharArquivosDS
! 
! ===================================================================================
!
      SUBROUTINE SETUPDX()

      use mGlobaisEscalares, only: NUMDX, NITGEO
!..
!...  PROGRAM TO SETUP MANAGER FILES FOR GRAPHICAL INTERFACE OPEN-DX
!
      CHARACTER*30 NIFEDX,NIFNOD,NIFMSH,NIRBAL
!
      CHARACTER*2 ASTEP
!
      IF (NUMDX.EQ.0) RETURN
!
      WRITE(ASTEP,'(I2.2)') NITGEO
!
      NIFEDX = 'dxcreep'//ASTEP//'/nodestoc.dx'
      NIFNOD = 'dxcreep'//ASTEP//'/fnodes.stoc'
      NIFMSH = 'dxcreep'//ASTEP//'/femesh.stoc'
      NIRBAL = 'dxcreep'//ASTEP//'/global.mass'
!
      OPEN(UNIT=IFEDX, FILE= NIFEDX)
      OPEN(UNIT=IFNOD, FILE= NIFNOD)
      OPEN(UNIT=IFMSH, FILE= NIFMSH)
      OPEN(UNIT=IRBAL, FILE= NIRBAL)
!
      END SUBROUTINE
! 
! !
! !**** new **********************************************************************
! !
!       subroutine prntel(mat,conectElem,nen,numel,tipo)
!       implicit none
! !
! !.... program to print data for element with "nen" nodes
! !
! !        note: presently the label formats are limited to
! !              elements with one to nine nodes
! !
!       integer*4:: nen, numel
!       integer*4:: mat(*),conectElem(nen,*)
!       integer*4:: tipo
! !
!       integer*4 n, i
! !
!       if(tipo==1) then
! 	write(iconects,*) "# Conectividades nodais"
! 	do n=1,numel
! 	  write(iconects,2000) n,mat(n),(conectElem(i,n),i=1,nen)
! 	end do
!       end if
! 
!       if(tipo==2) then
! 	write(iconectsL,*) "# Conectividades ladais"
! 	do  n=1,numel
! 	  write(iconectsL,3000) n,mat(n),(conectElem(i,n),i=1,nen)
! 	end do
!       end if
! !
!       return
! !
!  1000 format(///10x,&
!      ' e l e m e n t   d a t a',//1x,&
!      ' element   material ',9('node ',i1,1x),/1x,&
!      '  number    number',//)
!  2000 format(1x,i10,9(2x,i10))
!  3000 format(1x,i10,7(2x,i10))
!       end subroutine
! 

      END MODULE
