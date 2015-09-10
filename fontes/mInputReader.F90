!!>
!!         programa de elementos finitos em fortran 90
!!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!!
!!         Pilato Jr , Vinicius Pessoa  {pilato@deepsoft.com.br, pessoa@deepsoft.com.br}
!!
!!         Desenvolvido por DeepSoft para LNCC/MCT
!!         Rio de Janeiro, 11.2013-04.2014
!!
!!=================================================================================

!> Modulo responsavel por reunir subrotinas para leitura do arquivo de entrada.
module mInputReader

    !> Armazena as linhas do arquivo de input.
    character(len=200), allocatable :: file_lines(:)
    !> Armazena o numero de linhas no arquivo.
    integer*4 number_of_lines

    contains

    !> Leitura e geracao de valores de condicoes de contorno.
    !!
    !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
    !! @param f             f.
    !! @param ndof          ndof
    !! @param numnp         numnp
    !! @param j             j
    !! @param nlvect        nlvect
    !! @param iprtin        iprtin
    subroutine leituraValoresCondContornoDS(keyword_name, f,ndof,numnp,j, nlvect,iprtin)
!     Alterado por Diego para entrada de valor de cond de contorno adimensionalizada. (Ago/2015)
        use mleituraEscrita, only: iecho, printf, printd
        use mParametros,     only: condContAdim
!        use mInputReader,    only: findKeyword, genflDS

        implicit none

        integer*4 :: ndof, numnp, j, nlvect, iprtin, keyword_line
        character(len=50) :: keyword_name
        real*8 f(ndof,numnp,nlvect)

        logical lzero
        integer*4 nlv
        character(len=35) :: rotulo

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.-1) return

        !
        !call clear(f,nlvect*numnp*ndof)
         f=0.0
! 	Para a CC adim
!        f(1:nlvect,1:numnp,1:ndof)=condContAdim

        do 100 nlv=1,nlvect
            call genflDS(f(1,1,nlv),ndof,keyword_line)
            call ztest(f(1,1,nlv),ndof*numnp,lzero)
            if (iprtin.eq.0) then
                 if (lzero) then
                    if (j.eq.0) write(iecho,1000) nlv
                    if (j.eq.1) write(iecho,2000)
                 else
                    if (j.eq.0) call printf(f,ndof,numnp,nlv)
                    if (j.eq.1) then
                        rotulo=" n o d a l  b o d y  f o r c e s"
                        call printd (rotulo, f,ndof,numnp,iecho)
                    end if
                endif
            endif
        100 continue
!         write(*,*) f(1,1:60,1), " fim do vetor forca "
!         stop
        return
 1000 format('1'//,' there are no nonzero prescribed forces and ',&
         'kinematic boundary conditions for load vector number ',i10)
 2000 format('1'//,' there are no nonzero nodal body forces')
    end subroutine leituraValoresCondContornoDS !**********************************************************************


    !> Leitura e geracao de codigos de condicoes de contorno.
    !!
    !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
    !! @param id            id.
    !! @param ndof          numnp
    !! @param numnp         numnp
    !! @param neq           neq
    !! @param iecho         iecho
    !! @param iprtin        iprtin
    subroutine leituraCodigosCondContornoDS(keyword_name, id, ndof, numnp, neq, iecho, iprtin)
!        use mInputReader,     only: findKeyword, igenDS


        integer*4, intent(in) :: ndof, numnp, iecho, iprtin
        integer*4 ::  neq, keyword_line
        integer*4, intent(inout) :: id(ndof,numnp)
        character(len=50) keyword_name

        integer*4:: nn, n, i
        logical pflag

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.-1) return

        id = 0
        call igenDS(id,ndof, keyword_line)
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
    end subroutine leituraCodigosCondContornoDS !********************************************************************



    !> Leitura e geracao de coordenadas
    !!
    !! @param x          Matriz onde serao armazeadas as coordenadas.
    !! @param nsd        Corresponde ao nao de linhas da matriz x.
    !! @param numnp      Numero de xxx
    !! @param icoords    Handle para o arquivo de coordenadas.
    !! @param iprtin     iprtin
    subroutine leituraGeracaoCoordenadasDS(x, nsd, numnp, icoords, iprtin)

      implicit none

      integer*4, intent(in)   :: nsd, numnp, icoords, iprtin
      real*8, intent(inout) ::  x(nsd,*)
      integer*4:: i, n, keyword_line
      character(len=50) keyword_name

      keyword_name = 'coordenadas_nodais_bm'
      keyword_line = findKeyword(keyword_name)

      call genflDS(x,nsd,keyword_line)

      if (iprtin.eq.1) return
      write(icoords,*) "# Coordenadas ", nsd
      do n=1,numnp
        write(icoords,2000) n,(x(i,n),i=1,nsd)
      end do

      return
      2000 format(6x,i12,10x,3(1pe15.8,2x))
    end subroutine leituraGeracaoCoordenadasDS!***************************************************************************************************


    !---------------------------------------------------------------------------------
    !> Le arquivo de input e armazena seu conteudo em um array.
    !! @param file_name Nome do arquivo a ser lido.
    subroutine readInputFileDS()
        use mLeituraEscrita,   only: iin

        implicit none

        integer*4 success, lines_count
        character(len=200) file_line

        integer*4 :: main_number_of_lines, main_number_of_includes, i, include_index, inc_nlines, inc_inc, merge_lines
        character(len=200), allocatable :: main_file_lines(:)
        character(len=200), allocatable :: include_files(:)
        character(len=200) include_file
        integer*4, allocatable :: include_indexes(:), include_number_of_lines(:)

        main_number_of_lines = 0
        main_number_of_includes = 0

        call analyzeFileInput(main_number_of_lines, main_number_of_includes)

        if (main_number_of_includes.eq.0) then
            call createSimpleInputFile()
            return
        end if

        allocate(main_file_lines(main_number_of_lines))
        allocate(include_indexes(main_number_of_includes))
        allocate(include_number_of_lines(main_number_of_includes))
        allocate(include_files(main_number_of_includes))

        lines_count = 1
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            main_file_lines(lines_count) = file_line
            lines_count = lines_count + 1
        end do
        rewind(iin)

        number_of_lines = main_number_of_lines

        !Number of lines
        do i=1, main_number_of_includes
            include_index = findInclude(i, main_file_lines, main_number_of_lines)
            read(main_file_lines(include_index), '(A)') include_file

            include_file = adJustl(include_file)
            call analyzeFile(include_file, inc_nlines, inc_inc)
            number_of_lines = number_of_lines + inc_nlines
            include_indexes(i) = include_index
            include_number_of_lines(i) = inc_nlines
            include_files(i) = include_file
        end do

        !Prepare final struct.
        call prepareFileLines(include_indexes, include_number_of_lines, main_number_of_includes, main_file_lines)

        !Merge contensts.
        merge_lines = 0
        do i=1, main_number_of_includes
            call mergeIncludeContents(include_files(i), include_indexes(i) + merge_lines)
            merge_lines = merge_lines + include_number_of_lines(i)
        end do

        deallocate(main_file_lines)
        deallocate(include_indexes)
        deallocate(include_number_of_lines)
        return
    end subroutine readInputFileDS !***********************************************************************************


    !> Cria a estretura de input usando um arquivo de entrada sem includes
    !! @param file_name Nome do arquivo a ser lido.
    subroutine createSimpleInputFile()
        use mLeituraEscrita,   only: iin

        implicit none

        integer*4 success, lines_count
        character(len=200) file_line

        number_of_lines = 0

        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            number_of_lines = number_of_lines + 1
        end do
        rewind(iin)

        allocate(file_lines(number_of_lines))
        !TO-DO avoid two-times read
        lines_count = 1
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            file_lines(lines_count) = file_line
            lines_count = lines_count + 1
        end do
        rewind(iin)
    end subroutine createSimpleInputFile !*****************************************************************************

    !> Le o conteudo do arquivo de include e armazena no array principal.
    !! @param   include_index   O index do include.
    !! @param   include_files   Array com includes.
    !! @param   include_line    A linha do include.
    subroutine mergeIncludeContents(include_file, include_line)

        implicit none

        integer*4 include_line
        character(len=200) include_file

        character(len=200) file_line
        integer*4 file_channel, success, current_index

        file_channel = 1

        current_index = include_line

        open(unit=file_channel, file=include_file)
        do
            read(file_channel, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            file_lines(current_index) = file_line
            current_index = current_index + 1
        end do
        close(file_channel)
    end subroutine mergeIncludeContents !******************************************************************************

    !> Efetua a aloca��o da estrutura definitiva, preparando a linha dos arquivos originais para receber os includes
    !! @param   include_indexes             Array os indices de ocorrencias dos includes.
    !! @param   include_number_of_lines     Array com o numero de linhas de cada include
    !! @param   number_of_includes          Numero de includes.
    !! @param   original_file_lines         Linhas do arquivo de entrada original.
    subroutine prepareFileLines(include_indexes, include_number_of_lines, number_of_includes, original_file_lines)
        integer*4 number_of_includes, number_of_original_lines, line_index, shift_lines
        integer*4 include_indexes(:), include_number_of_lines(:)
        character(len=200) original_file_lines(:)

        integer*4 current_include_index, original_index

        allocate(file_lines(number_of_lines))

        current_include_index = 1
        original_index = 1
        line_index = 1
        shift_lines = 0
        do while ( line_index <= number_of_lines)
            if (original_index.eq.(include_indexes(current_include_index))) then
                line_index = line_index + include_number_of_lines(current_include_index)
                current_include_index = current_include_index + 1
            end if
            file_lines(line_index) = original_file_lines(original_index)
            line_index = line_index + 1
            original_index = original_index + 1
        end do


    end subroutine prepareFileLines !**********************************************************************************

    !> Efetua algumas an�lises no arquivo recebido.
    !! @param   number_of_lines     N�mero de linhas.
    !! @param   number_of_include   N�mero de ocorr�ncias da palavra include.
    subroutine analyzeFileInput(number_of_lines, number_of_includes)
        use mLeituraEscrita,   only: iin

        character(len=200) file_line
        integer*4 number_of_lines, number_of_includes

        character(len=50) include_keyword, formated_keyword
        integer*4 keyword_len, success

        include_keyword = "include"
        keyword_len = len(trim(include_keyword)) + 2
        formated_keyword = trim('*' // trim(include_keyword) // '{')

        number_of_lines = 0
        number_of_includes = 0

!        lunitInicial = 15
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            number_of_lines = number_of_lines + 1
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                number_of_includes = number_of_includes + 1
            end if
        end do
        rewind(iin)

    end subroutine analyzeFileInput !***************************************************************************************

    !> Efetua algumas an�lises no arquivo recebido.
    !! @param   file_name           O nome do arquivo.
    !! @param   number_of_lines     N�mero de linhas.
    !! @param   number_of_include   N�mero de ocorr�ncias da palavra include.
    subroutine analyzeFile(file_name, number_of_lines, number_of_includes)
        character(len=200) file_name, file_line
        integer*4 number_of_lines, number_of_includes

        character(len=50) include_keyword, formated_keyword
        integer*4 keyword_len, file_channel, success

        include_keyword = "include"
        keyword_len = len(trim(include_keyword)) + 2
        formated_keyword = trim('*' // trim(include_keyword) // '{')

        number_of_lines = 0
        number_of_includes = 0

        file_channel = 2
        lunitInicial = 15
        file_channel = lunitInicial

        open(unit=file_channel, file=file_name)
        do
            read(file_channel, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            number_of_lines = number_of_lines + 1
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                number_of_includes = number_of_includes + 1
            end if
        end do
        close(file_channel)

    end subroutine analyzeFile !***************************************************************************************

    !> Procura a n-esima palavra-chave include.
    !! @param  position         Corresponde a posicao desejada.
    !! @param  file_lines       Linhas do arquivo.
    !! @param  number_of_lines  Numero de linhas atuais.
    !! @return O indice da palavra-chave no array que contem as linhas do arquivo de entrada.
    integer*4 function findInclude(position, file_lines, number_of_lines)
        implicit none
        integer*4 position, number_of_lines, current_position
        character(len=200) file_lines(:)
        character(50) keyword, formated_keyword
        character(len=120) file_line
        integer*4 i, keyword_len

        keyword = "include"
        keyword_len = len(trim(keyword)) + 2
        formated_keyword = trim('*' // trim(keyword) // '{')
        current_position = 0

        do i=1, number_of_lines
            file_line = file_lines(i)
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                current_position = current_position + 1
                if (current_position.eq.position) then
                    findInclude = i + 1
                    return
                end if
            end if
        end do
        findInclude = 0
        return
    end function findInclude !*****************************************************************************************

    !> Procura uma palavra-chave.
    !! @param  keyword A palavra-chave.
    !! @return O indice da palavra-chave no array que contem as linhas do arquivo de entrada.
    integer*4 function findKeyword(keyword)
        implicit none
        character(50) keyword, formated_keyword
        character(len=120) file_line
        integer*4 i, keyword_len
        do i=1, number_of_lines, 1
            file_line = file_lines(i)
            keyword_len = len(trim(keyword)) + 2
            formated_keyword = trim('*' // trim(keyword) // '{')
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                findKeyword = i + 1
                return
            end if
        end do
        findKeyword = 0
        return
    end function findKeyword !*****************************************************************************************

    !> Efetua a leitura de uma palavra-chave to tipo inteiro. Se nao encontrado, associa o valor defualt fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde o valor inteiro sera atribuido.
    !! @param default_value Valor default.
    subroutine readIntegerKeywordValue(keyword, target, default_value)
        implicit none
        character(50) keyword
        character(120) file_line
        integer*4 target, default_value, keyword_line
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            target = default_value
            return
        end if
        file_line = adjustL(trim(file_lines(keyword_line)))
        read(file_line, *) target
        return
    end subroutine readIntegerKeywordValue !***************************************************************************

    !> Efetua a leitura de uma palavra-chave to tipo string. Se nao encontrado, associa o valor defualt fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde a string sera atribuido.
    !! @param default_value Valor default.
    subroutine readStringKeywordValue(keyword, target, default_value)
        implicit none
        character(50) keyword
        character(120) file_line
        character*4          :: target(20), default_value(20)
        integer*4 keyword_line
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            target = default_value
            return
        end if
        read(file_lines(keyword_line), '(20a4)') target
        return
    end subroutine readStringKeywordValue !****************************************************************************

    !> Efetua a leitura de uma palavra-chave to tipo real. Se nao encontrado, associa o valor defualt fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde o real sera atribuido.
    !! @param default_value Valor default.
    subroutine readRealKeywordValue(keyword, target, default_value)
        implicit none
        character(50) keyword
        character(120) file_line
        real(8) target, default_value
        integer*4 keyword_line
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            target = default_value
            return
        end if
        file_line = adjustL(trim(file_lines(keyword_line)))
        read(file_line, *) target
        return
    end subroutine readRealKeywordValue !****************************************************************************

    !> Efetua a leitura de valores definidos para propriedades de sa�da. Na pr�tica s�o lidas 3 vari�veis.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        A vari�vel destino.
    !! @param var_out       O valor da vari�vel out.
    !! @param var_n         Valor n.
    subroutine readOutFlagKeyword(keyword, target, var_out, var_n)
        implicit none
        character(50) keyword
        integer*4 target, var_n
        character(120) var_out
        character(120) file_line
        integer*4 keyword_line
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            return
        end if
        read(file_lines(keyword_line), *) target
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line), "(a)") var_out
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line), *) var_n
    end subroutine readOutFlagKeyword


    !> Efetua a geracao de coordeadas, de acordo com parametros.
    !!
    !! @param a      Matriz onde serao armazenados os dados.
    !! @param nra    Inteiro indicando nra
    !! @param nLinhaArqInput    Indice da linha onde as coordenadas estao posicionadas no array linhas no arquivo de entrada.
    subroutine genflDS(a,nra, nLinhaArqInput )
      use mMalha, only: genfl1   
      use mGlobaisEscalares

      implicit none

      integer*4:: nra, nLinhaArqInput
      real*8  :: a(nra,*)
      real*8  :: temp(6,20)
      integer*4:: n,numgp,ninc(3),inc(3)
      integer*4:: i, j, m, mgen

      100 continue
      read(file_lines(nLinhaArqInput:),1000) n,numgp,(temp(i,1),i=1,nra)
      nLinhaArqInput = nLinhaArqInput + 1

      if (n.eq.0) return
     ! call move(a(1,n),temp,nra)
      a(1:nra,n) = temp(1:nra,1)

      if (numgp.ne.0) then
         do 200 j=2,numgp
         read(file_lines(nLinhaArqInput:),1000) m,mgen,(temp(i,j),i=1,nra)
         nLinhaArqInput = nLinhaArqInput + 1

        if (mgen.ne.0) temp(1:nra,j)=a(1:nra,m)         !B
         200    continue
         read(file_lines(nLinhaArqInput:),2000) (ninc(i),inc(i),i=1,3)
         nLinhaArqInput = nLinhaArqInput + 1

         call genfl1(a,nra, temp, n, numgp, ninc, inc)
      endif
      go to 100

      1000 format(2i10,6f10.0)
      2000 format(16i10)

    end  subroutine !**************************************************************************************************


    !> Subrotina para ler e gerar dados nodais inteiros.
    !! @param ia              Array de entrada.
    !! @param m               Numero de linhas na matriz de entrada.
    !! @param nLinhaArqInput  Indice da linha no array de linhas do arquivo de entrada.
    subroutine igenDS(ia, m, nLinhaArqInput)
        use mGlobaisEscalares

        integer*4:: m, ia(m,*), nLinhaArqInput
        integer*4:: ib(m)
        integer*4:: n, ne, ng
        integer*4:: i
        !
        100 continue
        read(file_lines(nLinhaArqInput:),1000) n,ne,ng,(ib(i),i=1,m)
        nLinhaArqInput = nLinhaArqInput + 1

        if (n.eq.0) return

        if (ng.eq.0) then
            ne = n
            ng = 1
        else
            ne = ne - mod(ne-n,ng)
        endif
        !
        do 200 i=n,ne,ng
        ia(:,i)=ib
        200 continue
        !
        go to 100
        !
        1000 format(16i10)
    end subroutine igenDS !********************************************************************************************
    !> Subrotina respons�vel por ler e gerar conectividades nodais e ladais.
    !> @param keyword_name  O nome da keyword associada.
    !> @param conecElem     C�digo do elemento
    !> @param mat           C�digo do material
    !> @param nen           N�mero de elementos.
    subroutine leituraGeracaoConectividadesDs(keyword_name, conecElem, mat, nen)
        use mLeituraEscrita , only: genel1
        integer*4:: n,nel(3),incel(3),inc(3)

        integer*4:: nen
        integer*4:: conecElem(nen,*),mat(*) 
        integer*4:: itemp(27) !B
        character(len=50) keyword_name
        !
        integer*4:: m,ng,i, keyword_line

!        write(*,'(a)') " em subroutine leituraGeracaoConectividadesDS(keyword_name, conecElem, mat, nen)"
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) return

        100 continue
        read(file_lines(keyword_line:),1000) n,m,(itemp(i),i=1,nen),ng
        !write(*,1000) n,m,(itemp(i),i=1,nen),ng
        keyword_line = keyword_line + 1

        if (n.eq.0) return
        !call imove(conecElem(1,n),itemp,nen)
        conecElem(1:nen,n)=itemp(1:nen)
        mat(n)=m
        if (ng.ne.0) then
            !Generate data
            read(file_lines(keyword_line:),1000) (nel(i),incel(i),inc(i),i=1,3)
            keyword_line = keyword_line + 1
            !write(*,1000) (nel(i),incel(i),inc(i),i=1,3)
            call genel1(conecElem,mat,nen,n,nel,incel,inc)
        endif
        go to 100
        1000 format(16i10,10x,14i10)
    end subroutine !***************************************************************************************************


    !> Subrotina respons�vel por ler e gerar elementos de face.
    !> @param keyword_name     A plavra-chave associada.
    !> @param conecElem        C�digo dos n�s dos elementos.
    !> @param nen              N�mero de lementos.
    !> @param nelx             N�mero de elmentos em x.
    !> @param nely             N�mero de elmentos em y.
    !> @param nelz             N�mero de elmentos em z.
    subroutine genelFacesDS(keyword_name, conecElem, nen, nelx, nely, nelz)
        use mGlobaisEscalares
        use mMalha, only: numel_bm

        implicit none
        integer*4:: nen, nelx, nely, nelz
        integer*4:: conecElem(nen,*)
        character(len=50) keyword_name

        integer*4:: ng, n, m, nel, i
        integer*4:: condicao, condicao2, keyword_line

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) return

        read(file_lines(keyword_line:),1000) n,m,(conecElem(i,1),i=1,nen),ng
        keyword_line = keyword_line + 1
        !write(*,1000) n,m,(conecElem(i,1),i=1,nen),ng

        condicao=0
        condicao2=0

        do nel=2, numel_bm
            if(condicao==0.and.condicao2==0) then
                do i=1, nen
                    conecElem(i,nel)=conecElem(i,nel-1)+1
                end do
            else
                if(condicao==1.and.condicao2==0) then
                    do i=1, nen
                        if(i<=4) then
                            conecElem(i,nel)=conecElem(i,nel-1)+(nelx+1)+1
                        else
                            conecElem(i,nel)=conecElem(i,nel-1)+1
                        endif
                    end do
                else
                    if(condicao==1.and.condicao2==1) then
                        do i=1, nen
                            if(i<=4) then
                                conecElem(i,nel)=conecElem(i,nel-1)+(nelx+1)+nelx+(nelx*nely)+1
                            else
                                conecElem(i,nel)=conecElem(i,nel-1)+(nelx*(nely+1))+(nely*(nelx+1))+1
                            endif
                        enddo
                    end if
                end if
            end if
            if(mod(nel, nelx)==0) then
                condicao=1
            else
                condicao=0
            end if
            if(mod(nel, nelx*nely)==0) then
                condicao2=1
            else
                condicao2=0
            end if
        end do
        1000 format(16i10,10x,14i10)
    end subroutine !***************************************************************************************************


    subroutine readNodeElementsDS
        use mGlobaisEscalares
        use mMalha,          only: nen_BM
        implicit none
        character(len=50) keyword_name
       ! integer*4 :: ntype, numat

        keyword_name = "ntype_pc"
        call readIntegerKeywordValue(keyword_name, ntype, ntype)
        keyword_name = "numat_pc"
        call readIntegerKeywordValue(keyword_name, numat_BM, numat_BM)
        keyword_name = "nen_pc"
        call readIntegerKeywordValue(keyword_name, nen_BM, nen_BM)
        keyword_name = "nicode_pc"
        call readIntegerKeywordValue(keyword_name, nicode_BM, nicode_BM)
    end subroutine !***************************************************************************************************


    !> Efetua a leitura de propriedades de materiais.
    !> @param keyword_name  Keyword especifica das  propriedades de materiais.
    subroutine readMaterialPropertiesDS(keyword_name)
        use mGlobaisEscalares
        use mGlobaisArranjos
        use mleituraEscrita, only: iecho

        implicit none
        character(len=50) keyword_name
        integer*4 n, m, i, keyword_line

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) return

        do 400 n=1,numat_BM
        if (mod(n,50).eq.1) write(iecho,4000) numat_BM
        read (file_lines(keyword_line:),  5000) m,(c_BM(i,m),i=1,3)
        keyword_line = keyword_line + 1
        write(iecho,6000) m,(c_BM(i,m),i=1,3)
        400 continue
        5000  format(i10,5x,5f10.0)
        6000  format(2x,i3,1x,5(1x,1pe11.4))
        4000  format(///,&
                ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
                ' number of material sets . . . . . . . . . . (numat ) = ',i10///,2x,'set',4x,'Kx ',&
                10x,'Ky',10x,'Kz')
    end subroutine readMaterialPropertiesDS !**************************************************************************

    !> Faz a leitura de constant body forces.
    !> @param keyword_name  Keyword especifica para constant body forces.
    subroutine readConstantBodyForcesDS(keyword_name)
        use mGlobaisArranjos, only: grav_BM
        use mleituraEscrita, only: iecho

        implicit none
        character(len=50) keyword_name
        integer*4 i, keyword_line

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) return

        read  (file_lines(keyword_line:),  7000) (grav_BM(i),i=1,3)
        write (iecho,8000) (grav_BM(i),i=1,3)

        7000 format(8 f10.0)
        8000 format(///,&
         ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
         ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
         ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
         ' exemplo 3............................ = ',      1pe15.8,//)
    end subroutine readConstantBodyForcesDS !**************************************************************************



end module mInputReader

