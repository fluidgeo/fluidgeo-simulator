! 
!         programa de elementos finitos em fortran 90 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!
!         Eduardo Garcia e Tuane Lopes
!         bidu@lncc.br, tuane@lncc.br
!
!         LNCC/MCT
!         Petropolis, 07.2013
! 
!     ************************************************************
!     *                                                          *
!     *                                                          *
!     *         A LINEAR STATIC FINITE ELEMENT PROGRAM FOR       *
!     *                                                          *
!     *                 GALERKIN METHOD                          *
!     *                                                          *
!     *                                                          *
!     ************************************************************
!
      module mGlobaisArranjos

        integer, allocatable :: npar_B(:)
        integer, allocatable :: npar_BM(:)
        real*8               :: etime(6), phi_range(2)
        character*4          :: title_B(20)
        character*4          :: title_BM(20)
        character(len=80)    :: reservoir_case, coupling_mode
        character(len=80)    :: reservoir_depth

        real*8,  allocatable :: uTempoN(:)

        integer, allocatable :: mat_BM(:)
        real*8, allocatable  :: grav_BM(:), bf_BM(:,:), c_BM(:,:), celast(:)

        integer, allocatable :: mat_B(:)
        real*8, allocatable  :: grav_B(:), bf_B(:,:), c_B(:,:), K_BM(:)
        real*8, allocatable  :: Knp_BM(:)
        REAL*8, allocatable  :: phi_n(:), phi_n0(:), trEps(:), trEpsTmpAnt(:)
        
      end module mGlobaisArranjos


      module mGlobaisEscalares

        character (LEN=2) :: dimModelo
        integer :: ntype
        integer :: exec,iprtin
        integer :: numParElem=15
        real*8, parameter  :: zero=0.0d0, one=1.0d0, two=2.0d0, three=3.0d0
        real*8, parameter  :: four=4.0d0, five=5.0d0, six=6.0d0
        real*8, parameter  :: pt1667=0.1666666666666667d0, pt25=0.25d0, pt5=0.5d0
        real*8  :: coef
        integer :: nRK, ordemRK
        integer :: optCC
        character(len=10) :: optSolver
        logical :: simetriaVel, simetriaGeo, random_porosity

        integer :: numat_BM
        integer :: nrowsh_BM,nicode_BM,npint_BM

        integer :: numat_B
        integer :: nrowsh_B,nicode_B,npint_B
        
        integer*4:: ndofD, nlvectD

      integer :: nvel,nnp,ns
      integer :: NITGEO, NCREEP, NLOOPS, NUMDX
      real(8) :: tzero,tTransporte,tc,tt,dtBlocoTransp
      real(8) :: tempoNucTrans
      real*8  :: ttv, ttp, tts, tempoSolverVel
      real*8  :: tmVel, tmGeo, tsGeo
      real*8  :: tgeoFase1, tgeoFase2, tgeoFase3, tgeoFase4, ttgeo
      real*8  :: lambda, mu

      integer :: geomech

      end module mGlobaisEscalares

