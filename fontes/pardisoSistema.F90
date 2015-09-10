  module mPardisoSistema

        use mCSR

        integer*4 pt(64), iparm(64)
        REAL*8  dparm(64)

        logical :: primeiravezVel

      contains

#ifdef withPardiso

      subroutine solverPardisoPPD_Nodal(ia, ja, a, b, neq, nonzeros, simetria, label, parte)

      implicit none
!
        integer*4, INTENT(IN)  :: NEQ, NONZEROS
        integer*4, intent(in)  :: ia(neq+1), ja(nonzeros)
        REAL*8, INTENT(INOUT)   :: a(nonzeros)
        REAL*8, INTENT(INOUT):: b(neq)
        LOGICAL, INTENT(IN)  :: simetria
        character(LEN=*), INTENT(IN) :: label
        character(LEN=4), INTENT(IN) :: parte
!
        REAL*8, allocatable        :: x(:)
        REAL*8,  save :: ddum
        integer*4, save :: maxfct, mnum, mtype, phase, nrhs, error,  msglvl
        integer*4, save :: idum, solver
        integer*4 :: i
        REAL*8  :: t1, t2 , tt1, tt2
        integer*4:: omp_get_num_threads

        !integer*4:: ptAux(64), iparmAux(64)
! 
!
!  .. Setup Pardiso control parameters und initialize the solvers     
!     internal adress pointers. This is only necessary for the FIRST   
!     call of the PARDISO solver.                                     
!     
! $OMP PARALLEL

      allocate ( x(neq)) 
      !stop

!      write(*,*) "em, subroutine solverPardisoPPD_Nodal(ia, ja, a, b, neq, nonzeros, s"

      if(parte=='fact'.or. parte=='full') then

      iparm=0
      mtype     =  2   ! real and symmetric matrix, positive definite
      mtype     = -2   ! real and symmetric matrix, indefinite
      if(simetria.eqv..false.)  mtype     = 11   ! unsymmetric matrix, indefinite

      iparm(2)  = 1    !ou 0?  !Fill-In reduction reordering.
      iparm(11) = 0    ! Do not use (symmetric matrices).
      iparm(16) = 0    ! 

      solver     = 0    ! use 0 for sparse direct method or 1 multi-recursive iterative solver
      msglvl     = 0    ! with statistical information

      iparm(33) = 0    ! compute determinant 
      iparm(52) = 1    !For OpenMP-threaded solver
      iparm(27) = 1
      
      nrhs      = 1
      mnum      = 1
      pt        = 0
      idum      = 0
      ddum      = 0.0
      dparm    = 0.0
      maxfct    = 1
      error     = 0
!
!  .. Numbers of Processors ( value of OMP_NUM_THREADS )
!
       iparm(3) = 1
!   experimentando, bidu,  atribui valores default a todo vetor iparm
       !iparm(1) = 0
! 
#ifdef withOMP
!$OMP PARALLEL shared(iparm)  
       !write(*,*) " omp_get_num_threads() = ", omp_get_num_threads()
       iparm(3) =   omp_get_num_threads()
!$OMP END PARALLEL
#endif
       write(*,'(3a,i3)')  "++, em solverPardisoPPD_Nodal,",label,", com numThreads:", iparm(3) 
!        
       call timing(tt1)
!
!  .. PARDISO license check and initialize solver
!      call pardisoinit(pt, mtype, solver, iparm, dparm, error)
!
      IF (error .NE. 0) THEN
        IF (error.EQ.-10 ) WRITE(*,*) 'No license file found'
        IF (error.EQ.-11 ) WRITE(*,*) 'License is expired'
        IF (error.EQ.-12 ) WRITE(*,*) 'Wrong username or hostname'
        STOP ' (error .NE. 0) '
      ELSE
        WRITE(*,*) '[PARDISO]: License check was successful ... '
      END IF
!         phase     = -1     ! Release all internal memory for all matrices 
!       CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
!                    idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!..   Reordering and Symbolic Factorization, This step also allocates
!     all memory that is necessary for the factorization
!
#ifdef DEBUG
       WRITE(*,*) 'Begining reordering ... ', label
#endif
         call timing(t1)
         phase     = 11     ! only reordering and symbolic factorization

         CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                    idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm)
         call timing(t2)

#ifdef mostrarTempos
      write(*,*) "reordering: ", label, ", tempo de parede = ", t2 - t1
#endif

      IF (error .NE. 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', error
        STOP
      END IF

       WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
       WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

!
!.. Factorization.
!
#ifdef DEBUG
     WRITE(*,*) 'Begining factorization  ... '
#endif
       call timing(t1)
      phase     = 22  ! only factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, a, ia, ja,  &
                    idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm) 
      call timing(t2)
#ifdef mostrarTempos
      write(*,*) "factorization: ", label, ", tempo = ", t2 - t1
#endif

    IF (error .NE. 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
      STOP
    ENDIF 

      endif  !     if(parte=='fact'.or. parte=='full') then
! 
      if(parte=='back'.or. parte=='full') then

!.. Back substitution and iterative refinement
#ifdef DEBUG
      WRITE(*,*) 'Begining backsubstitution  ... '
#endif
       call timing(t1)
      iparm(8) = 1   ! max numbers of iterative refinement steps
      phase     = 33  ! only solve
      CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                   idum, nrhs, iparm, msglvl, b, x, error, dparm) 
       call timing(t2)
#ifdef mostrarTempos
      write(*,*) "backsubstitution: ",label, ", tempo = ", t2 - t1
#endif

       call timing(tt2)
#ifdef mostrarTempos
      WRITE(*,*) 'Solve completed ...  ',label, ", tempo =", tt2-tt1
#endif

      b(1:neq) = x(1:neq)

      endif !if(parte=='back'.or. parte=='full') then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if(parte=='full') then

!       desalocando memoria
#ifdef DEBUG
      WRITE(*,*) 'Begining memory desalocation  ... '
#endif
      call timing(t1)
      phase = -1

      CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                     idum, nrhs, iparm, msglvl, b, x, error, dparm)
         
      call timing(t2)
#ifdef mostrarTempos
      write(*,*) "memory desalocation: ", label, ", tempo = ", t2 - t1
#endif

      endif !if(parte=='full') then
      deallocate ( x) 


      end subroutine solverPardisoPPD_Nodal

#endif
!
!=======================================================================
!    
      subroutine solverPardisoEsparso(alhs, b, Ap, Ai, nonzeros, neq, simetria, parte, label)

      implicit none 
!
      integer*4, intent(in) :: nonzeros, neq
      real*8              :: alhs(nonZeros), b(neq)
      integer*4, intent(in) :: Ap(neq+1), Ai(nonZeros)
      logical, intent(in) :: simetria
      character(LEN=4), intent(in) :: parte
      character(LEN=*), intent(in) :: label
! 
#ifdef withPardiso
      call solverPardisoPPD_Nodal(Ap, Ai, alhs, b, neq, nonzeros, simetria, label, parte)
#endif
!      write(*,'(5f10.5)') b(1:5)
!      write(*,'(5f10.5)') b(neq-4: neq)

      end subroutine solverPardisoEsparso
!
! **** new *********************************************************************
!
       subroutine addlhsCSR(alhs,eleffm,lm,Ap,Ai,nee)
! 
! .... program to add element left-hand-side matrix to
!         global left-hand-side matrix
! 
! 
       implicit none
! 
! .... deactivate above card(s) for single-precision operation
! 
       integer*4:: nee
       real*8  :: alhs(*),eleffm(nee,*)
       integer*4:: lm(*)
       integer*4::  Ap(*), Ai(*)
!
       integer*4:: i, j, linha, coluna
       integer*4:: inicio, fim, jj
       integer*4, save ::  nel = 0
! 
!      nel = nel + 1 ! ? por que nao passar o parametro nel?
                     !  resposta: esta é uma variavel auxiliar
       do j=1,nee
          coluna = lm(j)
          if (coluna.gt.0) then
             do i=1,nee
                linha = lm(i)
                if (linha.gt.0) then
                   inicio=Ap(coluna)
                      fim=Ap(coluna+1)-1
                   do jj=inicio, fim
                      if(Ai(jj)==linha) then
              !          write(500 , '(3i3,e15.5)', advance='NO') nel, linha, jj, alhs(jj)
                         alhs(jj)= alhs(jj)+eleffm(i,j)
              !          write(500 ,'(2e15.5)')  eleffm(i,j), alhs(jj)
                      endif
                   enddo
                endif
             enddo
          endif
       enddo
! 
       return
       end subroutine
!     
     end module mPardisoSistema
