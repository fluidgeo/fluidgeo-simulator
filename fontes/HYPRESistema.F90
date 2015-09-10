 module mHYPREsistema
      integer*4 :: myid, num_procs
      integer*4 :: mpi_comm

#ifdef withHYPRE
      include 'mpif.h'
#endif

      integer*4 :: posPonteiro, contPonteiro, posColunas, posCoef
      integer*4 :: lda, nonzeros, nonzerosEst
      logical   :: primeiravezVel, primeiravezGeo
      integer*8, parameter :: HYPRE_PARCSR=5555

      contains

      subroutine inicializarMPI(myid_, num_procs_, mpi_comm_)
!
      integer*4, intent(inout) :: myid_, num_procs_
      integer*4  :: mpi_comm_
      integer*4 ::  ierr  
#ifdef withHYPRE
      print *, " em subroutine inicializarMPI(myid_, num_procs_, mpi_comm_)"
      !mpi_comm_ = 1000 ! MPI_COMM_WORLD
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid,      ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
      print *, " em subroutine inicializar, myid=", myid_,", numprocs=", num_procs_,", mpi_comm=", mpi_comm_
!      stop
!   Convert the Fortran MPI communicator to a C version that can be used in hypre.
!   Uncomment the line below if you are using MPICH
     mpi_comm = MPI_COMM_WORLD
       write(*,'(a,i12)') "i+++, mpi_comm_ =" ,mpi_comm_
       write(*,'(a,i12)') "i+++, MPI_COMM_WORLD =" ,MPI_COMM_WORLD
!   Uncomment the line below if you are using LAM-MPI
      !call HYPRE_MPI_Comm_f2c(mpi_comm, MPI_COMM_WORLD, ierr)
#endif
      end subroutine inicializarMPI
!
!=============================================================================
!

      subroutine criarSistemaAlgHYPRE(A_, parcsr_A_, b_, par_b_, u_, par_u_, &
                 solver_, Clower_, Cupper_, mpi_comm_)

      integer*8, intent(out) :: A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      integer*8, intent(in)  :: Clower_, Cupper_
      integer*4, intent(in)  :: mpi_comm_

      integer*8 :: ierr 
!     Create the matrix.

!     Note that this is a square matrix, so we indicate the row partition
!     size twice (since number of rows = number of cols)

#ifdef withHYPRE
      write(*,*) "  "
      write(*,*) " em subroutine criarSistemaAlgHYPRE(A_, parcsr_A_, b_, par_b_, ..."
      A_=0; parcsr_A_=0; b_=0; par_b_=0; u_=0; par_u_=0; solver_=0;
       write(*,'(a,7i15)') "c+++, ", A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_

       write(*,'(a,i12)') "c+++, mpi_comm =" ,mpi_comm
      call HYPRE_IJMatrixCreate     (mpi_comm, Clower_, Cupper_, Clower_, Cupper_, A_, ierr )
!     Choose a parallel csr format storage (see the User's Manual)
      call HYPRE_IJMatrixSetObjectType (A_, HYPRE_PARCSR,  ierr)
!      call HYPRE_StructMatrixSetSymmetric (A_, 1)
!     Initialize before setting coefficients
      call HYPRE_IJMatrixInitialize    (A_, ierr)
       
!     Create the rhs and solution
      call HYPRE_IJVectorCreate     (mpi_comm_world, Clower_, Cupper_, b_, ierr )
      call HYPRE_IJVectorSetObjectType (b_, HYPRE_PARCSR,  ierr)
      call HYPRE_IJVectorInitialize    (b_, ierr)
!
      call HYPRE_IJVectorCreate     (mpi_comm_world, Clower_, Cupper_, u_, ierr )
      call HYPRE_IJVectorSetObjectType (u_, HYPRE_PARCSR,  ierr)
      call HYPRE_IJVectorInitialize    (u_, ierr)
       write(*,'(a,7i15)') "c+++, ", A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_

#endif
      end subroutine criarSistemaAlgHYPRE

      subroutine criarB_HYPRE(b_, Clower_, Cupper_, mpi_comm_)

      integer*8, intent(out) :: b_
      integer*8, intent(in)  :: Clower_, Cupper_
      integer*4, intent(in)  :: mpi_comm_
      integer*8 :: ierr 

#ifdef withHYPRE
      write(*,*) "  "
      write(*,*) " em subroutine criarB_HYPRE"
      write(*,'(a,i10)') "cB+++, ", b_
      call HYPRE_IJVectorCreate        (mpi_comm_world, Clower_, Cupper_, b_, ierr )
      call HYPRE_IJVectorSetObjectType (b_, HYPRE_PARCSR, ierr)
      call HYPRE_IJVectorInitialize    (b_, ierr)
      write(*,'(a,i10)') "cB+++, ", b_
#endif
      end subroutine criarB_HYPRE

      subroutine criarU_HYPRE(u_, Clower_, Cupper_, mpi_comm_)

      integer*8, intent(out) :: u_
      integer*8, intent(in)  :: Clower_, Cupper_
      integer*4, intent(in)  :: mpi_comm_
      integer*8 :: ierr 

#ifdef withHYPRE
      write(*,*) "  "
      write(*,*) " em subroutine criarU_HYPRE"
      write(*,'(a,i10)') "cU+++, ", u_
      call HYPRE_IJVectorCreate        (mpi_comm_world, Clower_, Cupper_, u_, ierr )
      call HYPRE_IJVectorSetObjectType (u_, HYPRE_PARCSR, ierr)
      call HYPRE_IJVectorInitialize    (u_, ierr)
      write(*,'(a,i10)') "cU+++, ", u_
#endif
      end subroutine criarU_HYPRE

!
!=============================================================================
!
subroutine finalizarMontagemSistemaHYPRE(A_, parcsr_A_, b_, par_b_, u_, par_u_)

      implicit none
      integer*8, intent(in)  :: A_, b_, u_
      integer*8, intent(out) :: parcsr_A_, par_b_, par_u_
      integer*8 :: ierr

#ifdef withHYPRE
      write(*,*) " finalizarMontagemSistemaHYPRE(A_, parcsr_A_, b_, par_b_, u_, par_u_)"
      write(*,'(a,7i15)') "f+++, ", A_, parcsr_A_, b_, par_b_, u_, par_u_
!     Assemble after setting the coefficients
      call HYPRE_IJMatrixAssemble( A_, ierr)
!     Get parcsr matrix object
      call HYPRE_IJMatrixGetObject( A_, parcsr_A_, ierr)

      call HYPRE_IJVectorAssemble  (b_, ierr)
      call HYPRE_IJVectorGetObject (b_, par_b_, ierr)

      call HYPRE_IJVectorAssemble  (u_, ierr)
      call HYPRE_IJVectorGetObject (u_, par_u_, ierr)
      write(*,'(a,7i15)') "f+++, ", A_, parcsr_A_, b_, par_b_, u_, par_u_
#endif

end subroutine finalizarMontagemSistemaHYPRE
!
!=============================================================================
!
      subroutine resolverSistemaAlgHYPRE  (A_, parcsr_A_, b_, par_b_, u_, par_u_, &
                                          initialGuess_, solver_, &
                                          solver_id_, precond_id_, tol_, &
                                          num_iterations_, final_res_norm_, &
                                          myid_, mpi_comm_, brhs_, rows_,  neq_) 

      implicit none
      integer*8, intent(in)         :: A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      integer*4, intent(in)         :: solver_id_, precond_id_
      double precision,  intent(in) :: tol_ ! = 1.0e-08 
      integer*4, intent(out)        :: num_iterations_
      double precision, intent(out) :: final_res_norm_
      integer*4, intent(in)         :: myid_
      integer*4, intent(in)         :: mpi_comm_ 
      integer*4, intent(in)         :: neq_
      integer*4                     :: rows_(neq_)
      double precision              :: brhs_(neq_),  initialGuess_(neq_)

      integer*8 ::  precond
      integer*8        :: precond_id;
      integer*8        :: ierr
      integer*8, save  :: contador = 1
      integer*8        :: local_size , numMaxIterations = 10000
      integer*4        :: printLevel, solver_id
      integer*8        :: i, Flower, Fupper
      real*8 :: t1, t2, t3, t4, elapsedT 
      local_size = neq_
      printLevel = 0

#ifdef withHYPRE
      write(*,*) "em subroutine resolverSistemaAlgHYPRE  (A_, parcsr_A_, b_, par_b_, "
      write(*,'(a,7i15)') "r+++, ", A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      !block
      !character(len=9)  :: nomeU = "pot.out.x"
      !character(len=10) :: nomeB = "brhs.out.x"
      !character(len=10) :: nomeA = "alhs.out.x"
      !call HYPRE_IJVectorPrint( u_, nomeU, ierr)
      !call HYPRE_IJVectorPrint( b_, nomeB, ierr)
      !call HYPRE_IJMatrixPrint( A_, nomeA, ierr)
      !end block
      call HYPRE_IJVectorSetValues (u_, neq_, rows_, initialGuess_, ierr)
      call HYPRE_IJVectorSetValues (b_, neq_, rows_, brhs_, ierr)

      call timing(t1)
      if     (solver_id_==0) then
!     AMG
        write(*,*) " ... BoomerAMGCreate , BD"
        call solverOptionA ()

      elseif (solver_id_==1) then
        write(*,*) " ... PCG with AMG preconditioner, BD";
       ! if(precond_id_==2)tol=1.0e-10 !test and error experience chosen, better to potential
       ! if(precond_id_==1)tol=1.0e-08 !test and error experience chosen, better to post-processed flux
        write(*,*) "++ precond_id_ =", precond_id_,", tol =", tol_
        call solverOptionB ()
!     PCG with ParaSails
      elseif (solver_id_==2) then
        write(*,*) " ... PCG with ParaSails"
        call solverOptionC ()
      else
         if (myid_ .eq. 0) then 
           print *,'Invalid solver id specified'
           stop
         endif  
      endif
      call timing(t2)
      elapsedT = t2 - t1
      ! colocando valores no vetor de trabalho externo
      call HYPRE_IJVectorGetValues (u_, local_size, rows_, brhs_, ierr)


      !block
      !  integer*4 :: printSol=0
      call escreverResultadosHYPRE (u_, num_iterations_, final_res_norm_, elapsedT, myId, 0) 
      !end block
      write(*,'(a,7i15)') "r+++, ", A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_

      contains

      subroutine solverOptionA
!        Create solver
         call HYPRE_BoomerAMGCreate          (solver_, ierr)

!        Set some parameters (See Reference Manual for more parameters)

!        print solve info + parameters 
         call HYPRE_BoomerAMGSetPrintLevel   (solver_, 3,      ierr)  
!        Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType  (solver_, 6,      ierr) 
!        G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType    (solver_, 6,      ierr)     
!        Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps    (solver_, 1,      ierr)  
!        maximum number of levels 
         call HYPRE_BoomerAMGSetMaxLevels    (solver_, 200,     ierr) 
!        conv. tolerance
         call HYPRE_BoomerAMGSetTol          (solver_, tol_, ierr)    

!        Now setup and solve!
         call HYPRE_BoomerAMGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr )
         call HYPRE_BoomerAMGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr )

!        Run info - needed logging turned on 
         call HYPRE_BoomerAMGGetNumIterations(solver_, num_iterations_, ierr)
         call HYPRE_BoomerAMGGetFinalReltvRes(solver_, final_res_norm_, ierr)

!        Destroy solver_
         call HYPRE_BoomerAMGDestroy         (solver_, ierr )

!     PCG with AMG preconditioner
      end subroutine solverOptionA

      subroutine solverOptionB ()
!        Create solver_
         call HYPRE_ParCSRPCGCreate        (MPI_COMM_WORLD, solver_, ierr)

!        Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter    (solver_, numMaxIterations,   ierr)
         call HYPRE_ParCSRPCGSetTol        (solver_, tol_, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (solver_, 1,      ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (solver_, 2,      ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (solver_, printLevel,      ierr)

!        Now set up the AMG preconditioner and specify any parameters

         call HYPRE_BoomerAMGCreate          (precond, ierr)

!        Set some parameters (See Reference Manual for more parameters)

!  Relaxation Parameters:
!   Visiting Grid:                     down   up  coarse
!            Number of sweeps:            1    1     1 
!   Type 0=Jac, 3=hGS, 6=hSGS, 9=GE:      6    6     6 

!        print less solver_ info since a preconditioner
         call HYPRE_BoomerAMGSetPrintLevel   (precond, 1,     ierr); 
!        Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType  (precond, 6,     ierr) 
!        SYMMETRIC G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType    (precond, 6,     ierr)     
!        Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps    (precond, 1,     ierr)  
!        conv. tolerance
         call HYPRE_BoomerAMGSetTol          (precond, 0.0d-6, ierr)     
!        do only one iteration! 
         call HYPRE_BoomerAMGSetMaxIter      (precond, 1,     ierr)

!        set amg as the pcg preconditioner
         call HYPRE_ParCSRPCGSetPrecond      (solver_, precond_id_, precond, ierr)

!        Now setup and solve!
         call HYPRE_ParCSRPCGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr)
!        Run info - needed logging turned on 
        call HYPRE_ParCSRPCGGetNumIterations (solver_, num_iterations_, ierr)
        call HYPRE_ParCSRPCGGetFinalRelative (solver_, final_res_norm_, ierr)

!       Destroy precond and solver
        call HYPRE_BoomerAMGDestroy          (precond, ierr )
        call HYPRE_ParCSRPCGDestroy          (solver_, ierr)
      end subroutine solverOptionB

      subroutine solverOptionC()

!        Create solver
         call HYPRE_ParCSRPCGCreate          (MPI_COMM_WORLD, solver_, ierr)

!        Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter      (solver_, 1000, ierr)
         call HYPRE_ParCSRPCGSetTol          (solver_, 1.0d-7, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm      (solver_, 1, ierr)
         call HYPRE_ParCSRPCGSetPrintLevel   (solver_, printLevel, ierr)
         call HYPRE_ParCSRPCGSetLogging      (solver_, 1, ierr)

!        Now set up the Parasails preconditioner and specify any parameters
         call HYPRE_ParaSailsCreate          (MPI_COMM_WORLD, precond,ierr)
         call HYPRE_ParaSailsSetParams       (precond, 0.1d0, 1, ierr)
         call HYPRE_ParaSailsSetFilter       (precond, 0.05d0, ierr)
         call HYPRE_ParaSailsSetSym          (precond, 1)
         call HYPRE_ParaSailsSetLogging      (precond, 3, ierr)

!        set parsails as the pcg preconditioner
         !precond_id_ = 4
         call HYPRE_ParCSRPCGSetPrecond      (solver_, precond_id_, precond, ierr)


!        Now setup and solve!
         call HYPRE_ParCSRPCGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr)

!        Run info - needed logging turned on 

        call HYPRE_ParCSRPCGGetNumIterations (solver_, num_iterations_, ierr)
        call HYPRE_ParCSRPCGGetFinalRelative (solver_, final_res_norm_, ierr)

!       Destroy precond and solver
        call HYPRE_ParaSailsDestroy          (precond, ierr )
        call HYPRE_ParCSRPCGDestroy          (solver_, ierr)
      end subroutine solverOptionC

#endif

      end subroutine resolverSistemaAlgHYPRE
!
!
!=============================================================================
!  
      subroutine destruirB_HYPRE   ( b_)
      integer*8 , intent(in) ::  b_
      integer*8 ::  ierr  
#ifdef withHYPRE
      write(*,*) " ++, em destruirB_HYPRE   ( b_)"
      write(*,'(a,i10)') "dB+++, b=",  b_ 
         call HYPRE_IJVectorDestroy(b_, ierr)
      write(*,'(a,i10)') "dB+++, b=",  b_ 
#endif
      end subroutine destruirB_HYPRE       
!=============================================================================
!  
      subroutine destruirU_HYPRE   ( u_)
      integer*8 , intent(in) ::  u_
      integer*8 ::  ierr  
#ifdef withHYPRE
      write(*,*) " ++, em destruirU_HYPRE   (u_)"
      write(*,'(a,i10)') "dU+++, u=",  u_ 
         call HYPRE_IJVectorDestroy(u_, ierr)
      write(*,'(a,i10)') "dU+++, u=",  u_ 
#endif
      end subroutine destruirU_HYPRE       
!=============================================================================
      subroutine destruirSistemaAlgHYPRE   (A_, b_, u_)
      integer*8 , intent(inout) :: A_, b_, u_
      integer*8 ::  ierr  
#ifdef withHYPRE
      write(*,*) " ++, em destruirSistemaAlgHYPRE   (A_, b_, u_)"
       write(*,'(a,7i10)') "d+++, ", A_, b_, u_ 
         call HYPRE_IJMatrixDestroy(A_, ierr)
         call HYPRE_IJVectorDestroy(u_, ierr)
         call HYPRE_IJVectorDestroy(b_, ierr)
       write(*,'(a,i10,a,7i10)') "d+++, A=", A_,", b_=", b_, u_ 
!         A_=0; b_=0; u_=0;
#endif
      end subroutine destruirSistemaAlgHYPRE       
!
!=============================================================================
!
      subroutine escreverResultadosHYPRE (x_, num_iterations_, final_res_norm_, elapsedT_, myId_, print_solution_)

      implicit none
      integer*8   :: x_
      integer*4        num_iterations_
      integer*4, intent(in) :: myId_
      integer*4           :: print_solution_ 
!
      double precision final_res_norm_, elapsedT_

      integer*8 ::  ierr  
      character (LEN=16)  :: SolverStatus="incomplete"

#ifdef withHYPRE
    write(*,*) 'em subroutine escreverResultadosHYPRE'
    SolverStatus="complete"
    if(SolverStatus=="complete" .and. myid_ .eq. 0) then
        print *
        print *, "Final Relative Residual Norm = ", final_res_norm_
        print *, "Iterations                   = ", num_iterations_
        print *, 'Elapsed real time            = ', elapsedT_
        print *
    endif

!     Print the solution
      !print_solution_ = 0 
      if ( print_solution_ .ne. 0 ) then
         call HYPRE_IJVectorPrint( x_, "ij.out.x", ierr)
      endif
#endif
      end subroutine escreverResultadosHYPRE           

      subroutine finalizarMPI                 ()
      integer*8 :: ierr
#ifdef withHYPRE
         call MPI_Finalize(ierr)
#endif
      end subroutine finalizarMPI   

      subroutine addnslHYPRE(A_, eleffm,lm,nee,lsym)
!
!         program to add element left-hand-side matrix to          
!                global left-hand-side matrix                      
!                                                                  
!        ldiag = .true.,  add diagonal element matrix              
!                                                                  
!        ldiag = .false, then                                     
!        add full nonsymmetric element matrix                   
!                                                                  
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*8, intent(inout) :: A_
      integer*4:: nee
      real*8  :: eleffm(nee,*)
      integer*4:: lm(nee)
      logical :: lsym
!
      integer*4, parameter:: numCoefPorLinha=100 
      integer*4:: i,j,k,l, eq, c, nnz, ierr 
      integer*4::  cols(numCoefPorLinha)
      real*8  :: values(numCoefPorLinha) 
      integer, save :: contador = 0
#ifdef withHYPRE
!
!    write(*,*) '+++ em addnslHYPRE, A_ =', A_
   
      contador = contador + 1
     ! write(*,*) contador, "lm = ", lm(:)
         do 400 j=1,nee
            nnz = 0
            eq = abs(lm(j))
            if(eq > 0) then 
               do 200 i=1,nee
                   c = abs(lm(i))
                   if(c > 0) then 
                       nnz = nnz + 1
                       cols(nnz) = c
                     values(nnz) = eleffm(j,i)
                   end if
  200          continue
               cols = cols - 1
               eq = eq - 1
!              write(*,'(a,i3,a,10i3)') " eq = ", eq,  ", ", cols(1:nnz)
               call HYPRE_IJMatrixAddToValues(A_, 1, nnz, eq, cols, values, ierr)
            end if
  400    continue
#endif
!
      return
      end subroutine addnslHYPRE
!
!**** new **********************************************************************
      subroutine addrhsHYPREGeo (b_, brhs_, elresf,lm,nee)
!
!.... program to add element residual-force vector to
!        global right-hand-side vector
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: nee
      integer*8 , intent(inout) ::  b_
      real*8, intent(in)  :: brhs_(*)
      real*8  :: elresf(nee)
      integer*4:: lm(nee)
!
      integer*4:: k, j, nnz, ierr
      integer*4:: rows(nee+20), i
      real*8   :: values(nee+20)
      integer, save :: nel = 0
#ifdef withHYPRE
!
      !write(*,*) brhs_(20000:20100); write(*,*) " ++++"; stop
      values = 0.0
      nel = nel + 1
       write(1002,'(i5,i5,a,16i6)') nee, nel, "lm = ", lm(:)
      nnz = 0     
      do 100 j=1,nee
      k = lm(j)
      if (k.gt.0) then
            nnz = nnz + 1
            values(nnz) = elresf(j) + brhs_(k)  
            rows  (nnz) = k-1
      end if
  100 continue
!      rows = rows - 1
!     write(1002,'(a,i5,a,16i6)') " nel = ", nel,  ", ", rows(1:nnz)
!     write(1002,'(a,i5,a,16e15.4)') " nel = ", nel,", brhs_=",  brhs_(rows(1:nnz))
!     write(1002,'(a,i3,a,i9,16e15.4)') " nel = ", nel,  ", ", b_, values(1:nnz)
      !call HYPRE_IJVectorSetValues (b_, nnz, rows, elresf, ierr)
      call HYPRE_IJVectorAddToValues (b_, nnz, rows, values, ierr)
!

!     write(*,'(a,4i5)') "+++++  nee = ", nee
!     write(*,'(a,4i5)') "+++++",   lm(1:nee)
!     write(*,'(a,4f8.4)') "+++++", elresf(1:nee)
!     write(*,'(a,4i5)') "+++++",   rows(1:nnz)
!     write(*,'(a,4f8.4)') "+++++", values(1:nnz)
!     do i = 1, 12
!      rows (i) = i-1
!     enddo
!     call HYPRE_IJVectorGetValues (b_, 12, rows, values, ierr)
!      write(*,'(a,12i5)') "+++++",   rows(1:12)
!     write(*,'(a,12f8.4)') "+++++",   values(1:12)
      return
#endif
      end subroutine addrhsHYPREGeo
!**** new **********************************************************************
      subroutine addrhsHYPRE (b_, brhs_, elresf,lm,nee)
!
!.... program to add element residual-force vector to
!        global right-hand-side vector
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: nee
      integer*8 , intent(inout) ::  b_
      real*8, intent(in)  :: brhs_(*)
      real*8  :: elresf(nee)
      integer*4:: lm(nee)
!
      integer*4:: k, j, nnz, ierr
      integer*4:: rows(nee+20), i
      real*8   :: values(nee+20)
      integer, save :: nel = 0
#ifdef withHYPRE
!
      values = 0.0
      nel = nel + 1
!       write(*,'(i5,a,16i4)') nel, "lm = ", lm(:)
      nnz = 0     
      do 100 j=1,nee
      k = lm(j)
      if (k.gt.0) then
            nnz = nnz + 1
            values(nnz) = elresf(j) + brhs_(k)  
            rows  (nnz) = k-1
      end if
  100 continue
!      rows = rows - 1
!     write(*,'(a,i3,a,16i6)') " nel = ", nel,  ", ", rows(1:nnz)
!     write(*,'(a,i3,a,i9,16e15.4)') " nel = ", nel,  ", ", b_, values(1:nnz)
      !call HYPRE_IJVectorSetValues (b_, nnz, rows, elresf, ierr)
      call HYPRE_IJVectorAddToValues (b_, nnz, rows, values, ierr)
!

!     write(*,'(a,4i5)') "+++++  nee = ", nee
!     write(*,'(a,4i5)') "+++++",   lm(1:nee)
!     write(*,'(a,4f8.4)') "+++++", elresf(1:nee)
!     write(*,'(a,4i5)') "+++++",   rows(1:nnz)
!     write(*,'(a,4f8.4)') "+++++", values(1:nnz)
!     do i = 1, 12
!      rows (i) = i-1
!     enddo
!     call HYPRE_IJVectorGetValues (b_, 12, rows, values, ierr)
!      write(*,'(a,12i5)') "+++++",   rows(1:12)
!     write(*,'(a,12f8.4)') "+++++",   values(1:12)
      return
#endif
      end subroutine addrhsHYPRE

      subroutine lerValoresSistemaAlgHYPRE (alhs_ , Ap_, Ai_, rhs_, x_, rows_, &
                                   neq_, nonzerosT_, Clower_, Cupper_,  &
                                A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_,  &
                                 myid_, mpi_comm_) 

      implicit none
      integer*4, intent(in)  :: nonzerosT_
      integer*4, intent(in)  :: neq_
      integer*8, intent(in)  :: Clower_, Cupper_
      double precision       :: alhs_(nonZerosT_)  
      integer*4              :: Ap_(neq_), Ai_(nonZerosT_)  
      double precision       :: rhs_(neq_), x_(neq_)
      integer*4              :: rows_(neq_)  
      integer*8              :: A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      integer*4, intent(in)  :: myid_
      integer*4, intent(in)  :: mpi_comm_

      integer*8     :: nnz, i, j, k, eq, ierr, printLevel
      integer*8     :: local_size, Flower, Fupper
      integer*4     :: cols(90)
      integer*4     :: colsB(90), eqB
      double precision ::  values(90)
     
#ifdef withHYPRE
 !   write(*,*) " em subroutine lerValoresSistemaAlgCSR_HYPRE ( A_, parcsr_A_, b_, par_b_, u_, par_u_, .."
    print*, " atribuindo valores da matrix para o HYPRE_IJMatrix "
    Flower    = Clower_+1; Fupper    = Cupper_+1;
   !   HYPRE_IJMatrixRead( <filename>, MPI_COMM_WORLD, HYPRE_PARCSR, &A );
    local_size = Cupper_ - Clower_ + 1

       !HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD, HYPRE_PARCSR, &b ); 
    print*, " ++++ B atribuindo valores de RHS para o HYPRE_IJVector ", Flower, Fupper
      !call HYPRE_IJVectorSetValues (b_, local_size, rows_, rhs_, ierr )
   !   call HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD, par_b_, b_ ); 

      print*, " +++++ atribuindo valores de guess para o HYPRE_IJVector "
      !call HYPRE_IJVectorSetValues (u_, local_size, rows_, x_, ierr)
   !   call HYPRE_IJVectorRead( , MPI_COMM_WORLD, par_u_, u_ ); 

#endif
 end subroutine lerValoresSistemaAlgHYPRE

      subroutine atribuirValoresSistemaAlgHYPRE (alhs_ , Ap_, Ai_, rhs_, x_, rows_, &
                                   neq_, nonzerosT_, Clower_, Cupper_,  &
                                A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_,  &
                                 myid_, mpi_comm_) 

      implicit none
      integer*4, intent(in)   :: nonzerosT_
      integer*4, intent(in)   :: neq_
      integer*8, intent(in)   :: Clower_, Cupper_
      double precision        :: alhs_(nonZerosT_)  
      integer*4               :: Ap_(neq_), Ai_(nonZerosT_)  
      double precision        :: rhs_(neq_), x_(neq_)
      integer*4               :: rows_(neq_)  
      integer*8               :: A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      integer*4, intent(in)   :: myid_
      integer*4, intent(in)   :: mpi_comm_

      integer*8     :: nnz, i, j, k, eq, ierr, printLevel
      integer*8     :: local_size, Flower, Fupper
      integer*4     :: cols(90)
      integer*4     :: colsB(90), eqB
      double precision ::  values(90)
     
#ifdef withHYPRE
    write(*,*) " em subroutine atribuirValoresSistemaAlgCSR_HYPRE ( A_, parcsr_A_, b_, par_b_, u_, par_u_, .."
    print*, " atribuindo valores da matrix para o HYPRE_IJMatrix "
    write(*,*) " Clower_ = ",  Clower_, ", Cupper_=", Cupper_
    Flower    = Clower_+1; Fupper    = Cupper_+1;

    do i = Flower, Fupper
    !   nnz    =  Ap_(i+1) - Ap_(i) 
    !   cols   (1:nnz) = Ai_  (Ap_(i)   : Ap_(i+1)-1)-1 !
    !   values (1:nnz) = alhs_(Ap_(i)   : Ap_(i+1)-1) 
       j = 1
       eq = i - 1
       nnz=  Ap_(i+1) - Ap_(i)
       do k = Ap_(i), Ap_(i+1)-1
            cols(j)   = Ai_(Ap_(i) +  (k - Ap_(i)))   - 1
            values(j) = alhs_(k)
            j = j + 1
       end do
     !write(*,'(15e12.4)') values(1:nnz) 
!      write(*,'(15i4)') eq, cols(1:nnz) 
       call HYPRE_IJMatrixSetValues(A_, 1, nnz, eq, cols, values, ierr)

!     Note that here we are setting one row at a time, though
!     one could set all the rows together (see the User's Manual).
!        write(*,*) ' equacao ... ... ... ', i-1
!    incluindo a parte inferior da matriz : informando a matriz toda
      do j = 2, nnz
        eqB = cols(j)
        colsB(1) = eq
        !write(*,*)" ++++27",  eq, cols(j), values(j)
!        write(*,*)" ++++27",  eqB, colsB(1), values(j)
       call HYPRE_IJMatrixSetValues(A_, 1, 1, eqB, colsB(1), values(j), ierr)
      end do
    enddo

    local_size = Cupper_ - Clower_ + 1
    print*, " ++++ B atribuindo valores de RHS para o HYPRE_IJVector ", Flower, Fupper
     do i = Flower, Fupper
!       rows_(i)  = i - 1
     end do

      call HYPRE_IJVectorSetValues (b_, local_size, rows_, rhs_, ierr )
!      write(*,*) " ++++ local_size = ", local_size
!      write(*,*) " ++++ rhs_, ", rhs_(1:local_size)
      !write(*,*) " ++++ rows+, ", rows_(1:local_size)

      print*, " +++++ atribuindo valores de guess para o HYPRE_IJVector "
      call HYPRE_IJVectorSetValues (u_, local_size, rows_, x_, ierr)
#endif
 end subroutine atribuirValoresSistemaAlgHYPRE

!
!    
      subroutine solverHYPRE(alhs_, brhs_, initialGuess_,  Ap_, Ai_, neq_, nonzerosT_, myid_, num_procs_, mpi_comm_)

      implicit none 
!
      integer*4, intent(in)    :: nonzerosT_
      integer*4, intent(in)   :: neq_
      real*8 , intent(in)    :: alhs_(*)
      real*8 , intent(inout) :: brhs_(neq_), initialGuess_(neq_)
      integer*4, intent(in)  :: Ap_(*), Ai_(*)
      integer*4, intent(in)  :: myid_, num_procs_
      integer*4              :: mpi_comm_
!

      integer*4           :: solver_id, precond_id, print_solution 
      integer*8           :: A, parcsr_A, b, par_b, u, par_u, solver
      character (LEN=16)  :: SolverStatus="incomplete"
      integer*8           :: Clower, Cupper
      integer*4           :: rows(neq_)
      
      double precision, allocatable  :: x(:)


      integer*8        :: t1, t2, t0, t3, clock_rate, clock_max
      integer*4        :: num_iterations, ierr
      double precision :: final_res_norm
      double precision :: elapsedT, tol

#ifdef withHYPREEXCLUIDO
      call system_clock                   (t0, clock_rate, clock_max)
      allocate(x(neq_)); x=0.0d0
      x = initialGuess_

      write(*,*) ' +++  iniciando a execucao do  HYPRE '
      write(*,*)  myid_, num_procs_, mpi_comm_
!
      solver_id       = 2; ! PARSAIL
      solver_id       = 1; ! PCG+AMG
      print_solution  = 1
 
!-----------------------------------------------------------------------
!     Initialize MPI
!-----------------------------------------------------------------------
!     call inicializarMPI                 (myid, num_procs, mpi_comm)

       Clower = 0 
       Cupper = neq_-1

      call criarSistemaAlgHYPRE(A, parcsr_A, b, par_b, u, par_u, solver, Clower, Cupper, mpi_comm_)

      call atribuirValoresSistemaAlgHYPRE (alhs_ , Ap_, Ai_, brhs_, x, rows, &
                                         neq_, nonZerosT_, Clower, Cupper,  & 
                                         A, parcsr_A, b, par_b, u, par_u, solver,  &
                                         myid_, mpi_comm_) 

      call finalizarMontagemSistemaHYPRE  (A, parcsr_A, b, par_b, u, par_u)

      print* , 'solver_id =', solver_id
      call system_clock                   (t1, clock_rate, clock_max)

      call resolverSistemaAlgHYPRE    (A, parcsr_A, b, par_b, u, par_u, &
                                       solver, solver_id, precond_id, tol,  &
                                     num_iterations, final_res_norm, myid_, mpi_comm_, x, rows, neq_)   
      print*, "----------------depois de resolverSistemaAlgHYPRE"

      call system_clock                   (t2, clock_rate, clock_max)
      elapsedT  =  real ( t2 - t1 ) / real ( clock_rate )


      call escreverResultadosHYPRE        (u, num_iterations, final_res_norm, elapsedT,&
                                            myId_, print_solution)

!     Clean up
      call destruirSistemaAlgHYPRE        (A, b, u)

      brhs_ = x

      deallocate(x)
      call system_clock                   (t3, clock_rate, clock_max)
      write(*,*) " +++ em solverHYPRE, tempo = ", real ( t3 - t0 ) / real ( clock_rate )
#endif
      end subroutine solverHYPRE
     end module mHYPREsistema

