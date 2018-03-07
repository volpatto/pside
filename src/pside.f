C$Id: pside.f,v 1.23 1998/11/25 09:30:17 walter Exp $
      SUBROUTINE PSIDE(NEQN,Y,DY,GEVAL,
     +                 JNUM,NLJ,NUJ,JEVAL,
     +                 MNUM,NLM,NUM,MEVAL,
     +                 T,TEND,RTOL,ATOL,IND,
     +                 LRWORK,RWORK,LIWORK,IWORK,RPAR,IPAR,
     +                 IDID)
CF90  IMPLICIT NONE
      INTEGER NEQN,NLJ,NUJ,NLM,NUM,IND(*),LRWORK,LIWORK,
     +                 IWORK(LIWORK),IPAR(*),IDID
      DOUBLE PRECISION Y(NEQN),DY(NEQN),T,TEND,RTOL(*),ATOL(*),
     +                 RWORK(LRWORK),RPAR(*)
      LOGICAL JNUM,MNUM
      EXTERNAL GEVAL,JEVAL,MEVAL
CF90  INTENT(IN)    NEQN,JNUM,NLJ,NUJ,MNUM,NLM,NUM,TEND,RTOL,ATOL,IND,
CF90 +              LRWORK,LIWORK
CF90  INTENT(INOUT) Y,DY,T,RWORK,IWORK,RPAR,IPAR
CF90  INTENT(OUT)   IDID
C-----------------------------------------------------------------------
C
C Purpose
C =======
C
C PSIDE - Parallel Software for Implicit Differential Equations - is a
C Fortran 77 code for solving implicit differential equations on shared
C memory parallel computers.
C
C Authors
C =======
C
C Jacques J.B. de Swart, Walter M. Lioen, and Wolter A. van der Veen
C CWI, P.O. Box 94079, 1090 GB Amsterdam, The Netherlands
C
C Availability
C ============
C
C The latest version, including documentation [3,4] and example drivers
C can be found at the PSIDE home page [1].  The PSIDE Users' Guide
C contains information on installing PSIDE [3, Section 7].
C
C This is version 1.3, November 25, 1998.
C
C Description
C ===========
C
C PSIDE solves Implicit Differential Equations (IDEs) of the form
C
C                                            d
C          g(t,y,y') = 0,    g,y element of R ,
C                                                                    (1)
C          t  <= t <= t   ,  y(t ) = y ,         y'(t ) = y' ,
C           0          end      0     0              0      0
C
C where y  and y'  are such that g(t ,y ,y' ) = 0 (for higher-index
C        0       0                  0  0   0
C problems the initial values have to satisfy more conditions; see [3,
C Section 4]).  It uses the four-stage Radau IIA method.  The nonlinear
C systems are solved by a modified Newton process, in which every Newton
C iterate itself is computed by means of the Parallel Iterative Linear
C system Solver for Runge-Kutta (PILSRK) proposed in [2].  This process
C is constructed such that the four stage values can be computed
C simultaneously, thereby making PSIDE suitable for execution on four
C processors.  Full details about the algorithmic choices and the
C implementation of PSIDE can be found in [4].
C
C Arguments
C =========
C
C NEQN
C    On entry, this is the dimension d of the IDE (1), the number of
C    equations to be solved.
C
C Y(NEQN)
C    On entry, this array contains the initial value y .
C                                                     0
C    On exit, Y contains y(T), the computed solution approximation at T.
C    (After successful return, T = TEND.)
C
C DY(NEQN)
C    On entry, this array contains the initial value y' .
C                                                      0
C    On exit, DY contains y'(T), the computed derivative approximation
C    at T.
C    (After successful return, T = TEND.)
C
C GEVAL
C    This is the subroutine which you provide to define the IDE
C
C          SUBROUTINE GEVAL(NEQN,T,Y,DY,G,IERR,RPAR,IPAR)
C          INTEGER NEQN,IERR,IPAR(*)
C          DOUBLE PRECISION T,Y(NEQN),DY(NEQN),G(NEQN),RPAR(*)
C    C     INTENT(IN)    NEQN,T,Y,DY
C    C     INTENT(INOUT) IERR,RPAR,IPAR
C    C     INTENT(OUT)   G
C
C    For the given values of T, Y, and DY the subroutine should return
C    the residual of the IDE
C
C                             G = g(T,Y,DY).
C
C    You must declare the name GEVAL in an EXTERNAL statement in your
C    program that calls PSIDE.
C
C    IERR is an integer flag which is always equal to zero on input.
C    Subroutine GEVAL should set IERR = -1 if GEVAL can not be evaluated
C    for the current values of Y and DY.  PSIDE will then try to prevent
C    IERR = -1 by using a smaller stepsize.
C
C    All other parameters have the same meaning as within subroutine
C    PSIDE.
C
C JNUM
C    To solve the IDE it is necessary to use the partial derivatives
C    J = dg/dy.  The solution will be more reliable if you provide J via
C    the subroutine JEVAL, in this case set JNUM = .FALSE..  If you do
C    not provide a subroutine to evaluate J, provide a dummy JEVAL, set
C    JNUM = .TRUE. and PSIDE will approximate J by numerical
C    differencing.
C
C NLJ and NUJ
C    If J is a full matrix, set NLJ = NEQN, otherwise set NLJ and NUJ
C    equal to the lower bandwidth and upper bandwidth of J,
C    respectively.
C
C JEVAL
C    This is the subroutine which you provide to define J (if JNUM .EQ.
C    .FALSE.)
C
C          SUBROUTINE JEVAL(LDJ,NEQN,NLJ,NUJ,T,Y,DY,DGDY,RPAR,IPAR)
C          INTEGER LDJ,NEQN,NLJ,NUJ,IPAR(*)
C          DOUBLE PRECISION T,Y(NEQN),DY(NEQN),DGDY(LDJ,NEQN),RPAR(*)
C    C     INTENT(IN)    LDJ,NEQN,NLJ,NUJ,T,Y,DY
C    C     INTENT(INOUT) RPAR,IPAR
C    C     INTENT(OUT)   DGDY
C
C    For the given values of T, Y, and DY the subroutine should return
C    the partial derivatives, such that
C
C     - DGDY(I,J) contains dg (T,Y,DY)/dy  if J is a full matrix
C                            I           J
C       (NLJ = NEQN);
C
C     - DGDY(I-J+NUJ+1,J) contains dg (T,Y,DY)/dy  if J is a band matrix
C                                    I           J
C       (0 <= NLJ < NEQN) (LAPACK / LINPACK / BLAS storage).
C
C    You must declare the name JEVAL in an EXTERNAL statement in your
C    program that calls PSIDE.
C
C    LDJ denotes the leading dimension of J.
C
C    All other parameters have the same meaning as within subroutine
C    PSIDE.
C
C MNUM
C    To solve the IDE it is necessary to use the partial derivatives
C    M = dg/dy'.  The solution will be more reliable if you provide M
C    via MEVAL, in this case set MNUM = .FALSE..  If you do not provide
C    a subroutine to evaluate M, provide a dummy MEVAL, set MNUM =
C    .TRUE. and PSIDE will approximate M by numerical differencing.
C
C NLM and NUM
C    If M is a full matrix, set NLM = NEQN, otherwise set NLM and NUM
C    equal to the lower bandwidth and upper bandwidth of M,
C    respectively.  It is supposed that NLM .LE. NLJ and NUM .LE. NUJ.
C
C MEVAL
C    This is the subroutine which you provide to define M (if MNUM .EQ.
C    .FALSE.)
C
C          SUBROUTINE MEVAL(LDM,NEQN,NLM,NUM,T,Y,DY,DGDDY,RPAR,IPAR)
C          INTEGER LDM,NEQN,NLM,NUM,IPAR(*)
C          DOUBLE PRECISION T,Y(NEQN),DY(NEQN),DGDDY(LDM,NEQN),RPAR(*)
C    C     INTENT(IN)    LDM,NEQN,NLM,NUM,T,Y,DY
C    C     INTENT(INOUT) RPAR,IPAR
C    C     INTENT(OUT)   DGDDY
C
C    For the given values of T, Y, and DY the subroutine should return
C    the partial derivatives, such that
C
C     - DGDDY(I,J) contains dg (T,Y,DY)/dy'  if M is a full matrix
C                             I            J
C       (NLM = NEQN);
C
C     - DGDDY(I-J+NUM+1,J) contains dg (T,Y,DY)/dy'  if M is a band
C                                     I            J
C       matrix (0 <= NLM < NEQN) (LAPACK / LINPACK / BLAS storage).
C
C    You must declare the name MEVAL in an EXTERNAL statement in your
C    program that calls PSIDE.
C
C    LDM denotes the leading dimension of M.
C
C    All other parameters have the same meaning as within subroutine
C    PSIDE.
C
C T
C    On entry, T must specify t , the initial value of the independent
C                              0
C    variable.
C    On successful exit (IDID .EQ. 1), T contains TEND.
C    On an error return, T is the point reached.
C
C TEND
C    On entry, TEND must specify the value of the independent variable
C    at which the solution is desired.
C
C RTOL and ATOL
C    You must assign relative RTOL and absolute ATOL error tolerances to
C    tell the code how small you want the local errors to be.  You have
C    two choices
C
C     - both RTOL and ATOL are scalars (set IWORK(1) = 0): the code
C       keeps, roughly, the local error of Y(I) below
C       RTOL*ABS(Y(I))+ATOL;
C
C     - both RTOL and ATOL are vectors (set IWORK(1) = 1): the code
C       keeps the local error of Y(I) below RTOL(I)*ABS(Y(I))+ATOL(I).
C
C    In either case all components must be non-negative.
C
C IND
C    If IWORK(2) .EQ. 1 , then IND should be declared of length NEQN and
C    IND(I) must specify the index of variable I.  If IWORK(2) .EQ. 0 ,
C    then IND is not referenced and the problem is assumed to be of
C    index 1.
C    See [3, Section 4] for information how to determine the index of
C    variables of certain problem classes.
C
C LRWORK
C    On entry LRWORK must specify the length of the RWORK array.  You
C    must have for the full partial derivatives case (when NLJ = NEQN)
C
C       LRWORK .GE. 20 + 27*NEQN + 6*NEQN**2,
C
C    for the case where M is banded and J is full (when NLJ = NEQN and
C    NLM < NEQN)
C
C       LRWORK .GE. 20 + (27 + NLM+NUM+1 + 5*NEQN)*NEQN,
C
C    and for the case where both partial derivatives are banded (when
C    NLJ < NEQN)
C
C       LRWORK .GE. 20 + (27 + NLJ+NUJ+NLM+NUM+2 +
C                              4*(2*NLJ+NUJ+1))*NEQN.
C
C RWORK
C    Real work array of length LRWORK.  RWORK(1),...,RWORK(20) serve as
C    parameters for the code.  For standard use, set
C    RWORK(1),...,RWORK(20) to zero before calling.
C
C    On entry:
C
C     - if RWORK(1) .GT. 0D0 then PSIDE will use RWORK(1) as initial
C       stepsize instead of determining it internally.
C
C    On exit:
C
C     - RWORK(1) contains the stepsize used on the last successful step.
C
C LIWORK
C    On entry LIWORK must specify the length of the IWORK array.  You
C    must have
C
C       LIWORK .GE. 20 + 4*NEQN.
C
C IWORK
C    Integer work array of length LIWORK.  IWORK(1),...,IWORK(20) serve
C    as parameters for the code.  For standard use, set
C    IWORK(1),...,IWORK(20) to zero before calling.
C
C    On entry:
C
C     - if IWORK(1) .EQ. 1 then RTOL and ATOL are vectors instead of
C       scalars,
C
C     - if IWORK(2) .EQ. 1 then IND is a vector,
C
C     - set IWORK(10) = 0 if PSIDE is called for the first time; for
C       subsequent calls of PSIDE do not reinitialize the parameters
C       IWORK(10),...,IWORK(19) to zero.
C
C    On exit:
C
C     - IWORK(10) contains the number of successive PSIDE calls,
C
C     - IWORK(11) contains the number of g evaluations,
C
C     - IWORK(12) contains the number of J and M evaluations (J and M
C       are computed in tandem and count as 1),
C
C     - IWORK(13) contains the number of LU-decompositions.
C
C     - IWORK(14) contains the number of forward/backward solves,
C
C     - IWORK(15) contains the total number of steps (including rejected
C       steps),
C
C     - IWORK(16) contains the number of rejected steps due to error
C       control,
C
C     - IWORK(17) contains the number of rejected steps due to Newton
C       failure,
C
C     - IWORK(18) contains the number of rejected steps due to excessive
C       growth of the solution,
C
C     - IWORK(19) contains the number of rejected steps due to IERR .EQ.
C       -1 return of GEVAL.  The integration characteristics in
C       IWORK(11),...,IWORK(14) refer to an implementation on a one-
C       processor computer.  When implemented on a parallel computer
C       with four processors, one may divide these numbers by four to
C       obtain the number of sequential evaluations, decompositions and
C       solves.
C
C RPAR and IPAR
C    RPAR and IPAR are double precision and integer arrays which you can
C    use for communication between your calling program and the
C    subroutines GEVAL, and/or JEVAL, MEVAL.  They are not altered by
C    PSIDE.  If you do not need RPAR and IPAR, ignore these parameters
C    by treating them as dummy arguments.  If you choose to use them,
C    dimension them in GEVAL and/or JEVAL, MEVAL as arrays of
C    appropriate length.  Because of the parallel implementation of
C    PSIDE, GEVAL must not alter RPAR and IPAR to prevent concurrent
C    updating.  JEVAL and MEVAL may alter them.
C
C IDID
C    On exit:
C
C     - if IDID .EQ. 1 then the integration was successful,
C
C     - if IDID .EQ. -1 then PSIDE could not reach TEND because the
C       stepsize became too small,
C
C     - if IDID .EQ. -2 then something else went wrong.
C       For example this happens when the input was invalid.  An error
C       message will be printed.
C
C Acknowledgements
C ================
C
C This project was supported by STW (project number CWI22.2703).
C STW is the Dutch Technology Foundation.
C
C This work was sponsored by the Stichting Nationale
C Computerfaciliteiten (National Computing Facilities Foundation, NCF)
C for the use of supercomputer facilities, with financial support from
C the Nederlandse Organisatie voor Wetenschappelijk Onderzoek
C (Netherlands Organization for Scientific Research, NWO).
C
C References
C ==========
C
C 1. PSIDE Home Page,  http://www.cwi.nl/cwi/projects/PSIDE/
C
C 2. P.J. van der Houwen and J.J.B. de Swart, "Parallel Linear System
C    Solvers for Runge-Kutta Methods," Advances in Computational
C    Mathematics Vol. 7 pp. 157-181 (1997).
C
C 3. W.M. Lioen, J.J.B. de Swart, and W.A. van der Veen, PSIDE Users'
C    Guide, CWI, Amsterdam (1997).  Available at [1].
C
C 4. J.J.B. de Swart, W.M. Lioen, and W.A. van der Veen, Specification
C    of PSIDE, CWI, Amsterdam (1997).  Available at [1].
C
C-----------------------------------------------------------------------
C begin included $Id: disclm.h,v 1.1 1998/11/25 08:58:23 walter Exp $
C
C THE TEXT ON THIS PAGE IS ADDRESSED TO ALL PSIDE USERS.
C
C
C BY USING THE PSIDE SOFTWARE, YOU ARE CONSENTING TO BE BOUND BY THIS
C AGREEMENT.
C
C LICENSE
C
C SMC and STW grant you a non-exclusive license to use the PSIDE
C Software free of charge for research purposes.
C
C DISCLAIMER OF WARRANTY
C
C The PSIDE Software is provided on an "AS IS" basis, without warranty
C of any kind, including without limitation the warranties of
C merchantability, fitness for a particular purpose and
C non-infringement. The entire risk as to the quality and performance
C of the Software is borne by you. You must determine that the Software
C sufficiently meets your requirements.
C
C LIMITATION OF LIABILITY
C
C UNDER NO CIRCUMSTANCES AND UNDER NO LEGAL THEORY, TORT, CONTRACT, OR
C OTHERWISE, SHALL SMC OR STW BE LIABLE TO YOU OR ANY OTHER PERSON FOR
C ANY INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES OF ANY
C CHARACTER INCLUDING, WITHOUT LIMITATION, DAMAGES FOR LOSS OF
C GOODWILL, WORK STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, OR ANY AND
C ALL OTHER COMMERCIAL DAMAGES OR LOSSES. IN NO EVENT WILL SMC BE
C LIABLE FOR ANY DAMAGES, EVEN IF SMC SHALL HAVE BEEN INFORMED OF THE
C POSSIBILITY OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. YOU
C AGREE TO INDEMNIFY AND HOLD SMC AND STW HARMLESS WITH RESPECT TO ALL
C CLAIMS BY THIRD PARTIES ARISING OUT OF YOUR USE OF THE RESULTS OR
C OPERATION OF THE SOFTWARE.
C
C TITLE
C
C Title, ownership rights, and intellectual property rights in the
C Software shall remain with SMC and STW. The Software is protected by
C the copyright laws and treaties.
C
C MISCELLANEOUS
C
C You may not remove any proprietary notices or labels on the Software,
C nor remove this disclaimer. If any provision of this Agreement is
C held to be unenforceable, such provision shall be reformed only to
C the extent necessary to make it enforceable.
C
C This Agreement shall be governed by the Law of the Netherlands.
C
C SMC stands for Stichting Mathematisch Centrum and is the Dutch
C Foundation for promotion of mathematics and computer science and their
C applications. SMC administers CWI, which is the National Research
C Institute for Mathematics and Computer Science.
C
C STW stands for Stichting Technische Wetenschappen and is the Dutch
C Foundation for technical sciences. STW stimulates
C technical-scientific research and its applications.
C
C end included $Id: disclm.h,v 1.1 1998/11/25 08:58:23 walter Exp $
C-----------------------------------------------------------------------
C
C meaning of the identifiers used in PSIDE routines:
C
C-----------------------------------------------------------------------
C
C     A         : Runge-Kutta matrix (method parameter)
C     ALPHA     : estimated rate of convergence
C     ALPHA1    : strategy parameter in VERGEN:
C                 initial rate of convergence
C     ALPJAC    : strategy parameter in CTRL:
C                 threshold for update of Jacobians
C     ALPLU     : strategy parameter in CTRL:
C                 threshold for update of factorization of iteration
C                 matrix
C     ALPREF    : strategy parameter in CTRL and VERGEN:
C                 desired rate of convergence
C     ATOL      : absolute error tolerance
C     AUX       : auxiliary array of dimension NEQN
C     AUX2      : auxiliary array of dimension NEQN
C     AUX3      : auxiliary array of dimension NEQN
C     AUXS      : auxiliary array of dimension NEQN x S
C     AUXS2     : auxiliary array of dimension NEQN x S
C     AUXS3     : auxiliary array of dimension NEQN x S
C     B         : block diagonal matrix (method parameter)
C     B0        : free parameter in reference formula (method parameter)
C     C         : abscissa vector (method parameter)
C     COMPH0    : function that computes initial stepsize
C     D         : diagonal of diagonal matrix D (method parameter)
C     DECLUS    : array for containing S LU decompositions
C     DEFAUL    : strategy parameter in COMPH0:
C                 default value for H0
C     DELTA     : saves a difference (only in NJEVAL and NMEVAL)
C     DGDDY     : derivative of G w.r.t. DY
C     DGDY      : derivative of G w.r.t. Y
C     DIVER     : logical which tells whether Newton diverges
C     DSINV     : inverse of D(S) (only in ERROR)
C     DY        : steppoint derivative value
C     DYS       : stage derivative vector
C     DYSP      : previous stage derivative vector
C     EPS       : estimated error (only in ERROR)
C     EPSP      : previous estimated error (only in ERROR)
C     EPSREJ    : estimated error of the previous rejected step
C                 (only in ERROR)
C     EXACT     : logical which tells whether Newton has found the exact
C                 solution
C     FACNEW    : logical which is true if iteration matrix must be
C                 factorized
C     FH0       : strategy parameter in COMPH0:
C                 determines influence of integration interval on H0
C     FIRST     : logical which is only true in first step
C     FMAX      : strategy parameter in CTRL:
C                 maximal factor by which stepsize can change
C     FMIN      : strategy parameter in CTRL:
C                 minimal factor by which stepsize can change
C     FRIG      : strategy parameter in CTRL:
C                 rigid factor by which H is divided if analysis fails
C     GAM       : strategy parameter in NJEVAL and NMEVAL:
C                 if component is smaller than GAM, it is set equal to
C                 GAM
C     GAMMA     : strategy parameter in VERGEN:
C                 determines when process is considered to be diverging
C     GEVAL     : routine containing the function g
C     GFAC      : strategy parameter in VERGEN:
C                 the maximal allowed growfactor
C     GROWTH    : logical which tells whether iterates are growing
C                 excessively
C     H         : stepsize
C     HALPHA    : stepsize that should lead to desired rate of
C                 convergence (only in CTRL)
C     HLU       : stepsize for iteration matrix
C     HNEW      : new proposed stepsize (only in CTRL)
C     HP        : previous stepsize
C     HR        : stepsize suggested by error controller
C     HREJ      : stepsize of the previous rejected step (only in ERROR)
C     IDID      : indicator reporting what the code did
C     IERR      : scalar integer for communication between subroutine
C                 GEVAL and PSIDE
C     IERRS     : integer array of dimension S for communication
C                 between subroutine GEVAL and PSIDE (only in PILSRK)
C     IIPVTS    : pointer in IWORK
C     IND       : vector with index for components of Y
C     INDGT1    : switch for higher index problems
C     INFO      : linear algebra return code (only in JACFAC)
C     IPAR      : integer array for communication between
C                 subroutines GEVAL/JEVAL and PSIDE
C     IPVTS     : pivots for S LU decompositions
C     IWORK     : work array
C     JACNEW    : logical which is true if Jacobians must be updated
C     JACU2D    : logical which tells whether the Jacobians are
C                 up to date
C     JBND      : switch for banded form of DGDY
C     JEVAL     : routine containing DGDY
C     JNUM      : boolean switch for numerical or analytical DGDY
C     K         : index for Newton iteration
C     KAPPA     : strategy parameter in VERGEN:
C                 if relative update is smaller than KAPPA*UROUND
C                 we stop
C     KB        : index for loops over blocks in method parameter matrix
C                 (only in PILSRK)
C     KM        : index for loops over PILSRK iteration (only in PILSRK)
C     KMAX      : strategy parameter in VERGEN:
C                 maximal allowed number of iterations
C     KN,KN2    : index for loops over the problem dimension
C     KS,KS2,KS3: index for loops over the stages
C     LDJ       : leading dimension of DGDY
C     LDLU      : leading dimension of LU decompositions
C     LDM       : leading dimension of DGDDY
C     LIWORK    : dimension of IWORK
C     LRWORK    : dimension of RWORK
C     M         : number of PILSRK iterations (only in PILSRK)
C     MBND      : switch for banded form of DGDDY
C     MEVAL     : routine containing DGDDY
C     MNUM      : boolean switch for numerical or analytical DGDDY
C     NEQN      : dimension of the problem
C     NF        : number of function calls
C     NFB       : number of forward backward substitutions of problem
C                 dimension
C     NJAC      : number of Jacobian evaluation (1 Jacobian evaluation
C                 counts for the computation of both DGDY and DGDDY)
C     NLJ       : number of non-zero lower co-diagonals in DGDY
C     NLM       : number of non-zero lower co-diagonals in DGDDY
C     NLU       : number of LU decompositions of matrix of dimension
C                 NEQN x NEQN
C     NREJE     : number of rejected steps due to error test failure
C     NREJG     : number of rejected steps due to excessive growth of
C                 the solution
C     NREJI     : number of rejected steps due to IERR .EQ. -1 return
C                 of GEVAL
C     NREJN     : number of rejected steps due to convergence
C                 failure of the Newton process
C     NST       : total number of steps
C     NUJ       : number of non-zero upper co-diagonals in DGDY
C     NUM       : number of non-zero upper co-diagonals in DGDDY
C     OMEGA     : parameter to adjust stepsize such that the remaining
C                 integration interval is a multiple of the stepsize
C                 (only in CTRL)
C     PEST      : estimated order of the reference formula
C                 (only in ERROR)
C     PMIN      : strategy parameter in ERROR:
C                 minimal effective order of reference formula
C     Q         : transformation matrix (method parameter)
C     QINV      : inverse of Q (method parameter)
C     R         : error estimate of dimension NEQN (only in ERROR)
C     RAUX      : pointer in RWORK
C     RAUX2     : pointer in RWORK
C     RAUX3     : pointer in RWORK
C     RAUXS     : pointer in RWORK
C     RAUXS2    : pointer in RWORK
C     RAUXS3    : pointer in RWORK
C     RDECLU    : pointer in RWORK
C     RDGDDY    : pointer in RWORK
C     RDGDY     : pointer in RWORK
C     RDYS      : pointer in RWORK
C     RDYSP     : pointer in RWORK
C     READY     : logical which is true if Newton iteration
C                 is finished
C     REMAIN    : the estimated number of remaining steps (only in CTRL)
C     RH        : ratio of subsequent stepsizes (only in PREDIC)
C     RPAR      : double precision array for communication between
C                 subroutines GEVAL/JEVAL and PSIDE
C     RTOL      : relative error tolerance
C     RWORK     : work array
C     RYS       : pointer in RWORK
C     S         : the number of stages (method parameter)
C     SLOW      : logical which tells whether Newton process is
C                 converging too slowly
C     SNORM     : function that computes the scaled norm of a
C                 vector of dimension NEQN
C     SNORMS    : function that computes scaled norm of a vector of
C                 dimension NEQN x S
C     SNRM      : scalar to store the scaled norm of DY (only in COMPH0)
C     SOLVED    : logical which tells whether Newton has converged
C     SQUR      : square root of UROUND (only in NJEVAL and NMEVAL)
C     SUCREJ    : logical which is true if steps are rejected
C                 successively (only in ERROR)
C     T         : the current value of the independent variable
C     TAU       : strategy parameter in VERGEN:
C                 determines how small iteration error has to be w.r.t.
C                 local error
C     TBEGIN    : start of integration interval (only in COMPH0)
C     TEND      : end of integration interval
C     TGS       : transformed of -G(T+C*H,YS,DYS) (only in PILSRK)
C     THETA     : strategy parameter in VERGEN:
C                 determines influence of previous rates of convergence
C     TOLVEC    : switch for tolerance vector or scalar
C     TUDYS     : update for transformed stage vector (only in PILSRK)
C     U         : norm of update (only in VERGEN)
C     UDYS      : update for stage vector (only in NEWTON and PILSRK)
C     UINV      : auxiliary S x S matrix (only in PREDIC)
C     UP        : norm of previous update (only in VERGEN)
C     UROUND    : unit roundoff error
C     UYS       : update for stage vector (only in NEWTON and VERGEN)
C     V         : other parameters for reference formula
C                 (method parameter)
C     VM        : auxiliary S x S matrix (only in PREDIC)
C     W         : auxiliary S x S matrix (only in PREDIC)
C     XI        : strategy parameter in CTRL:
C                 prevents that the stepsize is reduced too little if
C                 the Newton process converges too slowly
C     Y         : steppoint value
C     YNORM     : vector of dimension NEQN of which we want to compute
C                 the norm (only in SNORM)
C     YS        : stage vector
C     YSNORM    : vector of dimension NEQN x S of which we want to
C                 compute the norm (only in SNORMS)
C     ZETA      : strategy parameter in COMPH0:
C                 determines influence of DY on H0
C     ZETA      : strategy parameter in ERROR:
C                 safety factor for stepsize selection
C-----------------------------------------------------------------------
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER LDJ,LDM,LDLU,
     +        NST,NREJE,NREJN,NREJG,NREJI,NF,NFB,NJAC,NLU,
     +        RYS,RDYS,RDYSP,RDECLU,RDGDY,RDGDDY,
     +        RAUX,RAUX2,RAUX3,RAUXS,RAUXS2,RAUXS3,
     +        IIPVTS
      DOUBLE PRECISION H,UROUND
      LOGICAL JBND,MBND,TOLVEC,INDGT1
      DOUBLE PRECISION DLAMCH

      UROUND = DLAMCH('Epsilon')

      JBND = NLJ.NE.NEQN
      MBND = NLM.NE.NEQN

      IF (NLM.GT.NLJ) THEN
         IDID = -2
         PRINT *, 'PSIDE: ERROR: on input NLM > NLJ ', NLM, NLJ
         RETURN
      ENDIF
      IF (JBND .AND. (.NOT.MBND .OR. MBND.AND.NUM.GT.NUJ)) THEN
         IDID = -2
         PRINT *, 'PSIDE: ERROR: on input NUM > NUJ', NUM, NUJ
         RETURN
      ENDIF

      IF (JBND) THEN
         LDJ  = NLJ+NUJ+1
         LDLU = 2*NLJ+NUJ+1
      ELSE
         LDJ  = NEQN
         LDLU = NEQN
      ENDIF
      IF (MBND) THEN
         LDM  = NLM+NUM+1
      ELSE
         LDM  = NEQN
      ENDIF

      IF (T.GE.TEND) THEN
         IDID = -2
         PRINT *, 'PSIDE: ERROR: on input T >= TEND', T, TEND
         RETURN
      ENDIF

      TOLVEC = IWORK(1).EQ.1
      INDGT1 = IWORK(2).EQ.1

      H=RWORK(1)
C
C     partition RWORK into the following arrays
C
C     double precision YS(NEQN,S)
C     double precision DYS(NEQN,S)
C     double precision DYSP(NEQN,S)
C     double precision DECLUS(LDLU,NEQN,S)
C     double precision DGDY(LDJ,NEQN)
C     double precision DGDDY(LDM,NEQN)
C     double precision AUX(NEQN)
C     double precision AUX2(NEQN)
C     double precision AUX3(NEQN)
C     double precision AUXS(NEQN,S)
C     double precision AUXS2(NEQN,S)
C     double precision AUXS3(NEQN,S)
C
C     format: pointer = previous pointer + previous length
C
      RYS    = 21
      RDYS   = RYS    + NEQN*S
      RDYSP  = RDYS   + NEQN*S
      RDECLU = RDYSP  + NEQN*S
      RDGDY  = RDECLU + LDLU*NEQN*S
      RDGDDY = RDGDY  + LDJ*NEQN
      RAUX   = RDGDDY + LDM*NEQN
      RAUX2  = RAUX   + NEQN
      RAUX3  = RAUX2  + NEQN
      RAUXS  = RAUX3  + NEQN
      RAUXS2 = RAUXS  + NEQN*S
      RAUXS3 = RAUXS2 + NEQN*S
C                     + NEQN*S

      IF (RAUXS3+NEQN*S-1 .GT. LRWORK) THEN
         IDID = -2
         PRINT *, 'PSIDE: ERROR: ',
     +            'insufficient storage for RWORK ',
     +            'minimal LRWORK = ', RAUXS3+NEQN*S-1
         RETURN
      ENDIF
C
C     partition IWORK into the following array(s)
C
C     integer IPVTS(NEQN,S)
C
      IIPVTS = 21
C                    + NEQN*S

      IF (IIPVTS+NEQN*S-1 .GT. LIWORK) THEN
         IDID = -2
         PRINT *,  'PSIDE: ERROR: ',
     +             'insufficient storage for IWORK ',
     +             'minimal LIWORK = ', IIPVTS+NEQN*S-1
         RETURN
      ENDIF

      IF (IWORK(10).EQ.0) THEN
         NF    = 0
         NJAC  = 0
         NLU   = 0
         NFB   = 0
         NST   = 0
         NREJE = 0
         NREJN = 0
         NREJG = 0
         NREJI = 0
      ELSE
         NF    = IWORK(11)
         NJAC  = IWORK(12)
         NLU   = IWORK(13)
         NFB   = IWORK(14)
         NST   = IWORK(15)
         NREJE = IWORK(16)
         NREJN = IWORK(17)
         NREJG = IWORK(18)
         NREJI = IWORK(19)
      ENDIF

      CALL TSTEPS(LDJ,LDM,LDLU,NEQN,Y,DY,T,TEND,H,GEVAL,JEVAL,
     +            MEVAL,
     +            JNUM,JBND,NLJ,NUJ,MNUM,MBND,NLM,NUM,
     +            TOLVEC,RTOL,ATOL,INDGT1,IND,
     +            RPAR,IPAR,
     +            UROUND,NST,NREJN,NREJE,NREJG,NREJI,NF,NFB,NJAC,NLU,
     +            RWORK(RYS),RWORK(RDYS),RWORK(RDYSP),RWORK(RDECLU),
     +            IWORK(IIPVTS),RWORK(RDGDY),RWORK(RDGDDY),
     +            RWORK(RAUX),RWORK(RAUX2),RWORK(RAUX3),
     +            RWORK(RAUXS),RWORK(RAUXS2),RWORK(RAUXS3),IDID)

      IWORK(10) = IWORK(10)+1
      IWORK(11) = NF
      IWORK(12) = NJAC
      IWORK(13) = NLU
      IWORK(14) = NFB
      IWORK(15) = NST
      IWORK(16) = NREJE
      IWORK(17) = NREJN
      IWORK(18) = NREJG
      IWORK(19) = NREJI

      RWORK(1)  = H

      RETURN
      END
C$Id: tsteps.f,v 1.14 1998/11/25 08:58:34 walter Exp $
      SUBROUTINE TSTEPS(LDJ,LDM,LDLU,NEQN,Y,DY,T,TEND,H,
     +                  GEVAL,JEVAL,MEVAL,
     +                  JNUM,JBND,NLJ,NUJ,MNUM,MBND,NLM,NUM,
     +                  TOLVEC,RTOL,ATOL,INDGT1,IND,
     +                  RPAR,IPAR,
     +                  UROUND,NST,NREJN,NREJE,NREJG,NREJI,NF,NFB,NJAC,
     +                  NLU,
     +                  YS,DYS,DYSP,DECLUS,IPVTS,DGDY,DGDDY,AUX,AUX2,
     +                  AUX3,
     +                  AUXS,AUXS2,AUXS3,IDID)
CF90  IMPLICIT NONE
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER LDJ,LDM,LDLU,NEQN,NLJ,NUJ,NLM,NUM,IND(*),
     +        IPAR(*),IPVTS(NEQN,S),NST,NREJN,NREJE,NREJG,NREJI,NF,NFB,
     +        NJAC,NLU,IDID
      DOUBLE PRECISION Y(NEQN),DY(NEQN),T,TEND,H,
     +                 RTOL(*),ATOL(*),RPAR(*),
     +                 UROUND,YS(NEQN,S),DYS(NEQN,S),DYSP(NEQN,S),
     +                 DECLUS(LDLU,NEQN,S),
     +                 DGDY(LDJ,NEQN),DGDDY(LDJ,NEQN),
     +                 AUX(NEQN),AUX2(NEQN),AUX3(NEQN),
     +                 AUXS(NEQN,S),AUXS2(NEQN,S),
     +                 AUXS3(NEQN,S)
      LOGICAL JNUM,MNUM,JBND,MBND,TOLVEC,INDGT1
      EXTERNAL GEVAL,JEVAL,MEVAL
CF90  INTENT(IN)    LDJ,LDM,LDLU,NEQN,TEND,JNUM,JBND,NLJ,
CF90 +              NUJ,MNUM,MBND,NLM,NUM,TOLVEC,RTOL,ATOL,INDGT1,IND,
CF90 +              UROUND
CF90  INTENT(INOUT) H,T,Y,DY,RPAR,IPAR
CF90  INTENT(OUT)   NST,NREJN,NREJE,NREJG,NREJI,NF,NFB,NJAC,NLU,
CF90 +              YS,DYS,DYSP,DECLUS,IPVTS,DGDY,DGDDY,AUX,AUX2,
CF90 +              AUX3,AUXS,
CF90 +              AUXS2,AUXS3,IDID
C-----------------------------------------------------------------------
C     this is the core integrator called by PSIDE
C-----------------------------------------------------------------------
      INTEGER KN,KS,IERR
      DOUBLE PRECISION HP,HLU,ALPHA,COMPH0
      LOGICAL FIRST,GROWTH,DIVER,SLOW,SOLVED,EXACT,JACNEW,FACNEW,JACU2D
      IF (H.LE.0D0) H=COMPH0(NEQN,DY,T,TEND,TOLVEC,RTOL,ATOL,INDGT1,IND)
      HP=H
      HLU=H
      FIRST=.TRUE.
      JACNEW=.TRUE.
      FACNEW=.TRUE.
      DO 20 KS=1,S
         DO 10 KN=1,NEQN
            DYSP(KN,KS)=DY(KN)
   10    CONTINUE
   20 CONTINUE
      IDID=1
   30 IF (ABS(TEND-T).GT.10D0*UROUND*ABS(T)) THEN
         CALL JACFAC(LDJ,LDM,LDLU,NEQN,Y,DY,T,GEVAL,JEVAL,MEVAL,
     +               JNUM,MNUM,JBND,MBND,NLJ,NUJ,NLM,NUM,HLU,
     +               RPAR,IPAR,NJAC,NLU,DECLUS,IPVTS,DGDY,DGDDY,
     +               UROUND,AUX,AUX2,AUX3,JACNEW,FACNEW,JACU2D,IDID)
         IF (IDID.NE.1) GOTO 40
         CALL NEWTON(LDM,LDLU,NEQN,Y,YS,DYS,DYSP,T,HP,H,
     +               GEVAL,JBND,MBND,NLJ,NUJ,NLM,NUM,IERR,RPAR,IPAR,
     +               TOLVEC,RTOL,ATOL,INDGT1,IND,
     +               UROUND,NF,NFB,DECLUS,IPVTS,DGDDY,AUXS,
     +               AUXS2,AUXS3,ALPHA,GROWTH,DIVER,SLOW,SOLVED,EXACT)
         CALL CTRL(LDLU,NEQN,Y,DY,T,TEND,GEVAL,
     +             JBND,NLJ,NUJ,HP,H,HLU,TOLVEC,RTOL,ATOL,
     +             INDGT1,IND,UROUND,
     +             IERR,RPAR,IPAR,
     +             NREJN,NREJE,NREJG,NREJI,NF,NFB,
     +             YS,DYS,DYSP,DECLUS,IPVTS,AUX,AUX2,
     +             ALPHA,GROWTH,DIVER,SLOW,SOLVED,EXACT,
     +             FIRST,JACNEW,FACNEW,JACU2D,IDID)
         IF (IDID.NE.1) GOTO 40
         NST = NST + 1
         GOTO 30
      ENDIF
   40 RETURN
      END
C$Id: comph0.f,v 1.6 1998/11/25 08:58:21 walter Exp $
      DOUBLE PRECISION FUNCTION COMPH0(NEQN,DY,TBEGIN,TEND,
     +                          TOLVEC,RTOL,ATOL,INDGT1,IND)
CF90  IMPLICIT NONE
      INTEGER NEQN,IND(*)
      DOUBLE PRECISION DY(NEQN),TBEGIN,TEND,RTOL(*),ATOL(*)
      LOGICAL TOLVEC,INDGT1
CF90  INTENT(IN)    NEQN,DY,TBEGIN,TEND,TOLVEC,RTOL,ATOL,INDGT1,IND
C-----------------------------------------------------------------------
C     computes the initial stepsize
C-----------------------------------------------------------------------
      DOUBLE PRECISION SNRM,SNORM,FH0,ZETA,DEFAUL
      PARAMETER(ZETA   = 5D-1,
     +          FH0    = 1D-5,
     +          DEFAUL = 1D-5)
      COMPH0=MIN(DEFAUL,FH0*ABS(TEND-TBEGIN))
      SNRM=SNORM(NEQN,DY,DY,1D0,TOLVEC,RTOL,ATOL,INDGT1,IND)
      IF (SNRM.GT.ZETA/COMPH0) COMPH0=ZETA/SNRM
      COMPH0=SIGN(COMPH0,TEND-TBEGIN)
      RETURN
      END
C$Id: jacfac.f,v 1.15 1998/11/25 08:58:24 walter Exp $
      SUBROUTINE JACFAC(LDJ,LDM,LDLU,NEQN,Y,DY,T,GEVAL,JEVAL,MEVAL,
     +                  JNUM,MNUM,JBND,MBND,NLJ,NUJ,NLM,NUM,HLU,
     +                  RPAR,IPAR,NJAC,NLU,DECLUS,IPVTS,DGDY,DGDDY,
     +                  UROUND,AUX,AUX2,AUX3,JACNEW,FACNEW,JACU2D,IDID)
CF90  IMPLICIT NONE
      use omp_lib
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER LDJ,LDM,LDLU,NEQN,NLJ,NUJ,NLM,NUM,IPAR(*),NJAC,NLU,
     +        IPVTS(NEQN,S),IDID
      DOUBLE PRECISION Y(NEQN),DY(NEQN),T,HLU,RPAR(*),
     +                 DECLUS(LDLU,NEQN,S),
     +                 DGDY(LDJ,NEQN),DGDDY(LDM,NEQN),
     +                 UROUND,AUX(NEQN),AUX2(NEQN),AUX3(NEQN)
      LOGICAL JNUM,MNUM,JBND,MBND,JACNEW,FACNEW,JACU2D
      EXTERNAL GEVAL,JEVAL,MEVAL
CF90  INTENT(IN)    LDJ,LDM,LDLU,NEQN,Y,DY,T,JNUM,
CF90 +              MNUM,JBND,MBND,NLJ,NUJ,NLM,NUM,HLU,UROUND,
CF90 +              JACNEW,FACNEW
CF90  INTENT(INOUT) NJAC,NLU,IPAR,RPAR
CF90  INTENT(OUT)   DECLUS,IPVTS,DGDY,DGDDY,IDID,AUX,AUX2,AUX3,JACU2D
C-----------------------------------------------------------------------
C     evaluates the Jacobians and factorizes the iteration matrix,
C     depending on JACNEW and FACNEW
C-----------------------------------------------------------------------
C begin included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
      DOUBLE PRECISION A(S,S),B(S,S),C(S),D(S),Q(S,S),QINV(S,S),B0,V(S)
C end included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
      INTEGER KS,KN,KN2,INFO(S),IERR
C begin included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
      DATA C(1) /8.85879595127039430D-02/
      DATA C(2) /4.09466864440734710D-01/
      DATA C(3) /7.87659461760847110D-01/
      DATA C(4) /1.00000000000000000D0/

      DATA A(1,1) /1.12999479323156150D-01/
      DATA A(1,2) /-4.03092207235221240D-02/
      DATA A(1,3) /2.58023774203363160D-02/
      DATA A(1,4) /-9.90467650726640190D-03/
      DATA A(2,1) /2.34383995747400210D-01/
      DATA A(2,2) /2.06892573935359120D-01/
      DATA A(2,3) /-4.78571280485409760D-02/
      DATA A(2,4) /1.60474228065163730D-02/
      DATA A(3,1) /2.16681784623249830D-01/
      DATA A(3,2) /4.06123263867374850D-01/
      DATA A(3,3) /1.89036518170054850D-01/
      DATA A(3,4) /-2.41821048998323020D-02/
      DATA A(4,1) /2.20462211176767560D-01/
      DATA A(4,2) /3.88193468843174740D-01/
      DATA A(4,3) /3.28844319980056810D-01/
      DATA A(4,4) /6.25000000000008880D-02/

      DATA D(1) /1.52077368976571730D-01/
      DATA D(2) /1.98631665602052860D-01/
      DATA D(3) /1.73704821245586140D-01/
      DATA D(4) /2.26879766524847200D-01/

      DATA B(1,1) / -3.36398745680207070D+00/
      DATA B(1,2) / -4.46547007540097850D-01/
      DATA B(1,3) /0D0/
      DATA B(1,4) /0D0/
      DATA B(2,1) /  2.53420388412422460D+01/
      DATA B(2,2) /  3.36398745680206800D+00/
      DATA B(2,3) /0D0/
      DATA B(2,4) /0D0/
      DATA B(3,1) /0D0/
      DATA B(3,2) /0D0/
      DATA B(3,3) / -4.37367276825310510D-01/
      DATA B(3,4) / -5.80576031184040710D-02/
      DATA B(4,1) /0D0/
      DATA B(4,2) /0D0/
      DATA B(4,3) /  3.29483348541735400D+00/
      DATA B(4,4) /  4.37367276825312510D-01/

      DATA Q(1,1) /2.95256291376448490D+00/
      DATA Q(1,2) /3.15939676526544260D-01/
      DATA Q(1,3) /1.53250361857358990D+00/
      DATA Q(1,4) /2.76001773070919800D-02/
      DATA Q(2,1) /-7.26638442522609210D+00/
      DATA Q(2,2) /-8.75577120872169210D-01/
      DATA Q(2,3) /-1.05525925554083820D+00/
      DATA Q(2,4) /-3.11277680445624430D-01/
      DATA Q(3,1) /3.42060134189704890D+00/
      DATA Q(3,2) /9.49331950091266920D-01/
      DATA Q(3,3) /-1.07997190626525940D+01/
      DATA Q(3,4) /-2.13491394363750460D+00/
      DATA Q(4,1) /3.48973092842014550D+01/
      DATA Q(4,2) /4.37528029636525950D+00/
      DATA Q(4,3) /-4.29039265780423160D+01/
      DATA Q(4,4) /-5.89600020104458620D+00/

      DATA QINV(1,1) /4.94042191453623210D-01/
      DATA QINV(1,2) /2.69406327540352320D-01/
      DATA QINV(1,3) /-2.07753732935469560D-01/
      DATA QINV(1,4) /6.33161132809504090D-02/
      DATA QINV(2,1) /-3.53358212942913720D+00/
      DATA QINV(2,2) /-2.98582140519829590D+00/
      DATA QINV(2,3) /1.75646748288234010D+00/
      DATA QINV(2,4) /-4.94914306872105080D-01/
      DATA QINV(3,1) /4.87641455081039950D-01/
      DATA QINV(3,2) /1.23938205146711190D-01/
      DATA QINV(3,3) /4.23770339324015460D-02/
      DATA QINV(3,4) /-1.96050751500893220D-02/
      DATA QINV(4,1) /-3.24650638473406870D+00/
      DATA QINV(4,2) /-1.52301305545598620D+00/
      DATA QINV(4,3) /-2.34591215977400570D-01/
      DATA QINV(4,4) /-1.94525303087971780D-02/

      DATA B0 /  1.00000000000000000D-02/

      DATA V(1) /  1.57753763977411530D-02/
      DATA V(2) / -9.73676595200762000D-03/
      DATA V(3) /  6.46138955426500680D-03/
      DATA V(4) /  2.24379766524848060D-01/
C end included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
C      call omp_get_max_threads()
      !call omp_set_dynamic(.true.)
      !print*, "Aqui"
      IF (JACNEW) THEN
         IF (JNUM) THEN
            CALL NJEVAL(LDJ,NEQN,T,Y,DY,DGDY,GEVAL,UROUND,
     +                  JBND,NLJ,NUJ,IERR,RPAR,IPAR,AUX,AUX2,AUX3,IDID)
            IF (IDID.EQ.-2) RETURN
         ELSE
            CALL JEVAL(LDJ,NEQN,NLJ,NUJ,T,Y,DY,DGDY,RPAR,IPAR)
         ENDIF
         IF (MNUM) THEN
            CALL NMEVAL(LDM,NEQN,T,Y,DY,DGDDY,GEVAL,UROUND,
     +                  MBND,NLM,NUM,IERR,RPAR,IPAR,AUX,AUX2,AUX3,IDID)
            IF (IDID.EQ.-2) RETURN
         ELSE
            CALL MEVAL(LDM,NEQN,NLM,NUM,T,Y,DY,DGDDY,RPAR,IPAR)
         ENDIF
         NJAC=NJAC+1
         JACU2D=.TRUE.
      ENDIF
      IF (FACNEW) THEN
         IF (JBND) THEN
C
C           J banded -> M banded
C           NLM <= NLJ and NUM <= NUJ assumed
C
C           the extra NLJ-offset in DECLUS is required by LAPACK
C           for the storage of an additional NLJ superdiagonals,
C           generated by fill-in as a result of row interchanges
C
C      omp_get_max_threads()
C            DO 50 KS=1,S
C$OMP PARALLEL PRIVATE(KS,KN,KN2)
C$OMP&      SHARED(HLU,D,DGDY,NLM,NUM,DGDDY,
C$OMP&             NEQN,NLJ,NUJ,DECLUS,LDLU,IPVTS,INFO)
C$OMP DO
            DO 50 KS=1,S
               DO 20 KN=1,NEQN
                  DO 10 KN2=MAX(1,         NUJ+1 + 1-KN),
     +                      MIN(NLJ+NUJ+1, NUJ+1 + NEQN-KN)
                     DECLUS(NLJ+KN2,KN,KS)=HLU*D(KS)*DGDY(KN2,KN)
   10             CONTINUE
   20          CONTINUE
               DO 40 KN=1,NEQN
                  DO 30 KN2=MAX(1,         NUM+1 + 1-KN),
     +                      MIN(NLM+NUM+1, NUM+1 + NEQN-KN)
                     DECLUS(NLJ+NUJ-NUM+KN2,KN,KS)=DGDDY(KN2,KN)+
     +               DECLUS(NLJ+NUJ-NUM+KN2,KN,KS)
   30             CONTINUE
   40          CONTINUE
               INFO(KS)=0
               CALL DGBTRF(NEQN,NEQN,NLJ,NUJ,DECLUS(1,1,KS),LDLU,
     +                     IPVTS(1,KS),INFO(KS))
   50       CONTINUE
C$OMP END PARALLEL
C            ENDDO
         ELSEIF (MBND) THEN
C
C           J full and M banded
C
C            DO 100 KS=1,S
C$OMP PARALLEL PRIVATE(KS,KN,KN2)
C$OMP&      SHARED(HLU,D,DGDY,NLM,NUM,DGDDY,
C$OMP&             NEQN,DECLUS,LDLU,IPVTS,INFO)
C$OMP DO
            DO 100 KS=1,S
               DO 70 KN=1,NEQN
                  DO 60 KN2=1,NEQN
                     DECLUS(KN2,KN,KS)=HLU*D(KS)*DGDY(KN2,KN)
   60             CONTINUE
   70          CONTINUE
               DO 90 KN=1,NEQN
                  DO 80 KN2=MAX(1,         NUM+1 + 1-KN),
     +                      MIN(NLM+NUM+1, NUM+1 + NEQN-KN)
                     DECLUS(KN-1+KN2-NUM,KN,KS)=DGDDY(KN2,KN)+
     +               DECLUS(KN-1+KN2-NUM,KN,KS)
   80             CONTINUE
   90          CONTINUE
               INFO(KS)=0
               CALL DGETRF(NEQN,NEQN,DECLUS(1,1,KS),LDLU,IPVTS(1,KS),
     +                     INFO(KS))
  100       CONTINUE
C$OMP END PARALLEL
         ELSE
C
C           J and M both full
C
C            DO 130 KS=1,S
C$OMP PARALLEL PRIVATE(KS,KN,KN2)
C$OMP&      SHARED(DGDDY,HLU,D,DGDY,
C$OMP&             NEQN,DECLUS,LDLU,IPVTS,INFO)
C$OMP DO
            DO 130 KS=1,S
               DO 120 KN=1,NEQN
                  DO 110 KN2=1,NEQN
                     DECLUS(KN2,KN,KS)=DGDDY(KN2,KN)+
     +                                 HLU*D(KS)*DGDY(KN2,KN)
  110             CONTINUE
  120          CONTINUE
               INFO(KS)=0
               CALL DGETRF(NEQN,NEQN,DECLUS(1,1,KS),LDLU,IPVTS(1,KS),
     +                     INFO(KS))
  130       CONTINUE
C$OMP END PARALLEL
         ENDIF
         DO 140 KS=1,S
            IF (INFO(KS).NE.0) THEN
               IDID = -2
               PRINT *, 'PSIDE: ERROR: ',
     +                  'the iteration matrix has become singular'
               RETURN
            ENDIF
  140    CONTINUE
         NLU=NLU+S
      ENDIF
      RETURN
      END
C$Id: njeval.f,v 1.13 1998/11/25 08:58:26 walter Exp $
      SUBROUTINE NJEVAL(LDJ,NEQN,T,Y,DY,DGDY,GEVAL,UROUND,
     +                  JBND,NLJ,NUJ,IERR,RPAR,IPAR,AUX,AUX2,AUX3,IDID)
CF90  IMPLICIT NONE
      INTEGER LDJ,NEQN,NLJ,NUJ,IERR,IPAR(*),IDID
      DOUBLE PRECISION T,Y(NEQN),DY(NEQN),DGDY(LDJ,NEQN),
     +                 UROUND,RPAR(*),AUX(NEQN),AUX2(NEQN),AUX3(NEQN)
      LOGICAL JBND
CF90  INTENT(IN)    LDJ,NEQN,T,Y,DY,UROUND,JBND,NLJ,NUJ
CF90  INTENT(INOUT) IERR,RPAR,IPAR
CF90  INTENT(OUT)   DGDY,AUX,AUX2,AUX3,IDID
C-----------------------------------------------------------------------
C     computation of numerical Jacobian DGDY by finite
C     differences using 2*NEQN+1 function evaluations.
C-----------------------------------------------------------------------
      INTEGER KN,KN2
      DOUBLE PRECISION SQUR,DELTA,GAM
      PARAMETER(GAM=1D-5)
C     GAM is a threshold parameter; if a component of Y is smaller
C         than GAM, then this component is set equal to GAM;
C         for the time-being is is set equal to 1D-5; in the future,
C         GAM could be made component-dependent, like in DASSL.
      EXTERNAL GEVAL
      SQUR=SQRT(UROUND)
      IERR=0
      CALL GEVAL(NEQN,T,Y,DY,AUX,IERR,RPAR,IPAR)
      IF (IERR.NE.0) THEN
         IDID=-2
         PRINT *, 'PSIDE: ERROR: ',
     +            'cannot handle IERR = -1 in numerical Jacobians'
         RETURN
      ENDIF
      DO 5 KN=1,NEQN
         AUX2(KN)=Y(KN)
   5  CONTINUE
      DO 30 KN=1,NEQN
         DELTA=SQUR*MAX(ABS(Y(KN)),GAM)
         AUX2(KN)=Y(KN)+DELTA
         IERR=0
         CALL GEVAL(NEQN,T,AUX2,DY,AUX3,IERR,RPAR,IPAR)
         IF (IERR.NE.0) THEN
            IDID=-2
            PRINT *, 'PSIDE: ERROR: ',
     +               'cannot handle IERR = -1 in numerical Jacobians'
            RETURN
         ENDIF
         IF (JBND) THEN
            DO 10 KN2=MAX(1,KN-NUJ),MIN(NEQN,KN+NLJ)
               DGDY(KN2-KN+NUJ+1,KN)=(AUX3(KN2)-AUX(KN2))/DELTA
   10       CONTINUE
         ELSE
            DO 20 KN2=1,NEQN
               DGDY(KN2,KN)=(AUX3(KN2)-AUX(KN2))/DELTA
   20       CONTINUE
         ENDIF
         AUX2(KN)=Y(KN)
   30 CONTINUE
      RETURN
      END
C$Id: nmeval.f,v 1.10 1998/11/25 08:58:27 walter Exp $
      SUBROUTINE NMEVAL(LDM,NEQN,T,Y,DY,DGDDY,GEVAL,UROUND,
     +                  MBND,NLM,NUM,IERR,RPAR,IPAR,AUX,AUX2,AUX3,IDID)
CF90  IMPLICIT NONE
      INTEGER LDM,NEQN,NLM,NUM,IERR,IPAR(*),IDID
      DOUBLE PRECISION T,Y(NEQN),DY(NEQN),DGDDY(LDM,NEQN),
     +                 UROUND,RPAR(*),AUX(NEQN),AUX2(NEQN),AUX3(NEQN)
      LOGICAL MBND
CF90  INTENT(IN)    LDM,NEQN,T,Y,DY,UROUND,MBND,NLM,NUM
CF90  INTENT(INOUT) IERR,RPAR,IPAR
CF90  INTENT(OUT)   DGDDY,AUX,AUX2,AUX3,IDID
C-----------------------------------------------------------------------
C     computation of numerical Jacobian DGDDY by finite
C     differences using 2*NEQN+1 function evaluations.
C-----------------------------------------------------------------------
      INTEGER KN,KN2
      DOUBLE PRECISION SQUR,DELTA,GAM
      PARAMETER(GAM=1D-5)
C     GAM is a threshold parameter; if a component of Y is smaller
C         than GAM, then this component is set equal to GAM;
C         for the time-being is is set equal to 1D-5; in the future,
C         GAM could be made component-dependent, like in DASSL.
      EXTERNAL GEVAL
C-----------------------------------------------------------------------
      SQUR=SQRT(UROUND)
      IERR=0
      CALL GEVAL(NEQN,T,Y,DY,AUX,IERR,RPAR,IPAR)
      IF (IERR.NE.0) THEN
         IDID=-2
         PRINT *, 'PSIDE: ERROR: ',
     +            'cannot handle IERR = -1 in numerical Jacobians'
         RETURN
      ENDIF
      DO 5 KN=1,NEQN
         AUX2(KN)=DY(KN)
    5 CONTINUE
      DO 30 KN=1,NEQN
         DELTA=SQUR*MAX(ABS(DY(KN)),GAM)
         AUX2(KN)=DY(KN)+DELTA
         IERR=0
         CALL GEVAL(NEQN,T,Y,AUX2,AUX3,IERR,RPAR,IPAR)
         IF (IERR.NE.0) THEN
            IDID=-2
            PRINT *, 'PSIDE: ERROR ',
     +               'cannot handle IERR = -1 in numerical Jacobians'
            RETURN
         ENDIF
         IF (MBND) THEN
            DO 10 KN2=MAX(1,KN-NUM),MIN(NEQN,KN+NLM)
               DGDDY(KN2-KN+NUM+1,KN)=(AUX3(KN2)-AUX(KN2))/DELTA
   10       CONTINUE
         ELSE
            DO 20 KN2=1,NEQN
               DGDDY(KN2,KN)=(AUX3(KN2)-AUX(KN2))/DELTA
   20       CONTINUE
         ENDIF
         AUX2(KN)=DY(KN)
   30 CONTINUE
      RETURN
      END
C$Id: newton.f,v 1.14 1998/11/25 08:58:25 walter Exp $
      SUBROUTINE NEWTON(LDM,LDLU,NEQN,Y,YS,DYS,DYSP,T,HP,H,
     +                  GEVAL,JBND,MBND,NLJ,NUJ,NLM,NUM,IERR,RPAR,IPAR,
     +                  TOLVEC,RTOL,ATOL,INDGT1,IND,
     +                  UROUND,NF,NFB,DECLUS,IPVTS,DGDDY,UYS,
     +                  UDYS,AUXS,ALPHA,GROWTH,DIVER,SLOW,SOLVED,EXACT)
CF90  IMPLICIT NONE
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER LDM,LDLU,NEQN,NLJ,NUJ,NLM,NUM,IERR,IND(*),
     +        IPAR(*),NF,NFB,IPVTS(NEQN,S)
      DOUBLE PRECISION Y(NEQN),YS(NEQN,S),DYS(NEQN,S),
     +                 DYSP(NEQN,S),T,HP,H,RPAR(*),RTOL(*),ATOL(*),
     +                 UROUND,DECLUS(LDLU,NEQN,S),
     +                 DGDDY(LDM,NEQN),UYS(NEQN,S),
     +                 UDYS(NEQN,S),AUXS(NEQN,S),ALPHA
      LOGICAL JBND,MBND,TOLVEC,INDGT1,GROWTH,DIVER,SLOW,SOLVED,EXACT
      EXTERNAL GEVAL
CF90  INTENT(IN)    LDM,LDLU,NEQN,Y,DYSP,T,HP,H,JBND,MBND,NLJ,NUJ,
CF90 +              NLM,NUM,RPAR,IPAR,TOLVEC,RTOL,ATOL,INDGT1,IND,
CF90 +              UROUND,DECLUS,IPVTS,DGDDY
CF90  INTENT(INOUT) IERR,NF,NFB
CF90  INTENT(OUT)   YS,DYS,UYS,UDYS,AUXS,ALPHA,GROWTH,DIVER,SLOW,SOLVED,
CF90 +              EXACT
C-----------------------------------------------------------------------
C     performs the Newton iteration
C-----------------------------------------------------------------------
C begin included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
      DOUBLE PRECISION A(S,S),B(S,S),C(S),D(S),Q(S,S),QINV(S,S),B0,V(S)
C end included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
      INTEGER K,KS,KS2,KN
      LOGICAL READY
C begin included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
      DATA C(1) /8.85879595127039430D-02/
      DATA C(2) /4.09466864440734710D-01/
      DATA C(3) /7.87659461760847110D-01/
      DATA C(4) /1.00000000000000000D0/

      DATA A(1,1) /1.12999479323156150D-01/
      DATA A(1,2) /-4.03092207235221240D-02/
      DATA A(1,3) /2.58023774203363160D-02/
      DATA A(1,4) /-9.90467650726640190D-03/
      DATA A(2,1) /2.34383995747400210D-01/
      DATA A(2,2) /2.06892573935359120D-01/
      DATA A(2,3) /-4.78571280485409760D-02/
      DATA A(2,4) /1.60474228065163730D-02/
      DATA A(3,1) /2.16681784623249830D-01/
      DATA A(3,2) /4.06123263867374850D-01/
      DATA A(3,3) /1.89036518170054850D-01/
      DATA A(3,4) /-2.41821048998323020D-02/
      DATA A(4,1) /2.20462211176767560D-01/
      DATA A(4,2) /3.88193468843174740D-01/
      DATA A(4,3) /3.28844319980056810D-01/
      DATA A(4,4) /6.25000000000008880D-02/

      DATA D(1) /1.52077368976571730D-01/
      DATA D(2) /1.98631665602052860D-01/
      DATA D(3) /1.73704821245586140D-01/
      DATA D(4) /2.26879766524847200D-01/

      DATA B(1,1) / -3.36398745680207070D+00/
      DATA B(1,2) / -4.46547007540097850D-01/
      DATA B(1,3) /0D0/
      DATA B(1,4) /0D0/
      DATA B(2,1) /  2.53420388412422460D+01/
      DATA B(2,2) /  3.36398745680206800D+00/
      DATA B(2,3) /0D0/
      DATA B(2,4) /0D0/
      DATA B(3,1) /0D0/
      DATA B(3,2) /0D0/
      DATA B(3,3) / -4.37367276825310510D-01/
      DATA B(3,4) / -5.80576031184040710D-02/
      DATA B(4,1) /0D0/
      DATA B(4,2) /0D0/
      DATA B(4,3) /  3.29483348541735400D+00/
      DATA B(4,4) /  4.37367276825312510D-01/

      DATA Q(1,1) /2.95256291376448490D+00/
      DATA Q(1,2) /3.15939676526544260D-01/
      DATA Q(1,3) /1.53250361857358990D+00/
      DATA Q(1,4) /2.76001773070919800D-02/
      DATA Q(2,1) /-7.26638442522609210D+00/
      DATA Q(2,2) /-8.75577120872169210D-01/
      DATA Q(2,3) /-1.05525925554083820D+00/
      DATA Q(2,4) /-3.11277680445624430D-01/
      DATA Q(3,1) /3.42060134189704890D+00/
      DATA Q(3,2) /9.49331950091266920D-01/
      DATA Q(3,3) /-1.07997190626525940D+01/
      DATA Q(3,4) /-2.13491394363750460D+00/
      DATA Q(4,1) /3.48973092842014550D+01/
      DATA Q(4,2) /4.37528029636525950D+00/
      DATA Q(4,3) /-4.29039265780423160D+01/
      DATA Q(4,4) /-5.89600020104458620D+00/

      DATA QINV(1,1) /4.94042191453623210D-01/
      DATA QINV(1,2) /2.69406327540352320D-01/
      DATA QINV(1,3) /-2.07753732935469560D-01/
      DATA QINV(1,4) /6.33161132809504090D-02/
      DATA QINV(2,1) /-3.53358212942913720D+00/
      DATA QINV(2,2) /-2.98582140519829590D+00/
      DATA QINV(2,3) /1.75646748288234010D+00/
      DATA QINV(2,4) /-4.94914306872105080D-01/
      DATA QINV(3,1) /4.87641455081039950D-01/
      DATA QINV(3,2) /1.23938205146711190D-01/
      DATA QINV(3,3) /4.23770339324015460D-02/
      DATA QINV(3,4) /-1.96050751500893220D-02/
      DATA QINV(4,1) /-3.24650638473406870D+00/
      DATA QINV(4,2) /-1.52301305545598620D+00/
      DATA QINV(4,3) /-2.34591215977400570D-01/
      DATA QINV(4,4) /-1.94525303087971780D-02/

      DATA B0 /  1.00000000000000000D-02/

      DATA V(1) /  1.57753763977411530D-02/
      DATA V(2) / -9.73676595200762000D-03/
      DATA V(3) /  6.46138955426500680D-03/
      DATA V(4) /  2.24379766524848060D-01/
C end included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
      READY=.FALSE.
      K=0
      CALL PREDIC(NEQN,DYS,DYSP,H,HP)
      DO 40 KS=1,S
         DO 10 KN=1,NEQN
            YS(KN,KS)=Y(KN)
   10    CONTINUE
         DO 30 KS2=1,S
            DO 20 KN=1,NEQN
               YS(KN,KS)=YS(KN,KS)+H*A(KS,KS2)*DYS(KN,KS2)
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE
      CALL VERGEN(NEQN,Y,YS,UYS,H,K,
     +            TOLVEC,RTOL,ATOL,INDGT1,IND,UROUND,
     +            ALPHA,READY,GROWTH,DIVER,SLOW,SOLVED,EXACT)
   50 IF (.NOT.READY) THEN
         K=K+1
C-----------------------------------------------------------------------
C we temporarily use UYS as work array for PILSRK
C-----------------------------------------------------------------------
         CALL PILSRK(LDM,LDLU,NEQN,YS,DYS,GEVAL,T,H,NF,NFB,
     +               DECLUS,IPVTS,DGDDY,JBND,NLJ,NUJ,MBND,NLM,NUM,
     +               UDYS,UYS,AUXS,IERR,RPAR,IPAR,INDGT1)
         IF (IERR.EQ.-1) RETURN
         DO 90 KS=1,S
            DO 60 KN=1,NEQN
               UYS(KN,KS)=0D0
   60       CONTINUE
            DO 80 KS2=1,S
               DO 70 KN=1,NEQN
                  UYS(KN,KS)=UYS(KN,KS)+H*A(KS,KS2)*UDYS(KN,KS2)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         DO 110 KS=1,S
            DO 100 KN=1,NEQN
               YS(KN,KS)=YS(KN,KS)+UYS(KN,KS)
               DYS(KN,KS)=DYS(KN,KS)+UDYS(KN,KS)
  100       CONTINUE
  110    CONTINUE
         CALL VERGEN(NEQN,Y,YS,UYS,H,K,
     +               TOLVEC,RTOL,ATOL,INDGT1,IND,UROUND,
     +               ALPHA,READY,GROWTH,DIVER,SLOW,SOLVED,EXACT)
         GOTO 50
      ENDIF
      RETURN
      END
C$Id: predic.f,v 1.7 1998/11/25 08:58:29 walter Exp $
      SUBROUTINE PREDIC(NEQN,DYS,DYSP,H,HP)
CF90  IMPLICIT NONE
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER NEQN
      DOUBLE PRECISION DYS(NEQN,S),DYSP(NEQN,S),H,HP
CF90  INTENT(IN)    NEQN,DYSP,H,HP
CF90  INTENT(OUT)   DYS
C-----------------------------------------------------------------------
C     extrapolation of order S on DYSP resulting in DYS
C-----------------------------------------------------------------------
C begin included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
      DOUBLE PRECISION A(S,S),B(S,S),C(S),D(S),Q(S,S),QINV(S,S),B0,V(S)
C end included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
      INTEGER KS,KS2,KS3,KN
      DOUBLE PRECISION RH,VM(S,S),W(S,S),UINV(S,S)
C begin included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
      DATA C(1) /8.85879595127039430D-02/
      DATA C(2) /4.09466864440734710D-01/
      DATA C(3) /7.87659461760847110D-01/
      DATA C(4) /1.00000000000000000D0/

      DATA A(1,1) /1.12999479323156150D-01/
      DATA A(1,2) /-4.03092207235221240D-02/
      DATA A(1,3) /2.58023774203363160D-02/
      DATA A(1,4) /-9.90467650726640190D-03/
      DATA A(2,1) /2.34383995747400210D-01/
      DATA A(2,2) /2.06892573935359120D-01/
      DATA A(2,3) /-4.78571280485409760D-02/
      DATA A(2,4) /1.60474228065163730D-02/
      DATA A(3,1) /2.16681784623249830D-01/
      DATA A(3,2) /4.06123263867374850D-01/
      DATA A(3,3) /1.89036518170054850D-01/
      DATA A(3,4) /-2.41821048998323020D-02/
      DATA A(4,1) /2.20462211176767560D-01/
      DATA A(4,2) /3.88193468843174740D-01/
      DATA A(4,3) /3.28844319980056810D-01/
      DATA A(4,4) /6.25000000000008880D-02/

      DATA D(1) /1.52077368976571730D-01/
      DATA D(2) /1.98631665602052860D-01/
      DATA D(3) /1.73704821245586140D-01/
      DATA D(4) /2.26879766524847200D-01/

      DATA B(1,1) / -3.36398745680207070D+00/
      DATA B(1,2) / -4.46547007540097850D-01/
      DATA B(1,3) /0D0/
      DATA B(1,4) /0D0/
      DATA B(2,1) /  2.53420388412422460D+01/
      DATA B(2,2) /  3.36398745680206800D+00/
      DATA B(2,3) /0D0/
      DATA B(2,4) /0D0/
      DATA B(3,1) /0D0/
      DATA B(3,2) /0D0/
      DATA B(3,3) / -4.37367276825310510D-01/
      DATA B(3,4) / -5.80576031184040710D-02/
      DATA B(4,1) /0D0/
      DATA B(4,2) /0D0/
      DATA B(4,3) /  3.29483348541735400D+00/
      DATA B(4,4) /  4.37367276825312510D-01/

      DATA Q(1,1) /2.95256291376448490D+00/
      DATA Q(1,2) /3.15939676526544260D-01/
      DATA Q(1,3) /1.53250361857358990D+00/
      DATA Q(1,4) /2.76001773070919800D-02/
      DATA Q(2,1) /-7.26638442522609210D+00/
      DATA Q(2,2) /-8.75577120872169210D-01/
      DATA Q(2,3) /-1.05525925554083820D+00/
      DATA Q(2,4) /-3.11277680445624430D-01/
      DATA Q(3,1) /3.42060134189704890D+00/
      DATA Q(3,2) /9.49331950091266920D-01/
      DATA Q(3,3) /-1.07997190626525940D+01/
      DATA Q(3,4) /-2.13491394363750460D+00/
      DATA Q(4,1) /3.48973092842014550D+01/
      DATA Q(4,2) /4.37528029636525950D+00/
      DATA Q(4,3) /-4.29039265780423160D+01/
      DATA Q(4,4) /-5.89600020104458620D+00/

      DATA QINV(1,1) /4.94042191453623210D-01/
      DATA QINV(1,2) /2.69406327540352320D-01/
      DATA QINV(1,3) /-2.07753732935469560D-01/
      DATA QINV(1,4) /6.33161132809504090D-02/
      DATA QINV(2,1) /-3.53358212942913720D+00/
      DATA QINV(2,2) /-2.98582140519829590D+00/
      DATA QINV(2,3) /1.75646748288234010D+00/
      DATA QINV(2,4) /-4.94914306872105080D-01/
      DATA QINV(3,1) /4.87641455081039950D-01/
      DATA QINV(3,2) /1.23938205146711190D-01/
      DATA QINV(3,3) /4.23770339324015460D-02/
      DATA QINV(3,4) /-1.96050751500893220D-02/
      DATA QINV(4,1) /-3.24650638473406870D+00/
      DATA QINV(4,2) /-1.52301305545598620D+00/
      DATA QINV(4,3) /-2.34591215977400570D-01/
      DATA QINV(4,4) /-1.94525303087971780D-02/

      DATA B0 /  1.00000000000000000D-02/

      DATA V(1) /  1.57753763977411530D-02/
      DATA V(2) / -9.73676595200762000D-03/
      DATA V(3) /  6.46138955426500680D-03/
      DATA V(4) /  2.24379766524848060D-01/
C end included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
C begin included $Id: predic.h,v 1.4 1997/11/25 17:11:50 walter Rel $
      DATA UINV(1,1) /   0D0                  /
      DATA UINV(1,2) /   0D0                  /
      DATA UINV(1,3) /   0D0                  /
      DATA UINV(1,4) /   1D0                  /
      DATA UINV(2,1) /  -0.6133376973486749D0 /
      DATA UINV(2,2) /   2.700531288824063D0  /
      DATA UINV(2,3) /  -9.587193591475394D0  /
      DATA UINV(2,4) /   7.5D0                /
      DATA UINV(3,1) /  -3.927079477247392D0  /
      DATA UINV(3,2) /  15.68094527819332D0   /
      DATA UINV(3,3) / -26.75386580094595D0   /
      DATA UINV(3,4) /  15D0                  /
      DATA UINV(4,1) /  -4.891279419672913D0  /
      DATA UINV(4,2) /  13.95409058457028D0   /
      DATA UINV(4,3) / -17.81281116489738D0   /
      DATA UINV(4,4) /   8.75D0               /
C end included $Id: predic.h,v 1.4 1997/11/25 17:11:50 walter Rel $
      RH=H/HP
      DO 20 KS=1,S
         DO 10 KS2=1,S
            VM(KS2,KS)=(RH*C(KS2))**(KS-1)
   10    CONTINUE
   20 CONTINUE
      DO 60 KS=1,S
         DO 30 KS2=1,S
            W(KS2,KS)=0D0
   30    CONTINUE
         DO 50 KS3=1,S
            DO 40 KS2=1,S
               W(KS2,KS)=W(KS2,KS)+VM(KS2,KS3)*UINV(KS3,KS)
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      DO 100 KS=1,S
         DO 70 KN=1,NEQN
            DYS(KN,KS)=DYSP(KN,S)
   70    CONTINUE
         DO 90 KS2=1,S
            DO 80 KN=1,NEQN
               DYS(KN,KS)=DYS(KN,KS)+W(KS,KS2)*(DYSP(KN,KS2)-DYSP(KN,S))
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
C$Id: pilsrk.f,v 1.18 1998/11/25 08:58:28 walter Exp $
      SUBROUTINE PILSRK(LDM,LDLU,NEQN,YS,DYS,GEVAL,T,H,NF,NFB,
     +                  DECLUS,IPVTS,DGDDY,JBND,NLJ,NUJ,MBND,NLM,NUM,
     +                  UDYS,TUDYS,TGS,IERR,RPAR,IPAR,INDGT1)
CF90  IMPLICIT NONE
C      include 'omp_lib.h'
      use omp_lib
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER LDM,LDLU,NEQN,NF,NFB,IPVTS(NEQN,S),NLJ,NUJ,NLM,NUM,
     +        IERR,IPAR(*)
      DOUBLE PRECISION YS(NEQN,S),DYS(NEQN,S),T,H,
     +                 DECLUS(LDLU,NEQN,S),DGDDY(LDM,NEQN),
     +                 UDYS(NEQN,S),TUDYS(NEQN,S),TGS(NEQN,S),RPAR(*)
      LOGICAL JBND,MBND,INDGT1
      EXTERNAL GEVAL
CF90  INTENT(IN)    LDM,LDLU,NEQN,YS,DYS,T,H,DECLUS,IPVTS,
CF90 +              DGDDY,JBND,NLJ,NUJ,MBND,NLM,NUM,IPAR,RPAR,INDGT1
CF90  INTENT(INOUT) NF,NFB,IERR
CF90  INTENT(OUT)   UDYS,TUDYS,TGS
C-----------------------------------------------------------------------
C     the Parallel Iterative Linear system Solver for
C     Runge-Kutta methods
C-----------------------------------------------------------------------
      INTEGER INFO(S),IERRS(S),M,KM,KB,KS,KS2,KN,KN2
C begin included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
      DOUBLE PRECISION A(S,S),B(S,S),C(S),D(S),Q(S,S),QINV(S,S),B0,V(S)
C end included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
C begin included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
      DATA C(1) /8.85879595127039430D-02/
      DATA C(2) /4.09466864440734710D-01/
      DATA C(3) /7.87659461760847110D-01/
      DATA C(4) /1.00000000000000000D0/

      DATA A(1,1) /1.12999479323156150D-01/
      DATA A(1,2) /-4.03092207235221240D-02/
      DATA A(1,3) /2.58023774203363160D-02/
      DATA A(1,4) /-9.90467650726640190D-03/
      DATA A(2,1) /2.34383995747400210D-01/
      DATA A(2,2) /2.06892573935359120D-01/
      DATA A(2,3) /-4.78571280485409760D-02/
      DATA A(2,4) /1.60474228065163730D-02/
      DATA A(3,1) /2.16681784623249830D-01/
      DATA A(3,2) /4.06123263867374850D-01/
      DATA A(3,3) /1.89036518170054850D-01/
      DATA A(3,4) /-2.41821048998323020D-02/
      DATA A(4,1) /2.20462211176767560D-01/
      DATA A(4,2) /3.88193468843174740D-01/
      DATA A(4,3) /3.28844319980056810D-01/
      DATA A(4,4) /6.25000000000008880D-02/

      DATA D(1) /1.52077368976571730D-01/
      DATA D(2) /1.98631665602052860D-01/
      DATA D(3) /1.73704821245586140D-01/
      DATA D(4) /2.26879766524847200D-01/

      DATA B(1,1) / -3.36398745680207070D+00/
      DATA B(1,2) / -4.46547007540097850D-01/
      DATA B(1,3) /0D0/
      DATA B(1,4) /0D0/
      DATA B(2,1) /  2.53420388412422460D+01/
      DATA B(2,2) /  3.36398745680206800D+00/
      DATA B(2,3) /0D0/
      DATA B(2,4) /0D0/
      DATA B(3,1) /0D0/
      DATA B(3,2) /0D0/
      DATA B(3,3) / -4.37367276825310510D-01/
      DATA B(3,4) / -5.80576031184040710D-02/
      DATA B(4,1) /0D0/
      DATA B(4,2) /0D0/
      DATA B(4,3) /  3.29483348541735400D+00/
      DATA B(4,4) /  4.37367276825312510D-01/

      DATA Q(1,1) /2.95256291376448490D+00/
      DATA Q(1,2) /3.15939676526544260D-01/
      DATA Q(1,3) /1.53250361857358990D+00/
      DATA Q(1,4) /2.76001773070919800D-02/
      DATA Q(2,1) /-7.26638442522609210D+00/
      DATA Q(2,2) /-8.75577120872169210D-01/
      DATA Q(2,3) /-1.05525925554083820D+00/
      DATA Q(2,4) /-3.11277680445624430D-01/
      DATA Q(3,1) /3.42060134189704890D+00/
      DATA Q(3,2) /9.49331950091266920D-01/
      DATA Q(3,3) /-1.07997190626525940D+01/
      DATA Q(3,4) /-2.13491394363750460D+00/
      DATA Q(4,1) /3.48973092842014550D+01/
      DATA Q(4,2) /4.37528029636525950D+00/
      DATA Q(4,3) /-4.29039265780423160D+01/
      DATA Q(4,4) /-5.89600020104458620D+00/

      DATA QINV(1,1) /4.94042191453623210D-01/
      DATA QINV(1,2) /2.69406327540352320D-01/
      DATA QINV(1,3) /-2.07753732935469560D-01/
      DATA QINV(1,4) /6.33161132809504090D-02/
      DATA QINV(2,1) /-3.53358212942913720D+00/
      DATA QINV(2,2) /-2.98582140519829590D+00/
      DATA QINV(2,3) /1.75646748288234010D+00/
      DATA QINV(2,4) /-4.94914306872105080D-01/
      DATA QINV(3,1) /4.87641455081039950D-01/
      DATA QINV(3,2) /1.23938205146711190D-01/
      DATA QINV(3,3) /4.23770339324015460D-02/
      DATA QINV(3,4) /-1.96050751500893220D-02/
      DATA QINV(4,1) /-3.24650638473406870D+00/
      DATA QINV(4,2) /-1.52301305545598620D+00/
      DATA QINV(4,3) /-2.34591215977400570D-01/
      DATA QINV(4,4) /-1.94525303087971780D-02/

      DATA B0 /  1.00000000000000000D-02/

      DATA V(1) /  1.57753763977411530D-02/
      DATA V(2) / -9.73676595200762000D-03/
      DATA V(3) /  6.46138955426500680D-03/
      DATA V(4) /  2.24379766524848060D-01/
      !call omp_set_dynamic(.true.)
C end included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
      IF (INDGT1) THEN
         M=2
      ELSE
         M=1
      ENDIF
C-----------------------------------------------------------------------
C     we temporarily use UDYS for storing G
C-----------------------------------------------------------------------
      IERR=0
C      DO 10 KS=1,S
C$OMP PARALLEL PRIVATE(KS)
C$OMP&      SHARED(NEQN,T,C,H,YS,DYS,UDYS,IERRS,RPAR,IPAR)
C$OMP DO
      DO 10 KS=1,S
         IERRS(KS)=0
         CALL GEVAL(NEQN,T+C(KS)*H,YS(1,KS),DYS(1,KS),UDYS(1,KS),
     +              IERRS(KS),RPAR,IPAR)
  10  CONTINUE
C$OMP END PARALLEL
      NF=NF+S
      DO 15 KS=1,S
         IF (IERRS(KS).EQ.-1) THEN
            IERR=-1
            RETURN
         ENDIF
  15  CONTINUE
C-----------------------------------------------------------------------
C     first PILSRK iteration:
C     transform GS / solve linear systems in one loop
C     to make the parallel loops as large as possible
C-----------------------------------------------------------------------
      IF (JBND) THEN
C         DO 60 KS=1,S
C$OMP PARALLEL PRIVATE(KS,KN,KS2)
C$OMP&      SHARED(TGS,QINV,UDYS,
C$OMP&             NEQN,NLJ,NUJ,DECLUS,LDLU,IPVTS,TUDYS,INFO)
C$OMP DO
         DO 60 KS=1,S
            DO 20 KN=1,NEQN
               TGS(KN,KS)=0D0
   20       CONTINUE
            DO 40 KS2=1,S
               DO 30 KN=1,NEQN
                  TGS(KN,KS)=TGS(KN,KS)-QINV(KS,KS2)*UDYS(KN,KS2)
   30          CONTINUE
   40       CONTINUE
            DO 50 KN=1,NEQN
               TUDYS(KN,KS)=TGS(KN,KS)
   50       CONTINUE
C           we do not have to check the INFO value on exit:
C           INFO returns 0 (no argument with illegal value)
            CALL DGBTRS('n',NEQN,NLJ,NUJ,1,DECLUS(1,1,KS),LDLU,
     +                  IPVTS(1,KS),TUDYS(1,KS),NEQN,INFO(KS))
   60    CONTINUE
C$omp END PARALLEL
      ELSE
C         DO 110 KS=1,S
C$OMP PARALLEL PRIVATE(KS,KN,KS2)
C$OMP&      SHARED(TGS,QINV,UDYS,
C$OMP&             NEQN,DECLUS,LDLU,IPVTS,TUDYS,INFO)
C$OMP DO
         DO 110 KS=1,S
C           DO 70 KN=1,NEQN
C              TGS(KN,KS)=0D0
C  70       CONTINUE
C           DO 90 KS2=1,S
C              DO 80 KN=1,NEQN
C                 TGS(KN,KS)=TGS(KN,KS)-QINV(KS,KS2)*UDYS(KN,KS2)
C  80          CONTINUE
C  90       CONTINUE
C           DO 100 KN=1,NEQN
C              TUDYS(KN,KS)=TGS(KN,KS)
C 100       CONTINUE
C
C loop 80 and 90 interchanged;
C unrolled inner loop (S=4);
C fused with 70 and 100
            DO 80 KN=1,NEQN
               TGS(KN,KS)=-QINV(KS,1)*UDYS(KN,1)
     +                    -QINV(KS,2)*UDYS(KN,2)
     +                    -QINV(KS,3)*UDYS(KN,3)
     +                    -QINV(KS,4)*UDYS(KN,4)
               TUDYS(KN,KS)=TGS(KN,KS)
   80       CONTINUE
C           we do not have to check the INFO value on exit:
C           INFO returns 0 (no argument with illegal value)
            CALL DGETRS('n',NEQN,1,DECLUS(1,1,KS),LDLU,IPVTS(1,KS),
     +                  TUDYS(1,KS),NEQN,INFO(KS))
  110    CONTINUE
C$OMP END PARALLEL
      ENDIF
      NFB=NFB+S
C-----------------------------------------------------------------------
C     second upto M-th PILSRK iterations
C-----------------------------------------------------------------------
      DO 280 KM=2,M
C-----------------------------------------------------------------------
C     we temporarily use UDYS for storing -(LU)^(-1)((B x M)TUDYS+TGS)
C-----------------------------------------------------------------------
         IF (JBND) THEN
C            DO 150 KS=1,S
C$OMP PARALLEL PRIVATE(KS,KN,KN2,KB)
C$OMP&      SHARED(TGS,NLM,NUM,DGDDY,B,TUDYS,
C$OMP&             NEQN,NLJ,NUJ,DECLUS,LDLU,IPVTS,UDYS,INFO)
C$OMP DO
            DO 150 KS=1,S
               KB=KS-MOD(KS+1,2)
               DO 120 KN=1,NEQN
                  UDYS(KN,KS)=TGS(KN,KS)
  120          CONTINUE
               DO 140 KN=1,NEQN
                  DO 130 KN2=MAX(1,KN-NUM),MIN(NEQN,KN+NLM)
                     UDYS(KN2,KS)=UDYS(KN2,KS)-DGDDY(NUM+1+KN2-KN,KN)*
     +                            (B(KS,KB)*TUDYS(KN,KB)+
     +                             B(KS,KB+1)*TUDYS(KN,KB+1))
  130             CONTINUE
  140          CONTINUE
C              we do not have to check the INFO value on exit:
C              INFO returns 0 (no argument with illegal value)
               CALL DGBTRS('n',NEQN,NLJ,NUJ,1,DECLUS(1,1,KS),LDLU,
     +                     IPVTS(1,KS),UDYS(1,KS),NEQN,INFO(KS))
  150       CONTINUE
C$OMP END PARALLEL 
         ELSEIF (MBND) THEN
C            DO 190 KS=1,S
C$OMP PARALLEL PRIVATE(KS,KN,KN2,KB)
C$OMP&      SHARED(TGS,NLM,NUM,DGDDY,B,TUDYS,
C$OMP&             NEQN,DECLUS,LDLU,IPVTS,UDYS,INFO)
C$OMP DO
            DO 190 KS=1,S
               KB=KS-MOD(KS+1,2)
               DO 160 KN=1,NEQN
                  UDYS(KN,KS)=TGS(KN,KS)
  160          CONTINUE
               DO 180 KN=1,NEQN
                  DO 170 KN2=MAX(1,KN-NUM),MIN(NEQN,KN+NLM)
                     UDYS(KN2,KS)=UDYS(KN2,KS)-DGDDY(NUM+1+KN2-KN,KN)*
     +                            (B(KS,KB)*TUDYS(KN,KB)+
     +                             B(KS,KB+1)*TUDYS(KN,KB+1))
  170             CONTINUE
  180          CONTINUE
C              we do not have to check the INFO value on exit:
C              INFO returns 0 (no argument with illegal value)
               CALL DGETRS('n',NEQN,1,DECLUS(1,1,KS),LDLU,IPVTS(1,KS),
     +                     UDYS(1,KS),NEQN,INFO(KS))
  190       CONTINUE
C$OMP END PARALLEL
         ELSE
C            DO 230 KS=1,S
C$OMP PARALLEL PRIVATE(KS,KN,KN2,KB)
C$OMP&      SHARED(TGS,DGDDY,B,TUDYS,
C$OMP&             NEQN,DECLUS,LDLU,IPVTS,UDYS,INFO)
C$OMP DO
            DO 230 KS=1,S
               KB=KS-MOD(KS+1,2)
               DO 200 KN=1,NEQN
                  UDYS(KN,KS)=TGS(KN,KS)
  200          CONTINUE
               DO 220 KN=1,NEQN
                  DO 210 KN2=1,NEQN
                     UDYS(KN2,KS)=UDYS(KN2,KS)-DGDDY(KN2,KN)*
     +                            (B(KS,KB)*TUDYS(KN,KB)+
     +                             B(KS,KB+1)*TUDYS(KN,KB+1))
  210             CONTINUE
  220          CONTINUE
C              we do not have to check the INFO value on exit:
C              INFO returns 0 (no argument with illegal value)
               CALL DGETRS('n',NEQN,1,DECLUS(1,1,KS),LDLU,IPVTS(1,KS),
     +                     UDYS(1,KS),NEQN,INFO(KS))
  230       CONTINUE
C$OMP END PARALLEL
         ENDIF
         DO 250 KS=1,S
            KB=KS-MOD(KS+1,2)
            DO 240 KN=1,NEQN
               UDYS(KN,KS)=UDYS(KN,KS)+
     +                     (B(KS,KB)*TUDYS(KN,KB)+
     +                      B(KS,KB+1)*TUDYS(KN,KB+1))
  240       CONTINUE
  250    CONTINUE
         DO 270 KS=1,S
            DO 260 KN=1,NEQN
               TUDYS(KN,KS)=UDYS(KN,KS)
  260       CONTINUE
  270    CONTINUE
         NFB=NFB+S
  280 CONTINUE
C-----------------------------------------------------------------------
C     transform TUDYS to UDYS
C-----------------------------------------------------------------------
      DO 320 KS=1,S
C        DO 290 KN=1,NEQN
C           UDYS(KN,KS)=0D0
C 290    CONTINUE
C        DO 310 KS2=1,S
C           DO 300 KN=1,NEQN
C              UDYS(KN,KS)=UDYS(KN,KS)+Q(KS,KS2)*TUDYS(KN,KS2)
C 300       CONTINUE
C 310    CONTINUE
C
C loop 300 and 310 interchanged;
C unrolled inner loop (S=4);
C fused with 290
         DO 300 KN=1,NEQN
            UDYS(KN,KS)=+Q(KS,1)*TUDYS(KN,1)
     +                  +Q(KS,2)*TUDYS(KN,2)
     +                  +Q(KS,3)*TUDYS(KN,3)
     +                  +Q(KS,4)*TUDYS(KN,4)
  300    CONTINUE
  320 CONTINUE
      RETURN
      END
C$Id: vergen.f,v 1.17 1998/11/25 08:58:35 walter Exp $
      SUBROUTINE VERGEN(NEQN,Y,YS,UYS,H,K,
     +                  TOLVEC,RTOL,ATOL,INDGT1,IND,UROUND,
     +                  ALPHA,READY,GROWTH,DIVER,SLOW,SOLVED,EXACT)
CF90  IMPLICIT NONE
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER NEQN,K,IND(*)
      DOUBLE PRECISION Y(NEQN),YS(NEQN,S),UYS(NEQN,S),H,
     +                 RTOL(*),ATOL(*),UROUND,ALPHA
      LOGICAL TOLVEC,INDGT1,READY,GROWTH,DIVER,SLOW,SOLVED,EXACT
CF90  INTENT(IN)    NEQN,Y,YS,UYS,H,K,TOLVEC,RTOL,ATOL,INDGT1,IND,UROUND
CF90  INTENT(INOUT) ALPHA
CF90  INTENT(OUT)   READY,GROWTH,DIVER,SLOW,SOLVED,EXACT
C-----------------------------------------------------------------------
C     checks whether the Newton process converges
C-----------------------------------------------------------------------
      INTEGER KMAX,KN
      DOUBLE PRECISION UP,U,SNORM,SNORMS,TAU,KAPPA,GAMMA,THETA,GFAC,
     +                 ALPHA1
      PARAMETER(TAU    = 0.01D0,
     +          KAPPA  = 100D0,
     +          KMAX   = 15,
     +          GAMMA  = 1D0,
     +          THETA  = 0.5D0,
     +          GFAC   = 100D0,
     +          ALPHA1 = 0.1D0)
      SAVE U
      GROWTH=.FALSE.
      DIVER =.FALSE.
      SLOW  =.FALSE.
      SOLVED=.FALSE.
      EXACT =.FALSE.
      IF (TOLVEC) THEN
         IF (INDGT1) THEN
            DO 10 KN=1,NEQN
               GROWTH = GROWTH .OR.
     +                  IND(KN).LE.1 .AND.
     +                  ABS(YS(KN,S)).GT.
     +                  GFAC*MAX(ABS(Y(KN)),ATOL(KN))
   10       CONTINUE
         ELSE
            DO 20 KN=1,NEQN
               GROWTH = GROWTH .OR.
     +                  ABS(YS(KN,S)).GT.
     +                  GFAC*MAX(ABS(Y(KN)),ATOL(KN))
   20       CONTINUE
         ENDIF
      ELSE
         IF (INDGT1) THEN
            DO 30 KN=1,NEQN
               GROWTH = GROWTH .OR.
     +                  IND(KN).LE.1 .AND.
     +                  ABS(YS(KN,S)).GT.
     +                  GFAC*MAX(ABS(Y(KN)),ATOL(1))
   30       CONTINUE
         ELSE
            DO 40 KN=1,NEQN
               GROWTH = GROWTH .OR.
     +                  ABS(YS(KN,S)).GT.
     +                  GFAC*MAX(ABS(Y(KN)),ATOL(1))
   40       CONTINUE
         ENDIF
      ENDIF
      IF (.NOT.GROWTH) THEN
         IF (K.EQ.1) THEN
            U=SNORMS(NEQN,Y,UYS,H,TOLVEC,RTOL,ATOL,INDGT1,IND)
            ALPHA=ALPHA1
            EXACT=U.EQ.0D0
            SOLVED=EXACT
         ELSEIF (K.GT.1) THEN
            UP=U
            U=SNORMS(NEQN,Y,UYS,H,TOLVEC,RTOL,ATOL,INDGT1,IND)
            ALPHA=ALPHA**THETA*(U/UP)**(1-THETA)
            IF (ALPHA.GE.GAMMA) THEN
               DIVER=.TRUE.
            ELSEIF ((U*ALPHA.LT.(1D0-ALPHA)*TAU).OR.(U.LT.KAPPA*UROUND*
     +            SNORM(NEQN,Y,Y,H,TOLVEC,RTOL,ATOL,INDGT1,IND))) THEN
               SOLVED=.TRUE.
            ELSEIF ((K.EQ.KMAX).OR.
     +            (U*ALPHA**(KMAX-K).GT.TAU*(1D0-ALPHA))) THEN
               SLOW=.TRUE.
            ENDIF
         ENDIF
      ENDIF
      READY = GROWTH .OR. DIVER .OR. SLOW .OR. SOLVED .OR. EXACT
      RETURN
      END
C$Id: ctrl.f,v 1.18 1998/11/25 08:58:22 walter Exp $
      SUBROUTINE CTRL(LDLU,NEQN,Y,DY,T,TEND,GEVAL,
     +                JBND,NLJ,NUJ,HP,H,HLU,TOLVEC,RTOL,ATOL,
     +                INDGT1,IND,UROUND,
     +                IERR,RPAR,IPAR,
     +                NREJN,NREJE,NREJG,NREJI,NF,NFB,
     +                YS,DYS,DYSP,DECLUS,IPVTS,AUX,AUX2,ALPHA,
     +                GROWTH,DIVER,SLOW,SOLVED,EXACT,
     +                FIRST,JACNEW,FACNEW,JACU2D,IDID)
CF90  IMPLICIT NONE
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER LDLU,NEQN,NLJ,NUJ,IND(*),IERR,IPAR(*),IPVTS(NEQN,S),
     +        NREJN,NREJE,NREJG,NREJI,NF,NFB,IDID
      DOUBLE PRECISION Y(NEQN),DY(NEQN),T,TEND,HP,H,HLU,
     +                 RTOL(*),ATOL(*),RPAR(*),UROUND,
     +                 YS(NEQN,S),DYS(NEQN,S),DYSP(NEQN,S),
     +                 DECLUS(LDLU,NEQN,S),AUX(NEQN),AUX2(NEQN),ALPHA
      LOGICAL TOLVEC,INDGT1,JBND,FIRST,JACNEW,FACNEW,JACU2D,
     +        GROWTH,DIVER,SLOW,SOLVED,EXACT
      EXTERNAL GEVAL
CF90  INTENT(IN)    LDLU,NEQN,TEND,JBND,NLJ,NUJ,TOLVEC,RTOL,ATOL,
CF90 +              INDGT1,IND,UROUND,RPAR,IPAR,YS,DYS,DECLUS,IPVTS,
CF90 +              ALPHA,GROWTH,DIVER,SLOW,SOLVED,EXACT
CF90  INTENT(INOUT) Y,DY,T,H,HLU,IERR,NREJN,NREJE,NREJG,NREJI,NF,NFB,
CF90 +              FIRST,JACNEW,JACU2D
CF90  INTENT(OUT)   HP,DYSP,AUX,AUX2,FACNEW,IDID
C-----------------------------------------------------------------------
C     decides whether the Jacobians should be updated
C     and whether the iteration matrix should be factorized
C-----------------------------------------------------------------------
      DOUBLE PRECISION HR,HNEW,HALPHA,ALPREF,ALPJAC,ALPLU,FMIN,FMAX,
     +                 FRIG,XI,REMAIN,OMEGA
      PARAMETER(ALPREF = 0.15D0,
     +          ALPJAC = 0.1D0,
     +          ALPLU  = 0.2D0,
     +          FMIN   = 0.2D0,
     +          FMAX   = 2D0,
     +          FRIG   = 2D0,
     +          XI     = 1.2D0,
     +          OMEGA  = 0.05D0)
      JACNEW = .FALSE.
      IF (.NOT. GROWTH) HALPHA=H*ALPREF/MAX(ALPHA,ALPREF/FMAX)
      IF (IERR.EQ.-1) THEN
         HNEW=H/FRIG
         NREJI=NREJI+1
      ELSEIF (SOLVED) THEN
         CALL ERROR(LDLU,NEQN,Y,DY,T,TEND,GEVAL,
     +              JBND,NLJ,NUJ,HP,H,HR,TOLVEC,RTOL,ATOL,
     +              INDGT1,IND,UROUND,
     +              IERR,RPAR,IPAR,NREJE,NF,NFB,
     +              YS,DYS,DYSP,DECLUS,IPVTS,AUX,AUX2,FIRST,JACU2D)
         IF (IERR.EQ.-1) THEN
            HNEW=H/FRIG
            NREJI=NREJI+1
         ELSEIF (ABS(TEND-T).GT.10D0*UROUND*ABS(T)) THEN
            IF (JACU2D.AND.ALPHA.GT.ALPREF) THEN
               HNEW=MIN(FMAX*H,MAX(FMIN*H,MIN(HR,HALPHA)))
            ELSE
               HNEW=MIN(FMAX*H,MAX(FMIN*H,HR))
            ENDIF
            IF (.NOT.EXACT.AND.ALPHA-ABS(H-HLU)/HLU.GT.ALPJAC) THEN
               IF (JACU2D) THEN
                  HNEW=H/FRIG
               ELSE
                  JACNEW=.TRUE.
               ENDIF
            ENDIF
         ENDIF
      ELSEIF (GROWTH) THEN
         HNEW=H/FRIG
         NREJG=NREJG+1
      ELSEIF (DIVER) THEN
         HNEW=MIN(FMAX*H,MAX(FMIN*H,HALPHA))
         JACNEW=.NOT.JACU2D
         NREJN=NREJN+1
      ELSEIF (SLOW) THEN
         IF (JACU2D) THEN
            IF (ALPHA.GT.XI*ALPREF) THEN
               HNEW=MIN(FMAX*H,MAX(FMIN*H,HALPHA))
            ELSE
               HNEW=H/FRIG
            ENDIF
         ELSE
            HNEW=H
            JACNEW=.TRUE.
         ENDIF
         NREJN=NREJN+1
      ELSE
         PRINT *, 'PSIDE: ERROR: impossible(?) error in CTRL'
         STOP
      ENDIF
      IF (ABS(TEND-T).GT.10D0*UROUND*ABS(T)) THEN
         IF (HNEW.LT.10D0*UROUND*ABS(T)) THEN
            IDID=-1
            RETURN
         ENDIF
         REMAIN=(TEND-T)/HNEW
         IF (REMAIN-AINT(REMAIN).GT.OMEGA .OR. AINT(REMAIN).EQ.0D0) THEN
            REMAIN=AINT(REMAIN)+1D0
         ELSE
            REMAIN=AINT(REMAIN)
         ENDIF
         H=(TEND-T)/REMAIN
      ENDIF
      FACNEW=JACNEW.OR.ABS(H-HLU).GT.ALPLU*HLU
      IF (FACNEW) HLU=H
      RETURN
      END
C$Id: error.f,v 1.21 1998/11/25 08:58:23 walter Exp $
      SUBROUTINE ERROR(LDLU,NEQN,Y,DY,T,TEND,GEVAL,
     +                 JBND,NLJ,NUJ,HP,H,HR,TOLVEC,RTOL,ATOL,
     +                 INDGT1,IND,UROUND,
     +                 IERR,RPAR,IPAR,NREJE,NF,NFB,
     +                 YS,DYS,DYSP,DECLUS,IPVTS,AUX,R,
     +                 FIRST,JACU2D)
CF90  IMPLICIT NONE
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER LDLU,NEQN,NLJ,NUJ,IND(*),IERR,IPAR(*),IPVTS(NEQN,S),
     +        NREJE,NF,NFB
      DOUBLE PRECISION Y(NEQN),DY(NEQN),T,TEND,HP,H,HR,
     +                 RTOL(*),ATOL(*),UROUND,RPAR(*),
     +                 YS(NEQN,S),DYS(NEQN,S),DYSP(NEQN,S),
     +                 DECLUS(LDLU,NEQN,S),AUX(NEQN),R(NEQN)
      LOGICAL JBND,TOLVEC,INDGT1,FIRST,JACU2D
      EXTERNAL GEVAL
CF90  INTENT(IN)    LDLU,NEQN,TEND,JBND,NLJ,NUJ,H,TOLVEC,RTOL,ATOL,
CF90 +              INDGT1,IND,UROUND,RPAR,IPAR,YS,DYS,DECLUS,IPVTS
CF90  INTENT(INOUT) Y,DY,T,HP,IERR,NREJE,NF,NFB,FIRST
CF90  INTENT(OUT)   HR,DYSP,AUX,R,JACU2D
C-----------------------------------------------------------------------
C     computes the local error estimate,
C     proposes a new step size
C     decides whether the current step can be accepted,
C     if so, it shifts iterates.
C-----------------------------------------------------------------------
C begin included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
      DOUBLE PRECISION A(S,S),B(S,S),C(S),D(S),Q(S,S),QINV(S,S),B0,V(S)
C end included $Id: metdec.h,v 1.3 1997/11/25 17:11:45 walter Rel $
      INTEGER KS,KN,INFO
      DOUBLE PRECISION EPS,EPSP,HREJ,EPSREJ,PEST,ZETA,PMIN,SNORM,
     +                 DSINV,FMAX
      LOGICAL SUCREJ
      PARAMETER(ZETA = 0.8D0,
     +          PMIN = 0.1D0,
     +          FMAX = 2D0)
C     FMAX is the same as in module CTRL
      SAVE EPSP,HREJ,EPSREJ,SUCREJ
C begin included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
      DATA C(1) /8.85879595127039430D-02/
      DATA C(2) /4.09466864440734710D-01/
      DATA C(3) /7.87659461760847110D-01/
      DATA C(4) /1.00000000000000000D0/

      DATA A(1,1) /1.12999479323156150D-01/
      DATA A(1,2) /-4.03092207235221240D-02/
      DATA A(1,3) /2.58023774203363160D-02/
      DATA A(1,4) /-9.90467650726640190D-03/
      DATA A(2,1) /2.34383995747400210D-01/
      DATA A(2,2) /2.06892573935359120D-01/
      DATA A(2,3) /-4.78571280485409760D-02/
      DATA A(2,4) /1.60474228065163730D-02/
      DATA A(3,1) /2.16681784623249830D-01/
      DATA A(3,2) /4.06123263867374850D-01/
      DATA A(3,3) /1.89036518170054850D-01/
      DATA A(3,4) /-2.41821048998323020D-02/
      DATA A(4,1) /2.20462211176767560D-01/
      DATA A(4,2) /3.88193468843174740D-01/
      DATA A(4,3) /3.28844319980056810D-01/
      DATA A(4,4) /6.25000000000008880D-02/

      DATA D(1) /1.52077368976571730D-01/
      DATA D(2) /1.98631665602052860D-01/
      DATA D(3) /1.73704821245586140D-01/
      DATA D(4) /2.26879766524847200D-01/

      DATA B(1,1) / -3.36398745680207070D+00/
      DATA B(1,2) / -4.46547007540097850D-01/
      DATA B(1,3) /0D0/
      DATA B(1,4) /0D0/
      DATA B(2,1) /  2.53420388412422460D+01/
      DATA B(2,2) /  3.36398745680206800D+00/
      DATA B(2,3) /0D0/
      DATA B(2,4) /0D0/
      DATA B(3,1) /0D0/
      DATA B(3,2) /0D0/
      DATA B(3,3) / -4.37367276825310510D-01/
      DATA B(3,4) / -5.80576031184040710D-02/
      DATA B(4,1) /0D0/
      DATA B(4,2) /0D0/
      DATA B(4,3) /  3.29483348541735400D+00/
      DATA B(4,4) /  4.37367276825312510D-01/

      DATA Q(1,1) /2.95256291376448490D+00/
      DATA Q(1,2) /3.15939676526544260D-01/
      DATA Q(1,3) /1.53250361857358990D+00/
      DATA Q(1,4) /2.76001773070919800D-02/
      DATA Q(2,1) /-7.26638442522609210D+00/
      DATA Q(2,2) /-8.75577120872169210D-01/
      DATA Q(2,3) /-1.05525925554083820D+00/
      DATA Q(2,4) /-3.11277680445624430D-01/
      DATA Q(3,1) /3.42060134189704890D+00/
      DATA Q(3,2) /9.49331950091266920D-01/
      DATA Q(3,3) /-1.07997190626525940D+01/
      DATA Q(3,4) /-2.13491394363750460D+00/
      DATA Q(4,1) /3.48973092842014550D+01/
      DATA Q(4,2) /4.37528029636525950D+00/
      DATA Q(4,3) /-4.29039265780423160D+01/
      DATA Q(4,4) /-5.89600020104458620D+00/

      DATA QINV(1,1) /4.94042191453623210D-01/
      DATA QINV(1,2) /2.69406327540352320D-01/
      DATA QINV(1,3) /-2.07753732935469560D-01/
      DATA QINV(1,4) /6.33161132809504090D-02/
      DATA QINV(2,1) /-3.53358212942913720D+00/
      DATA QINV(2,2) /-2.98582140519829590D+00/
      DATA QINV(2,3) /1.75646748288234010D+00/
      DATA QINV(2,4) /-4.94914306872105080D-01/
      DATA QINV(3,1) /4.87641455081039950D-01/
      DATA QINV(3,2) /1.23938205146711190D-01/
      DATA QINV(3,3) /4.23770339324015460D-02/
      DATA QINV(3,4) /-1.96050751500893220D-02/
      DATA QINV(4,1) /-3.24650638473406870D+00/
      DATA QINV(4,2) /-1.52301305545598620D+00/
      DATA QINV(4,3) /-2.34591215977400570D-01/
      DATA QINV(4,4) /-1.94525303087971780D-02/

      DATA B0 /  1.00000000000000000D-02/

      DATA V(1) /  1.57753763977411530D-02/
      DATA V(2) / -9.73676595200762000D-03/
      DATA V(3) /  6.46138955426500680D-03/
      DATA V(4) /  2.24379766524848060D-01/
C end included $Id: metpar.h,v 1.7 1997/11/25 17:11:46 walter Rel $
      DO 10 KN=1,NEQN
         AUX(KN)=-B0*DY(KN)
   10 CONTINUE
      DO 30 KS=1,S
         DO 20 KN=1,NEQN
            AUX(KN)=AUX(KN)+V(KS)*DYS(KN,KS)
   20    CONTINUE
   30 CONTINUE
      DSINV=1D0/D(S)
      DO 40 KN=1,NEQN
         AUX(KN)=DSINV*AUX(KN)
   40 CONTINUE
      IERR=0
      CALL GEVAL(NEQN,T+H,YS(1,S),AUX,R,IERR,RPAR,IPAR)
      NF=NF+1
      IF (IERR.EQ.-1) RETURN
C     we do not have to check the INFO value on exit:
C     INFO returns 0 (no argument with illegal value)
      IF (JBND) THEN
         CALL DGBTRS('n',NEQN,NLJ,NUJ,1,DECLUS(1,1,S),LDLU,
     +               IPVTS(1,S),R,NEQN,INFO)
      ELSE
         CALL DGETRS('n',NEQN,1,DECLUS(1,1,S),LDLU,IPVTS(1,S),
     +               R,NEQN,INFO)
      ENDIF
      NFB=NFB+1
      DO 50 KN=1,NEQN
         R(KN)=-H*D(S)*R(KN)
   50 CONTINUE
      EPS=SNORM(NEQN,Y,R,H,TOLVEC,RTOL,ATOL,INDGT1,IND)
      IF (EPS.LT.1D0) THEN
         IF (EPS.EQ.0D0) THEN
            HR=FMAX*H
         ELSEIF (FIRST.OR.SUCREJ) THEN
            FIRST=.FALSE.
            HR=ZETA*H*EPS**(-0.2D0)
         ELSE
            HR=ZETA*H**2/HP*(EPSP/EPS**2)**(0.2D0)
         ENDIF
         DO 60 KN=1,NEQN
            Y(KN)=YS(KN,S)
            DY(KN)=DYS(KN,S)
   60    CONTINUE
         DO 80 KS=1,S
            DO 70 KN=1,NEQN
               DYSP(KN,KS)=DYS(KN,KS)
   70       CONTINUE
   80    CONTINUE
         T=T+H
         IF (ABS(TEND-T).LT.10D0*UROUND*ABS(T)) T=TEND
         HP=H
         EPSP=EPS
         SUCREJ=.FALSE.
         JACU2D=.FALSE.
      ELSE
         IF (.NOT.FIRST.AND.SUCREJ) THEN
            PEST=MIN(5D0,MAX(PMIN,LOG10(EPS/EPSREJ)/LOG10(H/HREJ)))
            HR=ZETA*H*EPS**(-1D0/PEST)
         ELSE
            HR=ZETA*H*EPS**(-0.2D0)
         ENDIF
         HREJ=H
         EPSREJ=EPS
         SUCREJ=.TRUE.
         NREJE=NREJE+1
      ENDIF
      RETURN
      END
C$Id: snorm.f,v 1.10 1998/11/25 08:58:31 walter Exp $
      DOUBLE PRECISION FUNCTION SNORM(NEQN,Y,YNORM,H,TOLVEC,RTOL,ATOL,
     +                                INDGT1,IND)
CF90  IMPLICIT NONE
      INTEGER NEQN,IND(*)
      DOUBLE PRECISION Y(NEQN),YNORM(NEQN),H,RTOL(*),ATOL(*)
      LOGICAL TOLVEC,INDGT1
CF90  INTENT(IN)    NEQN,Y,YNORM,H,TOLVEC,RTOL,ATOL,INDGT1,IND
C-----------------------------------------------------------------------
C     function that computes the scaled norm of a vector
C     of dimension NEQN
C-----------------------------------------------------------------------
      INTEGER KN
      SNORM=0D0
      IF (INDGT1) THEN
         IF (TOLVEC) THEN
            DO 10 KN=1,NEQN
                  SNORM=SNORM+(H**(IND(KN)-1)*
     +                  YNORM(KN)/(ATOL(KN)+RTOL(KN)*ABS(Y(KN))))**2
   10       CONTINUE
         ELSE
            DO 20 KN=1,NEQN
               SNORM=SNORM+(H**(IND(KN)-1)*
     +               YNORM(KN)/(ATOL(1)+RTOL(1)*ABS(Y(KN))))**2
   20       CONTINUE
         ENDIF
      ELSE
         IF (TOLVEC) THEN
            DO 30 KN=1,NEQN
               SNORM=SNORM+
     +               (YNORM(KN)/(ATOL(KN)+RTOL(KN)*ABS(Y(KN))))**2
   30       CONTINUE
         ELSE
            DO 40 KN=1,NEQN
               SNORM=SNORM+
     +               (YNORM(KN)/(ATOL(1)+RTOL(1)*ABS(Y(KN))))**2
   40       CONTINUE
         ENDIF
      ENDIF
      SNORM=SQRT(SNORM/NEQN)
      RETURN
      END
C$Id: snorms.f,v 1.10 1998/11/25 08:58:32 walter Exp $
      DOUBLE PRECISION FUNCTION SNORMS(NEQN,Y,YSNORM,H,TOLVEC,RTOL,ATOL,
     +                                INDGT1,IND)
CF90  IMPLICIT NONE
C begin included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER S
      PARAMETER (S = 4)
C end included $Id: s.h,v 1.3 1997/11/25 17:11:51 walter Rel $
      INTEGER NEQN,IND(*)
      DOUBLE PRECISION Y(NEQN),YSNORM(NEQN,S),H,RTOL(*),ATOL(*)
      LOGICAL TOLVEC,INDGT1
CF90  INTENT(IN)    NEQN,Y,YSNORM,H,TOLVEC,RTOL,ATOL,INDGT1,IND
C-----------------------------------------------------------------------
C     function that computes the scaled norm of a vector
C     of dimension NEQN x S
C-----------------------------------------------------------------------
      INTEGER KS,KN
      SNORMS=0D0
      IF (INDGT1) THEN
         IF (TOLVEC) THEN
            DO 20 KS=1,S
               DO 10 KN=1,NEQN
                  SNORMS=SNORMS+(H**(IND(KN)-1)*
     +                 YSNORM(KN,KS)/(ATOL(KN)+RTOL(KN)*ABS(Y(KN))))**2
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 KS=1,S
               DO 30 KN=1,NEQN
                  SNORMS=SNORMS+(H**(IND(KN)-1)*
     +                   YSNORM(KN,KS)/(ATOL(1)+RTOL(1)*ABS(Y(KN))))**2
   30          CONTINUE
   40       CONTINUE
         ENDIF
      ELSE
         IF (TOLVEC) THEN
            DO 60 KS=1,S
               DO 50 KN=1,NEQN
                  SNORMS=SNORMS+
     +                (YSNORM(KN,KS)/(ATOL(KN)+RTOL(KN)*ABS(Y(KN))))**2
   50          CONTINUE
   60       CONTINUE
         ELSE
            DO 80 KS=1,S
               DO 70 KN=1,NEQN
                  SNORMS=SNORMS+
     +                  (YSNORM(KN,KS)/(ATOL(1)+RTOL(1)*ABS(Y(KN))))**2
   70          CONTINUE
   80       CONTINUE
         ENDIF
      ENDIF
      SNORMS=SQRT(SNORMS/(NEQN*S))
      RETURN
      END
