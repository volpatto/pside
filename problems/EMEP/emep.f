c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        EMEP problem
c        ODE of dimension 66
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/emep.f
c
c     This is revision
c     $Id: emep.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c-----------------------------------------------------------------------
      integer function pidate()
      pidate = 20060828
      return 
      end
c-----------------------------------------------------------------------
      subroutine prob(fullnm,problm,type,
     +                neqn,ndisc,t,
     +                numjac,mljac,mujac,
     +                nummas,mlmas,mumas,
     +                ind)
      character*(*) fullnm, problm, type
      integer neqn,ndisc,mljac,mujac,mlmas,mumas,ind(*)
      double precision t(0:*)
      logical numjac, nummas
C
C=======================================================================
C
C     EMEP MSC-W OZONE MODEL CHEMISTRY
C
C=======================================================================
C
      INTEGER I
C Parameters
      INTEGER NSPEC
      PARAMETER (NSPEC=66)
C     NSPEC : 66    number of species

      FULLNM = 'EMEP problem'
      PROBLM = 'emep'
      TYPE   = 'ODE'
      NEQN   = NSPEC
      NDISC  = 8
      T(0)   =  4D0*3600D0
      T(1)   = 20D0*3600D0
      DO 10 I=1,4
         T(2*I)   = T(0) + 24D0*3600D0*DBLE(I)
         T(2*I+1) = T(1) + 24D0*3600D0*DBLE(I)
   10 CONTINUE
      NUMJAC = .FALSE.
      MLJAC  = NEQN
      MUJAC  = NEQN

      RETURN
      END
c-----------------------------------------------------------------------
      subroutine init(neqn,t,y,yprime,consis)
      integer neqn
      double precision t,y(neqn),yprime(neqn)
      logical consis

      INTEGER I
C
C
C Establishment of initial conditions:
C    completely arbitrary at this stage !
      DO 10 I = 1, 13
          Y(I)=1.D7
 10   CONTINUE
C
      DO 20 I = 14, NEQN
         Y(I)=100.D0
 20   CONTINUE
      Y(1) = 1.0D9
      Y(2) = 5.0D9
      Y(3) = 5.0D9
      Y(4) =3.8D12
      Y(5) =3.5D13
      Y(14)=5.D11
      Y(38)=1.D-3
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine settolerances(neqn,rtol,atol,tolvec)
      integer neqn 
      logical tolvec
      double precision rtol(neqn), atol(neqn)
       
      tolvec  = .false.
      

      return
      end
c-----------------------------------------------------------------------
      subroutine setoutput(neqn,solref,printsolout,
     +                    nindsol,indsol)

      logical solref, printsolout
      integer neqn, nindsol
      integer indsol(neqn)

c the reference solution is available
      solref = .true.  

c output file is required
      printsolout = .true.

c default values if printsolout is .true.
      nindsol = 6
c only nindsol component of indsol are referenced
      indsol(1) =  1
      indsol(2) =  2
      indsol(3) =  3
      indsol(4) =  5
      indsol(5) = 14
      indsol(6) = 40

      return
      end


c-----------------------------------------------------------------------
      subroutine feval(neqn,time,y,yprime,dy,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision time,y(neqn),yprime(neqn),dy(neqn),rpar(*)

      DOUBLE PRECISION M, O2, XN2, RPATH3, RPATH4
      DOUBLE PRECISION S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11
C
C=======================================================================
C
C     EMEP MSC-W OZONE MODEL CHEMISTRY
C
C=======================================================================
C
C Parameters
      INTEGER NSPEC, NRC, NDJ
      DOUBLE PRECISION HMIX
      PARAMETER (NSPEC=66, NRC=266, NDJ=16)
      PARAMETER (HMIX=1.2D5)
C     NSPEC : 66    number of species
C     NRC   : 266   size of rate constant array
C     NDJ   : 16    number of photolysis reactions
C     HMIX  : mixing height in cm
C     EMIS1..EMIS13 : emitted species
C
C Listing of species
C     Y(1:NSPEC) =
C        NO,     NO2,    SO2,    CO,     CH4,     C2H6,
C        NC4H10, C2H4,   C3H6,   OXYL,   HCHO,    CH3CHO,
C        MEK,    O3,     HO2,    HNO3,   H2O2,    H2,
C        CH3O2,  C2H5OH, SA,     CH3O2H, C2H5O2,  CH3COO,
C        PAN,    SECC4H, MEKO2,  R2OOH,  ETRO2,   MGLYOX,
C        PRRO2,  GLYOX,  OXYO2,  MAL,    MALO2,   OP,
C        OH,     OD,     NO3,    N2O5,   ISOPRE,  NITRAT,
C        ISRO2,  MVK,    MVKO2,  CH3OH,  RCO3H,   OXYO2H,
C        BURO2H, ETRO2H, PRRO2H, MEKO2H, MALO2H,  MACR,
C        ISNI,   ISRO2H, MARO2,  MAPAN,  CH2CCH3, ISONO3,
C        ISNIR,  MVKO2H, CH2CHR, ISNO3H, ISNIRH,  MARO2H
C
C=======================================================================
C
C Emissions in molec/(cm**2*s), of NO, NO2, SO2, CO, CH4, C2H6, NC4H10,
C    C2H4, C3H6, O-XYLENE, C2H5OH and ISOPRENE
      DOUBLE PRECISION EMIS1, EMIS2, EMIS3, EMIS4, EMIS5, EMIS6,
     +   EMIS7, EMIS8, EMIS9, EMIS10, EMIS11, EMIS12, EMIS13
      DOUBLE PRECISION FRAC6, FRAC7, FRAC8, FRAC9, FRAC10, FRAC13
      DOUBLE PRECISION VEKT6, VEKT7, VEKT8, VEKT9, VEKT10, VEKT13
      DOUBLE PRECISION VMHC, EMNOX, EMHC, FACISO, FACHC, FACNOX
C
C   distribution of VOC emissions among species:
      PARAMETER (FRAC6=0.07689D0, FRAC7=0.41444D0, FRAC8=0.03642D0,
     +           FRAC9=0.03827D0, FRAC10=0.24537D0, FRAC13=0.13957D0)
      PARAMETER (VEKT6=30.D0, VEKT7=58.D0, VEKT8=28.D0,
     +           VEKT9=42.D0, VEKT10=106.D0, VEKT13=46.D0)
      PARAMETER (VMHC=1.D0/(FRAC6/VEKT6+FRAC7/VEKT7+FRAC8/VEKT8+
     +           FRAC9/VEKT9+FRAC10/VEKT10+FRAC13/VEKT13))
C
C   choose values for NOX and HC emissions in molecules/cm**2xs
      PARAMETER (EMNOX=2.5D11, EMHC=2.5D11)
C
C   rural case
      PARAMETER ( FACISO=1.0D0, FACHC=1.0D0, FACNOX=1.0D0)
C
      PARAMETER (EMIS1 = EMNOX * FACNOX,
     +           EMIS2 = 0.D0,
     +           EMIS3 = EMNOX * FACNOX,
     +           EMIS4 = EMHC*10.D0 * FACHC,
     +           EMIS5 = 0.D0,
     +           EMIS6 = EMHC*FRAC6/VEKT6*VMHC  * FACHC,
     +           EMIS7 = EMHC*FRAC7/VEKT7*VMHC  * FACHC,
     +           EMIS8 = EMHC*FRAC8/VEKT8*VMHC  * FACHC,
     +           EMIS9 = EMHC*FRAC9/VEKT9*VMHC  * FACHC,
     +           EMIS10= EMHC*FRAC10/VEKT10*VMHC* FACHC,
     +           EMIS11= 0.5D0 * FACISO * EMHC,
     +           EMIS12= 0.D0,
     +           EMIS13= EMHC*FRAC13/VEKT13*VMHC* FACHC)
C
C
      DOUBLE PRECISION RC(NRC), DJ(NDJ), H2O
C
C=======================================================================
C
C Compute time-dependent EMEP coefficients
      CALL EMEPCF (TIME, RC, DJ, H2O)
C
C
      M = 2.55D19
      O2 = 5.2D18
      XN2= 1.99D19
C..  pathways for decay of secc4h9o:
C..  Newer assumption, from Atkisnon , 1991
      RPATH3 = 0.65D0
      RPATH4 = 0.35D0
C=======================================================================
C
C
      DY(1) = DJ(3)*Y(2)+DJ(13)*Y(39)+RC(19)*Y(2)*Y(39)+EMIS1/HMIX-(RC
     &(5)*Y(36)+RC(11)*Y(14)+RC(17)*Y(15)+RC(72)*Y(23)+RC(79)*(Y(24)+Y(5
     &7))+RC(15)*Y(39)+RC(60)*(Y(19)+Y(26)+Y(27)+Y(29)+Y(31)+Y(33)+Y(35)
     &+Y(43)+Y(45)+Y(59)+Y(61)+Y(60)))*Y(1)
      DY(2) = Y(1)*(RC(5)*Y(36)+RC(11)*Y(14)+RC(17)*Y(15)+RC(72)*Y(23)+R
     &C(79)*(Y(24)+Y(57))+0.2D1*RC(15)*Y(39))+RC(60)*Y(1)*(Y(19)+Y(26)+Y
     &(27)+Y(29)+Y(31)+Y(33)+Y(35)+Y(59))+RC(60)*Y(1)*(0.86D0*Y(43)+0.19
     &D1*Y(61)+0.11D1*Y(60)+0.95D0*Y(45))+DJ(14)*Y(39)+DJ(5)*Y(16)+DJ(15
     &)*Y(40)+RC(29)*Y(40)+RC(78)*(Y(25)+Y(58))-(DJ(3)+RC(12)*Y(14)+RC(2
     &0)*Y(39)+RC(21)*Y(37)+RC(48)+RC(77)*(Y(24)+Y(57)))*Y(2)
      DY(3) = EMIS3/HMIX-(RC(39)*Y(37)+RC(40)*Y(19)+RC(47))*Y(3)
      DY(4) = EMIS4/HMIX+Y(37)*(RC(66)*Y(11)+2.D0*RC(221)*Y(32)+RC(222
     &)*Y(30))+Y(14)*(0.44D0*RC(112)*Y(8)+0.4D0*RC(123)*Y(9)+0.5D-1*RC(1
     &60)*Y(44)+0.5D-1*RC(150)*Y(41))+Y(11)*(DJ(6)+DJ(7)+RC(69)*Y(39))+D
     &J(8)*Y(12)+DJ(11)*Y(30)+2.D0*DJ(7)*Y(32)-RC(70)*Y(37)*Y(4)
      DY(5) = EMIS5/HMIX+0.7D-1*RC(123)*Y(14)*Y(9)-RC(59)*Y(37)*Y(5)
      DY(6) = EMIS6/HMIX-RC(71)*Y(37)*Y(6)
      DY(7) = EMIS7/HMIX-RC(81)*Y(37)*Y(7)
      DY(8) = EMIS8/HMIX-(RC(109)*Y(37)+RC(112)*Y(14))*Y(8)
      DY(9) = EMIS9/HMIX+0.7D-1*RC(150)*Y(14)*Y(41)-(RC(123)*Y(14)+RC(
     &125)*Y(37))*Y(9)
      DY(10) = EMIS10/HMIX-RC(234)*Y(37)*Y(10)
      DY(11) = Y(19)*(RC(60)*Y(1)+(2.D0*RC(61)+RC(62))*Y(19)+RC(80)*Y(24
     &)+RC(40)*Y(3))+Y(37)*(RC(63)*Y(46)+RC(67)*Y(22))+Y(1)*RC(60)*(2.D0
     &*Y(29)+Y(31)+0.74D0*Y(43)+0.266D0*Y(45)+0.15D0*Y(60))+Y(14)*(0.5D0
     &*RC(123)*Y(9)+RC(112)*Y(8)+0.7D0*RC(157)*Y(56)+0.8D0*RC(160)*Y(44)
     &+0.8D0*RC(150)*Y(41))+2.D0*DJ(7)*Y(32)+DJ(16)*(Y(22)+0.156D1*Y(50)
     &+Y(51))-(RC(66)*Y(37)+DJ(6)+DJ(7)+RC(69)*Y(39)+RC(53))*Y(11)
      DY(12) = Y(1)*(RC(72)*Y(23)+RC(83)*Y(26)*RPATH4+RC(105)*Y(27)+RC(1
     &26)*Y(31)+0.95D0*RC(162)*Y(61)+0.684D0*RC(154)*Y(45))+Y(37)*(RC(64
     &)*Y(20)+RC(76)*Y(28)+RC(76)*Y(50))+0.5D0*RC(123)*Y(14)*Y(9)+0.4D-1
     &*RC(160)*Y(14)*Y(44)+DJ(16)*(Y(28)+0.22D0*Y(50)+0.35D0*Y(49)+Y(51)
     &+Y(52))-(DJ(8)+RC(75)*Y(37)+RC(53))*Y(12)
      DY(13) = RC(83)*Y(1)*Y(26)*RPATH3+(0.65D0*DJ(16)+RC(76)*Y(37))*Y(4
     &9)+RC(76)*Y(37)*Y(51)+(RC(159)*Y(1)+DJ(16))*Y(59)+0.95D0*RC(162)*Y
     &(1)*Y(61)-(DJ(9)+RC(86)*Y(37)+RC(53))*Y(13)
      DY(14) = RC(1)*Y(36)+RC(89)*Y(15)*Y(24)-(RC(11)*Y(1)+RC(12)*Y(2)+R
     &C(13)*Y(37)+RC(14)*Y(15)+RC(49)+RC(112)*Y(8)+RC(123)*Y(9)+RC(157)*
     &Y(56)+RC(160)*Y(44)+RC(150)*Y(41)+DJ(1)+DJ(2))*Y(14)
      s1 = Y(37)*(RC(13)*Y(14)+RC(31)*Y(17)+RC(33)*Y(18)+RC(39)*Y(3)+RC(
     &63)*Y(46)+RC(64)*Y(20)+RC(66)*Y(11)+RC(70)*Y(4)+RC(221)*Y(32))+Y(1
     &9)*(RC(40)*Y(3)+2.D0*RC(61)*Y(19)+0.5D0*RC(80)*Y(24))+DJ(11)*Y(30)
     &+Y(1)*RC(60)*(Y(19)+Y(29)+Y(31)+Y(33)+Y(35)+0.95D0*Y(45)+Y(26)*RPA
     &TH3+0.78D0*Y(43)+Y(59)+0.5D-1*Y(61)+0.8D0*Y(60))+RC(72)*Y(1)*Y(23)
     &+DJ(8)*Y(12)
      DY(15) = s1+2.D0*DJ(6)*Y(11)+DJ(16)*(Y(22)+Y(28)+0.65D0*Y(49)+Y(50
     &)+Y(51)+Y(48)+Y(53))+Y(39)*(RC(26)*Y(17)+RC(69)*Y(11))+Y(14)*(0.12
     &D0*RC(112)*Y(8)+0.28D0*RC(123)*Y(9)+0.6D-1*RC(160)*Y(44))+0.6D-1*R
     &C(150)*Y(14)*Y(41)-(RC(14)*Y(14)+RC(17)*Y(1)+RC(30)*Y(37)+2.D0*RC(
     &36)*Y(15)+RC(65)*Y(19)+RC(74)*Y(23)+(RC(88)+RC(89))*Y(24)+RC(85)*(
     &Y(26)+Y(29)+Y(31)+Y(27)+Y(57)+Y(45)+Y(61)+Y(59)+Y(33)+Y(35)+Y(43)+
     &Y(60)))*Y(15)
      DY(16) = RC(21)*Y(2)*Y(37)+Y(39)*(RC(26)*Y(17)+RC(69)*Y(11))-(RC(3
     &5)*Y(37)+DJ(5)+RC(45))*Y(16)
      DY(17) = RC(36)*Y(15)**2-(RC(31)*Y(37)+DJ(4)+RC(43)+RC(26)*Y(39)+R
     &C(47))*Y(17)
      DY(18) = DJ(7)*Y(11)+Y(14)*(0.13D0*RC(112)*Y(8)+0.7D-1*RC(123)*Y(9
     &))-RC(33)*Y(37)*Y(18)
      DY(19) = Y(37)*(RC(59)*Y(5)+RC(68)*Y(22))+Y(24)*(RC(79)*Y(1)+2.D0*
     &RC(94)*Y(24))+DJ(8)*Y(12)+DJ(16)*Y(47)+0.31D0*RC(123)*Y(14)*Y(9)-(
     &RC(40)*Y(3)+RC(60)*Y(1)+2.D0*RC(61)*Y(19)+2.D0*RC(62)*Y(19)+RC(65)
     &*Y(15)+0.5D0*RC(80)*Y(24))*Y(19)
      DY(20) = EMIS13/HMIX-RC(64)*Y(37)*Y(20)
      DY(21) = (RC(40)*Y(19)+RC(39)*Y(37))*Y(3)+0.5D-1*EMIS3/HMIX-RC(5
     &1)*Y(21)
      DY(22) = RC(65)*Y(15)*Y(19)-(RC(43)+DJ(16)+(RC(67)+RC(68))*Y(37))*
     &Y(22)
      DY(23) = Y(37)*(RC(71)*Y(6)+RC(68)*Y(28))+0.35D0*DJ(16)*Y(49)+RC(8
     &3)*Y(1)*Y(26)*RPATH4+DJ(9)*Y(13)-(RC(72)*Y(1)+RC(74)*Y(15))*Y(23)
      DY(24) = Y(37)*(RC(75)*Y(12)+RC(222)*Y(30)+RC(68)*Y(47))+RC(105)*Y
     &(1)*Y(27)+RC(78)*Y(25)+DJ(11)*Y(30)+DJ(9)*Y(13)+DJ(16)*Y(52)+0.684
     &D0*RC(154)*Y(1)*Y(45)-(RC(77)*Y(2)+RC(79)*Y(1)+RC(80)*Y(19)+2.D0*R
     &C(94)*Y(24)+(RC(88)+RC(89))*Y(15))*Y(24)
      DY(25) = RC(77)*Y(24)*Y(2)-(RC(50)+RC(78))*Y(25)
      DY(26) = Y(37)*(RC(81)*Y(7)+RC(68)*Y(49))-(RC(83)*Y(1)+RC(85)*Y(15
     &))*Y(26)
      DY(27) = Y(37)*(RC(86)*Y(13)+RC(87)*Y(52))-(RC(105)*Y(1)+RC(85)*Y(
     &15))*Y(27)
      DY(28) = RC(74)*Y(15)*Y(23)-((RC(76)+RC(68))*Y(37)+DJ(16)+RC(52))*
     &Y(28)
      DY(29) = Y(37)*(RC(109)*Y(8)+RC(68)*Y(50))-(RC(110)*Y(1)+RC(85)*Y(
     &15))*Y(29)
      DY(30) = RC(236)*Y(1)*Y(33)+RC(220)*Y(1)*Y(35)+0.266D0*RC(154)*Y(1
     &)*Y(45)+0.82D0*RC(160)*Y(14)*Y(44)+DJ(16)*(Y(48)+Y(53))-(DJ(11)+RC
     &(222)*Y(37))*Y(30)
      DY(31) = Y(37)*(RC(125)*Y(9)+RC(68)*Y(51))-(RC(126)*Y(1)+RC(85)*Y(
     &15))*Y(31)
      DY(32) = RC(220)*Y(1)*Y(35)+DJ(16)*Y(53)-(2.D0*DJ(7)+RC(221)*Y(37)
     &)*Y(32)
      DY(33) = Y(37)*(RC(234)*Y(10)+RC(235)*Y(48))-(RC(236)*Y(1)+RC(85)*
     &Y(15))*Y(33)
      DY(34) = RC(236)*Y(1)*Y(33)+DJ(16)*Y(48)-RC(219)*Y(37)*Y(34)
      DY(35) = Y(37)*(RC(219)*Y(34)+RC(223)*Y(53))-(RC(220)*Y(1)+RC(85)*
     &Y(15))*Y(35)
      DY(36) = DJ(1)*Y(14)+DJ(3)*Y(2)+DJ(14)*Y(39)+RC(7)*Y(38)+0.2D0*RC(
     &160)*Y(14)*Y(44)+0.3D0*RC(150)*Y(14)*Y(41)-(RC(1)+RC(5)*Y(1))*Y(36
     &)
      s1 = 2.D0*RC(8)*H2O*Y(38)+Y(15)*(RC(14)*Y(14)+RC(17)*Y(1))+2.D0*DJ
     &(4)*Y(17)+DJ(5)*Y(16)
      s2 = s1+DJ(16)*(Y(22)+Y(28)+Y(47)+Y(49)+Y(50)+Y(52)+Y(48)+Y(53))
      s3 = s2+Y(14)*(0.15D0*RC(123)*Y(9)+0.8D-1*RC(160)*Y(44))
      s4 = s3
      s6 = 0.55D0*RC(150)*Y(14)*Y(41)
      s8 = -1
      s11 = RC(222)*Y(30)+RC(75)*Y(12)+RC(81)*Y(7)+RC(87)*Y(52)+RC(86)*Y
     &(13)+RC(235)*Y(48)+RC(109)*Y(8)+RC(125)*Y(9)+RC(234)*Y(10)+RC(223)
     &*Y(53)+RC(219)*Y(34)+RC(31)*Y(17)+RC(21)*Y(2)+RC(148)*Y(62)+RC(64)
     &*Y(20)+RC(39)*Y(3)+RC(71)*Y(6)
      s10 = s11+RC(30)*Y(15)+RC(59)*Y(5)+RC(70)*Y(4)+RC(35)*Y(16)+RC(13)
     &*Y(14)+RC(221)*Y(32)+RC(68)*(Y(22)+Y(28)+Y(47)+Y(50)+Y(51)+Y(49))+
     &RC(66)*Y(11)+RC(151)*Y(41)+RC(153)*Y(44)+RC(63)*Y(46)+RC(33)*Y(18)
     &+RC(158)*Y(54)+RC(146)*Y(63)+RC(149)*(Y(65)+Y(66))+RC(147)*Y(64)+R
     &C(161)*Y(55)
      s11 = Y(37)
      s9 = s10*s11
      s7 = s8*s9
      s5 = s6+s7
      DY(37) = s4+s5
      DY(38) = DJ(2)*Y(14)-(RC(7)+RC(8)*H2O)*Y(38)
      DY(39) = (RC(29)+DJ(15))*Y(40)+RC(12)*Y(14)*Y(2)+RC(35)*Y(37)*Y(16
     &)-(RC(15)*Y(1)+RC(26)*Y(17)+RC(163)*Y(41)+RC(19)*Y(2)+RC(20)*Y(2)+
     &DJ(13)+DJ(14)+RC(69)*Y(11))*Y(39)
      DY(40) = RC(20)*Y(39)*Y(2)-(RC(29)+DJ(15)+RC(45))*Y(40)
      DY(41) = EMIS11/HMIX-(RC(151)*Y(37)+RC(163)*Y(39)+RC(150)*Y(14))
     &*Y(41)
      DY(42) = RC(45)*Y(16)+2.D0*RC(44)*Y(40)-RC(51)*Y(42)
      DY(43) = Y(37)*(RC(151)*Y(41)+RC(156)*Y(56))+0.12D0*RC(152)*Y(1)*Y
     &(43)-(RC(152)*Y(1)+RC(155)*Y(15))*Y(43)
      DY(44) = RC(60)*Y(1)*(0.42D0*Y(43)+0.5D-1*Y(60))+0.26D0*RC(150)*Y(
     &14)*Y(41)-(RC(153)*Y(37)+RC(160)*Y(14))*Y(44)
      DY(45) = RC(153)*Y(44)*Y(37)+RC(148)*Y(37)*Y(62)-(RC(154)*Y(1)+RC(
     &85)*Y(15))*Y(45)
      DY(46) = RC(62)*Y(19)**2-RC(63)*Y(37)*Y(46)
      DY(47) = RC(88)*Y(15)*Y(24)-(RC(68)*Y(37)+DJ(16)+RC(52))*Y(47)
      DY(48) = RC(85)*Y(15)*Y(33)-(RC(235)*Y(37)+DJ(16)+RC(52))*Y(48)
      DY(49) = RC(85)*Y(15)*Y(26)-((RC(76)+RC(68))*Y(37)+DJ(16)+RC(52))*
     &Y(49)
      DY(50) = RC(85)*Y(15)*Y(29)-((RC(76)+RC(68))*Y(37)+DJ(16)+RC(52))*
     &Y(50)
      DY(51) = RC(85)*Y(15)*Y(31)-((RC(76)+RC(68))*Y(37)+DJ(16)+RC(52))*
     &Y(51)
      DY(52) = RC(85)*Y(15)*Y(27)-(RC(87)*Y(37)+DJ(16)+RC(52))*Y(52)
      DY(53) = RC(85)*Y(15)*Y(35)-(RC(223)*Y(37)+DJ(16)+RC(52))*Y(53)
      DY(54) = RC(60)*Y(1)*(0.32D0*Y(43)+0.1D0*Y(60))+0.67D0*RC(150)*Y(1
     &4)*Y(41)-RC(158)*Y(37)*Y(54)
      DY(55) = RC(60)*Y(1)*(0.14D0*Y(43)+0.5D-1*Y(45)+0.85D0*Y(60))-RC(1
     &61)*Y(37)*Y(55)
      DY(56) = RC(155)*Y(15)*Y(43)-(RC(156)*Y(37)+RC(157)*Y(14)+RC(52))*
     &Y(56)
      DY(57) = 0.5D0*RC(158)*Y(37)*Y(54)+RC(78)*Y(58)+RC(149)*Y(37)*Y(66
     &)-(RC(77)*Y(2)+RC(79)*Y(1)+RC(85)*Y(15))*Y(57)
      DY(58) = RC(77)*Y(57)*Y(2)-(RC(50)+RC(78))*Y(58)
      DY(59) = RC(79)*Y(1)*Y(57)+RC(146)*Y(37)*Y(63)-(RC(159)*Y(1)+RC(85
     &)*Y(15))*Y(59)
      DY(60) = RC(163)*Y(39)*Y(41)+RC(147)*Y(37)*Y(64)-(RC(164)*Y(1)+RC(
     &85)*Y(15))*Y(60)
      DY(61) = RC(161)*Y(37)*Y(55)+RC(149)*Y(37)*Y(65)-(RC(162)*Y(1)+RC(
     &85)*Y(15))*Y(61)
      DY(62) = RC(85)*Y(15)*Y(45)-(RC(148)*Y(37)+RC(52))*Y(62)
      DY(63) = RC(85)*Y(15)*Y(59)-(RC(146)*Y(37)+RC(52))*Y(63)
      DY(64) = RC(85)*Y(15)*Y(60)-(RC(147)*Y(37)+RC(52))*Y(64)
      DY(65) = RC(85)*Y(15)*Y(61)-(RC(149)*Y(37)+RC(52))*Y(65)
      DY(66) = RC(85)*Y(15)*Y(57)-(RC(149)*Y(37)+RC(52))*Y(66)
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine jeval(ldim,neqn,time,y,yprime,jac,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision time,y(neqn),yprime(neqn),jac(ldim,neqn),rpar(*)

      INTEGER I, J
      DOUBLE PRECISION M, O2, XN2, RPATH3, RPATH4
      DOUBLE PRECISION S1
C
C Parameters
      INTEGER NSPEC, NRC, NDJ
      DOUBLE PRECISION HMIX
      PARAMETER (NSPEC=66, NRC=266, NDJ=16)
      PARAMETER (HMIX=1.2D5)
C
      DOUBLE PRECISION EMIS1, EMIS2, EMIS3, EMIS4, EMIS5, EMIS6,
     +   EMIS7, EMIS8, EMIS9, EMIS10, EMIS11, EMIS12, EMIS13
      DOUBLE PRECISION FRAC6, FRAC7, FRAC8, FRAC9, FRAC10, FRAC13
      DOUBLE PRECISION VEKT6, VEKT7, VEKT8, VEKT9, VEKT10, VEKT13
      DOUBLE PRECISION VMHC, EMNOX, EMHC, FACISO, FACHC, FACNOX
C
C   distribution of VOC emissions among species:
      PARAMETER (FRAC6=0.07689D0, FRAC7=0.41444D0, FRAC8=0.03642D0,
     +           FRAC9=0.03827D0, FRAC10=0.24537D0, FRAC13=0.13957D0)
      PARAMETER (VEKT6=30.D0, VEKT7=58.D0, VEKT8=28.D0,
     +           VEKT9=42.D0, VEKT10=106.D0, VEKT13=46.D0)
      PARAMETER (VMHC=1.D0/(FRAC6/VEKT6+FRAC7/VEKT7+FRAC8/VEKT8+
     +           FRAC9/VEKT9+FRAC10/VEKT10+FRAC13/VEKT13))
C
C   choose values for NOX and HC emissions in molecules/cm**2xs
      PARAMETER (EMNOX=2.5D11, EMHC=2.5D11)
C
C   rural case
      PARAMETER ( FACISO=1.0D0, FACHC=1.0D0, FACNOX=1.0D0)
C
      PARAMETER (EMIS1 = EMNOX * FACNOX,
     +           EMIS2 = 0.D0,
     +           EMIS3 = EMNOX * FACNOX,
     +           EMIS4 = EMHC*10.D0 * FACHC,
     +           EMIS5 = 0.D0,
     +           EMIS6 = EMHC*FRAC6/VEKT6*VMHC  * FACHC,
     +           EMIS7 = EMHC*FRAC7/VEKT7*VMHC  * FACHC,
     +           EMIS8 = EMHC*FRAC8/VEKT8*VMHC  * FACHC,
     +           EMIS9 = EMHC*FRAC9/VEKT9*VMHC  * FACHC,
     +           EMIS10= EMHC*FRAC10/VEKT10*VMHC* FACHC,
     +           EMIS11= 0.5D0 * FACISO * EMHC,
     +           EMIS12= 0.D0,
     +           EMIS13= EMHC*FRAC13/VEKT13*VMHC* FACHC)
C
C
      DOUBLE PRECISION RC(NRC), DJ(NDJ), H2O
C
C=======================================================================
C
C Compute time-dependent EMEP coefficients
      CALL EMEPCF (TIME, RC, DJ, H2O)
C
C
      M = 2.55D19
      O2 = 5.2D18
      XN2= 1.99D19
C..  pathways for decay of secc4h9o:
C..  Newer assumption, from Atkisnon , 1991
      RPATH3 = 0.65D0
      RPATH4 = 0.35D0
C=======================================================================
C
      DO 10 I = 1, 66
      DO 10 J = 1, 66
         JAC(I,J) = 0.0D0
   10 CONTINUE
      JAC(1,1) = -RC(5)*Y(36)-RC(11)*Y(14)-RC(17)*Y(15)-RC(72)*Y(23)-RC(
     &79)*(Y(24)+Y(57))-RC(15)*Y(39)-RC(60)*(Y(19)+Y(26)+Y(27)+Y(29)+Y(3
     &1)+Y(33)+Y(35)+Y(43)+Y(45)+Y(59)+Y(61)+Y(60))
      JAC(1,2) = DJ(3)+RC(19)*Y(39)
      JAC(1,14) = -RC(11)*Y(1)
      JAC(1,15) = -RC(17)*Y(1)
      JAC(1,19) = -RC(60)*Y(1)
      JAC(1,23) = -RC(72)*Y(1)
      JAC(1,24) = -RC(79)*Y(1)
      JAC(1,26) = -RC(60)*Y(1)
      JAC(1,27) = -RC(60)*Y(1)
      JAC(1,29) = -RC(60)*Y(1)
      JAC(1,31) = -RC(60)*Y(1)
      JAC(1,33) = -RC(60)*Y(1)
      JAC(1,35) = -RC(60)*Y(1)
      JAC(1,36) = -RC(5)*Y(1)
      JAC(1,39) = DJ(13)+RC(19)*Y(2)-RC(15)*Y(1)
      JAC(1,43) = -RC(60)*Y(1)
      JAC(1,45) = -RC(60)*Y(1)
      JAC(1,57) = -RC(79)*Y(1)
      JAC(1,59) = -RC(60)*Y(1)
      JAC(1,60) = -RC(60)*Y(1)
      JAC(1,61) = -RC(60)*Y(1)
      JAC(2,1) = RC(5)*Y(36)+RC(11)*Y(14)+RC(17)*Y(15)+RC(72)*Y(23)+RC(7
     &9)*(Y(24)+Y(57))+0.2D1*RC(15)*Y(39)+RC(60)*(Y(19)+Y(26)+Y(27)+Y(29
     &)+Y(31)+Y(33)+Y(35)+Y(59))+RC(60)*(0.86D0*Y(43)+0.19D1*Y(61)+0.11D
     &1*Y(60)+0.95D0*Y(45))
      JAC(2,2) = -DJ(3)-RC(12)*Y(14)-RC(20)*Y(39)-RC(21)*Y(37)-RC(48)-RC
     &(77)*(Y(24)+Y(57))
      JAC(2,14) = RC(11)*Y(1)-RC(12)*Y(2)
      JAC(2,15) = RC(17)*Y(1)
      JAC(2,16) = DJ(5)
      JAC(2,19) = RC(60)*Y(1)
      JAC(2,23) = RC(72)*Y(1)
      JAC(2,24) = RC(79)*Y(1)-RC(77)*Y(2)
      JAC(2,25) = RC(78)
      JAC(2,26) = RC(60)*Y(1)
      JAC(2,27) = RC(60)*Y(1)
      JAC(2,29) = RC(60)*Y(1)
      JAC(2,31) = RC(60)*Y(1)
      JAC(2,33) = RC(60)*Y(1)
      JAC(2,35) = RC(60)*Y(1)
      JAC(2,36) = RC(5)*Y(1)
      JAC(2,37) = -RC(21)*Y(2)
      JAC(2,39) = 0.2D1*RC(15)*Y(1)+DJ(14)-RC(20)*Y(2)
      JAC(2,40) = RC(29)+DJ(15)
      JAC(2,43) = 0.86D0*RC(60)*Y(1)
      JAC(2,45) = 0.95D0*RC(60)*Y(1)
      JAC(2,57) = RC(79)*Y(1)-RC(77)*Y(2)
      JAC(2,58) = RC(78)
      JAC(2,59) = RC(60)*Y(1)
      JAC(2,60) = 0.11D1*RC(60)*Y(1)
      JAC(2,61) = 0.19D1*RC(60)*Y(1)
      JAC(3,3) = -RC(39)*Y(37)-RC(40)*Y(19)-RC(47)
      JAC(3,19) = -RC(40)*Y(3)
      JAC(3,37) = -RC(39)*Y(3)
      JAC(4,4) = -RC(70)*Y(37)
      JAC(4,8) = 0.44D0*RC(112)*Y(14)
      JAC(4,9) = 0.4D0*RC(123)*Y(14)
      JAC(4,11) = RC(66)*Y(37)+DJ(6)+DJ(7)+RC(69)*Y(39)
      JAC(4,12) = DJ(8)
      JAC(4,14) = 0.44D0*RC(112)*Y(8)+0.4D0*RC(123)*Y(9)+0.5D-1*RC(160)*
     &Y(44)+0.5D-1*RC(150)*Y(41)
      JAC(4,30) = DJ(11)+RC(222)*Y(37)
      JAC(4,32) = 2.D0*RC(221)*Y(37)+2.D0*DJ(7)
      JAC(4,37) = RC(66)*Y(11)+2.D0*RC(221)*Y(32)+RC(222)*Y(30)-RC(70)*Y
     &(4)
      JAC(4,39) = RC(69)*Y(11)
      JAC(4,41) = 0.5D-1*RC(150)*Y(14)
      JAC(4,44) = 0.5D-1*RC(160)*Y(14)
      JAC(5,5) = -RC(59)*Y(37)
      JAC(5,9) = 0.7D-1*RC(123)*Y(14)
      JAC(5,14) = 0.7D-1*RC(123)*Y(9)
      JAC(5,37) = -RC(59)*Y(5)
      JAC(6,6) = -RC(71)*Y(37)
      JAC(6,37) = -RC(71)*Y(6)
      JAC(7,7) = -RC(81)*Y(37)
      JAC(7,37) = -RC(81)*Y(7)
      JAC(8,8) = -RC(109)*Y(37)-RC(112)*Y(14)
      JAC(8,14) = -RC(112)*Y(8)
      JAC(8,37) = -RC(109)*Y(8)
      JAC(9,9) = -RC(123)*Y(14)-RC(125)*Y(37)
      JAC(9,14) = 0.7D-1*RC(150)*Y(41)-RC(123)*Y(9)
      JAC(9,37) = -RC(125)*Y(9)
      JAC(9,41) = 0.7D-1*RC(150)*Y(14)
      JAC(10,10) = -RC(234)*Y(37)
      JAC(10,37) = -RC(234)*Y(10)
      JAC(11,1) = Y(19)*RC(60)+RC(60)*(2.D0*Y(29)+Y(31)+0.74D0*Y(43)+0.2
     &66D0*Y(45)+0.15D0*Y(60))
      JAC(11,3) = RC(40)*Y(19)
      JAC(11,8) = RC(112)*Y(14)
      JAC(11,9) = 0.5D0*RC(123)*Y(14)
      JAC(11,11) = -RC(66)*Y(37)-DJ(6)-DJ(7)-RC(69)*Y(39)-RC(53)
      JAC(11,14) = 0.5D0*RC(123)*Y(9)+RC(112)*Y(8)+0.7D0*RC(157)*Y(56)+0
     &.8D0*RC(160)*Y(44)+0.8D0*RC(150)*Y(41)
      JAC(11,19) = RC(60)*Y(1)+2*(2.D0*RC(61)+RC(62))*Y(19)+RC(80)*Y(24)
     &+RC(40)*Y(3)
      JAC(11,22) = Y(37)*RC(67)+DJ(16)
      JAC(11,24) = RC(80)*Y(19)
      JAC(11,29) = 2.D0*RC(60)*Y(1)
      JAC(11,31) = RC(60)*Y(1)
      JAC(11,32) = 2.D0*DJ(7)
      JAC(11,37) = RC(63)*Y(46)+RC(67)*Y(22)-RC(66)*Y(11)
      JAC(11,39) = -RC(69)*Y(11)
      JAC(11,41) = 0.8D0*RC(150)*Y(14)
      JAC(11,43) = 0.74D0*RC(60)*Y(1)
      JAC(11,44) = 0.8D0*RC(160)*Y(14)
      JAC(11,45) = 0.266D0*RC(60)*Y(1)
      JAC(11,46) = Y(37)*RC(63)
      JAC(11,50) = 0.156D1*DJ(16)
      JAC(11,51) = DJ(16)
      JAC(11,56) = 0.7D0*RC(157)*Y(14)
      JAC(11,60) = 0.15D0*RC(60)*Y(1)
      JAC(12,1) = RC(72)*Y(23)+RC(83)*Y(26)*RPATH4+RC(105)*Y(27)+RC(126)
     &*Y(31)+0.95D0*RC(162)*Y(61)+0.684D0*RC(154)*Y(45)
      JAC(12,9) = 0.5D0*RC(123)*Y(14)
      JAC(12,12) = -DJ(8)-RC(75)*Y(37)-RC(53)
      JAC(12,14) = 0.5D0*RC(123)*Y(9)+0.4D-1*RC(160)*Y(44)
      JAC(12,20) = Y(37)*RC(64)
      JAC(12,23) = RC(72)*Y(1)
      JAC(12,26) = Y(1)*RC(83)*RPATH4
      JAC(12,27) = RC(105)*Y(1)
      JAC(12,28) = RC(76)*Y(37)+DJ(16)
      JAC(12,31) = RC(126)*Y(1)
      JAC(12,37) = RC(64)*Y(20)+RC(76)*Y(28)+RC(76)*Y(50)-RC(75)*Y(12)
      JAC(12,44) = 0.4D-1*RC(160)*Y(14)
      JAC(12,45) = 0.684D0*RC(154)*Y(1)
      JAC(12,49) = 0.35D0*DJ(16)
      JAC(12,50) = RC(76)*Y(37)+0.22D0*DJ(16)
      JAC(12,51) = DJ(16)
      JAC(12,52) = DJ(16)
      JAC(12,61) = 0.95D0*RC(162)*Y(1)
      JAC(13,1) = RC(83)*Y(26)*RPATH3+RC(159)*Y(59)+0.95D0*RC(162)*Y(61)
      JAC(13,13) = -DJ(9)-RC(86)*Y(37)-RC(53)
      JAC(13,26) = RC(83)*Y(1)*RPATH3
      JAC(13,37) = RC(76)*Y(49)+RC(76)*Y(51)-RC(86)*Y(13)
      JAC(13,49) = 0.65D0*DJ(16)+RC(76)*Y(37)
      JAC(13,51) = RC(76)*Y(37)
      JAC(13,59) = RC(159)*Y(1)+DJ(16)
      JAC(13,61) = 0.95D0*RC(162)*Y(1)
      JAC(14,1) = -RC(11)*Y(14)
      JAC(14,2) = -RC(12)*Y(14)
      JAC(14,8) = -RC(112)*Y(14)
      JAC(14,9) = -RC(123)*Y(14)
      JAC(14,14) = -RC(11)*Y(1)-RC(12)*Y(2)-RC(13)*Y(37)-RC(14)*Y(15)-RC
     &(49)-RC(112)*Y(8)-RC(123)*Y(9)-RC(157)*Y(56)-RC(160)*Y(44)-RC(150)
     &*Y(41)-DJ(1)-DJ(2)
      JAC(14,15) = RC(89)*Y(24)-RC(14)*Y(14)
      JAC(14,24) = RC(89)*Y(15)
      JAC(14,36) = RC(1)
      JAC(14,37) = -RC(13)*Y(14)
      JAC(14,41) = -RC(150)*Y(14)
      JAC(14,44) = -RC(160)*Y(14)
      JAC(14,56) = -RC(157)*Y(14)
      JAC(15,1) = RC(60)*(Y(19)+Y(29)+Y(31)+Y(33)+Y(35)+0.95D0*Y(45)+Y(2
     &6)*RPATH3+0.78D0*Y(43)+Y(59)+0.5D-1*Y(61)+0.8D0*Y(60))+RC(72)*Y(23
     &)-RC(17)*Y(15)
      JAC(15,3) = RC(40)*Y(19)+RC(39)*Y(37)
      JAC(15,4) = RC(70)*Y(37)
      JAC(15,8) = 0.12D0*RC(112)*Y(14)
      JAC(15,9) = 0.28D0*RC(123)*Y(14)
      JAC(15,11) = RC(66)*Y(37)+2.D0*DJ(6)+RC(69)*Y(39)
      JAC(15,12) = DJ(8)
      JAC(15,14) = RC(13)*Y(37)+0.12D0*RC(112)*Y(8)+0.28D0*RC(123)*Y(9)+
     &0.6D-1*RC(160)*Y(44)+0.6D-1*RC(150)*Y(41)-RC(14)*Y(15)
      JAC(15,15) = -4.D0*RC(36)*Y(15)-RC(14)*Y(14)-RC(17)*Y(1)-RC(30)*Y(
     &37)-RC(65)*Y(19)-RC(74)*Y(23)-(RC(88)+RC(89))*Y(24)-RC(85)*(Y(26)+
     &Y(29)+Y(31)+Y(27)+Y(57)+Y(45)+Y(61)+Y(59)+Y(33)+Y(35)+Y(43)+Y(60))
      JAC(15,17) = RC(31)*Y(37)+RC(26)*Y(39)
      JAC(15,18) = Y(37)*RC(33)
      JAC(15,19) = RC(40)*Y(3)+4.D0*RC(61)*Y(19)+0.5D0*RC(80)*Y(24)+RC(6
     &0)*Y(1)-RC(65)*Y(15)
      JAC(15,20) = Y(37)*RC(64)
      JAC(15,22) = DJ(16)
      JAC(15,23) = RC(72)*Y(1)-RC(74)*Y(15)
      JAC(15,24) = 0.5D0*RC(80)*Y(19)-(RC(88)+RC(89))*Y(15)
      JAC(15,26) = Y(1)*RC(60)*RPATH3-RC(85)*Y(15)
      JAC(15,27) = -RC(85)*Y(15)
      JAC(15,28) = DJ(16)
      JAC(15,29) = RC(60)*Y(1)-RC(85)*Y(15)
      JAC(15,30) = DJ(11)
      JAC(15,31) = RC(60)*Y(1)-RC(85)*Y(15)
      JAC(15,32) = RC(221)*Y(37)
      JAC(15,33) = RC(60)*Y(1)-RC(85)*Y(15)
      JAC(15,35) = RC(60)*Y(1)-RC(85)*Y(15)
      JAC(15,37) = RC(13)*Y(14)+RC(31)*Y(17)+RC(33)*Y(18)+RC(39)*Y(3)+RC
     &(63)*Y(46)+RC(64)*Y(20)+RC(66)*Y(11)+RC(70)*Y(4)+RC(221)*Y(32)-RC(
     &30)*Y(15)
      JAC(15,39) = RC(26)*Y(17)+RC(69)*Y(11)
      JAC(15,41) = 0.6D-1*RC(150)*Y(14)
      JAC(15,43) = 0.78D0*RC(60)*Y(1)-RC(85)*Y(15)
      JAC(15,44) = 0.6D-1*RC(160)*Y(14)
      JAC(15,45) = 0.95D0*RC(60)*Y(1)-RC(85)*Y(15)
      JAC(15,46) = Y(37)*RC(63)
      JAC(15,48) = DJ(16)
      JAC(15,49) = 0.65D0*DJ(16)
      JAC(15,50) = DJ(16)
      JAC(15,51) = DJ(16)
      JAC(15,53) = DJ(16)
      JAC(15,57) = -RC(85)*Y(15)
      JAC(15,59) = RC(60)*Y(1)-RC(85)*Y(15)
      JAC(15,60) = 0.8D0*RC(60)*Y(1)-RC(85)*Y(15)
      JAC(15,61) = 0.5D-1*RC(60)*Y(1)-RC(85)*Y(15)
      JAC(16,2) = RC(21)*Y(37)
      JAC(16,11) = RC(69)*Y(39)
      JAC(16,16) = -RC(35)*Y(37)-DJ(5)-RC(45)
      JAC(16,17) = RC(26)*Y(39)
      JAC(16,37) = RC(21)*Y(2)-RC(35)*Y(16)
      JAC(16,39) = RC(26)*Y(17)+RC(69)*Y(11)
      JAC(17,15) = 2*RC(36)*Y(15)
      JAC(17,17) = -RC(31)*Y(37)-DJ(4)-RC(43)-RC(26)*Y(39)-RC(47)
      JAC(17,37) = -RC(31)*Y(17)
      JAC(17,39) = -RC(26)*Y(17)
      JAC(18,8) = 0.13D0*RC(112)*Y(14)
      JAC(18,9) = 0.7D-1*RC(123)*Y(14)
      JAC(18,11) = DJ(7)
      JAC(18,14) = 0.13D0*RC(112)*Y(8)+0.7D-1*RC(123)*Y(9)
      JAC(18,18) = -Y(37)*RC(33)
      JAC(18,37) = -RC(33)*Y(18)
      JAC(19,1) = Y(24)*RC(79)-Y(19)*RC(60)
      JAC(19,3) = -RC(40)*Y(19)
      JAC(19,5) = RC(59)*Y(37)
      JAC(19,9) = 0.31D0*RC(123)*Y(14)
      JAC(19,12) = DJ(8)
      JAC(19,14) = 0.31D0*RC(123)*Y(9)
      JAC(19,15) = -RC(65)*Y(19)
      JAC(19,19) = -(2.D0*RC(61)+2.D0*RC(62))*Y(19)-RC(40)*Y(3)-RC(60)*Y
     &(1)-2.D0*RC(61)*Y(19)-2.D0*RC(62)*Y(19)-RC(65)*Y(15)-0.5D0*RC(80)*
     &Y(24)
      JAC(19,22) = RC(68)*Y(37)
      JAC(19,24) = RC(79)*Y(1)+4.D0*RC(94)*Y(24)-0.5D0*RC(80)*Y(19)
      JAC(19,37) = RC(59)*Y(5)+RC(68)*Y(22)
      JAC(19,47) = DJ(16)
      JAC(20,20) = -Y(37)*RC(64)
      JAC(20,37) = -RC(64)*Y(20)
      JAC(21,3) = RC(40)*Y(19)+RC(39)*Y(37)
      JAC(21,19) = RC(40)*Y(3)
      JAC(21,21) = -RC(51)
      JAC(21,37) = RC(39)*Y(3)
      JAC(22,15) = RC(65)*Y(19)
      JAC(22,19) = RC(65)*Y(15)
      JAC(22,22) = -RC(43)-DJ(16)-(RC(67)+RC(68))*Y(37)
      JAC(22,37) = -(RC(67)+RC(68))*Y(22)
      JAC(23,1) = RC(83)*Y(26)*RPATH4-RC(72)*Y(23)
      JAC(23,6) = RC(71)*Y(37)
      JAC(23,13) = DJ(9)
      JAC(23,15) = -RC(74)*Y(23)
      JAC(23,23) = -RC(72)*Y(1)-RC(74)*Y(15)
      JAC(23,26) = Y(1)*RC(83)*RPATH4
      JAC(23,28) = RC(68)*Y(37)
      JAC(23,37) = RC(71)*Y(6)+RC(68)*Y(28)
      JAC(23,49) = 0.35D0*DJ(16)
      JAC(24,1) = RC(105)*Y(27)+0.684D0*RC(154)*Y(45)-Y(24)*RC(79)
      JAC(24,2) = -RC(77)*Y(24)
      JAC(24,12) = RC(75)*Y(37)
      JAC(24,13) = DJ(9)
      JAC(24,15) = -(RC(88)+RC(89))*Y(24)
      JAC(24,19) = -RC(80)*Y(24)
      JAC(24,24) = -4.D0*RC(94)*Y(24)-RC(77)*Y(2)-RC(79)*Y(1)-RC(80)*Y(1
     &9)-(RC(88)+RC(89))*Y(15)
      JAC(24,25) = RC(78)
      JAC(24,27) = RC(105)*Y(1)
      JAC(24,30) = DJ(11)+RC(222)*Y(37)
      JAC(24,37) = RC(75)*Y(12)+RC(222)*Y(30)+RC(68)*Y(47)
      JAC(24,45) = 0.684D0*RC(154)*Y(1)
      JAC(24,47) = RC(68)*Y(37)
      JAC(24,52) = DJ(16)
      JAC(25,2) = RC(77)*Y(24)
      JAC(25,24) = RC(77)*Y(2)
      JAC(25,25) = -RC(50)-RC(78)
      JAC(26,1) = -RC(83)*Y(26)
      JAC(26,7) = RC(81)*Y(37)
      JAC(26,15) = -RC(85)*Y(26)
      JAC(26,26) = -RC(83)*Y(1)-RC(85)*Y(15)
      JAC(26,37) = RC(81)*Y(7)+RC(68)*Y(49)
      JAC(26,49) = RC(68)*Y(37)
      JAC(27,1) = -RC(105)*Y(27)
      JAC(27,13) = RC(86)*Y(37)
      JAC(27,15) = -RC(85)*Y(27)
      JAC(27,27) = -RC(105)*Y(1)-RC(85)*Y(15)
      JAC(27,37) = RC(86)*Y(13)+RC(87)*Y(52)
      JAC(27,52) = RC(87)*Y(37)
      JAC(28,15) = RC(74)*Y(23)
      JAC(28,23) = RC(74)*Y(15)
      JAC(28,28) = -(RC(76)+RC(68))*Y(37)-DJ(16)-RC(52)
      JAC(28,37) = -(RC(76)+RC(68))*Y(28)
      JAC(29,1) = -RC(110)*Y(29)
      JAC(29,8) = RC(109)*Y(37)
      JAC(29,15) = -RC(85)*Y(29)
      JAC(29,29) = -RC(110)*Y(1)-RC(85)*Y(15)
      JAC(29,37) = RC(109)*Y(8)+RC(68)*Y(50)
      JAC(29,50) = RC(68)*Y(37)
      JAC(30,1) = RC(236)*Y(33)+RC(220)*Y(35)+0.266D0*RC(154)*Y(45)
      JAC(30,14) = 0.82D0*RC(160)*Y(44)
      JAC(30,30) = -DJ(11)-RC(222)*Y(37)
      JAC(30,33) = RC(236)*Y(1)
      JAC(30,35) = RC(220)*Y(1)
      JAC(30,37) = -RC(222)*Y(30)
      JAC(30,44) = 0.82D0*RC(160)*Y(14)
      JAC(30,45) = 0.266D0*RC(154)*Y(1)
      JAC(30,48) = DJ(16)
      JAC(30,53) = DJ(16)
      JAC(31,1) = -RC(126)*Y(31)
      JAC(31,9) = RC(125)*Y(37)
      JAC(31,15) = -RC(85)*Y(31)
      JAC(31,31) = -RC(126)*Y(1)-RC(85)*Y(15)
      JAC(31,37) = RC(125)*Y(9)+RC(68)*Y(51)
      JAC(31,51) = RC(68)*Y(37)
      JAC(32,1) = RC(220)*Y(35)
      JAC(32,32) = -2.D0*DJ(7)-RC(221)*Y(37)
      JAC(32,35) = RC(220)*Y(1)
      JAC(32,37) = -RC(221)*Y(32)
      JAC(32,53) = DJ(16)
      JAC(33,1) = -RC(236)*Y(33)
      JAC(33,10) = RC(234)*Y(37)
      JAC(33,15) = -RC(85)*Y(33)
      JAC(33,33) = -RC(236)*Y(1)-RC(85)*Y(15)
      JAC(33,37) = RC(234)*Y(10)+RC(235)*Y(48)
      JAC(33,48) = RC(235)*Y(37)
      JAC(34,1) = RC(236)*Y(33)
      JAC(34,33) = RC(236)*Y(1)
      JAC(34,34) = -RC(219)*Y(37)
      JAC(34,37) = -RC(219)*Y(34)
      JAC(34,48) = DJ(16)
      JAC(35,1) = -RC(220)*Y(35)
      JAC(35,15) = -RC(85)*Y(35)
      JAC(35,34) = RC(219)*Y(37)
      JAC(35,35) = -RC(220)*Y(1)-RC(85)*Y(15)
      JAC(35,37) = RC(219)*Y(34)+RC(223)*Y(53)
      JAC(35,53) = RC(223)*Y(37)
      JAC(36,1) = -RC(5)*Y(36)
      JAC(36,2) = DJ(3)
      JAC(36,14) = DJ(1)+0.2D0*RC(160)*Y(44)+0.3D0*RC(150)*Y(41)
      JAC(36,36) = -RC(1)-RC(5)*Y(1)
      JAC(36,38) = RC(7)
      JAC(36,39) = DJ(14)
      JAC(36,41) = 0.3D0*RC(150)*Y(14)
      JAC(36,44) = 0.2D0*RC(160)*Y(14)
      JAC(37,1) = RC(17)*Y(15)
      JAC(37,2) = -RC(21)*Y(37)
      JAC(37,3) = -RC(39)*Y(37)
      JAC(37,4) = -RC(70)*Y(37)
      JAC(37,5) = -RC(59)*Y(37)
      JAC(37,6) = -RC(71)*Y(37)
      JAC(37,7) = -RC(81)*Y(37)
      JAC(37,8) = -RC(109)*Y(37)
      JAC(37,9) = 0.15D0*RC(123)*Y(14)-RC(125)*Y(37)
      JAC(37,10) = -RC(234)*Y(37)
      JAC(37,11) = -RC(66)*Y(37)
      JAC(37,12) = -RC(75)*Y(37)
      JAC(37,13) = -RC(86)*Y(37)
      JAC(37,14) = RC(14)*Y(15)+0.15D0*RC(123)*Y(9)+0.8D-1*RC(160)*Y(44)
     &+0.55D0*RC(150)*Y(41)-RC(13)*Y(37)
      JAC(37,15) = RC(14)*Y(14)+RC(17)*Y(1)-RC(30)*Y(37)
      JAC(37,16) = DJ(5)-RC(35)*Y(37)
      JAC(37,17) = 2.D0*DJ(4)-RC(31)*Y(37)
      JAC(37,18) = -Y(37)*RC(33)
      JAC(37,20) = -Y(37)*RC(64)
      JAC(37,22) = DJ(16)-RC(68)*Y(37)
      JAC(37,28) = DJ(16)-RC(68)*Y(37)
      JAC(37,30) = -RC(222)*Y(37)
      JAC(37,32) = -RC(221)*Y(37)
      JAC(37,34) = -RC(219)*Y(37)
      s1 = -RC(149)*(Y(65)+Y(66))-RC(21)*Y(2)-RC(35)*Y(16)-RC(87)*Y(52)-
     &RC(147)*Y(64)-RC(75)*Y(12)-RC(148)*Y(62)-RC(31)*Y(17)-RC(64)*Y(20)
     &-RC(70)*Y(4)-RC(153)*Y(44)-RC(86)*Y(13)-RC(13)*Y(14)-RC(109)*Y(8)-
     &RC(235)*Y(48)-RC(234)*Y(10)-RC(81)*Y(7)
      JAC(37,37) = s1-RC(59)*Y(5)-RC(219)*Y(34)-RC(30)*Y(15)-RC(221)*Y(3
     &2)-RC(63)*Y(46)-RC(66)*Y(11)-RC(222)*Y(30)-RC(33)*Y(18)-RC(39)*Y(3
     &)-RC(223)*Y(53)-RC(151)*Y(41)-RC(71)*Y(6)-RC(158)*Y(54)-RC(161)*Y(
     &55)-RC(125)*Y(9)-RC(68)*(Y(22)+Y(28)+Y(47)+Y(50)+Y(51)+Y(49))-RC(1
     &46)*Y(63)
      JAC(37,38) = 2.D0*RC(8)*H2O
      JAC(37,41) = 0.55D0*RC(150)*Y(14)-RC(151)*Y(37)
      JAC(37,44) = 0.8D-1*RC(160)*Y(14)-RC(153)*Y(37)
      JAC(37,46) = -Y(37)*RC(63)
      JAC(37,47) = DJ(16)-RC(68)*Y(37)
      JAC(37,48) = DJ(16)-RC(235)*Y(37)
      JAC(37,49) = DJ(16)-RC(68)*Y(37)
      JAC(37,50) = DJ(16)-RC(68)*Y(37)
      JAC(37,51) = -RC(68)*Y(37)
      JAC(37,52) = DJ(16)-RC(87)*Y(37)
      JAC(37,53) = DJ(16)-RC(223)*Y(37)
      JAC(37,54) = -RC(158)*Y(37)
      JAC(37,55) = -RC(161)*Y(37)
      JAC(37,62) = -RC(148)*Y(37)
      JAC(37,63) = -RC(146)*Y(37)
      JAC(37,64) = -RC(147)*Y(37)
      JAC(37,65) = -RC(149)*Y(37)
      JAC(37,66) = -RC(149)*Y(37)
      JAC(38,14) = DJ(2)
      JAC(38,38) = -RC(7)-RC(8)*H2O
      JAC(39,1) = -RC(15)*Y(39)
      JAC(39,2) = RC(12)*Y(14)-(RC(19)+RC(20))*Y(39)
      JAC(39,11) = -RC(69)*Y(39)
      JAC(39,14) = RC(12)*Y(2)
      JAC(39,16) = RC(35)*Y(37)
      JAC(39,17) = -RC(26)*Y(39)
      JAC(39,37) = RC(35)*Y(16)
      JAC(39,39) = -RC(15)*Y(1)-RC(26)*Y(17)-RC(163)*Y(41)-RC(19)*Y(2)-R
     &C(20)*Y(2)-DJ(13)-DJ(14)-RC(69)*Y(11)
      JAC(39,40) = RC(29)+DJ(15)
      JAC(39,41) = -RC(163)*Y(39)
      JAC(40,2) = RC(20)*Y(39)
      JAC(40,39) = RC(20)*Y(2)
      JAC(40,40) = -RC(29)-DJ(15)-RC(45)
      JAC(41,14) = -RC(150)*Y(41)
      JAC(41,37) = -RC(151)*Y(41)
      JAC(41,39) = -RC(163)*Y(41)
      JAC(41,41) = -RC(151)*Y(37)-RC(163)*Y(39)-RC(150)*Y(14)
      JAC(42,16) = RC(45)
      JAC(42,40) = 2.D0*RC(44)
      JAC(42,42) = -RC(51)
      JAC(43,1) = -0.88D0*RC(152)*Y(43)
      JAC(43,15) = -RC(155)*Y(43)
      JAC(43,37) = RC(151)*Y(41)+RC(156)*Y(56)
      JAC(43,41) = RC(151)*Y(37)
      JAC(43,43) = -0.88D0*RC(152)*Y(1)-RC(155)*Y(15)
      JAC(43,56) = RC(156)*Y(37)
      JAC(44,1) = RC(60)*(0.42D0*Y(43)+0.5D-1*Y(60))
      JAC(44,14) = 0.26D0*RC(150)*Y(41)-RC(160)*Y(44)
      JAC(44,37) = -RC(153)*Y(44)
      JAC(44,41) = 0.26D0*RC(150)*Y(14)
      JAC(44,43) = 0.42D0*RC(60)*Y(1)
      JAC(44,44) = -RC(153)*Y(37)-RC(160)*Y(14)
      JAC(44,60) = 0.5D-1*RC(60)*Y(1)
      JAC(45,1) = -RC(154)*Y(45)
      JAC(45,15) = -RC(85)*Y(45)
      JAC(45,37) = RC(153)*Y(44)+RC(148)*Y(62)
      JAC(45,44) = RC(153)*Y(37)
      JAC(45,45) = -RC(154)*Y(1)-RC(85)*Y(15)
      JAC(45,62) = RC(148)*Y(37)
      JAC(46,19) = 2*RC(62)*Y(19)
      JAC(46,37) = -RC(63)*Y(46)
      JAC(46,46) = -Y(37)*RC(63)
      JAC(47,15) = RC(88)*Y(24)
      JAC(47,24) = RC(88)*Y(15)
      JAC(47,37) = -RC(68)*Y(47)
      JAC(47,47) = -RC(68)*Y(37)-DJ(16)-RC(52)
      JAC(48,15) = RC(85)*Y(33)
      JAC(48,33) = RC(85)*Y(15)
      JAC(48,37) = -RC(235)*Y(48)
      JAC(48,48) = -RC(235)*Y(37)-DJ(16)-RC(52)
      JAC(49,15) = RC(85)*Y(26)
      JAC(49,26) = RC(85)*Y(15)
      JAC(49,37) = -(RC(76)+RC(68))*Y(49)
      JAC(49,49) = -(RC(76)+RC(68))*Y(37)-DJ(16)-RC(52)
      JAC(50,15) = RC(85)*Y(29)
      JAC(50,29) = RC(85)*Y(15)
      JAC(50,37) = -(RC(76)+RC(68))*Y(50)
      JAC(50,50) = -(RC(76)+RC(68))*Y(37)-DJ(16)-RC(52)
      JAC(51,15) = RC(85)*Y(31)
      JAC(51,31) = RC(85)*Y(15)
      JAC(51,37) = -(RC(76)+RC(68))*Y(51)
      JAC(51,51) = -(RC(76)+RC(68))*Y(37)-DJ(16)-RC(52)
      JAC(52,15) = RC(85)*Y(27)
      JAC(52,27) = RC(85)*Y(15)
      JAC(52,37) = -RC(87)*Y(52)
      JAC(52,52) = -RC(87)*Y(37)-DJ(16)-RC(52)
      JAC(53,15) = RC(85)*Y(35)
      JAC(53,35) = RC(85)*Y(15)
      JAC(53,37) = -RC(223)*Y(53)
      JAC(53,53) = -RC(223)*Y(37)-DJ(16)-RC(52)
      JAC(54,1) = RC(60)*(0.32D0*Y(43)+0.1D0*Y(60))
      JAC(54,14) = 0.67D0*RC(150)*Y(41)
      JAC(54,37) = -RC(158)*Y(54)
      JAC(54,41) = 0.67D0*RC(150)*Y(14)
      JAC(54,43) = 0.32D0*RC(60)*Y(1)
      JAC(54,54) = -RC(158)*Y(37)
      JAC(54,60) = 0.1D0*RC(60)*Y(1)
      JAC(55,1) = RC(60)*(0.14D0*Y(43)+0.5D-1*Y(45)+0.85D0*Y(60))
      JAC(55,37) = -RC(161)*Y(55)
      JAC(55,43) = 0.14D0*RC(60)*Y(1)
      JAC(55,45) = 0.5D-1*RC(60)*Y(1)
      JAC(55,55) = -RC(161)*Y(37)
      JAC(55,60) = 0.85D0*RC(60)*Y(1)
      JAC(56,14) = -RC(157)*Y(56)
      JAC(56,15) = RC(155)*Y(43)
      JAC(56,37) = -RC(156)*Y(56)
      JAC(56,43) = RC(155)*Y(15)
      JAC(56,56) = -RC(156)*Y(37)-RC(157)*Y(14)-RC(52)
      JAC(57,1) = -RC(79)*Y(57)
      JAC(57,2) = -RC(77)*Y(57)
      JAC(57,15) = -RC(85)*Y(57)
      JAC(57,37) = 0.5D0*RC(158)*Y(54)+RC(149)*Y(66)
      JAC(57,54) = 0.5D0*RC(158)*Y(37)
      JAC(57,57) = -RC(77)*Y(2)-RC(79)*Y(1)-RC(85)*Y(15)
      JAC(57,58) = RC(78)
      JAC(57,66) = RC(149)*Y(37)
      JAC(58,2) = RC(77)*Y(57)
      JAC(58,57) = RC(77)*Y(2)
      JAC(58,58) = -RC(50)-RC(78)
      JAC(59,1) = RC(79)*Y(57)-RC(159)*Y(59)
      JAC(59,15) = -RC(85)*Y(59)
      JAC(59,37) = RC(146)*Y(63)
      JAC(59,57) = RC(79)*Y(1)
      JAC(59,59) = -RC(159)*Y(1)-RC(85)*Y(15)
      JAC(59,63) = RC(146)*Y(37)
      JAC(60,1) = -RC(164)*Y(60)
      JAC(60,15) = -RC(85)*Y(60)
      JAC(60,37) = RC(147)*Y(64)
      JAC(60,39) = RC(163)*Y(41)
      JAC(60,41) = RC(163)*Y(39)
      JAC(60,60) = -RC(164)*Y(1)-RC(85)*Y(15)
      JAC(60,64) = RC(147)*Y(37)
      JAC(61,1) = -RC(162)*Y(61)
      JAC(61,15) = -RC(85)*Y(61)
      JAC(61,37) = RC(161)*Y(55)+RC(149)*Y(65)
      JAC(61,55) = RC(161)*Y(37)
      JAC(61,61) = -RC(162)*Y(1)-RC(85)*Y(15)
      JAC(61,65) = RC(149)*Y(37)
      JAC(62,15) = RC(85)*Y(45)
      JAC(62,37) = -RC(148)*Y(62)
      JAC(62,45) = RC(85)*Y(15)
      JAC(62,62) = -RC(148)*Y(37)-RC(52)
      JAC(63,15) = RC(85)*Y(59)
      JAC(63,37) = -RC(146)*Y(63)
      JAC(63,59) = RC(85)*Y(15)
      JAC(63,63) = -RC(146)*Y(37)-RC(52)
      JAC(64,15) = RC(85)*Y(60)
      JAC(64,37) = -RC(147)*Y(64)
      JAC(64,60) = RC(85)*Y(15)
      JAC(64,64) = -RC(147)*Y(37)-RC(52)
      JAC(65,15) = RC(85)*Y(61)
      JAC(65,37) = -RC(149)*Y(65)
      JAC(65,61) = RC(85)*Y(15)
      JAC(65,65) = -RC(149)*Y(37)-RC(52)
      JAC(66,15) = RC(85)*Y(57)
      JAC(66,37) = -RC(149)*Y(66)
      JAC(66,57) = RC(85)*Y(15)
      JAC(66,66) = -RC(149)*Y(37)-RC(52)
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine meval(ldim,neqn,t,y,yprime,dfddy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfddy(ldim,neqn),rpar(*)
c
c     dummy subroutine
c
      return
      end
c-----------------------------------------------------------------------
      subroutine solut(neqn,t,y)
      integer neqn
      double precision t,y(neqn)
C
C  RADAU5 applied to EMEP problem, Tend = 417600
C
C  relative error tolerance = 0.1000D-13,
C  absolute error tolerance = 0.1000D-07
C  and initial stepsize = 0.1000D-07
C  
C
C    number of integration steps      17817
C    number of accepted steps         17449
C    number of f evaluations         124335
C    number of Jacobian evaluations    5036
C    number of LU decompositions       9721
C
C
  
      y(  1) =  0.2564580511140732d+008      
      y(  2) =  0.5146134770952715d+011   
      y(  3) =  0.2315679957701715d+012     
      y(  4) =  0.1130936599472892d+014      
      y(  5) =  0.3459285326034955d+014      
      y(  6) =  0.1027236550974901d+012      
      y(  7) =  0.8508735586836855d+011   
      y(  8) =  0.4131285674624012d+010     
      y(  9) =  0.1270937854717943d+010     
      y( 10) =  0.5632880890404914d+010     
      y( 11) =  0.8263821552730888d+011    
      y( 12) =  0.3293552658908353d+011      
      y( 13) =  0.1054058186206315d+012      
      y( 14) =  0.3150308585365321d+013     
      y( 15) =  0.2488383936633755d+008    
      y( 16) =  0.1097565615411556d+011      
      y( 17) =  0.1365196557629180d+011      
      y( 18) =  0.3852048230503094d+012     
      y( 19) =  0.1137462721600089d+009    
      y( 20) =  0.2850982205905218d+011      
      y( 21) =  0.3786933619791445d+012    
      y( 22) =  0.3832384333863027d+010     
      y( 23) =  0.4804939902203071d+007      
      y( 24) =  0.3390546081092960d+008      
      y( 25) =  0.3494452929972591d+011     
      y( 26) =  0.1490576370801779d+008      
      y( 27) =  0.1000567871868853d+008     
      y( 28) =  0.1856753061984312d+010    
      y( 29) =  0.2272206156875160d+007     
      y( 30) =  0.5168375883757783d+010    
      y( 31) =  0.2317595952406701d+007    
      y( 32) =  0.3035139656996921d+010   
      y( 33) =  0.5928979479102226d+007      
      y( 34) =  0.2515297290687841d+010    
      y( 35) =  0.4569678276203798d+007      
      y( 36) =  0.3727926507368739d+001                
      y( 37) =  0.4063589349848112d+005    
      y( 38) =  0.8159151146269279d-037             
      y( 39) =  0.5577590724412284d+009      
      y( 40) =  0.7684596616753747d+009   
      y( 41) =  0.5338491777959816d+010      
      y( 42) =  0.6048604025328409d+012    
      y( 43) =  0.4222237216819787d+008      
      y( 44) =  0.3681852767009784d+010      
      y( 45) =  0.5835731567978018d+007   
      y( 46) =  0.1456507159081862d+010   
      y( 47) =  0.8862868684014137d+010    
      y( 48) =  0.9338753312528582d+009    
      y( 49) =  0.3595753586682656d+010    
      y( 50) =  0.4622427584057254d+009      
      y( 51) =  0.2395910064875511d+009     
      y( 52) =  0.6080936586253302d+010      
      y( 53) =  0.7704108297342240d+009    
      y( 54) =  0.3386220860753221d+010     
      y( 55) =  0.1838264387030593d+010      
      y( 56) =  0.1973252346091884d+010    
      y( 57) =  0.1313403940892729d+007      
      y( 58) =  0.1344471764571674d+010     
      y( 59) =  0.2568483117244701d+007      
      y( 60) =  0.1131854591900057d+010   
      y( 61) =  0.4591537937818946d+007      
      y( 62) =  0.8609382282952802d+009      
      y( 63) =  0.3509560971283513d+009   
      y( 64) =  0.6682350971104786d+009      
      y( 65) =  0.3207356403011160d+009     
      y( 66) =  0.1144987634882048d+009     


      RETURN
      END
c-----------------------------------------------------------------------
      subroutine solutold(neqn,t,y)
      integer neqn
      double precision t,y(neqn)
C
C  RADAU5 applied to EMEP problem, Tend = 417600
C
C  relative error tolerance = 0.1000D-11,
C  absolute error tolerance = 0.1000D+01
C  and initial stepsize = 0.1000D-09
C  Max stepsize    =   10.0
C
C  steps =    41401   (   41200 accepted and    201 rejected)
C  f-val =   169733
C   njac =      581
C    dec =     1069
C fb-sub =    42817
C

      Y( 1) =   0.25645805093601D+08
      Y( 2) =   0.51461347708556D+11
      Y( 3) =   0.23156799577319D+12
      Y( 4) =   0.11309365994733D+14
      Y( 5) =   0.34592853260350D+14
      Y( 6) =   0.10272365509751D+12
      Y( 7) =   0.85087355868443D+11
      Y( 8) =   0.41312856746209D+10
      Y( 9) =   0.12709378547180D+10
      Y(10) =   0.56328808904060D+10
      Y(11) =   0.82638215527479D+11
      Y(12) =   0.32935526589172D+11
      Y(13) =   0.10540581862144D+12
      Y(14) =   0.31503085853931D+13
      Y(15) =   0.24883839367388D+08
      Y(16) =   0.10975656154077D+11
      Y(17) =   0.13651965576432D+11
      Y(18) =   0.38520482305080D+12
      Y(19) =   0.11374627216097D+09
      Y(20) =   0.28509822059077D+11
      Y(21) =   0.37869336198310D+12
      Y(22) =   0.38323843339068D+10
      Y(23) =   0.48049399025208D+07
      Y(24) =   0.33905460811715D+08
      Y(25) =   0.34944529300020D+11
      Y(26) =   0.14905763708553D+08
      Y(27) =   0.10005678719131D+08
      Y(28) =   0.18567530620184D+10
      Y(29) =   0.22722061569502D+07
      Y(30) =   0.51683758837725D+10
      Y(31) =   0.23175959524799D+07
      Y(32) =   0.30351396570102D+10
      Y(33) =   0.59289794793138D+07
      Y(34) =   0.25152972906931D+10
      Y(35) =   0.45696782763907D+07
      Y(36) =   0.35165924336687D+01
      Y(37) =   0.40635893499335D+05
      Y(38) =   0.38043227330586D-36
      Y(39) =   0.55775907287099D+09
      Y(40) =   0.76845966195032D+09
      Y(41) =   0.53384917779330D+10
      Y(42) =   0.60486040253666D+12
      Y(43) =   0.42222372170072D+08
      Y(44) =   0.36818527669954D+10
      Y(45) =   0.58357315681736D+07
      Y(46) =   0.14565071591335D+10
      Y(47) =   0.88628686841955D+10
      Y(48) =   0.93387533127053D+09
      Y(49) =   0.35957535867357D+10
      Y(50) =   0.46224275841192D+09
      Y(51) =   0.23959100649045D+09
      Y(52) =   0.60809365864182D+10
      Y(53) =   0.77041082974937D+09
      Y(54) =   0.33862208607592D+10
      Y(55) =   0.18382643870322D+10
      Y(56) =   0.19732523461142D+10
      Y(57) =   0.13134039409286D+07
      Y(58) =   0.13444717645883D+10
      Y(59) =   0.25684831172265D+07
      Y(60) =   0.11318545919859D+10
      Y(61) =   0.45915379379973D+07
      Y(62) =   0.86093822830561D+09
      Y(63) =   0.35095609713139D+09
      Y(64) =   0.66823509701060D+09
      Y(65) =   0.32073564029578D+09
      Y(66) =   0.11449876348727D+09
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE EMEPCF (TIME, RC, DJ, H2O)
      INTEGER NSPEC, NRC, NDJ
      DOUBLE PRECISION TIME, HMIX
      PARAMETER (NSPEC=66, NRC=266, NDJ=16)
      PARAMETER (HMIX=1.2D5)
      DOUBLE PRECISION RC(NRC), DJ(NDJ), H2O
C
C Compute time-dependent EMEP coefficients
C   RC: reaction coefficients
C   DJ: dissociation rate coefficient
C   H2O water vapour concentrations
C
C A and B: DJ=A*exp(-B*SEC)
C    SEC = 1/cos(THETA) where THETA is solar zenith angle
C T temperature in K
C
      DOUBLE PRECISION A(NDJ), B(NDJ), SEC, T

      INTEGER ITIMEH, I24HRS, I
      DOUBLE PRECISION TIMEH, TIMEOD, PI, XLHA, FI, DEKL, XQ, RH, XZ
      DOUBLE PRECISION M, O2, XN2, DELTA

      DATA A/1.23D-3,2.00D-4,1.45D-2,2.20D-5,3.00D-6
     -      ,5.40D-5,6.65D-5,1.35D-5,2.43D-5,5.40D-4
     -      ,2.16D-4,5.40D-5,3.53D-2,8.94D-2,3.32D-5
     -      ,2.27D-5/
      DATA B/   0.60D0,   1.40D0,   0.40D0,   0.75D0,   1.25D0
     -      ,   0.79D0,   0.60D0,   0.94D0,   0.88D0,   0.79D0
     -      ,   0.79D0,   0.79D0,  0.081D0,  0.059D0,   0.57D0
     -      ,   0.62D0/

      TIMEH=TIME/3600.D0
      ITIMEH=int(TIMEH)
      I24HRS=ITIMEH/24+1
      TIMEOD=TIMEH-(I24HRS-1)*24.D0
C
C Meteorology
C
      PI=4.D0*ATAN(1.0D0)
C   XLHA local hour angle
      XLHA=(1.D0+TIMEOD*3600.D0/4.32D4)*PI
C   FI (Norwegian PHI!) latitude, dekl solar declination
C   here latitude = 50 deg. N
      FI=50.D0*PI/180.D0
      DEKL=23.5D0*PI/180.D0
      SEC=1.D0/(COS(XLHA)*COS(FI)*COS(DEKL)+SIN(FI)*SIN(DEKL))
C   def of temperature variation
C     XP=8.7D-5*TIMEOD*3600.D0-2.83
C     T=8.3D0*SIN(XP)+289.86
C   for simplicity
      T=298.D0
C
C   def of water vapor concentration
      XQ=-7.93D-5*TIMEOD*3600.D0+2.43D0
      RH=23.D0*SIN(XQ)+66.5D0
C
      XZ=(597.3D0-0.57D0*(T-273.16D0))*18.D0/1.986D0*
     1 (1.D0/T-1.D0/273.16D0)
      H2O=6.1078D0*EXP(-XZ)*10.D0/(1.38D-16*T)*RH
C
C Calculate  values of photolysis rates DJ(1..16), based
C   upon RGD A & B coefficients and correction factors from HOUGH (1988)
C
      IF (TIMEOD .LT. 4.0D0 .OR. TIMEOD .GE. 20.D0) THEN
C   in the dark:
         DO 100 I = 1,NDJ
            DJ(I)=1.D-40
 100     CONTINUE
      ELSE
C   daytime:
         DO 110 I = 1,NDJ
            DJ(I)=A(I) * EXP( -B(I) * SEC )
            IF( DJ(I) .LT. 0.0D0 ) STOP 'DJ'
 110     CONTINUE
      ENDIF
C
C Set up chemical reaction rate coefficients:
C     16/6/92: inclusion of M, N2, O2 values in rate-constants
C     reaction rate coefficient definition. units: 2-body reactions
C     cm**3/(molecule x s), unimolecular 1.D0/s, 3-body
C     cm**6/(molecule**2 x s)
C
      M  = 2.55D19
      O2 = 5.2D18
      XN2= 1.99D19
C
      DO 120 I = 1,NRC
         RC(I)=0.D0
 120  CONTINUE
c
c
c..A92, assuming 80% N2, 20% O2
      rc(1)  =5.7d-34*(t/300.0D0)**(-2.8D0) * m * o2
c..unchanged, as source of O+NO+O2 reaction rate unknown.
      rc(5)  =9.6d-32*(t/300.0D0)**(-1.6D0) * m
c..estimate from Demore(1990) (=A92) O(d)+N2 and DeMore O(D)+O2
      rc(7)  =2.0d-11*exp(100.D0/t) * m
c.. A92    >>>>>
      rc(11) =1.8d-12*exp(-1370.D0/t)
      rc(12) =1.2d-13*exp(-2450.D0/t)
      rc(13) =1.9d-12*exp(-1000.D0/t)
      rc(14) =1.4d-14*exp(-600.D0/t)
      rc(15) =1.8d-11*exp(+110.D0/t)
      rc(17) =3.7d-12*exp(240.0D0/t)
c.. <<<<<<    A92
c..M.J (from Wayne et al.)
      rc(19) =7.2d-14*exp(-1414.D0/t)
cfix  rc(19) =7.2d-13*exp(-1414.D0/t)
c..M.J suggests that rc(27) be omitted:
c     rc(27) =8.5d-13*exp(-2450.D0/t)
c..  My change to get similar to A92 troe.
      rc(29) =7.1d14*exp(-11080.D0/t)
c..A92
      rc(30) =4.8d-11*exp(+250.D0/t)
c..A92, De More,1990 .. no change : oh + h2o2
      rc(31) =2.9d-12*exp(-160.D0/t)
c..A92 : oh + h2
      rc(33) =7.7d-12*exp(-2100.D0/t)
c..My, similar to DeMore et al complex : oh+hno3
      rc(35) =1.0d-14*exp(785.0D0/t)
c.. Mike Jenkin`s suggestion for ho2 + ho2 reactions: (from DeMore et al.)
      rc(36) =2.3d-13*exp(600.D0/t)
      rc(36) = rc(36) + m * 1.7d-33*exp(1000.D0/t)
      rc(36) = rc(36) *
     &         (1.D0 + 1.4d-21 * h2o *exp(2200.D0/t))

c..A92
      rc(59) =3.9d-12*exp(-1885.D0/t)
c A92 : ch3o2 + no
      rc(60) =4.2d-12*exp(180.D0/t)
c A92 + A90 assumption that ka = 0.5D0 * k
      rc(61) =5.5d-14*exp(365.D0/t)
c A92 + A90 assumption that kb = 0.5D0 * k
      rc(62) =5.5d-14*exp(365.D0/t)
c A92
      rc(63) =3.3d-12*exp(-380.D0/t)
c A92 : ho2 + ch3o2
      rc(65) =3.8d-13*exp(780.D0/t)
c A92 new: ch3ooh + oh -> hcho + oh
      rc(67) =1.0d-12*exp(190.D0/t)
c A92 new: ch3ooh + oh -> ch3o2
      rc(68) =1.9d-12*exp(190.D0/t)
c.. A92
      rc(71) =7.8d-12*exp(-1020.D0/t)
c A92 new: c2h5o2 + ho2 -> c2h5ooh (r2ooh)
      rc(74) =6.5d-13*exp(650.D0/t)
c.. A92
      rc(75) =5.6d-12*exp(310.D0/t)
c TOR90 assumption w.r.t. rc(67) : c2h5ooh + oh -> ch3cho + oh
c     rc(76) = 5.8 * rc(67)
      rc(76) = 5.8d-12*exp(190.D0/t)
cA92 : approximation to troe expression
      rc(78) =1.34d16*exp(-13330.D0/t)
c additional reactions :-
c A92 : ho2 + ch3coo2 -> rco3h
      rc(88) =1.3d-13*exp(1040.D0/t)
c A92 : ho2 + ch3coo2 -> rco2h + o3
      rc(89) =3.0d-13*exp(1040.D0/t)
c.. A92
      rc(94) =2.8d-12*exp(530.D0/t)
c.. D & J, gives results very close to A92
      rc(81) =1.64d-11*exp(-559.D0/t)

c.. A90
      rc(83) =rc(60)
      rc(105)=rc(60)
cA90
      rc(110)=rc(60)
cA90/PS  isnir + no
      rc(162)=rc(60)
cA90/PS  isono3 + no
      rc(164)=rc(60)

c.. From RGD, but very similar to A92 Troe
      rc(109)=1.66d-12*exp(474.D0/t)
cA92
      rc(112)=1.2d-14*exp(-2630.D0/t)
cA92
      rc(123)=6.5d-15*exp(-1880.D0/t)
cA90
      rc(126) = rc(60)
c..A90
      rc(220) = rc(60)
c..A90
      rc(236) = rc(60)
cooooooooooooo   natural voc reactions 00000000000000000
c..A90  isoprene + o3 -> products
      rc(150) = 12.3d-15*exp(-2013.D0/t)
c..A90  isoprene + oh -> isro2
      rc(151) = 2.54d-11*exp(410.D0/t)
cA90  isoprene-RO2 + no
      rc(152) = rc(60)
cA90  methylvinylketone (mvk) + oh
      rc(153) = 4.13d-12*exp(452.D0/t)
cA90  mvko2 + no
      rc(154) = rc(60)
cA90  macr + oh
      rc(158) = 1.86d-11*exp(175.D0/t)
cA90  ch2cch3 + no
      rc(159) = rc(60)
cA90  mvk + o3
      rc(160) = 4.32d-15*exp(-2016.D0/t)
c
c     rc(255) = 9.9d-16*exp(-731/t)
c.......................................................................
c
coooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c     aerosol formation and depositions ..x
c     parameterization of heterogeneous loss in 1./s for humidities
c     less than 90%.  applied to hno3, h2o2, h2so4, ch3o2h.
      rc(43)=5.d-6
      if(rh.gt.0.90D0) rc(43)=1.d-4
      rc(44) =rc(43)
      rc(45) =rc(43)

c.. A92
      rc(8)  =2.2d-10
c.. My (?) , similar to A92 troe expression.
      rc(20) =1.4d-12
c.. My (?) , similar to A92 troe expression.
      rc(21) =1.4d-11
c.. RGD Note (=DeMore?)
      rc(26) =4.1d-16
c.. RGD
      rc(39) =1.35d-12
c.. RGD
      rc(40) =4.0d-17
c.. A92, assuming all products are ch3cho
      rc(64) =3.2d-12
c.. A92, with temperature dependance neglected.
      rc(66) =9.6d-12
c.. A92 >>>>>>>>
      rc(69) =5.8d-16
      rc(70) =2.4d-13
      rc(72) =8.9d-12
c.. <<<<<<<<<< A92
c.. A92 : approximation to troe expression
      rc(77) =1.0d-11
c.. A92 : ch3coo2 + no
      rc(79) =2.0d-11
c.. A92 : sum of two pathways
      rc(80) =1.1d-11
cmj   rc(84) =2.5d-14
cya   kho2ro2 : estimate of rate of ho2 + higher RO2, suggested by Yvonne
      rc(85) =1.0d-11
c..A90, ignoring slight temperature dependance
      rc(86) =1.15d-12
c..MJ suggestion.. combine oh + secc4h9o2, meko2 rates   =>
c     rc(87)= rc(68) + rc(86), approx. at 298
      rc(87)= 4.8d-12
cmj..new A92 : oh + ch3co2h -> ch3o2
      rc(90) =8.0d-13
c.. Approximates to A92 Troe ...
      rc(125)=2.86d-11
c.. rate for ch2chr + oh, = k68+k125 (propene), Yv.
      rc(146) = 3.2d-11
c.. rate for isono3h + oh, = k156 (isro2h), Yv.
      rc(147) = 2.0d-11
c.. MY GUESS rate of oh + mvko2h
      rc(148) = 2.2d-11
c.. MY GUESS rate of oh + other biogenic ooh
      rc(149) = 3.7d-11
cya   kho2ro2 : estimate of rate of ho2 + higher RO2, suggested by Yvonne
cA90  isro2 + ho2
      rc(155) =rc(85)
cPS   isro2h + oh
      rc(156) =2.0d-11
cPS   isro2h + o3
      rc(157) =8.0d-18
cPS   isni + oh
      rc(161) =3.35d-11
cPS   isopre + no3
      rc(163) =7.8d-13
cZZ   rc(163) =7.8d-16
c.. Unchanged, also in IVL scheme
      rc(219)=2.0d-11
c..A92
      rc(221)=1.1d-11
c..A92
      rc(222)=1.7d-11
c..MJ suggestion.. combine oh + malo2h rates   =>
c     rc(223)= rc(68) + rc(219)
      rc(223)= 2.4d-11
c..A90
      rc(234)=1.37d-11
c..MJ suggestion.. combine rc(68) with rc(234) =>
c     rc(235)= rc(68) + rc(234)
      rc(235)= 1.7d-11
c
c..............................................
c         deposition loss rate coefficients vd/hmix, vd in cm/s.
c         hno3     calculated     rc(46)
c         so2      0.5            rc(47)
c         h2o2     0.5            rc(47)
c         no2      0.2            rc(48)
c         o3       0.5            rc(49)
c         pan      0.2            rc(50)
c         h2so4    0.1            rc(51)
c... simple approx. for now - reduce all vg by 4 at night to allow
c    for surface inversion......
          delta=1.0D0
          if(timeod.ge.20.D0.or.timeod.lt.4.D0) delta=0.25D0
c         if(timeod.ge.20.D0.or.timeod.le.4.D0) delta=0.25D0
c         if (night) delta=0.25D0
          rc(46) =2.0D0 * delta /hmix
          rc(47) =0.5D0 * delta /hmix
          rc(48) =0.2D0 * delta /hmix
          rc(49) =0.5D0 * delta /hmix
          rc(50) =0.2D0 * delta /hmix
          rc(51) =0.1D0 * delta /hmix
c. dep. of organic peroxides = 0.5 cms-1
          rc(52) = 0.5D0 *delta /hmix
c. dep. of ketones, RCHO  = 0.3 cms-1
          rc(53) = 0.3D0 *delta /hmix
c
      RETURN
      END

 
