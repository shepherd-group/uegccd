! Change module name to dRPA
! Change the entering subroutine and figure out if I want to change anything
! else, perhaps just set things to private? 
! TODO - TNM: Add in the routines to print the structure factor for dRPA (Copy
! from CCD.f90 routine) 


Module dRPA

Use Precision
Use Constants
Use Pointers ! Has the abstract interface to ERI :)

Implicit None

Real (Kind=pr) :: DenomFactor = 2.0_pr

private
public :: DrvdRPA

Contains   

      Subroutine DrvdRPA(HEGData, UEGInfo, T2)

      Use Types, only: HEGDataType, UEGInfoType, CCScrDataType, DIISCCDataType
      Use HEG, only: FindTol
      Use CCScr, only: SetUpCCScr, ShutDownCCScr 
      Use DIISCC

      Type (HEGDataType), Intent(InOut) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo
      Type (CCScrDataType) :: CCScrData
      Type (DIISCCDataType) :: DIISCCData

      Real (Kind=pr), Intent(InOut) :: T2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

      Integer,        Parameter     :: MaxIter = 1000
      Real (Kind=pr)                :: TolMax = 1.0E-8_pr
      Real (Kind=pr)                :: FailRatio = 100
      Logical :: Fail
      Integer :: NIter
      Real (Kind=pr) :: dT, MaxRes

!========================================================!
!  This routine does the entire dRPA calculation for us. !
!  We take in the MBPT(2) results as input and write     !
!  them out as the first dRPA iteration even though we   !
!  don't actually repeat them.                           !
!                                                        !
!  Only the rings terms from the  CCD equations are      !
!  being accounted for to give rCCD, which is equivlaent !
!  to dRPA+SOSEX, and dRPA is obatined by subtracting    !
!  off the SOSEX terms                                   !
!                                                        !
!  We use DIIS to converge this better.                  !
!========================================================!

      Call FindTol(HEGData%rS,TolMax,FailRatio,DenomFactor)

      Write(6,*) 'Doing RPA...'
! Allocate space for everything.
      Call SetUpCCScr(UEGInfo%NOcc,UEGInfo%NOcc+1,UEGInfo%NAO,CCScrData)
      Call SetUpDIISCC(DIISCCData,UEGInfo%NOcc,UEGInfo%NOcc+1,UEGInfo%NAO,T2,T2)


! Initialize variables and write the MP2 energy out again
      NIter = 1
      dT = MaxVal(Abs(T2))
      Call CCEnergy(HEGData,UEGInfo,T2, Zero)
      Open(7,File='Output',Position='Append')
      Write(7,1010)
      Write(7,1020)
      If(UEGInfo%DoRing)   Write(7,1100)
      If(UEGInfo%DoMosaic) Write(7,1130)
      Write(7,1020)
      Write(7,1030)
      Write(7,1040) HEGData%ECorr,NIter,dT


! Start iterating!
      Do While(dT >= TolMax .and. NIter < MaxIter)
! Given OldT2 and R2, update T2 and R2
        Call GetTFromDIIS(T2,DIISCCData%OldT2abab,DIISCCData%R2abab, &
                          UEGInfo%NOcc,UEGInfo%NAO,DIISCCData%NDIIS,NIter)
! Given T2 from DIIS, calculate G2.
! Calculate the residual G[T]-HT
        Call TermsWrapper(HEGData,UEGInfo,CCScrData,T2,DIISCCData%R2abab,DIISCCData%NDIIS)
! Solve the CC equations, update OldT2, and get the energy
        Call SolveCC(UEGInfo,HEGData,T2,CCScrData%G2abab)
        Call CnvrgCC(T2,DIISCCData%OldT2abab,dT,UEGInfo%NOcc,UEGInfo%NAO,NIter,DIISCCData%NDIIS)
        Call CCEnergy(HEGData,UEGInfo,T2, Zero)
        Write(7,1040) HEGData%ECorr, NIter, dT
        Write(6,1040) HEGData%ECorr, NIter, dT
      End Do


! Check the resisdual
      Call TermsWrapper(HEGData,UEGInfo,CCScrData,T2,DIISCCData%R2abab,DIISCCData%NDIIS)
      MaxRes =  MaxVal(Abs(DIISCCData%R2abab(:,:,:,DIISCCData%NDIIS)))
      FailRatio=10000
      Fail = MaxRes > FailRatio*TolMax


! Finish the ouptut and shut down
      Write(7,1020)
      If(Fail) Then
        Write(7,2000)
      Else
        Write(7,1050) NIter
        Write(7,1060) HEGData%ECorr
        Write(7,1070) HEGData%ECorr+HEGData%EHF
        Call CCEnergy(HEGData,UEGInfo,T2, One)
        Write(7,1065) HEGData%ECorr
        Write(7,1075) HEGData%ECorr+HEGData%EHF
      End If
      Write(7,1090) MaxRes
      Write(7,1000)
      Close(7)
      Call ShutDownCCScr(CCScrData)
      Call ShutDownDIISCC(DIISCCData)


1000  Format(14x,'**************************************************')
1010  Format(14X,'*              dRPA summary follows              *')
1020  Format(14x,'*------------------------------------------------*')
1030  Format(14X,'* dRPA Energy      Iteration    Biggest change   *')
1040  Format(14X,'* ',F15.10,2x,I5,7X,F14.10,4x,'*')
1050  Format(14x,'* dRPA has converged in ',I3,' iterations',11x,'*')
1060  Format(14x,'* Ec(dRPA) is ',9x,F18.12,' a.u.   *')
1065  Format(14x,'* Ec(SOSEX) is ',8x,F18.12,' a.u.   *')
1070  Format(14x,'* Final dRPA Energy is ',F18.12,' a.u.   *')
1075  Format(14x,'* Final SOSEX Energy is',F18.12,' a.u.   *')
1090  Format(14x,'* Max CCD Residual is  ',4x,ES18.12,'    *')
1100  Format(14x,'* Including ring diagrams...                     *')
1110  Format(14x,'* Including crossed ring diagrams...             *')
1120  Format(14x,'* Including ladder diagrams...                   *')
1130  Format(14x,'* Including mosaic diagrams...                   *')
2000  Format(14x,'* Final dRPA Energy is   DID NOT CONVERGE        *')

      Return

      Contains

        Subroutine TermsWrapper(HEGData,UEGInfo,CCScrData,T2,R2abab,NDIIS)

          ! Term function calls, now with half the code.
          ! A wise being once said: "it's the little things in life" :)

          Use Types, only: HEGDataType, UEGInfoType, CCScrDataType

          Type (HEGDataType), Intent(InOut) :: HEGData
          Type (UEGInfoType), Intent(In) :: UEGInfo
          Type (CCScrDataType), Intent(InOut) :: CCScrData 

          Real (Kind=pr), Intent(InOut) :: T2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
          Real (Kind=pr), Intent(InOut) :: R2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO,NDIIS)

          Integer, Intent(In) :: NDIIS

          Associate(G2abab => CCScrData%G2abab, &
                    Joo => CCScrData%Joo,       &
                    Jvv => CCScrData%Jvv)

            Call GetG2Drivers(HEGData,UEGInfo,G2abab)
            If (UEGInfo%DoMosaic) Call GetG2Mosaics(HEGData,UEGInfo,G2abab,Joo,Jvv,T2)
            If (UEGInfo%DoRing) Call GetG2Rings(HEGData,UEGInfo,G2abab,T2)

            Call GetRes(HEGData,UEGInfo,T2,G2abab,R2abab,NDIIS)

          End Associate

        End Subroutine TermsWrapper

      End Subroutine DrvdRPA






      Subroutine CCEnergy(HEGData,UEGInfo,T2,CSosex)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType
      
      Type (HEGDataType), Intent(InOut) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(In) :: T2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(In) :: CSosex

      Real (Kind=pr) :: V_ijab, V_ijba, dE
      Integer        :: I, J, A, B

!===========================================!
!  This routine calculates the CCD energy.  !
!===========================================!

! THis needs to be fixed for momentum symmetry
      HEGData%ECorr = Zero
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A) 
        If(B <= UEGInfo%NOcc) Cycle
        V_ijab = ERI(HEGData,UEGInfo,A,B,I,J)
        V_ijba = ERI(HEGData,UEGInfo,B,A,I,J)
        dE = V_ijab*T2(I,J,A)*2  &   ! 2 from T2aaaa + T2abab = 2 T2
           - V_ijba*T2(I,J,A)*CSosex
        HEGData%ECorr = HEGData%ECorr + dE ! 1/4 * 2 from spin * 2 from dRPA
      End Do
      End Do
      End Do

      Return
      End Subroutine CCEnergy






      Subroutine SolveCC(UEGInfo,HEGData,T2,G2)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(Out)    :: T2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(In)     :: G2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Real (Kind=pr) :: Denom
      Integer        :: I, J, A, B

!=============================!
!  Solves H.T = G.  Trivial.  !
!=============================!

      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle
        Denom = HEGData%Eigen(I) + HEGData%Eigen(J) - HEGData%Eigen(A) - HEGData%Eigen(B)
!       T2(I,J,A) = G2(I,J,A)/Denom
        T2(I,J,A) = (G2(I,J,A)/Denom - (One-DenomFactor)*T2(I,J,A))/DenomFactor
      End Do
      End Do
      End Do

      Return
      End Subroutine SolveCC






      Subroutine GetRes(HEGData,UEGInfo,T2,G2,R2,NDIIS)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(In)    :: T2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                       G2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(InOut) :: R2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO,NDIIS)

      Integer,        Intent(In)    :: NDIIS

! Local variables
      Real (Kind=pr) :: Rijab, Denom
      Integer :: I, J, A, B

!===============================!
!  Computes the residual HT-G.  !
!===============================!

      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle
        Denom = HEGData%Eigen(I) + HEGData%Eigen(J) - HEGData%Eigen(A) - HEGData%Eigen(B)
        Rijab = -G2(I,J,A) + Denom*T2(I,J,A)
        R2(I,J,A,NDIIS) = Rijab
      End Do
      End Do
      End Do

      Return
      End Subroutine GetRes






      Subroutine CnvrgCC(T2,OldT2,dT,NOcc,NAO,NIter,NDIIS)

      Integer,        Intent(In)    :: NOcc, NAO, NDIIS
      Integer,        Intent(InOut) :: NIter
      Real (Kind=pr), Intent(Out)   :: dT
      Real (Kind=pr), Intent(In)    :: T2(NOcc,NOcc,NOcc+1:NAO)
      Real (Kind=pr), Intent(InOut) :: OldT2(NOcc,NOcc,NOcc+1:NAO,NDIIS)

      Integer, Parameter :: MaxIter = 10000

!======================================!
!  Update the DIIS amplitudes OldT2*.  !
!======================================!

! Update the iteration count
      NIter = NIter + 1

! Move the DIIS amplitudes
      OldT2(:,:,:,1:NDIIS-1) = OldT2(:,:,:,2:NDIIS)

! Get the change in T
      dT = MaxVal(Abs(T2 - OldT2(:,:,:,NDIIS)))

! Update the most recent DIIS amplitudes
      OldT2(:,:,:,NDIIS) = T2
      If(NIter == MaxIter) Stop 'CCD equations did not converge'

      Return
      End Subroutine CnvrgCC






      Subroutine GetG2Drivers(HEGData,UEGInfo,G2)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(Out) :: G2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Integer :: I, J, A, B
      Real (Kind=pr) :: V_ijab

!==================================!
!  Build the driving terms first.  !
!==================================!

      G2 = Zero
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle
        V_ijab = ERI(HEGData,UEGInfo,I,J,A,B,DummyFlag=UEGInfo%IRangeDriverDirect)
        If(V_ijab < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
        G2(I,J,A) =  V_ijab
      End Do
      End Do
      End Do

      Return
      End Subroutine GetG2Drivers


      Subroutine GetG2Mosaics(HEGData,UEGInfo,G2,Joo,Jvv,T2)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(In)    :: T2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(InOut) :: G2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(InOut) :: Joo(UEGInfo%NOcc),                                            &
                                       Jvv(UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Integer :: I, J, A, B, K, L, C, D
      Real (Kind=pr) :: V_cdkl

!===========================================================!
!  Handle the mosaic terms.  These are                      !
!   dG_ij^ab = -1/2 vbar^{kl}_{cd} t_{kl}^{ad} t_{ij}^{cb}  !
!              -1/2 vbar^{kl}_{cd} t_{kl}^{cb} t_{ij}^{ad}  !
!              -1/2 vbar^{kl}_{cd} t_{il}^{cd} t_{kj}^{ab}  !
!              -1/2 vbar^{kl}_{cd} t_{kj}^{cd} t_{il}^{ab}  !
!  They're particularly simple - we only have to build the  !
!  diagonal intermediates in this case.                     !
!                                                           !
!  Let's see that explicitly.  Consider first Jvv.          !
!  We would have                                            !
!     J^a_c = t_kl^ad v^kl_cd                               !
!  and then build                                           !
!     dG_ij^ab = t_ij^cb J^a_c                              !
!  But we know that i+j+a+b = 0, and we need i+j+c+b = 0,   !
!  which means that a=c => we need only the diagonal.       !
!-----------------------------------------------------------!
!  Start by building intermediates.                         !
!  It sure looks like the Joo and Jvv are identical!        !
!===========================================================!

      Joo = Zero
      Jvv = Zero
      Do K = 1,UEGInfo%NOcc
      Do L = 1,UEGInfo%NOcc
      Do C = UEGInfo%NOcc+1,UEGInfo%NAO
        D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,K,L,C)
        If(D <= UEGInfo%NOcc) Cycle
        V_cdkl = ERI(HEGData,UEGInfo,C,D,K,L,DummyFlag=UEGInfo%IRangeMosaic)
        If(V_cdkl < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
        Jvv(C) = Jvv(C) + 2*V_cdkl*T2(K,L,C)   ! 2 from T2aaaa + T2abab = 2 T2
        Joo(K) = Joo(K) + 2*V_cdkl*T2(K,L,C)
      End Do
      End Do
      End Do

! The intermediates are done, so contract with T2
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle
        G2(I,J,A) = G2(I,J,A) - F12*(Jvv(A)*T2(I,J,A) + Jvv(B)*T2(I,J,A))
        G2(I,J,A) = G2(I,J,A) - F12*(Joo(I)*T2(I,J,A) + Joo(J)*T2(I,J,A))
      End Do
      End Do
      End Do

      Return
      End Subroutine GetG2Mosaics






      Subroutine GetG2Rings(HEGData,UEGInfo,G2,T2)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(In)    :: T2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(InOut) :: G2(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Integer :: I, J, A, B, K, L, C, D, MaxDummyFlag
      Real (Kind=pr) :: V_cjkb, V_cika, V_cdkl, J_idal

!===========================================================!
!  Build the ring terms.  These are                         !
!   dG_ij^ab = t_{ik}^{ac} vbar^{kb}_{cj}                   !
!            + t_{kj}^{cb} vbar^{ak}_{ic}                   !
!            + t_{ik}^{ac} t_{lj}^{db} vbar^{kl}_{cd}       !
!===========================================================!

! These are the linear terms
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle

! Do the t_ik^ac V_jc^bk contraction
        Do K = 1,UEGInfo%NOcc
          C = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,K,A)
          If(C <= UEGInfo%NOcc) Cycle
          MaxDummyFlag = Max(UEGInfo%IRangeRing, &
                             UEGInfo%IRangeDirectRings, &
                             UEGInfo%IRangeLinRings)
          V_cjkb = ERI(HEGData,UEGInfo,C,J,K,B,DummyFlag=MaxDummyFlag)
          If(V_cjkb < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2(I,J,A) = G2(I,J,A) + 2*V_cjkb*T2(I,K,A)   ! 2 from T2aaaa + T2abab = 2 T2
        End Do

! Do the t_jk^bc V_ic^ak contraction
        Do K = 1,UEGInfo%NOcc
          C = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,J,K,B)
          If(C <= UEGInfo%NOcc) Cycle
          MaxDummyFlag = Max(UEGInfo%IRangeRing, &
                             UEGInfo%IRangeDirectRings, &
                             UEGInfo%IRangeLinRings)
          V_cika = ERI(HEGData,UEGInfo,C,I,K,A,DummyFlag=MaxDummyFlag)
          If(V_cika < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2(I,J,A) = G2(I,J,A) + 2*V_cika*T2(J,K,B)   ! 2 from T2aaaa + T2abab = 2 T2
        End Do
      End Do
      End Do
      End Do


! These are the quadratic terms.  They require intermediates.
      Do I = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        Do L = 1,UEGInfo%NOcc
          J_idal = Zero
          Do K = 1,UEGInfo%NOcc
            C = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,K,A)
            If(C <= UEGInfo%NOcc) Cycle
            D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,K,L,C)
            If(D <= UEGInfo%NOcc) Cycle
            MaxDummyFlag = Max(UEGInfo%IRangeRing, &
                               UEGInfo%IRangeDirectRings, &
                               UEGInfo%IRangeQuadRings)
            V_cdkl = ERI(HEGData,UEGInfo,C,D,K,L,DummyFlag=MaxDummyFlag)
            If(V_cdkl < -80.0_pr) Then
              Print *, "Disallowed excitation not trapped!"
              Cycle
            End If
            J_idal = J_idal + 2*T2(I,K,A)*V_cdkl     ! 2 from T2aaaa + T2abab = 2 T2
          End Do

! The intermediates are done, so contract with T2
! In principle, we've got the right D here.
          Do J = 1,UEGInfo%NOcc
            B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
            If(B <= UEGInfo%NOcc) Cycle
            D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,J,L,B)
            If(D <= UEGInfo%NOcc) Cycle
            G2(I,J,A) = G2(I,J,A) + 2*J_idal*T2(L,J,D)   ! 2 from T2aaaa + T2abab = 2 T2
          End Do
        End Do
      End Do
      End Do

      Return
      End Subroutine GetG2Rings

   End Module dRPA

