!TODO: TNM - Need to troubleshoot incorrect energy on second iteration by printing the vectors and check against test code 

Module CCD

Use Precision
Use Constants
Use Pointers ! Has the abstract interface to ERI :)

Implicit None

Real (Kind=pr) :: DenomFactor = 2.0_pr

private
public :: DrvCCD

Contains

      Subroutine DrvCCD(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)

      Use Types, only: HEGDataType, UEGInfoType, CCScrDataType, DIISCCDataType
      Use HEG, only: FindTol
      Use CCScr, only: SetUpCCScr, ShutDownCCScr 
      Use DIISCC

      Type (HEGDataType), Intent(InOut) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo
      Type (CCScrDataType) :: CCScrData
      Type (DIISCCDataType) :: DIISCCData

      Real (Kind=pr), Intent(InOut) :: T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                       T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                       T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

      Real (Kind=pr)                :: TolMax = 1.0E-8_pr
      Real (Kind=pr)                :: FailRatio = 100
      Integer,        Parameter     :: MaxIter = 1000
      Logical :: Fail
      Integer :: NIter
      Real (Kind=pr) :: dT, MaxRes

!========================================================!
!  This routine does the entire CCD calculation for us.  !
!  We take in the MBPT(2) results as input and write     !
!  them out as the first CCD iteration even though we    !
!  don't actually repeat them.                           !
!                                                        !
!  Certain terms of the CCD equations can be omitted in  !
!  what I think is an obvious way.                       !
!                                                        !
!  We use DIIS to converge this better.                  !
!========================================================!

      Call FindTol(HEGData%rS,TolMax,FailRatio,DenomFactor)

      Write(6,*) 'Doing CCD...'
! Allocate space for everything.
      Call SetUpCCScr(UEGInfo%NOcc,UEGInfo%NOcc+1,UEGInfo%NAO,CCScrData)
      Call SetUpDIISCC(DIISCCData,UEGInfo%NOcc,UEGInfo%NOcc+1,UEGInfo%NAO,T2abab,T2abba)


! Initialize variables and write the MP2 energy out again
      NIter = 1
      dT = Max(MaxVal(Abs(T2abab)),MaxVal(Abs(T2abba)))
      Call CCEnergy(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)
      if (UEGInfo%DoSFCalcMP2) Call CCEnergyAndSF(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)

      Open(7,File='Output',Position='Append')
      Write(7,1010)
      Write(7,1020)
      If(UEGInfo%DoRing)   Write(7,1100)
      If(UEGInfo%DoXRing)  Write(7,1110)
      If(UEGInfo%DoLadder) Write(7,1120)
      If(UEGInfo%DoMosaic) Write(7,1130)
      Write(7,1020)
      Write(7,1030)
      Write(7,1040) HEGData%ECorr,NIter,dT


! Start iterating!
      Do While(dT >= TolMax .and. NIter < MaxIter)
! Given OldT2 and R2, update T2 and R2
        Call GetTFromDIIS(T2abab,DIISCCData%OldT2abab,DIISCCData%R2abab, &
                          UEGInfo%NOcc,UEGInfo%NAO,DIISCCData%NDIIS,NIter)
        Call GetTFromDIIS(T2abba,DIISCCData%OldT2abba,DIISCCData%R2abba, &
                          UEGInfo%NOcc,UEGInfo%NAO,DIISCCData%NDIIS,NIter)
        T2aaaa = T2abab + T2abba
! Given T2 from DIIS, calculate G2 and calculate the residual G[T]-HT
        Call TermsWrapper(HEGData,UEGInfo,CCScrData,T2aaaa,T2abab,T2abba, &
                          DIISCCData%R2abab,DIISCCData%R2abba,DIISCCData%NDIIS)
! Solve the CC equations, update OldT2, and get the energy
        Call SolveCC(HEGData,UEGInfo,T2abab,CCScrData%G2abab)
        Call SolveCC(HEGData,UEGInfo,T2abba,CCScrData%G2abba)
        T2aaaa = T2abab + T2abba
        Call CnvrgCC(T2abab,DIISCCData%OldT2abab,T2abba,DIISCCData%OldT2abba, &
                     dT,UEGInfo%NOcc,UEGInfo%NAO,NIter,DIISCCData%NDIIS)
        Call CCEnergy(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)
        Write(7,1040) HEGData%ECorr, NIter, dT
        Write(6,1040) HEGData%ECorr, NIter, dT
      End Do


! One final check of the resisdual
      Call TermsWrapper(HEGData,UEGInfo,CCScrData,T2aaaa,T2abab,T2abba, &
                        DIISCCData%R2abab,DIISCCData%R2abba,DIISCCData%NDIIS)
      MaxRes =  Max(MaxVal(Abs(DIISCCData%R2abab(:,:,:,DIISCCData%NDIIS))), &
                    MaxVal(Abs(DIISCCData%R2abba(:,:,:,DIISCCData%NDIIS))))
      Fail = MaxRes > FailRatio*TolMax

      if (UEGInfo%DoSFCalcCCD) Call CCEnergyAndSF(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)

! Finish the ouptut and shut down
      Write(7,1020)
      If(Fail) Then
        Write(7,2000)
      Else
        Write(7,1050) NIter
        Write(7,1060) HEGData%ECorr
        Write(7,1070) HEGData%ECorr+HEGData%EHF
      End If
      Write(7,1090) MaxRes
      Write(7,1000)
      Close(7)
      Call ShutDownCCScr(CCScrData)
      Call ShutDownDIISCC(DIISCCData)


1000  Format(14x,'**************************************************')
1010  Format(14X,'*               CCD summary follows              *')
1020  Format(14x,'*------------------------------------------------*')
1030  Format(14X,'*  CCD Energy      Iteration    Biggest change   *')
1040  Format(14X,'* ',F15.10,2x,I5,7X,F14.10,4x,'*')
1050  Format(14x,'* CCD has converged in ',I3,' iterations',12x,'*')
1060  Format(14x,'* Ec(CCD) is ',10x,F18.12,' a.u.   *')
1070  Format(14x,'* Final CCD Energy is  ',F18.12,' a.u.   *')
1090  Format(14x,'* Max CCD Residual is  ',4x,ES18.12,'    *')
1100  Format(14x,'* Including ring diagrams...                     *')
1110  Format(14x,'* Including crossed ring diagrams...             *')
1120  Format(14x,'* Including ladder diagrams...                   *')
1130  Format(14x,'* Including mosaic diagrams...                   *')
2000  Format(14x,'* Final CCD Energy is   DID NOT CONVERGE         *')

      Return

      Contains

        Subroutine TermsWrapper(HEGData,UEGInfo,CCScrData,T2aaaa,T2abab,T2abba,R2abab,R2abba,NDIIS)
          ! CK: consider starting parallel block here instead of each indv subroutine
          ! Term function calls, now with half the code.
          ! A wise being once said: "it's the little things in life" :)

          Use Types, only: HEGDataType, UEGInfoType, CCScrDataType

          Type (HEGDataType), Intent(InOut) :: HEGData
          Type (UEGInfoType), Intent(In) :: UEGInfo
          Type (CCScrDataType), Intent(InOut) :: CCScrData 

          Real (Kind=pr), Intent(InOut) :: T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                           T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                           T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
          Real (Kind=pr), Intent(InOut) :: R2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO,NDIIS), &
                                           R2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO,NDIIS)

          Integer, Intent(In) :: NDIIS

          Associate(G2abab => CCScrData%G2abab, &
                    G2abba => CCScrData%G2abba, &
                    Joo => CCScrData%Joo,       &
                    Jvv => CCScrData%Jvv)

            Call GetG2Drivers(HEGData,UEGInfo,G2abab,G2abba)
            If (UEGInfo%DoMosaic) Call GetG2Mosaics(HEGData,UEGInfo,G2abab,G2abba,Joo,Jvv,T2aaaa,T2abab,T2abba)
            If (UEGInfo%DoLadder) Call GetG2Ladders(HEGData,UEGInfo,G2abab,G2abba,T2aaaa,T2abab,T2abba)
            If (UEGInfo%DoRing) Call GetG2Rings(HEGData,UEGInfo,G2abab,G2abba,T2aaaa,T2abab,T2abba)
            If (UEGInfo%DoXRing) Call GetG2XRings(HEGData,UEGInfo,G2abab,G2abba,T2aaaa,T2abab,T2abba)

            Call GetRes(HEGData,UEGInfo,T2abab,G2abab,R2abab,NDIIS)
            Call GetRes(HEGData,UEGInfo,T2abba,G2abba,R2abba,NDIIS)

          End Associate

        End Subroutine TermsWrapper

      End Subroutine DrvCCD


      Subroutine CCEnergy(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(InOut) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(In) :: T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Real (Kind=pr) :: V_ijab, V_ijba, Eaaaa, Eabab, Eabba, dE
      Integer        :: I, J, A, B

!===========================================!
!  This routine calculates the CCD energy.  !
!===========================================!

! THis needs to be fixed for momentum symmetry
      HEGData%ECorr = Zero
      Eaaaa = Zero
      Eabab = Zero
      Eabba = Zero
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle
        V_ijab = ERI(HEGData,UEGInfo,A,B,I,J,DummyFlag=UEGInfo%IRangeEnergy)
        V_ijba = ERI(HEGData,UEGInfo,A,B,J,I,DummyFlag=UEGInfo%IRangeEnergy)
        If(Min(V_ijab,V_ijba) < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
!       Eaaaa = Eaaaa + T2aaaa(I,J,A)*(V_ijab - V_ijba)
!       Eabab = Eabab + T2abab(I,J,A)*V_ijab
!       Eabba = Eabba - T2abba(I,J,A)*V_ijba
        dE = V_ijab*(T2aaaa(I,J,A) + T2abab(I,J,A))    &
           - V_ijba*(T2aaaa(I,J,A) + T2abba(I,J,A))
        HEGData%ECorr = HEGData%ECorr + F12*dE
      End Do
      End Do
      End Do
!     Eaaaa = Eaaaa*F12
!     Eabab = Eabab*F12
!     Eabba = Eabba*F12
!     ECorr = ECorr + Eaaaa + Eabab + Eabba

      Return
      End Subroutine CCEnergy 


      Subroutine CCEnergyAndSF(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)

      Use HEG, Only: FindIndex, SendForAllocation, dK, dK2, dK2scaled
      Use Types, Only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(InOut) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(In) :: T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Real (Kind=pr) :: V_ijab, V_ijba, Eaaaa, Eabab, Eabba, dE
      Integer        :: I, J, A, B
      Logical :: inCore
      Real (pr), Allocatable :: SfSum(:,:,:), SfG(:,:,:), SfV(:,:,:)
      Real (pr), Allocatable :: SfDSum(:,:,:), SfXSum(:,:,:)
      Integer :: SfLength
      Integer :: this_qvec(3)
      Real (pr) :: this_q2

!===========================================!
!  This routine calculates the CCD energy.  !
!===========================================!
      write(80,*) "--------"
      write(80,*) "Starting a new structure factor"
      write(80,*) "--------"
      
      call SendForAllocation(HEGData,SfSum,SfLength)
      call SendForAllocation(HEGData,SfG,SfLength)
      call SendForAllocation(HEGData,SfV,SfLength)
      call SendForAllocation(HEGData,SfXSum,SfLength)
      call SendForAllocation(HEGData,SfDSum,SfLength)
      SfSum=0.0_pr
      SfXSum=0.0_pr
      SfDSum=0.0_pr
      SfG=0.0_pr
      SfV=0.0_pr

! THis needs to be fixed for momentum symmetry
      HEGData%ECorr = Zero
      Eaaaa = Zero
      Eabab = Zero
      Eabba = Zero
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle
        V_ijab = ERI(HEGData,UEGInfo,A,B,I,J,DummyFlag=UEGInfo%IRangeEnergy)
        V_ijba = ERI(HEGData,UEGInfo,A,B,J,I,DummyFlag=UEGInfo%IRangeEnergy)
        If(Min(V_ijab,V_ijba) < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
!       Eaaaa = Eaaaa + T2aaaa(I,J,A)*(V_ijab - V_ijba)
!       Eabab = Eabab + T2abab(I,J,A)*V_ijab
!       Eabba = Eabba - T2abba(I,J,A)*V_ijba
        dE = V_ijab*(T2aaaa(I,J,A) + T2abab(I,J,A))    &
           - V_ijba*(T2aaaa(I,J,A) + T2abba(I,J,A))
        HEGData%ECorr = HEGData%ECorr + F12*dE
        
        ! Set up for the direct term in V
        this_qvec=dK(HEGData,UEGInfo%CoreN,a,b,i,j,inCore)
        SfSum(this_qvec(1),this_qvec(2),this_qvec(3))=SfSum(this_qvec(1),this_qvec(2),this_qvec(3))+(T2aaaa(I,J,A)+ T2abab(I,J,A))
        this_q2=dK2scaled(HEGData,UEGInfo%CoreN,a,b,i,j,inCore)
        SfG(this_qvec(1),this_qvec(2),this_qvec(3))=this_q2**0.5_pr
        SfV(this_qvec(1),this_qvec(2),this_qvec(3))=V_ijab
        
        ! Accumulate direct/exchange pieces
        SfDSum(this_qvec(1),this_qvec(2),this_qvec(3))=SfDSum(this_qvec(1),this_qvec(2),this_qvec(3))+(2.0_pr*T2abab(I,J,A))
        SfXSum(this_qvec(1),this_qvec(2),this_qvec(3))=SfXSum(this_qvec(1),this_qvec(2),this_qvec(3))+(T2aaaa(I,J,A)- T2abab(I,J,A))
        
        ! Set up for the exchange term in V
        this_qvec=dK(HEGData,UEGInfo%CoreN,a,b,j,i,inCore)
        SfSum(this_qvec(1),this_qvec(2),this_qvec(3))=SfSum(this_qvec(1),this_qvec(2),this_qvec(3))-(T2aaaa(I,J,A) + T2abba(I,J,A))
        this_q2=dK2scaled(HEGData,UEGInfo%CoreN,a,b,j,i,inCore)
        SfG(this_qvec(1),this_qvec(2),this_qvec(3))=this_q2**0.5_pr
        SfV(this_qvec(1),this_qvec(2),this_qvec(3))=V_ijba
        
        ! Accumulate direct/exchange pieces
        ! Assume that v_ijab t_ijba = v_ijba t_ijab and that this switches the spin labels
        SfDSum(this_qvec(1),this_qvec(2),this_qvec(3))=SfDSum(this_qvec(1),this_qvec(2),this_qvec(3))-(2.0_pr*T2abba(I,J,A))
        SfXSum(this_qvec(1),this_qvec(2),this_qvec(3))=SfXSum(this_qvec(1),this_qvec(2),this_qvec(3))-(T2aaaa(I,J,A)- T2abba(I,J,A)) 
        
        !write(80,*) SfSum(this_qvec(1),this_qvec(2),this_qvec(3)) 
        !write(80,*) dK(a,b,i,j,inCore), dK2(a,b,i,j,inCore), dK2scaled(a,b,i,j,inCore), (T2aaaa(I,J,A) + T2abab(I,J,A)) 
        !write(80,*) dK(a,b,i,j,inCore), dK2(a,b,j,i,inCore), dK2scaled(a,b,j,i,inCore), -1.0_pr*(T2aaaa(I,J,A) + T2abba(I,J,A))
        !write(80,*) SfSum(this_qvec(1),this_qvec(2),this_qvec(3)) 
      End Do
      End Do
      End Do
!     Eaaaa = Eaaaa*F12
!     Eabab = Eabab*F12
!     Eabba = Eabba*F12
!     ECorr = ECorr + Eaaaa + Eabab + Eabba

      Do I = -SfLength,SfLength
      Do J = -SfLength,SfLength
      Do A = -SfLength,SfLength
        if (SfG(I,J,A) > 1.0E-10_pr) then
            write(80,*) I, J, A, SfG(I,J,A), SfSum(I,J,A), SfV(I,J,A), SfDSum(I,J,A), SfXSum(I,J,A)
        endif
      End Do
      End Do
      End Do

      deallocate(SfSum,SfG,SfV,SfDSum,SfXSum)

      Return
      End Subroutine CCEnergyAndSF


      Subroutine SolveCC(HEGData,UEGInfo,T2,G2)

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


      Subroutine CnvrgCC(T2abab,OldT2abab,T2abba,OldT2abba,dT,NOcc,NAO,NIter,NDIIS)

      Integer,        Intent(In)    :: NOcc, NAO, NDIIS
      Integer,        Intent(InOut) :: NIter
      Real (Kind=pr), Intent(Out)   :: dT
      Real (Kind=pr), Intent(In)    :: T2abab(NOcc,NOcc,NOcc+1:NAO)
      Real (Kind=pr), Intent(In)    :: T2abba(NOcc,NOcc,NOcc+1:NAO)
      Real (Kind=pr), Intent(InOut) :: OldT2abab(NOcc,NOcc,NOcc+1:NAO,NDIIS)
      Real (Kind=pr), Intent(InOut) :: OldT2abba(NOcc,NOcc,NOcc+1:NAO,NDIIS)

      Integer, Parameter :: MaxIter = 10000
      Real (Kind=pr) :: dTabab, dTabba

!======================================!
!  Update the DIIS amplitudes OldT2*.  !
!======================================!

! Update the iteration count
      NIter = NIter + 1

! Move the DIIS amplitudes
      OldT2abab(:,:,:,1:NDIIS-1) = OldT2abab(:,:,:,2:NDIIS)
      OldT2abba(:,:,:,1:NDIIS-1) = OldT2abba(:,:,:,2:NDIIS)

! Get the change in T
      dTabab = MaxVal(Abs(T2abab - OldT2abab(:,:,:,NDIIS)))
      dTabba = MaxVal(Abs(T2abba - OldT2abba(:,:,:,NDIIS)))
      dT = Max(dTabab, dTabba)

! Update the most recent DIIS amplitudes
      OldT2abab(:,:,:,NDIIS) = T2abab
      OldT2abba(:,:,:,NDIIS) = T2abba
      If(NIter == MaxIter) Stop 'CCD equations did not converge'

      Return
      End Subroutine CnvrgCC


      Subroutine GetG2Drivers(HEGData,UEGInfo,G2abab,G2abba)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(Out) :: G2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                     G2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Integer :: I, J, A, B
      Real (Kind=pr) :: V_ijab, V_ijba

!==================================!
!  Build the driving terms first.  !
!==================================!

      G2abab = Zero
      G2abba = Zero
!$omp parallel
!$omp do collapse(3) private(B, V_ijab, V_ijba)
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle
        V_ijab = ERI(HEGData,UEGInfo,I,J,A,B,DummyFlag=UEGInfo%IRangeDriverDirect)
        V_ijba = ERI(HEGData,UEGInfo,I,J,B,A,DummyFlag=UEGInfo%IRangeDriverExchange)
        If(Min(V_ijab,V_ijba) < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
        G2abab(I,J,A) =  V_ijab
        G2abba(I,J,A) = -V_ijba
      End Do
      End Do
      End Do
!$omp end do
!$omp end parallel
      Return
      End Subroutine GetG2Drivers


      Subroutine GetG2Mosaics(HEGData,UEGInfo,G2abab,G2abba,Joo,Jvv,T2aaaa,T2abab,T2abba)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(InOut) :: G2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                       G2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(InOut) :: Joo(UEGInfo%NOcc),                                            &
                                       Jvv(UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(In) :: T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Integer :: I, J, A, B, K, L, C, D
      Real (Kind=pr) :: V_cdkl, V_cdlk

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
!$omp parallel 
!$omp do collapse(3) reduction(+:Jvv, Joo) private(D, V_cdkl, V_cdlk)
      Do K = 1,UEGInfo%NOcc
      Do L = 1,UEGInfo%NOcc
      Do C = UEGInfo%NOcc+1,UEGInfo%NAO
        D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,K,L,C)
        If(D <= UEGInfo%NOcc) Cycle
        V_cdkl = ERI(HEGData,UEGInfo,C,D,K,L,DummyFlag=UEGInfo%IRangeMosaic)
        V_cdlk = ERI(HEGData,UEGInfo,C,D,L,K,DummyFlag=UEGInfo%IRangeMosaic)
        If(Min(V_cdkl,V_cdlk) < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
        Jvv(C) = Jvv(C) + V_cdkl*(T2aaaa(K,L,C) + T2abab(K,L,C))               &
                        - V_cdlk*(T2aaaa(K,L,C) + T2abba(K,L,C))
        Joo(K) = Joo(K) + V_cdkl*(T2aaaa(K,L,C) + T2abab(K,L,C))               &
                        - V_cdlk*(T2aaaa(K,L,C) + T2abba(K,L,C))
      End Do
      End Do
      End Do
!$omp end do
      
! The intermediates are done, so contract with T2
!$omp do collapse(3) private(B)
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle
        G2abab(I,J,A) = G2abab(I,J,A) - F12*(Jvv(A)*T2abab(I,J,A) + Jvv(B)*T2abab(I,J,A))
        G2abba(I,J,A) = G2abba(I,J,A) - F12*(Jvv(A)*T2abba(I,J,A) + Jvv(B)*T2abba(I,J,A))
        G2abab(I,J,A) = G2abab(I,J,A) - F12*(Joo(I)*T2abab(I,J,A) + Joo(J)*T2abab(I,J,A))
        G2abba(I,J,A) = G2abba(I,J,A) - F12*(Joo(I)*T2abba(I,J,A) + Joo(J)*T2abba(I,J,A))
      End Do
      End Do
      End Do
!$omp end do
!$omp end parallel

      Return
      End Subroutine GetG2Mosaics


      Subroutine GetG2Ladders(HEGData,UEGInfo,G2abab,G2abba,T2aaaa,T2abab,T2abba)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(InOut) :: G2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                       G2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(In) :: T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Integer :: I, J, A, B, K, L, C, D, MaxDummyFlag
      Real (Kind=pr) :: J2_ijkl, J3_ijkl, V_cdkl, V_cdlk, V_cdab, V_cdba

!=================================================================!
!  Gets the ladder contributions here.  Specialized for the HEG.  !
!  As vectors, we have A + B = I + J = K + L = C + D.             !
!=================================================================!

! Build the intermediates for each IJKL => hole-hole ladder, ladder T.V.T
! We'll have J2_ijkl which goes with T2abab(K,L,A,B) in G2abab(I,J,A,B) and with T2abba(K,L,A,B) in G2abba(I,J,A,B)
! We'll have J3_ijkl which goes with T2abba(K,L,A,B) in G2abab(I,J,A,B) and with T2abab(K,L,A,B) in G2abba(I,J,A,B)
!$omp parallel
!$omp do collapse(3) &
!$omp&   private(C, A, L, D, B, J2_ijkl, J3_ijkl, V_cdkl, V_cdlk, MaxDummyFlag)
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
        Do K = 1,UEGInfo%NOcc
          L = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,K)
          If(L > UEGInfo%NOcc .or. L <= 0) Cycle
          MaxDummyFlag = Max(UEGInfo%IRangeLadder, &
                             UEGInfo%IRangeDirectLadders, &
                             UEGInfo%IRangeLinLadders)
          J2_ijkl =  ERI(HEGData,UEGInfo,I,J,K,L,DummyFlag=MaxDummyFlag)              ! For T2abab in G2abab and T2abba in G2abba
          MaxDummyFlag = Max(UEGInfo%IRangeLadder, &
                             UEGInfo%IRangeExchangeLadders, &
                             UEGInfo%IRangeLinLadders)
          J3_ijkl = -ERI(HEGData,UEGInfo,I,J,L,K,DummyFlag=MaxDummyFlag)              ! For T2abba in G2abab and T2abab in G2abba
          If(Min(J2_ijkl,-J3_ijkl) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
        ! CK: Time making this a reduction over J2 & J3. May need to ration threads.
          Do C = UEGInfo%NOcc+1,UEGInfo%NAO
            D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,C)
            If(D <= UEGInfo%NOcc) Cycle
            MaxDummyFlag = Max(UEGInfo%IRangeLadder, &
                               UEGInfo%IRangeDirectLadders, &
                               UEGInfo%IRangeQuadLadders)
            V_cdkl = ERI(HEGData,UEGInfo,C,D,K,L,DummyFlag=MaxDummyFlag)
            MaxDummyFlag = Max(UEGInfo%IRangeLadder, &
                               UEGInfo%IRangeExchangeLadders, &
                               UEGInfo%IRangeQuadLadders)
            V_cdlk = ERI(HEGData,UEGInfo,C,D,L,K,DummyFlag=MaxDummyFlag)
            If(Min(V_cdkl,V_cdlk) < -80.0_pr) Then
              Print *, "Disallowed excitation not trapped!"
             Cycle
            End If
            J2_ijkl = J2_ijkl  + F12*(V_cdkl*T2abab(I,J,C) - V_cdlk*T2abba(I,J,C))
            J3_ijkl = J3_ijkl  + F12*(V_cdkl*T2abba(I,J,C) - V_cdlk*T2abab(I,J,C))
          End Do
! The intermediates are done, so contract with T2
          Do A = UEGInfo%NOcc+1,UEGInfo%NAO
            B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
            If(B <= UEGInfo%NOcc) Cycle
            G2abab(I,J,A) = G2abab(I,J,A) + F12*(J2_ijkl*T2abab(K,L,A) + J3_ijkl*T2abba(K,L,A))
            G2abba(I,J,A) = G2abba(I,J,A) + F12*(J2_ijkl*T2abba(K,L,A) + J3_ijkl*T2abab(K,L,A))
          End Do
        End Do
      End Do
      End Do
!$omp end do

!$omp do collapse(3) private(A, C, B, D, V_cdab, V_cdba, MaxDummyFlag)
! Add the particle-particle ladder
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
        Do A = UEGInfo%NOcc+1,UEGInfo%NAO
          B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
          If(B <= UEGInfo%NOcc) Cycle
          Do C = UEGInfo%NOcc+1,UEGInfo%NAO
            D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,C)
            If(D <= UEGInfo%NOcc) Cycle
            MaxDummyFlag = Max(UEGInfo%IRangeLadder, &
                               UEGInfo%IRangeDirectLadders, &
                               UEGInfo%IRangeLinLadders)
            V_cdab = ERI(HEGData,UEGInfo,C,D,A,B,DummyFlag=MaxDummyFlag)
            MaxDummyFlag = Max(UEGInfo%IRangeLadder, &
                               UEGInfo%IRangeExchangeLadders, &
                               UEGInfo%IRangeLinLadders)
            V_cdba = ERI(HEGData,UEGInfo,C,D,B,A,DummyFlag=MaxDummyFlag)
            If(Min(V_cdab,V_cdba) < -80.0_pr) Then
             Print *, "Disallowed excitation not trapped!"
             Cycle
            End If
            G2abab(I,J,A) = G2abab(I,J,A) + F12*(V_cdab*T2abab(I,J,C) - V_cdba*T2abba(I,J,C))
            G2abba(I,J,A) = G2abba(I,J,A) + F12*(V_cdab*T2abba(I,J,C) - V_cdba*T2abab(I,J,C))
          End Do
        End Do
      End Do
      End Do
!$omp end do
!$omp end parallel

      Return
      End Subroutine GetG2Ladders


      Subroutine GetG2Rings(HEGData,UEGInfo,G2abab,G2abba,T2aaaa,T2abab,T2abba)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(InOut) :: G2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                       G2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(In) :: T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Integer :: I, J, A, B, K, L, C, D, MaxDummyFlag
      Real (Kind=pr) :: J1_idal, J2_idal, J3_idal
      Real (Kind=pr) :: V_cjkb, V_cjbk, V_cika, V_ciak, V_cdkl, V_cdlk

!===========================================================!
!  Build the ring terms.  These are                         !
!   dG_ij^ab = t_{ik}^{ac} vbar^{kb}_{cj}                   !
!            + t_{kj}^{cb} vbar^{ak}_{ic}                   !
!            + t_{ik}^{ac} t_{lj}^{db} vbar^{kl}_{cd}       !
!===========================================================!

!$omp parallel
!$omp do collapse(3) &
!$omp&   private(K, B, C, V_cjkb, V_cjbk, V_cika, V_ciak, MaxDummyFlag)
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
          MaxDummyFlag = Max(UEGInfo%IRangeRing, &
                             UEGInfo%IRangeExchangeRings, &
                             UEGInfo%IRangeLinRings)
          V_cjbk = ERI(HEGData,UEGInfo,C,J,B,K,DummyFlag=MaxDummyFlag)
          If(Min(V_cjkb, V_cjbk) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2abab(I,J,A) = G2abab(I,J,A) + (V_cjkb - V_cjbk)*T2abab(I,K,A) +  V_cjkb*T2aaaa(I,K,A)
          G2abba(I,J,A) = G2abba(I,J,A) - V_cjbk*T2abba(I,K,A)
        End Do

! Do the t_jk^bc V_ic^ak contraction
        Do K = 1,UEGInfo%NOcc
          C = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,J,K,B)
          If(C <= UEGInfo%NOcc) Cycle
          MaxDummyFlag = Max(UEGInfo%IRangeRing, &
                             UEGInfo%IRangeDirectRings, &
                             UEGInfo%IRangeLinRings)
          V_cika = ERI(HEGData,UEGInfo,C,I,K,A,DummyFlag=MaxDummyFlag)
          MaxDummyFlag = Max(UEGInfo%IRangeRing, &
                             UEGInfo%IRangeExchangeRings, &
                             UEGInfo%IRangeLinRings)
          V_ciak = ERI(HEGData,UEGInfo,C,I,A,K,DummyFlag=MaxDummyFlag)
          If(Min(V_cika, V_ciak) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2abab(I,J,A) = G2abab(I,J,A) + (V_cika - V_ciak)*T2abab(J,K,B) + V_cika*T2aaaa(J,K,B)
          G2abba(I,J,A) = G2abba(I,J,A) - V_ciak*T2abba(J,K,B)
        End Do
      End Do
      End Do
      End Do
!$omp end do

! These are the quadratic terms.  They require intermediates.
! We'll have J1_idal, which contracts with T2aaaa(L,J,D,B) in G2aaaa(I,J,A,B) and with T2abab(L,J,D,B) in G2abab(I,J,A,B)
! We'll have J2_idal, which contracts with T2abab(L,J,D,B) in G2aaaa(I,J,A,B) and with T2aaaa(L,J,D,B) in G2abab(I,J,A,B
! We'll have J3_idal, which contracts with T2abba(L,J,D,B) in G2abba(I,J,A,B)
! This one is a true pain in the ass.
!$omp do collapse(3) &
!$omp&   private(K, J, C, D, B, J1_idal, J2_idal, J3_idal, V_cdkl, V_cdlk, MaxDummyFlag)
      Do I = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        Do L = 1,UEGInfo%NOcc
          J1_idal = Zero   ! Intermediate for T2aaaa(L,J,D,B) in G2aaaa(I,J,A,B) and T2abab(L,J,D,B) in G2abab(I,J,A,B)
          J2_idal = Zero   ! Intermediate for T2abab(L,J,D,B) in G2aaaa(I,J,A,B) and T2aaaa(L,J,D,B) in G2abab(I,J,A,B)
          J3_idal = Zero   ! Intermediate for T2abba(L,J,D,B) in G2abba(I,J,A,B)
          Do K = 1,UEGInfo%NOcc
            C = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,K,A)
            If(C <= UEGInfo%NOcc) Cycle
            D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,K,L,C)
            If(D <= UEGInfo%NOcc) Cycle
            MaxDummyFlag = Max(UEGInfo%IRangeRing, &
                               UEGInfo%IRangeDirectRings, &
                               UEGInfo%IRangeQuadRings)
            V_cdkl = ERI(HEGData,UEGInfo,C,D,K,L,DummyFlag=MaxDummyFlag)
            MaxDummyFlag = Max(UEGInfo%IRangeRing, &
                               UEGInfo%IRangeExchangeRings, &
                               UEGInfo%IRangeQuadRings)
            V_cdlk = ERI(HEGData,UEGInfo,C,D,L,K,DummyFlag=MaxDummyFlag)
            If(Min(V_cdkl,V_cdlk) < -80.0_pr) Then
              Print *, "Disallowed excitation not trapped!"
              Cycle
            End If
            J1_idal = J1_idal + T2aaaa(I,K,A)*(V_cdkl - V_cdlk) + T2abab(I,K,A)*V_cdkl
            J2_idal = J2_idal + T2abab(I,K,A)*(V_cdkl - V_cdlk) + T2aaaa(I,K,A)*V_cdkl
            J3_idal = J3_idal - T2abba(I,K,A)*V_cdlk
          End Do

! The intermediates are done, so contract with T2
! In principle, we've got the right D here.
          Do J = 1,UEGInfo%NOcc
            B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
            If(B <= UEGInfo%NOcc) Cycle
            D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,J,L,B)
            If(D <= UEGInfo%NOcc) Cycle
            G2abab(I,J,A) = G2abab(I,J,A) + J2_idal*T2aaaa(L,J,D) + J1_idal*T2abab(L,J,D)
            G2abba(I,J,A) = G2abba(I,J,A) + J3_idal*T2abba(L,J,D)
          End Do
        End Do
      End Do
      End Do
!$omp end do
!$omp end parallel
      Return
      End Subroutine GetG2Rings


      Subroutine GetG2XRings(HEGData,UEGInfo,G2abab,G2abba,T2aaaa,T2abab,T2abba)

      Use HEG, Only: FindIndex
      Use Types, only: HEGDataType, UEGInfoType

      Type (HEGDataType), Intent(In) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo

      Real (Kind=pr), Intent(InOut) :: G2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                       G2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)
      Real (Kind=pr), Intent(In) :: T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                    T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

! Local variables
      Integer :: I, J, A, B, K, L, C, D, MaxDummyFlag
      Real (Kind=pr) :: J1_idlb, J2_idlb, J3_idlb
      Real (Kind=pr) :: V_jcka, V_ickb, V_jcak, V_icbk, V_cdkl, V_cdlk

!===========================================================!
!  Build the crossed-ring terms.  These are                 !
!   dG_ij^ab = - t_{ik}^{cb} vbar^{ka}_{jc}                 !
!              - t_{jk}^{ca} vbar^{kb}_{ic}                 !
!              - t_{ik}^{cb} t_{jl}^{da} vbar^{kl}_{cd}     !
!  Most of the signs drop out on using vbar!                !
!===========================================================!

!$omp parallel
!$omp do collapse(3) private(K, B, C, V_jcka, V_jcak, V_ickb, V_icbk, MaxDummyFlag)
! These are the linear terms
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle

! Do the t_ik^cb V_jc^ka contraction
        Do K = 1,UEGInfo%NOcc
          C = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,K,I,B)
          If(C <= UEGInfo%NOcc) Cycle
          MaxDummyFlag = Max(UEGInfo%IRangeXRing, &
                             UEGInfo%IRangeExchangeRings, &
                             UEGInfo%IRangeLinRings)
          V_jcka = ERI(HEGData,UEGInfo,J,C,K,A,DummyFlag=MaxDummyFlag)
          MaxDummyFlag = Max(UEGInfo%IRangeXRing, &
                             UEGInfo%IRangeDirectRings, &
                             UEGInfo%IRangeLinRings)
          V_jcak = ERI(HEGData,UEGInfo,J,C,A,K,DummyFlag=MaxDummyFlag)
          If(Min(V_jcka,V_jcak) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2abab(I,J,A) = G2abab(I,J,A) - V_jcka*T2abab(K,I,B)
          G2abba(I,J,A) = G2abba(I,J,A) + (V_jcak - V_jcka)*T2abba(K,I,B)  +  V_jcak*T2aaaa(K,I,B)
        End Do

! Do the t_jk^ca v_ic^kb contraction
        Do K = 1,UEGInfo%NOcc
          C = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,K,J,A) ! should return the same as GetIndex(I,K,B), but check!
          If(C <= UEGInfo%NOcc) Cycle
          MaxDummyFlag = Max(UEGInfo%IRangeXRing, &
                             UEGInfo%IRangeExchangeRings, &
                             UEGInfo%IRangeLinRings)
          V_ickb = ERI(HEGData,UEGInfo,I,C,K,B,DummyFlag=MaxDummyFlag)
          MaxDummyFlag = Max(UEGInfo%IRangeXRing, &
                             UEGInfo%IRangeDirectRings, &
                             UEGInfo%IRangeLinRings)
          V_icbk = ERI(HEGData,UEGInfo,I,C,B,K,DummyFlag=MaxDummyFlag)
          If(Min(V_ickb,V_icbk) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2abab(I,J,A) = G2abab(I,J,A) - V_ickb*T2abab(K,J,A)
          G2abba(I,J,A) = G2abba(I,J,A) + (V_icbk - V_ickb)*T2abba(K,J,A) + V_icbk*T2aaaa(K,J,A)
        End Do
      End Do
      End Do
      End Do
!$omp end do

! These are the quadratic terms.  They require intermediates.
! We'll have J1_idlb, which contracts with T2aaaa(L,J,A,D) in G2aaaa(I,J,A,B) and with T2abba(L,J,A,D) in G2abba(I,J,A,B)
! We'll have J2_idlb, which contracts with T2abba(L,J,A,D) in G2aaaa(I,J,A,B) and with T2aaaa(L,J,A,D) in G2abba(I,J,A,B
! We'll have J3_idlb, which contracts with T2abab(L,J,A,D) in G2abab(I,J,A,B)
!$omp do collapse(3) &
!$omp&   private(K, J, C, D, A, J1_idlb, J2_idlb, J3_idlb, V_cdkl, V_cdlk, MaxDummyFlag)
      Do I = 1,UEGInfo%NOcc
      Do B = UEGInfo%NOcc+1,UEGInfo%NAO
        Do L = 1,UEGInfo%NOcc
          J1_idlb = Zero   ! Intermediate for T2aaaa(L,J,A,D) in G2aaaa(I,J,A,B) and T2abba(L,J,A,D) in G2abba(I,J,A,B)
          J2_idlb = Zero   ! Intermediate for T2abba(L,J,A,D) in G2aaaa(I,J,A,B) and T2aaaa(L,J,A,D) in G2abba(I,J,A,B
          J3_idlb = Zero   ! Intermediate for T2abab(L,J,D,B) in G2abab(I,J,A,B)
          Do K = 1,UEGInfo%NOcc
            C = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,K,B)
            If(C <= UEGInfo%NOcc) Cycle
            D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,K,L,C)
            If(D <= UEGInfo%NOcc) Cycle
            MaxDummyFlag = Max(UEGInfo%IRangeXRing, &
                               UEGInfo%IRangeDirectRings, &
                               UEGInfo%IRangeQuadRings)
            V_cdkl = ERI(HEGData,UEGInfo,C,D,K,L,DummyFlag=MaxDummyFlag)
            MaxDummyFlag = Max(UEGInfo%IRangeXRing, &
                               UEGInfo%IRangeExchangeRings, &
                               UEGInfo%IRangeQuadRings)
            V_cdlk = ERI(HEGData,UEGInfo,C,D,L,K,DummyFlag=MaxDummyFlag)
            If(Min(V_cdkl,V_cdlk) < -80.0_pr) Then
              Print *, "Disallowed excitation not trapped!"
              Cycle
            End If
            J1_idlb = J1_idlb + T2aaaa(K,I,B)*(V_cdkl - V_cdlk) + T2abba(K,I,B)*V_cdkl
            J2_idlb = J2_idlb + T2abba(K,I,B)*(V_cdkl - V_cdlk) + T2aaaa(K,I,B)*V_cdkl
            J3_idlb = J3_idlb + T2abab(K,I,B)*V_cdlk
          End Do
! The intermediates are done, so contract with T2
          Do J = 1,UEGInfo%NOcc
            A = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,J,I,B)
            If(A <= UEGInfo%NOcc) Cycle
            D = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,L,J,A)
            If(D <= UEGInfo%NOcc) Cycle
            G2abba(I,J,A) = G2abba(I,J,A) - J1_idlb*T2abba(L,J,A) - J2_idlb*T2aaaa(L,J,A)
            G2abab(I,J,A) = G2abab(I,J,A) + J3_idlb*T2abab(L,J,A)
          End Do
        End Do
      End Do
      End Do
!$omp end do
!$omp end parallel
      Return
      End Subroutine GetG2XRings

   End Module CCD
