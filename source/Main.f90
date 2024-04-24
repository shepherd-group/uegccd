
      Program TermedCC
      Use Precision
      Use IO, only: ReadInput
      Use HEG, only: init_HEG, change_rs, GenerateTwist, ChangeTwist
      Use MP2, only: DrvMBPT
      Use CCD, only: DrvCCD
      Use dRPA, only: DrvdRPA
      Use Types, only: UEGInfoType, HEGDataType, TwistAngleCalculation, TwistAngleData
      Use Pointers
      Implicit None
! Dimensioning variables
      Type (UEGInfoType) :: UEGInfo
      Type (HEGDataType) :: HEGData
      Type (TwistAngleCalculation) :: TwistAverages
      Type (TwistAngleData), Allocatable :: TwistData(:)
      Real (Kind=pr) :: rS
      Integer :: iRSPoint
! Correlated stuff
      Real (Kind=pr), Allocatable :: T2aaaa(:,:,:), T2abab(:,:,:), T2abba(:,:,:)
      Real (Kind=pr), Allocatable :: X2aaaa(:,:,:), X2abab(:,:,:), X2abba(:,:,:)
! Error checking variables
      Integer, Parameter :: NAlloc = 6
      Integer :: IAlloc(NAlloc)
      !Logical, Parameter :: T = .true., F=.false.
! Twist angle stuff

      Integer :: iTwist
      Integer :: iConnPoints, iConnMax
      integer,parameter :: seed = 86456
      real(pr) :: ConnectivityDiff(400)
      real(pr) :: ConnectivityDiff2
      Real (kind=pr) :: TwistAveragedEnergies(2,4)

      TwistAveragedEnergies = 0.0_pr

!==========================================!
!  This code implements RHF-based CCD.     !
!  We give it the option of keeping the    !
!  ring, ladder, crossed-ring, and mosaic  !
!  terms each on a case-by-case basis.     !
!------------------------------------------!
!  The first thing we must do is read the  !
!  basic information for the calculation.  !
!==========================================!

      Call ReadInput(UEGInfo)

      Call Init_HEG(UEGInfo,HEGData)

      Call SetPointers(UEGInfo)
!==========================================!
!  Now we can allocate the memory and go!  !
!==========================================!

      IAlloc = 0
      Allocate(T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO),  Stat=IAlloc(1))
      Allocate(T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO),  Stat=IAlloc(2))
      Allocate(T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO),  Stat=IAlloc(3))
      Allocate(X2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO),  Stat=IAlloc(4))
      Allocate(X2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO),  Stat=IAlloc(5))
      Allocate(X2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO),  Stat=IAlloc(6))
      If(Any(IAlloc /= 0)) Stop "Could not allocate in main"
      Open(7,File='Output',Status="Replace")
      Close(7)

! TODO: This would be the place where we could loop over the twist
! angles. Inside the rs loop, and subject to an input parameter.
! - New input parameter NumTwists
! - ChangeTwist subroutine -- DONE
! - Store the CCD energy (etc) between DrvCCD calls. -- DONE
! - Average them over multiple iterations 
! - Consider other data to collect, perhaps in a derived type? -- DONE
! - "Characteristic momentum vector"


      ! Set the system rs
      UEGInfo%rS = UEGInfo%rSMin
      call change_rs(HEGData,UEGInfo)

      call srand(seed)

      ! If we are only doing TAHF, run that and stop.
      If (UEGInfo%DoCalcOnlyTAHF) then
          Do iTwist = 1, UEGInfo%TACalcN
              Call GenerateTwist(HEGData,UEGInfo%DoTwistx,UEGInfo%DoTwisty,UEGInfo%DoTwistz) ! TODO: this needs testing
              Call ChangeTwist(HEGData,UEGInfo) ! TODO : this needs testing
              Call WriteTwistData(iTwist,HEGData,UEGInfo,TwistAveragedEnergies)
          End Do
          ! We are done, stop the program.
          Call SuccessfulDeallocateAndStop(T2aaaa,T2abab,T2abba,X2aaaa,X2abab,X2abba)
      End If

      if (UEGInfo%DoSingleCalc) then 
          Call DrvMBPT(HEGData,UEGInfo,X2aaaa,X2abab,X2abba)

          T2aaaa = X2aaaa
          T2abab = X2abab
          T2abba = X2abba

          if (.not.UEGInfo%DodRPASOSEX) then
              Call DrvCCD(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)
          else
              Call DrvdRPA(HEGData,UEGInfo,T2aaaa)
          endif
      endif
      if (UEGInfo%DoSkipTA) then
          Stop "Finishing calculation after single point calculation"
      endif

      ! Ready for accumulation of TA stats
      TwistAverages%CcdSum=0.0_pr
      TwistAverages%Mp2Sum=0.0_pr
      TwistAverages%HfSum=0.0_pr
      TwistAverages%NTwistSum=0.0_pr
      TwistAverages%CcdSum2=0.0_pr
      TwistAverages%Mp2Sum2=0.0_pr
      TwistAverages%HfSum2=0.0_pr
      TwistAverages%NTwistSum2=0.0_pr
      TwistAverages%ConnectivitySum=0

      ! set in code, but could be set in input file    
      TwistAverages%NumTwists=UEGInfo%TACalcN
      allocate(TwistData(TwistAverages%NumTwists))

      ! Loop through each twist
      Do iTwist = 1,TwistAverages%NumTwists
          write(60,*) iTwist,"of", TwistAverages%NumTwists
          Call GenerateTwist(HEGData,UEGInfo%DoTwistx,UEGInfo%DoTwisty,UEGInfo%DoTwistz) ! TODO: this needs testing
          Call ChangeTwist(HEGData,UEGInfo) ! TODO : this needs testing

          Call WriteTwistData(iTwist,HEGData,UEGInfo,TwistAveragedEnergies)

          TwistData(iTwist)%ntwist = HEGData%ntwist
          If (iTwist == 1) Then
              write(60,*)"EIGEN",size(HEGData%Eigen)
              TwistAverages%EigenLength=size(HEGData%Eigen)
              Allocate(TwistAverages%EigenSum(TwistAverages%EigenLength))
              TwistAverages%EigenSum=0.0_pr
              Allocate(TwistAverages%EigenAverage(TwistAverages%EigenLength))
              TwistAverages%EigenAverage=0.0_pr
          End If

          Allocate(TwistData(iTwist)%Eigen(TwistAverages%EigenLength))
          Call DrvMBPT(HEGData,UEGInfo,X2aaaa,X2abab,X2abba,Conn=TwistData(iTwist)%Connectivity)
          TwistData(iTwist)%Eigen=HEGData%Eigen
          TwistData(iTwist)%mp2=HEGData%ECorr
          TwistData(iTwist)%hf=HEGData%EHF

          T2aaaa = X2aaaa
          T2abab = X2abab
          T2abba = X2abba

          if (UEGInfo%DoCalcTACCD) then
              if (.not.UEGInfo%DodRPASOSEX) then
                  Call DrvCCD(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)
              else
                  Call DrvdRPA(HEGData,UEGInfo,T2aaaa)
              endif
          endif
          TwistData(iTwist)%ccd=HEGData%ECorr
          write(59,*) TwistData(iTwist)%ntwist, TwistData(iTwist)%hf,TwistData(iTwist)%mp2, TwistData(iTwist)%ccd
          If (iTwist == 1) Then
              Do iConnPoints = size(TwistData(iTwist)%Connectivity),1,-1
                  If (TwistData(iTwist)%Connectivity(iConnPoints) > 0) Then
                      iConnMax=iConnPoints
                      Exit
                  End If
              End Do
          End If
          write(58,*) TwistData(iTwist)%Connectivity(:iConnMax)
          !write(59,*) TwistData(iTwist)%Connectivity
          ! TODO -- accumulate statistics 
          ! for <x>
          TwistAverages%CcdSum=TwistAverages%CcdSum+TwistData(iTwist)%ccd
          TwistAverages%Mp2Sum=TwistAverages%Mp2Sum+TwistData(iTwist)%mp2
          TwistAverages%HfSum=TwistAverages%HfSum+TwistData(iTwist)%hf
          TwistAverages%NTwistSum=TwistAverages%NTwistSum+TwistData(iTwist)%ntwist
          TwistAverages%ConnectivitySum=TwistAverages%ConnectivitySum+TwistData(iTwist)%Connectivity
          TwistAverages%EigenSum=TwistAverages%EigenSum+TwistData(iTwist)%Eigen
          ! for <x^2>
          TwistAverages%CcdSum2=TwistAverages%CcdSum2+(TwistData(iTwist)%ccd)**2.0_pr
          TwistAverages%Mp2Sum2=TwistAverages%Mp2Sum2+(TwistData(iTwist)%mp2)**2.0_pr
          TwistAverages%HfSum2=TwistAverages%HfSum2+(TwistData(iTwist)%hf)**2.0_pr
          TwistAverages%NTwistSum2=TwistAverages%NTwistSum2+(TwistData(iTwist)%ntwist)**2.0_pr
      End Do
      ! TODO -- calc averages
      ! for <x>
      TwistAverages%NTwist=TwistAverages%NTwistSum/TwistAverages%NumTwists
      TwistAverages%Hf=TwistAverages%HfSum/TwistAverages%NumTwists
      TwistAverages%Mp2=TwistAverages%Mp2Sum/TwistAverages%NumTwists
      TwistAverages%Ccd=TwistAverages%CcdSum/TwistAverages%NumTwists
      TwistAverages%ConnectivityAverage=real(TwistAverages%ConnectivitySum,pr)/TwistAverages%NumTwists
      TwistAverages%EigenAverage=TwistAverages%EigenSum/TwistAverages%NumTwists
      ! for <x^2>
      TwistAverages%NTwist2=TwistAverages%NTwistSum2/TwistAverages%NumTwists
      TwistAverages%Hf2=TwistAverages%HfSum2/TwistAverages%NumTwists
      TwistAverages%Mp22=TwistAverages%Mp2Sum2/TwistAverages%NumTwists
      TwistAverages%Ccd2=TwistAverages%CcdSum2/TwistAverages%NumTwists
      ! for sqrt((1/(N-1))*(<x^2>-<x>^2))
      TwistAverages%NTwistErr=((TwistAverages%NTwist2-TwistAverages%NTwist**2.0_pr)/(TwistAverages%NumTwists-1))**0.5_pr
      TwistAverages%HfErr=((TwistAverages%Hf2-TwistAverages%Hf**2.0_pr)/(TwistAverages%NumTwists-1))**0.5_pr
      TwistAverages%Mp2Err=((TwistAverages%Mp22-TwistAverages%Mp2**2.0_pr)/(TwistAverages%NumTwists-1))**0.5_pr
      TwistAverages%CcdErr=((TwistAverages%Ccd2-TwistAverages%Ccd**2.0_pr)/(TwistAverages%NumTwists-1))**0.5_pr
      write(59,*) TwistAverages%NTwist, TwistAverages%Hf, TwistAverages%Mp2, TwistAverages%Ccd, "averages"
      write(59,*) TwistAverages%NTwist2, TwistAverages%Hf2, TwistAverages%Mp22, TwistAverages%Ccd2, "squared averages"
      write(59,*) TwistAverages%NTwistErr, TwistAverages%HfErr, TwistAverages%Mp2Err, TwistAverages%CcdErr, "errors"
      write(59,*) TwistAverages%EigenAverage, "eigenvalue averages"
      Do iConnPoints = size(TwistAverages%ConnectivitySum),1,-1
          !write(59,*) iConnPoints,TwistAverages%ConnectivitySum(iConnPoints)
          If (TwistAverages%ConnectivitySum(iConnPoints) > 0) Then
              iConnMax=iConnPoints
              Exit
          End If
          !write(59,*) TwistAverages%ConnectivityAverage, "average connectivity"
      End Do
      write(59,*) TwistAverages%ConnectivityAverage(:iConnMax), "conn"

! -- Eigenvalue averages are stored in TwistAverages%EigenAverage
! -- Connectivity for each twist angle is: TwistData(iTwist)%Connectivity
! -- Eigenvalue averages are stored in TwistAverages%ConnectivityAverage
! ====================================
! TODO:
! -- Loop over twist angle, and find the connectivity that best matches
! the average connectivity
! -- Run the MP2 and CCD with the averaged eigenvalues and the matrix 
! elements for the point that best matches the connectivity
! -- Copy this same routine from a previous version of the code that was
! less heavily modified
      
      TwistAverages%MinConnectivityDiff2=1.0_pr*100000
      Do iTwist = 1,TwistAverages%NumTwists
        ConnectivityDiff=TwistData(iTwist)%Connectivity-TwistAverages%ConnectivityAverage
        Do iConnPoints = 1,size(ConnectivityDiff)
            ConnectivityDiff(iConnPoints)=ConnectivityDiff(iConnPoints)/iConnPoints**2.0_pr
        End Do
        ConnectivityDiff2=dot_product(ConnectivityDiff,ConnectivityDiff)
        If (iTwist == 1) Then
            TwistAverages%MinConnectivityDiff2=ConnectivityDiff2
            TwistAverages%SpecialTwistAngle=TwistData(iTwist)%ntwist
        End If
        If (ConnectivityDiff2 < TwistAverages%MinConnectivityDiff2) Then
            TwistAverages%MinConnectivityDiff2=ConnectivityDiff2
            TwistAverages%SpecialTwistAngle=TwistData(iTwist)%ntwist
        End If
        write(59,*) iTwist,ConnectivityDiff2
      End Do
      write(59,*) "special twist angle found",TwistAverages%SpecialTwistAngle
      write(59,*) "lowest connectivity diff squared", TwistAverages%MinConnectivityDiff2

      HEGData%this_ntwist = TwistAverages%SpecialTwistAngle
      Call ChangeTwist(HEGData,UEGInfo)

      HEGData%Eigen = TwistAverages%EigenAverage 
      Call DrvMBPT(HEGData,UEGInfo,X2aaaa,X2abab,X2abba)

      T2aaaa = X2aaaa
      T2abab = X2abab
      T2abba = X2abba

      if (.not.UEGInfo%DodRPASOSEX) then
          Call DrvCCD(HEGData,UEGInfo,T2aaaa,T2abab,T2abba)
      else
          Call DrvdRPA(HEGData,UEGInfo,T2aaaa)
      endif


!==============================================!
!  Lastly, deallocate memory and exit safely.  !
!==============================================!

      Call SuccessfulDeallocateAndStop(T2aaaa,T2abab,T2abba,X2aaaa,X2abab,X2abba)

      Contains

        Subroutine SetPointers(UEGInfo)

          Use HEG, only: CasedERI, NoCaseERI, CoreERI, GapERI, BaseERI, SphericallyTruncatedERI
          Use Types, Only: UEGInfoType
          Use Pointers

          Type (UEGInfoType), Intent(In) :: UEGInfo

          Integer :: SumRange

          Logical :: SetSphericallyTruncPointer

          SetSphericallyTruncPointer = .false.

          SumRange = 0
          SumRange = SumRange + Abs(UEGInfo%IRangeRing)
          SumRange = SumRange + Abs(UEGInfo%IRangeXRing)
          SumRange = SumRange + Abs(UEGInfo%IRangeLadder)
          SumRange = SumRange + Abs(UEGInfo%IRangeMosaic)
          SumRange = SumRange + Abs(UEGInfo%IRangeDriverDirect)
          SumRange = SumRange + Abs(UEGInfo%IRangeDriverExchange)
          SumRange = SumRange + Abs(UEGInfo%IRangeEnergy)
          SumRange = SumRange + Abs(UEGInfo%IRangeLinRings)
          SumRange = SumRange + Abs(UEGInfo%IRangeQuadRings)
          SumRange = SumRange + Abs(UEGInfo%IRangeDirectRings)
          SumRange = SumRange + Abs(UEGInfo%IRangeExchangeRings)
          SumRange = SumRange + Abs(UEGInfo%IRangeLinLadders)
          SumRange = SumRange + Abs(UEGInfo%IRangeQuadLadders)
          SumRange = SumRange + Abs(UEGInfo%IRangeDirectLadders)
          SumRange = SumRange + Abs(UEGInfo%IRangeExchangeLadders)

          If ((SumRange > 0 .And. UEGInfo%GapII > 0.0_pr .And. UEGInfo%CoreN > 0) .Or. UEGInfo%SafeERI) Then
            ! Include all the bells and whistles (all checks)
            ERI => CasedERI
          Else If (UEGInfo%GapII > 0.0_pr .And. UEGInfo%CoreN > 0) Then
            ! We can ignore the DummyFlag checks
            ERI => NoCaseERI
          Else If (UEGInfo%CoreN > 0) Then
            ! We only consider CoreN + Consistency Checks
            ERI => CoreERI
          Else If (UEGInfo%GapII > 0.0_pr) Then
            ! We only consider Gap + Consistency Checks
            ERI => GapERI
          Else If (UEGInfo%DoSphericalTruncV) Then
            ERI => SphericallyTruncatedERI
            SetSphericallyTruncPointer = .true.
          Else
            ! No checks at all, be warned!
            ERI => BaseEri
          End If

          If (UEGInfo%DoSphericalTruncV .and. .not. SetSphericallyTruncPointer) then
            Error Stop 'Failed to set the Spherical Truncation ERI pointer!'
          End If

        End Subroutine SetPointers


        Subroutine WriteTwistData(iTwist,HEGData,UEGInfo,TAE)

          ! Write the offset specific kinetic, exchange,
          ! and Hartree--Fock energy to fort.97. At the end
          ! write the twist averaged quantities as well.

          Use Types, only: HEGDataType, UEGInfoType

          Type (HEGDataType), Intent(In) :: HEGData
          Type (UEGInfoType), Intent(In) :: UEGInfo

          Integer, Intent(In) :: iTwist
          Real (Kind=pr), Intent(InOut) :: TAE(2,4)

          Real (Kind=pr) :: kinetic, exchange, total, Ncalcs

          kinetic = (HEGData%EHF - HEGData%XHF)/HEGData%Nel - HEGData%Madelung/2.0_pr
          exchange = HEGData%XHF/HEGData%Nel + HEGData%Madelung/2.0_pr
          total = HEGData%EHF - HEGData%Madelung*HEGData%Nel/2.0_pr
          total = total/HEGData%Nel + HEGData%Madelung/2.0_pr

          ! Write the header
          If (iTwist == 1) Then
              TAE = 0.0_pr
              Write(97, '(1X, "# Offset energy table")')
              Write(97, '(1X, "# -------------------")')
              Write(97, '(1X, A10, 4(A26))') '#    Twist', &
                  'Kinetic', 'Exchange', 'Total', 'Madelung'
          End If

          ! Update the running total for averaging at the end
          TAE(1,1) = TAE(1,1) + kinetic
          TAE(1,2) = TAE(1,2) + exchange
          TAE(1,3) = TAE(1,3) + total
          TAE(1,4) = TAE(1,4) + HEGData%Madelung
          TAE(2,1) = TAE(2,1) + kinetic**2.0_pr
          TAE(2,2) = TAE(2,2) + exchange**2.0_pr
          TAE(2,3) = TAE(2,3) + total**2.0_pr
          TAE(2,4) = TAE(2,4) + HEGData%Madelung**2.0_pr

          ! Write out the current off-set energies
          Write(97, '(1X, I10, 4(F26.16))') iTwist, &
              kinetic, exchange, total, HEGData%Madelung

          ! If this is the last off-set to run, report the averaged energies.
          If (iTwist == UEGInfo%TACalcN) Then
              Write(97, '()')
              Write(97, '(1X, "# Twist averaged energies and errors")')
              Write(97, '(1X, "# ----------------------------------")')
              Ncalcs = Real(UEGInfo%TACalcN, pr)
              TAE(1,1) = TAE(1,1)/Ncalcs
              TAE(1,2) = TAE(1,2)/Ncalcs
              TAE(1,3) = TAE(1,3)/Ncalcs
              TAE(1,4) = TAE(1,4)/Ncalcs
              TAE(2,1) = TAE(2,1) - Ncalcs*(TAE(1,1)**2.0_pr)
              TAE(2,2) = TAE(2,2) - Ncalcs*(TAE(1,2)**2.0_pr)
              TAE(2,3) = TAE(2,3) - Ncalcs*(TAE(1,3)**2.0_pr)
              TAE(2,4) = TAE(2,4) - Ncalcs*(TAE(1,4)**2.0_pr)
              TAE(2,1) = (TAE(2,1)/(Ncalcs - 1.0_pr))**0.5_pr
              TAE(2,2) = (TAE(2,2)/(Ncalcs - 1.0_pr))**0.5_pr
              TAE(2,3) = (TAE(2,3)/(Ncalcs - 1.0_pr))**0.5_pr
              TAE(2,4) = (TAE(2,4)/(Ncalcs - 1.0_pr))**0.5_pr
              TAE(2,1) = TAE(2,1)/(Ncalcs**0.5_pr)
              TAE(2,2) = TAE(2,2)/(Ncalcs**0.5_pr)
              TAE(2,3) = TAE(2,3)/(Ncalcs**0.5_pr)
              TAE(2,4) = TAE(2,4)/(Ncalcs**0.5_pr)
              Write(97, '(1X, A10, 4(A26))') '#         ', &
                  'Kinetic', 'Exchange', 'Total', 'Madelung'
              Write(97, '(1X, A10, 4(F26.16))') 'Average:', TAE(1,:)
              Write(97, '(1X, A10, 4(F26.16))') 'Error:', TAE(2,:)
          End If

        End Subroutine WriteTwistData


        Subroutine SuccessfulDeallocateAndStop(T2aaaa,T2abab,T2abba,&
                                               X2aaaa,X2abab,X2abba)

          Real (Kind=pr), Allocatable, Intent(InOut) :: T2aaaa(:,:,:), T2abab(:,:,:), &
                                                        T2abba(:,:,:), X2aaaa(:,:,:), &
                                                        X2abab(:,:,:), X2abba(:,:,:)

          Integer :: IAlloc(6)

          IAlloc = 0

          DeAllocate(T2aaaa, Stat=IAlloc(1))
          DeAllocate(T2abab, Stat=IAlloc(2))
          DeAllocate(T2abba, Stat=IAlloc(3))
          DeAllocate(X2aaaa, Stat=IAlloc(4))
          DeAllocate(X2abab, Stat=IAlloc(5))
          DeAllocate(X2abba, Stat=IAlloc(6))

          If (Any(IAlloc /= 0)) Then
              Error Stop "ERROR: Could not deallocate in main"
          Else
              Stop
          End If

        End Subroutine SuccessfulDeallocateAndStop

      End Program TermedCC

