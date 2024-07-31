Module Types

Use Precision

Implicit None

Type HEGDataType
    Integer :: Nel
    Integer :: Cell
    Real (Kind=pr) :: L
    Real (Kind=pr) :: rS
    Real (Kind=pr) :: Omega
    Real (Kind=pr) :: Madelung
    Real (Kind=pr) :: Density
    Real (Kind=pr) :: FermiWaveVector
    Real (Kind=pr) :: ScreeningDistance
    Real (Kind=pr) :: nTwist(3)
    Real (Kind=pr) :: this_nTwist(3)
    Real (Kind=pr) :: XHF
    Real (Kind=pr) :: EHF
    Real (Kind=pr) :: ECorr
    ! Allocatables
    Real (Kind=pr), Allocatable :: kVec(:,:)
    Real (Kind=pr), Allocatable :: Eigen(:)
    Integer, Allocatable :: nVec(:,:)
    Integer, Allocatable :: kPointToBasisFn(:,:,:)
End Type HEGDataType

Type UEGInfoType
    Integer :: NElectron
    Integer :: NOcc
    Integer :: NAO
    Integer :: MaxKPoint
    Integer :: CoreN
    Integer :: NumRSPoints
    Integer :: IRangeRing
    Integer :: IRangeXRing
    Integer :: IRangeLadder
    Integer :: IRangeMosaic
    Integer :: IRangeDriverDirect
    Integer :: IRangeDriverExchange
    Integer :: IRangeEnergy
    Integer :: IRangeLinRings
    Integer :: IRangeQuadRings
    Integer :: IRangeDirectRings
    Integer :: IRangeExchangeRings
    Integer :: IRangeLinLadders
    Integer :: IRangeQuadLadders
    Integer :: IRangeDirectLadders
    Integer :: IRangeExchangeLadders
    Integer :: TACalcN
    Real (Kind=pr) :: rS
    Real (Kind=pr) :: rSMin
    Real (Kind=pr) :: rSMax
    Real (Kind=pr) :: NTwist(3)
    Real (Kind=pr) :: SinglePointCalcTwist
    Real (Kind=pr) :: CoreFactor
    Real (Kind=pr) :: fcutMultiplier
    Real (Kind=pr) :: Nfermi
    Real (Kind=pr) :: exchangefactor
    Real (Kind=pr) :: madelungfactor
    Real (Kind=pr) :: Gap
    Real (Kind=pr) :: GapII
    Real (Kind=pr) :: fcut
    Real (Kind=pr) :: corefcut
    Logical :: SafeERI
    Logical :: DoOnlyMP2Grid
    Logical :: fcut2ecut
    Logical :: DoRing
    Logical :: DoXRing
    Logical :: DoLadder
    Logical :: DoMosaic
    Logical :: DoCalcTACCD
    Logical :: DodRPASOSEX
    Logical :: DoTruncCoulombHF
    Logical :: DoTruncCoulombAll
    Logical :: DoKS
    Logical :: DoMeanPotential
    Logical :: DoSingleCalc
    Logical :: DoSFCalcMP2
    Logical :: DoSFCalcCCD
    Logical :: DoSkipTA
End Type UEGInfoType

Type BasisSet
    Integer :: sorted_index
    Real (Kind=pr) :: k2
    Integer :: n2
    Integer :: n(3)
    Integer :: ms
    Real (Kind=pr) :: ntwisted(3)
    Real (Kind=pr) :: ntwisted2
End Type BasisSet

Type TwistAngleData
    Real (Kind=pr) :: ntwist(3)
    Real (Kind=pr) :: hf
    Real (Kind=pr) :: mp2
    Real (Kind=pr) :: ccd
    Integer, Allocatable :: connectivity(:)
    Real (Kind=pr), Allocatable :: Eigen(:)
    Real (Kind=pr) :: ConnectivityDiff(400)
End Type TwistAngleData

Type TwistAngleCalculation
    Real (Kind=pr) :: NTwistSum(3)
    Real (Kind=pr) :: HfSum
    Real (Kind=pr) :: Mp2Sum
    Real (Kind=pr) :: CcdSum
    Real (Kind=pr) :: NTwistSum2(3)
    Real (Kind=pr) :: HfSum2
    Real (Kind=pr) :: Mp2Sum2
    Real (Kind=pr) :: CcdSum2
    Real (Kind=pr) :: NTwist(3)
    Real (Kind=pr) :: Hf
    Real (Kind=pr) :: Mp2
    Real (Kind=pr) :: Ccd
    Real (Kind=pr) :: NTwist2(3)
    Real (Kind=pr) :: Hf2
    Real (Kind=pr) :: Mp22
    Real (Kind=pr) :: Ccd2
    Real (Kind=pr) :: NTwistErr(3)
    Real (Kind=pr) :: HfErr
    Real (Kind=pr) :: Mp2Err
    Real (Kind=pr) :: CcdErr
    Integer :: NumTwists
    Integer :: ConnectivitySum(400)
    Real (Kind=pr) :: ConnectivityAverage(400)
    Real (Kind=pr), allocatable :: EigenSum(:)
    Real (Kind=pr), allocatable :: EigenAverage(:)
    Integer :: EigenLength
    Real (Kind=pr) :: MinConnectivityDiff2
    Real (Kind=pr) :: SpecialTwistAngle(3)
End Type TwistAngleCalculation

Type CCScrDataType 
! Arrays for the subroutines in CCSrc.f90
   Real (Kind=pr), Allocatable :: G2aaaa(:,:,:)
   Real (Kind=pr), Allocatable :: G2abab(:,:,:)
   Real (Kind=pr), Allocatable :: G2abba(:,:,:)
   Real (Kind=pr), Allocatable :: Joo(:)
   Real (Kind=pr), Allocatable :: Jvv(:)
End Type CCScrDataType

Type DIISCCDataType
! Arrays for the subroutines in DIISCC.f90
   Integer :: NDIIS
   Real (Kind=pr), Allocatable :: OldT2abab(:,:,:,:)
   Real (Kind=pr), Allocatable :: R2abab(:,:,:,:)
   Real (Kind=pr), Allocatable :: OldT2abba(:,:,:,:)
   Real (Kind=pr), Allocatable :: R2abba(:,:,:,:)
End Type DIISCCDataType

End Module Types
