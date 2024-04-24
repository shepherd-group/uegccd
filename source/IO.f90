! TODO - TNM: Add a input to UEGCCD-Dev for controlling the convergence
! criterial of CCD calculation (e.g like in VASP)
! TODO - WZV: Add an rng seed read in parameter.

   Module IO

   Use Precision
   Implicit None
   Private
   Public  :: ReadInput

   Contains

      Subroutine ReadInput(UEGInfo)
      Use Types, only: UegInfoType
      Implicit None
      Type (UEGInfoType), intent(Out) :: UEGInfo
      Integer, Parameter    :: NParams = 51
      Integer, Parameter    :: LName   = 19
      Integer, Parameter    :: LLine   = 79
      Logical               :: Error, Exists
      Character (len=5)     :: FormatString
      Character (len=LLine) :: Line, KeyWord, Value
      Character (len=LName) :: ParamName(NParams)
      Logical               :: SetOnce(NParams), SetTwice(NParams)
      Integer               :: I, ExStatus, LineNumber
      Logical               :: EndInput

!=====================================================================!
!  This is a pretty complicated subroutine for me, as I know little   !
!  about string-handling.  But here's what everything does.           !
!                                                                     !
!  We open Input and read it one line at a time.  Lines are assumed   !
!  to be either of the general form                                   !
!       Keyword: Value                                                !
!       Keyword:                                                      !
!  We want each of the possible keywords to be set exactly once, so   !
!  we keep track of which are set and of which are set twice.  We     !
!  also check to see if the user tried to set an invalid keyword,     !
!  telling him on which line he tried to do so.                       !
!                                                                     !
!  The parameters NParams, LName, and LLine mean:                     !
!     NParams: The number of keywords we expect to be set.            !
!     LName:   The length of the string naming each keyword.          !
!     LLine:   The length of a line.                                  !
!  Other variables are:                                               !
!     Charge:        The total charge of the system, read in.         !
!     Multiplicity:  The multiplicity of the system, read in.         !
!     CorFunc:       The correlation functional to be used, read in.  !
!     CorPot:        The correlation potential to be used, read in.   !
!     ExFunc:        The exchange functional to be used, read in.     !
!     ExPot:         The exchange potential to be used, read in.      !
!     Exists:        .True. if Input file exists.                     !
!     Error:         .True. if there is an error in the routine.      !
!     FormatString:  If LLine = xyz, FormatString = (Axyz).           !
!     Line:          The line we read in.                             !
!     Keyword:       The keyword we read in.                          !
!     Value:         The value corresponding to that keyword.         !
!     ParamName:     Array of acceptable keyword names.               !
!     SetOnce:       Array, .True. when the right keyword is set.     !
!     SetTwice:      Array, .True. if the keyword is multiply set.    !
!                                                                     !
!  After creating FormatString, we set the defaults for charge,       !
!  multiplicity, and Kohn-Sham stuff.                                 !
!=====================================================================!

      ExStatus = 0

      Write(6,*) 'Reading Input...'
      Write(FormatString,'(a2,i2,a1)') '(a',LLine,')'

!============================================================!
!  Initialize the name list, and the Input checking arrays.  !
!============================================================!

      SetOnce  = .False.
      SetTwice = .False.
      ParamName = (/'# Electrons        ',    &
                    'Momentum Cutoff    ',    &
                    'rSMin              ',    &
                    'rSMax              ',    &
                    '# rS Points        ',    &
                    'Safe ERI           ',    &
                    'Do Rings           ',    &
                    'Do XRings          ',    &
                    'Do Ladders         ',    &
                    'Do Mosaics         ',    &
                    'Rings Range        ',    &
                    'XRings Range       ',    &
                    'Ladders Range      ',    &
                    'Mosaics Range      ',    &
                    'DriverDir Range    ',    &  
                    'DriverEx Range     ',    &   ! this implicitly sets the driver
                    'Energy Range       ',    &   ! for the energy expression
                    'LinRings Range     ',    &   ! for the linear rings and cross rings
                    'QuadRings Range    ',    &   ! for the quadratic rings and cross rings
                    'DRings Range       ',    &   ! for the dRPA terms
                    'ExRings Range      ',    &   ! for the RPAX terms
                    'LinLadd Range      ',    &   ! for the linear ladders
                    'QuadLadd Range     ',    &   ! for the quadratic ladders
                    'DLadders Range     ',    &   ! by analogy to rings 
                    'ExLadders Range    ',    &   ! by analogy to rings 
                    'CoreFactor         ',    &   ! new CAS parameter 
                    'TACalcN            ',    &   ! number of twist angles for TA and cTA
                    'DoCalcOnlyTAHF     ',    &   ! True if we run TAHF only and stop
                    'DoCalcTACCD        ',    &   ! calculates CCD during TA
                    'DodRPASOSEX        ',    &   ! calculates dRPA+SOSEX instead of CCD
                    'DoTruncCoulombHF   ',    &   ! truncates coulomb intergrals in the HF
                    'DoTruncCoulombAll  ',    &   ! truncates all coulomb integrals
                    'DoSphericalTrunc   ',    &   ! Use the spherical truncation for X(G) only.
                    'DoSphericalTruncV  ',    &   ! Use the spherical truncation for V(G) as well.
                    'DoKS               ',    &   ! only KE is included in eigenvalues
                    'MadelungFactor     ',    &   ! multiplies the madelung constant by this number
                    'DoSingleCalc       ',    &   ! 
                    'DoTwistx           ',    &
                    'DoTwisty           ',    &
                    'DoTwistz           ',    &
                    'Gap                ',    &   ! opens a gap in the Fock matrix
                    'fcut               ',    &   ! 
                    'corefcut           ',    &   ! 
                    'fcut2ecut          ',    &   ! this uses fcut to set ecut
                    'GapII              ',    &   ! opens a gap in the ERI energy expression
                    'DoSFCalcMP2        ',    &   ! Structure factor for mp2  
                    'DoSFCalcCCD        ',    &   ! Structure factor for CCD
                    'DoSkipTA           ',    &   ! Do not do TA at all
                    'DoOnlyMP2Grid      ',    &   ! This just calculates the MP2 grid for ML project
                    'UseFindTol         ',    &   ! If true use the if statements in FindTol, if F ignore.
                    'EndInput           '/)       ! end of input 

!================================================================!
!  Assuming the Input file exists, open it and get ready to go.  !
!================================================================!

      Inquire(File='Input',Exist=Exists)
      If(.not. Exists) Stop 'No Input file'
      Open(4,File='Input',Status='Old')
      LineNumber = 0

      Do
       LineNumber = LineNumber + 1
       Read(4,FormatString,End=10) Line
       Call ParseLine(Line,KeyWord,Value)

       Do I = 1,NParams        !  Input debugging checks...
        If(Trim(ParamName(I)) == Trim(AdjustL(KeyWord))) Then
         If(SetOnce(I))  SetTwice(I) = .True.
         If(SetTwice(I)) Write(6,1020) Trim(ParamName(I))
         SetOnce(I) = .True.
        End If
       End Do

       Select Case (Trim(AdjustL(Keyword)))     !  Setting the variables...
        Case ('SKIP')

        Case ('# Electrons')
         Read(Value,*) UEGInfo%NElectron
         UEGInfo%NOcc = UEGInfo%NElectron/2
         If(UEGInfo%NElectron < 0)         Stop 'Must have positive number of electrons'
         If(Mod(UEGInfo%NElectron,2) == 1) Stop 'Must have closed shell'

        Case ('Momentum Cutoff')
         Read(Value,*) UEGInfo%MaxKPoint

        Case ('rSMin')
         Read(Value,*) UEGInfo%rSMin

        Case ('rSMax')
         Read(Value,*) UEGInfo%rSMax

        Case ('# rS Points')
         Read(Value,*) UEGInfo%NumRSPoints
         If(UEGInfo%NumRSPoints > 1)         Stop 'Cannot have more than one rs point'

        Case ('Safe ERI')
         Read(Value,*) UEGInfo%SafeERI

        Case ('Do Rings')
         Read(Value,*) UEGInfo%DoRing

        Case ('Do XRings')
         Read(Value,*) UEGInfo%DoXRing

        Case ('Do Ladders')
         Read(Value,*) UEGInfo%DoLadder

        Case ('Do Mosaics')
         Read(Value,*) UEGInfo%DoMosaic

        Case ('Rings Range')
         Read(Value,*) UEGInfo%IRangeRing

        Case ('XRings Range')
         Read(Value,*) UEGInfo%IRangeXRing

        Case ('Ladders Range')
         Read(Value,*) UEGInfo%IRangeLadder

        Case ('Mosaics Range')
         Read(Value,*) UEGInfo%IRangeMosaic

        Case ('DriverDir Range') ! this implicitly sets the driver
         Read(Value,*) UEGInfo%IRangeDriverDirect

        Case ('DriverEx Range') ! this implicitly sets the driver
         Read(Value,*) UEGInfo%IRangeDriverExchange

        Case ('Energy Range') ! for the energy expression
         Read(Value,*) UEGInfo%IRangeEnergy

        Case ('LinRings Range') ! for the linear rings and cross rings
         Read(Value,*) UEGInfo%IRangeLinRings

        Case ('QuadRings Range') ! for the quadratic rings and cross rings
         Read(Value,*) UEGInfo%IRangeQuadRings
        
        Case ('DRings Range') ! for the dRPA terms
         Read(Value,*) UEGInfo%IRangeDirectRings

        Case ('ExRings Range') ! for the RPAX terms
         Read(Value,*) UEGInfo%IRangeExchangeRings

        Case ('LinLadd Range') ! for the linear ladders
         Read(Value,*) UEGInfo%IRangeLinLadders

        Case ('QuadLadd Range') ! for the quadratic ladders
         Read(Value,*) UEGInfo%IRangeQuadLadders
        
        Case ('DLadders Range') ! by analogy to rings 
         Read(Value,*) UEGInfo%IRangeDirectLadders

        Case ('ExLadders Range') ! by analogy to rings
         Read(Value,*) UEGInfo%IRangeExchangeLadders
        
        Case ('NTwist') ! twist angle
         Read(Value,*) UEGInfo%NTwist

        Case ('CoreFactor') ! for CAS
         Read(Value,*) UEGInfo%CoreFactor
         UEGInfo%CoreN=int(UEGInfo%NElectron/2*UEGInfo%CoreFactor)
         write(61,*) UEGInfo%CoreN*2,"core electrons"
         write(61,*) UEGInfo%NElectron-UEGInfo%CoreN*2,"active electrons"
        
        Case ('DoCalcTACCD') ! calculates CCD during TA
         Read(Value,*) UEGInfo%DoCalcTACCD
                    
        Case ('TACalcN') ! number of twists
         Read(Value,*) UEGInfo%TACalcN

        Case ('DoCalcOnlyTAHF')
         Read(Value,*) UEGInfo%DoCalcOnlyTAHF
                    
        Case ('DodRPASOSEX')
         Read(Value,*) UEGInfo%DodRPASOSEX
                    
        Case ('DoTruncCoulombHF') 
         Read(Value,*) UEGInfo%DoTruncCoulombHF
         If(UEGInfo%DoTruncCoulombHF)         Stop 'Not Implemented'
                    
        Case ('DoTruncCoulombAll') 
         Read(Value,*) UEGInfo%DoTruncCoulombAll
         If(UEGInfo%DoTruncCoulombAll)         Stop 'Not Implemented'
                    
        Case ('DoSphericalTrunc') 
         Read(Value,*) UEGInfo%DoSphericalTrunc
         If (UEGInfo%DoSphericalTrunc) Then
            Write(60,*) 'WARNING: Using spherical truncation on X(G) terms, &
                         this is an experimental feature!'
         End If

        Case ('DoSphericalTruncV') 
         Read(Value,*) UEGInfo%DoSphericalTruncV
         If (UEGInfo%DoSphericalTrunc) Then
            Write(60,*) 'WARNING: Spherical truncation on V(G) sets SafeERI to &
                        .false. and is an experimental feature!'
            UEGInfo%SafeERI = .false.
         EndIf

        Case ('DoKS')
         Read(Value,*) UEGInfo%DoKS
         !If(DoKS)         Stop 'Not Implemented'
         If(UEGInfo%DoKS) then
             UEGInfo%exchangefactor=0.0_pr
         else
             UEGInfo%exchangefactor=1.0_pr
         endif
        
        Case ('MadelungFactor')
         Read(Value,*) UEGInfo%madelungfactor
                    
        Case ('DoSingleCalc') 
         Read(Value,*) UEGInfo%DoSingleCalc
        
        Case ('DoTwistx')
         Read(Value,*) UEGInfo%DoTwistx

        Case ('DoTwisty')
         Read(Value,*) UEGInfo%DoTwisty

        Case ('DoTwistz')
         Read(Value,*) UEGInfo%DoTwistz

        Case ('DoSFCalcMP2') 
            Read(Value,*) UEGInfo%DoSFCalcMP2
        
        Case ('DoSFCalcCCD') 
            Read(Value,*) UEGInfo%DoSFCalcCCD
        
        Case ('DoSkipTA') 
         Read(Value,*) UEGInfo%DoSkipTA
        
        Case ('DoOnlyMP2Grid') 
            Read(Value,*) UEGInfo%DoOnlyMP2Grid ! kills calc after one MP2 calc
        
        Case ('UseFindTol') 
            Read(Value,*) UEGInfo%UseFindTol ! kills calc after one MP2 calc
        
        Case ('Gap') ! for an insulating gas
         Read(Value,*) UEGInfo%Gap
                    
        Case ('fcut') ! for a basis set that scales with the electron number
         Read(Value,*) UEGInfo%fcut
        
        Case ('corefcut') ! for an insulating gas
         Read(Value,*) UEGInfo%corefcut
        
        Case ('fcut2ecut') ! to automatically set the ecut from the fcut
         Read(Value,*) UEGInfo%fcut2ecut
         if (UEGInfo%fcut2ecut) then
            UEGInfo%Nfermi=UEGInfo%NElectron**(2.0_pr/3.0_pr)/4.0_pr ! approximately divided by 4
            UEGInfo%fcutMultiplier=UEGInfo%fcut
            if (UEGInfo%fcut < 2.0) then
               UEGInfo%fcutMultiplier=2.0
            endif
            UEGInfo%MaxKPoint=UEGInfo%Nfermi*UEGInfo%fcutMultiplier*2.0_pr
            if (UEGInfo%Nfermi < 10) then
               UEGInfo%MaxKPoint=10*UEGInfo%fcutMultiplier*2.0_pr
            endif
        endif
        
        Case ('GapII') ! for an insulating gas
         Read(Value,*) UEGInfo%GapII
                    
        Case ('EndInput')
         Read(Value,*) EndInput
         If(.not.EndInput)         Stop 'EndInput must be T'
                    
        Case Default   !  Unknown keyword
         ExStatus=1
         write(6,*) Value
         Exit

       End Select
      End Do
10    Close(4,Status='Keep')


!==================================================================!
!  If the user tried to set an undefined keyword, tell him where.  !
!  Then check if the required keywords are set.                    !
!  Then check if any are multiply set.                             !
!  Finally, if there was an error, we stop.                        !
!                                                                  !
!  Since parameters 1:2 have defaults, we tell the code as much.   !
!==================================================================!

      SetOnce(1:2) = .True.
      If(ExStatus == 1) Write(6,1000) LineNumber
      Do I = 1,NParams
       If(.not. SetOnce(I)) Write(6,1010) Trim(ParamName(I))
      End Do
      Error = ExStatus == 1 .or. Any(SetTwice) .or. .not. All(SetOnce)
      If(Error) Stop 'Error in Input'


!==================!
!  Over to James.  !
!==================!

1000  Format('Error: Line number ',I4,' of Input not recognized; ',  &
             'subsequent lines not read.')
1010  Format('Error: Parameter ',A,' not set.')
1020  Format('Error: Parameter ',A,' multiply set.')

      Return
      End Subroutine ReadInput






      Subroutine ParseLine(Line,KeyWord,Value)
      Implicit None
      Character (Len=*),         Intent(In)  :: Line
      Character (Len=Len(Line)), Intent(Out) :: KeyWord, Value
      Character (Len=Len(Line))              :: WorkingLine
      Integer :: I, CommentPos, ColonPos

!==================================================!
!  Swiped from Paul, reads keywords and values.    !
!--------------------------------------------------!
!  Check for comments and convert them to spaces,  !
!  then turns control characters into spaces       !
!==================================================!

      WorkingLine = Line
      CommentPos  = Index(WorkingLine,'!')
      If(CommentPos /= 0) WorkingLine = WorkingLine(1:CommentPos-1)
      Do I = 1,Len(Line)
       If(IAChar(WorkingLine(I:I)) <= 31) WorkingLine(I:I) = ' '
      End Do


!=================================================================!
!  If we have a blank line at this point, we skip this line.      !
!  Otherwise, we presumably have Keyword: Value and return them.  !
!=================================================================!

      If(Len_Trim(WorkingLine) == 0) Then
       KeyWord = 'SKIP'
       Return
      End If
      ColonPos = Index(WorkingLine,':',back=.False.)
      KeyWord = WorkingLine(1:ColonPos-1)
      Value   = WorkingLine(ColonPos+1:)

      Return
      End Subroutine ParseLine

   End Module IO

