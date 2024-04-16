
   Module MP2

   Use Precision
   Use Constants

   Implicit None

   Contains

      Subroutine DrvMBPT(HEGData,UEGInfo,T2aaaa,T2abab,T2abba,Conn)

      Use HEG, Only: FindIndex, SendForAllocation, dK, dK2, dK2scaled
      Use Types, Only: HEGDataType, UEGInfoType
      Use Pointers

      Type (HEGDataType), Intent(InOut) :: HEGData
      Type (UEGInfoType), Intent(In) :: UEGInfo
      Real (Kind=pr), Intent(Out) :: T2aaaa(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                     T2abab(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO), &
                                     T2abba(UEGInfo%NOcc,UEGInfo%NOcc,UEGInfo%NOcc+1:UEGInfo%NAO)

      Integer :: I, J, A, B
      Real (Kind=pr) :: V_ijab, V_ijba, Denom
      ! Connectivity calculation
      Integer :: ConnectivityLimit=400
      Integer :: Connectivity(400)
      Integer :: dK2_ijab, dK2_ijba
      Logical :: InCore
      Integer, Allocatable, Optional :: Conn(:)
      Real (pr), Allocatable :: SfSum(:,:,:), SfG(:,:,:), SfV(:,:,:)
      integer :: SfLength
      Integer :: this_qvec(3)
      Real (pr) :: this_q2

!==========================================!
!  This does the MP2 calculation.  Glory!  !
!==========================================!

! TODO:
! - Calculate the connectivity <--- this first
! - Pass back the connectivity to the Main.f90 algorithm? 
! - Maybe calculate the connectivity in a separate MBPT loop? 
! - Needs to have a "mode" call, where the vector passed out is different
      
      call SendForAllocation(HEGData,SfG,SfLength)
      call SendForAllocation(HEGData,SfV,SfLength)
      call SendForAllocation(HEGData,SfSum,SfLength)

      Write(6,*) "Doing MBPT2..."
      T2abab = Zero
      T2abba = Zero
      T2aaaa = Zero
      HEGData%ECorr = Zero
      Connectivity = Zero ! connectivity calculation
      Do I = 1,UEGInfo%NOcc
      Do J = 1,UEGInfo%NOcc
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        B = FindIndex(HEGData,UEGInfo%MaxKPoint,UEGInfo%NAO,I,J,A)
        If(B <= UEGInfo%NOcc) Cycle
        V_ijab = ERI(HEGData,UEGInfo,I,J,A,B)
        V_ijba = ERI(HEGData,UEGInfo,I,J,B,A)
        If(Min(V_ijab,V_ijba) < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
        Denom = HEGData%Eigen(I) + HEGData%Eigen(J) - HEGData%Eigen(A) - HEGData%Eigen(B)
        T2abab(I,J,A) = V_ijab/Denom
        T2abba(I,J,A) = -V_ijba/Denom
        T2aaaa(I,J,A) = T2abab(I,J,A) + T2abba(I,J,A)
        HEGData%ECorr = HEGData%ECorr &
                      + F12*T2aaaa(I,J,A)*(V_ijab - V_ijba) &
                      + F12*T2abab(I,J,A)*(V_ijab) &
                      - F12*T2abba(I,J,A)*(V_ijba)
        dK2_ijab = dK2(HEGData,UEGInfo%CoreN,I,J,A,B,InCore) ! TODO need this function
        dK2_ijba = dK2(HEGData,UEGInfo%CoreN,I,J,B,A,InCore) 
        If ((dK2_ijab > ConnectivityLimit) .or. (dK2_ijba > ConnectivityLimit)) Then
          Print *, "Connectivity greater than current ConnectivityLimit"
          Print *, "Increase ConnectivityLimit"
          Stop "Program stopped in MP2.f90"
        End If
        If (.not. InCore) then
            Connectivity(dK2_ijba)=Connectivity(dK2_ijba)+1
            Connectivity(dK2_ijab)=Connectivity(dK2_ijab)+1
        End If
        this_qvec=dK(HEGData,UEGInfo%CoreN,a,b,i,j,inCore)
        SfSum(this_qvec(1),this_qvec(2),this_qvec(3))=SfSum(this_qvec(1),this_qvec(2),this_qvec(3))+(T2aaaa(I,J,A)+ T2abab(I,J,A))
        this_q2=dK2scaled(HEGData,UEGInfo%CoreN,a,b,i,j,inCore)
        SfG(this_qvec(1),this_qvec(2),this_qvec(3))=this_q2**0.5_pr
        SfV(this_qvec(1),this_qvec(2),this_qvec(3))=V_ijab
        this_qvec=dK(HEGData,UEGInfo%CoreN,a,b,j,i,inCore)
        SfSum(this_qvec(1),this_qvec(2),this_qvec(3))=SfSum(this_qvec(1),this_qvec(2),this_qvec(3))-(T2aaaa(I,J,A) + T2abba(I,J,A))
        this_q2=dK2scaled(HEGData,UEGInfo%CoreN,a,b,j,i,inCore)
        SfG(this_qvec(1),this_qvec(2),this_qvec(3))=this_q2**0.5_pr
        SfV(this_qvec(1),this_qvec(2),this_qvec(3))=V_ijba
      End Do
      End Do
      End Do
      
      If (Present(Conn)) Then
          Allocate(Conn(ConnectivityLimit))
          Conn=Connectivity
      End If

      write(81,*) "--------"
      write(81,*) "Starting a new structure factor"
      write(81,*) "--------" 

      Do I = -SfLength,SfLength
      Do J = -SfLength,SfLength
      Do A = -SfLength,SfLength
        if (SfG(I,J,A) > 0.0_pr) then
            write(81,*) I, J, A, SfG(I,J,A), SfSum(I,J,A), SfV(I,J,A)
        endif
      End Do
      End Do
      End Do

      deallocate(SfSum,SfG,SfV)


! Write to Output and we're done!
      Open(7,File='Output',Position="Append")
      Write(7,1000)
      Write(7,1005) HEGData%EHF
      Write(7,1040)
      Write(7,1030)
      Write(7,1040)
      Do I = 1,UEGInfo%NOcc
        Write(7,1050) HEGData%Eigen(I)
      End Do
      Write(7,1060)
      Do A = UEGInfo%NOcc+1,UEGInfo%NAO
        Write(7,1050) HEGData%Eigen(A)
      End Do
      Write(7,1000)
      Write(7,1010) HEGData%ECorr
      Write(7,1020) HEGData%ECorr + HEGData%EHF
      Write(7,1000)
      Close(7)


1000  Format(14x,'**************************************************')
1005  Format(14x,'* E(SCF) is ',12x,F18.12,' a.u.',2x,'*')
1010  Format(14x,'* E(2) is ',14x,F18.12,' a.u.',2x,'*')
1020  Format(14x,'* E(SCF) + E(2) is ',5x,F18.12,' a.u.',2x,'*')
1030  Format(14X,'*                RHF Eigenvalues:                *')
1040  Format(14x,'*------------------------------------------------*')
1050  Format(14x,'*',12x,ES24.16,12x,'*')
1060  Format(14x,'*',12x,24("-"),12x,'*')

      if (UEGInfo%DoOnlyMP2Grid) then
          stop "DoOnlyMP2Grid flag is true, so calculation stopping here"
      endif 

      Return
      End Subroutine DrvMBPT

   End Module MP2

