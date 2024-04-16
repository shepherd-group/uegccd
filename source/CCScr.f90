
   Module CCScr

   Use Precision
   Implicit None

   Contains

      Subroutine SetUpCCScr(O,V,N,CCScrData)

      Use Types, only: CCScrDataType

      Type (CCScrDataType), Intent(InOut) :: CCScrData
      
      Integer, Intent(In) :: O, V, N

      Integer :: IAlloc

      Allocate(CCScrData%G2aaaa(1:O,1:O,V:N),     &
               CCScrData%G2abab(1:O,1:O,V:N),     &
               CCScrData%G2abba(1:O,1:O,V:N),     &
               CCScrData%Joo(1:O),                &
               CCScrData%Jvv(V:N),                &
               Stat=IAlloc)

      If(IAlloc /= 0) Stop 'Could not allocate in CCScr'

      Return

      End Subroutine SetUpCCScr


      Subroutine ShutDownCCScr(CCScrData)

      Use Types, only: CCScrDataType

      Type (CCScrDataType), Intent(InOut) :: CCScrData

      Integer :: IAlloc

      DeAllocate(CCScrData%G2aaaa, &
                 CCScrData%G2abab, &
                 CCScrData%G2abba, &
                 CCScrData%Joo,    &
                 CCScrData%Jvv,    &
                 Stat=IAlloc)

      If(IAlloc /= 0) Stop 'Could not deallocate in CCScr'

      Return

      End Subroutine ShutDownCCScr

   End Module CCScr

