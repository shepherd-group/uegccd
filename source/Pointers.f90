Module Pointers

Use Precision, Only: pr

Implicit None

Abstract Interface

    Pure Function AIERI(HEGData,UEGInfo,A,B,I,J,DummyFlag) Result(Intgrl)
        Use Types, Only: HEGDataType, UEGInfoType
        Import :: pr
        Type (HEGDataType), Intent(In) :: HEGData
        Type (UEGInfoType), Intent(In) :: UEGInfo
        Integer, Intent(In) :: I, J, A, B
        Integer, Intent(In), Optional :: DummyFlag
        Real (Kind=pr) :: Intgrl
    End Function AIERI

End Interface

Procedure(AIERI), Pointer :: ERI => Null()

End Module Pointers
