module heg 

use Precision
use Constants

implicit none

private

public :: init_HEG, FindIndex, change_rs, findtol, ChangeTwist, &
          GenerateTwist, dK2, dK2scaled, dK, SendForAllocation, &
          CasedERI, NoCaseERI, CoreERI, GapERI, BaseERI

contains

    subroutine init_HEG(UEGInfo,HEGData)
        
        Use Types, only: UEGInfoType, HEGDataType

        Type (UEGInfoType), Intent(InOut) :: UEGInfo
        Type (HEGDataType), Intent(InOut) :: HEGData

        ! prints a table
        call PrintShells(10) 

        HEGData%nTwist = UEGInfo%NTwist
        HEGData%Nel = UEGInfo%NElectron
        HEGData%rS = 1.0_pr
        
        HEGData%L=HEGData%rS*(4.0_pr/3.0_pr*pi* HEGData%Nel)**(1.0_pr/3.0_pr)
        HEGData%Omega=HEGData%L**3
        HEGData%Madelung=CalcMadelung(HEGData%Omega,HEGData%L)*UEGInfo%madelungfactor
        ! The condition on cell is: cell > sqrt(ecut)
        HEGData%Cell=int(sqrt(real(UEGInfo%MaxKPoint,dp)))+3 ! this looks really dodgy to me 
        write(60,*) "Madelung constant calculated as", HEGData%Madelung
        HEGData%Density=HEGData%Nel/HEGData%L**3.0_pr
        HEGData%FermiWaveVector=(3.0_pr*pi**2.0_pr*HEGData%Density)**(1.0_pr/3.0_pr)
        HEGData%ScreeningDistance=(4.0_dp/pi*HEGData%FermiWaveVector)**0.5_dp

        ! makes and deallocates G1, makes eigenvector array
        call SetupBasis(HEGData,UEGInfo)

    end subroutine


    subroutine PrintShells(kc) 
    ! kc is an interger input
    ! genuinely kc in integer coordinates
    ! corresponds to "cell" keyword

        integer :: l1,l2,l3
        integer :: kc
        integer :: bins
        integer :: energy
        integer, allocatable :: histogram(:)
        integer :: culm

        bins=kc**2
        allocate(histogram(0:bins))
        histogram = 0

        do l1=-kc,kc
            do l2=-kc,kc
                do l3=-kc,kc
                    energy=l1**2+l2**2+l3**2
                    if (energy.lt.kc**2) then
                        histogram(energy)=histogram(energy)+1
                    endif
                enddo
            enddo
        enddo
        culm=0
        write(60,*) ""
        write(60,*) "--------------------"
        write(60,*) "Basis sets available"
        write(60,*) "--------------------"
        write(60,*) "# ecut k-points M"

        do l1=0,kc**2
            if (histogram(l1).ne.0) then
                culm=culm+histogram(l1)
                write(60,*) l1, culm, culm*2
            endif
        enddo

    end subroutine


    function CalcMadelung(Omega,L) result(Madelung)

        Real (Kind=pr), Intent(In) :: Omega, L

        Integer :: l1,l2,l3,l4,n2
        Real (Kind=pr) :: Madelung,kappa,k2,ek2,recipsum2,modr,er2,realsum2,term2,term4

        write(60,*) ""
        write(60,*) "-----------------------------------------------------------"
        write(60,*) "Calculating Madelung Constant - Fraser et al. PRB 53 4 1814"
        write(60,*) "-----------------------------------------------------------"
        kappa=2.8_pr/Omega**(1.0_pr/3.0_pr)
        write(60,*) "kappa taken from CASINO manual to be", kappa

        term2=-pi/(kappa**2.0_pr*L**3.0_pr)
        write(60,*) term2, "term2"
        term4=-2.0_pr*kappa/sqrt(pi)
        write(60,*) term4, "term4"

        recipsum2=0.0_pr
        do l4=1,10
            recipsum2=0.0_pr
            do l1=-l4,l4
                do l2=-l4,l4
                    do l3=-l4,l4
                        n2=l1**2+l2**2+l3**2
                        k2=(1.0_pr/L**2.0_pr)*(l1**2+l2**2+l3**2)
                        ek2=(1.0_pr/L**3.0_pr)*(1.0_pr/(pi*k2))*exp(-pi**2.0_pr*k2/kappa**2.0_pr)
                        if (n2.ne.0) then
                            !write(60,*) k2,ek2 ! for testing
                            recipsum2=recipsum2+ek2
                        endif
                    enddo
                enddo
            enddo
            write(60,*) l4,recipsum2
        enddo
        write(60,*) "reciprocal space", recipsum2
        
        realsum2=0.0_pr
        do l4=1,10
            realsum2=0.0_pr
            do l1=-l4,l4
                do l2=-l4,l4
                    do l3=-l4,l4
                        modr=L*sqrt(real((l1**2+l2**2+l3**2),dp))
                        if (modr.ne.0.0_pr) then
                            er2=erfc(kappa*modr)/modr
                            realsum2=realsum2+er2
                            !write(60,*) modr*0.5_pr,er2 ! for testing
                        endif
                    enddo
                enddo
            enddo
            write(60,*) l4,realsum2
        enddo
        write(60,*) "real space", realsum2
        Madelung=realsum2+recipsum2+term2+term4

    end function

    subroutine SetupBasis(HEGData,UEGInfo)

        Use Types, only: HEGDataType, UEGInfoType, BasisSet

        Type (HEGDataType), Intent(InOut) :: HEGData
        Type (UEGInfoType), Intent(InOut) :: UEGInfo

        Type (BasisSet), allocatable :: G1(:)
        Type (BasisSet) :: G1_temp

        Logical :: swapped
        Integer :: l1, l2, l3, M, Mfcut, offset, Mcorefcut, nBasis
        Real (Kind=pr) :: Mtarget, qvec(3)

        write(60,*) ""
        write(60,*) "-----------------------------------------------------------"
        write(60,*) "Setting up HEG basis set"
        write(60,*) "-----------------------------------------------------------"
        allocate(G1((2*HEGData%cell+1)**3))
        write(60,*) "G1 memory allocated for", (2*HEGData%cell+1)**3, "basis functions"

        ! Generates basis set: all plane waves which have an less than the energy 
        ! cutoff ecut where E=(1/2)(n^2) i.e. scaled units

        ! Plane wave bais set:
        ! \psi_j(r,\sigma) = \sqrt(1/\Omega)*exp(i*k_i.r) \delta_{\sigma_i,\sigma}
        ! where real-space cell volume \Omega = L**3 and 
        ! the reciprocal lattice vectors are k = ((2*pi)/L)(n,m,l) 
        ! n,m,l=0,+/-1,+/-2,... are integers

        M=0
        do l1=-HEGData%cell,HEGData%cell
            do l2=-HEGData%cell,HEGData%cell
                do l3=-HEGData%cell,HEGData%cell
                    M=M+1
                    G1(M)%n(1)=l1
                    G1(M)%n(2)=l2
                    G1(M)%n(3)=l3
                    G1(M)%n2=l1**2+l2**2+l3**2
                    G1(M)%k2=(l1**2+l2**2+l3**2)*(2*pi/HEGData%L)**2
                    G1(M)%ntwisted(1)=l1+HEGData%ntwist(1)
                    G1(M)%ntwisted(2)=l2+HEGData%ntwist(2)
                    G1(M)%ntwisted(3)=l3+HEGData%ntwist(3)
                    G1(M)%ntwisted2=G1(M)%ntwisted(1)**2+G1(M)%ntwisted(2)**2+G1(M)%ntwisted(3)**2
                    if (G1(M)%n2 .gt. UEGInfo%MaxKPoint) then 
                        M=M-1 ! By construction G1(M+1) exists but is never used (caution!)
                    endif
                enddo
            enddo
        enddo

        nBasis=M

        write(60,*) "Actual basis set size", nBasis

        ! Sorting basis set according to ascending energy
        do
            swapped=.false.
            do l1=2,nBasis
                if (G1(l1)%ntwisted2.lt.G1(l1-1)%ntwisted2) then
                    G1_temp=G1(l1-1)
                    G1(l1-1)=G1(l1)
                    G1(l1)=G1_temp
                    swapped=.true.
                endif
            enddo
            if (.not.swapped) exit
        enddo
        
        ! Create a look-up table from the k-points to the orbital indices
        allocate(HEGData%kPointToBasisFn(-HEGData%cell:HEGData%cell, &
                                         -HEGData%cell:HEGData%cell, &
                                         -HEGData%cell:HEGData%cell))

        Mfcut=int(UEGInfo%fcut**(3.0_pr/2.0_pr)*HEGData%Nel/2)
        Mtarget=UEGInfo%fcut**(3.0_pr/2.0_pr)*HEGData%Nel/2
        ! now, correct if it's artificially too low.
        if (abs(Mfcut-Mtarget) > abs(Mfcut+1-Mtarget)) then
            Mfcut=Mfcut+1
        endif

        if (Mfcut < nBasis) then
            nbasis=Mfcut
            M=Mfcut
        endif

        if (abs(FcutToBasisFn(UEGInfo%fcut,HEGData%Nel)-Mfcut) > 0) then
            write(6,*) FcutToBasisFn(UEGInfo%fcut, HEGData%Nel), UEGInfo%fcut, Mfcut
            stop "Error in FcutToBasisFn during setup_basis" 
        endif 

        write(60,*) "Basis set size after fcut", nBasis
        
        ! "Insulating" gas code below
        ! ===========================
        ! corefcut should remove orbitals that are too close to the Fermi sphere
        ! do this with a factor with respect to N, I think
        ! e.g. M=gamma N
        ! but for consistency do we apply this with an fcut cutoff?
        ! I think probably yes
        !
        ! so that means that corefcut minimum is 1
        ! 
        Mcorefcut=FcutToBasisFn(UEGInfo%corefcut,HEGData%Nel)
        if ((UEGInfo%corefcut-1.0_pr) < 1e-6) then
            ! corecut is 1.0 and this is the default
            if (Mcorefcut .ne. HEGData%Nel/2) then
                write(6,*) Mcorefcut,HEGData%Nel/2
                stop "ERROR in setup_basis: corefcut is less than 1 rather than equal to 1"
            endif
        endif
        nbasis=nbasis-Mcorefcut+HEGData%Nel/2 ! Mcorefcut contains the occupied manifold too
        M=nbasis
        UEGInfo%NAO=M
        
        write(60,*) "Basis set size after corefcut", nBasis

        ! Initiatlise arrays
        allocate(HEGData%kvec(M,3))
        allocate(HEGData%eigen(M))
        allocate(HEGData%nvec(M,3))
        HEGData%kvec=0.0_pr
        HEGData%eigen=0.0_pr
        HEGData%nvec=0

        ! Form the array of k vectors (scaled)
        ! and n vectors (integers)
        do l1=1,nbasis
            offset=0
            if (l1 > HEGData%Nel/2) then
                offset=Mcorefcut-HEGData%Nel/2 ! this is the relative offset between G1 and k
            endif
            HEGData%kvec(l1,1)=G1(l1+offset)%ntwisted(1)*(2.0_pr*pi/HEGData%L)
            HEGData%kvec(l1,2)=G1(l1+offset)%ntwisted(2)*(2.0_pr*pi/HEGData%L)
            HEGData%kvec(l1,3)=G1(l1+offset)%ntwisted(3)*(2.0_pr*pi/HEGData%L)
            HEGData%nvec(l1,1)=G1(l1+offset)%n(1) ! When writing twist angle: I think I need to leave these alone for momentum symmetry checking
            HEGData%nvec(l1,2)=G1(l1+offset)%n(2)
            HEGData%nvec(l1,3)=G1(l1+offset)%n(3)
        enddo

        ! Form eigenvector array
        do l1=1,nbasis
            HEGData%eigen(l1)=0.0_pr
            HEGData%eigen(l1)=0.5_pr*dot_product(HEGData%kvec(l1,:),HEGData%kvec(l1,:)) ! kinetic
            do l2=1,HEGData%Nel/2
                if (l2.ne.l1) then
                    ! find "q" momentum transfer
                    qvec=HEGData%kvec(l1,:)-HEGData%kvec(l2,:)
                    HEGData%eigen(l1)=HEGData%eigen(l1)-(1/HEGData%Omega)*(4.0_pr*pi/dot_product(qvec,qvec))*UEGInfo%exchangefactor
                else 
                    HEGData%eigen(l1)=HEGData%eigen(l1)+HEGData%Madelung-UEGInfo%gap ! self-interaction, madelung is a negative constant
                endif
            enddo
            write(60,*) l1, HEGData%nvec(l1,1), HEGData%nvec(l1,2), HEGData%nvec(l1,3), HEGData%eigen(l1) !TODO: edit to include new output
        enddo
        
        HEGData%kPointToBasisFn=0 ! there are some elements which aren't in the basis, these will
                          ! be returned as zero
        do l1=1,nBasis
            offset=0
            if (l1 > HEGData%Nel/2) then
                offset=Mcorefcut-HEGData%Nel/2 ! this is the relative offset between G1 and k
            endif
            HEGData%kPointToBasisFn(  G1(l1+offset)%n(1), &
                                      G1(l1+offset)%n(2), &
                                      G1(l1+offset)%n(3)     ) = l1
        enddo

        ! Print out E_HF
        call print_HF(HEGData)

        ! G1 is now an unnecessary array
        deallocate(G1)

    end subroutine SetupBasis

    function FcutToBasisFn(ThisFcut,Nel) result(FcutM)

        Real (Kind=pr), Intent(In) :: ThisFcut
        Integer, Intent(In) :: Nel

        Real (Kind=pr) :: MTarget
        Integer :: FcutM

        FcutM=int(ThisFcut**(3.0_pr/2.0_pr)*Nel/2)

        ! now, correct if it's artificially too low.
        MTarget=ThisFcut**(3.0_pr/2.0_pr)*Nel/2
        if (abs(FcutM-MTarget) > abs(FcutM+1-MTarget)) then
            FcutM=FcutM+1
        endif

    end function 

    ! This function takes integer inputs as orbital indices and finds the fourth 
    ! momentum-allowed index and returns the value in the function
    !
    ! Returns FindIndex=0 if there isn't a momentum allowed orbital in the WHOLE BASIS
    !
    ! The momentum conservation law used here is that for an excitation 
    ! ij -> ab
    ! k_i+k_j = k_a+k_b
    ! i.e. the total momentum is the same before and after the excitation
    !
    ! Inputs: m,n,l - indecies of orbitals
    ! Global data requirements: kPointToBasisFn, nvec, ecut, nBasis
    !
    Function FindIndex(HEGData,EnergyCutoff,nBasis,M,N,L) Result(K)

        Use Types, Only: HEGDataType

        Type (HEGDataType), Intent(In) :: HEGData

        Integer, Intent(In) :: M, N, L, EnergyCutoff, nBasis

        Integer :: K, eTarget, kTarget(3)
        logical :: InBasisSet

        kTarget=HEGData%nVec(M,:)+HEGData%nVec(N,:)-HEGData%nVec(L,:)

        ! Is this in the basis set at all?
        InBasisSet=.true.
        ! First, check range
        if (abs(kTarget(1)).gt.HEGData%cell) InBasisSet=.false. 
        if (abs(kTarget(2)).gt.HEGData%cell) InBasisSet=.false. 
        if (abs(kTarget(3)).gt.HEGData%cell) InBasisSet=.false. 
        ! Second, check energy
        eTarget=dot_product(kTarget,kTarget)
        if (eTarget.gt.EnergyCutoff) InBasisSet=.false.
        
        ! If it is out then return zero
        if (.not.InBasisSet) then
            K=0
            return
        endif

        ! Now it's safe to use the look-up table
        K=HEGData%kPointToBasisFn(kTarget(1),kTarget(2),kTarget(3))
        if (K.gt.nBasis) then
            K=0
            return
        endif

        !! corefcut breaks this because there are things within corefcut which
        !! now return zero because they are not in the basis
        !! I'm not sure there's a way to check this. 
        !
        !if ((FindIndex.eq.0) .or. (FindIndex.gt.nBasis))  then
        !    write (6,*) kTarget,FindIndex
        !    stop "Error in FindIndex" ! This is now an error since it should 
        !    ! have been filtered out before it hit the lookup table. We check 
        !    ! this because the table might have array bound errors
        !endif

    End Function FindIndex

    ! As FindIndex(m,n,l) but only returns non-zero index if it's also in the 
    ! occupied manifold
    !
    function FindIndexOcc(m,n,l)

        integer :: FindIndexOcc
        integer :: m,n,l

        !!!FindIndexOcc=FindIndex(m,n,l)
        !!!if (FindIndexOcc.gt.nEl/2) then 
        !!!    FindIndexOcc=0
        !!!endif

    end function

    ! As FindIndex(m,n,l) but only returns non-zero index if it's also in the 
    ! virtual manifold
    !
    function FindIndexVirt(m,n,l)

        integer :: FindIndexVirt
        integer :: m,n,l

        !!!FindIndexVirt=FindIndex(m,n,l)
        !!!if (.not.(FindIndexVirt.le.nEl/2)) then ! i.e. if it's not in the OCCUPIED manifold
        !!!    FindIndexVirt=0
        !!!endif

    end function

    ! This routine, principally a debug routine, checks for momentum symmetry of four indices
    ! 
    Pure Function CheckMomentumSymmetry(I,J,A,B,nVec) Result(SymAllowed)

        Integer, Intent(In) :: I, J, A, B, nVec(:,:)

        Integer :: G1(3), G2(3)
        Logical :: SymAllowed

        SymAllowed=.False.
        G1=-nVec(i,:)+nVec(a,:) ! i.e. g1 goes from i to a
        G2=-nVec(j,:)+nVec(b,:) ! i.e. g2 goes from j to b
        if (G1(1).eq.-G2(1).and.G1(2).eq.-G2(2).and.G1(3).eq.-G2(3)) SymAllowed=.True.

    End Function CheckMomentumSymmetry

    ! This takes an INTEGER VECTOR input and returns the electron repulsion integral
    !
    function ERI1(qvec)

        real(pr) :: ERI1
        integer :: qvec(3)

        !!!ERI1=(1/om)*(4.0_pr*pi/dot_product(qvec,qvec))

    end function

    ! This takes FOUR INTEGER ORBITAL INDICES and returns the electron repulsion integral
    ! for an excitation ij->ab
    ! Strangely, ERI is returned simply as a negative number if it's momentum disallowed
    ! (when the integral is actually zero)
    !
!    function ERI(a,b,i,j)
!
!        integer :: i,j,a,b
!        real(pr) :: ERI
!        integer :: qvec(3)
!
!        if (min(a,b,i,j) <= 0) then
!            stop "ERROR: this is not an allowed excitation"
!        endif
!        if (.not.check_momentum_symmetry(i,j,a,b)) then
!            ERI=-100.0_pr
!!           write(60,*) "WARNING: ERI called as zero, returned as -100!"
!            return
!        endif
!        qvec=-nvec(i,:)+nvec(a,:)
!        if ((qvec(1).eq.0).and.(qvec(2).eq.0).and.(qvec(3).eq.0)) then
!            ERI=-madelung
!            !stop "ERROR: zero momentum error"
!        else
!            !ERI=(1/om)*(4.0_pr*pi/dot_product(qvec,qvec))
!            ERI= One/(L*Pi*Dot_Product(qVec,qVec))
!        endif
!
!    end function

    
    subroutine change_rs(HEGData,UEGInfo)

        Use Types, Only: HEGDataType, UEGInfoType

        Type (HEGDataType), Intent(InOut) :: HEGData
        Type (UEGInfoType), Intent(In) :: UEGInfo

        real(pr) :: rescale, qvec(3)
        integer :: l1, l2

        rescale=UEGInfo%rS/HEGData%rS
        ! imagine old_rs as 1 and everything makes sense...
        HEGData%L=HEGData%L*rescale
        HEGData%Omega=HEGData%Omega*rescale**3.0_pr
        HEGData%kvec=HEGData%kvec*(1.0_pr/rescale)
        HEGData%madelung=CalcMadelung(HEGData%Omega,HEGData%L)*UEGInfo%MadelungFactor
        HEGData%Density=HEGData%Nel/HEGData%L**3.0_pr
        HEGData%FermiWaveVector=(3.0_pr*pi**2.0_pr*HEGData%Density)**(1.0_pr/3.0_pr)
        HEGData%ScreeningDistance=(4.0_dp/pi*HEGData%FermiWaveVector)**0.5_dp

        ! Form eigenvector array
        HEGData%Eigen=0.0_pr
        do l1=1,UEGInfo%NAO
            HEGData%eigen(l1)=0.0_pr
            HEGData%eigen(l1)=0.5_pr*dot_product(HEGData%kvec(l1,:),HEGData%kvec(l1,:)) ! kinetic
            do l2=1,HEGData%Nel/2
                if (l2.ne.l1) then
                    ! find "q" momentum transfer
                    qvec=HEGData%kvec(l1,:)-HEGData%kvec(l2,:)
                    HEGData%eigen(l1)=HEGData%eigen(l1)-(1/HEGData%Omega)*(4.0_pr*pi/dot_product(qvec,qvec))*UEGInfo%ExchangeFactor
                else 
                    HEGData%eigen(l1)=HEGData%eigen(l1)+HEGData%Madelung-UEGInfo%Gap ! self-interaction, madelung is a negative constant
                endif
            enddo
            write(60,*) l1, HEGData%nvec(l1,1), HEGData%nvec(l1,2), HEGData%nvec(l1,3), HEGData%eigen(l1)
        enddo
        
        call print_HF(HEGData)

    end subroutine
    
    subroutine GenerateTwist(HEGData)

        Use Types, only: HEGDataType

        Type (HEGDataType), Intent(InOut) :: HEGData
    
        write (60,*) "randomizing twist angle"
       
        HEGData%this_ntwist(1)=rand()-0.5_pr
        HEGData%this_ntwist(2)=rand()-0.5_pr
        HEGData%this_ntwist(3)=rand()-0.5_pr

        write (60,*) "generated",HEGData%this_ntwist

    end subroutine
    
    subroutine ChangeTwist(HEGData,UEGInfo)

        Use Types, only: HEGDataType, UEGInfoType

        Type (HEGDataType), Intent(InOut) :: HEGData
        Type (UEGInfoType), Intent(InOut) :: UEGInfo

        ! Flow through the program is:
        !   ReadInput looks at the input file
        !   Calls initHEG to set up constants etc.
        !   Setup basis is used to construct G1 and the eigenvector array etc.
        ! 
        ! So what we need to do is:
        !   Change global twist angle
        write (60,*) "changing twist angle from",HEGData% ntwist,"to",HEGData%this_ntwist
        HEGData%ntwist=HEGData%this_ntwist
        
        !   Deallocate arrays that were previously set up
        deallocate(HEGData%kVec)
        deallocate(HEGData%Eigen)
        deallocate(HEGData%nVec)
        deallocate(HEGData%kPointToBasisFn)
        
        !   Call setup basis to allocate new arrays (includes print_HF)
        call SetupBasis(HEGData,UEGInfo) 
        ! print_HF now includes print_ex_SF

    end subroutine


    subroutine print_HF(HEGData)

        Use Types, only: HEGDataType

        Type (HEGDataType), Intent(InOut) :: HEGData

        Integer :: l1, l2
        Real (Kind=pr) :: qvec(3)

        ! Print out E_HF
        HEGData%EHF=0.0_pr
        HEGData%XHF=0.0_pr
        do l1=1,HEGData%Nel/2
            HEGData%EHF=HEGData%EHF+0.5_pr*dot_product(HEGData%kvec(l1,:),HEGData%kvec(l1,:))
            do l2=1,HEGData%Nel/2
                if (l2.ne.l1) then
                    qvec=HEGData%kvec(l1,:)-HEGData%kvec(l2,:)
                    HEGData%EHF=HEGData%EHF-0.5_pr*(1/HEGData%Omega)*(4.0_pr*pi/dot_product(qvec,qvec))
                    HEGData%XHF=HEGData%XHF-0.5_pr*(1/HEGData%Omega)*(4.0_pr*pi/dot_product(qvec,qvec))
                endif
            enddo
        enddo
        HEGData%EHF=HEGData%EHF*2.0_pr
        HEGData%XHF=HEGData%XHF*2.0_pr
        write(60,*) "L:", HEGData%L
        write(60,*) "Omega:", HEGData%Omega
        write(60,*) "cell:", HEGData%Cell
        write(60,*) "HF energy:", HEGData%EHF/HEGData%Nel+HEGData%Madelung/2.0
        write(60,*) "Exchange energy:", HEGData%XHF/HEGData%Nel+HEGData%Madelung/2.0
        write(60,*) "Madelung:", HEGData%Madelung
        HEGData%EHF = HEGData%EHF + HEGData%Madelung*F12*HEGData%Nel

        call print_ex_SF(HEGData)

    end subroutine print_HF
    
    subroutine print_ex_SF(HEGData)

        Use Types, Only: HEGDataType

        Type (HEGDataType), Intent(InOut) :: HEGData

        integer :: l1, l2, I, J, A
        integer :: qvec(3)
        integer :: qvec_max(3), maxq
        integer, allocatable :: ExSf(:,:,:)
        real(pr), allocatable :: ExV(:,:,:)

        qvec_max=0
        
        ! Print out E_HF
        do l1=1,HEGData%Nel/2
            do l2=1,HEGData%Nel/2
                if (l2.ne.l1) then
                    qvec=abs(HEGData%nVec(l1,:)-HEGData%nVec(l2,:))
                    if (qvec(1).gt.qvec_max(1)) qvec_max(1)=qvec(1)
                    if (qvec(2).gt.qvec_max(2)) qvec_max(2)=qvec(2)
                    if (qvec(3).gt.qvec_max(3)) qvec_max(3)=qvec(3)
                endif
            enddo
        enddo
        maxq=maxval(qvec_max)+1
        allocate(ExSf(-maxq:maxq,-maxq:maxq,-maxq:maxq))
        allocate(ExV(-maxq:maxq,-maxq:maxq,-maxq:maxq))
        ExSf=0
        ExV=0.0_pr
        do l1=1,HEGData%Nel/2
            do l2=1,HEGData%Nel/2
                if (l2.ne.l1) then
                    qvec=HEGData%nVec(l1,:)-HEGData%nVec(l2,:)
                    ExSf(qvec(1),qvec(2),qvec(3))=ExSf(qvec(1),qvec(2),qvec(3))+1
                    ExV(qvec(1),qvec(2),qvec(3))=(1/HEGData%Omega)*&
                        (4.0_pr*pi/dot_product(HEGData%kvec(l1,:)-HEGData%kvec(l2,:), &
                                               HEGData%kvec(l1,:)-HEGData%kvec(l2,:)))
                endif
            enddo
        enddo
      
        write(82,*) "--------"
        write(82,*) "Starting a new structure factor"
        write(82,*) "--------"
        do I = -maxq,maxq
        do J = -maxq,maxq
        do A = -maxq,maxq
            if (ExSf(I,J,A) > 1.0E-10_pr) then
                write(82,*) I, J, A, ExSf(I,J,A),ExV(I,J,A)
            endif
        end do
        end do
        end do

        deallocate(ExSf,ExV)

    end subroutine


    ! This takes FOUR INTEGER ORBITAL INDICES and returns the electron repulsion integral
    ! for an excitation ij->ab
    ! Strangely, ERI is returned simply as a negative number if it's momentum disallowed
    ! (when the integral is actually zero)
    ! NEW: now has a flag argument for calculation of other kernels
    !
    !real(pr) function ERI(a,b,i,j,dummy_flag) result(intgrl)

    ! TODO - WZV
    ! Write two versions, called with a pointer one with and without
    ! the flags.

    Pure Function CasedERI(HEGData,UEGInfo,A,B,I,J,DummyFlag) result(Intgrl)

        ! Take in four orbital indices and calculate the
        ! two-particle integral for the UEG
        !
        ! In:
        !    i, j, a, b : int
        !        Orbital indices for the two-particle integral.
        !    dummy_flag : int
        !        Not used in the plain integral routine!
        ! Out:
        !    intgrl : float
        !        The two-particle integrals.

        Use Types, Only: HEGDataType, UEGInfoType

        Type (HEGDataType), Intent(In) :: HEGData
        Type (UEGInfoType), Intent(In) :: UEGInfo

        Integer, Intent(In) :: I, J, A, B
        Integer, Intent(In), Optional :: DummyFlag

        Integer :: Flag, qvec(3)
        Real (Kind=pr) :: Intgrl, qvec2, qvec2_scaled

        if (.not.present(DummyFlag)) then
            flag=0
        else
            flag=DummyFlag
        endif

        if (min(a,b,i,j) <= 0) then
            intgrl=-100.0_pr
            return
        endif
        If (.Not. CheckMomentumSymmetry(I,J,A,B,HEGData%nVec)) then
            intgrl=-100.0_pr
            return
        endif
        qvec=-HEGData%nvec(i,:)+HEGData%nvec(a,:)
        if ((qvec(1).eq.0).and.(qvec(2).eq.0).and.(qvec(3).eq.0)) then
            select case (flag)
                case (0) ! normal
                    intgrl=-HEGData%madelung
                case (1) ! short-range
                    intgrl=(1.0_dp/HEGData%L**3.0_dp)*4.0_pr*pi/HEGData%ScreeningDistance**2.0_pr
                case (2) ! long-range = normal-short range
                    intgrl=-HEGData%madelung-(1.0_dp/HEGData%L**3.0_dp)*4.0_pr*pi/HEGData%ScreeningDistance**2.0_pr
                case (3) ! completely screened, zero
                    intgrl=0.0_pr
                case default
                    intgrl=-100.0_pr
                    return
            end select
        else
            qvec2=dot_product(qvec,qvec)
            qvec2_scaled=qvec2*(2.0_pr*pi/HEGData%L)**2.0_pr
            select case (flag)
                case (0) ! normal
                    intgrl= One/(HEGData%L*Pi*qvec2)
                case (1) ! short-range
                    intgrl= (1.0_dp/HEGData%L**3.0_dp)*4.0_dp*pi/(HEGData%ScreeningDistance**2.0_dp+qvec2_scaled)
                case (2) ! long-range = normal-short range
                    intgrl= (1.0_dp/HEGData%L**3.0_dp)*4.0_dp*pi/(HEGData%ScreeningDistance**2.0_dp+qvec2_scaled) &
                        *HEGData%ScreeningDistance**2.0_dp/qvec2_scaled
                case (3) ! completely screened, zero
                    intgrl=0.0_pr
                case default
            end select
        endif

        if ((a.le.UEGInfo%CoreN).or.(b.le.UEGInfo%CoreN)&
            .or.(i.le.UEGInfo%CoreN).or.(j.le.UEGInfo%CoreN)) then ! CAS
            Intgrl=0.0_pr
        endif

        if (Intgrl.lt.0.0_pr) then
            intgrl=-100.0_pr
            return
        endif

        if ((Intgrl > UEGInfo%gapII).and.(UEGInfo%gapII > 0.0_pr)) then ! insulator II
            Intgrl=0.0_pr
        endif

    End Function CasedERI
    
    Pure Function NoCaseERI(HEGData,UEGInfo,A,B,I,J,DummyFlag) result(Intgrl)

        ! Take in four orbital indices and calculate the
        ! two-particle integral for the UEG
        !
        ! In:
        !    i, j, a, b : int
        !        Orbital indices for the two-particle integral.
        !    dummy_flag : int
        !        Not used in the plain integral routine!
        ! Out:
        !    intgrl : float
        !        The two-particle integrals.

        Use Types, Only: HEGDataType, UEGInfoType

        Type (HEGDataType), Intent(In) :: HEGData
        Type (UEGInfoType), Intent(In) :: UEGInfo

        Integer, Intent(In) :: I, J, A, B
        Integer, Intent(In), Optional :: DummyFlag

        Integer :: qvec(3)
        Real (Kind=pr) :: Intgrl, qvec2, qvec2_scaled

        if (min(a,b,i,j) <= 0) then
            Intgrl=-100.0_pr
            return
        endif
        If (.Not. CheckMomentumSymmetry(I,J,A,B,HEGData%nVec)) then
            Intgrl=-100.0_pr
            return
        endif

        qvec=-HEGData%nvec(i,:)+HEGData%nvec(a,:)
        if ((qvec(1).eq.0).and.(qvec(2).eq.0).and.(qvec(3).eq.0)) then
            Intgrl=-HEGData%Madelung
        else
            qvec2=dot_product(qvec,qvec)
            qvec2_scaled=qvec2*(2.0_pr*pi/HEGData%L)**2.0_pr
            Intgrl= One/(HEGData%L*Pi*qvec2)
        endif

        if ((a.le.UEGInfo%CoreN).or.(b.le.UEGInfo%CoreN)&
            .or.(i.le.UEGInfo%CoreN).or.(j.le.UEGInfo%CoreN)) then ! CAS
            Intgrl=0.0_pr
        endif

        if (Intgrl.lt.0.0_pr) then
            Intgrl=-100.0_pr
            return
        endif

        if ((Intgrl > UEGInfo%gapII).and.(UEGInfo%gapII > 0.0_pr)) then ! insulator II
            Intgrl=0.0_pr
        endif

    End Function NoCaseERI

    Pure Function CoreERI(HEGData,UEGInfo,A,B,I,J,DummyFlag) result(Intgrl)

        ! Take in four orbital indices and calculate the
        ! two-particle integral for the UEG
        !
        ! In:
        !    i, j, a, b : int
        !        Orbital indices for the two-particle integral.
        !    dummy_flag : int
        !        Not used in the plain integral routine!
        ! Out:
        !    intgrl : float
        !        The two-particle integrals.

        Use Types, Only: HEGDataType, UEGInfoType

        Type (HEGDataType), Intent(In) :: HEGData
        Type (UEGInfoType), Intent(In) :: UEGInfo

        Integer, Intent(In) :: I, J, A, B
        Integer, Intent(In), Optional :: DummyFlag

        Integer :: qvec(3)
        Real (Kind=pr) :: Intgrl, qvec2, qvec2_scaled

        if (min(a,b,i,j) <= 0) then
            Intgrl=-100.0_pr
            return
        endif
        If (.Not. CheckMomentumSymmetry(I,J,A,B,HEGData%nVec)) then
            Intgrl=-100.0_pr
            return
        endif

        qvec=-HEGData%nvec(i,:)+HEGData%nvec(a,:)
        if ((qvec(1).eq.0).and.(qvec(2).eq.0).and.(qvec(3).eq.0)) then
            Intgrl=-HEGData%Madelung
        else
            qvec2=dot_product(qvec,qvec)
            qvec2_scaled=qvec2*(2.0_pr*pi/HEGData%L)**2.0_pr
            Intgrl= One/(HEGData%L*Pi*qvec2)
        endif

        if (Intgrl.lt.0.0_pr) then
            Intgrl=-100.0_pr
            return
        endif

        if ((a.le.UEGInfo%CoreN).or.(b.le.UEGInfo%CoreN)&
            .or.(i.le.UEGInfo%CoreN).or.(j.le.UEGInfo%CoreN)) then ! CAS
            Intgrl=0.0_pr
        endif

    End Function CoreERI

    Pure Function GapERI(HEGData,UEGInfo,A,B,I,J,DummyFlag) result(Intgrl)

        ! Take in four orbital indices and calculate the
        ! two-particle integral for the UEG
        !
        ! In:
        !    i, j, a, b : int
        !        Orbital indices for the two-particle integral.
        !    dummy_flag : int
        !        Not used in the plain integral routine!
        ! Out:
        !    intgrl : float
        !        The two-particle integrals.

        Use Types, Only: HEGDataType, UEGInfoType

        Type (HEGDataType), Intent(In) :: HEGData
        Type (UEGInfoType), Intent(In) :: UEGInfo

        Integer, Intent(In) :: I, J, A, B
        Integer, Intent(In), Optional :: DummyFlag

        Integer :: qvec(3)
        Real (Kind=pr) :: Intgrl, qvec2, qvec2_scaled

        if (min(a,b,i,j) <= 0) then
            Intgrl=-100.0_pr
            return
        endif
        If (.Not. CheckMomentumSymmetry(I,J,A,B,HEGData%nVec)) then
            Intgrl=-100.0_pr
            return
        endif

        qvec=-HEGData%nvec(i,:)+HEGData%nvec(a,:)
        if ((qvec(1).eq.0).and.(qvec(2).eq.0).and.(qvec(3).eq.0)) then
            Intgrl=-HEGData%Madelung
        else
            qvec2=dot_product(qvec,qvec)
            qvec2_scaled=qvec2*(2.0_pr*pi/HEGData%L)**2.0_pr
            Intgrl= One/(HEGData%L*Pi*qvec2)
        endif

        if (Intgrl.lt.0.0_pr) then
            Intgrl=-100.0_pr
            return
        endif

        if ((Intgrl > UEGInfo%gapII).and.(UEGInfo%gapII > 0.0_pr)) then ! insulator II
            Intgrl=0.0_pr
        endif

    End Function GapERI

    Pure Function BaseERI(HEGData,UEGInfo,A,B,I,J,DummyFlag) Result(Intgrl)

        ! Take in four orbital indices and calculate the
        ! two-particle integral for the UEG
        !
        ! In:
        !    i, j, a, b : int
        !        Orbital indices for the two-particle integral.
        !    dummy_flag : int
        !        Not used in the plain integral routine!
        ! Out:
        !    intgrl : float
        !        The two-particle integrals.

        Use Types, Only: HEGDataType, UEGInfoType

        Type (HEGDataType), Intent(In) :: HEGData
        Type (UEGInfoType), Intent(In) :: UEGInfo

        Integer, Intent(In) :: I, J, A, B
        Integer, Intent(In), Optional :: DummyFlag

        Integer :: qvec(3)
        Real (Kind=pr) :: Intgrl, qvec2, qvec2_scaled

        qvec=-HEGData%nvec(i,:)+HEGData%nvec(a,:)
        if ((qvec(1).eq.0).and.(qvec(2).eq.0).and.(qvec(3).eq.0)) then
            Intgrl=-HEGData%Madelung
        else
            qvec2=dot_product(qvec,qvec)
            qvec2_scaled=qvec2*(2.0_pr*pi/HEGData%L)**2.0_pr
            Intgrl= One/(HEGData%L*Pi*qvec2)
        endif

    End Function BaseERI

    ! Modified form of ERI function above
    ! ===================================
    ! This takes FOUR INTEGER ORBITAL INDICES and returns the squared momentum value
    ! for an excitation ij->ab
    !
    Function dK2(HEGData,CoreN,A,B,I,J,InCore) Result(this_qvec2)

        Use Types, Only: HEGDataType

        Type (HEGDataType), Intent(In) :: HEGData

        Integer, Intent(In) :: CoreN, A, B, I, J

        Logical, Intent(InOut) :: InCore

        Integer :: this_qvec(3), this_qvec2

        this_qvec=-HEGData%nvec(i,:)+HEGData%nvec(a,:)
        this_qvec2=dot_product(this_qvec,this_qvec)
        inCore=.false.
        if ((a.le.CoreN).or.(b.le.CoreN).or.(i.le.CoreN).or.(j.le.CoreN)) then ! CAS
            inCore=.true.
        endif


    End Function dK2
    
    ! Modified form of dK2scaled function above
    ! =========================================
    ! This takes FOUR INTEGER ORBITAL INDICES and returns the squared momentum value
    ! for an excitation ij->ab
    !
    Function dK2scaled(HEGData,CoreN,A,B,I,J,InCore) Result(this_qvec2)

        Use Types, Only: HEGDataType

        Type (HEGDataType), Intent(In) :: HEGData

        Integer, Intent(In) :: CoreN, A, B, I, J

        Logical, Intent(InOut) :: InCore

        Integer :: this_qvec(3), this_qvec2

        this_qvec=-HEGData%nvec(i,:)+HEGData%nvec(a,:)
        this_qvec2=dot_product(this_qvec,this_qvec)
        this_qvec2=this_qvec2*(2.0_pr*pi/HEGData%L)**2.0_pr
        inCore=.false.
        if ((a.le.CoreN).or.(b.le.CoreN).or.(i.le.CoreN).or.(j.le.CoreN)) then ! CAS
            inCore=.true.
        endif


    end function
    
    Function dK(HEGData,CoreN,A,B,I,J,InCore) Result(this_qvec)

        Use Types, Only: HEGDataType

        Type (HEGDataType), Intent(In) :: HEGData

        Integer, Intent(In) :: CoreN, A, B, I, J

        Logical, Intent(InOut) :: InCore

        Integer :: this_qvec(3)

        this_qvec=-HEGData%nvec(i,:)+HEGData%nvec(a,:)
        inCore=.false.
        if ((a.le.CoreN).or.(b.le.CoreN).or.(i.le.CoreN).or.(j.le.CoreN)) then ! CAS
            inCore=.true.
        endif

    end function
   
    Subroutine SendForAllocation(HEGData,this_Sf,SfLength)

        Use Types, Only: HEGDataType

        Type (HEGDataType), Intent(In) :: HEGData

        Real (Kind=pr), Intent(InOut), Allocatable :: this_Sf(:,:,:)
        Integer, Intent(InOut) :: SfLength

        SfLength=maxval(HEGData%nvec(:,:))*2+1
        allocate(this_Sf(-SfLength:SfLength,-SfLength:SfLength,-SfLength:SfLength))
        !!!write(80,*) "SfLength=", SfLength
        this_Sf=0.0_pr

    End Subroutine SendForAllocation


    Subroutine FindTol(rs,TolMax,FailRatio,DenomFactor)

        Real (Kind=pr), Intent(In) :: rs
        Real (Kind=pr), Intent(InOut) :: TolMax, FailRatio, DenomFactor

        TolMax = 1.0E-8_pr
        FailRatio = 100
        if (rs.gt.5.1_pr) TolMax = 1.0E-6_pr
        if (rs.gt.5.1_pr) denomfactor=5.0_pr
        if (rs.gt.9.9_pr) TolMax = 1.0E-5_pr
        if (rs.gt.80.0_pr) TolMax = 1.0E-3_pr
        if (rs.lt.0.95_pr) FailRatio = 1E6

    End Subroutine FindTol
    
end module 
