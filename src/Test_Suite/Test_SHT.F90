!
!  Copyright (C) 2018 by the authors of the RAYLEIGH code.
!
!  This file is part of RAYLEIGH.
!
!  RAYLEIGH is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  RAYLEIGH is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!


Module Test_SHT
    !/////////////////////////////////////////////////
    !  Spherical Harmonic Transform Testing Module
    Use ProblemSize
    Use Fields
    Use Parallel_Framework
    Use Fourier_Transform
    Use Legendre_Transforms, Only : Legendre_Transform
    Use SendReceive
    Implicit None
    Integer :: ntest_legendre =1
Contains

    Subroutine Test_Spherical_Transforms()
        !Call Amp_Test_Parallel()
        Call test_LT()
    End Subroutine Test_Spherical_Transforms

    Subroutine Test_LT()
        Integer :: colrank, rowrank, nf=3, lval, mval
        Integer :: mp, l, m, r, i, f
        Integer :: fcount(3,2)
        Real*8 :: mxdiff, diff, diff1, diff2, ans, norm
        Real*8, allocatable :: to_phys_norm(:), to_spec_norm(:)
        type(SphericalBuffer) :: mytest

        fcount(:,:) = nf

        ! allocate fields in spectral space
        Call mytest%init(field_count = fcount, config = 's2a')
        Call mytest%construct('s2a')

        ! initialize specific (l,m) values
        lval = 2; mval = 1
        do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            do l=m,l_max
                mytest%s2a(mp)%data(l,:,:,:) = 0.0d0
                if ((l .eq. lval) .and. (m .eq. mval)) then
                    do r=my_r%min,my_r%max
                        do f=1,nf
                            mytest%s2a(mp)%data(l,r,1,f) = (2.0d0*f+1.0d0)*radius(r)
                            mytest%s2a(mp)%data(l,r,2,f) = -(3.0d0*f*f+2.0d0)*radius(r)**2
                        enddo
                    enddo
                endif
            enddo
        enddo

        ! build FFT normalizations so they can be accounted for later
        !     if using Rayleigh: defined in Legendre_Polynomials.F90 (search PTS or STP)
        !     if using SHTns: defined in Legendre_Transforms.F90 (search PTS or STP)
        allocate(to_phys_norm(0:l_max), to_spec_norm(0:l_max))
        to_phys_norm(0) = 1.0d0  ! m=0
        to_phys_norm(1:) = 0.5d0 ! m/=0
        to_spec_norm(0) = 1.0d0/(n_phi)    ! m=0
        to_spec_norm(1:) = 1.0d0/(n_theta) ! m/=0

        ! build physical space
        Call mytest%construct('p2a')

        ! run tests
        do i=1, ntest_legendre
            Call Legendre_Transform(mytest%s2a,mytest%p2a) ! to physical

            Call Legendre_Transform(mytest%p2a,mytest%s2a) ! back to spectral

            ! zero out l_max modes (this is done in production runs too)
            do mp=my_mp%min,my_mp%max
                mytest%s2a(mp)%data(l_max,:,:,:) = 0.0d0
            enddo
        enddo

        ! compute error
        diff = -1.0d0
        mxdiff = -1.0d0
        do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            norm = (to_phys_norm(m)*to_spec_norm(m))**ntest_legendre ! undo FFT normalization
            do l=m,l_max
                if ((l .eq. lval) .and. (m .eq. mval)) then
                    do r=my_r%min,my_r%max
                        do f=1,nf
                           ans = (2.0d0*f + 1.0d0)*radius(r)*norm
                           diff1 = abs(ans - mytest%s2a(mp)%data(l,r,1,f))

                           ans = -(3.0d0*f*f + 2.0d0)*radius(r)**2*norm
                           diff2 = abs(ans - mytest%s2a(mp)%data(l,r,2,f))

                           diff = max(diff1, diff2)
                        enddo
                    enddo
                else
                    diff = maxval(abs(mytest%s2a(mp)%data(l,:,:,:)))
                endif
                if (diff .gt. mxdiff) mxdiff = diff
            enddo
        enddo
        rowrank = pfi%rcomm%rank
        colrank = pfi%ccomm%rank
        write(*,*) 'rowrank=',rowrank,'colrank=',colrank,'max error=',mxdiff

        ! cleanup
        Call mytest%deconstruct('p2a')
        Call mytest%deconstruct('s2a')

        deallocate(to_phys_norm, to_spec_norm)

    End Subroutine Test_LT

    Subroutine Amp_Test()
        Implicit None
        Integer :: i,m,mp,l, mxl,mxm
        Integer :: fcount(3,2)

        Integer :: colrank, rowrank, nmodes, this_mode, mcount, offset
        Integer :: p, np, testing_tag = 90
        Real*8  :: mxdiff
        Real*8, Allocatable :: buff(:),reldiff(:), reldiffs(:)
        Integer, Allocatable :: all_l_values(:), all_m_values(:)
        Logical :: report
        Character*4 :: lstring
        type(SphericalBuffer) :: test
        fcount(:,:) = 1


        !//////////////////////////////////////////////////////////
        ! Test 1:   Legendre Transform White Noise Test.
        !                 All modes given unit power.
        Call test%init(field_count = fcount, config = 's2a')
        Call test%construct('s2a')

        nmodes = 0
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            Do l = m, l_max
                nmodes = nmodes+1
                test%s2a(mp)%data(l,:,:,:) = 1.0d0
            Enddo
        Enddo

        Allocate(reldiff(1:nmodes))
        Call test%construct('p2a')

        Do i = 1, ntest_legendre

            Call Legendre_Transform(test%s2a,test%p2a)
            !Do mp = my_mp%min, my_mp%max
                    !test%s2a(mp)%data(:,:) = 0.0d0
            !Enddo
            Call Legendre_Transform(test%p2a,test%s2a)
        Enddo

        mxdiff = -1.0d0
        mxl = -1
        mxm = -1
        this_mode = 1
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            Do l = m, l_max
                reldiff(this_mode) = 1.0d0-test%s2a(mp)%data(l,1,1,1)
                this_mode = this_mode+1
                !If (diff .gt. mxdiff) then
                !    mxdiff = diff
                !    mxl = l
                !    mxm = m
                !Endif

            Enddo
        Enddo

        !Write(6,*)'My maxdiff was at: ', mxl,mxm,mxdiff
        Call test%deconstruct('p2a')
        Call test%deconstruct('s2a')
        colrank = pfi%ccomm%rank
        report = .false.

        if (colrank .eq. 0) report = .true.

        If (report) Then
            rowrank = pfi%rcomm%rank
            np = pfi%rcomm%np
            If (rowrank .eq. 0) Then
                !/////////////
                !  First pass - book keeping
                mcount = 0
                Do p = 0, np -1
                    Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
                        m = m_values(mp)
                        Do l = m, l_max
                            mcount = mcount+1
                        Enddo
                    Enddo
                enddo
                offset = 1
                Allocate(all_l_values(1:mcount))
                Allocate(all_m_values(1:mcount))
                Allocate(reldiffs(1:mcount))
                reldiffs(1:nmodes) = reldiff(1:nmodes)
                Do mp = my_mp%min, my_mp%max
                    m = m_values(mp)
                    Do l = m, l_max
                        all_l_values(offset) = l
                        all_m_values(offset) = m
                        offset = offset+1
                    enddo
                enddo
                !///////////////////
                ! now receive from everyone
                !offset1 = nmodes+1
                Do p = 1, np -1
                    mcount = 0
                    Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
                        m = m_values(mp)
                        Do l = m, l_max
                            all_l_values(offset) = l
                            all_m_values(offset) = m
                            offset = offset+1
                            mcount = mcount+1
                        Enddo
                    Enddo
                    Allocate(buff(1:mcount))
                    buff(:) = -20.0d0
                    Call receive(buff, source= p,tag=testing_tag,grp = pfi%rcomm)
                    offset = offset-mcount
                    reldiffs(offset:offset+mcount-1) = buff(1:mcount)
                    DeAllocate(buff)
                    offset = offset+mcount
                Enddo
                Write(lstring,'(i4.4)') l_max
            Open(unit=15,file='legendre_accuracy_'//lstring,form='unformatted', status='replace')
            Write(15)ntest_legendre
                Write(15)(offset-1)
            Write(15)(reldiffs(i),i=1,offset-1)
            Write(15)(all_l_values(i),i=1,offset-1)
                Write(15)(all_m_values(i),i=1,offset-1)
            Close(15)

                DeAllocate(reldiffs)
                DeAllocate(all_l_values)
                DeAllocate(all_m_values)
            Else
                Call send(reldiff, dest = 0,tag=testing_tag, grp=pfi%rcomm)

            Endif
        Endif
        DeAllocate(reldiff)
    End Subroutine Amp_Test

    Subroutine Amp_Test_Parallel()
        Implicit None
        Integer :: i,m,mp,l, mxl,mxm
        Integer :: fcount(3,2)
        Integer :: nrl
        Integer :: nf = 3
        Integer :: colrank, rowrank, nmodes, this_mode, mcount, offset
        Integer :: p, np, testing_tag = 90, testing_tag2 = 91
        Real*8  :: mxdiff, ans
        Real*8, Allocatable :: buff(:),reldiff_real(:), reldiff_imag(:), reldiffs_real(:), reldiffs_imag(:)
        Integer, Allocatable :: all_l_values(:), all_m_values(:)
        Logical :: report
        Character*4 :: lstring
        type(SphericalBuffer) :: test
        fcount(:,:) = nf


        !//////////////////////////////////////////////////////////
        ! Tests The Full Spherical Harmonic Transform
        ! (Good for checking transposes and normalization of FFT routine)
        ! We also test multiple fields in the transpose... Do this tomorrow
        Call test%init(field_count = fcount, config = 's2a')
        Call test%construct('s2a')
        nrl = my_r%delta
        nmodes = 0

        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            Do l = m, l_max
                nmodes = nmodes+1
                test%s2a(mp)%data(l,1:nrl,1,1) = (l+1)*1.0d0
                test%s2a(mp)%data(l,1:nrl,2,1) = (m+1)*1.0d0
                if ((m .eq. 0) ) then
                    ! Be careful not to give m = 0 any imaginary part
                    ! That would be absolutely foolish
                    test%s2a(mp)%data(l,1:nrl,1,1) = 14.1d0
                    test%s2a(mp)%data(l,1:nrl,2,1) = 0.0d0
                Endif
            Enddo
        Enddo

        Allocate(reldiff_real(1:nmodes))
        Allocate(reldiff_imag(1:nmodes))
        Call test%construct('p2a')

        Do i = 1, ntest_legendre

            Call Legendre_Transform(test%s2a,test%p2a)
            test%config = 'p2a'
            Call test%reform()
            Call fft_to_physical(test%p3a,rsc = .true.)

            Call test%construct('p3b')
            test%p3b = test%p3a
            test%config = 'p3b'
            Call test%deconstruct('p3a')
            Call fft_to_spectral(test%p3b, rsc = .true.)
            Call test%reform()    ! Move to p2b

            !Call test%deconstruct('p3b')
            Call Legendre_Transform(test%p2b,test%s2a)
            Call test%deconstruct('p2b')
            Call test%construct('p2a')
        Enddo

        mxdiff = -1.0d0
        mxl = -1
        mxm = -1
        this_mode = 1
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            Do l = m, l_max
                if ( (m .eq. 0) ) then
                    reldiff_real(this_mode) = (14.1-test%s2a(mp)%data(l,1,1,1))/14.1
                    reldiff_imag(this_mode) = test%s2a(mp)%data(l,1,2,1)
                    write(6,*)test%s2a(mp)%data(l,1,1,1)
                else
                    ans = (l+1)*1.0d0
                    reldiff_real(this_mode) = (ans-test%s2a(mp)%data(l,1,1,1))/ans
                    ans = (m+1)*1.0d0
                    reldiff_imag(this_mode) = (ans-test%s2a(mp)%data(l,1,2,1))/ans
                endif


                this_mode = this_mode+1
                !If (diff .gt. mxdiff) then
                !    mxdiff = diff
                !    mxl = l
                !    mxm = m
                !Endif

            Enddo
        Enddo

        !Write(6,*)'My maxdiff was at: ', mxl,mxm,mxdiff
        Call test%deconstruct('p2a')
        Call test%deconstruct('s2a')
        colrank = pfi%ccomm%rank
        report = .true.

        if (colrank .eq. 0) report = .true.

        If (report) Then
            rowrank = pfi%rcomm%rank
            np = pfi%rcomm%np
            If (rowrank .eq. 0) Then
                !/////////////
                !  First pass - book keeping
                mcount = 0
                Do p = 0, np -1
                    Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
                        m = m_values(mp)
                        Do l = m, l_max
                            mcount = mcount+1
                        Enddo
                    Enddo
                enddo
                offset = 1
                Allocate(all_l_values(1:mcount))
                Allocate(all_m_values(1:mcount))
                Allocate(reldiffs_real(1:mcount))
                Allocate(reldiffs_imag(1:mcount))

                reldiffs_real(1:nmodes) = reldiff_real(1:nmodes)
                reldiffs_imag(1:nmodes) = reldiff_imag(1:nmodes)
                Do mp = my_mp%min, my_mp%max
                    m = m_values(mp)
                    Do l = m, l_max
                        all_l_values(offset) = l
                        all_m_values(offset) = m
                        offset = offset+1
                    enddo
                enddo
                !///////////////////
                ! now receive from everyone
                !offset1 = nmodes+1
                Do p = 1, np -1
                    mcount = 0
                    Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
                        m = m_values(mp)
                        Do l = m, l_max
                            all_l_values(offset) = l
                            all_m_values(offset) = m
                            offset = offset+1
                            mcount = mcount+1
                        Enddo
                    Enddo
                    Allocate(buff(1:mcount))
                    buff(:) = -20.0d0
                    Call receive(buff, source= p,tag=testing_tag,grp = pfi%rcomm)
                    offset = offset-mcount
                    reldiffs_real(offset:offset+mcount-1) = buff(1:mcount)
                    Call receive(buff, source= p,tag=testing_tag2,grp = pfi%rcomm)
                    reldiffs_imag(offset:offset+mcount-1) = buff(1:mcount)
                    DeAllocate(buff)
                    offset = offset+mcount
                Enddo
                Write(lstring,'(i4.4)') l_max
            Open(unit=15,file='SHT_accuracy_'//lstring,form='unformatted', status='replace')
            Write(15)ntest_legendre
                Write(15)(offset-1)
            Write(15)(reldiffs_real(i),i=1,offset-1)
                Write(15)(reldiffs_imag(i),i=1,offset-1)
            Write(15)(all_l_values(i),i=1,offset-1)
                Write(15)(all_m_values(i),i=1,offset-1)
            Close(15)

                DeAllocate(reldiffs_real, reldiffs_imag)
                DeAllocate(all_l_values)
                DeAllocate(all_m_values)
            Else
                Call send(reldiff_real, dest = 0,tag=testing_tag, grp=pfi%rcomm)
                Call send(reldiff_imag, dest = 0,tag=testing_tag2, grp=pfi%rcomm)
            Endif
        Endif
        DeAllocate(reldiff_real)
        DeAllocate(reldiff_imag)
    End Subroutine Amp_Test_Parallel

End Module Test_SHT
