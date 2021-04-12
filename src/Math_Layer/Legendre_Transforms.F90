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

! This module contains routines for performing the Legendre transforms.
! These are really just matrix-matrix multiplies that we offload to dgemm.
! The polynomials themselves are computed in Legendre_Polynomials.F90 at startup.
! Parity is accounted for, as well as "Legendre" data structures.
! Notes:
!            (Dec 5, 2013 - 9:15 p.m.) : Verified that this module sees the correct lower bounds of
!                        data_in/out(i)%data, which is allocated data(m_value(i):l_max) in external modules.
!                        Was concerned that it might index array data starting at 1, not m_value(i) [Fortran standard?]

Module Legendre_Transforms
    Use Legendre_Polynomials
    Use Structures
#ifdef USE_SHTns
    Use iso_c_binding
    !Type, Public :: rmcontainer
    !    Real*8, Allocatable :: data(:,:)
    !End Type rmcontainer
    Implicit None

    ! include the SHTns interfaces
    Include 'shtns.f03'

    Type(c_ptr) :: SHTns_c             ! main SHTns structure built on initialization
    Type(SHTns_info), Pointer :: SHTns ! the Fortran friendly main SHTns structure

    Integer :: SHTns_nthreads          ! keep track of how many threads SHTns is using

    ! include FFT normalizations in the Legendre transform
    Real*8, Allocatable :: PTS_normalization(:), STP_normalization(:)

    Interface Legendre_Transform
        Module Procedure SHTns_ToSpectral, SHTns_ToPhysical
    End Interface
#else
    Interface Legendre_Transform
        Module Procedure PtS_2d_dgpv2, StP_2d_dgp
        Module Procedure PtS_4d_dgpv3, StP_4d_dgp2  ! <<<<<< These are the two legendre transforms that are used in practice (at bottom of file)
    End Interface
#endif
Contains

#ifndef USE_SHTns
Subroutine Test_Legendre
    Implicit None
    Real*8, Allocatable :: theta(:),tmp(:)
    Real*8, Allocatable :: ylm(:,:)
    Integer :: i,l,m, nrhs
    Integer, Allocatable :: l_test(:), m_test(:)
    Type(p_lm_array), Allocatable :: ans(:), ans2(:)
    Real*8 :: st, ct, chk
    ! Test the physical to spectral legendre transform
    !  using tabulated associated legendre transforms
    nrhs = 6
    Allocate(theta(1:n_theta))
    Allocate(tmp(1:n_theta))
    Allocate(ylm(1:n_theta,1:nrhs)) ! Test the first six spherical harmonics (modulu the exp(i m phi) piece)
    Allocate(ans(1:n_m))
    Allocate(ans2(1:n_m))
    Allocate(l_test(1:nrhs))
    Allocate(m_test(1:nrhs))

    Do m = 1,n_m
        Allocate( ans(m)%data(m_values(m):l_max,nrhs))
        Allocate(ans2(m)%data(m_values(m):l_max,nrhs))
         ans(m)%data(:,:) = 0.0d0
        ans2(m)%data(:,:) = 53.0d0
    Enddo

    Do i = 1, n_theta
        theta(i) = acos(coloc(i))
        ct = cos(theta(i))
        st = sin(theta(i))
        ylm(i,1) = (0.25d0/PiQuad)**0.5d0    !00
        ylm(i,2) = -st*(3.0d0/8.0d0/PiQuad)**0.5d0 !11
        ylm(i,3) = ct*(3.0d0/4.0/PiQuad)**0.5d0 ! 10
        ylm(i,4) = st*st*0.25d0*(15.0d0/2.0d0/PiQuad)**0.5d0 !22
        ylm(i,5) = -st*ct*(15.0d0/8.0d0/PiQuad)**0.5d0 !21
        ylm(i,6) = (1.5d0*ct*ct-0.5d0)*(5.0d0/4.0/PiQuad)**0.5d0    !20
        tmp(i) = 1.0d0/(1-coloc(i)*coloc(i))
    Enddo
    l_test(1) = 0
    l_test(2:3) = 1
    l_test(4:6) = 2
    m_test(1) = 0
    m_test(2) = 1
    m_test(3) = 0
    m_test(4) = 2
    m_test(5) = 1
    m_test(6) = 0

    !Next perform the transform on our test harmonics

    If (parity) Then
        write(6,*)'testing parity transform'
        Call Pts_2D_drp(ylm,ans,nrhs)
        Call PtS_2D_dgp(ylm, ans2,nrhs)
    Else
        Call Pts_2D_dr(ylm, ans,nrhs)
        Call PtS_2D_dg(ylm, ans2,nrhs)
    Endif


    !Verify the orthogonality of m1 = m2, but l1 ne l2
    Do i = 1, nrhs
        Write(6,*)'============================='
        Write(6,*)'     Test mode  '
        Write(6,*)'l = ', l_test(i), 'm = ', m_test(i)
        Do m = 1, n_m
            if (m_test(i) .eq. m_values(m)) then
                Do l = m_values(m), l_max
                    Write(6,*)l, m_values(m), ans(m)%data(l,i), ans2(m)%data(l,i)
                Enddo
            Endif
        Enddo
    Enddo

    ! Next verify the orthogonality of m1 ne m2, but l1 eq l2
    !Next perform the transform on our test harmonics
    Do i = 1, nrhs
        ylm(:,i) = ylm(:,i)*tmp(:)
    Enddo


    If (parity) Then
        Call Pts_2D_drp(ylm, ans,nrhs)
        Call PtS_2D_dgp(ylm, ans2,nrhs)
    Else
        Call Pts_2D_dr(ylm, ans,nrhs)
        Call PtS_2D_dg(ylm, ans2,nrhs)
    Endif


    Do i = 1, nrhs
        Write(6,*)'============================='
        Write(6,*)'     Test 2   '
        Write(6,*)'l = ', l_test(i), 'm = ', m_test(i)
        Do m = 1, n_m
            Do l = m_values(m), l_max
                if (l_test(i) .eq. l) then
                    Write(6,*)l, m_values(m), ans(m)%data(l,i), ans2(m)%data(l,i)
                    if(m_test(i) .eq. m_values(m)) then
                        chk = (2.0d0*l+1.0d0)/2.0d0 !4.0d0/pi
                        chk = chk/m_values(m)
                        write(6,*)'check: ', chk
                    endif
                Endif
            Enddo
        Enddo
    Enddo



103 format(d11.2)

    DeAllocate(theta)
    DeAllocate(ylm)
    DeAllocate(tmp)
    Do m = 1,n_m
        DeAllocate(ans(m)%data)
        DeAllocate(ans2(m)%data)
    Enddo
    DeAllocate(ans,ans2)
    DeAllocate(l_test,m_test)
End Subroutine Test_Legendre

Subroutine Test_Simple_Dgemm
    Real*8, Allocatable :: a(:,:), b(:,:), c(:,:)
    Integer :: m, n, k, i, ntimes
    Real*8 :: alpha,beta

    ! Just to make sure I've done dgemm correctly
    ! a = | 1 2 |  b = |5|  a#b = c = | 17 |
    !     = | 3 4 |      |6|              = | 39 |

    m = 2
    n = 2
    k = 2
    Allocate(a(1:2,1:2))
    Allocate(b(1:2,1:2))
    Allocate(c(1:2,1:2))

    a(1,1) = 1.0d0
    a(2,1) = 2.0d0
    a(1,2) = 3.0d0
    a(2,2) = 4.0d0
    b(1,1) = 5.0d0
    b(2,1) = 6.0d0
    b(1,2) = 7.0d0
    b(2,2) = 8.0d0
    c(:,:) = 0.0d0
    write(6,*)'calling dgemm'
    alpha = 1.0d0
    beta = 0.0d0
    ntimes = 100
    do i = 1, ntimes
        CALL DGEMM('T','N',m,n,k, alpha, a, m, b, k, beta,c,m)
    enddo
    Write(6,*)'c is : ', c

    a(1,1) = 1.0d0
    a(1,2) = 2.0d0
    a(2,1) = 3.0d0
    a(2,2) = 4.0d0
    CALL DGEMM('N','N',m,n,k, alpha, a, m, b, k, beta,c,m)
    Write(6,*)'c is : ', c

    b(1,1) = 5.0d0
    b(1,2) = 6.0d0
    b(2,1) = 7.0d0
    b(2,2) = 8.0d0
    CALL DGEMM('N','T',m,n,k, alpha, a, m, b, k, beta,c,m)
    Write(6,*)'c is : ', c

    a(1,1) = 1.0d0
    a(2,1) = 2.0d0
    a(1,2) = 3.0d0
    a(2,2) = 4.0d0
    CALL DGEMM('T','T',m,n,k, alpha, a, m, b, k, beta,c,m)
    Write(6,*)'c is : ', c

    DeAllocate(a,b,c)



End Subroutine Test_simple_dgemm


Subroutine Test_Simple_Dgemm2
    Real*8, Allocatable :: a(:,:), b(:,:), c(:,:)
    Integer :: m, n, k, i, ntimes
    Real*8 :: alpha, beta

    !//////////////////////////////////////////////
    ! Just to make sure I've done dgemm correctly
    !
    ! a = | 1 2 |  b = |7  9|  a#b = c = | 23  29 |
    !       | 3 4 |      |8 10|                | 53  67 |
    !        | 5 6 |                        | 83 105 |


    !//////////////////////////////////
    !  Test 1.  Columns run fastest.
    m = 3
    n = 2
    k = 2
    Allocate(a(1:m,1:k))

    Allocate(b(1:k,1:n))

    Allocate(c(1:m,1:n))



    a(1,1) = 1.0d0
    a(1,2) = 2.0d0
    a(2,1) = 3.0d0
    a(2,2) = 4.0d0
    a(3,1) = 5.0d0
    a(3,2) = 6.0d0




    b(1,1) = 7.0d0
    b(1,2) = 9.0d0
    b(2,1) = 8.0d0
    b(2,2) = 10.0d0

    c(:,:) = 0.0d0
    write(6,*)'calling dgemm'
    alpha = 1.0d0
    beta = 0.0d0
    ntimes = 1
    do i = 1, ntimes
        CALL DGEMM('N','N',m,n,k, alpha, a, m, b, k, beta,c,m)
    enddo
    Write(6,*)'c is : ', c



    DeAllocate(a,b,c)

    !//////////////////////////////////
    !  Test 2.  Rows of A run fastest.
    !  B and C remain the same
    m = 3
    n = 2
    k = 2
    Allocate(a(1:k,1:m))

    Allocate(b(1:k,1:n))

    Allocate(c(1:m,1:n))



    a(1,1) = 1.0d0
    a(2,1) = 2.0d0
    a(1,2) = 3.0d0
    a(2,2) = 4.0d0
    a(1,3) = 5.0d0
    a(2,3) = 6.0d0




    b(1,1) = 7.0d0
    b(1,2) = 9.0d0
    b(2,1) = 8.0d0
    b(2,2) = 10.0d0

    c(:,:) = 0.0d0
    write(6,*)'calling dgemm'
    alpha = 1.0d0
    beta = 0.0d0
    ntimes = 1
    do i = 1, ntimes
        CALL DGEMM('T','N',m,n,k, alpha, a, k, b, k, beta,c,m)
    enddo
    Write(6,*)'c is : ', c



    DeAllocate(a,b,c)



End Subroutine Test_simple_dgemm2
#endif




#ifdef USE_SHTns
   Subroutine SHTns_Initialize(n_threads, &
                               on_the_fly, information, polar_threshold, theta_contiguous)
       Integer, intent(in) :: information, n_threads
       Logical, intent(in) :: on_the_fly, theta_contiguous
       Real*8, intent(in) :: polar_threshold

       Integer :: norm, layout, m_res, m_max, n_phi

       n_phi = 2*n_theta
       m_max = l_max

       m_res = 1 ! 2*pi/m_res is the azimuthal periodicity, m_max*m_res is max azimuthal order

       ! choose grid, data layout, and Y_l^m normalization --- be consistent with Rayleigh
       If (on_the_fly) Then
           layout = SHT_gauss_fly
       Else
           layout = SHT_gauss
       Endif
       ! we only need the scalar transforms
       ! Rayleigh has x in (-1,1) & theta in (pi,0) ---> so south pole is first
       layout = layout + SHT_scalar_only + SHT_south_pole_first
       if (theta_contiguous) then
          layout = layout + SHT_theta_contiguous
       else
          layout = layout + SHT_phi_contiguous
       endif

       ! Rayleigh uses the very sane choice of orthonormal Y_l^m
       norm = SHT_orthonormal

       Call SHTns_verbose(information) ! set how much information SHTns will display

       SHTns_nthreads = SHTns_use_threads(n_threads) ! set OpenMP threads

       ! initialize/allocate transforms and build useful arrays
       SHTns_c = SHTns_create(l_max, m_max, m_res, norm)

       ! attach a grid to the SHT object & determine optimal algorithm
       Call SHTns_set_grid(SHTns_c, layout, polar_threshold, n_theta, n_phi)

       ! map the C SHTns structure to the Fortran one
       Call C_F_pointer(cptr=SHTns_C, fptr=SHTns)
       !Call C_F_pointer(cptr=SHTns%ct, fptr=costheta, shape=[SHTns%nlat])
       !Call C_F_pointer(cptr=SHTns%st, fptr=sintheta, shape=[SHTns%nlat])

       ! apply some m-dependent FFT normalizations during the Legendre Transforms
       Allocate(PTS_normalization(0:l_max), STP_normalization(0:l_max))

       PTS_normalization(0) = 1.0d0/n_phi    ! m=0
       STP_normalization(0) = 1.0d0
       PTS_normalization(1:) = 1.0d0/n_theta ! m/=0
       STP_normalization(1:) = 0.5d0

   End Subroutine SHTns_Initialize

   Subroutine SHTns_Finalize()
       Call SHTns_unset_grid(SHTns_c)
       Call SHTns_destroy(SHTns_c)
       DeAllocate(PTS_normalization, STP_normalization)
   End Subroutine SHTns_Finalize

   Subroutine SHTns_ToSpectral(data_in, data_out)
       Real*8, Intent(In) :: data_in(:,:,:)
       Type(rmcontainer4d), Intent(InOut) :: data_out(1:)
       ! ingoing data has shape:
       !   data_out(th,nrhs,lm)
       !       th = theta (in-processor)
       !     nrhs = number of RHS elements, stacked over: radius/real/imag/nfields
       !       lm = mode index (distributed)
       ! nrhs axis is 1-based indexing, except for radius, ordered: radius-real-imag-field
       !     index = (r-rlo+1) + (imi-1)*Nr + (f-1)*Nr*2
       !
       ! outgoing data has shape:
       !   data_out(lm)%data(l,r,imi,nf)
       !       lm = mode index (distributed)
       !        l = spherical harmonic (in-processor)
       !        r = radius (distributed)
       !      imi = real/imag parts
       !       nf = number of fields
       Complex*16, Allocatable :: temp_spec(:), temp_phys(:)
       Complex*16 :: ai, ar
       Integer :: m, f, r, mode
       Integer :: ddims(3), oddims(4)
       Integer :: n_m, nrhs, nfield, rmn, rmx, my_Nr, indr, indi, lind
       Real*8 :: norm

       ai = (0.0d0, 1.0d0) ! imaginary unit
       ar = (1.0d0, 0.0d0) ! regular one

       ! find lower/upper bounds of incoming/outgoing data
       oddims = shape(data_out(1)%data)
       nfield = oddims(4)
       rmn = LBOUND(data_out(1)%data,2)
       rmx = UBOUND(data_out(1)%data,2)

       my_Nr = rmx - rmn + 1

       ddims = shape(data_in)
       n_m = ddims(3)
       nrhs = ddims(2)

       Allocate(temp_spec(1:l_max+1), temp_phys(1:n_theta)) ! SHTns expects complex arrays

       Do mode = 1, n_m ! loop over lm modes

           m = m_values(mode) ! extract actual m value

           norm = PTS_normalization(m) ! extract FFT normalizations, based on m

           lind = l_max - m + 1 ! number of modes for this m-value

           Do f = 1, nfield ! number of fields
               Do r = rmn, rmx ! radius

                   ! package incoming data into complex array at this radius/field
                   indr = (r-rmn+1) + (f-1)*my_Nr*2
                   indi = (r-rmn+1) + my_Nr + (f-1)*my_Nr*2
                   temp_phys(:) = ar*data_in(:,indr,mode) + ai*data_in(:,indi,mode)

                   temp_spec(:) = 0.0d0 ! spectral output is of size (lmax-m+1)
                   Call spat_to_sh_ml(SHTns_c, m, temp_phys, temp_spec(1:lind), l_max)

                   ! extract results and store in output array
                   data_out(mode)%data(m:l_max,r,1,f) = norm*real(temp_spec(1:lind))
                   data_out(mode)%data(m:l_max,r,2,f) = norm*aimag(temp_spec(1:lind))

               Enddo
           Enddo
       Enddo

       DeAllocate(temp_spec, temp_phys)

   End Subroutine SHTns_ToSpectral

   Subroutine SHTns_ToPhysical(data_in, data_out)
       Type(rmcontainer4D), Intent(In) :: data_in(1:)
       Real*8, Intent(InOut) :: data_out(:,:,:)
       ! ingoing data has shape:
       !   data_out(lm)%data(l,r,imi,nf)
       !       lm = mode index (distributed)
       !        l = spherical harmonic (in-processor)
       !        r = radius (distributed)
       !      imi = real/imag parts
       !       nf = number of fields
       !
       ! outgoing data has shape:
       !   data_out(th,nrhs,lm)
       !       th = theta (in-processor)
       !     nrhs = number of RHS elements, stacked over: radius/real/imag/nfields
       !       lm = mode index (distributed)
       Complex*16, Allocatable :: temp_phys(:), temp_spec(:)
       Complex*16 :: ai, ar
       Integer :: odims(3), nm, nrhs, my_Nr
       Integer :: idims(4), nfield, rmn, rmx, mode, m, f, r, ind
       Real*8 :: norm

       ai = (0.0d0, 1.0d0) ! imaginary unit
       ar = (1.0d0, 0.0d0) ! regular one

       odims = shape(data_out)
       nm = odims(3)
       nrhs = odims(2)

       idims = shape(data_in(1)%data)
       nfield = idims(4)
       rmn = LBOUND(data_in(1)%data,2)
       rmx = UBOUND(data_in(1)%data,2)

       my_Nr = rmx - rmn + 1

       ! build temporary storage spaces
       allocate(temp_phys(1:n_theta), temp_spec(0:l_max))

       Do mode = 1, nm ! loop over lm modes

           m = m_values(mode) ! extract actual m value

           norm = STP_normalization(m) ! extract FFT normalizations, based on m

           Do f = 1, nfield ! number of fields
               Do r = rmn, rmx ! radius

                  ! package the input data into a complex array of size lmax-m+1
                  temp_spec(:) = (0.0d0, 0.0d0)
                  temp_spec(m:l_max) = ar*data_in(mode)%data(m:l_max,r,1,f) &
                                     + ai*data_in(mode)%data(m:l_max,r,2,f)

                  temp_phys(:) = (0.0d0, 0.0d0)
                  Call sh_to_spat_ml(SHTns_c, m, temp_spec(m:l_max), temp_phys, l_max)

                  ! 1-based indexing, except for radius, ordered: radius-real-imag-field
                  !     index = (r-rlo+1) + (imi-1)*Nr + (f-1)*Nr*2
                  ind = (r-rmn+1) + (f-1)*my_Nr*2         ! real part for this r/field
                  data_out(1:n_theta,ind,mode) = norm*real(temp_phys(1:n_theta))

                  ind = (r-rmn+1) + my_Nr + (f-1)*my_Nr*2 ! imag part for this r/field
                  data_out(1:n_theta,ind,mode) = norm*aimag(temp_phys(1:n_theta))

               Enddo
           Enddo
       Enddo
       deallocate(temp_phys, temp_spec)

   End Subroutine SHTns_ToPhysical

#else

Subroutine PtS_2d_dg(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...DGEMM
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Real*8 :: alpha, beta
    Integer :: m,nl

    alpha = 1.0d0
    beta = 0.0d0
    !Do m = 0, l_max,m_mod
    Do m = 1, n_m
            nl = l_max-m_values(m)+1
            CALL DGEMM('T','N',nl,nrhs,n_theta, alpha, p_lm(m)%data, n_theta,data_in , n_theta, beta,data_out(m)%data,nl)
            !CALL DGEMM('T','N',nl(m),nrhs,n_theta, alpha, p_lm(m)%data, n_theta,data_in , n_theta, beta,data_out(m)%data,nl(m))
    Enddo

End Subroutine PtS_2d_dg

Subroutine PtS_2d_dr(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...Direct
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Integer :: m,i,l
    ! Need to try a version 2 of this where the sum function is not invoked
    !Do m = 0, l_max,m_mod
    Do m = 1, n_m
        Do i = 1, nrhs
            Do l = m_values(m), l_max
                data_out(m)%data(l,i) = Sum(data_in(:,i)*p_lm(m)%data(:,l))
            Enddo
        Enddo
    Enddo

End Subroutine PtS_2d_dr

Subroutine PtS_2d_drp(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...Direct...Parity
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Real*8, Allocatable :: feven(:,:), fodd(:,:)
    Integer :: m,i,j,l,nt1

    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs))
    Allocate(fodd(1:n_theta/2,1:nrhs))
    feven(1:n_theta/2,:) = data_in(1:n_theta/2,:)
    fodd(1:n_theta/2,:) = data_in(1:n_theta/2,:)
    nt1 = n_theta+1
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i) = feven(j,i)+data_in(nt1-j,i)
             fodd(j,i) =  fodd(j,i)-data_in(nt1-j,i)
        Enddo
    enddo

    ! One possible way - involves a lot of indexing.
    Do m = 1, n_m
        Do i = 1, nrhs
            Do j = 1, n_l_odd(m)
                l = lvals(m)%odd(j)
                data_out(m)%data(l,i) = Sum(fodd(:,i)*p_lm_odd(m)%data(:,j))
            Enddo
        Enddo
    Enddo
    Do m = 1, n_m
        Do i = 1, nrhs
            Do j = 1, n_l_even(m)
                l = lvals(m)%even(j)
                data_out(m)%data(l,i) = Sum(feven(:,i)*p_lm_even(m)%data(:,j))
            Enddo
        Enddo
    Enddo

    ! A possible alternative is to simply store the odds up front and evens later
    ! This would require re-interleaving in the end
    !Do m = 1, n_m
    !    Do i = 1, nrhs
    !        Do j = 1, n_ell_odd(m)
    !            !l = lvals(m)%odd(j)
    !            data_out(m)%data(j,i) = Sum(fodd(:,i)*p_lm_odd(m)%data(:,j))
    !        Enddo
    !        offset = n_ell_odd(m)
    !        Do j = 1, n_ell_even(m)
    !            data_out(m)%data(j+offset,i) = Sum(feven(:,i)*p_lm_even(m)%data(:,j))
    !        Enddo
    !    Enddo
    !Enddo

    DeAllocate(feven,fodd)

End Subroutine PtS_2d_drp

Subroutine PtS_2d_dgp(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...DGEMM.... Parity
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),fodd(:,:), feven(:,:)
    Integer :: m,nt1,i,j,l, nt2



    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs))
    Allocate(fodd(1:n_theta/2,1:nrhs))
    feven(1:n_theta/2,:) = data_in(1:n_theta/2,:)
    fodd(1:n_theta/2,:) = data_in(1:n_theta/2,:)
    nt1 = n_theta+1
    nt2 = n_theta/2
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i) = feven(j,i)+data_in(nt1-j,i)
             fodd(j,i) =  fodd(j,i)-data_in(nt1-j,i)
        Enddo
    enddo



    alpha = 1.0d0
    beta = 0.0d0
    Do m = 1, n_m

            If (n_l_even(m) .gt. 0) then
                Allocate(temp(1:n_l_even(m),1:nrhs))
                CALL DGEMM('T','N',n_l_even(m),nrhs,nt2, alpha, p_lm_even(m)%data, nt2,feven , nt2, beta,temp,n_l_even(m))
                Do i =1, nrhs
                Do j = 1, n_l_even(m)
                    l = lvals(m)%even(j)
                    data_out(m)%data(l,i) = temp(j,i)
                Enddo
                Enddo
                DeAllocate(temp)
            Endif

            If (n_l_odd(m) .gt. 0) then
                Allocate(temp(1:n_l_odd(m),1:nrhs))
                CALL DGEMM('T','N',n_l_odd(m),nrhs,nt2, alpha, p_lm_odd(m)%data, nt2,fodd , nt2, beta,temp,n_l_odd(m))
                Do i = 1, nrhs
                Do j = 1, n_l_odd(m)
                    l = lvals(m)%odd(j)
                    data_out(m)%data(l,i) = temp(j,i)
                Enddo
                Enddo
                DeAllocate(temp)
            Endif


    Enddo

    DeAllocate(feven)
    DeAllocate(fodd)

End Subroutine PtS_2d_dgp

Subroutine PtS_2d_dgpv2(data_in, data_out)
    ! Physical-to_Spectral ..2D...DGEMM.... Parity
    ! Data_in is a simple 2-d array (theta is first index)
    ! version 2 --- data_in is a 3D array (ntheta,nrhs,n_m)
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(rmcontainer), Intent(InOut) :: data_out(1:)
    Real*8, Intent(In) :: data_in(:,:,:)
    Integer  :: nrhs
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),fodd(:,:,:), feven(:,:,:)
    Integer :: m,nl,nt1,i,j,l, nt2,ddims(3),k

    ddims = shape(data_in)
    n_m = ddims(3)
    nrhs = ddims(2)
    alpha = 1.0d0
    beta = 0.0d0
    if (.not. parity) then
        Do m =1, n_m
            nl = l_max-m_values(m)+1

            CALL DGEMM('T','N',nl,nrhs,n_theta, alpha, ip_lm(m)%data, &
                n_theta,data_in(:,:,m) , n_theta, beta,data_out(m)%data,nl)
        Enddo
    else
    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs,1:n_m))
    Allocate(fodd(1:n_theta/2,1:nrhs,1:n_m))
    feven(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    fodd(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    nt1 = n_theta+1
    nt2 = n_theta/2
    Do k = 1, n_m
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i,k) = feven(j,i,k)+data_in(nt1-j,i,k)
             fodd(j,i,k) =  fodd(j,i,k)-data_in(nt1-j,i,k)
        Enddo
    enddo
    Enddo



    Do m = 1, n_m

            If (n_l_even(m) .gt. 0) then
                Allocate(temp(1:n_l_even(m),1:nrhs))
                CALL DGEMM('T','N',n_l_even(m),nrhs,nt2, alpha, ip_lm_even(m)%data, nt2,feven(:,:,m) , nt2, beta,temp,n_l_even(m))
                Do i =1, nrhs
                Do j = 1, n_l_even(m)
                    l = lvals(m)%even(j)
                    data_out(m)%data(l,i) = temp(j,i)
                Enddo
                Enddo
                DeAllocate(temp)
            Endif

            If (n_l_odd(m) .gt. 0) then
                Allocate(temp(1:n_l_odd(m),1:nrhs))
                CALL DGEMM('T','N',n_l_odd(m),nrhs,nt2, alpha, ip_lm_odd(m)%data, nt2,fodd(:,:,m) , nt2, beta,temp,n_l_odd(m))
                Do i = 1, nrhs
                Do j = 1, n_l_odd(m)
                    l = lvals(m)%odd(j)
                    data_out(m)%data(l,i) = temp(j,i)
                Enddo
                Enddo
                DeAllocate(temp)
            Endif


    Enddo

    DeAllocate(feven)
    DeAllocate(fodd)
    endif
End Subroutine PtS_2d_dgpv2


Subroutine StP_2d_dgp(data_in, data_out)
    ! Physical-to_Spectral ..2D...DGEMM.... Parity
    ! Data_in is a spectral structure data_in(m)%data(l,i) ! i is radius or what have you
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(rmcontainer), Intent(In) :: data_in(:)
    Real*8, Intent(InOut) :: data_out(:,:,:)
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),temp2(:,:)
    Integer :: m,nl,nt1,i,j,l, nt2, ddims(3), nrhs
    ddims = shape(data_out)
    n_m = ddims(3)
    nrhs = ddims(2)

    alpha = 1.0d0
    beta = 0.0d0
    if (.not. parity) then

        Do m = 1, n_m
            nl = l_max-m_values(m)+1
            CALL DGEMM('T','N',n_theta,nrhs,nl, alpha, p_lm(m)%data,  &
                nl,data_in(m)%data , nl, beta,data_out,n_theta)
        Enddo

    else
    !////////////////////////////////////
    ! In progress
    nt1 = n_theta+1
    nt2 = n_theta/2
    data_out(:,:,:) = 0.0d0
    Allocate(temp(1:nt2,1:nrhs))
    ! Solve for odd and even functions
    Do m = 1, n_m

        If (n_l_even(m) .gt. 0) then
            ! This feels unnecessarily clunky.  Might want to consider storing spectral data as even/odd modes.
            ! Just get it running for now
            Allocate(temp2(1:n_l_even(m),1:nrhs))
            Do i =1, nrhs
            Do j = 1, n_l_even(m)
                l = lvals(m)%even(j)
                temp2(j,i) = data_in(m)%data(l,i)
            Enddo
            Enddo
            CALL DGEMM('T','N',nt2,nrhs,n_l_even(m), alpha, p_lm_even(m)%data, n_l_even(m),temp2 , n_l_even(m), beta,temp,nt2)
            data_out(1:nt2,:,m) = temp    ! store symmetric part in data_out
            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(nt1-j,i,m) = temp(j,i)    ! reflect even modes about equator
                Enddo
            Enddo
            DeAllocate(temp2)
        Endif

        If (n_l_odd(m) .gt. 0) then
            Allocate(temp2(1:n_l_odd(m),1:nrhs))
            Do i =1, nrhs
            Do j = 1, n_l_odd(m)
                l = lvals(m)%odd(j)
                temp2(j,i) = data_in(m)%data(l,i)
            Enddo
            Enddo
            CALL DGEMM('T','N',nt2,nrhs,n_l_odd(m), alpha, p_lm_odd(m)%data, n_l_odd(m),temp2 , n_l_odd(m), beta,temp,nt2)
            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(j,i,m) = data_out(j,i,m)+temp(j,i)
                    data_out(nt1-j,i,m) = data_out(nt1-j,i,m)-temp(j,i)    ! antisymmetric about equator
                Enddo
            Enddo
            DeAllocate(temp2)
        Endif
    Enddo


    ! Note - not sure if it's faster to make a variable named nt2j1 = nt2-j+1 or just let it compute on the fly
    DeAllocate(temp)
    Endif

End Subroutine StP_2d_dgp



Subroutine PtS_2d_dr2(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...Direct
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Integer :: m,i,l,th
    ! Does not use the sum function
    !Do m = 0, l_max
    Do m = 1,n_m
        Do i = 1, nrhs
            Do l = m_values(m), l_max
                data_out(m)%data(l,i) =   data_in(1,i)*p_lm(m)%data(1,l)
                Do th = 2, n_theta
                    data_out(m)%data(l,i) = data_out(m)%data(l,i)+data_in(th,i)*p_lm(m)%data(th,l)
                Enddo
            Enddo
        Enddo
    Enddo

End Subroutine PtS_2d_dr2

Subroutine PtS_2d_dr3(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...Direct
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Integer :: m,i,l,th
    ! Does not use the sum function
    ! Exhanged l and i loops
    !Do m = 0, l_max
    Do m = 1,n_m
        Do l = m_values(m), l_max
        Do i = 1, nrhs
            !Do l = m_values(m), l_max
                data_out(m)%data(l,i) =   data_in(1,i)*p_lm(m)%data(1,l)
                Do th = 2, n_theta
                    data_out(m)%data(l,i) = data_out(m)%data(l,i)+data_in(th,i)*p_lm(m)%data(th,l)
                Enddo
            Enddo
        Enddo
    Enddo

End Subroutine PtS_2d_dr3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!////////////////////
!//////   4-D routines
Subroutine PtS_4d_dgpv2(data_in, data_out)
    ! Physical-to_Spectral ..2D...DGEMM.... Parity
    ! Data_in is a simple 2-d array (theta is first index)
    ! version 2 --- data_in is a 3D array (ntheta,nrhs,n_m)
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(rmcontainer4D), Intent(InOut) :: data_out(1:)
    Real*8, Intent(In) :: data_in(:,:,:)
    Integer  :: nrhs
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),fodd(:,:,:), feven(:,:,:)
    Integer :: m,nl,nt1,i,j,l, nt2,ddims(3),k
    Integer :: oddims(4), nfield
    Integer :: rmn, rmx, nr, f, imi, istart, iend

    oddims = shape(data_out(1)%data)
    nfield = oddims(4)
    rmn = LBOUND(data_out(1)%data,2)
    rmx = UBOUND(data_out(1)%data,2)
    nr = rmx-rmn+1

    ddims = shape(data_in)
    n_m = ddims(3)
    nrhs = ddims(2)
    alpha = 1.0d0
    beta = 0.0d0
    if (.not. parity) then
        Do m =1, n_m
            nl = l_max-m_values(m)+1

            CALL DGEMM('T','N',nl,nrhs,n_theta, alpha, ip_lm(m)%data, &
                n_theta,data_in(:,:,m) , n_theta, beta,data_out(m)%data,nl)
        Enddo
    else
    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs,1:n_m))
    Allocate(fodd(1:n_theta/2,1:nrhs,1:n_m))
    feven(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    fodd(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    nt1 = n_theta+1
    nt2 = n_theta/2
    Do k = 1, n_m
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i,k) = feven(j,i,k)+data_in(nt1-j,i,k)
             fodd(j,i,k) =  fodd(j,i,k)-data_in(nt1-j,i,k)
        Enddo
    enddo
    Enddo



    Do m = 1, n_m

            If (n_l_even(m) .gt. 0) then
                Allocate(temp(1:n_l_even(m),1:nrhs))
                CALL DGEMM('T','N',n_l_even(m),nrhs,nt2, alpha, ip_lm_even(m)%data, nt2,feven(:,:,m) , nt2, beta,temp,n_l_even(m))
                !Do i =1, nrhs
                !Do j = 1, n_l_even(m)
                !    l = lvals(m)%even(j)
                !    data_out(m)%data(l,i) = temp(j,i)
                !Enddo
                !Enddo

                istart = 1
                iend = nr
                Do f = 1, nfield
                Do imi =1, 2
                    Do j = 1, n_l_even(m)
                        l = lvals(m)%even(j)
                        data_out(m)%data(l,rmn:rmx,imi,f) = temp(j,istart:iend)
                    Enddo
                    istart = istart+nr
                    iend = iend+nr
                Enddo
                Enddo


                DeAllocate(temp)
            Endif

            If (n_l_odd(m) .gt. 0) then
                Allocate(temp(1:n_l_odd(m),1:nrhs))
                CALL DGEMM('T','N',n_l_odd(m),nrhs,nt2, alpha, ip_lm_odd(m)%data, nt2,fodd(:,:,m) , nt2, beta,temp,n_l_odd(m))
                !Do i = 1, nrhs
                !Do j = 1, n_l_odd(m)
                !    l = lvals(m)%odd(j)
                !    data_out(m)%data(l,i) = temp(j,i)
                !Enddo
                !Enddo
                istart = 1
                iend = nr
                Do f = 1, nfield
                Do imi =1, 2
                    Do j = 1, n_l_odd(m)
                        l = lvals(m)%odd(j)
                        data_out(m)%data(l,rmn:rmx,imi,f) = temp(j,istart:iend)
                    Enddo
                    istart = istart+nr
                    iend = iend+nr
                Enddo
                Enddo
                DeAllocate(temp)
            Endif


    Enddo

    DeAllocate(feven)
    DeAllocate(fodd)
    endif
End Subroutine PtS_4d_dgpv2


Subroutine StP_4d_dgp(data_in, data_out)
    ! Physical-to_Spectral ..2D...DGEMM.... Parity
    ! Data_in is a spectral structure data_in(m)%data(l,i) ! i is radius or what have you
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(rmcontainer4d), Intent(In) :: data_in(:)
    Real*8, Intent(InOut) :: data_out(:,:,:)
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),temp2(:,:)
    Integer :: m,nl,nt1,i,j,l, nt2, ddims(3), nrhs
    Integer :: nfield, rmn, rmx, nr, oddims(4),imi,f,iend,istart

    ddims = shape(data_out)
    n_m = ddims(3)
    nrhs = ddims(2)

    oddims = shape(data_in(1)%data)
    nfield = oddims(4)

    rmn = LBOUND(data_in(1)%data,2)
    rmx = UBOUND(data_in(1)%data,2)
    nr = rmx-rmn+1


    alpha = 1.0d0
    beta = 0.0d0
    if (.not. parity) then
        ! unsure if this works with the 4d layout, but it should
        Do m = 1, n_m
            nl = l_max-m_values(m)+1
            CALL DGEMM('T','N',n_theta,nrhs,nl, alpha, p_lm(m)%data,  &
                nl,data_in(m)%data , nl, beta,data_out,n_theta)
        Enddo

    else
    !////////////////////////////////////
    ! In progress
    nt1 = n_theta+1
    nt2 = n_theta/2
    data_out(:,:,:) = 0.0d0
    Allocate(temp(1:nt2,1:nrhs))
    ! Solve for odd and even functions
    Do m = 1, n_m

        If (n_l_even(m) .gt. 0) then
            ! This feels unnecessarily clunky.  Might want to consider storing spectral data as even/odd modes.
            ! Just get it running for now
            Allocate(temp2(1:n_l_even(m),1:nrhs))
            !Do i =1, nrhs
            !Do j = 1, n_l_even(m)
            !    l = lvals(m)%even(j)
            !    temp2(j,i) = data_in(m)%data(l,i)
            !Enddo
            !Enddo

            istart = 1
            iend = nr
            Do f = 1, nfield
            Do imi =1, 2
                Do j = 1, n_l_even(m)
                    l = lvals(m)%even(j)
                    !!data_out(m)%data(l,rmn:rmx,imi,f) = temp(j,istart:iend)
                    temp2(j,istart:iend) = data_in(m)%data(l,rmn:rmx,imi,f)
                Enddo
                istart = istart+nr
                iend = iend+nr
            Enddo
            Enddo


            CALL DGEMM('T','N',nt2,nrhs,n_l_even(m), alpha, p_lm_even(m)%data, n_l_even(m),temp2 , n_l_even(m), beta,temp,nt2)
            data_out(1:nt2,:,m) = temp    ! store symmetric part in data_out
            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(nt1-j,i,m) = temp(j,i)    ! reflect even modes about equator
                Enddo
            Enddo
            DeAllocate(temp2)
        Endif

        If (n_l_odd(m) .gt. 0) then
            Allocate(temp2(1:n_l_odd(m),1:nrhs))
            !Do i =1, nrhs
            !Do j = 1, n_l_odd(m)
            !    l = lvals(m)%odd(j)
            !    temp2(j,i) = data_in(m)%data(l,i)
            !Enddo
            !Enddo


            istart = 1
            iend = nr
            Do f = 1, nfield
            Do imi =1, 2
                Do j = 1, n_l_odd(m)
                    l = lvals(m)%odd(j)
                    !!data_out(m)%data(l,rmn:rmx,imi,f) = temp(j,istart:iend)
                    temp2(j,istart:iend) = data_in(m)%data(l,rmn:rmx,imi,f)
                Enddo
                istart = istart+nr
                iend = iend+nr
            Enddo
            Enddo


            CALL DGEMM('T','N',nt2,nrhs,n_l_odd(m), alpha, p_lm_odd(m)%data, n_l_odd(m),temp2 , n_l_odd(m), beta,temp,nt2)
            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(j,i,m) = data_out(j,i,m)+temp(j,i)
                    data_out(nt1-j,i,m) = data_out(nt1-j,i,m)-temp(j,i)    ! antisymmetric about equator
                Enddo
            Enddo
            DeAllocate(temp2)
        Endif
    Enddo


    ! Note - not sure if it's faster to make a variable named nt2j1 = nt2-j+1 or just let it compute on the fly
    DeAllocate(temp)
    Endif

End Subroutine StP_4d_dgp


!/////////////////////////////////////////////////////////////////////////////
! These 4-D routines use a higher dimension of temp array
Subroutine StP_4d_dgp2(data_in, data_out)
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(rmcontainer4d), Intent(In) :: data_in(:)
    Real*8, Intent(InOut) :: data_out(:,:,:)
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),temp2(:,:), temp3(:,:,:,:)
    Integer :: m,nl,nt1,i,j,l, nt2, ddims(3), nrhs, jstart, jend,r,lstart
    Integer :: nfield, rmn, rmx, oddims(4),imi,f,iend,istart

    ddims = shape(data_out)
    n_m = ddims(3)
    nrhs = ddims(2)

    oddims = shape(data_in(1)%data)
    nfield = oddims(4)

    rmn = LBOUND(data_in(1)%data,2)
    rmx = UBOUND(data_in(1)%data,2)



    alpha = 1.0d0
    beta = 0.0d0
    if (.not. parity) then
        Do m = 1, n_m
            nl = l_max-m_values(m)+1
            CALL DGEMM('T','N',n_theta,nrhs,nl, alpha, p_lm(m)%data,  &
                nl,data_in(m)%data , nl, beta,data_out,n_theta)
        Enddo

    else
    !////////////////////////////////////
    ! In progress
    nt1 = n_theta+1
    nt2 = n_theta/2
    data_out(:,:,:) = 0.0d0
    Allocate(temp(1:nt2,1:nrhs))
    ! Solve for odd and even functions
    Do m = 1, n_m

        If (n_l_even(m) .gt. 0) then

            jend = n_l_even(m)
            !write(6,*)jend,rmn,rmx,nfield
            Allocate(temp3(1:jend,rmn:rmx,1:2,1:nfield))
            Do f = 1, nfield
            Do imi = 1, 2
            Do r = rmn, rmx
            Do j = 1,jend
                l =  lvalsi(m)%even(j)
                temp3(j,r,imi,f) = data_in(m)%data(l,r,imi,f)
            Enddo
            Enddo
            Enddo
            Enddo


            CALL DGEMM('T','N',nt2,nrhs,n_l_even(m), alpha, p_lm_even(m)%data, n_l_even(m),temp3 , n_l_even(m), beta,temp,nt2)

            data_out(1:nt2,1:nrhs,m) = temp(1:nt2,1:nrhs)    ! store symmetric part in data_out

            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(nt1-j,i,m) = temp(j,i)    ! reflect even modes about equator
                Enddo
            Enddo
            DeAllocate(temp3)
        Endif

        If (n_l_odd(m) .gt. 0) then

            jend = n_l_odd(m)

            Allocate(temp3(1:jend,rmn:rmx,1:2,1:nfield))
            Do f = 1, nfield
            Do imi = 1, 2
            Do r = rmn, rmx
            Do j = 1,jend
                l =  lvalsi(m)%odd(j)
                temp3(j,r,imi,f) = data_in(m)%data(l,r,imi,f)
            Enddo
            Enddo
            Enddo
            Enddo



            CALL DGEMM('T','N',nt2,nrhs,n_l_odd(m), alpha, p_lm_odd(m)%data, n_l_odd(m),temp3 , n_l_odd(m), beta,temp,nt2)
            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(j,i,m) = data_out(j,i,m)+temp(j,i)
                    data_out(nt1-j,i,m) = data_out(nt1-j,i,m)-temp(j,i)    ! antisymmetric about equator
                Enddo
            Enddo
            DeAllocate(temp3)
        Endif
    Enddo


    ! Note - not sure if it's faster to make a variable named nt2j1 = nt2-j+1 or just let it compute on the fly
    DeAllocate(temp)
    Endif

End Subroutine StP_4d_dgp2


Subroutine PtS_4d_dgpv3(data_in, data_out)
    Implicit None
    Type(rmcontainer4D), Intent(InOut) :: data_out(1:)
    Real*8, Intent(In) :: data_in(:,:,:)
    Integer  :: nrhs
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:,:,:),fodd(:,:,:), feven(:,:,:)
    Integer :: m,nl,nt1,i,j,l, nt2,ddims(3),k
    Integer :: oddims(4), nfield
    Integer :: rmn, rmx, f, imi, istart, iend, jend,r

    oddims = shape(data_out(1)%data)
    nfield = oddims(4)
    rmn = LBOUND(data_out(1)%data,2)
    rmx = UBOUND(data_out(1)%data,2)


    ddims = shape(data_in)
    n_m = ddims(3)
    nrhs = ddims(2)
    alpha = 1.0d0
    beta = 0.0d0
    if (.not. parity) then
        Do m =1, n_m
            nl = l_max-m_values(m)+1

            CALL DGEMM('T','N',nl,nrhs,n_theta, alpha, ip_lm(m)%data, &
                n_theta,data_in(:,:,m) , n_theta, beta,data_out(m)%data,nl)
        Enddo
    else
    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs,1:n_m))
    Allocate(fodd(1:n_theta/2,1:nrhs,1:n_m))
    feven(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    fodd(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    nt1 = n_theta+1
    nt2 = n_theta/2
    Do k = 1, n_m
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i,k) = feven(j,i,k)+data_in(nt1-j,i,k)
             fodd(j,i,k) =  fodd(j,i,k)-data_in(nt1-j,i,k)
        Enddo
    enddo
    Enddo



    Do m = 1, n_m

            If (n_l_even(m) .gt. 0) then
                Allocate(temp(1:n_l_even(m),rmn:rmx,1:2,1:nfield))
                CALL DGEMM('T','N',n_l_even(m),nrhs,nt2, alpha, ip_lm_even(m)%data, nt2,feven(:,:,m) , nt2, beta,temp,n_l_even(m))


                jend = n_l_even(m)
                Do f = 1, nfield
                Do imi = 1, 2
                Do r = rmn, rmx
                Do j = 1,jend
                    l =  lvalsi(m)%even(j)
                    data_out(m)%data(l,r,imi,f)=temp(j,r,imi,f)
                Enddo
                Enddo
                Enddo
                Enddo

                DeAllocate(temp)
            Endif

            If (n_l_odd(m) .gt. 0) then
                Allocate(temp(1:n_l_odd(m),rmn:rmx,1:2,1:nfield))
                CALL DGEMM('T','N',n_l_odd(m),nrhs,nt2, alpha, ip_lm_odd(m)%data, nt2,fodd(:,:,m) , nt2, beta,temp,n_l_odd(m))


                jend = n_l_odd(m)
                Do f = 1, nfield
                Do imi = 1, 2
                Do r = rmn, rmx
                Do j = 1,jend
                    l =  lvalsi(m)%odd(j)
                    data_out(m)%data(l,r,imi,f)=temp(j,r,imi,f)
                Enddo
                Enddo
                Enddo
                Enddo


                DeAllocate(temp)
            Endif


    Enddo

    DeAllocate(feven)
    DeAllocate(fodd)
    endif
End Subroutine PtS_4d_dgpv3
#endif

End Module Legendre_Transforms
