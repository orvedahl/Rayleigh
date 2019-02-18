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

Module Sorting
    ! This module performs basic sorting, these are not optimized for speed/efficiency

    Implicit None

    Interface Sort
        Module Procedure Wrap_Bubble_Sort_i1d
        Module Procedure Wrap_Bubble_Sort_r1d
    End Interface

Contains

    ! wrapper to bubble sort for 1D integer array
    Subroutine Wrap_Bubble_Sort_i1d(x, ignore, ignore_sign, tolerance)
        ! x is the array to be sorted
        ! instances where "x==ignore" will not be included in the sort
        !    for example, x = [9, 2, 3, 5, -1, -1, -1, -1] with ignore=-1
        !    will produce x = [2, 3, 5, 9, -1, -1, -1, -1]
        ! ignore_sign indicates if sorting is done on the abolute value or the actual value
        !    for example, x = [9, 2, 3, -5, -1, -1, -1, -1] with ignore=-1, ignore_sign=.true.
        !    will produce x = [2, 3, -5, 9, -1, -1, -1, -1]
        ! tolerance indicates the degree to which "x==ignore", i.e., "x - ignore < tol"
        Integer, Intent(InOut) :: x(:)
        Integer, Intent(In) :: ignore
        Logical, Intent(In), Optional :: ignore_sign
        Real*8, Intent(In), Optional :: tolerance
        Integer, Allocatable :: y(:)
        Logical, Allocatable :: mask(:)
        Real*8 :: tol
        Logical :: ign_sign
        Integer :: i, j, npts, N_valid
        If (present(ignore_sign)) Then
            ign_sign = ignore_sign
        Else
            ign_sign = .false.
        EndIf
        If (present(tolerance)) Then
            tol = tolerance
        Else
            tol = 1.0d-8
        EndIf

        npts = Size(x, dim=1)
        Allocate(mask(1:npts))

        Where (abs(x-ignore) .lt. tol) ! find location of values /= ignore_value
            mask = .false.
        ElseWhere
            mask = .true.
        EndWhere

        N_valid = Count(mask) ! number of valid locations
        If (N_valid .lt. 1) Then
            DeAllocate(mask)
            Return
        EndIf

        Allocate(y(1:N_valid))

        j = 1
        Do i=1,npts
           If (mask(i)) Then ! extract valid numbers
               y(j) = x(i)
               j = j + 1
           EndIf
        EndDo

        Call Bubble_Sort_i1d(y, ignore_sign=ign_sign) ! sort valid numbers

        x(1:N_valid) = y(1:N_valid) ! store valid numbers and append ignored values
        x(N_valid+1:) = ignore

        DeAllocate(mask, y)

    End Subroutine Wrap_Bubble_Sort_i1d
   
    ! wrapper to bubble sort for 1D real array
    Subroutine Wrap_Bubble_Sort_r1d(x, ignore, ignore_sign, tolerance)
        ! x is the array to be sorted
        ! instances where "x==ignore" will not be included in the sort
        !    for example, x = [9.1, 2.4, 3.1, 5., -3, -3, -3] with ignore=-3d0
        !    will produce x = [2.4, 3.1, 5., 9.1, -3, -3, -3]
        ! ignore_sign indicates if sorting is done on the abolute value or the actual value
        !    for example, x = [9.1, 2.1, 3.1, -5., -3, -3] with ignore=-3d0, ignore_sign=.true.
        !    will produce x = [2.1, 3.1, -5., 9.1, -3, -3]
        ! tolerance indicates the degree to which "x==ignore", i.e., "x - ignore < tol"
        Real*8, Intent(InOut) :: x(:)
        Real*8, Intent(In) :: ignore
        Logical, Intent(In), Optional :: ignore_sign
        Real*8, Intent(In), Optional :: tolerance
        Real*8, Allocatable :: y(:)
        Logical, Allocatable :: mask(:)
        Real*8 :: tol
        Logical :: ign_sign
        Integer :: i, j, npts, N_valid
        If (present(ignore_sign)) Then
            ign_sign = ignore_sign
        Else
            ign_sign = .false.
        EndIf
        If (present(tolerance)) Then
            tol = tolerance
        Else
            tol = 1.0d-8
        EndIf

        npts = Size(x, dim=1)
        Allocate(mask(1:npts))

        Where (abs(x-ignore) .lt. tol) ! find location of values /= ignore_value
            mask = .false.
        ElseWhere
            mask = .true.
        EndWhere

        N_valid = Count(mask) ! number of valid locations
        If (N_valid .lt. 1) Then
            DeAllocate(mask)
            Return
        EndIf

        Allocate(y(1:N_valid))

        j = 1
        Do i=1,npts
           If (mask(i)) Then ! extract valid numbers
               y(j) = x(i)
               j = j + 1
           EndIf
        EndDo

        Call Bubble_Sort_r1d(y, ignore_sign=ign_sign) ! sort valid numbers

        x(1:N_valid) = y(1:N_valid) ! store valid numbers and append ignored values
        x(N_valid+1:) = ignore

        DeAllocate(mask, y)

    End Subroutine Wrap_Bubble_Sort_r1d

    ! bubble sort for 1D integer array
    Subroutine Bubble_Sort_i1d(x, ignore_sign)
        ! x is the array to be sorted
        ! ignore_sign indicates if sorting is done on the abolute value or the actual value
        !    for example, x = [9, 2, 3, -5, -1] ignore_sign=.true.
        !    will produce x = [-1, 2, 3, -5, 9]
        Integer, Intent(InOut) :: x(:)
        Logical, Intent(In), Optional :: ignore_sign
        Logical :: ign_sign
        Integer :: npts, i, j
        Real*8 :: tmpv
        If (present(ignore_sign)) Then
            ign_sign = ignore_sign
        Else
            ign_sign = .false.
        EndIf
        npts = Size(x, dim=1)

        If (ign_sign) Then
            Do i=1,npts
               Do j=1,npts-1,1
                  If (abs(x(j)) > abs(x(j+1))) Then ! compare magnitudes
                     tmpv = x(j)          ! swap entries
                     x(j) = x(j+1)
                     x(j+1) = tmpv
                  EndIf
               EndDo
            EndDo
        Else
            Do i=1,npts
               Do j=1,npts-1,1
                  If (x(j) > x(j+1)) Then ! compare actual values
                     tmpv = x(j)          ! swap entries
                     x(j) = x(j+1)
                     x(j+1) = tmpv
                  EndIf
               EndDo
            EndDo
        EndIf
    End Subroutine Bubble_Sort_i1d

    ! bubble sort for 1D real array
    Subroutine Bubble_Sort_r1d(x, ignore_sign)
        ! x is the array to be sorted
        ! ignore_sign indicates if sorting is done on the abolute value or the actual value
        !    for example, x = [9.1, 2.1, 3.1, -5., -3] with ignore_sign=.true.
        !    will produce x = [2.1, -3, 3.1, -5., 9.1]
        Real*8, Intent(InOut) :: x(:)
        Logical, Intent(In), Optional :: ignore_sign
        Logical :: ign_sign
        Integer :: npts, i, j
        Real*8 :: tmpv
        If (present(ignore_sign)) Then
            ign_sign = ignore_sign
        Else
            ign_sign = .false.
        EndIf
        npts = Size(x, dim=1)

        If (ign_sign) Then
            Do i=1,npts
               Do j=1,npts-1,1
                  If (abs(x(j)) > abs(x(j+1))) Then ! compare magnitudes
                     tmpv = x(j)          ! swap entries
                     x(j) = x(j+1)
                     x(j+1) = tmpv
                  EndIf
               EndDo
            EndDo
        Else
            Do i=1,npts
               Do j=1,npts-1,1
                  If (x(j) > x(j+1)) Then ! compare actual values
                     tmpv = x(j)          ! swap entries
                     x(j) = x(j+1)
                     x(j+1) = tmpv
                  EndIf
               EndDo
            EndDo
        EndIf
    End Subroutine Bubble_Sort_r1d

End Module Sorting
