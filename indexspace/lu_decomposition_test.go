/*
 File		:	$File: //depot/stillwater-sc/godfa/indexspace/lu_decomposition_test.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	29 January 2017

 Source Control Information:
 Version	:	$Revision: #1 $
 Latest		:	$Date: 2017/01/29 $
 Location	:	$Id: //depot/stillwater-sc/godfa/indexspace/lu_decomposition_test.go#1 $

 Organization:
		Stillwater Supercomputing, Inc.
		P.O Box 720
		South Freeport, ME 04078-0720

Copyright (c) 2006-2017 E. Theodore L. Omtzigt.  All rights reserved.

Licence      : Stillwater license as defined in this directory

 */
package indexspace

import (
	"testing"
)

func TestLU_with_pivoting(t *testing.T) {
	A := NewMatrixShell()
	A, _ = A.AddRow(Vector{0,0,1})
	A, _ = A.AddRow(Vector{0,1,0})
	A, _ = A.AddRow(Vector{1,0,0})

	LU, pivots, fullrank := LU_with_pivoting(A,false)
	if !fullrank {
		t.Error("Result should have been fullrank")
	}
	pivots_ref := "  2   1   2 "
	if pivots.String() != pivots_ref {
		t.Error("Pivoting Result is incorrect")
		t.Logf("A:\n%s\n", A)
		t.Logf("LU:\n%s\n", LU)
		t.Logf("Pivots: %s\n", pivots)
		t.Logf("Pivots: %x\n", pivots)
		t.Logf("Pivots: %x\n", pivots_ref)
	}

	// start testing LU solve
	b := NewVector(3)
	x := NewVector(3)

	A = NewMatrix(3,3)
	A[0][0] = 1.0; A[0][1] = 2.0; A[0][2] = 3.0
	A[1][0] = 1.0; A[1][1] = 2.0; A[1][2] = 3.0
	A[2][0] = 1.0; A[2][1] = 2.0; A[2][2] = 3.0
	LU, pivots, fullrank = LU_with_pivoting(A,false)
	if fullrank {
		t.Error("Result should have been singular")
	}
	pivots_ref = "  0   1   2 "
	if pivots.String() != pivots_ref {
		t.Error("Pivoting Result is incorrect")
		t.Logf("A:\n%s\n", A)
		t.Logf("LU:\n%s\n", LU)
		t.Logf("Pivots: %s\n", pivots)
		t.Logf("Pivots: %x\n", pivots)
		t.Logf("Pivots: %x\n", pivots_ref)
	}

	A[0][0] = 1.0; A[0][1] = 2.0; A[0][2] = 3.0
	A[1][0] = 0.0; A[1][1] = 1.0; A[1][2] = 2.0
	A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 1.0
	LU, pivots, fullrank = LU_with_pivoting(A,false)
	if !fullrank {
		t.Error("Result should have been fullrank")
	}
	pivots_ref = "  0   1   2 "
	if pivots.String() != pivots_ref {
		t.Error("Pivoting Result is incorrect")
		t.Logf("A:\n%s\n", A)
		t.Logf("LU:\n%s\n", LU)
		t.Logf("Pivots: %s\n", pivots)
		t.Logf("Pivots: %x\n", pivots)
		t.Logf("Pivots: %x\n", pivots_ref)
	}
	b[0] = 1.0; b[1] = 1.0; b[2] = 1.0
	x = LUsolve(LU, pivots, b)
	x_ref := "  0.00  -1.00   1.00 "
	if x.String() != x_ref {
		t.Error("Solution Result is incorrect")
		t.Logf("b: %s\n", b)
		t.Logf("x: %s\n", x)
		t.Logf("x_ref: %x\n", x_ref)
		t.Logf("x    : %x\n", x)
	}

	// Scramble the matrix so that pivoting is necessary: should still yield the same result but with pivots
	A[0][0] = 0.0; A[0][1] = 0.0; A[0][2] = 1.0
	A[1][0] = 0.0; A[1][1] = 1.0; A[1][2] = 2.0
	A[2][0] = 1.0; A[2][1] = 2.0; A[2][2] = 3.0
	LU, pivots, fullrank = LU_with_pivoting(A,false)
	if !fullrank {
		t.Error("Result should have been fullrank")
	}
	b[0] = 1.0; b[1] = 1.0; b[2] = 1.0
	x = LUsolve(LU, pivots, b)
	x_ref = "  0.00  -1.00   1.00 "
	if x.String() != x_ref {
		t.Error("Solution Result is incorrect")
		t.Logf("b: %s\n", b)
		t.Logf("x: %s\n", x)
		t.Logf("x_ref: %x\n", x_ref)
		t.Logf("x    : %x\n", x)
	}

	// a simple stencil matrix
	A[0][0] = 1.0; A[0][1] = 2.0; A[0][2] = 0.0
	A[1][0] = 2.0; A[1][1] = 1.0; A[1][2] = 2.0
	A[2][0] = 0.0; A[2][1] = 2.0; A[2][2] = 1.0
	LU, pivots, fullrank = LU_with_pivoting(A,false)
	if !fullrank {
		t.Error("Result should have been fullrank")
	}
	x = LUsolve(LU, pivots, b)
	x_ref = "  0.14   0.43   0.14 "
	if x.String() != x_ref {
		t.Error("Solution Result is incorrect")
		t.Logf("b: %s\n", b)
		t.Logf("x: %s\n", x)
		t.Logf("x_ref: %x\n", x_ref)
		t.Logf("x    : %x\n", x)
	}
}