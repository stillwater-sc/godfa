/*
 File		:	$File: //depot/stillwater-sc/godfa/linalg/lu_decomposition.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	29 January 2017

 Source Control Information:
 Version	:	$Revision: #1 $
 Latest		:	$Date: 2017/01/29 $
 Location	:	$Id: //depot/stillwater-sc/godfa/linalg/lu_decomposition.go#1 $

 Organization:
		Stillwater Supercomputing, Inc.
		P.O Box 720
		South Freeport, ME 04078-0720

Copyright (c) 2006-2017 E. Theodore L. Omtzigt.  All rights reserved.

Licence      : Stillwater license as defined in this directory

 */
package la

import (
	"math"
)


/*
Straight forward LU decomposition with pivoting
 * @param A
 * @param bInplace
 * @returns {{LU: *, P: Array, fullRank: boolean}}
 * @constructor
*/
func LU_with_pivoting (A Matrix, bInplace bool) (Matrix, IntVector, bool) {
	var absAjk, Akk float64
	var Ak, Ai Vector
	var max float64
	var epsilon float64 = 1e-12
	var fullRank bool = true
	var n int = len(A)
	var n1 int = n-1
	var P IntVector = NewIntVector(n)
	var Pk int
	if !bInplace {
		A = A.Clone()
	}

	for k := 0; k < n; k++ {
		Pk = k
		Ak = A[k]
		max = math.Abs(Ak[k])
		for j := k + 1; j < n; j++ {
			absAjk = math.Abs(A[j][k])
			if max < absAjk {
				max = absAjk
				Pk = j
			}
		}
		P[k] = Pk
		if (Pk != k) {
			A[k] = A[Pk]
			A[Pk] = Ak
			Ak = A[k]
		}

		Akk = Ak[k]

		if math.Abs(Akk) > epsilon {
			for i := k + 1; i < n; i++ {
				A[i][k] /= Akk
			}
		} else {
			fullRank = false
		}

		var j int
		for i := k + 1; i < n; i++ {
			Ai = A[i]
			for j = k + 1; j < n1; j++ {
				Ai[j] -= Ai[k] * Ak[j]
				j++
				Ai[j] -= Ai[k] * Ak[j]
			}
			if j == n1 {
				Ai[j] -= Ai[k] * Ak[j]
			}
		}
	}

	return A, P, fullRank
}

/*
LUSolve computes the x in Ax = b
 */
func LUsolve(LU Matrix, P IntVector, b Vector) (x Vector) {
	var i, j int
	var n int = LU.Rows()
	var LUi Vector
	var Pi int

	x = b.Clone()
	for i = n-1; i != -1; i-- {
		x[i] = b[i];
	}
	for i = 0; i < n; i++ {
		Pi = P[i];
		if Pi != i {
			tmp := x[i]
			x[i] = x[Pi]
			x[Pi] = tmp
		}

		LUi = LU[i]
		for j = 0; j < i; j++ {
			x[i] = x[i] - (x[j] * LUi[j])
		}
	}

	for i = n - 1; i >= 0; i-- {
		LUi = LU[i]
		for j = i + 1; j < n; j++ {
			x[i] = x[i] - (x[j] * LUi[j])
		}
		x[i] = x[i] / LUi[i]
	}

	return x;
}
