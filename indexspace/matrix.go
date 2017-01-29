/*
 File		:	$File: //depot/stillwater-sc/godfa/indexspace/matrix.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	29 January 2017

 Source Control Information:
 Version	:	$Revision: #1 $
 Latest		:	$Date: 2017/01/29 $
 Location	:	$Id: //depot/stillwater-sc/godfa/indexspace/matrix.go#1 $

 Organization:
		Stillwater Supercomputing, Inc.
		P.O Box 720
		South Freeport, ME 04078-0720

Copyright (c) 2006-2017 E. Theodore L. Omtzigt.  All rights reserved.

Licence      : Stillwater license as defined in this directory

 */
package indexspace

import (
	"fmt"
	"errors"
)

/*
Matrix represents an [m][n] matrix

The Matrix type is intended for small matrix manipulations, and providing a tailored interface to geometric questions
that arise in the manipulation of index spaces and domains of computation
*/
type Matrix [][]float64

// Clone creates a deep copy of the receiver matrix
func (m Matrix) Clone() Matrix {
	newMatrix := make([][]float64, m.Rows())
	for i := 0; i < m.Rows(); i++ {
		if newMatrix[i] == nil {
			newMatrix[i] = make([]float64, m.Rows())
		}
		for j := 0; j < m.Columns(); j++ {
			newMatrix[i][j] = m[i][j]
		}
	}
	return newMatrix
}

// Rows returns the number of rows in the Matrix
func (m Matrix) Rows() int {
	return len(m)
}

// Columns returns the number of columns in the Matrix
func (m Matrix) Columns() int {
	if m.Rows() == 0 {
		return 0
	}
	return len(m[0])
}

// Dimensionality returns the dimension of the Matrix
func (mat Matrix) Dimensionality() int {
	return mat.Columns()
}

// NewMatrix creates a ready to use Matrix of size rows x columns
// This is an alternative to using a Matrix shell with AddRow to create matrices.
func NewMatrix(rows, columns int) Matrix {
	mat := make(Matrix, rows)
	for i := 0; i < rows; i++ {
		mat[i] = make(Vector, columns)
	}
	return mat
}

// NewMatrixShell creates an empty Matrix slice that can be used with AddRow to compose a Matrix
func NewMatrixShell() Matrix {
	return make(Matrix, 0)        // since we are using append to add to the slice
}

// AddRow checks the dimensionality of the argument and adds it if consistent, or returns old set and error
func (mat Matrix) AddRow(vec Vector) (Matrix, error) {
	if len(mat) != 0 && len(vec) != len(mat[0]) {
		return mat, errors.New("Dimensionality of the new vector is not consistent with the existing matrix")
	}
	mat = append(mat, vec)
	return mat, nil
}

// GetRow returns the row vector at index 'index'
func (mat Matrix) GetRow(index int) (Vector, error) {
	if index >= len(mat) {
		return nil, errors.New("Index out of bound")
	}
	return mat[index], nil
}


/*
String represents the Matrix in row-order form
 */
func (mat Matrix) String() string {
	if len(mat) == 0 {
		return "Empty Matrix"
	}
	var str string
	dim := mat.Columns()
	for _, v := range mat {
		for i := 0; i < dim; i++ {
			str = str + fmt.Sprintf("%6.2f ", v[i])
		}
		str = str + "\n"
	}
	return str
}
