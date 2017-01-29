/*
 File		:	$File: //depot/stillwater-sc/godfa/indexspace/vector.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	29 January 2017

 Source Control Information:
 Version	:	$Revision: #1 $
 Latest		:	$Date: 2017/01/29 $
 Location	:	$Id: //depot/stillwater-sc/godfa/indexspace/vector.go#1 $

 Organization:
		Stillwater Supercomputing, Inc.
		P.O Box 720
		South Freeport, ME 04078-0720

Copyright (c) 2006-2017 E. Theodore L. Omtzigt.  All rights reserved.

Licence      : Stillwater license as defined in this directory

 */
package indexspace

import "fmt"

// IntVector represents an integer vector
type IntVector []int

// Dimensionality returns the length of the vector
func (v IntVector) Dimensionality() int {
	return len(v)
}

// Vector represents a vector
type Vector []float64

// Clone creates a deep copy of the receiver matrix
func (v Vector) Clone() Vector {
	n := len(v)
	newVector := make([]float64, n)
	for i := 0; i < n; i++ {
		newVector[i] = v[i]
	}
	return newVector
}

// Dimensionality returns the length of the vector
func (v Vector) Dimensionality() int {
	return len(v)
}

// NewIntVector creates a new slice of size 'size'
func NewIntVector(size int) IntVector {
	return make(IntVector, size)
}

// NewVector creates a new slice of size 'size'
func NewVector(size int) Vector {
	return make(Vector, size)
}

// Stringer for a IntVector
func (v IntVector) String() string {
	var str string
	for i := 0; i < len(v); i++ {
		str = str + fmt.Sprintf("%3d ", v[i])
	}
	return str
}

// Stringer for a Vector
func (v Vector) String() string {
	var str string
	for i := 0; i < len(v); i++ {
		str = str + fmt.Sprintf("%6.2f ", v[i])
	}
	return str
}