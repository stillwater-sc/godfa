/*
 File		:	$File: //depot/stillwater-sc/godfa/indexspace/constraintset.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	6 April 2014

 Source Control Information:
 Version	:	$Revision: #1 $
 Latest		:	$Date: 2014/04/06 $
 Location	:	$Id: //depot/stillwater-sc/godfa/indexspace/constraintset.go#1 $

 Organization:
		Stillwater Supercomputing, Inc.
		P.O Box 720
		South Freeport, ME 04078-0720

Copyright (c) 2006-2017 E. Theodore L. Omtzigt.  All rights reserved.

Licence      : Stillwater license as defined in this directory

 */
package indexspace

import (
	"errors"
	"fmt"
)

/*
DoC = Domain of Computation.

DoC is the set of constraints defining the domain of computation
of a system of recurrence equations. The domain flow language defines
the domains for the system of recurrences, but individual recurrences
interact with the domain through more specific constraints. It is the
role of the compiler to optimize the DoC for each recurrence equation.

Here is an example of a domain flow algorithm:
  compute( (i,j,k) | 1 <= i,j,k <= N) {
      a(i,j,k) = a(i-1,j,k)                            <- recurrence 1
      b(i,j,k) = b(i,j-1,k)                            <- recurrence 2
      c(i,j,k) = c(i,j,k-1) + a(i-1,j,k) * b(i,j-1,k)  <- recurrence 3
  }

For this example, we have a DoC that is a 3D cube defined by
the constraint set { (i,j,k) | 1 <= i,j,k <= N }
This translates into the following inequalities
  i <= N
  j <= N
  k <= N
  i => 1 => -i <= -1
  -j <= -1
  -k <= -1
The C^T matrix is represented by:
 (i,j,k)[  1  0  0 -1  0  0 ]    [  N ]
        [  0  1  0  0 -1  0 ] <= [  N ]
        [  0  0  1  0  0 -1 ]    [  N ]

*/

// A Constraint is a row vector
type Constraint []int

// A ConstraintSet is a slice of Constraints
type ConstraintSet []Constraint

// NewConstraint creates a new slice of size 'size'
func NewConstraint(size int) Constraint {
	return make(Constraint, size)
}

// Clone a constraint
func (c Constraint) Clone() Constraint {
	newC := NewConstraint(c.Dimensionality())
	for i := range c {
		newC[i] = c[i]
	}
	return newC
}

// NewConstraintSet creates an empty ConstraintSet slice
func NewConstraintSet() ConstraintSet {
	return make(ConstraintSet, 0)        // since we are using append to add to the slice
}

// AddConstraint checks the dimensionality of the argument and adds it if consistent, or returns old set and error
func (cSet ConstraintSet) AddConstraint(constraint Constraint) (ConstraintSet, error) {
	if len(cSet) != 0 && len(constraint) != len(cSet[0]) {
		return cSet, errors.New("Dimensionality of the new constraint is not consistent with the set")
	}
	cSet = append(cSet, constraint.Clone())
	return cSet, nil
}

// GetConstraint returns the Constraint at index 'index'
func (cSet ConstraintSet) GetConstraint(index int) (Constraint, error) {
	if index >= len(cSet) {
		return nil, errors.New("Index out of bound")
	}
	return cSet[index], nil
}

// NrOfConstraints returns the number of Constraints in the Set
func (cSet ConstraintSet) NrOfConstraints() int {
	return len(cSet)
}

// Dimensionality returns the dimension of the Constraint
func (c Constraint) Dimensionality() int {
	return len(c)
}

// Dimensionality returns the dimension of the Constraints in the Set
func (cSet ConstraintSet) Dimensionality() int {
	if len(cSet) == 0 {
		return 0
	}
	return cSet[0].Dimensionality()
}

// Stringer for a Constraint
func (c Constraint) String() string {
	var str string
	for i := 0; i < len(c); i++ {
		str = str + fmt.Sprintf("%3d ", c[i])
	}
	return str
}

/*
String represents the ConstraintSet in its Transposed form

The C^T matrix is represented by:
 (i,j,k)[  1  0  0 -1  0  0 ]    [  N ]
        [  0  1  0  0 -1  0 ] <= [  N ]
        [  0  0  1  0  0 -1 ]    [  N ]
  */
func (cSet ConstraintSet) String() string {
	if len(cSet) == 0 {
		return "Empty Constraint Set"
	}
	var str string
	var dim int = len(cSet[0])
	for i := 0; i < dim; i++ {
		for _, c := range cSet {
			str = str + fmt.Sprintf("%3d ", c[i])
		}
		str = str + "\n"
	}
	return str
}