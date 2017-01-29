/*
 File		:	$File: //depot/stillwater-sc/godfa/linalg/matrix_test.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	28 January 2017

 Source Control Information:
 Version	:	$Revision: #1 $
 Latest		:	$Date: 2017/01/28 $
 Location	:	$Id: //depot/stillwater-sc/godfa/linalg/matrix_test.go#1 $

 Organization:
		Stillwater Supercomputing, Inc.
		P.O Box 720
		South Freeport, ME 04078-0720

Copyright (c) 2006-2017 E. Theodore L. Omtzigt.  All rights reserved.

Licence      : Stillwater license as defined in this directory

 */
package la

import (
	"testing"
	"fmt"
)

func TestMatrix_Creation(t *testing.T) {
	c1 := NewVector(3)
	c1[0] = 1
	cs := NewMatrixShell()
	cs, _ = cs.AddRow(c1)	// on the first it can't really throw any errors

	c2 := NewVector(3)
	c2[1] = 1
	cs, err := cs.AddRow(c2)	// on subsequent calls you may submit a constraint with a different dimensionality
	if err != nil {
		t.Error(err)
	}
	c3 := NewVector(3)
	c3[2] = 1
	cs, err = cs.AddRow(c3)
	if err != nil {
		t.Error(err)
	}

	cc, err := cs.GetRow(0)
	if err != nil {
		t.Error(err)
	}
	if cc[0] != c1[0] {
		t.Error("Constraint didn't persist properly")
	}

}

func TestMatrix_SpecialConditions(t *testing.T) {
	c1 := Vector{1,0,0}
	cs := NewMatrixShell()
	if cs.Dimensionality() != 0 {
		t.Error("Dimensionality of empty set is incorrect")
	}
	cs, _ = cs.AddRow(c1)
	if cs.Dimensionality() != c1.Dimensionality() {
		t.Error("Dimensionality of the Set is incorrect")
	}
	c2 := Vector{1,0,0,0}
	cs, err := cs.AddRow(c2)
	if err == nil {
		t.Error("ConstraintSet.AddRow didn't recognized incompatible constraint dimensions")
	}
	c2, err = cs.GetRow(2)
	if err == nil {
		t.Error("ConstraintSet.GetRow didn't recognize out of bounds index")
	}
}

func TestMatrix_String(t *testing.T) {
	c1 := Vector{1,0,0,-1,0,0}
	c1_str := c1.String()
	c1_ref := "  1.00   0.00   0.00  -1.00   0.00   0.00 "   // notice the trailing space
	if c1_str != c1_ref {
		t.Error("Constraint stringer failed")
		t.Logf("%x\n", c1_str)
		t.Logf("%x\n", c1_ref)
	}
	c2 := Vector{0,1,0,0,-1,0}
	c3 := Vector{0,0,1,0,0,-1}
	cs := NewMatrixShell()
	cs_str := cs.String()		// empty constraint set
	cs_ref := "Empty Matrix"
	if cs_str != cs_ref {
		t.Error("ConstraintSet stringer didn't recognize an empty matrix")
		t.Logf("%x\n", cs_str)
		t.Logf("%x\n", cs_ref)
	}
	cs, _ = cs.AddRow(c1)
	cs, _ = cs.AddRow(c2)
	cs, _ = cs.AddRow(c3)
	cs_str = cs.String()
	cs_ref = "  1.00   0.00   0.00  -1.00   0.00   0.00 \n  0.00   1.00   0.00   0.00  -1.00   0.00 \n  0.00   0.00   1.00   0.00   0.00  -1.00 \n"
	if cs_str != cs_ref {
		t.Error("ConstraintSet stringer failed")
		t.Logf("%x\n", cs_str)
		t.Logf("%x\n", cs_ref)
	}
}

func ExampleConstraintSet() {
	c1 := NewVector(3)
	c1[0] = 1
	cs := NewMatrixShell()
	cs, _ = cs.AddRow(c1)	// on the first it can't really throw any errors
	c1 = Vector{0,1,0}
	cs, _ = cs.AddRow(c1)
	c2 := Vector{0,0,1}
	cs, _ = cs.AddRow(c2)
	fmt.Print(cs)

	// OUTPUT
	//   1   0   0
	//   0   1   0
	//   0   0   1
}
