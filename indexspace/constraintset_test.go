/*
 File		:	$File: //depot/stillwater-sc/godfa/indexspace/constraintset_test.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	28 January 2017

 Source Control Information:
 Version	:	$Revision: #1 $
 Latest		:	$Date: 2017/01/28 $
 Location	:	$Id: //depot/stillwater-sc/godfa/indexspace/constraintset_test.go#1 $

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
	"fmt"
)

func TestConstraintSet(t *testing.T) {
	c1 := NewConstraint(3)
	c1[0] = 1
	cs := NewConstraintSet()
	cs, _ = cs.AddConstraint(c1)	// on the first it can't really throw any errors

	c2 := NewConstraint(3)
	c2[1] = 1
	cs, err := cs.AddConstraint(c2)	// on subsequent calls you may submit a constraint with a different dimensionality
	if err != nil {
		t.Error(err)
	}
	c3 := NewConstraint(3)
	c3[2] = 1
	cs, err = cs.AddConstraint(c3)
	if err != nil {
		t.Error(err)
	}

	cc, err := cs.GetConstraint(0)
	if err != nil {
		t.Error(err)
	}
	if cc[0] != c1[0] {
		t.Error("Constraint didn't persist properly")
	}

}

func TestConstraintSet_SpecialConditions(t *testing.T) {
	c1 := Constraint{1,0,0}
	cs := NewConstraintSet()
	cs, _ = cs.AddConstraint(c1)
	c2 := Constraint{1,0,0,0}
	cs, err := cs.AddConstraint(c2)
	if err == nil {
		t.Error("ConstraintSet.AddConstraint didn't recognized incompatible constraint dimensions")
	}
	c2, err = cs.GetConstraint(2)
	if err == nil {
		t.Error("ConstraintSet.GetConstraint didn't recognize out of bounds index")
	}
}

func TestConstraintSet_String(t *testing.T) {
	c1 := Constraint{1,0,0,-1,0,0}
	c1_str := c1.String()
	c1_ref := "  1   0   0  -1   0   0 "   // notice the trailing space
	if c1_str != c1_ref {
		t.Error("Constraint stringer failed")
	}
	c2 := Constraint{0,1,0,0,-1,0}
	c3 := Constraint{0,0,1,0,0,-1}
	cs := NewConstraintSet()
	cs_str := cs.String()		// empty constraint set
	cs_ref := "Empty Constraint Set"
	if cs_str != cs_ref {
		t.Error("ConstraintSet stringer didn't recognize an empty set")
	}
	cs, _ = cs.AddConstraint(c1)
	cs, _ = cs.AddConstraint(c2)
	cs, _ = cs.AddConstraint(c3)
	cs_str = cs.String()
	cs_ref = "  1   0   0 \n  0   1   0 \n  0   0   1 \n -1   0   0 \n  0  -1   0 \n  0   0  -1 \n"
	if cs_str != cs_ref {
		t.Error("ConstraintSet stringer failed")
		fmt.Printf("%x\n", cs_str)
		fmt.Printf("%x\n", cs_ref)
	}
}

func ExampleConstraintSet() {
	c1 := NewConstraint(3)
	c1[0] = 1
	cs := NewConstraintSet()
	cs, _ = cs.AddConstraint(c1)	// on the first it can't really throw any errors
	c1 = Constraint{0,1,0}
	cs, _ = cs.AddConstraint(c1)
	c2 := Constraint{0,0,1}
	cs, _ = cs.AddConstraint(c2)
	fmt.Print(cs)

	// OUTPUT
	//   1   0   0
	//   0   1   0
	//   0   0   1
}