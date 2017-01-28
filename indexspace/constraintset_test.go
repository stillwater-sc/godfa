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

func TestNewConstraintSet(t *testing.T) {
	c1 := NewConstraint(3)
	c1[0] = 1
	cs := NewConstraintSet()
	cs, err := cs.AddConstraint(c1)
	if err != nil {
		t.Error(err)
	}
	c2 := NewConstraint(3)
	c2[1] = 1
	cs, err = cs.AddConstraint(c2)
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
	fmt.Print(cs)
}