/*
 File		:	$File: //depot/stillwater-sc/godfa/indexspace/convexhull_test.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	6 April 2014

 Source Control Information:
 Version	:	$Revision: #1 $
 Latest		:	$Date: 2014/04/06 $
 Location	:	$Id: //depot/stillwater-sc/godfa/indexspace/convexhull_test.go#1 $

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
	"time"
)

/*
The constraint set { (i,j,k) | 1 <= i,j,k <= N } translates into the following inequalities
  i <= N
  j <= N
  k <= N
  i => 1 => -i <= -1
  -j <= -1
  -k <= -1
The C^T matrix is represented by:
 (i,j,k)[  1  0  0 -1  0  0 ]    [  N  N  N -1 -1 -1 ]
        [  0  1  0  0 -1  0 ] <=
        [  0  0  1  0  0 -1 ]
 */
func TestConvexHull_GenerateVertices(t *testing.T) {
	cs := NewConstraintSet()
	c := NewConstraint(3)
	c[0] = 1; c[1] = 0; c[2] = 0
	cs, _ = cs.AddConstraint(c)
	c[0] = 0; c[1] = 1; c[2] = 0
	cs, _ = cs.AddConstraint(c)
	c[0] = 0; c[1] = 0; c[2] = 1
	cs, _ = cs.AddConstraint(c)
	c[0] = -1; c[1] = 0; c[2] = 0
	cs, _ = cs.AddConstraint(c)
	c[0] = 0; c[1] = -1; c[2] = 0
	cs, _ = cs.AddConstraint(c)
	c[0] = 0; c[1] = 0; c[2] = -1
	cs, _ = cs.AddConstraint(c)

	N := 5.0
	b := NewVector(cs.NrOfConstraints())
	b[0] = N; b[1] = N; b[2] = N; b[3] = -1; b[4] = -1; b[5] = -1
	fmt.Printf("Nr of Constraints: %d\n", cs.NrOfConstraints())
	fmt.Printf("index * A^T:\n%s\n", cs.String())
	fmt.Printf("b: %s\n", b.String())

	time.Sleep(time.Duration(1 * time.Second))
	ch := NewConvexHull()
	ch.GenerateVertices(cs, b)
}


