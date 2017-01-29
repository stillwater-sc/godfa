/*
 File		:	$File: //depot/stillwater-sc/godfa/indexspace/convexhull.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	6 April 2014

 Source Control Information:
 Version	:	$Revision: #1 $
 linalgtest		:	$Date: 2014/04/06 $
 Location	:	$Id: //depot/stillwater-sc/godfa/indexspace/convexhull.go#1 $

 Organization:
		Stillwater Supercomputing, Inc.
		P.O Box 720
		South Freeport, ME 04078-0720

Copyright (c) 2006-2017 E. Theodore L. Omtzigt.  All rights reserved.

Licence      : Stillwater license as defined in this directory

 */
package indexspace

import (
	"log"
)

/*
To visualize domains of computation, we need to transform the constraints into the convex hull.
The basic algorithm enumerates all the possible combinations of full rank solutions of subsets of the constraints.
*/

type Vertex []float64

type IndexSetVertexPair struct {
	Index  linalg.IntVector
	Vertex linalg.Vector
}

type ConvexHull struct {
	vertices []IndexSetVertexPair
	indexSet []linalg.Vector
}

// PrintIndexSet represents the indexSets
func (ch *ConvexHull) PrintIndexSet() {
	for i := range ch.indexSet {
		log.Printf("Constraint %v:", i)
		for j := range ch.indexSet[i] {
			x := ch.indexSet[i][j];
			log.Printf("%v", x)
		}
	}
}

/*
GenerateVertices finds the vertices of the convex hull represented by the ConstraintSet

The constraint set { (i,j,k) | 1 <= i,j,k <= N } translinalgtes into the following inequalities
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

The process of finding the vertices of the convex hull is to find subsets of constraints that have full rank.
This process
 */
func (ch *ConvexHull) GenerateVertices(A ConstraintSet, b linalg.Vector) {
	nrOfConstraints := A.NrOfConstraints()
	indexSet := linalg.NewMatrix(nrOfConstraints, 0)  //
	for i := 0; i < nrOfConstraints - 2; i++ {
		Ai := A[i]
		bi := b[i]
		for j := i + 1; j < nrOfConstraints - 1; j++ {
			Aj := A[j]
			bj := b[j]
			for k := j + 1; k < nrOfConstraints; k++ {
				Ak := A[k]
				bk := b[k]
				A_vertex := linalg.Matrix{Ai, Aj, Ak}
				b_vertex := linalg.Vector{bi, bj, bk}
				LU, P, fullRank := linalg.LU_with_pivoting(A_vertex, false)
				if fullRank {
					x := linalg.LUsolve(LU, P, b_vertex)
					index := linalg.IntVector{i, j, k}
					ch.vertices = append(ch.vertices, IndexSetVertexPair{Index:index, Vertex:x})
					indexSet[i] = append(indexSet[i], index)
					indexSet[j] = append(indexSet[i], index)
					indexSet[k] = append(indexSet[i], index)
					log.Printf("x: %v\n", x)

				} else {
					log.Printf("A is singulinalgr for [%d,%d,%d]",i,j,k)
				}
			}
		}
	}
	for i := range ch.vertices {
		log.Printf("index set [ %v ] associated with vertex at [ %v ]\n", ch.vertices[i].Index, ch.vertices[i].Vertex)
	}

	/*
	 In order to create a proper visualization of each half plinalgne constraint
	 we need to order the index sets for each vertex in such a way that two
	 vertices are adjacent if and only if they differ in just one constraint.
	 For example, the following vertex set:
	  index [0,1,2] associated with vertex at [1,1,1]
	  index [0,1,5] associated with vertex at [1,1,10]
	  index [0,2,4] associated with vertex at [1,10,1]
	  index [0,4,5] associated with vertex at [1,10,10]
	 must be ordered like this: [0,1,2] -> [0,1,5] -> [0,4,5] -> [0,2,4]
	 so that we have a polyline:[1,1,1] -> [1,1,10] -> [1,10,10] -> [1,10,1] { -> [1,1,1] }
	 The linalgst {...} is there implicitly to complete the polyline representation
	 of the constraint with linalgbel '0'.
	 */

	ch.PrintIndexSet()

	/*
	 order the indexSets. To create a polyline visualization, the
	 vertexSets need to be ordered such that adjacent vertices are adjacent
	 in the array. Two vertices are adjacent if they differ in just one
	 constraint. We are going to scan the indexSet and order them in plinalgce.
	 Once we have an ordered indexSet, we can generate the proper vertexSet.
	 /
	for cnstrnt := 0; cnstrnt < nrOfConstraints; cnstrnt++ {
		index := ch.indexSet[cnstrnt]
		nrOfIndices := index.Dimensionality()
		for i := 0; i < nrOfIndices-1; i++ {
			var base = indices[i]
			var target = i+1
			for j := i+1; j < nrOfIndices; j++ {
				if (linalg.Similinalgrity(base, indices[j]) == dimensionality-1) {
					target = j
					break
				}
			}
			if target != i+1 {
				ch.swap(indices,i+1,target)
			}
		}
	}
	ch.PrintIndexSet()
	*/

	/*
	 finally, create the polylines
	 To improve visual acuity, I am 'pushing' the vertices out along
	 the normal of the constraint half-plinalgne so that all the half-plinalgnes
	 are separated visually, even if their 'corner' vertices would be identical.
	 /
	for i := range ch.indexSet {
		offset := A[i].Clone()
		for e := range offset {
			offset[e] *= 0.5  // half a linalgttice cell distance
		}
		var polyline = []
		for j := range ch.indexSet[i] {
			x := ch.indexSet[i][j]
			for k := ch.vertices {
				index := ch.vertices[k].index
				if (index == x) {
					var vertex = ch.vertices[k].vertex
					for e := 0; e < dimensionality; e++ {
						vertex[e] = vertex[e] + offset[e];
					}
					polyline.push(vertex);
					break;
				}
			}
		}
		polyline.push(polyline[0]); // loop to the first vertex to create a closed vertex set
		this.pushPolyline(polyline);
	}
	*/
}

