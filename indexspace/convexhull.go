/*
 File		:	$File: //depot/stillwater-sc/godfa/indexspace/convexhull.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	6 April 2014

 Source Control Information:
 Version	:	$Revision: #1 $
 latest		:	$Date: 2014/04/06 $
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
	Index  IntVector
	Vertex Vector
}

type ConvexHull struct {
	vertices []IndexSetVertexPair
	indexSet [][]IntVector
}

/*
TransformConstraintIntoVector transforms a Constraint into a Vector
 */
func TransformConstraintIntoVector(c Constraint) Vector {
	vec := make(Vector, c.Dimensionality())
	for i := range c {
		vec[i] = float64(c[i])
	}
	return vec
}

func NewConvexHull() *ConvexHull {
	return &ConvexHull{}
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
 Similarity calculates the similarity between two index vectors

 This function returns the number of indices that are shared among the two input vectors.
   @param i1: vector containing a set of indices
   @param i2: second vector to compare to
 @return: number of element matches between the two vectors
 Example: i1 = [ 0, 1, 2], i2 = [ 1, 2, 3], similarity will return 2
 for the sum of the matches of the elements, [1, 2].
 Otherwise stated, similarity calculates the cardinality of the intersection of the two input sets.

 I am using the O(n^2) brute force method of enumerating both sets.
 We can get away with that as the index sets are relatively small, typically
 less than 3 element triples and sets less than 5 or 6 elements.
  */
func (ch *ConvexHull) similarity(i1, i2 IntVector) int {
	// count the number of similarities
	var similarity int = 0
	if i1.Dimensionality() == i2.Dimensionality() {
		for i := range i1 {
			for j := range i2 {
				if i1[i] == i2[j] {
					similarity++
					break;
				}
			}
		}
	}
	return similarity
}

func (ch *ConvexHull) swap(indexSet []IntVector, i1, i2 int) {
	tmp := indexSet[i1]
	indexSet[i1] = indexSet[i2]
	indexSet[i2] = tmp
}
/*
GenerateVertices finds the vertices of the convex hull represented by the ConstraintSet

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

The process of finding the vertices of the convex hull is to find subsets of constraints that have full rank.
This process
 */
func (ch *ConvexHull) GenerateVertices(A ConstraintSet, b Vector) {
	nrOfConstraints := A.NrOfConstraints()
	ch.indexSet = make([][]IntVector, nrOfConstraints)
	for i := 0; i < nrOfConstraints; i++ {
		ch.indexSet[i] = make([]IntVector,0) // we'll use append to add elements to the vector
	}
	for i := 0; i < nrOfConstraints - 2; i++ {
		Ai := TransformConstraintIntoVector(A[i])
		bi := b[i]
		for j := i + 1; j < nrOfConstraints - 1; j++ {
			Aj := TransformConstraintIntoVector(A[j])
			bj := b[j]
			for k := j + 1; k < nrOfConstraints; k++ {
				Ak := TransformConstraintIntoVector(A[k])
				bk := b[k]
				A_vertex := Matrix{Ai, Aj, Ak}
				b_vertex := Vector{bi, bj, bk}
				LU, P, fullRank := LU_with_pivoting(A_vertex, false)
				if fullRank {
					x := LUsolve(LU, P, b_vertex)
					index := IntVector{i, j, k}
					ch.vertices = append(ch.vertices, IndexSetVertexPair{Index:index, Vertex:x})
					ch.indexSet[i] = append(ch.indexSet[i], index)
					ch.indexSet[j] = append(ch.indexSet[j], index)
					ch.indexSet[k] = append(ch.indexSet[k], index)
					log.Printf("x: %v", x)

				} else {
					log.Printf("A is singular for [%d,%d,%d]",i,j,k)
				}
			}
		}
	}
	for i := range ch.vertices {
		log.Printf("index set [ %v ] associated with vertex at [ %v ]\n", ch.vertices[i].Index, ch.vertices[i].Vertex)
	}

	/*
	 In order to create a proper visualization of each half plane constraint
	 we need to order the index sets for each vertex in such a way that two
	 vertices are adjacent if and only if they differ in just one constraint.
	 For example, the following vertex set:
	  index [0,1,2] associated with vertex at [1,1,1]
	  index [0,1,5] associated with vertex at [1,1,10]
	  index [0,2,4] associated with vertex at [1,10,1]
	  index [0,4,5] associated with vertex at [1,10,10]
	 must be ordered like this: [0,1,2] -> [0,1,5] -> [0,4,5] -> [0,2,4]
	 so that we have a polyline:[1,1,1] -> [1,1,10] -> [1,10,10] -> [1,10,1] { -> [1,1,1] }
	 The last {...} is there implicitly to complete the polyline representation
	 of the constraint with label '0'.
	 */

	ch.PrintIndexSet()

	/*
	 order the indexSets. To create a polyline visualization, the
	 vertexSets need to be ordered such that adjacent vertices are adjacent
	 in the polyline array. Two vertices are adjacent if they differ in just one
	 constraint. We are going to scan the indexSet and order them in place.
	 Once we have an ordered indexSet, we can generate the proper vertexSet.
	 */
	for cnstrnt := 0; cnstrnt < nrOfConstraints; cnstrnt++ {
		indexSet := ch.indexSet[cnstrnt]
		nrOfIndices := indexSet[0].Dimensionality()
		for i := 0; i < nrOfIndices-1; i++ {
			var base = indexSet[i]
			var target = i+1
			for j := i+1; j < nrOfIndices; j++ {
				log.Printf("base: %s  cmp: %s", base, indexSet[j])
				if ch.similarity(base, indexSet[j]) == nrOfIndices-1 {
					target = j
					log.Printf("Similar to within one: %d", nrOfIndices)
					break
				}
			}
			if target != i+1 {
				ch.swap(indexSet, i+1,target)
			}
		}
	}
	ch.PrintIndexSet()


	/*
	 finally, create the polylines
	 To improve visual acuity, I am 'pushing' the vertices out along
	 the normal of the constraint half-plane so that all the half-planes
	 are separated visually, even if their 'corner' vertices would be identical.
	 /
	for i := range ch.indexSet {
		offset := A[i].Clone()
		for e := range offset {
			offset[e] *= 0.5  // half a lattice cell distance
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

