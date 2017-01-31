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
	"sort"
)

/*
To visualize domains of computation, we need to transform the constraints into the convex hull.
The basic algorithm enumerates all the possible combinations of full rank solutions of subsets of the constraints.
*/

type Vertex []float64

type ConvexHull struct {
	/*
	 The ordered set of vertices of the convex hull
	 */
	vertices []Vector
	/*
	 for each vertex, we store the indices of all the half-plane constraints that create the vertex

	 If we have this constraint set:
	 The C^T matrix is represented by:
	  index    0  1  2  3  4  5
	 (i,j,k)[  1  0  0 -1  0  0 ]    [  N  N  N -1 -1 -1 ]
		[  0  1  0  0 -1  0 ] <=
		[  0  0  1  0  0 -1 ]
	then, for the vertex (N,N,N), which is defined by constraints at index 0, 1, and 2,
	we store these indices, that is 0,1,2 at the index for that vertex.
	If we assume that this vertex ends up on index 0, we store at indexSet[0] the array [0, 1, 2]
	  */
	indexSet          [][]int

	/*
	 Helper to keep the constraint set polygon sets ordered
	 */
	sortedHalfplanes  []int

	/*
	 Polygon vertex sets are used to organize the vertices that define a half-plane polygon
	 Each polygon is a slice of indexSets. Each indexSet identifies the indices of the constraints
	 that define the vertex. And a polygon is a sequence of adjacent vertices in the same half-plane.
	 */
	polygonIndexSet   map[int][][]int

	/*
	 The visual presentation of the polygons
	 */
	polylineVertexSet map[int][]Vector
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

/*
GenerateVertices finds the vertices of the convex hull represented by a ConstraintSet

The constraint set { (i,j,k) | 1 <= i,j,k <= N } translates into the following inequalities
  1 <= i -> i => 1 -> -i <= -1
 -i <= -1
 -j <= -1
 -k <= -1
  i <= N
  j <= N
  k <= N

The constraint set is represented by its transposed form index^T * C^T <= b^T:
 (i,j,k)[ -1  0  0  1  0  0 ] <= [ -1 -1 -1  N  N  N ]
        [  0 -1  0  0  1  0 ]
        [  0  0 -1  0  0  1 ]

The process of finding the vertices of the convex hull is to find subsets of constraints that have full rank.

Each constraint, such as i <= N or 1 <= j, represents a half space.
A full rank solution to a set of constraints, represents a vertex of the convex hull.
Thus the process of generating the vertices of a (closed) convex hull consists of enumerating all the combinations
of half-space constraints, and labeling these solutions as vertices.

We would also like to visualize the convex hull, and for that we need to organize the vertices in sets that
trace out the half-space polygon that represents the edge of the convex hull.

OpenGL has polylines, which can be used to draw such closed polygons, and the input requires that we use
a list/sequence of adjacent vertices. Two vertices are adjacent when they share dim-1 half-space constraints,
where dim is the dimensionality of the vertices. For example, a vertex in 3D Euclidean space is defined by
the full rank solution of three half-space constraints. Two vertices in 3D Euclidean space are adjacent,
when they share two half-space constraints.
 */
func (ch *ConvexHull) GenerateVertices(A ConstraintSet, b Vector) {

	nrOfHalfspaceConstraints := A.NrOfConstraints()
//	ch.indexSet = make([][]int, 0)		// we'll need to use append to push indexSets as they become known
	ch.polygonIndexSet = make(map[int][][]int) // for each nrOfHalfspaceConstraints we'll store the array of indexSets
	vertexIndex := 0
	for i := 0; i < nrOfHalfspaceConstraints - 2; i++ {
		Ai := TransformConstraintIntoVector(A[i])
		bi := b[i]
		for j := i + 1; j < nrOfHalfspaceConstraints - 1; j++ {
			Aj := TransformConstraintIntoVector(A[j])
			bj := b[j]
			for k := j + 1; k < nrOfHalfspaceConstraints; k++ {
				Ak := TransformConstraintIntoVector(A[k])
				bk := b[k]
				A_vertex := Matrix{Ai, Aj, Ak}
				b_vertex := Vector{bi, bj, bk}
				LU, P, fullRank := LU_with_pivoting(A_vertex, false)
				if fullRank {
					index := []int{i,j,k}
					x := LUsolve(LU, P, b_vertex)
					ch.vertices = append(ch.vertices, x)
					ch.indexSet = append(ch.indexSet, index)
					//log.Printf("x[%d]: %v", vertexIndex, x)
					vertexIndex = vertexIndex + 1
					// mark off the vertices that belong to each half plane
					ch.polygonIndexSet[i] = append(ch.polygonIndexSet[i], index)
					ch.polygonIndexSet[j] = append(ch.polygonIndexSet[j], index)
					ch.polygonIndexSet[k] = append(ch.polygonIndexSet[k], index)
				} else {
					//log.Printf("A is singular for [%d,%d,%d]",i,j,k)
				}
			}
		}
	}
	for index := range ch.polygonIndexSet {
		ch.sortedHalfplanes = append(ch.sortedHalfplanes, index)
	}
	sort.Ints(ch.sortedHalfplanes)

	for i := range ch.vertices {
		log.Printf("index set [ %v ] associated with vertex [ %v ]\n", ch.indexSet[i], ch.vertices[i])
	}
	ch.PrintHalfplanePolygonSets()


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

	 order the indexSets. To create a polyline visualization, the
	 vertexSets need to be ordered such that adjacent vertices are adjacent
	 in the polyline array. Two vertices are adjacent if they differ in just one
	 constraint. We are going to scan the indexSet and order them in place.
	 Once we have an ordered indexSet, we can generate the proper vertexSet.

	 As indicated before, each constraint represents a half-space. By construction,
	 we know that the visualization polygons for the convex hull will lie in the
	 half-plane dividing that half-space. Thus we need to create a polygon definition
	 for each constraint (-> half-plane containing the polygon).
	*/

	dim := A[0].Dimensionality()
	for halfplaneIndex := 0; halfplaneIndex < nrOfHalfspaceConstraints; halfplaneIndex++ {
		polygonIndexSet := ch.polygonIndexSet[halfplaneIndex]

		// it doesn't really matter what the starting vertex of the polygon is
		// so we simply take the first indexSet defining the vertex as base
		// Now use a bubble sort to organize the list
		nrOfVerticesInPolygon := len(polygonIndexSet)
		for i := 0; i < nrOfVerticesInPolygon-1; i++ {
			base := polygonIndexSet[i]
			target := i+1
			j := i+1
			for j = i+1; j < nrOfVerticesInPolygon; j++ {
				// log.Printf("base: %v  cmp: %v", base, polygonIndexSet[j])
				if ch.similarity(base, polygonIndexSet[j]) == dim - 1 {
					target = j
					//log.Printf("Similar to within one: %d", nrOfIndices)
					break
				}
			}
			if target != i+1 {
				ch.swap(polygonIndexSet, i+1,target)
			}
		}
	}
	ch.PrintHalfplanePolygonSets()

	/*
	 finally, create the polylines that we can send to an OpenGL renderer

	 To improve visual acuity, I am 'pushing' the vertices out along
	 the normal of the constraint half-plane so that all the half-planes
	 are separated visually, even if their 'corner' vertices would be identical.
	 */
	ch.polylineVertexSet = make(map[int][]Vector)
	for i, halfplaneIndex := range ch.sortedHalfplanes {
		// take the constraint vector, which is the normal of the halfplane
		// and use it to 'push' the visual representation out by a little bit
		offset := A[i].Clone().CreateVector()
		for e := range offset {
			offset[e] = 0.5 * offset[e]  // half a lattice cell distance
		}
		
		for _, vertexSet := range ch.polygonIndexSet[halfplaneIndex] {
			// find the index where the vertex is stored
			for k := range ch.indexSet {
				if ch.VertexIndexSetMatch(vertexSet, ch.indexSet[k]) {
					// found the index at 'k'
					var vertex = ch.vertices[k]	// copy the vertex
					for e := 0; e < dim; e++ {
						vertex[e] = vertex[e] + offset[e]; // push it out
					}
					ch.polylineVertexSet[halfplaneIndex] = append(ch.polylineVertexSet[halfplaneIndex], vertex) // store it
					break;
				}
			}
		}
		ch.polylineVertexSet[halfplaneIndex] = append(ch.polylineVertexSet[halfplaneIndex], ch.polylineVertexSet[halfplaneIndex][0]) // loop to the first vertex to create a closed vertex set
	}
	ch.PrintHalfplanePolylineSets()
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
func (ch *ConvexHull) similarity(i1, i2 []int) int {
	// count the number of similarities
	var similarity int = 0
	if len(i1) == len(i2) {
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

func (ch *ConvexHull) swap(indexSet [][]int, i1, i2 int) [][]int {
	log.Printf("Swapping index %d with %d", i1, i2)
	tmp := indexSet[i1]
	indexSet[i1] = indexSet[i2]
	indexSet[i2] = tmp
	return indexSet
}

// PrintHalfplanePolygonSets represents the Polygon indexsets into text
func (ch *ConvexHull) PrintHalfplanePolygonSets() {
	for _, halfplaneIndex := range ch.sortedHalfplanes {
		log.Printf("Half-plane %d:\n", halfplaneIndex)
		for i, indexSet := range ch.polygonIndexSet[halfplaneIndex] {
			log.Printf("\t[%3d] {%v}", i, indexSet)
		}
	}
}

// PrintHalfplanePolylineSets represents the Polylines into text
func (ch *ConvexHull) PrintHalfplanePolylineSets() {
	for _, halfplaneIndex := range ch.sortedHalfplanes {
		log.Printf("Half-plane %d:\n", halfplaneIndex)
		for i, vertex := range ch.polylineVertexSet[halfplaneIndex] {
			log.Printf("\t[%3d] {%v}", i, vertex)
		}
	}
}

func (ch *ConvexHull) VertexIndexSetMatch(i1, i2 []int) bool {
	dim := len(i1)
	for i := 0; i < dim; i++ {
		if i1[i] != i2[i] {
			return false
		}
	}
	return true
}
