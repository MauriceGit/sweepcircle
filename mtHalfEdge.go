package mtSweepCircle

import (
	"fmt"
)

// Just to make the code more readable.
//type *HEEdge int
//type *HEFace int
//type *HEVertex int

var EmptyFace *HEFace = nil
var EmptyEdge *HEEdge = nil
var EmptyVertex *HEVertex = nil
var InfiniteVertex *HEVertex = nil

/**
 * All edges, vertices and faces are kept in static, not-resizable slices
 * and therefore tight data representation.
 *
 * References to vertices, edges or faces are kept as indices into the corresponding slices.
 *
 * It is guaranteed, that the initial lists are long enough for the to-be constructed voronoi.
 * Data does not need to be relocated for normal voronoi creation.
 */

type HEVertex struct {
	Pos Vector
}

func (v *HEVertex) Valid() bool {
	return v != EmptyVertex && v != InfiniteVertex
}

type HEFace struct {
	// Especially for Voronoi
	ReferencePoint Vector

	// Points to an arbitrary edge of its polygon
	// Only arbitrary for closed faces!!!!
	// Otherwise it HAS to point to the FIRST edge (counter clockwise)!!!
	EEdge         *HEEdge
	EEdgeParallel *HEEdge
}

type HEEdge struct {
	VOrigin *HEVertex
	ETwin   *HEEdge
	// Starts from: this->ETwin->VOrigin. With this->FFace == this->ENext->FFace
	// Following ENext will traverse the polygon around FFace!
	ENext *HEEdge
	EPrev *HEEdge
	FFace *HEFace

	// Iff VOrigin is NOT defined, an edge representation is required
	// (general direction + some kind of starting point that is NOT an official HEVertex!)
	TmpEdge Edge
}

// Todo but not important right now.
func (f *HEFace) CalcNormal() Vector {
	return Vector{}
}

func (v HEVertex) String() string {
	return fmt.Sprintf("(%.2f, %.2f)", v.Pos.X, v.Pos.Y)
}

func (f HEFace) String() string {
	return fmt.Sprintf("(%3d)", f.EEdge)
}

func (e HEEdge) String() string {
	return fmt.Sprintf("(o: %3d, t: %3d, n: %3d, p: %3d, f: %3d)", e.VOrigin, e.ETwin, e.ENext, e.EPrev, e.FFace)
}
