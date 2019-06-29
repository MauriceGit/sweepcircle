package mtSweepCircle

import (
	"fmt"
	//"os"
	"errors"
	"math"

	"sort"

	//v "mtVector"

	//s "github.com/MauriceGit/skiplist"

	//"time"
	"math/rand"
	"sync"
)

const (
	EPS float64 = 0.0000001
)

type Delaunay struct {
	Vertices           []HEVertex
	firstFreeVertexPos VertexIndex
	Edges              []HEEdge
	firstFreeEdgePos   EdgeIndex
	Faces              []HEFace
	firstFreeFacePos   FaceIndex
	frontier           *ArrayMap

	edgeMemPool []HEEdge
}

type Voronoi Delaunay
type ConvexHull []Vector

type DelaunayPoint struct {
	Point      Vector
	Distance   float64
	PolarAngle float64
}
type DelaunayPointList struct {
	Points []DelaunayPoint
	Origin Vector
}

type FrontElement struct {
	//// Index into Delaunay.vertices data structure!
	//Index       VertexIndex
	// Index into Delaunay.edges data structure!
	EdgeIndex EdgeIndex
	// Polar angle of the vertex of EdgeI
	PolarAngle float64
	Radius     float64
}

// For optimized printing or drawing only!!!
type SimpleEdge struct {
	// Start vertex of the edge
	V1 Vector
	// End vertex of the edge
	V2 Vector
}

var EmptyF = HEFace{}
var EmptyE = HEEdge{}
var EmptyV = HEVertex{}

// We will need this matrix millions of times for bigger Delaunay. So we can cache the memory
// and just reuse it.
var g_matrixElements = [4][4]float64{{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}}

func (d Delaunay) String() string {
	s := "Vertices: \n"
	for i, v := range d.Vertices {
		if v != EmptyV {
			s += fmt.Sprintf("    %2d: %v\n", i, v)
		}
	}
	s += "\nEdges: \n"
	for i, e := range d.Edges {
		if e != EmptyE {
			s += fmt.Sprintf("    %2d: %v\n", i, e)
		}
	}
	s += "\nFaces: \n"
	for i, f := range d.Faces {
		if f != EmptyF {
			s += fmt.Sprintf("    %2d: %v\n", i, f)
		}
	}
	return s
}

func (dp DelaunayPoint) String() string {
	return fmt.Sprintf("{Point: %v, Distance: %v, PolarAngle: %v}", dp.Point, dp.Distance, dp.PolarAngle)
}
func (dpl DelaunayPointList) String() string {
	return fmt.Sprintf("\nDelaunayPointList:\n    Origin: %v\n    Points: %v\n", dpl.Origin, dpl.Points)
}

// By is the type of a "less" function that defines the ordering of its Planet arguments.
type By func(i, j *DelaunayPoint) bool

// planetSorter joins a By function and a slice of Planets to be sorted.
type pointSorter struct {
	points []DelaunayPoint
	by     func(i, j *DelaunayPoint) bool // Closure used in the Less method.
}

// Sort is a method on the function type, By, that sorts the argument slice according to the function.
func (by By) Sort(points []DelaunayPoint) {
	ps := &pointSorter{
		points: points,
		by:     by, // The Sort method's receiver is the function (closure) that defines the sort order.
	}
	sort.Sort(ps)
}

// Len is part of sort.Interface.
func (s *pointSorter) Len() int {
	return len(s.points)
}

// Swap is part of sort.Interface.
func (s *pointSorter) Swap(i, j int) {
	s.points[i], s.points[j] = s.points[j], s.points[i]
}

// Less is part of sort.Interface. It is implemented by calling the "by" closure in the sorter.
func (s *pointSorter) Less(i, j int) bool {
	return s.by(&s.points[i], &s.points[j])
}

func smaller(v1, v2 *DelaunayPoint) bool {
	// sort by radius
	if math.Abs(v1.Distance-v2.Distance) >= EPS {
		return v1.Distance < v2.Distance
	}
	// if radius is equal, sort by polar angle
	return v1.PolarAngle < v2.PolarAngle
}

// Custom  and short quicksort from Stackoverflow... https://stackoverflow.com/questions/23276417/golang-custom-sort-is-faster-than-native-sort
func qsort(a []DelaunayPoint) []DelaunayPoint {
	if len(a) < 2 {
		return a
	}

	left, right := 0, len(a)-1

	// Pick a pivot
	pivotIndex := rand.Int() % len(a)

	// Move the pivot to the right
	a[pivotIndex], a[right] = a[right], a[pivotIndex]

	// Pile elements smaller than the pivot on the left
	for i := range a {
		if smaller(&a[i], &a[right]) {
			a[i], a[left] = a[left], a[i]
			left++
		}
	}

	// Place the pivot after the last smaller element
	a[left], a[right] = a[right], a[left]

	// Go down the rabbit hole
	qsort(a[:left])
	qsort(a[left+1:])

	return a
}

////////////////////////////////////////////////////////////////////////
//  Pretty Print the  Voronoi Attributes
////////////////////////////////////////////////////////////////////////

func (d *Delaunay) pprint() {

	fmt.Println("Voronoi:")
	fmt.Printf("   Vertices (%v):\n", d.firstFreeVertexPos)

	for i, ve := range d.Vertices {
		if ve != EmptyV {
			fmt.Printf("\t%v:\tPos: %v\n", i, ve.Pos)
		}
	}
	fmt.Printf("\n")

	fmt.Printf("   Edges    (%v):\n", d.firstFreeEdgePos)
	for i, e := range d.Edges {
		if e != EmptyE {
			fmt.Printf("\t%v:\tOrigin: %v,\tTwin: %v,\tPrev: %v,\tNext: %v,\tFace: %v\n", i, e.VOrigin, e.ETwin, e.EPrev, e.ENext, e.FFace)
		}
	}
	fmt.Printf("\n")

	fmt.Printf("   Faces    (%v):\n", d.firstFreeFacePos)
	for i, f := range d.Faces {
		if f != EmptyF {
			fmt.Printf("\t%v:\tRefPoint: %v,\tEdge:%v\n", i, f.ReferencePoint, f.EEdge)
		}
	}
	fmt.Printf("\n")
}

// Verifies that the Delaunay is not corrupted or has miscalculated edges/faces
// or any other unvalid stuff.
func (d *Delaunay) Verify() error {

	for i, e := range d.Edges {
		if e != EmptyE {

			// Every valid edge MUST have a valid Twin edge!
			if e.ETwin == EmptyEdge || d.Edges[e.ETwin] == EmptyE {
				return errors.New(fmt.Sprintf("Edge %v: %v has an invalid twin edge", i, e))
			}

			// Twins must refer to each other!
			if i != int(d.Edges[e.ETwin].ETwin) {
				return errors.New(fmt.Sprintf("Edge %v and his Twin %v don't refer to each other", i, e.ETwin))
			}

			// Check if the origin vertex is valid (if there is one)
			if e.VOrigin == EmptyVertex || d.Vertices[e.VOrigin] == EmptyV {
				return errors.New(fmt.Sprintf("Origin vertex of edge %v: %v is invalid", i, e))
			}

			// Check, if the face is valid (if defined)
			if e.FFace != EmptyFace && d.Faces[e.FFace] == EmptyF {
				return errors.New(fmt.Sprintf("Face of edge %v: %v is invalid", i, e.FFace))
			}

			// Check, if the next edge is valid
			if e.ENext != EmptyEdge && d.Edges[e.ENext] == EmptyE {
				return errors.New(fmt.Sprintf("Next edge for edge %v: %v is invalid", i, e))
			}

			// Check, if the prev edge is valid
			if e.EPrev != EmptyEdge && d.Edges[e.EPrev] == EmptyE {
				return errors.New(fmt.Sprintf("Prev edge for edge %v: %v is invalid", i, e))
			}

			// If this->next is not empty, the next one must have a previous edge.
			if e.ENext != EmptyEdge && d.Edges[e.ENext].EPrev == EmptyEdge {
				return errors.New(fmt.Sprintf("Next edge for edge %v must have a non-empty prev", i))
			}

			// The this edge must correspond to the next->prev edge
			if e.ENext != EmptyEdge && d.Edges[e.ENext].EPrev != EmptyEdge && d.Edges[e.ENext].EPrev != EdgeIndex(i) {
				return errors.New(fmt.Sprintf("The edge %v has %v as next edge, which has %v as his previous. They must reference each other!", i, e.ENext, d.Edges[e.ENext].EPrev))
			}
		}
	}

	for i, f := range d.Faces {
		if f != EmptyF {

			// Check for valid EEdge
			if f.EEdge == EmptyEdge || d.Edges[f.EEdge] == EmptyE {
				return errors.New(fmt.Sprintf("Face %v: %v points to invalid edge", i, f))
			}

			// If the edge is the first one (no/infinite start-vertex), it's all good.
			// Otherwise check, that the edges actually go all the way around the face!
			startEdge := f.EEdge
			edgeCount := 0
			e := d.Edges[startEdge].ENext

			for e != EmptyEdge && e != startEdge {

				if int(d.Edges[e].FFace) != i {
					return errors.New(fmt.Sprintf("Edge %v of the face %v: %v does not point to the right face!", e, i, f))
				}

				e = d.Edges[e].ENext
				edgeCount += 1

				if edgeCount > 3 {
					return errors.New(fmt.Sprintf("Looping around the edges of the face %v: %v are more than 3", i, f))
				}
			}

		}
	}

	return nil
}

func (d *Delaunay) createFace(refPoint Vector, eEdge EdgeIndex) FaceIndex {

	d.Faces[d.firstFreeFacePos].ReferencePoint = refPoint
	d.Faces[d.firstFreeFacePos].EEdge = eEdge

	d.firstFreeFacePos += 1
	return d.firstFreeFacePos - 1
}

func (d *Delaunay) createEdge(vOrigin VertexIndex, eTwin, ePrev, eNext EdgeIndex, fFace FaceIndex, tmpEdge Edge) EdgeIndex {

	d.Edges[d.firstFreeEdgePos].VOrigin = vOrigin
	d.Edges[d.firstFreeEdgePos].ETwin = eTwin
	d.Edges[d.firstFreeEdgePos].ENext = eNext
	d.Edges[d.firstFreeEdgePos].EPrev = ePrev
	d.Edges[d.firstFreeEdgePos].FFace = fFace
	d.Edges[d.firstFreeEdgePos].TmpEdge = tmpEdge

	d.firstFreeEdgePos += 1
	return d.firstFreeEdgePos - 1
}

func (d *Delaunay) createVertex(pos Vector) VertexIndex {

	d.Vertices[d.firstFreeVertexPos].Pos = pos

	d.firstFreeVertexPos += 1
	return d.firstFreeVertexPos - 1
}

func calcPolarAngle(p Vector, origin Vector) float64 {

	if Equal(p, origin) {
		return 0.0
	}

	diff := Sub(p, origin)
	xAxis := Vector{1, 0}

	angle := AngleRad(diff, xAxis)
	if diff.Y < 0 {
		angle = math.Pi + (math.Pi - angle)
	}

	return angle
}

func preparePointList(points []Vector, threadCount int) DelaunayPointList {
	var dPointList DelaunayPointList

	dPointList.Points = make([]DelaunayPoint, len(points), len(points))

	// We cannot calculate the origin in the same loop as the distance, as distance depends on the final origin...
	origin := Vector{}
	for _, p := range points {
		origin.Add(p)
	}
	origin.Div(float64(len(points)))
	dPointList.Origin = origin

	//for i,p := range *points {
	//    dPointList.Points[i] = DelaunayPoint{
	//                                         Point: p,
	//                                         Distance: v.Length(v.Sub(p, dPointList.Origin)),
	//                                         PolarAngle: calcPolarAngle(p, dPointList.Origin),
	//                                        }
	//}

	var wg sync.WaitGroup
	wg.Add(threadCount)
	groupSize := len(dPointList.Points) / threadCount

	for t := 0; t < threadCount; t++ {
		upperCount := t*groupSize + groupSize
		if t == threadCount-1 {
			upperCount = len(dPointList.Points)
		}
		go func(tMin, tMax int) {
			defer wg.Done()
			for i := tMin; i < tMax; i++ {
				p := points[i]
				dPointList.Points[i] = DelaunayPoint{
					Point:      p,
					Distance:   LengthSquared(Sub(p, dPointList.Origin)),
					PolarAngle: calcPolarAngle(p, dPointList.Origin),
				}
			}
		}(t*groupSize, upperCount)
	}

	wg.Wait()

	/*
		distance := func(v1, v2 *DelaunayPoint) bool {
			// sort by radius
			if math.Abs(v1.Distance-v2.Distance) >= EPS {
				return v1.Distance < v2.Distance
			}
			// if radius is equal, sort by polar angle
			return v1.PolarAngle < v2.PolarAngle
		}
		By(distance).Sort(dPointList.Points)
	*/

	dPointList.Points = qsort(dPointList.Points)

	// ToDo:
	// This is very hacky and should be improved first thing, when the algorithms works...
	// There should be a way to determine the first triangle where Origin is in without sorting/recalculating the list twice!!!
	dPointList.Origin = Div(Add(Add(dPointList.Points[0].Point, dPointList.Points[1].Point), dPointList.Points[2].Point), 3.0)
	// This is useless.
	//for i,p := range dPointList.Points {
	//    dPointList.Points[i].Distance   = v.Length(v.Sub(p.Point, dPointList.Origin))
	//    dPointList.Points[i].PolarAngle = calcPolarAngle(p.Point, dPointList.Origin)
	//}

	wg.Add(threadCount)

	for t := 0; t < threadCount; t++ {
		upperCount := t*groupSize + groupSize
		if t == threadCount-1 {
			upperCount = len(dPointList.Points)
		}
		go func(tMin, tMax int) {
			defer wg.Done()
			for i := tMin; i < tMax; i++ {
				p := dPointList.Points[i]
				dPointList.Points[i].Distance = LengthSquared(Sub(p.Point, dPointList.Origin))
				dPointList.Points[i].PolarAngle = calcPolarAngle(p.Point, dPointList.Origin)
			}
		}(t*groupSize, upperCount)
	}

	wg.Wait()

	// Like ... Really.

	//start := time.Now()
	//By(distance).Sort(dPointList.Points)
	//binTime := time.Since(start).Nanoseconds()

	//fmt.Printf("Sort in Seconds: %.8f\n", float64(binTime)/1000000000.0)

	return dPointList
}

func (d *Delaunay) initializeTriangulation(pl *DelaunayPointList) *ArrayMap {

	//var frontier Frontier = make([]FrontElement, 3, len((*pl).Points))

	f0 := FrontElement{}
	f1 := FrontElement{}
	f2 := FrontElement{}

	// The closest three points are (per definition) the triangle surrounding origin.
	if pl.Points[0].PolarAngle < pl.Points[1].PolarAngle {
		// 0 < 1
		d.Vertices[0] = HEVertex{pl.Points[0].Point}
		f0.PolarAngle = pl.Points[0].PolarAngle
		if pl.Points[2].PolarAngle < pl.Points[0].PolarAngle || pl.Points[2].PolarAngle > pl.Points[1].PolarAngle {
			// 2 < 0 < 1
			d.Vertices[1] = HEVertex{pl.Points[1].Point}
			d.Vertices[2] = HEVertex{pl.Points[2].Point}
			f1.PolarAngle = pl.Points[1].PolarAngle
			f2.PolarAngle = pl.Points[2].PolarAngle
		} else {
			// 0 < 2 < 1
			d.Vertices[1] = HEVertex{pl.Points[2].Point}
			d.Vertices[2] = HEVertex{pl.Points[1].Point}
			f1.PolarAngle = pl.Points[2].PolarAngle
			f2.PolarAngle = pl.Points[1].PolarAngle
		}
	} else {
		d.Vertices[0] = HEVertex{pl.Points[2].Point}
		f0.PolarAngle = pl.Points[2].PolarAngle
		// 1 < 0
		if pl.Points[2].PolarAngle < pl.Points[1].PolarAngle || pl.Points[2].PolarAngle > pl.Points[0].PolarAngle {
			// 2 < 1 < 0
			d.Vertices[1] = HEVertex{pl.Points[1].Point}
			d.Vertices[2] = HEVertex{pl.Points[0].Point}
			f1.PolarAngle = pl.Points[1].PolarAngle
			f2.PolarAngle = pl.Points[0].PolarAngle
		} else {
			// 1 < 2 < 0
			d.Vertices[1] = HEVertex{pl.Points[0].Point}
			d.Vertices[2] = HEVertex{pl.Points[1].Point}
			f1.PolarAngle = pl.Points[0].PolarAngle
			f2.PolarAngle = pl.Points[1].PolarAngle
		}
	}
	d.firstFreeVertexPos = 3

	// We can write edge == 0 because e will be the very first edge that will always be at index 0.
	// We write a -1 vector so it is different to Vector{}.
	f := d.createFace(Vector{-1, -1}, 0)

	// Edge 0-->1
	e := d.createEdge(0, EmptyEdge, EmptyEdge, EmptyEdge, f, Edge{})
	// Edge 1-->0
	eT := d.createEdge(1, e, EmptyEdge, EmptyEdge, EmptyFace, Edge{})
	d.Edges[e].ETwin = eT
	d.Faces[f].EEdge = e
	//feT := d.createFace(Vector{}, eT)
	//d.Edges[eT].FFace = feT

	// Edge 1-->2
	e2 := d.createEdge(1, EmptyEdge, e, EmptyEdge, f, Edge{})
	// Edge 2-->1
	e2T := d.createEdge(2, e2, EmptyEdge, EmptyEdge, EmptyFace, Edge{})
	d.Edges[e2].ETwin = e2T
	//fe2T := d.createFace(Vector{}, e2T)
	//d.Edges[e2T].FFace = fe2T

	// Edge 2-->0
	e3 := d.createEdge(2, EmptyEdge, e2, e, f, Edge{})
	// Edge 0-->2
	e3T := d.createEdge(0, e3, EmptyEdge, EmptyEdge, EmptyFace, Edge{})
	d.Edges[e3].ETwin = e3T
	//fe3T := d.createFace(Vector{}, e3T)
	//d.Edges[e3T].FFace = fe3T

	d.Edges[e].EPrev = e3
	d.Edges[e].ENext = e2
	d.Edges[e2].ENext = e3

	f0.EdgeIndex = eT
	f0.PolarAngle = calcPolarAngle(d.Vertices[1].Pos, (*pl).Origin)
	f1.EdgeIndex = e2T
	f1.PolarAngle = calcPolarAngle(d.Vertices[2].Pos, (*pl).Origin)
	f2.EdgeIndex = e3T
	f2.PolarAngle = calcPolarAngle(d.Vertices[0].Pos, (*pl).Origin)

	// Pop first three points, because they are already triangulated by default
	pl.Points = pl.Points[3:]

	//frontier := s.NewEps(0.000000001)
	frontier := NewArrayMap(int(math.Sqrt(float64(len(pl.Points)))))

	//fmt.Printf("First three inserts: %v, %v, %v\n", f0, f1, f2)
	n1 := frontier.InsertAfter(f0, nil)
	n2 := frontier.InsertAfter(f1, n1)
	frontier.InsertAfter(f2, n2)

	return frontier
}

// Checks, if a given triangle is valid according to a given outer point.
// For details see: https://en.wikipedia.org/wiki/Delaunay_triangulation
// The forth element is omitted because it is always 1 (and already initialized as 1).
func triangleIsValid_mine(a, b, c, d Vector) bool {

	g_matrixElements[0][0] = a.X
	g_matrixElements[1][0] = a.Y
	g_matrixElements[2][0] = a.X*a.X + a.Y*a.Y

	g_matrixElements[0][1] = b.X
	g_matrixElements[1][1] = b.Y
	g_matrixElements[2][1] = b.X*b.X + b.Y*b.Y

	g_matrixElements[0][2] = c.X
	g_matrixElements[1][2] = c.Y
	g_matrixElements[2][2] = c.X*c.X + c.Y*c.Y

	g_matrixElements[0][3] = d.X
	g_matrixElements[1][3] = d.Y
	g_matrixElements[2][3] = d.X*d.X + d.Y*d.Y

	// Making the eps smaller does not seem to create problems.
	return Fast4x4Determinant(&g_matrixElements) <= 0.001
}

func triangleIsValid(a, b, c, p Vector) bool {
	dx := a.X - p.X
	dy := a.Y - p.Y
	ex := b.X - p.X
	ey := b.Y - p.Y
	fx := c.X - p.X
	fy := c.Y - p.Y

	ap := dx*dx + dy*dy
	bp := ex*ex + ey*ey
	cp := fx*fx + fy*fy

	return dx*(ey*cp-bp*fy)-dy*(ex*cp-bp*fx)+ap*(ex*fy-ey*fx) < 0.001
}

func (d *Delaunay) flipEdge(e EdgeIndex) {

	e1 := d.Edges[e].ENext
	e2 := d.Edges[e1].ENext
	eT := d.Edges[e].ETwin
	eT1 := d.Edges[eT].ENext
	eT2 := d.Edges[eT1].ENext

	// Let's make sure that the face references are not affected by flipping the middle edge
	// by just assigning them to an edge that will stay with the same face reference anyway!
	// There could be an extra check for this but this should be not really slower and is correct every time.
	d.Faces[d.Edges[e].FFace].EEdge = e2
	d.Faces[d.Edges[eT].FFace].EEdge = eT2

	d.Edges[e2].ENext = eT1
	d.Edges[eT1].EPrev = e2

	d.Edges[eT2].ENext = e1
	d.Edges[e1].EPrev = eT2

	// Instead of creating new edges, we can just substitute the to-be-deleted edges
	d.Edges[e].VOrigin = d.Edges[eT2].VOrigin
	d.Edges[e].ENext = e2
	d.Edges[e].EPrev = eT1

	d.Edges[eT].VOrigin = d.Edges[e2].VOrigin
	d.Edges[eT].ENext = eT2
	d.Edges[eT].EPrev = e1

	// Reconnect other edges correctly
	d.Edges[e1].ENext = eT
	d.Edges[eT2].EPrev = eT
	d.Edges[eT1].ENext = e
	d.Edges[e2].EPrev = e

	// Faces rotate 90° counter clockwise. So have to be redefined for some edges.
	d.Edges[e1].FFace = d.Edges[eT].FFace
	d.Edges[eT1].FFace = d.Edges[e].FFace

}

// Recursively flips edges until all triangles are correct delaunay triangles.
// This should never affect the frontier, as only internal edges can be flipped.
func (d *Delaunay) legalizeTriangle(e EdgeIndex, v Vector, firstE EdgeIndex, firstTime bool) {

	// We might recursively itereate the same triangles. So if we arive at the very first triangle again,
	// we exit the recursion! Things are fine.
	if e == firstE && !firstTime { // && v.Equal(t.v, firstCheck.v) && !firstTime {
		return
	}

	// recursion ends at outer bounds.
	e1 := d.Edges[e].ENext
	if e1 == EmptyEdge {
		return
	}
	e2 := d.Edges[e1].ENext
	if e2 == EmptyEdge {
		return
	}

	a := d.Vertices[d.Edges[e].VOrigin].Pos
	b := d.Vertices[d.Edges[e1].VOrigin].Pos
	c := d.Vertices[d.Edges[e2].VOrigin].Pos

	hadRecursion := false

	if !triangleIsValid(a, b, c, v) {

		//fmt.Println(e)
		hadRecursion = true

		// Flip edge here
		d.flipEdge(e)

		// Recursively check neighboring triangles
		d.legalizeTriangle(d.Edges[e1].ETwin, v, firstE, false)
		d.legalizeTriangle(d.Edges[e2].ETwin, v, firstE, false)
		d.legalizeTriangle(d.Edges[d.Edges[e].EPrev].ETwin, c, firstE, false)
		d.legalizeTriangle(d.Edges[d.Edges[e1].EPrev].ETwin, c, firstE, false)

		//return
	}

	if firstTime && hadRecursion {
		//fmt.Println("")
	}

}

// Recursive function that walks right as long as there is a triangle to be created (to fill holes)
// that match the pi/2 angle criteria.
// TODO: Pass vector position directly without going through the frontier and delaunay every time.
// Returns, how many consecutive triangles were created to the right.
func (d *Delaunay) createConsecutiveTrianglesRight(baseVertex VertexIndex, leafNode *List, maxFillAngle float64) {

	//nextFrontier := d.frontier.Next(leafNode)
	nextFrontier := leafNode.Next

	leafNodeValue := leafNode.Value
	nextFrontierValue := nextFrontier.Value

	basePos := d.Vertices[baseVertex].Pos

	e0 := leafNodeValue.EdgeIndex
	e1 := nextFrontierValue.EdgeIndex

	v0 := d.Vertices[d.Edges[e0].VOrigin].Pos
	v1 := d.Vertices[d.Edges[e1].VOrigin].Pos

	isRight := IsRight2D(basePos, v0, v1)
	angle := AngleRad(Sub(basePos, v0), Sub(v1, v0))

	if isRight && angle <= maxFillAngle {

		// v1
		// |\
		// | \
		// |  \
		// |   \ eNew2 // This triangle is newly created here with already existing edges e0 and e1
		// |e1  \
		// |     \
		// v0-e0--baseVertex
		// |     /
		// |    /
		// |   /

		// Create one new triangle before calling recursively.
		f := d.createFace(Vector{-1, -1}, e0)
		d.Edges[e0].FFace = f
		d.Edges[e1].FFace = f

		eNew1 := d.createEdge(baseVertex, EmptyEdge, e0, e1, f, Edge{})
		eNew2 := d.createEdge(d.Edges[e1].VOrigin, eNew1, EmptyEdge, EmptyEdge, EmptyFace, Edge{})

		//f := d.createFace(Vector{}, eNew2)
		d.Edges[eNew1].ETwin = eNew2
		//d.Edges[eNew2].FFace = f

		d.Edges[e0].ENext = eNew1
		d.Edges[e1].EPrev = eNew1

		d.Edges[e0].EPrev = e1
		d.Edges[e1].ENext = e0

		//d.frontier.ChangeValue(nextFrontier, FrontElement{eNew2, nextFrontierValue.PolarAngle, nextFrontierValue.Radius})
		nextFrontier.Value.EdgeIndex = eNew2

		//fmt.Printf("%v --> %v\n", leafNode, leafNodeValue)
		//d.frontier.Delete(leafNodeValue)
		d.frontier.Delete(leafNode)

		// Validate towards the previous triangle created consecutively
		tp := d.Vertices[d.Edges[d.Edges[d.Edges[e0].ETwin].EPrev].VOrigin].Pos
		d.legalizeTriangle(e0, tp, e0, true)
		// Validate towards inside the triangulation
		tp = d.Vertices[d.Edges[d.Edges[d.Edges[e1].ETwin].EPrev].VOrigin].Pos
		d.legalizeTriangle(e1, tp, e1, true)

		d.createConsecutiveTrianglesRight(baseVertex, nextFrontier, maxFillAngle)

	}

}

// Recursive function that walks right as long as there is a triangle to be created (to fill holes)
// that match the pi/2 angle criteria.
// TODO: Pass vector position directly without going through the frontier and delaunay every time.
func (d *Delaunay) createConsecutiveTrianglesLeft(baseVertex VertexIndex, leafNode *List, maxFillAngle float64) {

	leafNodeValue := leafNode.Value

	//prevFrontier := d.frontier.Prev(leafNode)
	prevFrontier := leafNode.Prev
	prevFrontierValue := prevFrontier.Value

	basePos := d.Vertices[baseVertex].Pos

	e0 := leafNodeValue.EdgeIndex
	e1 := prevFrontierValue.EdgeIndex

	v0 := d.Vertices[d.Edges[e1].VOrigin].Pos
	v1Vertex := d.Edges[d.Edges[e1].ETwin].VOrigin
	v1 := d.Vertices[v1Vertex].Pos

	isLeft := IsLeft2D(basePos, v0, v1)
	angle := AngleRad(Sub(basePos, v0), Sub(v1, v0))

	if isLeft && angle <= maxFillAngle {

		f := d.createFace(Vector{-1, -1}, e0)
		d.Edges[e0].FFace = f
		d.Edges[e1].FFace = f

		eNew1 := d.createEdge(v1Vertex, EmptyEdge, e1, e0, f, Edge{})
		eNew2 := d.createEdge(baseVertex, eNew1, EmptyEdge, EmptyEdge, EmptyFace, Edge{})

		//f := d.createFace(Vector{}, eNew2)
		d.Edges[eNew1].ETwin = eNew2
		//d.Edges[eNew2].FFace = f

		d.Edges[e1].ENext = eNew1
		d.Edges[e0].EPrev = eNew1

		d.Edges[e0].ENext = e1
		d.Edges[e1].EPrev = e0

		//d.frontier.ChangeValue(leafNode, FrontElement{eNew2, leafNodeValue.PolarAngle, leafNodeValue.Radius})
		leafNode.Value.EdgeIndex = eNew2

		//d.frontier.Delete(prevFrontierValue)
		d.frontier.Delete(prevFrontier)

		// Validate towards the previous triangle created consecutively
		tp := d.Vertices[d.Edges[d.Edges[d.Edges[e0].ETwin].EPrev].VOrigin].Pos
		d.legalizeTriangle(e0, tp, e0, true)
		// Validate towards inside the triangulation
		tp = d.Vertices[d.Edges[d.Edges[d.Edges[e1].ETwin].EPrev].VOrigin].Pos
		d.legalizeTriangle(e1, tp, e1, true)

		d.createConsecutiveTrianglesLeft(baseVertex, leafNode, maxFillAngle)
	}
}

func extendByPoint(p DelaunayPoint, d *Delaunay, center Vector) {

	frontierItem := d.frontier.FindGreaterOrEqual(p.PolarAngle)

	frontierItemValue := frontierItem.Value

	vi := d.createVertex(p.Point)

	existingE := frontierItemValue.EdgeIndex

	//nextLeaf := d.frontier.Next(frontierItem)
	nextLeaf := frontierItem.Next
	nextLeafV := nextLeaf.Value
	hitVertex := false
	// If we have all vertices on a line, we don't want to create two triangles to somewhere.
	// We only create two triangles when the second one is well defined (has some angle)
	if math.Abs(frontierItemValue.PolarAngle-p.PolarAngle) <= EPS {
		nextE := nextLeafV.EdgeIndex

		p0 := d.Vertices[d.Edges[existingE].VOrigin].Pos
		p1 := d.Vertices[d.Edges[nextE].VOrigin].Pos

		v1 := Sub(p0, p.Point)
		v2 := Sub(p1, p.Point)

		if AngleRad(v1, v2) > EPS {
			hitVertex = true
		}
	}

	fi1 := d.createFace(Vector{-1, -1}, existingE)
	d.Edges[existingE].FFace = fi1

	twinV := d.Edges[d.Edges[existingE].ETwin].VOrigin
	ei1 := d.createEdge(twinV, EmptyEdge, existingE, EmptyEdge, fi1, Edge{})
	ei2 := d.createEdge(vi, ei1, EmptyEdge, EmptyEdge, EmptyFace, Edge{})
	d.Edges[ei1].ETwin = ei2

	d.Edges[existingE].ENext = ei1

	ej1 := d.createEdge(vi, EmptyEdge, ei1, existingE, fi1, Edge{})
	ej2 := d.createEdge(d.Edges[existingE].VOrigin, ej1, EmptyEdge, EmptyEdge, EmptyFace, Edge{})
	//fj2 := d.createFace(Vector{}, ej2)
	//d.Edges[ej2].FFace = fj2
	d.Edges[ej1].ETwin = ej2

	d.Edges[ei1].ENext = ej1
	d.Edges[existingE].EPrev = ej1

	// If we hit a vertex directly, we have to create a second triangle to the left as well so we don't have duplicate
	// keys in our d.frontier data structure!
	if hitVertex {

		prevVertex := d.Edges[nextLeafV.EdgeIndex].VOrigin

		fj2 := d.createFace(Vector{-1, -1}, ej2)
		d.Edges[ej2].FFace = fj2

		ek1 := d.createEdge(vi, EmptyEdge, ej2, nextLeafV.EdgeIndex, fj2, Edge{})
		ek2 := d.createEdge(prevVertex, ek1, EmptyEdge, EmptyEdge, EmptyFace, Edge{})
		//fk2 := d.createFace(Vector{}, ek2)
		//d.Edges[ek2].FFace = fk2
		d.Edges[ek1].ETwin = ek2

		d.Edges[ej2].ENext = ek1
		d.Edges[nextLeafV.EdgeIndex].EPrev = ek1
		d.Edges[nextLeafV.EdgeIndex].ENext = ej2
		d.Edges[nextLeafV.EdgeIndex].FFace = fj2
		d.Edges[ej2].EPrev = nextLeafV.EdgeIndex

		//ok := d.frontier.ChangeValue(nextLeaf, FrontElement{ek2, nextLeafV.PolarAngle, nextLeafV.Radius})
		nextLeaf.Value.EdgeIndex = ek2

		//ok = d.frontier.ChangeValue(frontierItem, FrontElement{ei2, frontierItemValue.PolarAngle, p.Distance})
		frontierItem.Value.EdgeIndex = ei2
		frontierItem.Value.Radius = p.Distance

		// There must exist a triangle behind the previous d.frontier edge. So this should be save.
		d2 := d.Vertices[d.Edges[d.Edges[d.Edges[nextLeafV.EdgeIndex].ETwin].EPrev].VOrigin].Pos
		d.legalizeTriangle(nextLeafV.EdgeIndex, d2, nextLeafV.EdgeIndex, true)

	} else {
		//ok := d.frontier.ChangeValue(frontierItem, FrontElement{ej2, frontierItemValue.PolarAngle, frontierItemValue.Radius})
		frontierItem.Value.EdgeIndex = ej2

		//fmt.Printf("Insert: %v\n", FrontElement{ei2, p.PolarAngle, p.Distance})
		//d.frontier.Insert(FrontElement{ei2, p.PolarAngle, p.Distance})

		d.frontier.InsertAfter(FrontElement{ei2, p.PolarAngle, p.Distance}, frontierItem.Prev)

	}

	// There must exist a triangle behind the previous d.frontier edge. So this should be save.
	tp := d.Vertices[d.Edges[d.Edges[d.Edges[existingE].ETwin].EPrev].VOrigin].Pos
	d.legalizeTriangle(existingE, tp, existingE, true)

	//   x
	//   |\
	// e | \
	// x |  \ ej1/ej2
	// i |   \
	// s |    \
	// t |     vi
	// i |    /
	// n |   /
	// g |  / ei1/ei2
	// E | /
	//   |/
	//   x

	maxFillAngle := DegToRad(120)

	//previousLeaf := d.frontier.Prev(frontierItem)
	previousLeaf := frontierItem.Prev

	d.createConsecutiveTrianglesRight(vi, frontierItem, maxFillAngle)

	// x                    // createConsecutiveTrianglesRight will create more triangles on this side
	// |\
	// | \
	// |  \ nextFrontier    // nextFrontier will not change during createConsecutiveTrianglesLeft()
	// |   \                // AFTER createConsecutiveTrianglesRight is done!
	// |    \
	// |     vi             // The initialy created triangle.
	// |    /
	// |   /
	// |  / f1              // f1 will not change during createConsecutiveTrianglesRight()
	// | /
	// |/                   // createConsecutiveTrianglesLeft will create triangles on this side!
	// x

	d.createConsecutiveTrianglesLeft(vi, previousLeaf, maxFillAngle)

}

func (d *Delaunay) triangulatePoints(pl *DelaunayPointList) {

	var lastP DelaunayPoint
	for i, p := range pl.Points {
		if i == 0 || !Equal(p.Point, lastP.Point) {
			extendByPoint(p, d, (*pl).Origin)
		}
		lastP = p
	}

	// Finalization step so we have a valid convex hull afterwards!
	f := d.frontier.GetLargestNode()
	fV := f.Value

	//prev := d.frontier.Prev(f)
	prev := f.Prev
	prevV := prev.Value

	frontierVertex := d.Edges[prevV.EdgeIndex].VOrigin
	d.createConsecutiveTrianglesRight(frontierVertex, f, DegToRad(180))

	for fV.PolarAngle >= prevV.PolarAngle {

		f = prev
		fV = prevV
		//prev = d.frontier.Prev(f)
		prev = f.Prev
		prevV = prev.Value

		frontierVertex = d.Edges[prevV.EdgeIndex].VOrigin
		d.createConsecutiveTrianglesRight(frontierVertex, f, DegToRad(180))
	}
	// Just to fill the gap from the last to the first vertex.
	// Todo: This could be make more general by changing the condition of the loop.
	// but the first attempts didn't quite work so I'll leave it here for the moment.
	f = prev
	fV = prevV
	//prev = d.frontier.Prev(f)
	prev = f.Prev
	prevV = prev.Value
	frontierVertex = d.Edges[prevV.EdgeIndex].VOrigin
	d.createConsecutiveTrianglesRight(frontierVertex, f, DegToRad(180))

}

func TriangulateMultithreaded(pointList []Vector, threadCount int) Delaunay {

	if len(pointList) < 3 {
		s := NewArrayMap(0)
		return Delaunay{
			Vertices:           []HEVertex{},
			firstFreeVertexPos: 0,
			Edges:              []HEEdge{},
			firstFreeEdgePos:   0,
			Faces:              []HEFace{},
			firstFreeFacePos:   0,
			frontier:           s,
		}
	}

	dPointList := preparePointList(pointList, threadCount)
	//defer fmt.Printf("point to be inserted last: %v\n", dPointList.Points[52])
	//dPointList.Points = dPointList.Points[:7]

	//fmt.Printf("Center: %v\n", dPointList.Origin)

	// Any planar triangulation: total degree == 3f + k = 2e
	// With k == points on convex hull of all points. Let k = p
	// https://en.wikipedia.org/wiki/Planar_graph#Euler.27s_formula
	// Eulers formula states: v - e + f = 2
	// For a finite, sparse, connected graph, we have:
	// e <= 3v - 6
	edgeCount := (3*len(dPointList.Points) - 6) * 2
	delaunay := Delaunay{
		Vertices:           make([]HEVertex, len(dPointList.Points)),
		firstFreeVertexPos: 0,
		Edges:              make([]HEEdge, edgeCount),
		firstFreeEdgePos:   0,
		Faces:              make([]HEFace, edgeCount/3, edgeCount/3),
		firstFreeFacePos:   0,
	}
	delaunay.frontier = delaunay.initializeTriangulation(&dPointList)

	//frontier := delaunay.initializeTriangulation(&dPointList)

	delaunay.triangulatePoints(&dPointList)

	//delaunay.Edges = append([]HEEdge{}, delaunay.Edges[:delaunay.firstFreeEdgePos]...)
	//delaunay.Vertices = append([]HEVertex{}, delaunay.Vertices[:delaunay.firstFreeVertexPos]...)
	//delaunay.Faces = append([]HEFace{}, delaunay.Faces[:delaunay.firstFreeFacePos]...)

	delaunay.Edges = delaunay.Edges[:delaunay.firstFreeEdgePos]
	delaunay.Vertices = delaunay.Vertices[:delaunay.firstFreeVertexPos]
	delaunay.Faces = delaunay.Faces[:delaunay.firstFreeFacePos]

	//fmt.Printf("edges len: %d, cap: %v\n", cap(delaunay.Edges), len(delaunay.Edges))

	//newEdges := make([]HEEdge, delaunay.firstFreeEdgePos, delaunay.firstFreeEdgePos)
	//copied := copy(newEdges, delaunay.Edges)
	//fmt.Printf("edges copied: %d. Planned to copy: %d\n", copied, delaunay.firstFreeEdgePos)
	//delaunay.Edges = newEdges

	return delaunay

}

func Triangulate(pointList []Vector) Delaunay {
	return TriangulateMultithreaded(pointList, 4)
}

// Extracts a list of simple edge representations of all existing edges from the Delaunay triangulation.
// This list does not contain any duplicates from the original half-edge data structure and should only
// be used for drawing or operations that need all edges anyway.
func (d *Delaunay) ExtractEdgeList() []SimpleEdge {

	edges := []SimpleEdge{}

	for i, e := range d.Edges {
		// if the index of the twin edge is smaller than i, we already included this edge.
		if int(e.ETwin) > i {

			var p1 Vector
			var p2 Vector
			if e.VOrigin != EmptyVertex {
				p1 = d.Vertices[e.VOrigin].Pos
			} else {
				// Only one vertex can ever be missing. Never both.
				p1 = d.Vertices[d.Edges[e.ETwin].VOrigin].Pos
				p1.Add(Mult(e.TmpEdge.Dir, -10))
			}

			if d.Edges[e.ETwin].VOrigin != EmptyVertex {
				p2 = d.Vertices[d.Edges[e.ETwin].VOrigin].Pos
			} else {
				p2 = d.Vertices[e.VOrigin].Pos
				p2.Add(Mult(e.TmpEdge.Dir, -10))
			}

			edges = append(edges, SimpleEdge{p1, p2})
		}
	}

	return edges
}

// List of points that represent the convex hull.
func (d *Delaunay) ExtractConvexHull() ConvexHull {

	ch := ConvexHull{}

	largest := d.frontier.GetLargestNode()

	node := d.frontier.GetSmallestNode()
	for node != largest {
		//e := node.GetValue().(FrontElement).EdgeIndex
		e := node.Value.EdgeIndex
		ch = append(ch, d.Vertices[d.Edges[e].VOrigin].Pos)

		//node = d.frontier.Next(node)
		node = node.Next
	}

	//e := largest.GetValue().(FrontElement).EdgeIndex
	e := largest.Value.EdgeIndex
	ch = append(ch, d.Vertices[d.Edges[e].VOrigin].Pos)

	return ch

}

// Expects an Inner edge! Not the outward facing one!
func calcPerpendicularBisector(d *Delaunay, e EdgeIndex) Edge {
	p1 := d.Vertices[d.Edges[e].VOrigin].Pos
	p2 := d.Vertices[d.Edges[d.Edges[e].ETwin].VOrigin].Pos

	diff := Sub(p2, p1)
	//z := Vector{0, 0}
	//cross := v.Cross(diff, z)
	cross := Perpendicular(diff)
	diff = Mult(diff, 0.5)
	middle := Add(p1, diff)
	newPoint := Add(middle, cross)
	dir := Sub(middle, newPoint)
	return Edge{
		Pos: middle,
		Dir: dir,
	}
}

// Thanks to Paul Draper at
// http://stackoverflow.com/questions/20677795/find-the-point-of-intersecting-lines
func calculateVertexPosition(e1, e2 Edge) Vector {
	p12 := Add(e1.Pos, e1.Dir)
	p22 := Add(e2.Pos, e2.Dir)

	xdiff := Vector{e1.Pos.X - p12.X, e2.Pos.X - p22.X}
	ydiff := Vector{e1.Pos.Y - p12.Y, e2.Pos.Y - p22.Y}

	det2D := func(a, b Vector) float64 {
		return a.X*b.Y - a.Y*b.X
	}

	div := det2D(xdiff, ydiff)
	if math.Abs(div) <= EPS {
		fmt.Printf("Lines do not intersect!\n")
		return Vector{}
	}

	d := Vector{det2D(e1.Pos, p12), det2D(e2.Pos, p22)}
	x := det2D(d, xdiff) / div
	y := det2D(d, ydiff) / div
	return Vector{x, y}

}

// Verifies that the Voronoi is not corrupted or has miscalculated edges/faces
// or any other unvalid stuff.
func (v *Voronoi) Verify() error {

	for i, e := range v.Edges {
		if e != EmptyE {

			// Every valid edge MUST have a valid Twin edge!
			if e.ETwin == EmptyEdge || v.Edges[e.ETwin] == EmptyE {
				return errors.New(fmt.Sprintf("Edge %v: %v has an invalid twin edge", i, e))
			}

			// Twins must refer to each other!
			if i != int(v.Edges[e.ETwin].ETwin) {
				return errors.New(fmt.Sprintf("Edge %v and his Twin %v don't refer to each other", i, e.ETwin))
			}

			// Check if the origin vertex is valid (if there is one)
			if e.VOrigin != EmptyVertex && v.Vertices[e.VOrigin] == EmptyV {
				return errors.New(fmt.Sprintf("Origin vertex of edge %v: %v is invalid", i, e))
			}

			// Check, if the face is valid (if defined)
			if e.FFace == EmptyFace || v.Faces[e.FFace] == EmptyF {
				return errors.New(fmt.Sprintf("Face of edge %v: %v is invalid", i, e.FFace))
			}

			// Check, if the next edge is valid
			if e.ENext != EmptyEdge && v.Edges[e.ENext] == EmptyE {
				return errors.New(fmt.Sprintf("Next edge for edge %v: %v is invalid", i, e))
			}

			// Check, if the prev edge is valid
			if e.EPrev != EmptyEdge && v.Edges[e.EPrev] == EmptyE {
				return errors.New(fmt.Sprintf("Prev edge for edge %v: %v is invalid", i, e))
			}

			// If this->next is not empty, the next one must have a previous edge.
			if e.ENext != EmptyEdge && v.Edges[e.ENext].EPrev == EmptyEdge {
				return errors.New(fmt.Sprintf("Next edge for edge %v must have a non-empty prev", i))
			}

			// The this edge must correspond to the next->prev edge
			if e.ENext != EmptyEdge && v.Edges[e.ENext].EPrev != EmptyEdge && v.Edges[e.ENext].EPrev != EdgeIndex(i) {
				return errors.New(fmt.Sprintf("The edge %v has %v as next edge, which has %v as his previous. They must reference each other!", i, e.ENext, v.Edges[e.ENext].EPrev))
			}
		}
	}

	for i, f := range v.Faces {
		if f != EmptyF {

			// Check for valid EEdge
			if f.EEdge == EmptyEdge || v.Edges[f.EEdge] == EmptyE {
				return errors.New(fmt.Sprintf("Face %v: %v points to invalid edge", i, f))
			}

			// If the edge is the first one (no/infinite start-vertex), it's all good.
			// Otherwise check, that the edges actually go all the way around the face!
			startEdge := f.EEdge
			e := v.Edges[startEdge].ENext

			for e != EmptyEdge && e != startEdge {

				if v.Edges[e].FFace != FaceIndex(i) {
					return errors.New(fmt.Sprintf("Edge %v of the face %v: %v does not point to the right face!", e, i, f))
				}

				e = v.Edges[e].ENext
			}

			if e != startEdge && v.Edges[startEdge].EPrev != EmptyEdge {
				return errors.New(fmt.Sprintf("Edge %v does not correctly loop around face %v (All around, it should point to edge: %v)", e, f, startEdge))
			}

		}
	}

	return nil
}

// The v is actually the Voronoi not Delaunay. But before casting.
func (v *Delaunay) connectToExistingVoronoi(d *Delaunay, outgoingEdges [][]EdgeIndex, de1, e1, e2, e3 EdgeIndex, b1 Edge, df FaceIndex) {

	// If there is no neighboring triangle/face, we are at the outer bounds
	f1 := d.Edges[d.Edges[de1].ETwin].FFace
	if f1 == EmptyFace {
		// Again: Because we create all the Voronoi faces in a loop at the beginning with the Delaunay Vertices as reference (position),
		// the Delaunay vertex ids correspond to the Voronoi face ids! So we can just use them where we know it fits.
		e1T := v.createEdge(EmptyVertex, e1, EmptyEdge, e3, FaceIndex(d.Edges[de1].VOrigin), b1)
		v.Edges[e1].ETwin = e1T
		v.Edges[e3].EPrev = e1T

		// The incoming edge (just created) will be the reference for the face per definition (first edge in incomplete polygon!)
		v.Faces[d.Edges[de1].VOrigin].EEdge = e1T

	} else {
		// If f1 is smaller than the current face we are iterating, we already created all outgoing edges for this triangle. Meaning, we can just connect!
		if f1 < df {

			var incomingEdge EdgeIndex
			adjacentEdge := d.Faces[f1].EEdge
			var i int
			for i = 0; i < 3; i++ {

				// Found the adjacent delaunay edge. So we also know, which is the incoming voronoi edge because they are always! created in the same
				// order, starting with the edge the face points to and going to .next edges.
				if d.Edges[adjacentEdge].ETwin == de1 {
					incomingEdge = outgoingEdges[f1][i]
					break
				}

				adjacentEdge = d.Edges[adjacentEdge].ENext
			}

			v.Edges[e1].ETwin = incomingEdge
			v.Edges[incomingEdge].ETwin = e1

			v.Edges[incomingEdge].ENext = e3
			v.Edges[e3].EPrev = incomingEdge

			// Find the previous outgoing edge from the adjacent triangle to connect e1.next to.
			prevOutgoingEdge := outgoingEdges[f1][(i+2)%3]

			v.Edges[e1].ENext = prevOutgoingEdge
			v.Edges[prevOutgoingEdge].EPrev = e1
		}
	}

}

func (d *Delaunay) CreateVoronoi() Voronoi {

	// The data structure is identical. So we create this Delaunay struct and cast it to Voronoi on return.
	// So we can use the corresponding methods.
	edgeCount := (3*len(d.Vertices) - 6) * 2
	v := Delaunay{
		Vertices:           make([]HEVertex, (edgeCount+6)/3),
		firstFreeVertexPos: 0,
		Edges:              make([]HEEdge, edgeCount),
		firstFreeEdgePos:   0,
		Faces:              make([]HEFace, len(d.Vertices)),
		firstFreeFacePos:   0,
	}

	// Temporary slice of the three outgoing edges for each delaunay face for fast voronoi edge lookup
	outgoingEdges := make([][]EdgeIndex, len(d.Faces))

	for i, _ := range outgoingEdges {
		outgoingEdges[i] = make([]EdgeIndex, 3, 3)
	}

	// Create all Voronoi faces beforehand at the position of the Delaunay vertex. The connection to the edges will be added later.
	// The index of the Voronoi faces will correspond to the index Delaunay vertices.
	for _, dv := range d.Vertices {
		v.createFace(dv.Pos, EmptyEdge)
	}

	// Iterate over all Delaunay Faces. All Delaunay triangles have exactly three edges!
	for i, df := range d.Faces {

		if df == EmptyF {
			break
		}

		de1 := df.EEdge
		de2 := d.Edges[de1].ENext
		de3 := d.Edges[de2].ENext
		b1 := calcPerpendicularBisector(d, de1)
		b2 := calcPerpendicularBisector(d, de2)
		b3 := calcPerpendicularBisector(d, de3)

		newVertexPos := calculateVertexPosition(b1, b2)
		newVertex := v.createVertex(newVertexPos)

		// Create three outgoing voronoi edges perpendicular to the delaunay edges of the triangle
		// Face: Instead of calculating d1.twin.origin we can just take d2.origin for the same result!
		e1 := v.createEdge(newVertex, EmptyEdge, EmptyEdge, EmptyEdge, FaceIndex(d.Edges[de2].VOrigin), b1)
		e2 := v.createEdge(newVertex, EmptyEdge, EmptyEdge, EmptyEdge, FaceIndex(d.Edges[de3].VOrigin), b2)
		e3 := v.createEdge(newVertex, EmptyEdge, EmptyEdge, EmptyEdge, FaceIndex(d.Edges[de1].VOrigin), b3)

		// Assign the edge reference for a face (if not already done!)
		if v.Faces[d.Edges[de1].VOrigin].EEdge == EmptyEdge {
			v.Faces[d.Edges[de1].VOrigin].EEdge = e3
		}
		if v.Faces[d.Edges[de2].VOrigin].EEdge == EmptyEdge {
			v.Faces[d.Edges[de2].VOrigin].EEdge = e1
		}
		if v.Faces[d.Edges[de3].VOrigin].EEdge == EmptyEdge {
			v.Faces[d.Edges[de3].VOrigin].EEdge = e2
		}

		v.connectToExistingVoronoi(d, outgoingEdges, de1, e1, e2, e3, b1, FaceIndex(i))
		v.connectToExistingVoronoi(d, outgoingEdges, de2, e2, e3, e1, b2, FaceIndex(i))
		v.connectToExistingVoronoi(d, outgoingEdges, de3, e3, e1, e2, b3, FaceIndex(i))

		// For fast lookup later
		outgoingEdges[i][0] = e1
		outgoingEdges[i][1] = e2
		outgoingEdges[i][2] = e3

	}

	realV := Voronoi(v)

	realV.Edges = append([]HEEdge{}, realV.Edges[:realV.firstFreeEdgePos]...)
	realV.Vertices = append([]HEVertex{}, realV.Vertices[:realV.firstFreeVertexPos]...)
	realV.Faces = append([]HEFace{}, realV.Faces[:realV.firstFreeFacePos]...)

	//if e := realV.Verify(); e != nil {
	//	fmt.Println(e)
	//}

	return realV
}
