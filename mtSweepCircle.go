package mtSweepCircle

import (
    "fmt"
    //"os"
    he "mtHalfEdge"
    v "mtVector"
    "errors"
    "sort"
    "math"
    //t "github.com/MauriceGit/tree23"
    t "tree23"
    //"time"
)

const (
    EPS float64 = 0.0000001
)

type Delaunay struct {
    Vertices           []he.HEVertex
    FirstFreeVertexPos he.VertexIndex
    Edges              []he.HEEdge
    FirstFreeEdgePos   he.EdgeIndex
    Faces              []he.HEFace
    FirstFreeFacePos   he.FaceIndex
}

type DelaunayPoint struct {
    Point       v.Vector
    Distance    float64
    PolarAngle  float64
}
type DelaunayPointList struct {
    Points      []DelaunayPoint
    Origin      v.Vector
}

type FrontElement struct {
    //// Index into Delaunay.vertices data structure!
    //Index       he.VertexIndex
    // Index into Delaunay.edges data structure!
    EdgeIndex   he.EdgeIndex
    // Polar angle of the vertex of EdgeI
    PolarAngle  float64
    Radius      float64
}

// Used for building a stack of triangle checks so we don't infinitely flip edges.
type CheckTriangle struct {
    e   he.EdgeIndex
    v   v.Vector
}

// The frontier is a list of indices into the Delaunay.vertices data structure!
// The frontier should always be sorted according to their polar angle to allow optimized search in the list
// when including new points!!!
// type Frontier []FrontElement
//type Frontier *t.Tree23

var EmptyF = he.HEFace{}
var EmptyE = he.HEEdge{}
var EmptyV = he.HEVertex{}

// We will need this matrix millions of times for bigger Delaunay. So we can cache the memory
// and just reuse it.
var g_matrixElements = [4][4]float64{{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}}


func (dp DelaunayPoint) String() string {
    return fmt.Sprintf("{Point: %v, Distance: %v, PolarAngle: %v}", dp.Point, dp.Distance, dp.PolarAngle)
}
func (dpl DelaunayPointList) String() string {
    return fmt.Sprintf("\nDelaunayPointList:\n    Origin: %v\n    Points: %v\n", dpl.Origin, dpl.Points)
}

//func (f *Frontier) insertAt(fp FrontElement, i int) {
//    *f = append(*f, FrontElement{})
//    copy((*f)[i+1:], (*f)[i:])
//    (*f)[i] = fp
//}

// Len is part of sort.Interface.
//func (f *Frontier) Len() int {
//    return len(*f)
//}
// Smaller is part of the Searchable interface
//func (f *Frontier) Smaller(e interface{}, i int) bool {
//    return e.(float64) < (*f)[i].PolarAngle
//}
// Match is part of the Searchable interface
//func (f *Frontier) Match(e interface{}, i int) bool {
//
//    if i > 0 {
//        return  (*f)[i].PolarAngle   > e.(float64) &&
//                (*f)[i-1].PolarAngle < e.(float64)
//    }
//
//    return (*f)[i].PolarAngle > e.(float64)
//}
//// Match is part of the SearchableInterpolation interface
//func (f *Frontier) GetValue(i int) float64 {
//    return float64((*f)[i].PolarAngle)
//}

// Equal is part of the Tree23 interface
func (f FrontElement) Equal(f2 t.TreeElement) bool {

    fe := f2.(FrontElement)
    return math.Abs(f.PolarAngle - fe.PolarAngle) <= EPS && math.Abs(f.Radius - fe.Radius) <= EPS
    //return math.Abs(f.PolarAngle - f2.(FrontElement).PolarAngle) <= EPS
}
// ExtractValue is part of the Tree23 interface
func (f FrontElement) ExtractValue() float64 {
    return f.PolarAngle
}

// By is the type of a "less" function that defines the ordering of its Planet arguments.
type By func(i, j *DelaunayPoint) bool
// planetSorter joins a By function and a slice of Planets to be sorted.
type pointSorter struct {
    points []DelaunayPoint
    by     func(i,j *DelaunayPoint) bool // Closure used in the Less method.
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

////////////////////////////////////////////////////////////////////////
//  Pretty Print the  Voronoi Attributes
////////////////////////////////////////////////////////////////////////

func (v *Delaunay) pprint() {

    fmt.Println("Voronoi:")
    fmt.Printf("   Vertices (%v):\n", v.FirstFreeVertexPos)

    for i,ve := range v.Vertices {
        if ve != EmptyV {
            fmt.Printf("\t%v:\tPos: %v\n", i, ve.Pos)
        }
    }
    fmt.Printf("\n")

    fmt.Printf("   Edges    (%v):\n", v.FirstFreeEdgePos)
    for i,e := range v.Edges {
        if e != EmptyE {
            fmt.Printf("\t%v:\tOrigin: %v,\tTwin: %v,\tPrev: %v,\tNext: %v,\tFace: %v\n", i, e.VOrigin, e.ETwin, e.EPrev, e.ENext, e.FFace)
        }
    }
    fmt.Printf("\n")

    fmt.Printf("   Faces    (%v):\n", v.FirstFreeFacePos)
    for i,f := range v.Faces {
        if f != EmptyF {
            fmt.Printf("\t%v:\tRefPoint: %v,\tEdge:%v\n", i, f.ReferencePoint, f.EEdge)
        }
    }
    fmt.Printf("\n")
}

func print(frontier *t.Tree23, delaunay *Delaunay, title string) {

    //frontier.Pprint()

    //fmt.Printf("%s\n\t", title)
    //for _,f := range *frontier {
    //    //fmt.Printf("(%d - %.1f), ", f.EdgeIndex, f.PolarAngle)
    //    fmt.Printf("%d, ", delaunay.Edges[f.EdgeIndex].VOrigin)
    //}
    //fmt.Printf("\n")
}

// Verifies that the Delaunay is not corrupted or has miscalculated edges/faces
// or any other unvalid stuff.
func (de *Delaunay) Verify() error {

    for i,e := range de.Edges {
        if e != EmptyE {

            // Every valid edge MUST have a valid Twin edge!
            if e.ETwin == he.EmptyEdge || de.Edges[e.ETwin] == EmptyE {
                return errors.New(fmt.Sprintf("Edge %de: %de has an invalid twin edge", i, e))
            }

            // Twins must refer to each other!
            if i != int(de.Edges[e.ETwin].ETwin) {
                return errors.New(fmt.Sprintf("Edge %de and his Twin %de don't refer to each other", i, e.ETwin))
            }

            // Check if the origin vertex is valid (if there is one)
            if e.VOrigin != he.EmptyVertex && e.VOrigin != he.InfiniteVertex && de.Vertices[e.VOrigin] == EmptyV {
                return errors.New(fmt.Sprintf("Origin vertex of edge %de: %de is invalid", i, e))
            }

            // Check, if the face is valid
            if e.FFace == he.EmptyFace || de.Faces[e.FFace] == EmptyF {
                return errors.New(fmt.Sprintf("Face of edge %de: %de is invalid", i, e))
            }

            // if VOrigin is Empty, it MUST be the referenced edge for the face of the edge!
            if (e.VOrigin == he.EmptyVertex || e.VOrigin == he.InfiniteVertex) && i != int(de.Faces[e.FFace].EEdge) {
                return errors.New(fmt.Sprintf("A first edge %de: %de must be referenced as first edge by the corresponding face %de!", e, i, e.FFace))
            }

            // Check, if the next edge is valid
            if e.ENext != he.EmptyEdge && de.Edges[e.ENext] == EmptyE {
                return errors.New(fmt.Sprintf("Next edge for edge %de: %de is invalid", i, e))
            }

            // Check, if the prev edge is valid
            if e.EPrev != he.EmptyEdge && de.Edges[e.EPrev] == EmptyE {
                return errors.New(fmt.Sprintf("Prev edge for edge %de: %de is invalid", i, e))
            }

            // If the prev edge is valid
            if e.VOrigin != he.EmptyVertex && e.VOrigin != he.InfiniteVertex && e.EPrev == he.EmptyEdge {
                return errors.New(fmt.Sprintf("Prev edge for edge %de: %de is not defined but should be", i, e))
            }

            // If this->next is not empty, the next one must have a previous edge.
            if e.ENext != he.EmptyEdge && de.Edges[e.ENext].EPrev == he.EmptyEdge {
                return errors.New(fmt.Sprintf("Next edge for edge %de must have a non-empty prev", i))
            }

            // The this edge must correspond to the next->prev edge
            if e.ENext != he.EmptyEdge && de.Edges[e.ENext].EPrev != he.EmptyEdge && de.Edges[e.ENext].EPrev != he.EdgeIndex(i) {
                return errors.New(fmt.Sprintf("The edge %de has %de as next edge, which has %de as his previous. They must reference each other!", i, e.ENext, de.Edges[e.ENext].EPrev))
            }



        }
    }

    for i,f := range de.Faces {
        if f != EmptyF {

            // Check for valid reference point
            if f.ReferencePoint == v.InfinitePoint {
                return errors.New(fmt.Sprintf("Face %de: %de has infinite reference point", i, f))
            }

            // Check for valid EEdge
            if f.EEdge == he.EmptyEdge || de.Edges[f.EEdge] == EmptyE {
                return errors.New(fmt.Sprintf("Face %de: %de points to invalid edge", i, f))
            }

            // If the edge is the first one (no/infinite start-vertex), it's all good.
            // Otherwise check, that the edges actually go all the way around the face!
            if de.Edges[f.EEdge].VOrigin != he.EmptyVertex && de.Edges[f.EEdge].VOrigin != he.InfiniteVertex {
                startEdge := f.EEdge
                edgeCount := 0
                e := de.Edges[startEdge].ENext

                for e != he.EmptyEdge && e != startEdge {

                    if int(de.Edges[e].FFace) != i {
                        return errors.New(fmt.Sprintf("Edge %de of the face %de: %de does not point to the right face!", e, i, f))
                    }

                    if edgeCount > 50 {
                        return errors.New(fmt.Sprintf("Looping around the edges of the face %de: %de most likely does not stop (infinite loop)", i, f))
                    }

                    e = de.Edges[e].ENext
                    edgeCount += 1
                }
            }
        }
    }

    return nil
}

// So we can get a pointer of some data structure? What about scope issues?
func (v *Delaunay)createFace(refPoint v.Vector, eEdge he.EdgeIndex) he.FaceIndex {
    v.Faces[v.FirstFreeFacePos] = he.HEFace {
        ReferencePoint: refPoint,
        EEdge:          eEdge,
    }
    v.FirstFreeFacePos += 1
    return v.FirstFreeFacePos-1
}

// So we can get a pointer of some data structure? What about scope issues?
func (v *Delaunay)createEdge(vOrigin he.VertexIndex, eTwin, ePrev, eNext he.EdgeIndex, fFace he.FaceIndex, tmpEdge v.Edge) he.EdgeIndex {

    index := v.FirstFreeEdgePos
    v.FirstFreeEdgePos += 1

    v.Edges[index] = he.HEEdge {
        VOrigin:    vOrigin,
        ETwin:      eTwin,
        ENext:      eNext,
        EPrev:      ePrev,
        FFace:      fFace,
        TmpEdge:    tmpEdge,
    }
    return index
}

// So we can get a pointer of some data structure? What about scope issues?
func (v *Delaunay)createVertex(pos v.Vector) he.VertexIndex {
    index := v.FirstFreeVertexPos
    v.FirstFreeVertexPos += 1

    v.Vertices[index] = he.HEVertex {
        Pos:        pos,
    }

    return index
}

func calcPolarAngle(p v.Vector, origin v.Vector) float64 {

    if v.Equal(p, origin) {
        return 0.0
    }

    diff  := v.Sub(p, origin)
    xAxis := v.Vector{1,0,0}

    angle := v.Angle(diff, xAxis)
    if diff.Y < 0 {
        angle = 180 + (180-angle)
    }

    return angle
}

func preparePointList(points *v.PointList) DelaunayPointList {
    var dPointList DelaunayPointList

    dPointList.Points = make([]DelaunayPoint, points.Len(), points.Len())

    // We cannot calculate the origin in the same loop as the distance, as distance depends on the final origin...
    origin := v.Vector{}
    for _,p := range *points {
        origin.Add(p)
    }
    origin.Div(float64(points.Len()))
    dPointList.Origin = origin

    for i,p := range *points {
        dPointList.Points[i] = DelaunayPoint{
                                             Point: p,
                                             Distance: v.Length(v.Sub(p, dPointList.Origin)),
                                             PolarAngle: calcPolarAngle(p, dPointList.Origin),
                                            }
    }

    distance := func(v1,v2 *DelaunayPoint) bool {
        // sort by radius
        if math.Abs(v1.Distance - v2.Distance) >= EPS {
            return v1.Distance < v2.Distance
        }
        // if radius is equal, sort by polar angle
        return v1.PolarAngle < v2.PolarAngle
    }

    By(distance).Sort(dPointList.Points)

    // This is very hacky and should be improved first thing, when the algorithms works...
    // There should be a way to determine the first triangle where Origin is in without sorting/recalculating the list twice!!!
    dPointList.Origin = v.Div(v.Add(v.Add(dPointList.Points[0].Point, dPointList.Points[1].Point), dPointList.Points[2].Point), 3.0)
    // This is useless.
    for i,p := range dPointList.Points {
        dPointList.Points[i].Distance   = v.Length(v.Sub(p.Point, dPointList.Origin))
        dPointList.Points[i].PolarAngle = calcPolarAngle(p.Point, dPointList.Origin)

        //fmt.Printf("PointList: angle:%v, point:%v\n", dPointList.Points[i].PolarAngle, dPointList.Points[i].Point)
    }
    // Like ... Really.

    //start := time.Now()
    By(distance).Sort(dPointList.Points)
    //binTime := time.Since(start).Nanoseconds()

    //fmt.Printf("Sort in Seconds: %.8f\n", float64(binTime)/1000000000.0)


    return dPointList
}

func (delaunay *Delaunay)initializeTriangulation(pl *DelaunayPointList) *t.Tree23 {

    //var frontier Frontier = make([]FrontElement, 3, len((*pl).Points))

    f0 := FrontElement{}
    f1 := FrontElement{}
    f2 := FrontElement{}

    d := delaunay

    // The closest three points are (per definition) the triangle surrounding origin.
    if pl.Points[0].PolarAngle < pl.Points[1].PolarAngle {
        // 0 < 1
        delaunay.Vertices[0] = he.HEVertex{pl.Points[0].Point}
        f0.PolarAngle = pl.Points[0].PolarAngle
        if pl.Points[2].PolarAngle < pl.Points[0].PolarAngle || pl.Points[2].PolarAngle > pl.Points[1].PolarAngle {
            // 2 < 0 < 1
            delaunay.Vertices[1] = he.HEVertex{pl.Points[1].Point}
            delaunay.Vertices[2] = he.HEVertex{pl.Points[2].Point}
            f1.PolarAngle = pl.Points[1].PolarAngle
            f2.PolarAngle = pl.Points[2].PolarAngle
        } else {
            // 0 < 2 < 1
            delaunay.Vertices[1] = he.HEVertex{pl.Points[2].Point}
            delaunay.Vertices[2] = he.HEVertex{pl.Points[1].Point}
            f1.PolarAngle = pl.Points[2].PolarAngle
            f2.PolarAngle = pl.Points[1].PolarAngle
        }
    } else {
        delaunay.Vertices[0] = he.HEVertex{pl.Points[2].Point}
        f0.PolarAngle = pl.Points[2].PolarAngle
        // 1 < 0
        if pl.Points[2].PolarAngle < pl.Points[1].PolarAngle || pl.Points[2].PolarAngle > pl.Points[0].PolarAngle {
            // 2 < 1 < 0
            delaunay.Vertices[1] = he.HEVertex{pl.Points[1].Point}
            delaunay.Vertices[2] = he.HEVertex{pl.Points[0].Point}
            f1.PolarAngle = pl.Points[1].PolarAngle
            f2.PolarAngle = pl.Points[0].PolarAngle
        } else {
            // 1 < 2 < 0
            delaunay.Vertices[1] = he.HEVertex{pl.Points[0].Point}
            delaunay.Vertices[2] = he.HEVertex{pl.Points[1].Point}
            f1.PolarAngle = pl.Points[0].PolarAngle
            f2.PolarAngle = pl.Points[1].PolarAngle
        }
    }
    delaunay.FirstFreeVertexPos = 3

    f := delaunay.createFace(v.Vector{}, 0)

    // Edge 0-->1
    e   := delaunay.createEdge(0, he.EmptyEdge, he.EmptyEdge, he.EmptyEdge, f, v.Edge{})
    // Edge 1-->0
    eT  := delaunay.createEdge(1, e, he.EmptyEdge, he.EmptyEdge, he.EmptyFace, v.Edge{})
    d.Edges[e].ETwin = eT
    d.Faces[f].EEdge = e
    feT := d.createFace(v.Vector{}, eT)
    d.Edges[eT].FFace = feT

    // Edge 1-->2
    e2  := delaunay.createEdge(1, he.EmptyEdge, e, he.EmptyEdge, f, v.Edge{})
    // Edge 2-->1
    e2T := delaunay.createEdge(2, e2, he.EmptyEdge, he.EmptyEdge, he.EmptyFace, v.Edge{})
    d.Edges[e2].ETwin = e2T
    fe2T := d.createFace(v.Vector{}, e2T)
    d.Edges[e2T].FFace = fe2T

    // Edge 2-->0
    e3  := delaunay.createEdge(2, he.EmptyEdge, e2, e, f, v.Edge{})
    // Edge 0-->2
    e3T := delaunay.createEdge(0, e3, he.EmptyEdge, he.EmptyEdge, he.EmptyFace, v.Edge{})
    d.Edges[e3].ETwin = e3T
    fe3T := d.createFace(v.Vector{}, e3T)
    d.Edges[e3T].FFace = fe3T

    d.Edges[e].EPrev = e3
    d.Edges[e].ENext = e2
    d.Edges[e2].ENext = e3

    // TODO: Use precalculated polar angles instead of recalculating them. Should be straight forward.
    f0.EdgeIndex  = eT
    f0.PolarAngle = calcPolarAngle(delaunay.Vertices[1].Pos, (*pl).Origin)
    f1.EdgeIndex = e2T
    f1.PolarAngle = calcPolarAngle(delaunay.Vertices[2].Pos, (*pl).Origin)
    f2.EdgeIndex = e3T
    f2.PolarAngle = calcPolarAngle(delaunay.Vertices[0].Pos, (*pl).Origin)

    // Pop first three points, because they are already triangulated by default
    pl.Points = pl.Points[3:]

    frontier := t.New()

    fmt.Printf("pos: %v, origin: %v\n", delaunay.Vertices[1].Pos, (*pl).Origin)

    frontier.Insert(f0)
    fmt.Printf("Insert1: %v\n", f0)
    frontier.Insert(f1)
    fmt.Printf("Insert1: %v\n", f1)
    frontier.Insert(f2)
    fmt.Printf("Insert1: %v\n", f2)

    return frontier
}

// To project a point onto the frontier, we only need to know the polar angle, as the frontier
// is sorted by polar angle, starting with 0.
// The function returns both vertex indices in between the new point is projected (counter clockwise!)
// TODO: Use or implement quadratic binary search/binary search for optimized location finding. Otherwise we might result in O(n²).
//func (frontier Frontier)findFrontierPosition(polarAngle float64) int {
//
//    index, err := advsearch.InterpolationSearch(frontier, polarAngle)
//    if err != nil {
//        return 0
//    }
//    return index
//}



// Checks, if a given triangle is valid according to a given outer point.
// For details see: https://en.wikipedia.org/wiki/Delaunay_triangulation
// The forth element is omitted because it is always 1 (and already initialized as 1).
func triangleIsValid(a, b, c, d v.Vector) bool {

    g_matrixElements[0][0] = a.X
    g_matrixElements[1][0] = a.Y
    g_matrixElements[2][0] = a.X*a.X+a.Y*a.Y

    g_matrixElements[0][1] = b.X
    g_matrixElements[1][1] = b.Y
    g_matrixElements[2][1] = b.X*b.X+b.Y*b.Y

    g_matrixElements[0][2] = c.X
    g_matrixElements[1][2] = c.Y
    g_matrixElements[2][2] = c.X*c.X+c.Y*c.Y

    g_matrixElements[0][3] = d.X
    g_matrixElements[1][3] = d.Y
    g_matrixElements[2][3] = d.X*d.X+d.Y*d.Y

    return v.Fast4x4Determinant(&g_matrixElements) <= 0.0
}

func (delaunay *Delaunay) flipEdge(e he.EdgeIndex) {

    e1  := delaunay.Edges[e].ENext
    e2  := delaunay.Edges[e1].ENext
    eT  := delaunay.Edges[e].ETwin
    eT1 := delaunay.Edges[eT].ENext
    eT2 := delaunay.Edges[eT1].ENext

    // Check face reference to e or eT
    if delaunay.Faces[delaunay.Edges[e].FFace].EEdge == e {
        delaunay.Faces[delaunay.Edges[e].FFace].EEdge = e2
    }
    if delaunay.Faces[delaunay.Edges[eT].FFace].EEdge == eT {
        delaunay.Faces[delaunay.Edges[eT].FFace].EEdge = eT2
    }

    delaunay.Edges[e2].ENext  = eT1
    delaunay.Edges[eT1].EPrev = e2

    delaunay.Edges[eT2].ENext = e1
    delaunay.Edges[e1].EPrev  = eT2

    // Instead of creating new edges, we can just substitute the to-be-deleted edges
    delaunay.Edges[e].VOrigin = delaunay.Edges[eT2].VOrigin
    delaunay.Edges[e].ENext   = e2
    delaunay.Edges[e].EPrev   = eT1

    delaunay.Edges[eT].VOrigin = delaunay.Edges[e2].VOrigin
    delaunay.Edges[eT].ENext   = eT2
    delaunay.Edges[eT].EPrev   = e1

    // Reconnect other edges correctly
    delaunay.Edges[e1].ENext  = eT
    delaunay.Edges[eT2].EPrev = eT
    delaunay.Edges[eT1].ENext = e
    delaunay.Edges[e2].EPrev  = e

    // Faces rotate 90° counter clockwise. So have to be redefined for some edges.
    delaunay.Edges[e1].FFace  = delaunay.Edges[eT].FFace
    delaunay.Edges[eT1].FFace = delaunay.Edges[e].FFace


}

// Recursively flips edges until all triangles are correct delaunay triangles.
// This should never affect the frontier, as only internal edges can be flipped.
func (delaunay *Delaunay) legalizeTriangle(t CheckTriangle, firstCheck CheckTriangle) {

    // We might recursively itereate the same triangles. So if we arive at the very first triangle again,
    // we exit the recursion! Things are fine.
    if t.e == firstCheck.e && v.Equal(t.v, firstCheck.v) {
        return
    }

    // recursion ends at outer bounds.
    e1 := delaunay.Edges[t.e].ENext
    if e1 == he.EmptyEdge {
        return
    }
    e2 := delaunay.Edges[e1].ENext
    if e2 == he.EmptyEdge {
        return
    }

    a := delaunay.Vertices[delaunay.Edges[t.e].VOrigin].Pos
    b := delaunay.Vertices[delaunay.Edges[e1].VOrigin].Pos
    c := delaunay.Vertices[delaunay.Edges[e2].VOrigin].Pos

    if !triangleIsValid(a,b,c,t.v) {

        // Flip edge here
        delaunay.flipEdge(t.e)

        // Recursively check neighboring triangles
        delaunay.legalizeTriangle(CheckTriangle{delaunay.Edges[e1].ETwin, t.v}, firstCheck)
        delaunay.legalizeTriangle(CheckTriangle{delaunay.Edges[e2].ETwin, t.v}, firstCheck)
        delaunay.legalizeTriangle(CheckTriangle{delaunay.Edges[delaunay.Edges[t.e].EPrev].ETwin, c}, firstCheck)
        delaunay.legalizeTriangle(CheckTriangle{delaunay.Edges[delaunay.Edges[e1].EPrev].ETwin, c}, firstCheck)

        return
    }

}

// Recursive function that walks right as long as there is a triangle to be created (to fill holes)
// that match the pi/2 angle criteria.
// TODO: Pass vector position directly without going through the frontier and delaunay every time.
// Returns, how many consecutive triangles were created to the right.
func createConsecutiveTrianglesRight(frontier *t.Tree23, delaunay *Delaunay, baseVertex he.VertexIndex, leafNode t.TreeNodeIndex, maxFillAngle float64) {

    nextFrontier,_ := frontier.Next(leafNode)

    leafNodeValue := frontier.GetValue(leafNode).(FrontElement)
    nextFrontierValue := frontier.GetValue(nextFrontier).(FrontElement)

    basePos := delaunay.Vertices[baseVertex].Pos

    e0 := leafNodeValue.EdgeIndex
    e1 := nextFrontierValue.EdgeIndex

    v0 := delaunay.Vertices[delaunay.Edges[e0].VOrigin].Pos
    v1 := delaunay.Vertices[delaunay.Edges[e1].VOrigin].Pos

    isRight := v.IsRight2D(basePos, v0, v1)
    angle   := v.Angle(v.Sub(basePos, v0), v.Sub(v1, v0))

    if isRight && angle < maxFillAngle {

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
        eNew1 := delaunay.createEdge(baseVertex, he.EmptyEdge, e0, e1, delaunay.Edges[e0].FFace, v.Edge{})
        eNew2 := delaunay.createEdge(delaunay.Edges[e1].VOrigin, eNew1, he.EmptyEdge, he.EmptyEdge, he.EmptyFace, v.Edge{})

        f := delaunay.createFace(v.Vector{}, eNew2)
        delaunay.Edges[eNew1].ETwin = eNew2
        delaunay.Edges[eNew2].FFace = f

        delaunay.Edges[e0].ENext = eNew1
        delaunay.Edges[e1].EPrev = eNew1

        delaunay.Edges[e0].EPrev = e1
        delaunay.Edges[e1].ENext = e0

        frontier.ChangeValue(nextFrontier, FrontElement{eNew2, nextFrontierValue.PolarAngle, nextFrontierValue.Radius})

        //fmt.Printf("%v --> %v\n", leafNode, leafNodeValue)
        frontier.Delete(leafNodeValue)

        // Validate towards the previous triangle created consecutively
        d := delaunay.Vertices[delaunay.Edges[delaunay.Edges[delaunay.Edges[e0].ETwin].EPrev].VOrigin].Pos
        delaunay.legalizeTriangle(CheckTriangle{e0, d}, CheckTriangle{e0, d})
        // Validate towards inside the triangulation
        d = delaunay.Vertices[delaunay.Edges[delaunay.Edges[delaunay.Edges[e1].ETwin].EPrev].VOrigin].Pos
        delaunay.legalizeTriangle(CheckTriangle{e1, d}, CheckTriangle{e1, d})

        createConsecutiveTrianglesRight(frontier, delaunay, baseVertex, nextFrontier, maxFillAngle)

    }

}

// Similar to createConsecutiveTrianglesRight, this function checks to the right side, detects basins and
// removes them! This function should be called after createConsecutiveTrianglesRight so we already have
// the basic triangles in place before this optimization.
func removeBasinRight(frontier *t.Tree23, delaunay *Delaunay, leafNode t.TreeNodeIndex) {

    // Variable names are according to the paper!
    //vi := (*frontier)[fIndex].EdgeIndex

    //fVRPlus := fIndex+1
    //if fVRPlus >= len(*frontier) {
    //    fVRPlus = fVRPlus % len(*frontier)
    //}

    fVRPlus,_ := frontier.Next(leafNode)

    leafNodeValue := frontier.GetValue(leafNode).(FrontElement)
    vRPlus := frontier.GetValue(fVRPlus).(FrontElement)


    //vRPlus := (*frontier)[fVRPlus].EdgeIndex
    vRPlusVertex := delaunay.Edges[vRPlus.EdgeIndex].VOrigin

    r2 := vRPlus.Radius
    dR := leafNodeValue.Radius - r2
    d0 := vRPlus.PolarAngle - leafNodeValue.PolarAngle

    if vRPlus.PolarAngle < leafNodeValue.PolarAngle {
        d0 = vRPlus.PolarAngle+360 - leafNodeValue.PolarAngle
    }

    // basin detected!
    if dR / (r2*d0) > 2 {
        //fmt.Printf("Basin to the right detected!\n")

        foundMinimum := false
        fMin := fVRPlus
        vMin,_ := frontier.GetValue(fMin).(FrontElement)
        // Find local minimum radius
        for {
            fNext,_ := frontier.Next(fMin)
            vNext := frontier.GetValue(fNext).(FrontElement)
            //fNext := fMin+1
            //if fNext >= len(*frontier) {
            //    fNext = 0
            //}
            if vNext.Radius < vMin.Radius {
                fMin = fNext
                foundMinimum = true
            } else {
                break
            }
        }
        // No vertex after fVRPlus has a smaller radius than fVRPlus itself.
        if !foundMinimum {
            return
        }

        fCurrent := fMin
        // Find local maximum on the other side of the basin depending on
        // positive triangle area or vertex position (which side!)
        for {
            fLast := fCurrent
            //fCurrent = fCurrent+1
            //if fCurrent >= len(*frontier) {
            //    fCurrent = 0
            //}
            fCurrent,_ = frontier.Next(fCurrent)

            vLast := frontier.GetValue(fLast).(FrontElement)
            vCurrent := frontier.GetValue(fCurrent).(FrontElement)


            v0 := delaunay.Vertices[vRPlusVertex].Pos
            v1 := delaunay.Vertices[delaunay.Edges[vLast.EdgeIndex].VOrigin].Pos
            vNext := delaunay.Vertices[delaunay.Edges[vCurrent.EdgeIndex].VOrigin].Pos

            if v.IsRight2D(v0, v1, vNext) {
                //fmt.Printf("Basin filling triangle created.\n")

                // re-set fLast to the newly created outer edge!
                // don't forget to remove the edge(s) from the frontier!


                lastEdge    := vLast.EdgeIndex
                currentEdge := vCurrent.EdgeIndex

                eNew1 := delaunay.createEdge(vRPlusVertex, he.EmptyEdge, lastEdge, currentEdge, delaunay.Edges[lastEdge].FFace, v.Edge{})
                eNew2 := delaunay.createEdge(delaunay.Edges[currentEdge].VOrigin, eNew1, he.EmptyEdge, he.EmptyEdge, he.EmptyFace, v.Edge{})

                face  := delaunay.createFace(v.Vector{}, eNew2)
                delaunay.Edges[eNew1].ETwin = eNew2
                delaunay.Edges[eNew2].FFace = face

                delaunay.Edges[lastEdge].ENext = eNew1
                delaunay.Edges[currentEdge].EPrev = eNew1

                delaunay.Edges[lastEdge].EPrev = currentEdge
                delaunay.Edges[currentEdge].ENext = lastEdge

                //(*frontier)[fCurrent].EdgeIndex = eNew2
                frontier.ChangeValue(fCurrent, FrontElement{eNew2, vCurrent.PolarAngle, vCurrent.Radius})
                frontier.Delete(vLast)

                //(*frontier) = append((*frontier)[:fLast], (*frontier)[fLast+1:]...)
                //(*frontier) = append((*frontier)[:fLast], (*frontier)[fCurrent:]...)


                // It should stay the same. But in some circumstances, we delete one element underneath it. So the
                // list gets shorter and the frontier position shifts.
                //newFCurrent := fCurrent
                //if fLast < fCurrent {
                //    newFCurrent -= 1
                //    if newFCurrent < 0 {
                //        newFCurrent = newFCurrent+len(*frontier)
                //    }
                //}

                //fCurrent = newFCurrent

                d := delaunay.Vertices[delaunay.Edges[delaunay.Edges[delaunay.Edges[lastEdge].ETwin].EPrev].VOrigin].Pos
                delaunay.legalizeTriangle(CheckTriangle{lastEdge, d}, CheckTriangle{lastEdge, d})
                d = delaunay.Vertices[delaunay.Edges[delaunay.Edges[delaunay.Edges[currentEdge].ETwin].EPrev].VOrigin].Pos
                delaunay.legalizeTriangle(CheckTriangle{currentEdge, d}, CheckTriangle{currentEdge, d})

            } else {
                break
            }

        }




    }


}

// Recursive function that walks right as long as there is a triangle to be created (to fill holes)
// that match the pi/2 angle criteria.
// TODO: Pass vector position directly without going through the frontier and delaunay every time.
func createConsecutiveTrianglesLeft(frontier *t.Tree23, delaunay *Delaunay, baseVertex he.VertexIndex, leafNode t.TreeNodeIndex, maxFillAngle float64) {

    leafNodeValue := frontier.GetValue(leafNode).(FrontElement)

    prevFrontier,_ := frontier.Previous(leafNode)
    prevFrontierValue := frontier.GetValue(prevFrontier).(FrontElement)

    basePos := delaunay.Vertices[baseVertex].Pos

    e0 := leafNodeValue.EdgeIndex
    e1 := prevFrontierValue.EdgeIndex

    v0 := delaunay.Vertices[delaunay.Edges[e1].VOrigin].Pos
    v1Vertex := delaunay.Edges[delaunay.Edges[e1].ETwin].VOrigin
    v1 := delaunay.Vertices[v1Vertex].Pos

    isLeft := v.IsLeft2D(basePos, v0, v1)
    angle   := v.Angle(v.Sub(basePos, v0), v.Sub(v1, v0))

    if isLeft && angle < maxFillAngle {

        eNew1 := delaunay.createEdge(v1Vertex, he.EmptyEdge, e1, e0, delaunay.Edges[e0].FFace, v.Edge{})
        eNew2 := delaunay.createEdge(baseVertex, eNew1, he.EmptyEdge, he.EmptyEdge, he.EmptyFace, v.Edge{})

        f := delaunay.createFace(v.Vector{}, eNew2)
        delaunay.Edges[eNew1].ETwin = eNew2
        delaunay.Edges[eNew2].FFace = f

        delaunay.Edges[e1].ENext = eNew1
        delaunay.Edges[e0].EPrev = eNew1

        delaunay.Edges[e0].ENext = e1
        delaunay.Edges[e1].EPrev = e0

        frontier.ChangeValue(leafNode, FrontElement{eNew2, leafNodeValue.PolarAngle, leafNodeValue.Radius})

        frontier.Delete(prevFrontierValue)

        // Validate towards the previous triangle created consecutively
        d := delaunay.Vertices[delaunay.Edges[delaunay.Edges[delaunay.Edges[e0].ETwin].EPrev].VOrigin].Pos
        delaunay.legalizeTriangle(CheckTriangle{e0, d}, CheckTriangle{e0, d})
        // Validate towards inside the triangulation
        d = delaunay.Vertices[delaunay.Edges[delaunay.Edges[delaunay.Edges[e1].ETwin].EPrev].VOrigin].Pos
        delaunay.legalizeTriangle(CheckTriangle{e1, d}, CheckTriangle{e1, d})

        createConsecutiveTrianglesLeft(frontier, delaunay, baseVertex, leafNode, maxFillAngle)
    }
}

// Similar to createConsecutiveTrianglesRight, this function checks to the right side, detects basins and
// removes them! This function should be called after createConsecutiveTrianglesRight so we already have
// the basic triangles in place before this optimization.
func removeBasinLeft(frontier *t.Tree23, delaunay *Delaunay, leafNode t.TreeNodeIndex) {

    // Variable names are according to the paper!
    //vi := (*frontier)[fIndex].EdgeIndex

    //fVRPlus := fIndex-1
    //if fVRPlus < 0 {
    //    fVRPlus = fVRPlus+len(*frontier)
    //}

    fVRPlus,_ := frontier.Previous(leafNode)

    leafNodeValue := frontier.GetValue(leafNode).(FrontElement)
    vRPlus := frontier.GetValue(fVRPlus).(FrontElement)
    vRPlusVertex := delaunay.Edges[vRPlus.EdgeIndex].VOrigin




    //vRPlus := (*frontier)[fVRPlus].EdgeIndex
    //vRPlusVertex := delaunay.Edges[vRPlus].VOrigin

    r2 := vRPlus.Radius
    dR := leafNodeValue.Radius - r2
    d0 := leafNodeValue.PolarAngle - vRPlus.PolarAngle

    if vRPlus.PolarAngle > leafNodeValue.PolarAngle {
        d0 = leafNodeValue.PolarAngle+360 - vRPlus.PolarAngle
    }

    // basin detected!
    if dR / (r2*d0) > 2 {
        //fmt.Printf("Basin to the left detected!\n")

        foundMinimum := false
        fMin := fVRPlus
        vMin := frontier.GetValue(fMin).(FrontElement)
        // Find local minimum radius
        for {

            fNext,_ := frontier.Previous(fMin)
            vNext := frontier.GetValue(fNext).(FrontElement)

            //fNext := fMin-1
            //if fNext < 0 {
            //    fNext = len(*frontier)-1
            //}
            if vNext.Radius < vMin.Radius {
                fMin = fNext
                foundMinimum = true
            } else {
                break
            }
        }
        // No vertex after fVRPlus has a smaller radius than fVRPlus itself.
        if !foundMinimum {
            return
        }

        fCurrent := fMin
        // Find local maximum on the other side of the basin depending on
        // positive triangle area or vertex position (which side!)
        for {
            fLast := fCurrent
            vLast := frontier.GetValue(fLast).(FrontElement)
            //fCurrent = fCurrent-1

            fCurrent,_ = frontier.Previous(fCurrent)
            vCurrentE := frontier.GetValue(fCurrent).(FrontElement)

            //if fCurrent < 0 {
            //    fCurrent = len(*frontier)-1
            //}

            v0 := delaunay.Vertices[vRPlusVertex].Pos
            v1 := delaunay.Vertices[delaunay.Edges[vLast.EdgeIndex].VOrigin].Pos
            vCurrent := delaunay.Edges[vCurrentE.EdgeIndex].VOrigin

            if v.IsLeft2D(v0, v1, delaunay.Vertices[vCurrent].Pos) {
                //fmt.Printf("Basin filling triangle created.\n")

                lastEdge    := vLast.EdgeIndex
                //currentEdge := (*frontier)[fCurrent].EdgeIndex



                eNew1 := delaunay.createEdge(vCurrent, he.EmptyEdge, lastEdge, vRPlus.EdgeIndex, delaunay.Edges[lastEdge].FFace, v.Edge{})
                eNew2 := delaunay.createEdge(vRPlusVertex, eNew1, he.EmptyEdge, he.EmptyEdge, he.EmptyFace, v.Edge{})

                face  := delaunay.createFace(v.Vector{}, eNew2)
                delaunay.Edges[eNew1].ETwin = eNew2
                delaunay.Edges[eNew2].FFace = face

                delaunay.Edges[lastEdge].ENext = eNew1
                delaunay.Edges[vRPlus.EdgeIndex].EPrev = eNew1

                delaunay.Edges[lastEdge].EPrev = vRPlus.EdgeIndex
                delaunay.Edges[vRPlus.EdgeIndex].ENext = lastEdge

                //(*frontier)[fVRPlus].EdgeIndex = eNew2
                //(*frontier) = append((*frontier)[:fLast], (*frontier)[fVRPlus:]...)
                //frontier.removeFrontierAt(fLast)

                frontier.ChangeValue(fVRPlus, FrontElement{eNew2, vRPlus.PolarAngle, vRPlus.Radius})
                frontier.Delete(vLast)

                // It should stay the same. But in some circumstances, we delete one element underneath it. So the
                // list gets shorter and the frontier position shifts.
                //newFCurrent := fCurrent
                //if fLast < fCurrent {
                //    newFCurrent -= 1
                //    if newFCurrent < 0 {
                //        newFCurrent = newFCurrent+len(*frontier)
                //    }
                //}
                //
                //fCurrent = newFCurrent

                d := delaunay.Vertices[delaunay.Edges[delaunay.Edges[delaunay.Edges[lastEdge].ETwin].EPrev].VOrigin].Pos
                delaunay.legalizeTriangle(CheckTriangle{lastEdge, d}, CheckTriangle{lastEdge, d})
                d = delaunay.Vertices[delaunay.Edges[delaunay.Edges[delaunay.Edges[vRPlus.EdgeIndex].ETwin].EPrev].VOrigin].Pos
                delaunay.legalizeTriangle(CheckTriangle{vRPlus.EdgeIndex, d}, CheckTriangle{vRPlus.EdgeIndex, d})

            } else {
                break
            }

        }




    }


}

func extendByPoint(frontier *t.Tree23, p DelaunayPoint, delaunay *Delaunay, center v.Vector) {

    frontierItem, err := frontier.FindFirstLargerLeaf(p.PolarAngle)
    // In this case, p.PolarAngle is after the very last node
    if err != nil {
        frontierItem,_ = frontier.GetSmallestLeaf()
    }

    frontierItemValue := frontier.GetValue(frontierItem).(FrontElement)

    vi := delaunay.createVertex(p.Point)

    existingE := frontierItemValue.EdgeIndex

    fi1 := delaunay.Edges[existingE].FFace
    twinV := delaunay.Edges[delaunay.Edges[existingE].ETwin].VOrigin
    ei1 := delaunay.createEdge(twinV, he.EmptyEdge, existingE, he.EmptyEdge, fi1, v.Edge{})
    ei2 := delaunay.createEdge(vi, ei1, he.EmptyEdge, he.EmptyEdge, he.EmptyFace, v.Edge{})
    fi2 := delaunay.createFace(v.Vector{}, ei2)
    delaunay.Edges[ei2].FFace = fi2
    delaunay.Edges[ei1].ETwin = ei2

    delaunay.Edges[existingE].ENext = ei1

    ej1 := delaunay.createEdge(vi, he.EmptyEdge, ei1, existingE, fi1, v.Edge{})
    ej2 := delaunay.createEdge(delaunay.Edges[existingE].VOrigin, ej1, he.EmptyEdge, he.EmptyEdge, he.EmptyFace, v.Edge{})
    fj2 := delaunay.createFace(v.Vector{}, ej2)

    delaunay.Edges[ej2].FFace = fj2
    delaunay.Edges[ej1].ETwin = ej2

    delaunay.Edges[ei1].ENext = ej1

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


    frontier.ChangeValue(frontierItem, FrontElement{ej2, frontierItemValue.PolarAngle, frontierItemValue.Radius})
    fmt.Printf("Insert0: %v\n", FrontElement{ei2, p.PolarAngle, p.Distance})
    frontier.Insert(FrontElement{ei2, p.PolarAngle, p.Distance})

    maxFillAngle := 125.0

    // There must exist a triangle behind the previous frontier edge. So this should be save.
    d := delaunay.Vertices[delaunay.Edges[delaunay.Edges[delaunay.Edges[existingE].ETwin].EPrev].VOrigin].Pos
    delaunay.legalizeTriangle(CheckTriangle{existingE, d}, CheckTriangle{existingE, d})


    previousLeaf,err := frontier.Previous(frontierItem)

    createConsecutiveTrianglesRight(frontier, delaunay, vi, frontierItem, maxFillAngle)
    //removeBasinRight(frontier, delaunay, previousLeaf)



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

    createConsecutiveTrianglesLeft(frontier, delaunay, vi, previousLeaf, maxFillAngle)
    //removeBasinLeft(frontier, delaunay, previousLeaf)

}

func (delaunay *Delaunay)triangulatePoints(pl *DelaunayPointList, frontier *t.Tree23) {

    for _,p := range pl.Points {
        extendByPoint(frontier, p, delaunay, (*pl).Origin)
    }

}

func Triangulate(pointList v.PointList) Delaunay {

    dPointList := preparePointList(&pointList)
    //defer fmt.Printf("point to be inserted last: %v\n", dPointList.Points[52])
    //dPointList.Points = dPointList.Points[:95]

    // Any planar triangulation: total degree == 3f + k = 2e
    // With k == points on convex hull of all points. Let k = p
    // https://en.wikipedia.org/wiki/Planar_graph#Euler.27s_formula
    // Eulers formula states: v - e + f = 2
    // For a finite, sparse, connected graph, we have:
    // e <= 3v - 6
    edgeCount := (3*len(dPointList.Points)-6)*2
    delaunay := Delaunay {
        Vertices:           make([]he.HEVertex, len(dPointList.Points)),
        FirstFreeVertexPos: 0,
        Edges:              make([]he.HEEdge, edgeCount + 1000),
        FirstFreeEdgePos:   0,
        Faces:              make([]he.HEFace, edgeCount/3 + len(dPointList.Points)+1000),
        FirstFreeFacePos:   0,
    }

    frontier := delaunay.initializeTriangulation(&dPointList)

    delaunay.triangulatePoints(&dPointList, frontier)

    return delaunay

}
