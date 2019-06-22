package mtSweepCircle

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"
	//"github.com/pkg/profile"
)

const (
	maxN = 10
)

// timeTrack will print out the number of nanoseconds since the start time divided by n
// Useful for printing out how long each iteration took in a benchmark
func timeTrack(start time.Time, n int, name string) {
	loopNS := time.Since(start).Nanoseconds() / int64(n)
	fmt.Printf("%s: %d\n", name, loopNS)
}

func TestInsertAndFind(t *testing.T) {
	
	toKey := func(e int) float64 {
		return float64(e)/float64(maxN)*2.*math.Pi
	}
	
	list := NewArrayMap(int(math.Sqrt(float64(maxN))))

	if !list.IsEmpty() {
		t.Fail()
		fmt.Println("list not empty.")
	}

	// Test at the beginning of the list.
	var last *List = nil
	for i := 0; i < maxN; i++ {
		last = list.InsertAfter(FrontElement{PolarAngle: toKey(i)}, last)
	}
	
	fmt.Println(list)
	
	for i := 0; i < maxN; i++ {
		if e := list.FindGreaterOrEqual(toKey(i)); math.Abs(e.Value.PolarAngle - toKey(i)) > 0.00001 {
			fmt.Printf("Fail to find an element. %.3f != %.3f\n", e.Value.PolarAngle, toKey(i))
			t.Fail()
		}
	}
	
	//fmt.Println(list)
	
	

	list = NewArrayMap(int(math.Sqrt(float64(maxN))))
	// Test at random positions in the list.
	rList := rand.Perm(maxN)
	for _, e := range rList {
		n := list.FindGreaterOrEqual(toKey(e))
		var prev *List = nil
		if n != nil {
			prev = n.Prev
		}
		list.InsertAfter(FrontElement{PolarAngle: toKey(e)}, prev)
	}
		
	fmt.Println(list)
	
	for i, e := range rList {
		
		if v := list.FindGreaterOrEqual(toKey(e)); math.Abs(toKey(e) - v.Value.PolarAngle) >= 0.000001 {
			fmt.Printf("Fail to find an exact element at %v. %.3f != %.3f\n", i, toKey(e), v.Value.PolarAngle)
			t.Fail()
		}
	}
	
	//rList = rand.Perm(maxN)
	//for _, e := range rList {
	//	n := list.FindGreaterOrEqual(toKey(e))
	//	var prev *List = nil
	//	if n != nil {
	//		prev = n.Prev
	//	}
	//	inserted := list.InsertAfter(FrontElement{PolarAngle: toKey(e)}, prev)
	//	
	//	// Find the one we just inserted.
	//	if v := list.FindGreaterOrEqual(toKey(e)); math.Abs(toKey(e) - v.Value.PolarAngle) >= 0.000001 {
	//		fmt.Printf("Inserted element not found %v != %v\n", toKey(e), v.Value.PolarAngle)
	//		t.Fail()
	//	}
	//	
	//	list.Delete(inserted)
	//	
	//	// Find the one we just inserted.
	//	if v := list.FindGreaterOrEqual(toKey(e)); math.Abs(toKey(e) - v.Value.PolarAngle) <= 0.00001 {
	//		fmt.Printf("Element not correctly removed exact element %v != %v\n", toKey(e), v.Value.PolarAngle)
	//		t.Fail()
	//	}
	//}
	
}

//func TestPrev(t *testing.T) {
//	list := New()
//
//	for i := 0; i < maxN; i++ {
//		list.Insert(FrontElement(i))
//	}
//
//	smallest := list.GetSmallestNode()
//	largest := list.GetLargestNode()
//
//	lastNode := largest
//	node := lastNode
//	for node != smallest {
//		node = list.Prev(node)
//		// Must always be incrementing here!
//		if node.value.(FrontElement) >= lastNode.value.(FrontElement) {
//			t.Fail()
//		}
//		// Next.Prev must always point to itself!
//		if list.Prev(list.Next(node)) != node {
//			t.Fail()
//		}
//		lastNode = node
//	}
//
//	if list.Prev(smallest) != largest {
//		t.Fail()
//	}
//}
//
//func TestNext(t *testing.T) {
//	list := New()
//
//	for i := 0; i < maxN; i++ {
//		list.Insert(FrontElement(i))
//	}
//
//	smallest := list.GetSmallestNode()
//	largest := list.GetLargestNode()
//
//	lastNode := smallest
//	node := lastNode
//	for node != largest {
//		node = list.Next(node)
//		// Must always be incrementing here!
//		if node.value.(FrontElement) <= lastNode.value.(FrontElement) {
//			t.Fail()
//		}
//		// Next.Prev must always point to itself!
//		if list.Next(list.Prev(node)) != node {
//			t.Fail()
//		}
//		lastNode = node
//	}
//
//	if list.Next(largest) != smallest {
//		t.Fail()
//	}
//}

