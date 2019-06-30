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
	maxN = 2000000
)

// timeTrack will print out the number of nanoseconds since the start time divided by n
// Useful for printing out how long each iteration took in a benchmark
func timeTrack(start time.Time, n int, name string) {
	loopNS := time.Since(start).Nanoseconds() / int64(n)
	fmt.Printf("%s: %d\n", name, loopNS)
}

func TestInsertAndFind(t *testing.T) {

	epsilon := math.Nextafter(1, 2) - 1
	toKey := func(e int) float64 {
		return float64(e) / float64(maxN) * 2. * math.Pi
	}
	toKeyF := func(e float64) float64 {
		return e * 2. * math.Pi
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

	//fmt.Println(list)

	for i := 0; i < maxN; i++ {
		if e := list.FindGreaterOrEqual(toKey(i)); math.Abs(e.Value.PolarAngle-toKey(i)) > epsilon {
			fmt.Printf("Fail to find an element. %.3f != %.3f\n", e.Value.PolarAngle, toKey(i))
			t.Fail()
		}
	}

	//fmt.Println(list)
	r := rand.New(rand.NewSource(9))

	list = NewArrayMap(maxN)
	// Test at random positions in the list.
	rList := r.Perm(maxN)
	for i, e := range rList {
		n := list.FindGreaterOrEqual(toKey(e))
		var prev *List = nil
		if i != 0 {
			prev = n.Prev
		}

		list.InsertAfter(FrontElement{PolarAngle: toKey(e)}, prev)

		//fmt.Println(list)
	}

	for i, e := range rList {

		if v := list.FindGreaterOrEqual(toKey(e)); math.Abs(toKey(e)-v.Value.PolarAngle) >= epsilon {
			fmt.Printf("Fail to find an exact element at %v. %.3f != %.3f\n", i, toKey(e), v.Value.PolarAngle)
			t.Fail()
		}
	}

	//fmt.Println(list)

	//fmt.Printf("find larger: %.3f\n", list.FindGreaterOrEqual(6).Value.PolarAngle)

	sum := 0
	for i := 0; i < len(list.lookup); i++ {
		if list.lookup[i] != nil {
			sum++
		}
	}
	fmt.Printf("Lookup table: %v\n", sum)

	r = rand.New(rand.NewSource(99))

	for i := 0; i < maxN; i++ {
		e := r.Float64()
		n := list.FindGreaterOrEqual(toKeyF(e))
		prev := n.Prev

		//fmt.Printf("About to insert %.3f --> %.3f\n", e, toKeyF(e))
		inserted := list.InsertAfter(FrontElement{PolarAngle: toKeyF(e)}, prev)

		// Find the one we just inserted.
		if v := list.FindGreaterOrEqual(toKeyF(e)); math.Abs(toKeyF(e)-v.Value.PolarAngle) >= epsilon {
			fmt.Printf("Inserted element not found %v != %v\n", toKeyF(e), v.Value.PolarAngle)
			t.Fail()
		}

		//fmt.Printf("Inserted: %.3f\n", inserted.Value.PolarAngle)
		//fmt.Println(list)
		list.Delete(inserted)
		//fmt.Println(list)

		// Find the one we just inserted.
		if v := list.FindGreaterOrEqual(toKeyF(e)); math.Abs(toKeyF(e)-v.Value.PolarAngle) <= epsilon {
			fmt.Printf("Element not correctly removed exact element %v != %v\n", toKeyF(e), v.Value.PolarAngle)
			t.Fail()
		}
	}

	sum = 0
	for i := 0; i < len(list.lookup); i++ {
		if list.lookup[i] != nil {
			sum++
		}
	}
	fmt.Printf("Lookup table: %v\n", sum)

	//fmt.Println(list)

}
