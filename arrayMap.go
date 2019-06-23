package mtSweepCircle

import (
	"fmt"
	"math"
)

type List struct {
	Next  *List
	Prev  *List
	Value FrontElement
}

// SkipList is the actual skiplist representation.
// It saves all nodes accessible from the start and end and keeps track of element count, eps and levels.
type ArrayMap struct {
	lookup []*List
	list   *List
}

// New returns a new empty, initialized Skiplist.
func NewArrayMap(size int) *ArrayMap {
	//fmt.Printf("New ArrayMap --> %v\n", size)
	return &ArrayMap{
		lookup: make([]*List, size, size),
		list:   nil,
	}
}

// IsEmpty checks, if the skiplist is empty.
func (t *ArrayMap) IsEmpty() bool {
	return t.list == nil
}

// FindGreaterOrEqual finds the first element, that is greater or equal to the given FrontElement e.
// The comparison is done on the keys (So on ExtractKey()).
// FindGreaterOrEqual runs in approx. O(log(n))
func (t *ArrayMap) FindGreaterOrEqual(key float64) *List {

	if t.list == nil {
		return nil
	}

	//fmt.Printf("    find %.3f\n", key)

	i := int(key * 0.159154943 * float64(len(t.lookup)))
	n := t.lookup[i]

	//fmt.Printf("    initial i: %v\n", i)
	//fmt.Printf("    list: %v\n", t.list)

	i0 := i
	// Get the first valid reference to the linked list!
	i++
	for n == nil && i != i0 {
		if i >= len(t.lookup) {
			i = 0
		}
		n = t.lookup[i]
		i++

	}
	if n == nil {
		n = t.list
	}

	//fmt.Printf("    initial: %.3f\n", n.Value.PolarAngle)

	// Default.
	if n.Prev.Value.PolarAngle < key && n.Value.PolarAngle > key {
		return n
	}

	// Move along the linked list until we actually have the first greater-or-equal element!
	if n.Value.PolarAngle < key {
		//fmt.Println(t)
		//fmt.Println("Smaller.")
		//fmt.Printf()
		for n.Value.PolarAngle < key {
			// For when we are all around the list and no element matches, the one we look for is just bigger than the biggest one in here.
			if n == t.list.Prev {
				return n.Next
			}

			//if key > t.list.Prev.Value.PolarAngle {
			//return t.list
			//}

			n = n.Next
		}
	} else {
		for n.Prev.Value.PolarAngle >= key {
			if n == t.list {
				return n
			}
			n = n.Prev
		}
	}

	return n
}

// Delete removes an element equal to e from the skiplist, if there is one.
// If there are multiple entries with the same value, Delete will remove one of them
// (Which one will change based on the actual skiplist layout)
// Delete runs in approx. O(log(n))
func (t *ArrayMap) Delete(n *List) {
	epsilon := math.Nextafter(1, 2) - 1
	n.Prev.Next = n.Next
	n.Next.Prev = n.Prev

	if n == t.list {
		t.list = n.Next
	}

	// Remove element from the lookup table if we are sure, it really is the right one.
	key := int(n.Value.PolarAngle * 0.159154943 * float64(len(t.lookup)))
	n2 := t.lookup[key]
	if n2 != nil && math.Abs(n.Value.PolarAngle-n2.Value.PolarAngle) <= epsilon {
		t.lookup[key] = nil
	}

	// Create new link, if the next array position is empty. This way, we try to sustain as many references as we had before
	if key < len(t.lookup)-1 && t.lookup[key+1] == nil {
		t.lookup[int(n.Next.Value.PolarAngle*0.159154943*float64(len(t.lookup)))] = n.Next
	}
	if key > 0 && t.lookup[key-1] == nil {
		t.lookup[int(n.Prev.Value.PolarAngle*0.159154943*float64(len(t.lookup)))] = n.Prev
	}

	n.Next = nil
	n.Prev = nil
}

// Insert inserts the given FrontElement into the skiplist.
// Insert runs in approx. O(log(n))
func (t *ArrayMap) InsertAfter(e FrontElement, prev *List) *List {

	if prev != nil {
		//fmt.Printf("    insert %.3f after %.3f\n", e.PolarAngle, prev.Value.PolarAngle)
	}

	n := &List{nil, nil, e}

	// Very first one
	if prev == nil {
		n.Prev = n
		n.Next = n
		t.list = n
	} else {
		n.Prev = prev
		n.Next = prev.Next
		prev.Next = n
		n.Next.Prev = n

		// New first element
		if e.PolarAngle < prev.Value.PolarAngle {
			t.list = n
		}
	}

	//fmt.Printf("%v\n", t.list)
	//fmt.Printf("Angle: %.2f, len: %v\n", e.PolarAngle, len(t.lookup))
	//fmt.Printf("    Insert %.2f at index %v\n", e.PolarAngle, int(e.PolarAngle * 0.159154943 * float64(len(t.lookup))))

	// 1 / (2*pi)
	t.lookup[int(e.PolarAngle*0.159154943*float64(len(t.lookup)))] = n

	return n
}

// GetSmallestNode returns the very first/smallest node in the skiplist.
// GetSmallestNode runs in O(1)
func (t *ArrayMap) GetSmallestNode() *List {
	if t.list == nil {
		return nil
	}

	return t.list
}

// GetLargestNode returns the very last/largest node in the skiplist.
// GetLargestNode runs in O(1)
func (t *ArrayMap) GetLargestNode() *List {
	if t.list == nil {
		return nil
	}
	return t.list.Prev
}

func (t *ArrayMap) String() string {
	s := ""
	first := t.list
	current := first.Next

	if first != nil {
		s += fmt.Sprintf("%.3f --> ", first.Value.PolarAngle)
	}

	for current != first {
		s += fmt.Sprintf("%.3f --> ", current.Value.PolarAngle)
		current = current.Next
	}

	s += "\n"
	for i := 0; i < len(t.lookup); i++ {
		v := "     "
		if t.lookup[i] != nil {
			v = "<ptr>"
		}
		s += fmt.Sprintf("%v ... ", v)
	}

	return s
}
