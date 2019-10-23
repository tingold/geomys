package geomys

import (
	"github.com/reconditematter/mym"
	"math"
)

// Spheroid -- represents the figure of the Earth as an oblate ellipsoid of revolution.
type Spheroid struct {
	a, f float64
}

// Clarke1866 -- returns the spheroid used in the North American Datum 1927.
func Clarke1866() Spheroid {
	const a = 6378206.4
	const b = 6356583.8
	return Spheroid{a, (a - b) / a}
}

// International1924 -- returns the spheroid adopted by IUGG in 1924 (Madrid).
func International1924() Spheroid {
	return Spheroid{6378388, 1.0 / 297.0}
}

// WGS1972 -- returns the World Geodetic System 1972 spheroid.
func WGS1972() Spheroid {
	return Spheroid{6378135, 1.0 / 298.26}
}

// GRS1967 -- returns the spheroid adopted by IUGG in 1967 (Lucerne).
func GRS1967() Spheroid {
	return Spheroid{6378160, 1.0 / 298.247167427}
}

// GRS1980 -- returns the spheroid adopted by IUGG in 1979 (Canberra).
func GRS1980() Spheroid {
	return Spheroid{6378137, 1.0 / 298.257222101}
}

// WGS1984 -- returns the World Geodetic System 1984 spheroid.
func WGS1984() Spheroid {
	return Spheroid{6378137, 1.0 / 298.257223563}
}

// IERS2003 -- returns the IERS Technical Note No.32 spheroid.
func IERS2003() Spheroid {
	return Spheroid{6378136.6, 1.0 / 298.25642}
}

// SRMmax -- returns a spheroid with an exaggerated flattening (f=1/150)
// for the purposes of algorithm testing.
func SRMmax() Spheroid {
	return Spheroid{6400000, 1.0 / 150}
}

// NewSpheroid -- returns a spheroid with the equatorial (major) axis `a`
// and the (first) flattening `f`.
// This function causes a runtime panic when a∉[1,10²²] or f∉[0,1/150].
func NewSpheroid(a, f float64) Spheroid {
	if !(1 <= a && a <= 1e22) {
		panic("geomys.NewSpheroid: domain error: `a`")
	}
	if !(0 <= f && f <= 1.0/150.0) {
		panic("geomys.NewSpheroid: domain error: `f`")
	}
	if f < mym.SqrtEps {
		f = 0
	}
	return Spheroid{a, f}
}

// NewSphere -- returns a sphere with the radius `r`.
// This function causes a runtime panic when r∉[1,10²²].
func NewSphere(r float64) Spheroid {
	if !(1 <= r && r <= 1e22) {
		panic("geomys.NewSphere: domain error: `r`")
	}
	return Spheroid{r, 0}
}

// A -- returns the equatorial (major) axis (a) of `s`.
func (s Spheroid) A() float64 {
	if s.a > 1 {
		return s.a
	}
	return 1
}

// B -- returns the polar (minor) axis (b) of `s`.
func (s Spheroid) B() float64 {
	return s.A() * (1 - s.F())
}

// F -- returns the (first) flattening (f) of `s`:
//
//    f = (a-b)/a.
func (s Spheroid) F() float64 {
	return s.f
}

// Fp -- returns the second flattening (f') of `s`:
//
//    f' = (a-b)/b.
func (s Spheroid) Fp() float64 {
	return s.F() / (1 - s.F())
}

// Fpp -- returns the third flattening (f") of `s`:
//
//    f" = (a-b)/(a+b).
func (s Spheroid) Fpp() float64 {
	return s.F() / (2 - s.F())
}

// E2 -- returns the (first) eccentricity squared (e²) of `s`:
//
//    e² = (a²-b²)/a².
func (s Spheroid) E2() float64 {
	return s.F() * (2 - s.F())
}

// Ep2 -- returns the second eccentricity squared (e'²) of `s`:
//
//    e'² = (a²-b²)/b².
func (s Spheroid) Ep2() float64 {
	return s.F() * (2 - s.F()) / mym.Sq(1-s.F())
}

// Epp2 -- returns the third eccentricity squared (e"²) of `s`:
//
//    e"² = (a²-b²)/(a²+b²).
func (s Spheroid) Epp2() float64 {
	return s.F() * (2 - s.F()) / (1 + mym.Sq(1-s.F()))
}

// Rm -- returns the radius of a sphere having equal meridian length to that of `s`.
func (s Spheroid) Rm() float64 {
	// using Ramanujan's celebrated formula for the perimeter of an ellipse
	f := s.F()
	t := 3 * mym.Sq(f/(2-f))
	return s.A() * (1 - f/2) * (1 + t/(10+math.Sqrt(4-t)))
}

// Rs -- returns the radius of a sphere having equal surface area to that of `s`.
func (s Spheroid) Rs() float64 {
	f := s.F()
	e2 := f * (2 - f)
	series := 1 + e2*(1.0/3+e2*(1.0/5+e2*(1.0/7+e2*(1.0/9+e2*(1.0/11+e2*(1.0/13+e2*(1.0/15)))))))
	return s.A() * math.Sqrt((1+(1-e2)*series)/2)
}
