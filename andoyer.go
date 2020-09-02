package geomys

import (
	"github.com/reconditematter/mym"
	"math"
)

// Andoyer -- computes an approximation to the geodesic distance
// between `p1` and `p2` on the spheroid `s`.
//
// See: Andoyer, H. Bull. Géodésique (1932) 34: 77. https://doi.org/10.1007/BF03030136
func Andoyer(s Spheroid, p1, p2 Point) float64 {
	φ1, λ1,_ := p1.Geo()
	φ2, λ2,_ := p2.Geo()
	//
	F := (φ1 + φ2) / 2
	G := (φ1 - φ2) / 2
	L := (λ1 - λ2) / 2
	//
	sinF, cosF := mym.SinCosD(F)
	sinG, cosG := mym.SinCosD(G)
	sinL, cosL := mym.SinCosD(L)
	//
	S := math.Hypot(sinG*cosL, cosF*sinL)
	C := math.Hypot(cosG*cosL, sinF*sinL)
	ω := math.Atan2(S, C)
	D := 2 * s.A() * ω
	//
	R := S * C / ω
	H1 := (3*R - 1) / (2 * C * C)
	H2 := (3*R + 1) / (2 * S * S)
	d := D * (1 + s.F()*(H1*mym.Sq(sinF*cosG)-H2*mym.Sq(cosF*sinG)))
	//
	if mym.FiniteIs(d) {
		return d
	}
	// special case: antipodal p1&p2
	if mym.FiniteIs(R) {
		return math.Pi * s.Rm()
	}
	// special case: p1=p2
	return 0
}
