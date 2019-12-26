package geomys

import (
	"github.com/reconditematter/mym"
	"math"
)

// Geocentric -- a coordinate converter between
// the geographic and geocentric coordinate systems.
type Geocentric struct {
	s Spheroid
}

// NewGeocentric -- returns a geographic/geocentric
// coordinate converter for the spheroid `s`.
func NewGeocentric(s Spheroid) Geocentric {
	return Geocentric{s}
}

// Spheroid -- returns the spheroid of the converter `geocen`.
func (geocen Geocentric) Spheroid() Spheroid {
	return geocen.s
}

// Forward -- converts the geographic coordinates of `p`
// to the geocentric coordinates `xyz`.
func (geocen Geocentric) Forward(p Point) (xyz [3]float64) {
	a := geocen.s.A()
	e2 := geocen.s.E2()
	φ, λ := p.Geo()
	sinφ, cosφ := mym.SinCosD(φ)
	sinλ, cosλ := mym.SinCosD(λ)
	N := a / math.Sqrt(1-e2*sinφ*sinφ)
	xyz[0] = N * cosφ * cosλ
	xyz[1] = N * cosφ * sinλ
	xyz[2] = N * (1 - e2) * sinφ
	return
}
