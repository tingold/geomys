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
	φ, λ, _ := p.Geo()
	sinφ, cosφ := mym.SinCosD(φ)
	sinλ, cosλ := mym.SinCosD(λ)
	N := a / math.Sqrt(1-e2*sinφ*sinφ)
	xyz[0] = N * cosφ * cosλ
	xyz[1] = N * cosφ * sinλ
	xyz[2] = N * (1 - e2) * sinφ
	return
}

// Inverse -- converts the geocentric coordinates `xyz`
// to the pair of corresponding geographic coordinates.
// This function causes a runtime panic when any of `xyz` ∉ [-10²³,10²³].
//
// Reference: Fukushima, T. Transformation from Cartesian to Geodetic Coordinates
// Accelerated by Halley's Method. J Geodesy 79, 689–693 (2006).
//
// DOI: https://doi.org/10.1007/s00190-006-0023-2
func (geocen Geocentric) Inverse(xyz [3]float64) Point {
	if !(math.Abs(xyz[0]) <= 1e23 && math.Abs(xyz[1]) <= 1e23 && math.Abs(xyz[2]) <= 1e23) {
		panic("geomys.Geocentric.Inverse: domain error: `xyz`")
	}
	a := geocen.s.A()
	e2 := geocen.s.E2()
	fc := 1 - geocen.s.F()
	//
	P := math.Hypot(xyz[0], xyz[1]) / a
	Z := fc * math.Abs(xyz[2]) / a
	//
	S0 := Z
	C0 := fc * P
	A0 := math.Hypot(S0, C0)
	D0 := Z*mym.Cb(A0) + e2*mym.Cb(S0)
	F0 := P*mym.Cb(A0) - e2*mym.Cb(C0)
	B0 := 1.5 * P * mym.Sq(e2*S0*C0) * (A0 - fc)
	//
	S1 := D0*F0 - B0*S0
	C1 := F0*F0 - B0*C0
	A1 := math.Hypot(S1, C1)
	D1 := Z*mym.Cb(A1) + e2*mym.Cb(S1)
	F1 := P*mym.Cb(A1) - e2*mym.Cb(C1)
	B1 := 1.5 * e2 * S1 * C1 * ((P*S1-Z*C1)*A1 - e2*S1*C1)
	//
	S2 := D1*F1 - B1*S1
	C2 := F1*F1 - B1*C1
	//
	φ := math.Atan2(S2, fc*C2)
	if xyz[2] < 0 {
		φ = -φ
	}
	λ := math.Atan2(xyz[1], xyz[0])

	//todo add altitude
	alt := 0.0

	return Geo(φ*(180/math.Pi), λ*(180/math.Pi), alt)
}
