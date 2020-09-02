package geomys

import (
	"github.com/reconditematter/mym"
	"math"
)

// GreatEllipse -- great ellipse solver for a spheroidal model of the Earth.
type GreatEllipse struct {
	sph Spheroid
}

// NewGreatEllipse -- returns a great ellipse solver for the spheroid `sph`.
func NewGreatEllipse(sph Spheroid) GreatEllipse {
	return GreatEllipse{sph}
}

// Spheroid -- returns the spheroid of `g`.
func (g GreatEllipse) Spheroid() Spheroid {
	return g.sph
}

// Inverse -- solves the inverse problem: given two points `p1` and `p2`, find
// the great ellipse distance `s12` (meters) between the points, also find the azimuths
// `α1` (degrees) at `p1` and `α2` (degrees) at `p2`.
func (g GreatEllipse) Inverse(p1, p2 Point) (s12 float64, α1, α2 float64) {
	a, f := g.sph.A(), g.sph.F()
	f1, e2 := 1-f, g.sph.E2()
	//
	lat1, lon1, _ := p1.Geo()
	if lat1 > 90*(1-mym.Epsilon) {
		lat1 = 90 * (1 - mym.Epsilon)
	}
	if lat1 < -90*(1-mym.Epsilon) {
		lat1 = -90 * (1 - mym.Epsilon)
	}
	sβ1, cβ1 := mym.SinCosD(lat1)
	sβ1 *= f1
	sβ1, cβ1 = hat(sβ1, cβ1)
	//
	lat2, lon2, _ := p2.Geo()
	if lat2 > 90*(1-mym.Epsilon) {
		lat2 = 90 * (1 - mym.Epsilon)
	}
	if lat2 < -90*(1-mym.Epsilon) {
		lat2 = -90 * (1 - mym.Epsilon)
	}
	sβ2, cβ2 := mym.SinCosD(lat2)
	sβ2 *= f1
	sβ2, cβ2 = hat(sβ2, cβ2)
	//
	λ12 := lon2 - lon1
	if λ12 <= -180 {
		λ12 += 2 * 180
	} else if λ12 > 180 {
		λ12 -= 2 * 180
	}
	sλ12, cλ12 := mym.SinCosD(λ12)
	//
	sγ1 := cβ2 * sλ12
	cγ1 := +cβ1*sβ2 - sβ1*cβ2*cλ12
	//
	sγ2 := cβ1 * sλ12
	cγ2 := -sβ1*cβ2 + cβ1*sβ2*cλ12
	//
	sσ12 := math.Hypot(sγ1, cγ1)
	cσ12 := sβ1*sβ2 + cβ1*cβ2*cλ12
	//
	sγ1, cγ1 = hat(sγ1, cγ1)
	sγ2, cγ2 = hat(sγ2, cγ2)
	cγ0 := math.Hypot(cγ1, sγ1*sβ1)
	//
	sσ1, cσ1 := sβ1, cβ1*cγ1
	if sσ1 == 0 && cσ1 == 0 {
		cσ1 = 1
	}
	sσ1, cσ1 = hat(sσ1, cσ1)
	//
	sσ2 := sσ1*cσ12 + cσ1*sσ12
	cσ2 := cσ1*cσ12 - sσ1*sσ12
	//
	k2 := f * (2 - f) * cγ0 * cγ0
	eps := k2 / (2*(1+math.Sqrt(1-k2)) - k2)
	A1 := a * (1 + ellA1m1f(eps)) * (1 - eps) / (1 + eps)
	C1 := ellC1f(eps)
	//
	s12 = A1 * (math.Atan2(sσ12, cσ12) + (ellSinSeries(sσ2, cσ2, C1[:]) - ellSinSeries(sσ1, cσ1, C1[:])))
	//
	α1 = math.Atan2(sγ1, cγ1*math.Sqrt(1-e2*cβ1*cβ1)) * (180 / math.Pi)
	α2 = math.Atan2(sγ2, cγ2*math.Sqrt(1-e2*cβ2*cβ2)) * (180 / math.Pi)
	return
}

// Direct -- solves the direct problem: given the source point `p1`, the azimuth `α1` (degrees),
// and the great ellipse distance `s12` (meters), find the target point `p2`, also find the azimuth
// `α2` (degrees) at `p2`.
func (g GreatEllipse) Direct(p1 Point, α1 float64, s12 float64) (p2 Point, α2 float64) {
	a, f := g.sph.A(), g.sph.F()
	f1, e2 := 1-f, g.sph.E2()
	//
	lat1, lon1,_ := p1.Geo()
	sγ1, cγ1 := mym.SinCosD(α1)
	sβ1, cβ1 := mym.SinCosD(lat1)
	sβ1 *= f1
	sβ1, cβ1 = hat(sβ1, cβ1)
	sγ1, cγ1 = hat(sγ1*math.Sqrt(1-e2*cβ1*cβ1), cγ1)
	//
	sγ0 := sγ1 * cβ1
	cγ0 := math.Hypot(cγ1, sγ1*sβ1)
	sσ1 := sβ1
	sλ1 := sγ0 * sβ1
	cσ1 := cβ1 * cγ1
	if sβ1 == 0 && cγ1 == 0 {
		cσ1 = 1
	}
	cλ1 := cσ1
	sσ1, cσ1 = hat(sσ1, cσ1)
	//
	k2 := e2 * cγ0 * cγ0
	eps := k2 / (2*(1+math.Sqrt(1-k2)) - k2)
	A1 := a * (1 + ellA1m1f(eps)) * (1 - eps) / (1 + eps)
	C1a := ellC1f(eps)
	B11 := ellSinSeries(sσ1, cσ1, C1a[:])
	s, c := math.Sincos(B11)
	sτ1 := sσ1*c + cσ1*s
	cτ1 := cσ1*c - sσ1*s
	//
	C1pa := ellC1pf(eps)
	τ12 := s12 / A1
	s, c = math.Sincos(τ12)
	B12 := -ellSinSeries(sτ1*c+cτ1*s, cτ1*c-sτ1*s, C1pa[:])
	σ12 := τ12 - (B12 - B11)
	sσ12, cσ12 := math.Sincos(σ12)
	//
	sσ2 := sσ1*cσ12 + cσ1*sσ12
	cσ2 := cσ1*cσ12 - sσ1*sσ12
	sβ2 := cγ0 * sσ2
	cβ2 := math.Hypot(sγ0, cγ0*cσ2)
	//
	sλ2 := sγ0 * sσ2
	cλ2 := cσ2
	sγ2 := sγ0
	cγ2 := cγ0 * cσ2
	//
	lon12 := math.Atan2(sλ2*cλ1-cλ2*sλ1, cλ2*cλ1+sλ2*sλ1) * (180 / math.Pi)
	lon2 := lon1 + lon12
	if lon2 <= -180 {
		lon2 += 2 * 180
	} else if lon2 > 180 {
		lon2 -= 2 * 180
	}
	//
	lat2 := math.Atan2(sβ2, f1*cβ2) * (180 / math.Pi)
	p2 = Geo(lat2, lon2, 0.0)
	α2 = math.Atan2(sγ2, cγ2*math.Sqrt(1-e2*cβ2*cβ2)) * (180 / math.Pi)
	return
}
