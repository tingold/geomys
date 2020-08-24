package geomys

import (
	"github.com/reconditematter/mym"
	"math"
)

// Albers -- Albers conical equal-area map projection.
type Albers struct {
	sph          Spheroid
	par          map[string]float64
	n, c, ρ0, λ0 float64
}

// NewAlbers -- returns a new Albers map projection based on the spheroid `sph`.
//
// The projection has the following parameters:
//	lat1 -- latitude of the 1st standard parallel
//	lat2 -- latitude of the 2nd standard parallel
//	lat0 -- latitude of the center
//	lon0 -- longitude of the center
func NewAlbers(sph Spheroid, lat1, lat2, lat0, lon0 float64) Albers {
	par := map[string]float64{"lat1": lat1, "lat2": lat2, "lat0": lat0, "lon0": lon0}
	alb := Albers{sph: sph, par: par}
	alb.inialb()
	return alb
}

// Spheroid -- returns the spheroid of the map projection.
func (prj Albers) Spheroid() Spheroid {
	if prj.par == nil {
		panic("geomys.Albers.Spheroid: uninitialized structure")
	}
	//
	return prj.sph
}

// Params -- returns the parameters of the map projection.
func (prj Albers) Params() map[string]float64 {
	if prj.par == nil {
		panic("geomys.Albers.Params: uninitialized structure")
	}
	//
	par := make(map[string]float64)
	for k, v := range prj.par {
		par[k] = v
	}
	return par
}

// Project -- transforms a geographic point from
// the spheroid into a location on the plane.
func (prj Albers) Project(p Point) (xy [2]float64) {
	if prj.par == nil {
		panic("geomys.Albers.Project: uninitialized structure")
	}
	//
	lat, lon := p.Geo()
	sinφ, _ := mym.SinCosD(lat)
	q := prj.qalb(sinφ)
	ρ := prj.sph.A() * math.Sqrt(prj.c-prj.n*q) / prj.n
	θ := prj.n * (lon - prj.λ0)
	sinθ, cosθ := mym.SinCosD(θ)
	xy[0] = ρ * sinθ
	xy[1] = prj.ρ0 - ρ*cosθ
	return
}

func (prj *Albers) inialb() {
	lat1 := prj.par["lat1"]
	lat2 := prj.par["lat2"]
	lat0 := prj.par["lat0"]
	lon0 := prj.par["lon0"]
	//
	sinφ1, cosφ1 := mym.SinCosD(lat1)
	sinφ2, cosφ2 := mym.SinCosD(lat2)
	sinφ0, _ := mym.SinCosD(lat0)
	//
	m1 := prj.malb(sinφ1, cosφ1)
	m2 := prj.malb(sinφ2, cosφ2)
	//
	q1 := prj.qalb(sinφ1)
	q2 := prj.qalb(sinφ2)
	q0 := prj.qalb(sinφ0)
	//
	n := (m1 - m2) * (m1 + m2) / (q2 - q1)
	c := m1*m1 + n*q1
	ρ0 := prj.sph.A() * math.Sqrt(c-n*q0) / n
	//
	prj.n, prj.c, prj.ρ0, prj.λ0 = n, c, ρ0, lon0
}

func (prj Albers) malb(sinφ, cosφ float64) float64 {
	return cosφ / math.Sqrt(1-prj.sph.E2()*sinφ*sinφ)
}

func (prj Albers) qalb(sinφ float64) float64 {
	e2 := prj.sph.E2()
	e := math.Sqrt(e2)
	return (1 - e2) * (sinφ/(1-e2*sinφ*sinφ) - (1/(2*e))*math.Log((1-e*sinφ)/(1+e*sinφ)))
}
