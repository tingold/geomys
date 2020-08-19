package geomys

import (
	"github.com/reconditematter/mym"
	"math"
)

const (
	DistAndoyer  = iota // Andoyer's approximate distance
	DistEllipse         // great ellipse distance
	DistGeodesic        // geodesic distance
)

// GeoMatrix -- computes an n-by-n symmetric matrix of pairwise distances
// between the points p[0],...,p[n-1]. The distances are computed on the
// spheroid `sph` using a predefined method specified by `dist`
// (DistAndoyer,DistEllipse,DistGeodesic).
func GeoMatrix(sph Spheroid, p []Point, dist int) mym.Sym0 {
	n := len(p)
	M := mym.NewSym0(n)
	//
	switch dist {
	case DistAndoyer:
		for i, pi := range p {
			for j := i + 1; j < n; j++ {
				pj := p[j]
				geodist := Andoyer(sph, pi, pj)
				M.Set(i, j, geodist)
			}
		}
	case DistEllipse:
		grell := NewGreatEllipse(sph)
		for i, pi := range p {
			for j := i + 1; j < n; j++ {
				pj := p[j]
				geodist, _, _ := grell.Inverse(pi, pj)
				M.Set(i, j, geodist)
			}
		}
	case DistGeodesic:
		panic("geomys.GeoMatrix: not implemented: `DistGeodesic`")
	default:
		panic("geomys.Geomatrix: domain error: `dist`")
	}
	//
	return M
}

// PrjMatrix -- computes an n-by-n symmetric matrix of pairwise distances
// between the points p[0],...,p[n-1]. The distances are computed in the
// plane using the map projection transformation defined by `prj`.
func PrjMatrix(prj MapProjection, p []Point) mym.Sym0 {
	n := len(p)
	M := mym.NewSym0(n)
	//
	for i, pi := range p {
		ci := prj.Project(pi)
		for j := i + 1; j < n; j++ {
			pj := p[j]
			cj := prj.Project(pj)
			dx := ci[0] - cj[0]
			dy := ci[1] - cj[1]
			prjdist := math.Hypot(dx, dy)
			M.Set(i, j, prjdist)
		}
	}
	//
	return M
}
