package geomys

import (
	"math"
	"strings"
)

// GeoHash -- computes the geohash of length `n` of the given
// geographic latitude and longitude stored in `p`. Also computes
// the angular resolution of the geohash measured in degrees.
//
// See: https://en.wikipedia.org/wiki/Geohash
//
// The given length `n` is assumed to be n∈{5,7,9,11,13,15}.
// When n<5, it is set to 5; when n>15, it is set to 15.
// When n is an even number, it is set to the next odd number.
func GeoHash(n int, p Point) (hash string, res float64) {
	if n < 5 {
		n = 5
	}
	if n > 15 {
		n = 15
	}
	if n&1 == 0 {
		n++
	}
	//
	const gh = "0123456789bcdefghjkmnpqrstuvwxyz"
	const two45 = 1 << 45
	const lateps = 90.0 / two45
	const loneps = 180.0 / two45
	//
	abs := func(b bool) uint64 {
		if b {
			return 1
		}
		return 0
	}
	//
	lat, lon := p.Geo()
	if lat == 90 {
		lat -= lateps / 2
	}
	if lon == 180 {
		lon = -180
	}
	//
	ulat := uint64(math.Floor(lat/lateps) + two45)
	ulon := uint64(math.Floor(lon/loneps) + two45)
	//
	var (
		b uint64
		B strings.Builder
	)
	for i := 0; i < 5*n; i++ {
		if i&1 == 0 {
			b = (b << 1) + abs((ulon&two45) != 0)
			ulon <<= 1
		} else {
			b = (b << 1) + abs((ulat&two45) != 0)
			ulat <<= 1
		}
		if (i+1)%5 == 0 {
			B.WriteByte(gh[int(b)])
			b = 0
		}
	}
	//
	return B.String(), math.Ldexp(180, -(5 * n / 2))
}
