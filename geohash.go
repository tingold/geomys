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
// The given length `n` is assumed to be n∈{3,5,7,9,11,13,15}.
// When n<3, it is set to 3; when n>15, it is set to 15.
// When n is an even number, it is set to the next odd number.
func GeoHash(n int, p Point) (hash string, res float64) {
	if n < 3 {
		n = 3
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

// HashGeo -- decodes the given geohash `hash` into a pair of geographic
// coordinates returned in `p`.
//
// See: https://en.wikipedia.org/wiki/Geohash
//
// When `hash` is decoded without errors, sets `ok` to true;
// otherwise sets `ok` to false and returns (0,0) in `p`.
// When the length of `hash` is greater than 15, only the first 15 characters
// are decoded; when the length is less than 3, ((0,0),false) is returned.
func HashGeo(hash string) (p Point, ok bool) {
	const gh = "0123456789bcdefghjkmnpqrstuvwxyz"
	const two45 = 1 << 45
	const lateps = 90.0 / two45
	const loneps = 180.0 / two45
	//
	hash = strings.ToLower(hash)
	n := len(hash)
	if n < 3 {
		return Point{}, false
	}
	if n > 15 {
		n = 15
	}
	//
	abs := func(b bool) uint64 {
		if b {
			return 1
		}
		return 0
	}
	//
	var ulon, ulat uint64
	for k, j := 0, 0; k < n; k++ {
		b := strings.IndexByte(gh, hash[k])
		if b < 0 {
			return Point{}, false
		}
		for m := 16; m > 0; m >>= 1 {
			if j == 0 {
				ulon = (ulon << 1) + abs((b&m) != 0)
			} else {
				ulat = (ulat << 1) + abs((b&m) != 0)
			}
			j ^= 1
		}
	}
	ulon <<= 1
	ulat <<= 1
	ulon += 1
	ulat += 1
	s := 5 * (18 - n)
	ulon <<= (s / 2)
	ulat <<= s - (s / 2)
	//
	lat, lon := float64(ulat)*lateps-90, float64(ulon)*loneps-180
	//
	var scale float64
	switch {
	case n == 3 || n == 4:
		scale = 1
	case n == 5 || n == 6:
		scale = 1e2
	case n == 7 || n == 8:
		scale = 1e3
	case n == 9 || n == 10:
		scale = 1e5
	case n == 11 || n == 12:
		scale = 1e6
	case n == 13 || n == 14:
		scale = 1e8
	default:
		scale = 1e9
	}
	lat = math.Round(scale*lat) / scale
	lon = math.Round(scale*lon) / scale
	return Geo(lat, lon), true
}
