package geomys

// Point -- represents a pair of geographic coordinates
// (latitude,longitude).
type Point struct {
	lat, lon, alt float64
}

// Geo -- returns a point with the geographic latitude `lat`
// and the geographic longitude `lon`.
// This function causes a runtime panic when lat∉[-90,90]
// or lon∉[-180,180].
func Geo(lat, lon, alt float64) Point {
	if !(-90 <= lat && lat <= 90) {
		panic("geomys.Geo: lat not in [-90,90]")
	}
	if !(-180 <= lon && lon <= 180) {
		panic("geomys.Geo: lon not in [-180,180]")
	}
	return Point{lat, lon, alt}
}

// Geo -- returns the geographic coordinates of `p`.
func (p Point) Geo() (lat, lon, alt float64) {
	lat, lon = p.lat, p.lon
	return
}
