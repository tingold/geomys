package geomys

import (
	"github.com/reconditematter/mym"
	"math"
)

func hat(y, x float64) (s, c float64) {
	norm := math.Hypot(y, x)
	if norm < mym.Epsilon {
		return 0, 1
	}
	return y / norm, x / norm
}

func polyval(n int, p []float64, s int, x float64) float64 {
	var y float64
	if n >= 0 {
		y = p[s]
		s++
		n--
	}
	for ; n >= 0; n-- {
		y = y*x + p[s]
		s++
	}
	return y
}

func ellA1m1f(eps float64) float64 {
	const n = 8
	const m = n / 2
	coeff := [...]float64{25, 64, 256, 4096, 0, 16384}
	t := polyval(m, coeff[:], 0, eps*eps) / coeff[m+1]
	return (t + eps) / (1 - eps)
}

func ellC1f(eps float64) (C [9]float64) {
	const n = 8
	coeff := [...]float64{19, -64, 384, -1024, 2048, 7, -18, 128, -256, 4096, -9, 72, -128, 6144, -11, 96, -160, 16384, 35, -56, 10240, 9, -14, 4096, -33, 14336, -429, 262144}
	eps2, d := eps*eps, eps
	oo := 0
	for L := 1; L <= n; L++ {
		m := (n - L) / 2
		C[L] = d * polyval(m, coeff[:], oo, eps2) / coeff[oo+m+1]
		oo += m + 2
		d *= eps
	}
	return
}

func ellC1pf(eps float64) (C [9]float64) {
	const n = 8
	coeff := [...]float64{-4879, 9840, -20736, 36864, 73728, -86171, 120150, -142080, 115200, 368640, 8703, -7200, 3712, 12288, 1082857, -688608, 258720, 737280, -141115, 41604, 92160, -2200311, 533134, 860160, 459485, 516096, 109167851, 82575360}
	eps2, d := eps*eps, eps
	oo := 0
	for L := 1; L <= n; L++ {
		m := (n - L) / 2
		C[L] = d * polyval(m, coeff[:], oo, eps2) / coeff[oo+m+1]
		oo += m + 2
		d *= eps
	}
	return
}

func ellSinSeries(sin, cos float64, c []float64) float64 {
	k := len(c)
	n := k - 1
	ar := 2 * (cos - sin) * (cos + sin)
	var y0, y1 float64
	if (n & 1) != 0 {
		k--
		y0 = c[k]
	}
	n /= 2
	for ; n != 0; n-- {
		k--
		y1 = ar*y0 - y1 + c[k]
		k--
		y0 = ar*y1 - y0 + c[k]
	}
	return 2 * sin * cos * y0
}
