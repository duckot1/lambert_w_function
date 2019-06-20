package lambert

import (
	"math"
)

var gslDblEpsilon = 2.2204460492503131e-16
var oneOverE = 1 / math.E

type result struct {
	Val     float64
	err     float64
	iters   int
	success bool
}

func halleyInteration(x, wInitial float64, maxIters int) result {

	w := wInitial
	r := result{}

	var i int

	for i = 1; i <= maxIters; i++ {

		var tol float64

		e := math.Exp(w)
		p := w + 1.0
		t := w*e - x

		if w > 0 {
			/* Newton iteration */
			t = (t / p) / e
		} else {
			/* Halley iteration */
			t = t / (e*p - 0.5*(p+1.0)*t/p)
		}

		w = w - t

		tol = gslDblEpsilon * math.Max(math.Abs(w), 1.0/(math.Abs(p)*e))

		if math.Abs(t) < tol {
			r.Val = w
			r.err = 2.0 * tol
			r.iters = i
			r.success = true
			return r
		}
	}

	/* should never get here */
	r.Val = w
	r.err = math.Abs(w)
	r.iters = i
	r.success = true
	return r
}

func seriesEval(r float64) float64 {
	var c [12]float64
	c[0] = -1.0
	c[1] = 2.331643981597124203363536062168
	c[2] = -1.812187885639363490240191647568
	c[3] = 1.936631114492359755363277457668
	c[4] = -2.353551201881614516821543561516
	c[5] = 3.066858901050631912893148922704
	c[6] = -4.175335600258177138854984177460
	c[7] = 5.858023729874774148815053846119
	c[8] = -8.401032217523977370984161688514
	c[9] = 12.250753501314460424
	c[10] = -18.100697012472442755
	c[11] = 27.029044799010561650

	t8 := c[8] + r*(c[9]+r*(c[10]+r*c[11]))
	t5 := c[5] + r*(c[6]+r*(c[7]+r*t8))
	t1 := c[1] + r*(c[2]+r*(c[3]+r*(c[4]+r*t5)))
	return c[0] + r*t1
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

func GslSfLambertW0E(x float64) result {

	q := x + oneOverE

	r := result{}

	if x == 0.0 {
		r.Val = 0.0
		r.err = 0.0
		r.success = true
		return r
	} else if q < 0.0 {
		/* Strictly speaking this is an error. But because of the
		 * arithmetic operation connecting x and q, I am a little
		 * lenient in case of some epsilon overshoot. The following
		 * answer is quite accurate in that case. Anyway, we have
		 * to return GSL_EDOM.
		 */
		r.Val = -1.0
		r.err = math.Sqrt(-q)
		r.success = false // GSL_EDOM
		return r
	} else if q < 1.0e-03 {
		/* series near -1/E in sqrt(q) */
		root := math.Sqrt(q)
		r.Val = seriesEval(root)
		r.err = 2.0 * gslDblEpsilon * math.Abs(r.Val)
		r.success = true
		return r
	}

	maxIters := 100
	var w float64

	if x < 1.0 {
		/* obtain initial approximation from series near x=0;
		 * no need for extra care, since the Halley iteration
		 * converges nicely on this branch
		 */
		p := math.Sqrt(2.0 * math.E * q)
		w = -1.0 + p*(1.0+p*(-1.0/3.0+p*11.0/72.0))
	} else {
		/* obtain initial approximation from rough asymptotic */
		w = math.Log(x)
		if x > 3.0 {
			w = w - math.Log(w)
		}
	}
	return halleyInteration(x, w, maxIters)
}
