import graph;
import graph3;

// Рисуем эллипс, растягивая и поворачивая единичную окружность.
path ellipse(pair c, real a, real b, real rot=0, real sc=1) {
    return(shift(c)*rotate(rot)*scale(sc*a, sc*b)*unitcircle);
}

using path_real = path(real);
real findr(path_real fn, path target, real start=1, real tol=0.000001, int iter=100) {
    real lb=0, ub=start*2, r = start;
    real err=100;
    int i=0;

    path candidate = fn(ub);
    pair[] isects = intersectionpoints(target, candidate);
    if (isects.length == 0) abort("Wrong UB");

    while (err > tol) {
        ++i;
        if (i > iter) {
            write("Failed to converge!");
            return r;
        }
        r = (ub + lb) / 2;
        path candidate = fn(r);
        pair[] isects = intersectionpoints(target, candidate);
        if (isects.length == 0) {
            lb = r;
        } else if (isects.length == 2) {
            pair dp = isects[0] - isects[1];
            err = abs(dp);
            ub = r;
        } else {
            return r;
        }
    }

    return r;
}

// Функция распределения для N(0, 1).
real cdfGaussian(real x) { return 0.5 * (1 + erf(x)); }

real max0(real x) {
    return max(x, 0);
}

bool isclose(real a, real b, real relTol = 1e-6, real absTol = 1e-7) {
	return abs(a - b) <= max(relTol * max(abs(a), abs(b)), absTol);
}

bool isclose(triple a, triple b, real relTol = 1e-6, real absTol = 1e-7) {
	bool rx = isclose(a.x, b.x, relTol, absTol),
	     ry = isclose(a.y, b.y, relTol, absTol),
	     rz = isclose(a.z, b.z, relTol, absTol);
	return rx && ry && rz;
}

// Символ перпендикулярности.
// https://sourceforge.net/p/asymptote/discussion/409349/thread/5f2a11af/
void perpmark3d(triple ptI, triple ptO, triple ptJ, real size, pen p=black) {
   triple imI = ptO + size * unit(ptI - ptO),
          imJ = ptO + size * unit(ptJ - ptO),
          imK = imI + imJ - ptO;
   draw(imI--imK--imJ, p);
}

// Расширение оси, заданной ортом, если координата выходит за пределы.
void incaxis(triple vec, real lmb, real delta, pen p = dotted+black) {
    if (lmb == 0) return;
    if (lmb > 1 || lmb < 0) draw(O--((lmb + sgn(lmb) * delta) * vec), dotted);
}

void incaxis(triple vec, real beta, real lmb, real delta, pen p = dotted+black) {
    if (beta * lmb > 0) {
        real lmb = sgn(beta) * max(abs(lmb), abs(beta));
        incaxis(vec, lmb, delta, p);
    } else if (beta * lmb < 0) {
        incaxis(vec, lmb, delta, p);
        incaxis(vec, beta, delta, p);
    } else if (beta == 0 && lmb != 0) {
        incaxis(vec, lmb, delta, p);
    } else if (lmb == 0 && beta != 0) {
        incaxis(vec, beta, delta, p);
    } else abort("Wrong values of lmb and beta!");
}

// Координата проекции ветора на орт.
real perp2orth(triple vec, triple orth) {
    real lv = length(vec), lo = length(orth);
    real cosA = dot(vec, orth) / (lv * lo);
    return lv * cosA;
}

//
triple lintersection2(triple v1, triple p1, triple v2, triple p2) {
    triple vec1 = cross(p2 - p1, v2),
           vec2 = cross(p2 - p1, v1),
           vecV = cross(v1, v2);
    real vecV1 = dot(vecV, vec1),
    	 vecV2 = dot(vecV, vec2),
         vecVV = dot(vecV, vecV);
    if (vecVV == 0) {return (infinity, infinity, infinity);}
    real x1 = vecV1 / vecVV, x2 = vecV2 / vecVV;
    triple i1 = p1 + x1 * v1, i2 = p2 + x2 * v2;
    return isclose(i1, i2) ? i1 : (infinity, infinity, infinity);
}
