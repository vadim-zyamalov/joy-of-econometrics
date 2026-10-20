texpreamble("\usepackage{fontspec}\usepackage{unicode-math}\setmainfont{Noto Sans}\setmathfont{Noto Sans Math}");

settings.tex = "xelatex";
settings.outformat = "pdf";
settings.prc = false;
settings.render = 0;

import graph;
import stats;

picture pPool, pFE;
unitsize(1cm, pic=pPool);
unitsize(1cm, pic=pFE);

int nt = 100;
int n = 4;
real a1 = .33;
pair[] dxPool = {(0, 0), (1, .33), (2, .66), (3, .99)};
pair[] dxFE = {(0, 0), (1, 1), (2, 2), (3, 3)};
pen[] pens = {orange, red, green, blue};

real[][][] gendata(int n, int nt, real sigma = .1) {
	real[][] x, y, w, z;
	for (int i : sequence(n)) {
		x[i] = new real[nt];
		y[i] = new real[nt];
		w[i] = new real[nt];
		z[i] = new real[nt];
		for (int t : sequence(nt)) {
			real xit = unitrand();
			real yit = a1 * xit + Gaussrand() * sigma;
    		x[i][t] = xit + dxPool[i].x;
    		y[i][t] = yit + dxPool[i].y;
    		w[i][t] = xit + dxFE[i].x;
    		z[i][t] = yit + dxFE[i].y;
    	}
	}
	return new real[][][] {x, y, w, z};
}

real[][][] data = gendata(n, nt);
real[][] x = data[0], y = data[1], w = data[2], z = data[3];

xaxis(pPool, "$x$", xmin=0, arrow=Arrow(TeXHead));
xaxis(pFE, "$x$", xmin=0, arrow=Arrow(TeXHead));
yaxis(pPool, "$y$", arrow=Arrow(TeXHead));
yaxis(pFE, "$y$", arrow=Arrow(TeXHead));

real[] allx, ally, allw, allz;
for (int i : sequence(n)) {
	allx = concat(allx, x[i]);
	ally = concat(ally, y[i]);
	allw = concat(allw, w[i]);
	allz = concat(allz, z[i]);
	for (int t : sequence(nt)) {
		real minx = min(x[i]) - .1, maxx = max(x[i]) + .1;
		real minw = min(w[i]) - .1, maxw = max(w[i]) + .1;

		// Pool
		filldraw(pPool, Circle((x[i][t], y[i][t]), .025), white, pens[i]);
		linefit reg = leastsquares(x[i], y[i]);
		draw(pPool, (minx, reg.fit(minx))--(maxx, reg.fit(maxx)), black+dotted);

		// FE
		filldraw(pFE, Circle((w[i][t], z[i][t]), .025), white, pens[i]);
		linefit reg = leastsquares(w[i], z[i]);
		draw(pFE, (minw, reg.fit(minw))--(maxw, reg.fit(maxw)), black+dotted);
	}
}

linefit reg = leastsquares(allx, ally);
draw(pPool, (min(allx), reg.fit(min(allx)))--(max(allx), reg.fit(max(allx))), gray(.5)+dotted);

linefit reg = leastsquares(allw, allz);
draw(pFE, (min(allw), reg.fit(min(allw)))--(max(allw), reg.fit(max(allw))), gray(.5)+dotted);

shipout("panel-pool", pic = pPool);
shipout("panel-FE", pic = pFE);
