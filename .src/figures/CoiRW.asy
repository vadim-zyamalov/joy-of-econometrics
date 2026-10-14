texpreamble("\usepackage{fontspec}\usepackage{unicode-math}\setmainfont{Noto Sans}\setmathfont{Noto Sans Math}");

settings.tex = "xelatex";
settings.outformat = "pdf";
settings.prc = false;
settings.render = 0;

unitsize(1mm);

import graph;
import stats;

int n = 300;
real mlt = 3, sgm = 1, xmax = 60, a0 = -10, a1 = .5;

real cdfGaussian(real x) {
	static real b0 = 0.2316419,
				b1 = 0.319381530,
				b2 = -0.356563782,
				b3 = 1.781477937,
				b4 = -1.821255978,
				b5 = 1.330274429;
	real t = 1 / (1 + b0 * abs(x));
	real poly = b1 * t + b2 * t^2 + b3 * t^3 + b4 * t^4 + b5 * t^5;
	return (x <= 0) ? Gaussian(x) * poly : 1 - Gaussian(x) * poly;
}

real[][] gendata(int n, real mlt = 1.0) {
	real[] x, y;
	x[0] = 0;
	y[0] = 0;
	for (int i=1; i <= n; ++i) {
		x[i] = Gaussrand() + x[i-1];
		y[i] = x[i] + Gaussrand();
	}
	x.delete(0);
	y.delete(0);
    return new real[][] {mlt * x, mlt * y};
}

real[][] data = gendata(n, mlt);
real[] x = data[0], y = data[1];


xaxis("$t$", xmin=0, arrow=Arrow(TeXHead));
yaxis("$y$", arrow=Arrow(TeXHead));

draw(graph(sequence(n), x), black);
shipout("RW");

draw(graph(sequence(n), y), red);
draw(graph(sequence(n), y - x), gray(.5));
