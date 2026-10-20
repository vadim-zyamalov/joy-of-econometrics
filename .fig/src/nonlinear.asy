texpreamble("\usepackage{fontspec}\usepackage{unicode-math}\setmainfont{Noto Sans}\setmathfont{Noto Sans Math}");

settings.tex = "xelatex";
settings.outformat = "pdf";
settings.prc = false;
settings.render = 0;

unitsize(1cm);

import graph;
import stats;
import myfunc;

int n = 300;
real sgm = 0.25, xmax = 10, a0 = 1, a1 = 0.7, a2 = -.05;

real[][] gendata(int n) {
	real[] x, y, y0;
	for (int i=0; i < n; ++i) x[i] = unitrand() * xmax;
	x = sort(x);
    for (int i=0; i < n; ++i) y[i] = a0 + a1 * x[i] + a2 * x[i]^2 + Gaussrand()*sgm;
    return new real[][] {x, y};
}

real[][] data = gendata(n);
real[] x = data[0], y = data[1];

xaxis("$x$", xmin=0, arrow=Arrow(TeXHead));
yaxis("$y$", ymax = 1.1 * max(y), arrow=Arrow(TeXHead));

for (int i: sequence(n)) {
    filldraw(Circle((x[i], y[i]), .025), white, black);
}

linefit reg = leastsquares(x, y);
draw((min(x), reg.fit(min(x)))--(max(x), reg.fit(max(x))), blue);

file fout = output("nonlinear.txt");
write(fout, transpose(data));
close(fout);

shipout("nonlinear");
