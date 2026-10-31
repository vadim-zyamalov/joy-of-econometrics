texpreamble("\usepackage{fontspec}\usepackage{unicode-math}\setmainfont{Noto Sans}\setmathfont{Noto Sans Math}");

settings.tex = "xelatex";
settings.outformat = "pdf";
settings.prc = false;
settings.render = 0;

unitsize(1cm);

import graph;
import stats;
import myfunc;

file fin = input("../nonlinear.txt").line();
real[][] data = fin;
data = transpose(data);
real[] x = data[0], y = data[1];
close(fin);

file fin = input("../nonlinear_betas.txt").line();
real[] betas = fin;
real b1 = betas[0], b2 = betas[1], b0 = betas[2];
close(fin);

xaxis("$x$", xmin=0, arrow=Arrow(TeXHead));
yaxis("$y$", ymax = 1.1 * max(y), arrow=Arrow(TeXHead));

for (int i: sequence(x.length)) {
    filldraw(Circle((x[i], y[i]), .025), white, black);
}

linefit reg = leastsquares(x, y);
draw((min(x), reg.fit(min(x)))--(max(x), reg.fit(max(x))), blue);

real[] ynl;
for (int i : sequence(x.length)) ynl[i] = b0 + b1 * x[i] + b2 * x[i]^2;

draw(graph(x, ynl), red);
