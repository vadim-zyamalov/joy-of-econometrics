texpreamble("\usepackage{fontspec}\usepackage{unicode-math}\setmainfont{Noto Sans}\setmathfont{Noto Sans Math}");

settings.tex = "xelatex";
settings.outformat = "pdf";
settings.prc = false;
settings.render = 0;

unitsize(5mm);

import graph;
import stats;

scale(Linear, Linear(3));

int n = 25;
real lm1 = .5, lm2 = .7;
real lm3 = .5;
pair lm4 = (.5, .7);

real f1(real t) { return lm1^t + lm2^t; }
real f2(real t) { return lm3^t + t * lm3^t; }
real f3(real t) { return length(lm4)^t * cos(angle(lm4) * t); }

xaxis("$t$", xmin=0, arrow=Arrow(TeXHead));
yaxis("$y$", arrow=Arrow(TeXHead));

draw(graph(f1, 0, n, operator ..), black);
draw(graph(f2, 0, n, operator ..), blue);
draw(graph(f3, 0, n, operator ..), green);
