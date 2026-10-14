texpreamble("\usepackage{fontspec}\usepackage{unicode-math}\setmainfont{Noto Sans}\setmathfont{Noto Sans Math}");

settings.tex = "xelatex";
settings.outformat = "pdf";
settings.prc = false;
settings.render = 0;

import graph;
import stats;

real xmin=0, xmax=12;
real se=.5;

real pm(real se) { return se * 3; }

real F(real x) { return 2 + .1 * x; }

using realf = pair(real);
realf dens(real dx, real se) {
	return new pair(real x) {return (Gaussian(x, se) + dx, F(dx) + x);};
}

void sceleton(real x, real se = .5, real pm = 1.5, picture pic = currentpicture) {
	draw((x, F(x) - pm)--(x, F(x) + pm), dashed+blue, pic=pic);
	draw((x, F(x))--(x + Gaussian(0, se), F(x)), dashed+blue, pic=pic);
}

picture pic0, pic1;
unitsize(1cm, pic=pic0);
unitsize(1cm, pic=pic1);

path f0 = graph(F, 0, 12);
draw(pic = pic0, f0);
draw(pic = pic1, f0);

for (int i = 2; i <= 10; i += 2) {
	real se2 = se + (i - 2) / 20;
	realf innerf0 = dens(i, se);
	realf innerf1 = dens(i, se2);
	draw(pic = pic0, graph(innerf0, -pm(se), pm(se)), blue);
	draw(pic = pic1, graph(innerf1, -pm(se2), pm(se2)), blue);
	sceleton(i, pic = pic0);
	sceleton(i, pm = pm(se2), se = se2, pic = pic1);
}

xaxis("$x$", xmin=xmin-0.2, xmax=xmax+0.2, arrow=Arrow(TeXHead), pic = pic0);
xaxis("$x$", xmin=xmin-0.2, xmax=xmax+0.2, arrow=Arrow(TeXHead), pic = pic1);
yaxis("$y$", arrow=Arrow(TeXHead), pic = pic0);
yaxis("$y$", arrow=Arrow(TeXHead), pic = pic1);

shipout("OLS_homo_errors", pic = pic0);
shipout("OLS_hetero_errors", pic = pic1);
