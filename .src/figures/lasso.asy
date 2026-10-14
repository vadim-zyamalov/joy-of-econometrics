texpreamble("\usepackage{fontspec}\usepackage{unicode-math}\setmainfont{Noto Sans}\setmathfont{Noto Sans Math}");

settings.tex = "xelatex";
settings.outformat = "pdf";
settings.prc = false;
settings.render = 0;

import graph;
import myfunc;

using path_real = path(real);
path_real ell(pair c, real a, real b, real rot) {
	return new path(real sc) { return ellipse(c, a, b, rot, sc); };
}

picture pLasso, pRidge;

currentpicture=pLasso;

unitsize(1cm);

pair bhat=(2, 3);
real rot=30;
real a=1, b=.5, t=1.5;

path sq=(t,0)--(0,t)--(-t,0)--(0,-t)--(cycle);
filldraw(sq, green+opacity(.5), green*.5);

real sc = findr(ell(bhat, a, b, rot), sq, start=3);
real[] scs = {0.25sc, 0.5sc, 0.75sc, sc};

path[] ells;
for (int i: sequence(scs.length)) {
	ells[i] = ellipse(bhat, a, b, rot, scs[i]);
	draw(ells[i], red);
}

pair[] isect = intersectionpoints(sq, ells[3]);

dot("$\hat\beta_{Lasso}$", isect[0], green*.5, align=W);
dot("$\hat\beta_{OLS}$", bhat, red, align=W);

xaxis("$\beta_1$", xmin=-1.1*t, arrow=Arrow(TeXHead));
yaxis("$\beta_2$", ymin=-1.1*t, arrow=Arrow(TeXHead));

shipout("lasso");

currentpicture=pRidge;

unitsize(1cm);

pair bhat=(2, 3);
real rot=30;
real a=1, b=.5, t=1.5;

path sq=Circle((0, 0), t);
filldraw(sq, green+opacity(.5), green*.5);

real sc = findr(ell(bhat, a, b, rot), sq, start=3);
real[] scs = {0.25sc, 0.5sc, 0.75sc, sc};

path[] ells;
for (int i: sequence(scs.length)) {
	ells[i] = ellipse(bhat, a, b, rot, scs[i]);
	draw(ells[i], red);
}

pair[] isect = intersectionpoints(sq, ells[3]);

dot("$\hat\beta_{Ridge}$", isect[0], green*.5, align=W);
dot("$\hat\beta_{OLS}$", bhat, red, align=W);

xaxis("$\beta_1$", xmin=-1.1*t, arrow=Arrow(TeXHead));
yaxis("$\beta_2$", ymin=-1.1*t, arrow=Arrow(TeXHead));

shipout("ridge");

currentpicture=pLasso;
