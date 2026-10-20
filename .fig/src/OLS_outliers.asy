texpreamble("\usepackage{fontspec}\usepackage{unicode-math}\setmainfont{Noto Sans}\setmathfont{Noto Sans Math}");

settings.tex = "xelatex";
settings.outformat = "pdf";
settings.prc = false;
settings.render = 0;

unitsize(1cm);

import graph;
import stats;

int n = 25;
real a0 = 2, a1=.2;

pair[] gendata(int n) {
    pair[] res;
    for (int i=0; i < n; ++i) {
        res[i] = (i/2, a0 + a1 * i/2 + Gaussrand()*.5);
    }
    return res;
}

real[] getX(pair[] data) {
    int L = data.length;
    real[] x;
    for (int i; i < L; ++i) x[i] = data[i].x;
    return x;
}

real[] getY(pair[] data) {
    int L = data.length;
    real[] y;
    for (int i; i < L; ++i) y[i] = data[i].y;
    return y;
}

pair[] data = gendata(n);

for (int i=0; i < n; ++i) {
    filldraw(Circle(data[i], .1), black, black);
}

linefit reg0 = leastsquares(getX(data), getY(data));
draw((0, reg0.fit(0))--(24/2, reg0.fit(24/2)), black+1bp);

data[n] = (20/2, data[20].y + 5);
filldraw(Circle(data[n], .1), red, red);

linefit reg1 = leastsquares(getX(data), getY(data));
draw((0, reg1.fit(0))--(24/2, reg1.fit(24/2)), red+1bp);

xaxis("$x$", xmin=0, arrow=Arrow(TeXHead));
yaxis("$y$", ymin=0, arrow=Arrow(TeXHead));
