texpreamble("\usepackage{fontspec}\usepackage{unicode-math}\setmainfont{Noto Sans}\setmathfont{Noto Sans Math}");

settings.tex = "xelatex";
settings.outformat = "pdf";
settings.prc = false;
settings.render = 0;

import graph3;
import three;
import myfunc;

currentprojection = orthographic((-0.5, 1, 0.5), up=Z);

// Orths
triple vec1 = unit((1, 1, 0)),
       vecX = unit((-1, 3, 0)),
       vecZ = (0, 0, 4);

triple nv = cross(vec1, vecX);

real   b0 = 1.5, b1 = 2.5;
triple Yr = (0.5, 3.5, 3);
//triple Yr = b0 * vec1 + b1 * vecX + vecZ;

// Projection on a plane
triple Yh = planeproject(nv) * Yr;

// Projection Yh on 1
real   coord1 = perp2orth(Yh, vec1);
triple proj1  = coord1 * vec1;

// Projection Yh on X
real   coordX = perp2orth(Yh, vecX);
triple projX  = coordX * vecX;

// LS coords
triple p1 = lintersection2(vecX, Yh, vec1, O);
real   b1 = length(p1);
triple pX = lintersection2(vec1, Yh, vecX, O);
real   bX = length(pX);
write(p1);
write(pX);

picture pic0, pic1;
picture[] pics  = {pic0, pic1};
string[]  names = {"OLS", "OLS_nocons"};

for (picture pic : pics) {
    currentpicture = pic;
    unitsize(2.5cm);

    // Pseudo-Axis
    draw(O--vec1, Arrow3(TeXHead2));
    label("$\mathbf{1}$", vec1, align=W);

    draw(O--vecX, Arrow3(TeXHead2));
    label("$\mathbf{x}$", vecX, align=SE);
}

currentpicture=pics[0];

// Projection on 1
draw(p1--Yh, dashed+blue);
label("$\beta_0\mathbf{1}$", p1, blue, align=NE);

// Projection on X
draw(pX--Yh, dashed+blue);
label("$\beta_1\mathbf{x}$", pX, blue, align=E);

// Fitted values
draw(O--Yh, red+linewidth(.5pt), Arrow3(TeXHead2));
label("$\hat{\mathbf{y}}$", (Yh+O)/2, red, align=SE);

// Residuals
draw(Yh--Yr, deepgreen+linewidth(.5pt), Arrow3(TeXHead2));
label("$\hat{\mathbf{u}}$", (Yh+Yr)/2, deepgreen, align=W);
perpmark3d(O, Yh, Yr, .1, deepgreen);

// Projection on 1 ONLY
label("$\bar{\mathbf{y}}=\bar{y}\mathbf{1}$", proj1, align=NW, black);

draw(Yr--proj1, dashed+black);
label("$\mathbf{y}-\bar{\mathbf{y}}$", (proj1+Yr)/2, align=NW, black);
perpmark3d(Yr, proj1, O, .1, black);

draw(Yh--proj1, dashed+red);
label("$\hat{\mathbf{y}}-\bar{\mathbf{y}}$", (proj1+Yh)/2, align=SW, red);
perpmark3d(Yh, proj1, O, .1, red);

incaxis(vec1, b1, coord1, 0.1);
incaxis(vecX, bX, 0.1);

currentpicture=pics[1];

// Fitted
draw(O--projX, red+linewidth(.5pt), Arrow3(TeXHead2));
label("$\hat{\mathbf{y}}$", (projX+O)/2, red, align=E);

// Residuals w/o constant
b1 = perp2orth((Yr-projX), vec1);
triple ort1P = vec1 * b1;
draw((Yr-projX)--ort1P, dashed+heavygreen);
label("$\langle\mathbf{1}, \tilde{\mathbf{u}}\rangle\neq0$", ort1P, align=NW, deepgreen);
perpmark3d((Yr-projX), ort1P, O, .1, deepgreen);

draw(projX--Yr, deepgreen, Arrow3(TeXHead2));
label("$\tilde{\mathbf{u}}$", (projX+Yr)/2, align=NE, deepgreen);
perpmark3d(Yr, projX, O, .1, deepgreen);

draw(O--(Yr-projX), deepgreen, Arrow3(TeXHead2));
label("$\tilde{\mathbf{u}}$", (O+Yr-projX)/2, align=NE, deepgreen);

draw(Yr--(Yr-projX), dotted+heavygreen);

incaxis(vec1, b1, 0.1);
incaxis(vecX, bX, coordX, 0.1);

for (int i : sequence(2)) {
    currentpicture = pics[i];

    // Real data
    draw(O--Yr, linewidth(1pt), Arrow3(TeXHead2));
    label("$\mathbf{y}$", (Yr+O)/2, align=NE);

    shipout(names[i]);
}

exit();
