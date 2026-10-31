import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math

plt.rcParams["figure.figsize"] = [6, 3]

x = np.linspace(-3, 3, 1000)
xc = np.linspace(0, 25, 1000)
y = stats.norm.pdf(x, 0, 1)
y2 = stats.norm.pdf(x, 0.5, 1)
yc = stats.chi2.pdf(xc, df=5)

fig = plt.figure()
ax = fig.add_subplot()
ax.axhline(0, color="black")
ax.plot(x, y)
ax.vlines(x=-1.96, ymin=0, ymax=np.max(y))
ax.text(
    -1.9,
    0,
    r"$t^{crit}_{\frac{\alpha}{2}}$",
    ha="left",
    va="bottom",
    size="xx-large",
)
ax.vlines(x=1.96, ymin=0, ymax=np.max(y))
ax.text(
    1.9,
    0,
    r"$t^{crit}_{1-\frac{\alpha}{2}}$",
    ha="right",
    va="bottom",
    size="xx-large",
)
ax.fill_between(x, y, where=list(np.abs(x) > 1.96), color="blue", alpha=0.3)
ax.annotate(
    r"$\frac{\alpha}{2}$",
    xy=(-2.25, 0.01),
    xytext=(-2.5, 0.15),
    arrowprops=dict(arrowstyle="->"),
    size="xx-large",
)
ax.annotate(
    r"$\frac{\alpha}{2}$",
    xy=(2.25, 0.01),
    xytext=(2.5, 0.15),
    arrowprops=dict(arrowstyle="->"),
    size="xx-large",
)
plt.savefig("fig_twoside.png")

fig = plt.figure()
ax = fig.add_subplot()
ax.axhline(0, color="black")
ax.plot(x, y)
ax.vlines(x=-1.64, ymin=0, ymax=np.max(y))
ax.text(
    -1.59, 0, r"$t^{crit}_{\alpha}$", ha="left", va="bottom", size="xx-large"
)
ax.fill_between(x, y, where=list(x < -1.64), color="blue", alpha=0.3)
ax.annotate(
    r"$\alpha$",
    xy=(-1.96, 0.01),
    xytext=(-1.96, 0.15),
    arrowprops=dict(arrowstyle="->"),
    ha="right",
    size="xx-large",
)
plt.savefig("fig_left.png")

fig = plt.figure()
ax = fig.add_subplot()
ax.axhline(0, color="black")
ax.plot(x, y)
ax.vlines(x=1.64, ymin=0, ymax=np.max(y))
ax.text(
    1.59, 0, r"$t^{crit}_{1-\alpha}$", ha="right", va="bottom", size="xx-large"
)
ax.fill_between(x, y, where=list(x > 1.64), color="blue", alpha=0.3)
ax.annotate(
    r"$\alpha$",
    xy=(1.96, 0.01),
    xytext=(1.96, 0.15),
    arrowprops=dict(arrowstyle="->"),
    ha="left",
    size="xx-large",
)
plt.savefig("fig_right.png")

fig = plt.figure()
ax = fig.add_subplot()
ax.axhline(0, color="black")
ax.plot(x, y)
ax.vlines(x=-1, ymin=0, ymax=np.max(y))
ax.text(-0.95, 0, r"$-|t|$", ha="left", va="bottom", size="xx-large")
ax.vlines(x=1, ymin=0, ymax=np.max(y))
ax.text(0.95, 0, r"$|t|$", ha="right", va="bottom", size="xx-large")
ax.fill_between(x, y, where=list(x < -1), color="blue", alpha=0.3)
ax.fill_between(x, y, where=list(x > 1), color="blue", alpha=0.3)
plt.savefig("fig_p_twoside.png")

fig = plt.figure()
ax = fig.add_subplot()
ax.axhline(0, color="black")
ax.plot(x, y)
ax.vlines(x=1, ymin=0, ymax=np.max(y))
ax.text(1.05, 0, r"$t$", ha="left", va="bottom", size="xx-large")
ax.fill_between(x, y, where=list(x < 1), color="blue", alpha=0.3)
plt.savefig("fig_p_left.png")

fig = plt.figure()
ax = fig.add_subplot()
ax.axhline(0, color="black")
ax.plot(x, y)
ax.vlines(x=1, ymin=0, ymax=np.max(y))
ax.text(0.95, 0, r"$t$", ha="right", va="bottom", size="xx-large")
ax.fill_between(x, y, where=list(x > 1), color="blue", alpha=0.3)
plt.savefig("fig_p_right.png")

fig = plt.figure()
ax = fig.add_subplot()
ax.axhline(0, color="black")
ax.plot(x, y)
ax.vlines(x=-1.5, ymin=0, ymax=np.max(y2))
ax.vlines(x=1.5, ymin=0, ymax=np.max(y2))
ax.plot(x, y2, color="red")
ax.fill_between(x, y2, where=list(np.abs(x) < -1.5), color="red", alpha=0.3)
ax.fill_between(x, y2, where=list(np.abs(x) > 1.5), color="red", alpha=0.3)
plt.savefig("fig_power.png")

fig = plt.figure()
ax = fig.add_subplot()
ax.axhline(0, color="black")
ax.plot(xc, yc)
plt.savefig("fig_chi.png")
