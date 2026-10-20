import matplotlib.pyplot as plt
import numpy as np


def param_plot(D="all"):
    """This function creates the graph on page 189 of
    Sargent Macroeconomic Theory, second edition, 1987.
    """
    if D not in ("g", "e", "l", "all"):
        raise ValueError()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_aspect("equal")

    # Стили стрелок
    plot_arrowsl = {"arrowstyle": "->", "connectionstyle": "arc3, rad=-0.2"}
    plot_arrowsr = {"arrowstyle": "->", "connectionstyle": "arc3, rad=0.2"}
    plot_text = {
        "horizontalalignment": "center",
        "verticalalignment": "bottom",
        "transform": ax.transData,
    }

    # ~~~~~~~~~~~~~ ОСИ ~~~~~~~~~~~~~
    # Диапазоны по осям
    xmin, ymin = -2.5, -1.5
    xmax, ymax = -xmin, -ymin
    plt.axis([xmin, xmax, ymin, ymax])

    # Названия осей
    ax.set(xticks=[], yticks=[])
    ax.set_xlabel(r"$a_2$", fontsize=16)
    ax.xaxis.set_label_position("top")
    ax.set_ylabel(r"$a_1$", rotation=0, fontsize=16)
    ax.yaxis.set_label_position("right")

    # Отключаем стандартные оси
    for spine in ["left", "bottom", "top", "right"]:
        ax.spines[spine].set_visible(False)

    # Делаем собственные оси
    axes_arrows = {"arrowstyle": "<-", "lw": 1.3}
    ax.annotate("", xy=(xmin, 0), xytext=(xmax, 0), arrowprops=axes_arrows)
    ax.annotate("", xy=(0, ymin), xytext=(0, ymax), arrowprops=axes_arrows)

    # ~~~~~~~~~~~~~ ЛИНИИ ~~~~~~~~~~~~~
    a1 = np.linspace(-2, 2, 100)
    zz = -abs(a1) + 1
    hl = np.full_like(a1, -1)
    arc = -(a1**2 / 4)

    match D:
        case "g":
            ax.plot(a1, zz, "k:")
            ax.plot(a1, arc, "k:")
            ax.fill_between(a1, zz, arc, color="red", alpha=0.4)

        case "l":
            ax.plot(a1, hl, "k:")
            ax.plot(a1, arc, "k:")
            ax.fill_between(a1, arc, hl, color="blue", alpha=0.4)

        case "e":
            ax.plot(a1, arc, c="orange", linewidth=6, alpha=0.4)
            ax.plot(a1, arc, c="black")

        case _:
            ax.plot(a1, zz, "k:")
            ax.plot(a1, hl, "k:")
            ax.plot(a1, arc, c="black")
            ax.fill_between(a1, zz, arc, color="red", alpha=0.4)
            ax.fill_between(a1, arc, hl, color="blue", alpha=0.4)

    # ~~~~~~~~~~~~~ ПОДПИСИ ~~~~~~~~~~~~~
    match D:
        case "g":
            ax.annotate(
                r"$a_1 + a_2 < 1$",
                xy=(0.5, 0.3),
                xytext=(0.8, 0.6),
                arrowprops=plot_arrowsr,
                fontsize="12",
            )
            ax.text(
                1.5,
                -0.5,
                r"$a_1 + a_2 = 1$",
                rotation=-45,
                fontsize="12",
                **plot_text,
            )
            ax.annotate(
                r"$a_2 < 1 + a_1$",
                xy=(-0.5, 0.3),
                xytext=(-1.3, 0.6),
                arrowprops=plot_arrowsl,
                fontsize="12",
            )
            ax.text(
                -1.5,
                -0.5,
                r"$a_2 - a_1 = 1$",
                rotation=45,
                fontsize="12",
                **plot_text,
            )
            ax.annotate(
                r"${a_1}^2 + 4a_2 > 0$",
                xy=(-0.25, 0.25),
                xytext=(-1, -0.6),
                arrowprops=plot_arrowsl,
                fontsize="12",
            )
        case "l":
            ax.text(0.75, -1.15, r"$a_2 = -1$", fontsize="12", **plot_text)
            ax.annotate(
                r"${a_1}^2 + 4a_2 < 0$",
                xy=(-0.5, -0.5),
                xytext=(-1, 0.25),
                arrowprops=plot_arrowsl,
                fontsize="12",
            )
    if D in ("g", "e"):
        ax.plot([-2, -2], [0, -1], "k:")
        ax.plot([2, 2], [0, -1], "k:")
    if D in ("g", "e", "l"):
        ax.text(
            0.75,
            -0.45,
            r"$a_1^2 + 4 a_2 = 0$",
            rotation=-20,
            fontsize="12",
            **plot_text,
        )

    if D == "all":
        # Label categories of solutions
        ax.text(1.5, 1, "Взрывной\nрост", ha="center", fontsize=16)
        ax.text(-1.5, 1, "Взрывные\nосцилляции", ha="center", fontsize=16)
        ax.text(0.05, -1.5, "Взрывные\nосцилляции", ha="center", fontsize=16)
        ax.text(0.09, -0.75, "Ограниченные\nосцилляции", ha="center", fontsize=16)

    # Add small marker to y-axis
    ax.axhline(y=1, xmin=0.495, xmax=0.505, c="black")
    ax.text(-0.12, 0.98, "1", fontsize=10)

    ax.axhline(y=-1, xmin=0.495, xmax=0.505, c="black")
    ax.text(-0.15, -1.1, "-1", fontsize=10)

    ax.axvline(x=1, ymin=0.49, ymax=0.51, c="black")
    ax.text(1.05, 0.05, "1", fontsize=10)

    ax.axvline(x=2, ymin=0.49, ymax=0.51, c="black")
    ax.text(2.05, 0.05, "2", fontsize=10)

    ax.axvline(x=-1, ymin=0.49, ymax=0.51, c="black")
    ax.text(-1.15, 0.05, "-1", fontsize=10)

    ax.axvline(x=-2, ymin=0.49, ymax=0.51, c="black")
    ax.text(-2.15, 0.05, "-2", fontsize=10)

    fig.tight_layout()

    return fig


for m in ("g", "e", "l", "all"):
    param_plot(m).savefig(f"../AR2.sols.{m}.pdf", format="pdf")
