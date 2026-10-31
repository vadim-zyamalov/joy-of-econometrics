library(ggplot2)
library(ggbrace)
library(latex2exp)

s <- seq(-10, 10, by = .1)

p <- ggplot(mapping = aes(x = s)) +
  ylab("G") +
  xlab(TeX("$s_t$")) +
  geom_line(aes(y = plogis(s, -3, 2)), colour = "red", linewidth = 1) +
  annotate(
    "text",
    x = -5.5,
    y = 0.375,
    label = TeX("$G^2,\\ c = -3$"),
  ) +
  geom_line(aes(y = plogis(s, +3, 2)), colour = "blue", linewidth = 1) +
  annotate(
    "text",
    x = 5.5,
    y = 0.625,
    label = TeX("$G^3,\\ c = 3$")
  ) +
  geom_line(aes(y = 1), colour = "black", linewidth = 1) +
  annotate(
    "text",
    x = -7.5,
    y = 0.95,
    label = TeX("$G^1,\\ c \\rightarrow -\\infty$")
  ) +
  geom_line(aes(y = 0), colour = "black", linewidth = 1) +
  annotate(
    "text",
    x = +7.5,
    y = 0.05,
    label = TeX("$G^4,\\ c \\rightarrow \\infty$")
  )

## Вертикальная линия
p <- p +
  geom_vline(aes(xintercept = -2.5), linetype = "dashed") +
  annotate("text", y = 0, x = -2.5, label = "-2.5", hjust = 1.1, vjust = 1.1)

## Скобки
p <- p +
  stat_brace(aes(x = -2.5, y = c(plogis(-2.5, -3, 2), 1)),
    width = .75,
    linewidth = .5,
    rotate = 270
  ) +
  annotate(
    "text",
    y = (1 + plogis(-2.5, -3, 2)) / 2,
    x = -4,
    label = TeX("$G^1 - G^2$")
  ) +
  stat_brace(aes(x = -2.5, y = c(plogis(-2.5, -3, 2), plogis(-2.5, 3, 2))),
    width = .75,
    linewidth = .5,
    rotate = 90
  ) +
  annotate(
    "text",
    y = (plogis(-2.5, -3, 2) + plogis(-2.5, 3, 2)) / 2,
    x = -1,
    label = TeX("$G^2 - G^3$")
  ) +
  stat_brace(aes(x = -2.5, y = c(0, plogis(-2.5, 3, 2))),
    width = .25,
    linewidth = .5,
    rotate = 90
  ) +
  annotate(
    "text",
    y = (0 + plogis(-2.5, 3, 2)) / 2,
    x = -1.5,
    label = TeX("$G^3 - G^4$")
  )

p <- p +
  geom_point(
    aes(x = -2.5, y = c(0, 1, plogis(-2.5, -3, 2), plogis(-2.5, 3, 2))),
    size = 2,
    shape = 21,
    fill = "white"
  )

print(p)

ggsave(
  "STAR.pdf",
  plot = p,
  device = "pdf",
  width = 15,
  height = 6,
  units = "cm",
  scale = 1.75
)

ggsave(
  "STAR.png",
  plot = p,
  device = "png",
  width = 15,
  height = 6,
  units = "cm",
  scale = 1.75
)
