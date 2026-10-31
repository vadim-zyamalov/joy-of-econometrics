library(Cairo)
library(ggplot2)
library(ggarrow)
library(ggbrace)

setwd("d:\\YandexDisk\\05. Преподавание\\Эконометрика (МГУ)\\data\\")

## -----------------------------------------------------------------------------
plot.theme <- theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = "none",
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank(),
  plot.margin = margin(
    t = 1,
    b = 1,
    l = 1,
    r = 1
  )
)

arrows <- function(x0,
                   x1,
                   y0,
                   y1,
                   ...) {
  .l <- max(length(x0), length(x1), length(y0), length(y1))
  .x0 <- rep_len(x0, .l)
  .x1 <- rep_len(x1, .l)
  .y0 <- rep_len(y0, .l)
  .y1 <- rep_len(y1, .l)
  res <- list()
  for (i in 1:(.l)) {
    res[[i]] <- lapply(i, function(z)
      geom_arrow(
        aes(x = c(.x0[z], .x1[z]), y = c(.y0[z], .y1[z])),
        ...
      ))
  }
  print(res)
  res
}

points <- function(x, y, labels) {
  .l <- max(length(x), length(y), length(labels))
  .x <- rep_len(x, .l)
  .y <- rep_len(y, .l)
  .labels <- rep_len(labels, .l)
  
  list(
    geom_point(
      aes(x = .x, y = .y),
      size = 4,
      shape = 21,
      fill = "black"
    ),
    annotate(
      "text",
      x = .x,
      y = .y + 0.2,
      label = .labels,
      vjust = "inward",
      angle = 45,
    )
  )
}

brace <- function(x, y, label, distance = 0) {
  list(
    stat_brace(
      aes(x = x, y = y),
      width = .2,
      size = .2,
      rotate = 180,
      distance = distance
    ),
    annotate(
      "text",
      x = mean(x),
      y = y - distance - .3,
      label = label,
      vjust = "inward",
      hjust = "center"
    )
  )
}

## ---- Fixed ------------------------------------------------------------------

p <- ggplot() +
  geom_arrow(aes(x = 1:35, y = 0), length = 5) +
  brace(c(1, 12), 0, "Базовый период") +
  brace(c(13, 25), 0, "Прогнозный период") +
  points(x = c(1, 12, 13, 25), y = 0, labels = c("1", "T", "T+1", "T+h")) +
  plot.theme

ggsave(
  "F.fixed.pdf",
  plot = p,
  width = 17.5,
  height = 2,
  units = "cm",
  scale = 1,
  device = cairo_pdf
)

## ---- Рекурсивный ------------------------------------------------------------

p <- ggplot() +
  geom_arrow(aes(x = 1:35, y = 1), length = 5) +
  geom_arrow(aes(x = 1:35, y = 0), length = 5) +
  arrows(x0 = c(1, 12, 13, 25),
         x1 = c(1, 13, 14, 26),
         y0 = 1,
         y1 = 0,
         linetype = "dashed",
         linewidth = .5,
         resect = 1.5,
         colour = "red") +
  brace(c(1, 12), 1, "Базовый период") +
  brace(c(13, 25), 1, "Прогнозный период") +
  brace(c(1, 13), 0, "Базовый период") +
  brace(c(14, 26), 0, "Прогнозный период") +
  points(c(1, 12, 13, 25), 1, c("1", "T", "T+1", "T+h")) +
  points(c(1, 13, 14, 26), 0, c("1", "T+1", "T+2", "T+h+1")) +
  plot.theme

ggsave(
  "F.recursive.pdf",
  plot = p,
  width = 17.5,
  height = 5,
  units = "cm",
  scale = 1,
  device = cairo_pdf
)

## ---- Скользящий -------------------------------------------------------------

p <- ggplot() +
  geom_arrow(aes(x = 1:35, y = 1), length = 5) +
  geom_arrow(aes(x = 1:35, y = 0), length = 5) +
  arrows(x0 = c(1, 12, 13, 25),
         x1 = c(2, 13, 14, 26),
         y0 = 1,
         y1 = 0,
         linetype = "dashed",
         linewidth = .5,
         resect = 1.5,
         colour = "red") +
  brace(c(1, 12), 1, "Базовый период") +
  brace(c(13, 25), 1, "Прогнозный период") +
  brace(c(2, 13), 0, "Базовый период") +
  brace(c(14, 26), 0, "Прогнозный период") +
  points(c(1, 12, 13, 25), 1, c("1", "T", "T+1", "T+h")) +
  points(c(2, 13, 14, 26), 0, c("2", "T+1", "T+2", "T+h+1")) +
  plot.theme

ggsave(
  "F.moving.pdf",
  plot = p,
  width = 17.5,
  height = 5,
  units = "cm",
  scale = 1,
  device = cairo_pdf
)

## ---- Наукастинг -------------------------------------------------------------

p <- ggplot() +
  geom_arrow(aes(x = 1:35, y = 0), length = 5) +
  annotate(
    "text",
    x = c(12, 25),
    y = 0.5,
    label = c("Момент\nпрогнозирования", "Текущий\nмомент"),
    vjust = "inward",
    hjust = "center"
  ) +
  brace(c(1, 12), 0, "Базовый период") +
  brace(c(13, 25), 0, "Наукастинг") +
  brace(c(26, 34), 0, "Прогноз") +
  points(c(1, 12, 13, 25, 26), 0, c("1", "T", "T+1", "T+H", "T+H+1")) +
  plot.theme

ggsave(
  "F.nowcasting.pdf",
  plot = p,
  width = 17.5,
  height = 3,
  units = "cm",
  scale = 1.5,
  device = cairo_pdf
)

## ---- Псевдовневыборочный-----------------------------------------------------

p <- ggplot() +
  geom_arrow(aes(x = 1:35, y = 0), length = 5) +
  annotate(
    "text",
    x = c(12, 25),
    y = 0.5,
    label = c("Момент\nпрогнозирования", "Текущий\nмомент"),
    vjust = "inward",
    hjust = "center"
  ) +
  brace(c(1, 12), 0, "Базовый период") +
  brace(c(13, 20), 0, "Прогнозный период") +
  brace(c(1, 25), -.4, "Известные величины") +
  points(c(1, 12, 13, 20, 25), 0, c("1", "T", "T+1", "T+h", "T+H")) +
  plot.theme

ggsave(
  "F.pseudo.pdf",
  plot = p,
  width = 17.5,
  height = 4,
  units = "cm",
  scale = 1.5,
  device = cairo_pdf
)
