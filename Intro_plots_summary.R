library(ggpubr)
library(gridExtra)
library(moments)
library(MASS)
library(ggplot2)
library(tidyverse)

# Specify the assumptions
set.seed(1076354)
s1 <- 1
s2 <- 1
compeffect <- 0.5
treateffect <- 0.7
exp_effect <- (treateffect - compeffect) * 100 / compeffect
size <- 800
rho1 <- c(-0.9, -0.6, -0.2, 0.3, 0.7, 0.9)
bvn <- NULL

# Simulate dataset for different correlations
for (rho in rho1) {
  bvn1 <-
    data.frame(cbind(mvrnorm(
      ss,
      mu = c(compeffect, treateffect),
      Sigma = matrix(c(s1, s1 * s2 * rho, s1 *
        s2 * rho, s2), 2)
    ), rho))
  colnames(bvn1) <- c("X1", "X2", "rho")
  bvn1$pch <- ((bvn1$X2 / bvn1$X1) - 1) * 100
  bvn <- rbind(bvn, bvn1)
}

# Add delta estimates to the data
delta0 <-
  ((mean((
    bvn %>% filter(rho == -0.9)
  )$X2) / mean((
    bvn %>% filter(rho == -0.9)
  )$X1)) - 1) * 100
delta1 <-
  ((mean((
    bvn %>% filter(rho == -0.6)
  )$X2) / mean((
    bvn %>% filter(rho == -0.6)
  )$X1)) - 1) * 100
delta2 <-
  ((mean((
    bvn %>% filter(rho == -0.2)
  )$X2) / mean((
    bvn %>% filter(rho == -0.2)
  )$X1)) - 1) * 100
delta3 <-
  ((mean((
    bvn %>% filter(rho == 0.3)
  )$X2) / mean((
    bvn %>% filter(rho == 0.3)
  )$X1)) - 1) * 100
delta4 <-
  ((mean((
    bvn %>% filter(rho == 0.7)
  )$X2) / mean((
    bvn %>% filter(rho == 0.7)
  )$X1)) - 1) * 100
delta5 <-
  ((mean((
    bvn %>% filter(rho == 0.9)
  )$X2) / mean((
    bvn %>% filter(rho == 0.9)
  )$X1)) - 1) * 100

bvn <- bvn %>% mutate(
  delta = case_when(
    rho == -0.9 ~ delta0,
    rho == -0.6 ~ delta1,
    rho == -0.2 ~ delta2,
    rho == 0.3 ~ delta3,
    rho == 0.7 ~ delta4,
    rho == 0.9 ~ delta5
  )
)

# plot for percent change from baseline with varying correlations
pch_plot <- ggplot(bvn, aes(x = pch, colour = factor(rho))) +
  geom_density() +
  geom_vline(
    xintercept = exp_effect,
    linetype = "solid",
    color = "black",
    size = 0.4
  ) +
  annotate(
    geom = "text",
    x = c(240),
    y = c(0.0075),
    label = c("Expected Effect"),
    col = c("black")
  ) +
  xlab("Percent change from baseline") +
  ylab("density") +
  scale_x_continuous(
    limits = c(-700, 700),
    breaks = c(-700, -400, 0, exp_effect, 400, 700)
  ) +
  theme(axis.ticks.length = unit(.25, "cm")) +
  scale_color_manual(
    name = "Correlation",
    labels = c("-0.9", "-0.6", "-0.2", "0.3", "0.7", "0.9"),
    values = c("red", "green", "blue", "black", "orange", "grey")
  )

tiff(
  "Figure_1.tiff",
  units = "in",
  width = 6.1,
  height = 4,
  res = 310,
  compression = "lzw"
)
ggarrange(pch_plot,
  ncol = 1,
  nrow = 2,
  heights = c(1, 0.2)
)
dev.off()

# QQ-plot
qqplot <- ggplot(bvn, aes(sample = pch, colour = factor(rho))) +
  geom_qq() +
  geom_qq_line() +
  theme(legend.position = "none") +
  ylim(-10000, 10000) + # ylim(-30,30) +
  facet_wrap(~rho,
    nrow = 2,
    ncol = 3,
    scales = "free"
  ) +
  ggtitle("QQ-plot for Normality check") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(
    name = "Correlation",
    labels = c("-0.9", "-0.6", "-0.2", "0.3", "0.7", "0.9"),
    values = c("red", "green", "blue", "black", "orange", "grey")
  )
tiff(
  "Figure_2.tiff",
  units = "in",
  width = 6.1,
  height = 4,
  res = 310,
  compression = "lzw"
)
ggarrange(qqplot,
  ncol = 1,
  nrow = 1,
  heights = c(2)
)
dev.off()

# summary from the simulated data with skewness and Kurtosis
stat <- bvn %>%
  group_by(rho) %>%
  summarize(
    Mean = round(mean(pch), 3),
    SD = round(sqrt(var(pch)), 3),
    Median = round(median(pch), 3),
    Min = round(min(pch), 3),
    Max = round(max(pch), 3),
    skewness = round(skewness(pch), 3),
    Kurtosis = round(kurtosis(pch), 3),
    # Shapiro-Wilk's normality test
    # "Test statistic" = round(shapiro.test(pch)$statistic, 5),
    # "KS" = ks.test(pch, "pnorm")$statistic,    #Kolmogorov-Smirnov statistic.
    "Shapiro–Wilk\n normality test,\n P-value" = case_when(
      shapiro.test(bvn$pch)$p.value < 0.0001 ~ "<0.0001",
      shapiro.test(bvn$pch)$p.value >= 0.0001
      ~ paste0(round(shapiro.test(bvn$pch)$p.value, 4))
    )
  )
table1 <- ggtexttable(stat,
  rows = NULL,
  theme = ttheme("mOrange")
)
tiff(
  "/cloud/project/Percent-change-from-baseline/Table_1.tiff",
  units = "in",
  width = 9,
  height = 4.5,
  res = 310,
  compression = "lzw"
)
ggarrange(table1,
  ncol = 1,
  nrow = 2,
  heights = c(0.5)
)
dev.off()

# Shapiro-Wilk's normality test
shapiro <- bvn %>%
  group_by(rho) %>%
  summarize(
    "Test statistic" = round(shapiro.test(pch)$statistic, 5),
    # "KS" = ks.test(pch, "pnorm")$statistic,    #Kolmogorov-Smirnov statistic.
    "P-value" = case_when(
      shapiro.test(bvn$pch)$p.value < 0.0001 ~ "<0.0001",
      shapiro.test(bvn$pch)$p.value >= 0.0001
      ~ paste0(round(shapiro.test(bvn$pch)$p.value, 4))
    )
  )

table2 <- ggtexttable(shapiro,
  rows = NULL,
  theme = ttheme("mOrange")
)
tiff(
  "Table_2_old.tiff",
  units = "in",
  width = 6.1,
  height = 4,
  res = 310,
  compression = "lzw"
)
ggarrange(table2,
  ncol = 1,
  nrow = 2,
  heights = c(0.5)
)
dev.off()

# facet plots with Delta estimate for the density plots in Introduction
density <-
  bvn %>% ggplot(aes(x = pch, colour = factor(rho))) +
  geom_density() +
  geom_vline(
    xintercept = 40,
    linetype = "solid",
    color = "yellow",
    size = 0.4
  ) +
  geom_vline(
    data = filter(bvn, rho == -0.9),
    aes(xintercept = delta0),
    colour = "red",
    linetype = "solid",
    size = 0.4
  ) +
  geom_vline(
    data = filter(bvn, rho == -0.6),
    aes(xintercept = delta1),
    colour = "green",
    linetype = "solid",
    size = 0.4
  ) +
  geom_vline(
    data = filter(bvn, rho == -0.2),
    aes(xintercept = delta2),
    colour = "blue",
    linetype = "solid",
    size = 0.4
  ) +
  geom_vline(
    data = filter(bvn, rho == 0.3),
    aes(xintercept = delta3),
    colour = "black",
    linetype = "solid",
    size = 0.4
  ) +
  geom_vline(
    data = filter(bvn, rho == 0.7),
    aes(xintercept = delta4),
    colour = "orange",
    linetype = "solid",
    size = 0.4
  ) +
  geom_vline(
    data = filter(bvn, rho == 0.9),
    aes(xintercept = delta5),
    colour = "grey",
    linetype = "solid",
    size = 0.4
  ) +
  facet_wrap(~rho,
    nrow = 2,
    ncol = 3,
    scales = "free"
  ) +
  xlab("Percent change from baseline") +
  ylab("density") +
  scale_x_continuous(
    limits = c(-400, 400),
    breaks = c(-400, 0, 400)
  ) +
  scale_color_manual(
    name = "Correlation",
    labels = c("-0.9", "-0.6", "-0.2", "0.3", "0.7", "0.9"),
    values = c("red", "green", "blue", "black", "orange", "grey")
  )

tiff(
  "Figure_3.tiff",
  units = "in",
  width = 6.1,
  height = 4,
  res = 310,
  compression = "lzw"
)
ggarrange(density,
  ncol = 1, nrow = 1
)
dev.off()
