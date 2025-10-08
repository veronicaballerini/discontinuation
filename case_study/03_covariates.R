################################################################################ 
# Evaluating causal effects on time-to-event outcomes in an RCT in Oncology 
# with treatment discontinuation
#
# Authors: V. Ballerini, B., Bornkamp, F. Mealli, C. Wang, Y. Zhang, A. Mattei
#
# Replication Package for results in Section 6 
# In particular: Figures 16, 17, 18.
#
# Code authors: Veronica Ballerini
# Last modified: October 6, 2025 
################################################################################ 

mycolors_x1 <- c("#A6CEE3", "#1F78B4", "#B2DF8A")
mycolors_x2 <- c("#FB9A99", "#E31A1C", "#FF7F00")
mycolors_x3 <- c("#CAB2D6", "#6A3D9A", "#B15928")

### Load results load('results.RData')

m1nd <- rep(NA, length(cov_distr))
for (i in 1:length(cov_distr)) {
    m1nd[i] <- cov_distr[[i]]$m_nd[1]
}
m1ed <- rep(NA, length(cov_distr))
for (i in 1:length(cov_distr)) {
    m1ed[i] <- cov_distr[[i]]$m_ed[1]
}
m1ld <- rep(NA, length(cov_distr))
for (i in 1:length(cov_distr)) {
    m1ld[i] <- cov_distr[[i]]$m_ld[1]
}

m2nd <- rep(NA, length(cov_distr))
for (i in 1:length(cov_distr)) {
    m2nd[i] <- cov_distr[[i]]$m_nd[2]
}
m2ed <- rep(NA, length(cov_distr))
for (i in 1:length(cov_distr)) {
    m2ed[i] <- cov_distr[[i]]$m_ed[2]
}
m2ld <- rep(NA, length(cov_distr))
for (i in 1:length(cov_distr)) {
    m2ld[i] <- cov_distr[[i]]$m_ld[2]
}

m3nd <- rep(NA, length(cov_distr))
for (i in 1:length(cov_distr)) {
    m3nd[i] <- cov_distr[[i]]$m_nd[3]
}
m3ed <- rep(NA, length(cov_distr))
for (i in 1:length(cov_distr)) {
    m3ed[i] <- cov_distr[[i]]$m_ed[3]
}
m3ld <- rep(NA, length(cov_distr))
for (i in 1:length(cov_distr)) {
    m3ld[i] <- cov_distr[[i]]$m_ld[3]
}

df1 <- data.frame(value = c(m1nd, m1ld, m1ed), group = factor(rep(c("ND",
    "Late D", "Early D"), times = c(length(m1nd), length(m1ld),
    length(m1ed)))))

df1$group <- factor(df1$group, levels = c("ND", "Late D",
    "Early D"))

jpeg("figures/Figure16.jpeg", width = 4500, height = 2864,
    units = "px", res = 500)
ggplot(df1, aes(x = value, color = group, fill = group)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
        alpha = 0.4, position = "identity") + scale_color_manual(values = mycolors_x1) +
    scale_fill_manual(values = mycolors_x1) + coord_cartesian(xlim = c(-0.4,
    0.7), ylim = c(0, 10)) + labs(x = expression(X[1]),
    y = "", color = NULL, fill = NULL) + theme_minimal(base_size = 16) +
    theme(legend.position = "top")
dev.off()

df2 <- data.frame(value = c(m2nd, m2ld, m2ed), group = factor(rep(c("ND",
    "Late D", "Early D"), times = c(length(m1nd), length(m1ld),
    length(m1ed)))))
df2$group <- factor(df2$group, levels = c("ND", "Late D",
    "Early D"))

jpeg("figures/Figure17.jpeg", width = 4500, height = 2864,
    units = "px", res = 500)
ggplot(df2, aes(x = value, color = group, fill = group)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
        alpha = 0.4, position = "identity") + scale_color_manual(values = mycolors_x2) +
    scale_fill_manual(values = mycolors_x2) + coord_cartesian(xlim = c(0.2,
    0.8), ylim = c(0, 20)) + labs(x = expression(X[2]),
    y = "", color = NULL, fill = NULL) + theme_minimal(base_size = 16) +
    theme(legend.position = "top")
dev.off()

df3 <- data.frame(value = c(m3nd, m3ld, m3ed), group = factor(rep(c("ND",
    "Late D", "Early D"), times = c(length(m1nd), length(m1ld),
    length(m1ed)))))
df3$group <- factor(df3$group, levels = c("ND", "Late D",
    "Early D"))

jpeg("figures/Figure18.jpeg", width = 4500, height = 2864,
    units = "px", res = 500)
ggplot(df3, aes(x = value, color = group, fill = group)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
        alpha = 0.4, position = "identity") + scale_color_manual(values = mycolors_x3) +
    scale_fill_manual(values = mycolors_x3) + coord_cartesian(xlim = c(0.05,
    0.6), ylim = c(0, 25)) + labs(x = expression(X[3]),
    y = "", color = NULL, fill = NULL) + theme_minimal(base_size = 16) +
    theme(legend.position = "top")
dev.off()
