## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, echo=FALSE-----------------------------------------
library(ggplot2)
library(magrittr)
library(CurveAware)

## ---- results='asis', echo=FALSE, fig.align='center', fig.width=6, fig.height=3----
set.seed(48)
infxn <- data.frame(time = 1:100, It = sapply(1:100, function(x, I0, r){I0*exp(r*x)}, I0 = 50, r = 0.10))
tod <- data.frame(time = 1:15, todt = sapply(1:15, function(x, m_od, s_od){pgamma(x, shape = 1/s_od^2, scale = m_od*s_od^2)}, m_od = 10, s_od = 0.5))


#..................
# plot
#..................

infxnplot <- infxn %>% 
  ggplot() + 
  geom_line(aes(x = time, y = It), color = "#3288bd", size = 3) +
  geom_area(aes(x = ifelse(time >= 0 & time <= 100, time, 0), y = It), fill = "#3288bd", alpha = 0.5) +
  theme_minimal() +
  xlab("Time") + ylab("Num. Infected") +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "#000000", size = 2)
  )
  
Xplot <- ggplot() +
  geom_text(aes(x = 1, y = 1, label = paste('bolditalic("X")')), parse = T, size = 15) + 
  xlim(0.9, 1.1) + ylim(0.9, 1.1) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_blank(), 
    strip.background = element_blank(),
  )

todplot <- tod %>% 
  ggplot() + 
  geom_line(aes(x = time, y = todt), color = "#d53e4f", size = 3) +
  geom_area(aes(x = ifelse(time >= 0 & time <= 100, time, 0), y = todt), fill = "#d53e4f", alpha = 0.5) +
  theme_minimal() +
  xlab("Time") + ylab("Prob. Onset to Death") +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "#000000", size = 2)
  )
  

cowplot::plot_grid(infxnplot, Xplot, todplot, 
                   rel_heights = c(2, 0.75, 2),
                   rel_widths = c(2, 0.75, 2), 
                   nrow = 1)



