library(ggplot2)
library(ggh4x)
library(colorBlindness)

# colors for different levels
siteOverallColors <- c("#990F0F", "#B22C2C", "#CC5151", "#E57E7E", "#FFB2B2")
contextColors <- c("#6B990F", "#85B22C", "#A3CC51", "#C3E57E", "#E5FFB2")
tfContextColors <- c("#260F99", "#422CB2", "#6551CC", "#8F7EE5", "#BFB2FF")
tfColors <- c("#0F6B99", "#2C85B2", "#51A3CC", "#7EC3E5", "#B2E5FF")

# colors for metrics
# AUCPR: viridis option B
# PPV 0.05: viridis option A
# PPV 0.1: viridis option D

# colors for methods other discrete groups
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
methodColors <- c("maxATAC"="#CC79A7",
                  "maxATAC.default"="#CC79A7",
                  "maxATAC.100"="#CC79A7",
                  "TFBlearner.full"="#56B4E9",
                  "TFBlearner.basic"="#0072B2",
                  "BMO"="#009E73",
                  "TOP"="#F0E442")



theme_sleek <- function() {
  theme_minimal() +
    theme(
      plot.tag = element_text(size = 16, face = "bold"),
      
      plot.title = element_text(size = 12, face = "bold", margin = margin(b = 8)),
      plot.subtitle = element_text(size = 10, color = "#666666", margin = margin(b = 6)),
      axis.title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 9, color = "#333333"),
      
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      legend.background = element_rect(fill = NA, color = NA),
      
      panel.grid.major = element_line(color = "#e0e0e0", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), #"#fafafa"
      
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 8, 8, 8),
      
      strip.text = element_text(size = 9, face = "bold"),
      strip.background = element_rect(fill = "#f0f0f0", color = NA)
    )
}
