#*****************************************************************************#
#------------------------------Ploting Function-------------------------------#

theme_CP <- function() {
  theme_bw()+
    theme(
      text = element_text(size = 7, lineheight = unit(0.5, "lines")), # lineheight is adjusting the space between lines
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      axis.line = element_line(linewidth = 0.05),
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
      panel.border = element_blank(), 
      panel.background = element_blank(),
      panel.grid = element_line(linewidth = 0.1), 
      strip.background = element_blank(),
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0, "cm"),
      legend.key.width = unit(0.2, "cm"),  # Make legend key width span the whole plot
      legend.key.height = unit(0.2, "cm"),  # Adjust legend key height
  )
}

# Chrysa_colors <- c("purple3",
#                    "aquamarine4",
#                    "royalblue2",
#                    "yellow3",
#                    "orange",
#                      "violetred4",
#                      "palegreen3",
#                      "lightskyblue1",
#                      "plum4")
# 
# Chrysa_colorsc <- c("mediumpurple1",
#                     "purple1",
#                     "mediumpurple4",
#                     "purple4")


# Chrysa_morecolors <- c("lavenderblush4",
#                    "violetred4",
#                    "lightgoldenrod",
#                    "aquamarine4",
#                    "royalblue2",
#                    "plum4",
#                    "yellow3",
#                    "palegreen3",
#                    "darkorchid4", 
#                    "lightskyblue2", 
#                    "coral")