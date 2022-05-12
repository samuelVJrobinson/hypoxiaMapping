

# library(remotes)
# install_github("r-tmap/tmaptools")
# install_github("r-tmap/tmap", force = T)
library(tmap)
library(raster)

library(tmaptools)
# palette_explorer()

getwd()

# f   <- "./figures/na_pixel_percent/ava_pixel_percent_IC_2014_chlor_a_90days.tif"
# f   <- "./figures/na_pixel_percent/ava_pixel_percent_IC_2014_sst_90days.tif"

pat <- "_80days.tif"

f_ls <- list.files("./figures/na_pixel_percent", pattern = paste0(".*", pat, "$"), full.names = T)
f_ls


for (f in f_ls) {
  
  
  ras <- raster(x = f)
  
  band <- basename(f) %>% 
    gsub('ava_pixel_pct_Y14_', '', .) %>% 
    gsub(pat, '', .)
  
  
  map <- 
    tm_shape(ras) + 
    tm_raster(style = "pretty", n = 5, ## pretty, equal, quantile, jenks, cont, cat, fixed, sd, kmeans, hclust, bclust, and fisher
              title = paste0("Pixel availability (%)"),
              # palette = rev(colorRampPalette( c("darkolivegreen4","yellow", "brown"))(12)),
              palette = "-viridis",
              legend.reverse = T, 
              legend.hist = TRUE)+
    tm_style("col_blind") +
    tm_layout(frame = FALSE,
              main.title = band, main.title.size = 1, main.title.fontface = "bold",
              main.title.position = c(0, 0.5),
              outer.margins = 0, inner.margins = 0,
              bg.color = NA) +
    tm_legend(outside = TRUE, 
              # legend.position = c("left","bottom"),
              # legend.position = "bottom",
              legend.format = list(fun = function(x) formatC(x, digits = 0, format = "f")),
              # legend.format = list(fun = function(x) scales::percent(x, accuracy = .01)),
              legend.outside.position = "bottom", 
              legend.stack = 'horizontal', 
              legend.hist.height = .7,
              legend.hist.width = 1)
  map
  ## save an image ("plot" mode)
  fname <- gsub('tif', 'png', basename(f)) %>% paste0(dirname(f), '/map_', .); fname
  tmap_save(tm = map, filename = fname, width = 6, height = 4, units = 'in', dpi = 200)
  
  
}


