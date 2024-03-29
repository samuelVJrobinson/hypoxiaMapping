


## This function aims to plot multiple maps in one column 


func_domap_1col <- function(n_plots, h, postfix = "") {  #@ the total number of plots
  
  fname <- paste0("./figures/do_maps_", which_dat, '_', which_do,  "_1col", postfix, ".png"); fname
  png(filename = fname, pointsize = 12, width = 7/3, height = h/3*n_plots, units="in", res = 600)
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = n_plots+1, ncol = 1, heights = c(0.03, rep(0.98/n_plots, n_plots))
  )))
  
  
  if (which_dat == 'FR') {
    text_header <- 'Functional Data Analysis'
  } else if (which_dat == 'LL') {
    # text_header <- paste0(which_dat, 'R')
    text_header <- 'Lagged Linear Regression'
  } else {
    text_header <- 'Random Forest Regression'
  }
  
  grid.text(label = text_header, vp = viewport(layout.pos.row = 1), gp = gpar(fontsize = 9, fontface = "bold"))
  
  z = 0
  for (i in 1:n_plots) {
    # print(i)
    # for (j in 1:3) {
    z = z + length(i)
    print(z)
    
    if (z>n_plots) {
      break
    }
    print(domap_ls[z], vp = viewport(layout.pos.row = z+1))
    # }
  }
  
  dev.off()
}
