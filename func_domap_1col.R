


## This function aims to plot multiple maps in one column 


func_domap_1col <- function(n_plots, h) {  #@ the total number of plots
  
  fname <- paste0("./figures/do_maps_", which_dat, '_', which_do,  "_1col.png"); fname
  png(filename = fname, pointsize = 12, width = 7/3, height = h/3*n_plots, units="in", res = 300)
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = n_plots+1, ncol = 1, heights = c(0.01, rep(0.99/n_plots, n_plots))
  )))
  
  grid.text(label = paste0(which_dat, 'R'), vp = viewport(layout.pos.row = 1), gp = gpar(fontsize = 10, fontface = "bold"))
  
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
