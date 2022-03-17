library(readr)
library(tidyverse)
library(cowplot)

getwd()



## temporal change of dz area in 2014
csv <- read_csv('./data/data_from_gee/area_dz_over_time.csv') %>%
  rename_at(vars(starts_with("system")), ~'Date') %>%
  dplyr::mutate(Date = as.Date(Date, "%B %d, %Y"))


ggplot(csv, aes(x = Date, y = Area)) +
  geom_point() +
  geom_smooth(method = 'loess', formula = 'y ~ x', color = 'darkorange') +
  ylab('Area (sq km)')+
  theme_bw()

fname <- paste0('./figures/', 'area_dz_over_time_2014.png'); fname
ggsave(filename = fname, plot = last_plot(), width = 6, height = 4, units = 'in', dpi = 300)



summary(csv$Area)


# July 26 - August 3
csv %>% 
  dplyr::filter(Date >= as.Date('2014-07-19'), 
                Date <= as.Date('2014-08-15')) %>%
  summary()
  