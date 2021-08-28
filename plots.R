


getwd()

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)


source('helperFunctions.R') ## mae() function



### to combine sub-plots for Figure 2

## 1. load llModErr and fdaErr --------------------------------------------------------------------
load(file = './data/errDat.Rdata')





## 2. load RF model error data --------------------------------------------------------------------
root <- dirname(getwd())

fname <- paste0(root, '/Dead_Zone_telecoupling/data/results_RF/err_rf_Day_DObottom_2014_2014.Rdata'); fname
load(fname)

acc_set <- 'testing'

### input data for plotting 
acc0 <- err_rf %>%
  gather(key = 'errType', value = 'value', RMSE:R2) %>%
  dplyr::mutate(errType = factor(errType, levels = c('R2', 'RMSE', 'MAE')))

### option 1. use mean, se OR sd,for plotting ----
###   `summarySE` provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
library(Rmisc)
rfErr <- summarySE(data = acc0, measurevar="value", groupvars=c("year", "nw", "errType")) %>%
  dplyr::select(-year, -N) %>%
  dplyr::rename(lag = nw)



### option 2. use median, max, min for plotting -----
rfErr <- acc0 %>%
  dplyr::select(-year) %>%
  dplyr::rename(lag = nw) %>%
  group_by(lag, errType) %>%
  dplyr::summarise(med = median(value, na.rm = T),
                   max = max(value, na.rm = T), 
                   min = min(value, na.rm = T), 
                   iqr=IQR(value, na.rm = T)) %>%
  dplyr::mutate(withinSampErr = NA)
  






## 3. bind three data sets -------------------------------------------------------------------------

names(fdaErr)
names(llModErr)
names(rfErr)

err_2models <- rbind(llModErr %>% dplyr::mutate(mod = 'llr'),
                     rfErr    %>% dplyr::mutate(mod = 'rf'))


(p1 <- ggplot(err_2models, aes(x=lag,y=med))+
    geom_ribbon(aes(ymax=max,ymin=min),alpha=0.3)+
    geom_line()+
    geom_line(aes(y=withinSampErr),col='red')+
    facet_grid(errType~mod, scales = 'free_y')+
    labs(x='Time lag',y='Model Accuracy') +
    theme_bw()
)




(p2 <- fdaErr %>% 
  dplyr::mutate(mod = 'fda') %>% 
  ggplot(aes(x=value))+#geom_freqpoly()+
  geom_histogram(fill='black', alpha=0.3, binwidth = 0.02)+
  geom_vline(aes(xintercept=med),col='black')+
  geom_vline(data=fdaWithin,aes(xintercept=value),col='red')+
    facet_grid(name~mod,scales='free_y')+
  
  coord_flip() +
    labs(x='',y='Count') +
    theme_bw()
)

(p <- ggarrange(p1,p2,nrow=1, widths = c(2,1)))

ggsave('./figures/outOfSamp_comparison_3models.png', p, width=8, height=4)
