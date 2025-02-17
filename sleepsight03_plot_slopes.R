
install.packages('RColorBrewer')
install.packages('beeswarm')
install.packages('effsize')

library('RColorBrewer')
library('beeswarm')
library('effsize')



# define paths
datdir = file.path(.Platform$file.sep, 'Users','petzka', 'Documents',
                   'projects', 'MPI', 'sleepsight_upload', 'eeg')

folder_dat = 'eeg4stats'


### define colour schemes
sleep_cols = c(brewer.pal(n = 9, name = "YlGn")[6],brewer.pal(n = 11, name = "Spectral")[9], brewer.pal(n = 11, name = "BrBG")[10])
insight_cols = c(brewer.pal(n = 9, name = "Blues")[4], brewer.pal(n = 9, name = "YlGnBu")[7])

########################
### plot figure 3A
########################

title_size = 2.5
label_size = 2

datOI = '1fslope_C4.csv'
dat = read.csv(file.path(datdir, folder_dat, datOI), 
               header = TRUE)
dat$sleep = as.factor(dat$sleep)
#dat$slope = dat$slope*-1

# barplot of 1/f slopes
mean0 = mean(dat$slope[dat$sleep == 0])
mean1 = mean(dat$slope[dat$sleep == 1])
mean2 = mean(dat$slope[dat$sleep == 2])

mean_ = c(mean0, mean1, mean2)

t.test(dat$slope[dat$sleep == 0], dat$slope[dat$sleep == 1])
cohen.d(dat$slope[dat$sleep == 0], dat$slope[dat$sleep == 1])
t.test(dat$slope[dat$sleep == 1], dat$slope[dat$sleep == 2])


se0 = sd(dat$slope[dat$sleep == 0])/sqrt(length(dat$slope[dat$sleep == 0]))
se1 = sd(dat$slope[dat$sleep == 1])/sqrt(length(dat$slope[dat$sleep == 1]))
se2 = sd(dat$slope[dat$sleep == 2])/sqrt(length(dat$slope[dat$sleep == 2]))

se_ = c(se0, se1, se2)

par(mar=c(5,1,4,7), lwd = 4)
k = barplot(c(mean0, mean1, mean2), col = 'white', border = sleep_cols, xlab = 'spectral slope', 
            xlim = c(-2.5,0), xaxt = 'n', ylab = '', main = '', width = 0.8, 
            space = 0.5, cex.lab = title_size, cex.names = label_size, lwd = 1.5,
            horiz = TRUE, names = c('', '', ''), las = 2)
axis(side = 4, at = k, labels = c('Wake', 'N1', 'N2'), las = 1, cex.axis = title_size, tick = FALSE)
beeswarm(dat$slope ~ dat$sleep, at = k,  
         pch = 21, col = sleep_cols, method = 'swarm', 
         horiz = TRUE, xlim = c(0,2.5), axes = FALSE, add = TRUE)
#axis(1, at = c(0,1,2), labels = c(0,1,2), cex.axis = 1.75, lwd = 1.5) 
axis(1, at = c(-2.5,-2,-1,0), labels = c('',-2,-1,0), cex.axis = label_size, lwd = 1.5) 

for (i in 1:length(mean_)){
  arrows(mean_[i] - se_[i], k[i], mean_[i] + se_[i], k[i], lwd = 2.5, length=0.1, 
         angle=90, code=3, 
         col = c(c(brewer.pal(n = 9, name = "YlGn")[6],brewer.pal(n = 11, name = "Spectral")[9], 
                   brewer.pal(n = 11, name = "BrBG")[10]))[i])
}


########################
### plot figure 3C
########################

# barplot of 1/f slopes
mean0 = mean(dat$slope[dat$insight == 0])
mean1 = mean(dat$slope[dat$insight == 1])

mean_ = c(mean0, mean1)

se0 = sd(dat$slope[dat$insight == 0])/sqrt(length(dat$slope[dat$insight == 0]))
se1 = sd(dat$slope[dat$insight == 1])/sqrt(length(dat$slope[dat$insight == 1]))

se_ = c(se0, se1)

title_size = 2.25
label_size = 1.75

par(mar=c(5,1,10,10), lwd = 3)
k = barplot(mean_, col = 'white', border = insight_cols, xlab = 'spectral slope', 
            xlim = c(-2.7,0), xaxt = 'n', ylab = '', main = '', width = 0.8, 
            space = 0.5, cex.lab = title_size, cex.names = label_size, lwd = 1.5,
            horiz = TRUE, names = c('', ''), las = 2)
axis(side = 4, at = k, labels = c('No insight', 'Insight'), las = 1, cex.axis = title_size, tick = FALSE)
beeswarm(dat$slope ~ dat$insight, at = k,  
         pch = 21, col = insight_cols, method = 'swarm', 
         horiz = TRUE, xlim = c(0,2.5), axes = FALSE, add = TRUE)
#axis(1, at = c(0,1,2), labels = c(0,1,2), cex.axis = 1.75, lwd = 1.5) 
axis(1, at = c(-2.5,-2,-1,0), labels = c('',-2,-1,0), cex.axis = label_size, lwd = 1.5) 

for (i in 1:length(mean_)){
  arrows(mean_[i] - se_[i], k[i], mean_[i] + se_[i], k[i], lwd = 2.5, length=0.1, 
         angle=90, code=3, 
         col = insight_cols[i])
}


