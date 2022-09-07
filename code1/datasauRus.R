# Generating the plots of the "datasauRus" dataset

# Load packages
library('datasauRus')
library('dplyr')
library('ggplot2')

# Generate scatterplots (and save them)
pdf(file = "datasauRusPlots.pdf", width = 7, height = 9)
ggplot(datasaurus_dozen, aes(x=x, y=y, colour=dataset))+
  geom_point()+
  coord_fixed(ratio=0.6)+
  theme(legend.position = "none")+
  facet_wrap(~dataset, ncol=3)
dev.off()

# Check the summary statistics (they are very similar)
  datasaurus_dozen %>% 
    group_by(dataset) %>% 
    summarize(
      mean_x    = mean(x),
      mean_y    = mean(y),
      std_dev_x = sd(x),
      std_dev_y = sd(y),
      corr_x_y  = cor(x, y)
    )


