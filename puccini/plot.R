library(ggplot2)
library(scales)

cores <- 1:11
comptime <- c(48506,31300,23763,19785,17065,15127,13504,13154,12142,11196,10228)
writing <- c(588,425,447,423,422,423,422,426,425,430,424)

speedup <- comptime[1]/comptime
plot(cores,comptime)
plot(cores,speedup)


##ggplot dataframe
data <- data.frame(Time = c(writing,comptime), Legend = c(rep("Writing",11),rep("Computation",11)),
                   Cores = cores)
data$Legend <- factor(data$Legend,levels = c("Writing","Computation"))

plot_time <- ggplot(data, aes(x=Cores, y=Time, fill=Legend)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("#51CBD1","#F27C5F"))


data2 <- data.frame(Speedup = speedup, Cores = cores)
plot_speedup <- ggplot(data2, aes(x=Cores,y=Speedup)) + geom_line(linetype = "dashed",color = "#5FA9F2") + geom_point(color = "red") + 
  scale_y_continuous(breaks=cores, labels=cores, limits = c(1,11)) +
  scale_x_continuous(breaks=cores, labels=cores, limits = c(1,11)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.3, color = "#62C658")

ggsave(plot_time, file="plot_time.png")
ggsave(plot_speedup, file="plot_speedup.png")                     
