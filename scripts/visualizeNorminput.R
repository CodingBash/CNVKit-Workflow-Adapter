#
# This script is for visualizing the normal segments lengthXratio to see how often small fragments can be artifacts
# Practically, this script is for determining a rough cutoff of removing artifacts based on segmedian length
#

length_threshold=1.1e7
# Determine cutoff
plot_ylim <- 0.5
plot_xlim <- 150000000
oNI <- norminput[order(norminput$length), ]
plot(oNI$length, oNI$segmedian, ylim = c(-plot_ylim, plot_ylim), xlim = c(0, plot_xlim))
dev.off()
nNI <- norminput[norminput$length > length_threshold,]
plot(nNI$length, nNI$segmedian, ylim = c(-plot_ylim, plot_ylim), xlim = c(0, plot_xlim))

