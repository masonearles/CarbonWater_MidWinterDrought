d.all$ID.unq = paste(d.all$ID, d.all$day, sep = "_")
d.insol$ID.unq = paste(d.insol$ID, d.insol$DAY, sep = "_")


d = merge(d.all, d.insol, by = 'ID.unq')


d.short = d[d$day %in% c(54, 84,129, 420, 466, 498),]
d.short$year = NA
d.short[d.short$day %in% c(54, 84, 129),]$year = 2014
d.short[d.short$day %in% c(420, 466, 498),]$year = 2015
p0 <- ggplot(d.short) + stat_summary(fun.data = 'mean_cl_boot', aes(x = VPD.zscore5, y = CONCAVG, color = SPECIES)) + facet_grid(. ~ year) + theme_bw()
p1 <- ggplot(d.short) + stat_summary(fun.data = 'mean_cl_boot', aes(x = P60, y = CONCAVG, color = species)) + facet_grid(. ~ year) + theme_bw()

grid.arrange(p0, p1)
