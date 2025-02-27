library(ggh4x)

#### plot 

## mean Var Coeff per hpi per window
ggvarcoeff <- ggplot(cov.hpi.sum.mean,
                     aes(x=as.factor(Time), y=varcoeff_coverage)) +
  geom_col(aes(fill=cell_line),
           position = 'dodge2', 
           alpha = 0.7) +
  scale_color_npg(alpha = 0.7) +
  scale_fill_npg(alpha = 0.7) +
  theme_bw() + 
  facet_nested(rows=vars(strand))

ggsave(paste0(outdir, '/coverage_sum_varcoeff.jpg'), ggvarcoeff, width=9, height = 7, limitsize = F)

## mean Var Coeff per hpi per genome average coverage
ggvarcoeff <- ggplot(norm_cov_summary.hpi,
                     aes(x=as.factor(Time), y=varcoeff_coverage)) +
  geom_col(aes(fill=cell_line),
           position = 'dodge2', 
           alpha = 0.7) +
  scale_color_npg(alpha = 0.7) +
  scale_fill_npg(alpha = 0.7) +
  theme_bw() 
# + facet_nested(rows=vars(strand))

ggsave(paste0(outdir, '/coverage_sum_varcoeff.hpi.jpg'), ggvarcoeff, width=9, height = 6, limitsize = F)


## mean and sd per hpi per genome average coverage
ggmeansd <- ggplot(norm_cov_summary.hpi,
                   aes(x=cell_line, y=mean_coverage.hpi, group=cell_line)) +
  geom_col(aes(fill=cell_line),
           position = position_dodge2(), 
           alpha = 0.7) +
  geom_errorbar(aes(fill = cell_line, 
                    ymin = mean_coverage.hpi + sd_coverage.hpi, 
                    ymax = mean_coverage.hpi - sd_coverage.hpi),
                position = position_dodge2(width = 1, padding = 0.5)) +
  scale_color_npg(alpha = 0.7) +
  scale_fill_npg(alpha = 0.7) +
  theme_bw() + 
  facet_nested_wrap(~as.factor(Time), nrow=1, scales = 'free')

ggsave(paste0(outdir, '/coverage_sum_meanSD.hpi.2.jpg'), ggmeansd, width=12, height = 5, limitsize = F)

## mean and sd per hpi per genome average coverage
ggcov <- ggplot(norm_cov_summary,
                aes(x=cell_line, y=average_coverage, group=cell_line)) +
  geom_col(aes(fill=cell_line),
           position = position_dodge2(), 
           alpha = 0.7) +
  scale_color_npg(alpha = 0.7) +
  scale_fill_npg(alpha = 0.7) +
  theme_bw() + 
  facet_nested_wrap(~as.factor(Time), nrow=1, scales = 'free')

ggsave(paste0(outdir, '/coverage_sum.hpi.2.jpg'), ggcov, width=17, height = 5, limitsize = F)


plot.data <- win.cov.hpi.sum
plot.data[,mean_coverage.hpi := fifelse(strand == '+', mean_coverage.hpi, mean_coverage.hpi * -1)]

ggp <- ggplot(
  plot.data, # [cell_line == 'PK-15',]
  aes(x=window_start, y=mean_coverage.hpi, fill=strand)) +
  geom_area() +
  coord_flip() + 
  theme_bw() +
  facet_nested(cols=vars(as.factor(Time), cell_line))

ggsave(paste0(outdir, '/win.cov.sum.jpg'), ggp, width=30, height = 12, limitsize = F)



