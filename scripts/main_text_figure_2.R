###
# Reproducing main text Figure 2: Sensitivity of mortality to temperature and tropical cyclones 
# Reproducing supplementary Figure S4: Changes in point sensitivities over time for temperature-mortality relationships shown in Fig 2
# add two new countries: Mex and Colombia
# add new layout for fig 2 and point sensitivities in appendix

# part 1: output: main/figure_2.pdf
# part 2: output: supplementary/figure_s4.pdf 

rm(list = ls())
gc()

# remove anything in environment and load library
unique_packages <- c(
  "tidyverse",
  "fixest",
  "MetBrewer",
  "DescTools",
  "plotrix",
  "arrow",
  "data.table",
  "sf",
  "magrittr",
  "haven",
  "cowplot"
)

# Identify packages that are not installed
not_installed <- setdiff(unique_packages, rownames(installed.packages()))

# Install missing packages
if (length(not_installed) > 0) {
  install.packages(not_installed)
}

# Load all packages
lapply(unique_packages, library, character.only = TRUE)

source('scripts/functions/mkheatmap.R')
source('scripts/functions/polynomial_evaluation_functions.R')

# CONTINUE FROM HERE: 
# 1. Load pop weighted temperature distribution
load('data/mortality/usa/USTemperatureQuantiles_popweighted_v2.RDA') # pop weighted temperature distributions for evaluation 
load('data/mortality/mexico/MEXTemperatureQuantiles_popweighted.RDA')
load('data/mortality/colombia/COTemperatureQuantiles_popweighted.RDA')
load('data/mortality/eu/EUTemperatureQuantiles_popweighted.RDA')

# Now make main text plot: response functions, sensitivity at different temperatures, total sensitivity
# pdf(file="plots/combined/Figure_Mortality_Combined.pdf",height=9,width=10)
par(mar=c(3,5,2,2),mfrow=c(5,3))
# plot response functions, plotting confidence interval on first one
# US data
coef <- read_rds(file='data/mortality/usa/USDecade_polynomial_bootstrap_popweight.rds')
decades = seq(1970,2010,10)

us_response <- plotresponsefnc(coef=coef,xx=-15:40,center=us_mmt,xname="pw_tmean_poly",ylim=c(-0.004,0.01),ylab="log monthly mortality rate", decades=decades,
                               addtext="yes",textattx=-5,textatty=seq(0.002,0.0005,length.out=length(decades)),textlab=paste0(decades,"s"),
                               addhist="yes",histlb=-0.004,histh=0.002,histbreaks=us_mids,histcounts=us_pwdall,title="US mortality - temperature")

us_ps <- plotpointsensitivities(coef=coef,base_mmt = us_decade_mmt,xvals=c(-10,0,30,35),decades = decades,xname="pw_tmean_poly",ylim=c(0,0.007),title="point sensitivities")
us_ts <- plottotalsensitivity(coef=coef,pwd=us_pwd,xname="pw_tmean_poly",locations=us_mids,decades=decades,ylim=c(-0.00008,0.00002),addyrs="yes")

# Europe
coef <- read_rds('data/mortality/eu/EUDecade_polynomial_bootstrap_popweight.rds')
decades = seq(1990,2010,10)

eu_response <- plotresponsefnc(coef=coef,xx=-10:33, center = eu_mmt, xname="tmean_poly",ylim=c(-0.0008,0.002),ylab="log annual mortality rate", decades=decades,
                               addtext="yes",textattx=-5,textatty=seq(0.0008,0.00095,length.out=length(decades)),textlab=paste0(decades,"s"),
                               addhist="yes",histlb=-0.0008,histh=0.0004,histbreaks=eu_mids, histcounts=eu_pwdall,title="EU mortality - temperature")

eu_ps <- plotpointsensitivities(coef=coef,base_mmt = eu_decade_mmt, xvals=c(-10,0,25,30),decades = decades,xname="tmean_poly",ylim=c(-0.0005,0.0015),title="point sensitivities")
eu_ts <- plottotalsensitivity(coef=coef,pwd=eu_pwd,xname="tmean_poly",locations=eu_mids,decades=decades,ylim=c(-0.00001,0.00001),addyrs="yes")

# MEX data
coef <- read_rds(file='data/mortality/mexico/MEXDecade_polynomial_bootstrap_popweight.rds')
decades = seq(1990,2010,10)
mex_response <- plotresponsefnc(coef=coef,xx=0:40,center=mex_mmt,xname="pw_tmean_poly",ylim=c(-0.01,0.03),ylab="log monthly mortality rate", decades=decades,
                                addtext="yes",textattx=5,textatty=seq(0.03,0.025,length.out=length(decades)),textlab=paste0(decades,"s"),
                                addhist="yes",histlb=-0.01,histh=0.005,histbreaks=mex_mids,histcounts=mex_pwdall,title="Mexico mortality - temperature")
mex_ps <- plotpointsensitivities(coef=coef,base_mmt = mex_decade_mmt,xvals=c(5,15,30,35),decades = decades,xname="pw_tmean_poly",ylim=c(0,0.02),title="point sensitivities")
mex_ts <- plottotalsensitivity(coef=coef,pwd=mex_pwd,xname="pw_tmean_poly",locations=mex_mids,decades=decades,ylim=c(-0.0005, 0.0001),addyrs="yes")

# COL data
coef <- read_rds(file='data/mortality/colombia/CODecade_polynomial_bootstrap_popweight.rds')
decades = seq(1980,2010,10)
# plot response functions, plotting confidence interval on first one
co_response <- plotresponsefnc(coef=coef, xx=9:33, center=co_mmt, xname="pw_tmean_poly", ylim=c(-0.01,0.03), ylab="log monthly mortality rate",  decades=decades,
                               addtext="yes", textattx=30, textatty=seq(0.03,0.025,length.out=length(decades)), textlab=paste0(decades,"s"),
                               addhist="yes", histlb=-0.01, histh=0.005, histbreaks=co_mids, histcounts=co_pwdall, title="Colombia monthly mortality")
co_ps <- plotpointsensitivities(coef=coef, base_mmt = co_decade_mmt, xvals=c(10,15,20,30), decades = decades, xname="pw_tmean_poly", ylim=c(-0.01,0.02), title="point sensitivities")
co_ts <- plottotalsensitivity(coef=coef, pwd=co_pwd, xname="pw_tmean_poly", locations=co_mids, decades=decades, ylim=c(-0.0009, 0.0001),addyrs="yes")

# TC- US mortality last panel
dt <- read_rds('data/mortality/tropical_cyclone/boots_n1000_period_model_est.rds')
pernames <- c("reference_period","time_4bins3","time_4bins4") #our three periods
colz=rev(met.brewer("Homer1",3))
att=seq(190,220,10 ) #where to put the dot and whisker plots
plot(1,type="n",xlim=c(0,220),ylim=c(0,10),las=1,axes=F,xlab="years since storm",ylab="deaths per 100k",main="US mortality - cyclones")
axis(1,at=c(0,seq(1,15,2))*12,c(0,seq(1,15,2)))
axis(2,las=1)
for (i in 1:length(pernames)) {
  df <- dt %>% filter(interacted==pernames[i]) %>% filter(maxs_lag<=180)
  toplot <- df %>% group_by(maxs_lag) %>%   
    summarise(median=median(cumulative_est),
              mean=mean(cumulative_est),
              cihi=quantile(cumulative_est,probs = 0.975),
              cilo=quantile(cumulative_est,probs = 0.025))
  if (i==1) {polygon(c(toplot$maxs_lag,rev(toplot$maxs_lag)),c(toplot$cilo,rev(toplot$cihi)),col=alpha(colz[i],0.3),border=NA)}
  lines(toplot$maxs_lag,toplot$median,col=colz[i],lwd=2)
  #box and whisker at the end
  segments(att[i],toplot$cilo[toplot$maxs_lag==180],att[i],toplot$cihi[toplot$maxs_lag==180],col=colz[i])
  points(att[i],toplot$median[toplot$maxs_lag==180],pch=19,col=colz[i])
}

dev.off()

# Main text Figure 2 (Version 2): 
# In the second version, we include both mex and colombia results
misc_dict <- list(
  "US" = list(
    decade_mmt = us_decade_mmt,
    pwd = us_pwd,
    pwdall = us_pwdall,
    mids = us_mids,
    mmt = us_mmt,
    temps = us_temps,
    decades = seq(1970,2010,10),
    response = us_response, 
    ps = us_ps,
    ts = us_ts
  ),
  "EU" = list(
    decade_mmt = eu_decade_mmt,
    pwd = eu_pwd,
    pwdall = t(eu_pwdall),
    mids = eu_mids,
    mmt = eu_mmt,
    temps = eu_temps,
    decades = seq(1990,2010,10),
    response = eu_response, 
    ps = eu_ps,
    ts = eu_ts
  ),
  "MEX" = list(
    decade_mmt = mex_decade_mmt,
    pwd = mex_pwd,
    pwdall = mex_pwdall,
    mids = mex_mids,
    mmt = mex_mmt,
    temps = mex_temps,
    decades = seq(1990,2010,10),
    response = mex_response, 
    ps = mex_ps,
    ts = mex_ts
  ),
  "CO" = list(
    decade_mmt = co_decade_mmt,
    pwd = co_pwd,
    pwdall = co_pwdall,
    mids = co_mids,
    mmt = co_mmt,
    temps = co_temps,
    decades = seq(1980,2010,10),
    response = co_response, 
    ps = co_ps,
    ts = co_ts
  )
)

# combine figure 2 and point sensitivities
plot_figure_2 <- function(){
  country_names <- c("US", "MEX", "CO", "EU")
  
  out_plots <- lapply(country_names, function(country_name){
    
    # country_name <- "MEX"
    needed <- misc_dict[[country_name]]
    
    decades <- needed$decades
    color_palette <- rev(met.brewer("Homer1",length(decades)))
    color_palette_names <- paste0(decades, "s")
    
    # prepare response data
    response <- needed$response %>% 
      mutate(decade = paste0(decade, "s"))
    
    if (country_name == "EU") {
      response <- response %>% mutate(across(starts_with("pred"), ~.x*12))
      needed$ts <- needed$ts %>% mutate(across(starts_with("p"), ~.x*12))
      needed$ps <- needed$ps %>% mutate(across(starts_with("p"), ~.x*12))
    }
    
    # create breaks for response
    original_breaks <- pretty(c(0, response$pred_975), n = 4)
    
    # prepare exposure data
    exposure <- as.data.frame(t(needed$pwdall))
    names(exposure) <- c("pwd")
    exposure$bins <- needed$mids
    max_response <- max(response$pred_50, na.rm = TRUE)
    max_exposure <- max(exposure$pwd, na.rm = TRUE)
    
    # adjust height of the exposure to mean of ci higher end
    exposure$scaled_pwd <- exposure$pwd / max_exposure * mean(response$pred_975) * 0.5
    exposure %<>% mutate(bins = ceiling(bins)) %>% filter(bins %in% unique(response$bins)) 
    
    # yoffset
    yoffset <- max(exposure$scaled_pwd) + max(exposure$scaled_pwd)*0.2 
    response %<>% mutate(across(starts_with("pred"), ~.x+yoffset))
    
    
    # new y breaks
    response_exposure_plot <- ggplot(response, aes(x = bins, y = pred_50, color = decade, fill = decade)) +
      geom_ribbon(
        data = response %>% filter(decade == color_palette_names[1]), 
        aes(ymin = pred_025, ymax = pred_975), alpha = 0.3, color = NA
      ) + 
      geom_line() + 
      geom_col(
        data = exposure, 
        aes(x = bins, y = scaled_pwd), 
        fill = "grey70", 
        color = "grey70",
        width = 1
      ) +
      scale_y_continuous(
        breaks = original_breaks + yoffset,
        labels = function(x) x-yoffset
      ) +
      geom_hline(yintercept = 0 + yoffset, linetype = "dashed") +
      scale_color_manual(values = setNames(color_palette, color_palette_names)) +
      scale_fill_manual(values = setNames(color_palette, color_palette_names)) +
      theme_classic() +
      labs(x = "daily temperature (°C)", y = "log monthly mortality rate", title = paste0(country_name, " mortality-temperature")) +
      # guides(fill = "none") + 
      theme(
        # general design
        plot.title = element_text(size = 10, color = "black", face = "bold"),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10),
        axis.line = element_line(color = 'black',linewidth  = 0.35),
        axis.ticks = element_line(colour = "black",linewidth = 0.35),
        # axis title and text
        axis.title.y = if (country_name == "US") element_text(margin=margin(t=0,r=10,b=0,l=0)) else element_blank(),
        axis.text.y = element_text(color="black"), 
        axis.title.x = element_text(margin=margin(t=0,r=10,b=0,l=0)),
        axis.text.x = element_text(color="black"),
        # axis.line.x = element_blank(),
        # axis.ticks.x = element_blank(),
        # legend settings
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.8),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'lines'), 
        legend.key.spacing.y = unit(0.1, 'lines'),
        legend.background = element_blank(),
        # plot margins
        plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
      )
    
    # total sensitivities
    ts_plot <- needed$ts %>% 
      mutate(decade = paste0(decade, "s")) %>% 
      ggplot(aes(x = factor(1), y = p50, color = decade)) +  
      geom_errorbar(aes(ymin = p025, ymax = p975), width = 0.1, position = position_dodge(width = 1)) + 
      geom_point(position = position_dodge(width = 1)) +  
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_manual(values = setNames(color_palette, color_palette_names)) + 
      scale_y_continuous(position = "left") + 
      theme_classic() +
      labs(y = "E[dy/dT]", title = paste0(country_name, " total sensitivities")) +
      theme(
        # general design
        plot.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10),
        axis.line = element_line(color = 'black',linewidth  = 0.35),
        axis.ticks = element_line(colour = "black",linewidth = 0.35),
        # axis title and text
        axis.title.y = if (country_name == "US") element_text(margin=margin(t=0,r=5,b=0,l=0)) else element_blank(),
        axis.text.y = element_text(color="black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        # legend settings
        legend.position ="none",
        # legend.title = element_blank(),
        # legend.key.size = unit(0.5, 'lines'), 
        # legend.key.spacing.y = unit(0.1, 'lines'),
        # legend.background = element_blank(),
        # plot margins
        plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
      )
    
    if (country_name %in% c("US","CO")) {
      ylab.text <- "effect of additional day, \n relative to period-specific MMT (°C day)"
    } else {
      ylab.text <- " \n  \n "
    }
    
    # plot point sensitivities
    ps_plot <- needed$ps %>% 
      mutate(eval_temp = factor(eval_temp, levels = unique(needed$ps$eval_temp))) %>% 
      mutate(decade = paste0(decade, "s")) %>% 
      ggplot(aes(x = eval_temp, y = p50, color = decade)) +  
      geom_errorbar(aes(ymin = p025, ymax = p975), width = 0.1, position = position_dodge(width = 0.3)) + 
      geom_point(position = position_dodge(width = 0.3)) +  
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(
        xintercept = seq(1.5, length(unique(needed$ps$eval_temp)) - 0.5), 
        linetype = "dashed", 
        color = "gray50", 
        linewidth = 0.3
      ) +
      geom_text(
        data = ~ distinct(., eval_temp) %>% mutate(y_pos = max(needed$ps$p975, na.rm = TRUE)),
        aes(x = eval_temp, y = y_pos * 1.05, label = eval_temp),
        color = "black",
        size = 3.5,
        vjust = 0
      ) +
      scale_color_manual(values = setNames(color_palette, color_palette_names)) +
      labs(y = ylab.text, title = paste0(country_name, " point sensitivities")) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 10, color = "black", face = "bold"),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10),
        axis.line = element_line(color = 'black', linewidth = 0.35),
        axis.ticks = element_line(colour = "black", linewidth = 0.35),
        axis.title.y = element_text(color="black"),
        axis.text.y = element_text(color="black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),  # Hide default x-axis text
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'lines'), 
        legend.key.spacing.y = unit(0.1, 'lines'),
        legend.background = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
      ) +
      coord_cartesian(clip = "off")
    
    return(list(response=response_exposure_plot, ts = ts_plot, ps = ps_plot))
  })
  
  # Extract all response plots
  response_plots <- lapply(out_plots, function(x) x$response)
  ts_plots <- lapply(out_plots, function(x) x$ts)
  ps_plots <- lapply(out_plots, function(x) x$ps)
  
  # plot_grid(plotlist = response_plots, ncol = 5, nrow = 1, align = "v")
  # tropical cyclone 
  dt <- read_rds('data/mortality/tropical_cyclone/boots_n1000_period_model_est.rds')
  pernames <- c("reference_period", "time_4bins3", "time_4bins4")
  colz <- rev(met.brewer("Homer1", 3))
  
  # Prepare summary data for ggplot
  df_all <- lapply(seq_along(pernames), function(i) {
    df <- dt %>% filter(interacted == pernames[i], maxs_lag <= 180)
    df %>%
      group_by(maxs_lag) %>%
      summarise(
        median = median(cumulative_est),
        cihi = quantile(cumulative_est, probs = 0.975),
        cilo = quantile(cumulative_est, probs = 0.025)
      ) %>%
      mutate(period = pernames[i], color = colz[i])
  }) %>% bind_rows()
  
  # tropical cyclone separate
  tc_response <- ggplot(df_all, aes(x = maxs_lag, y = median, color = period, fill = period)) +
    geom_ribbon(
      data = df_all %>% filter(period == "reference_period"), 
      aes(ymin = cilo, ymax = cihi), alpha = 0.3, color = NA) +
    geom_line() +
    scale_color_manual(
      values = setNames(colz, pernames),
      labels = c("reference_period" = "1952-1973", "time_4bins3" = "1975-1994", "time_4bins4" = "1995-2015")
    ) +
    scale_fill_manual(
      values = setNames(colz, pernames),
      labels = c("reference_period" = "1952-1973", "time_4bins3" = "1975-1994", "time_4bins4" = "1995-2015")
    ) +
    scale_x_continuous(
      breaks = c(0, seq(12, 180, by = 24)),
      labels = c(0, seq(1, 15, by = 2)),
      # expand = c(0, 0),
      name = "years since storm"
    ) +
    labs(ylab = "deaths per 100K", title = "US mortality-cyclones") +
    ylab("deaths per 100k") +
    theme_classic() +
    theme(
      # general design
      axis.text = element_text(size=10),
      axis.title = element_text(size=10),
      plot.title = element_text(size=10, color = "black", face = "bold"),
      axis.line = element_line(color = 'black',linewidth  = 0.35),
      axis.ticks = element_line(colour = "black",linewidth = 0.35),
      # axis title and text
      axis.title.y = element_text(margin=margin(t=0,r=10,b=0,l=10)),
      axis.text.y = element_text(color="black"), 
      axis.title.x = element_text(margin=margin(t=0,r=10,b=0,l=0)),
      axis.text.x = element_text(color="black"),
      # axis.line.x = element_blank(),
      # axis.ticks.x = element_blank(),
      # legend settings
      legend.position = "inside",
      legend.position.inside = c(0.3, 0.8),
      legend.title = element_blank(),
      legend.key.size = unit(0.5, 'lines'), 
      legend.key.spacing.y = unit(0.1, 'lines'),
      legend.background = element_blank(),
      # plot margins
      plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
    )
  
  # final whisker + dot
  tc_ts <- df_all %>% filter(maxs_lag == 180) %>% 
    ggplot(aes(x = factor(maxs_lag), y = median, color = period)) +  
    geom_errorbar(aes(ymin = cilo, ymax = cihi), width = 0.1, position = position_dodge(width = 1)) + 
    geom_point(position = position_dodge(width = 1)) +  
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(
      values = setNames(colz, pernames),
      labels = c("reference_period" = "1952-1973", "time_4bins3" = "1975-1994", "time_4bins4" = "1995-2015")
    ) +
    scale_y_continuous(position = "left") + 
    theme_classic() +
    theme(
      # general design
      plot.title = element_blank(),
      axis.text = element_text(size=10),
      axis.title = element_text(size=10),
      axis.line = element_line(color = 'black',linewidth  = 0.35),
      axis.ticks = element_line(colour = "black",linewidth = 0.35),
      # axis title and text
      axis.title.y = element_blank(),
      axis.text.y = element_text(color="black"), 
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      # axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      # legend settings
      legend.position ="none",
      # legend.title = element_blank(),
      # legend.key.size = unit(0.5, 'lines'), 
      # legend.key.spacing.y = unit(0.1, 'lines'),
      # legend.background = element_blank(),
      # plot margins
      plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
    )
  
  response_plots[[5]] <- tc_response
  ts_plots[[5]] <- tc_ts
  final_response_ts <- plot_grid(plotlist = c(response_plots, ts_plots), ncol = 5, nrow = 2, align = "hv", axis = "tlr", rel_heights = c(1, 0.5))
  
  # plot point sensitivities together
  final_ps <- plot_grid(plotlist = ps_plots, ncol=2, nrow = 2, axis = "tblr", align = "v")
  
  return(list(response_ts = final_response_ts, ps = final_ps))
}

final_plots <- plot_figure_2()

ggsave(
  filename = "fig/main/figure_2.pdf",
  plot = final_plots$response_ts,
  height = 6,
  width = 15
)

ggsave(
  filename = "fig/supplementary/figure_s4.pdf",
  plot = final_plots$ps,
  height = 8,
  width = 10
)
