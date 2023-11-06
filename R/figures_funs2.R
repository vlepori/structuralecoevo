# require(ggplot2)
# require(patchwork)


ssq <- scales::trans_new("safesq",
                         function(x) x^2,
                         function(y) sign(y)*sqrt(abs(y)))

slog <- scales::trans_new("safelog",
                         function(x) exp(x),
                         function(y) log(abs(y)))
                        # function(y) sign(y)*sqrt(abs(y)))


#*************************************************************************
# FIGURE 1 ---------------------------------------------------------------
#*************************************************************************

# args...log_scale
figure1 = function(df_sim){
  df_sim=filter(df_sim, time > 0.2)
  
  
  p1_1=
    df_sim %>% 
    select(time, x, Nstars) %>% 
    tidyr::unnest(cols=c(x, Nstars)) %>% 
    # filter(time>0.2) %>%
    ggplot()+
    geom_vline(xintercept = pull(filter(df_sim, branching==1),time), linetype = "dotted", color = "grey50")+
    geom_point(aes(x=time,y = x), size = .3)+
    theme_bw()+
    scale_x_log10()+
    # scale_y_continuous(name = expression(paste("Niche position (",mu,")")),sec.axis = sec_axis(~(.), name= " "))+
    scale_y_continuous(name = "Niche position")+
    # scale_color_viridis(name = expression(paste("Biomass")))+
    guides(color = guide_colorbar(title.vjust = .8))+ # manually adjust alignment
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.x = element_blank(),
          legend.direction="horizontal",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size=14),
          axis.ticks.y.right = element_blank(),
          axis.title.y.right = element_text(size=14),
          axis.text =element_text(size=12))  
  
  
  colorkey = c("eta" = "red", "theta" = "orange")
  # names(colorkey) = c(expression(mu),expression(theta))
  
  p1_2_2 = ggplot(df_sim, aes(x = time))+
    geom_vline(xintercept = pull(filter(df_sim, branching==1),time), linetype = "dotted", color = "grey50")+
    geom_line(aes(group = S, y = eta,col="eta")
    )+
    geom_line(aes(group = S, y = theta,col="theta")
    )+
    theme_bw()+
    # labs(y = expression(paste("Fitness difference (", theta,")")))+
    labs(y = "Structural resistance\n & fitness difference")+
    scale_x_log10(name = "Evolutionary time")+
    # scale_y_continuous(sec.axis = sec_axis(~(.), name=expression(paste("Robustness (",eta,")"))))+
    scale_color_manual(values = colorkey, labels = c(expression(eta),expression(theta)))+
    # theme(legend.position = c(log(mean(range(df_sim$time,na.rm = T))),mean(range(df_sim$eta, na.rm = T))), 
    theme(legend.position = c(.1,.5), 
                
          strip.background = element_rect(fill = "white"),
          strip.text.y = element_text(colour = "black",angle = 0),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=14),
          axis.title.y = element_text(size=14),
          # axis.text.y = element_text(color="orange"),
          # axis.text.y.right = element_text(color="red"),
          axis.title.y.right = element_text(size=14),
          axis.title=element_text(size=14),
          axis.text =element_text(size=10))
  
  # return(p1_2_2)
  
  # browser()
  
  p1_2_3 = ggplot(df_sim, aes(x = time))+
    geom_vline(xintercept = pull(filter(df_sim, branching==1),time), linetype = "dotted", color = "grey50")+
    geom_line(aes(group = S, y = omega_rel),col="darkgreen")+
    # geom_line(aes(group = S, y = sum_Y),col="blue"
    # )+
    theme_bw()+
    labs(y = "Structural\nniche difference")+
    # labs(y = expression(paste("     Structural\nniche difference (", hat(Omega),")")))+
    scale_x_log10(name = "Evolutionary time")+
    # scale_y_continuous(sec.axis = sec_axis(~(.), name=expression(paste("Total relative yield"))))+
    theme(legend.position = "none", 
          strip.background = element_rect(fill = "white"),
          strip.text.y = element_text(colour = "black",angle = 0),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.text.y = element_text(color="darkgreen"),
          # axis.text.y.right = element_text(color="blue"),
          axis.title.y = element_text(size=14),
          axis.title.y.right = element_text(size=13),
          axis.text =element_text(size=10))
  

  p1_2 <- p1_1  / p1_2_3 / p1_2_2 + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold')) # + p1_2_4 + plot_layout(ncol=2)
  # fig1 <- / p1_2
  # + plot_layout(heights = c(2,6)) 
  # + plot_annotation(tag_levels = (c('A',"1")))
  
  return(p1_2)
  
}


#*************************************************************************
# FIGURE 2 ---------------------------------------------------------------
#*************************************************************************

getfront = function(df){
  
  lseq <- function(from, to, length.out){exp(seq(log(from), log(to), length.out = length.out))}
  
  out = df[FALSE,]
  
  grd = lseq(from = min(df$theta,na.rm=T), to = max(df$theta,na.rm=T), length.out = 20)
  for (i in 2:20) {
    tmp = filter(df, between(theta,grd[i-1],grd[i]))%>% arrange(desc(omega_rel))  %>% head(n=1L) 
    out = rbind(out,tmp)
  }
  out
}

# getfront(filter(mcA,S==7)) %>% ggplot(aes(x=omega_rel,y=theta))+geom_line()+facet_wrap(~S)



figure2 = function(sim_df,sim_mc,fronts = FALSE,facet = TRUE){

  fronts_df = sim_mc[FALSE,]
  
  for(s in min(sim_mc$S):max(sim_mc$S)){
    if(s>1){
    fronts_df=rbind(fronts_df,getfront(filter(sim_mc,S==s)))
    }
  }
  
  # browser()
 sim_dt =
  sim_df %>% 
    filter(S > 1, !is.na(theta)) %>%
    group_by(S) %>% 
    mutate(rtime =  1:n()) %>% 
    ungroup() 
    
 endpoints = sim_dt %>%  group_by(S) %>% filter(time == max(time)) %>% ungroup()
 
  fig2 = 
    ggplot(sim_dt)+
    # {if(3 %in% sim_mc$id){geom_point(data = filter(sim_mc,id==3 & S>1), aes(x = omega_rel, y = theta),size = 1, col = "grey90", shape=16, alpha = .3)}}+
    geom_point(data = filter(sim_mc,id==2 & S>1), aes(x = omega_rel, y = theta),size = 1, col = "grey80", shape=16, alpha = .3)+ #light
    geom_point(data = filter(sim_mc,id==1 & S>1), aes(x = omega_rel, y = theta), size = 1, col = "grey50", shape=16, alpha = .3)+ #dark
    # {if(fronts){geom_line(data = fronts_df, aes(x = omega_rel, y = theta), col = "black")}}+
    geom_line(data = fronts_df, aes(x = omega_rel, y = theta), col = "black",linetype = 2,alpha=.8)+
    # geom_point(aes(x = omega_rel, y = theta, col = rtime), size = .5) +
    geom_line(aes(x = omega_rel, y = theta), col="darkgreen", size = .9) +
    geom_point(data=endpoints, aes(x = omega_rel, y = theta), col="darkgreen", size = 2.3) +
    {if(facet)facet_wrap( ~ S, scales = "fixed",labeller = as_labeller(function(x)paste("S =",x)), nrow = 1) }+
    scale_x_log10()+
    theme_bw() +
    # theme(strip.background = element_rect(fill = "white"),
    #       legend.position = c(.9, 0.1), # hacky
    #       legend.justification = c(1, 0))+
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          legend.position = "none",
          axis.title = element_blank(),
          # axis.title=element_text(size=15),
          axis.text =element_text(size=10))+
    labs(color = "Time",
         x = expression(paste("Standardized niche difference (", hat(Omega),")")), 
         y = expression(paste("Fitness difference (", theta,")")))+
    # scale_colour_continuous() +
    # scale_color_viridis_c(trans = "log", breaks = round(log_seq(min(simA$time[simA$time>0]),max(simA$time),length.out = 5)))
    scale_color_viridis_c(guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                                 label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                                 label.theme = element_text(angle = 90)),
                          option="plasma", trans = "log", breaks = 10^c(1:9))
  
  return(fig2)
}


# figure2b = function(sim_df,sim_mc){
#   fig2 = 
#     sim_df %>% 
#     filter(S > 1) %>%
#     ggplot()+
#     {if(3 %in% sim_mc$id){geom_point(data = filter(sim_mc,id==3 & S>1), aes(x = omega_rel, y = eta),size = 1, col = "grey90", shape=16, alpha = .3)}}+
#     geom_point(data = filter(sim_mc,id==2 & S>1), aes(x = omega_rel, y = eta),size = 1, col = "grey80", shape=16, alpha = .3)+ #light
#     geom_point(data = filter(sim_mc,id==1 & S>1), aes(x = omega_rel, y = eta), size = 1, col = "grey50", shape=16, alpha = .3)+ #dark
#     geom_point(aes(x = omega_rel, y = eta, col = time), size = .5) +
#     facet_wrap( ~ S, scales = "free_y",labeller = as_labeller(function(x)paste("S =",x))) +
#     scale_x_log10()+
#     theme_bw() +
#     # theme(strip.background = element_rect(fill = "white"),
#     #       legend.position = c(.9, 0.1), # hacky
#     #       legend.justification = c(1, 0))+
#     theme(strip.background = element_rect(fill = "white"),
#           legend.position = "none",
#           axis.title=element_text(size=15),
#           axis.text =element_text(size=10))+
#     labs(color = "Time",
#          x = expression(paste("Standardized niche difference (", hat(Omega),")")), 
#          y = expression(paste("Angle to feasibility domain border (", eta,")")))+
#     # scale_colour_continuous() +
#     # scale_color_viridis_c(trans = "log", breaks = round(log_seq(min(simA$time[simA$time>0]),max(simA$time),length.out = 5)))
#     scale_color_viridis_c(guide = guide_colorbar(direction = "horizontal", title.position = "top",
#                                                  label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
#                                                  label.theme = element_text(angle = 90)),
#                           option="plasma", trans = "log", breaks = 10^c(1:9))
#   
#   return(fig2)
# }

#*************************************************************************************

figure3 = function(f3){
  f3$S <- sapply(f3$rs, FUN = length)
  f3$omega <- sapply(FUN = tb$Omega.L1, X = f3$alpha) # log big omega
  f3$omega_rel <- 10^(f3$omega/f3$S)   # linear smalle omega
  
  out=
  f3 %>% 
    filter(!S == 1) %>% 
    filter(!(S==6 & sigmar == 3.1)) %>% 
    ggplot()+
    geom_tile(aes(x = sigmar, y = S, fill = exp(omega_rel)))+
    geom_text(data = filter(f3, S>1 & stability == "css"), aes(x = sigmar, y = S), size = 3, color = "white", label = "ESC")+
    geom_text(data = filter(f3, S>1 & stability == "branching"), aes(x = sigmar, y = S), size = 3, color = "white", label = " ")+
    scale_x_continuous(breaks = unique(f3$sigmar)[seq(1,24,by=2)])+
    scale_fill_viridis_c(option="plasma" )+
    labs(x = expression(paste("Resource width (", sigma[R],")")),
         fill = expression(paste("Standardized\nniche difference (", hat(Omega),")")), 
         y = "Morphs at the singular strategy")+
    coord_fixed(ratio = (max(f3$sigmar)-min(f3$sigmar))/(max(f3$S,na.rm = T)-min(f3$S,na.rm = T)))+ # force square
    cowplot::theme_minimal_hgrid()+
    guides(fill = guide_colorbar(title.vjust = .8, title.hjust = 1))+ # manually adjust alignement
    theme(legend.position = "top", 
          axis.title=element_text(size=15),
          axis.text =element_text(size=10),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()
    )
  
  return(out)
}

figure3b = function(f3){
  f3$S <- sapply(f3$rs, FUN = length)
  f3$omega <- sapply(FUN = tb$Omega.L1, X = f3$alpha) # log big omega
  f3$omega_rel <- 10^(f3$omega/f3$S)   # linear smalle omega
  
  out=
    f3 %>% 
    filter(!S == 1) %>% 
    # filter(!(S==6 & sigmar == 3.1)) %>% 
    ggplot()+
    geom_tile(aes(x = sigma, y = S, fill = exp(omega_rel)))+
    geom_text(data = filter(f3, S>1 & stability == "css"), aes(x = sigma, y = S), size = 3, color = "white", label = "ESC")+
    geom_text(data = filter(f3, S>1 & stability == "branching"), aes(x = sigma, y = S), size = 3, color = "white", label = " ")+
    scale_x_continuous(breaks = unique(f3$sigma))+
    scale_fill_viridis_c(option="plasma")+
    labs(x = expression(paste("Niche width (", sigma[alpha],")")),
         fill = expression(paste("Standardized niche difference (", hat(Omega),")")),
         y = "Morphs at the singular strategy")+
    coord_fixed(ratio = (max(f3$sigma)-min(f3$sigma))/(max(f3$S,na.rm = T)-min(f3$S,na.rm = T)))+ # force square
    cowplot::theme_minimal_hgrid()+
    # guides(fill = guide_colorbar(title.vjust = .8, title.hjust = 1))+ # manually adjust alignement
    theme(legend.position = "right", 
          axis.title=element_text(size=15),
          axis.text =element_text(size=10),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()
    )
  
  return(out)
}



#*************************************************************************************


# figure4 = function(sim_df,mc_df,thin=1){
#   
#   # dfN <-
#   #   left_join(
#   #     sim_df %>%
#   #       group_by(S) %>%
#   #       summarise(maxT = max(time),
#   #                 minT = min(time)),
#   #     mc_df %>%
#   #       filter(id == 2) %>%
#   #       group_by(S) %>%
#   #       summarise(maxN = max(sum_n_stars)),
#   #     by = "S"
#   #   )
#   
#   
#   dfN2 <-
#     mc_df %>% 
#     filter(id == 2) %>% 
#     left_join(sim_df %>% group_by(S) %>% mutate(logt = ifelse(time>0,log(time),NA)) %>% summarise(meanT = exp(mean(range(logt,na.rm=T)))), # midrange
#               by = "S")
#   
#   
#   ess = sim_df %>%  group_by(S) %>% filter(time == max(time)) %>% ungroup
#   perc = rep(NA,max(ess$S))
#   for (i in 1:max(ess$S)) {
#     tmp = filter(mc_df, S==i & id == 2) %>% pull(sum_n_stars)> ess$sum_n_stars[i] 
#     perc[i] = (1-round(mean(tmp,na.rm = T),2))
#   }
#   
#   # print(EPRC)
#   df3 = data.frame(x = dfN2 %>% with(unique(meanT)), l = perc)
#   
#   # browser()
#   out= 
#     ggplot(filter(dfN2, row_number() %% thin == 0), aes(x = log10(meanT), y = sum_n_stars))+
#     # geom_segment(data = dfN, aes(x = log10(minT), xend = log10(maxT), y = maxN, yend = maxN), color = "grey80")+
#     geom_vline(xintercept = log10(pull(filter(sim_df, branching==1),time)), linetype = "solid", color = "grey80")+
#     geom_quasirandom(aes(group = S), alpha = .5, color = "grey80", shape=16)+
#     geom_line(data=sim_df,aes(x = log10(time),y=sum_n_stars), size = 1, color = "red")+
#     geom_text(data=dfN2,aes(x=log10(meanT), label = paste("S = ",S)), y = exp(4.3), size = 4.5)+
#     geom_text(data=df3,aes(x=log10(x), label = l), y = exp(4.2), size = 4.5)+
#     labs(x="Evolutionary time", y = expression(paste("Total abundance")))+
#     # labs(x="Evolutionary time", y = expression(paste("Total abundance (", sum(N[i]),")")))+
#     theme_bw() +
#     scale_y_continuous(trans = slog, limits = c(0,4.3), breaks = 0:5)+
#     theme(panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.title=element_text(size=15),
#           axis.text =element_text(size=10))
#   return(out)
# }
# 
# figure4b = function(sim_df,mc_df){
#   
#   # dfN <-
#   #   left_join(
#   #     sim_df %>%
#   #       group_by(S) %>%
#   #       summarise(maxT = max(time),
#   #                 minT = min(time)),
#   #     mc_df %>%
#   #       filter(id == 1) %>%
#   #       group_by(S) %>%
#   #       summarise(maxN = max(sum_n_stars)),
#   #     by = "S"
#   #   )
#   
#   
#   dfN2 <-
#     mc_df %>% 
#     filter(id == 1) %>% 
#     left_join(sim_df %>% group_by(S) %>% mutate(logt = ifelse(time>0,log(time),NA)) %>% summarise(meanT = exp(mean(range(logt,na.rm=T)))), # midrange
#               by = "S")
#   
#   
#   ess = sim_df %>%  group_by(S) %>% filter(time == max(time)) %>% ungroup
#   perc = rep(NA,max(ess$S))
#   for (i in 1:max(ess$S)) {
#     tmp = filter(mc_df, S==i & id == 1) %>% pull(sum_n_stars)> ess$sum_n_stars[i] 
#     perc[i] = (1-round(mean(tmp,na.rm = T),2))
#   }
#   
#   
#   df3 = data.frame(x = dfN2 %>% with(unique(meanT)), l = perc)
#   
#   # browser()
#   out= 
#     ggplot(filter(dfN2, row_number()%%1 == 0), aes(x = log10(meanT), y = sum_n_stars))+
#     # geom_segment(data = dfN, aes(x = log10(minT), xend = log10(maxT), y = maxN, yend = maxN), color = "grey80")+
#     geom_vline(xintercept = log10(pull(filter(sim_df, branching==1),time)), linetype = "solid", color = "grey80")+
#     geom_quasirandom(aes(group = S), alpha = .5, color = "grey50", shape=16)+
#     geom_line(data=sim_df,aes(x = log10(time),y=sum_n_stars), size = 1, color = "red")+
#     geom_text(data=dfN2,aes(x=log10(meanT), label = paste("S = ",S)), y = exp(4.3))+
#     geom_text(data=df3,aes(x=log10(x), label = l), y = exp(4.2))+
#     # geom_label(data=df3,aes(x=log10(x), label = l), y = 15.5, label.size = 0)+
#     labs(x="Evolutionary time", y = expression(paste("Total abundance")))+
#     # labs(x="Evolutionary time", y = expression(paste("Total abundance (", sum(N[i]),")")))+
#     theme_bw() +
#   scale_y_continuous(trans = slog, limits = c(0,4.3), breaks = 0:5)+
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title=element_text(size=15),
#         axis.text =element_text(size=10))
#   return(out)
# }


## Figure 4 for relative yield


# figure4c = function(sim_df,mc_df,thin=1){
#   dfN2 <-
#     mc_df %>% 
#     filter(id == 2) %>% 
#     left_join(sim_df %>% group_by(S) %>% mutate(logt = ifelse(time>0,log(time),NA)) %>% summarise(meanT = exp(mean(range(logt,na.rm=T)))), # midrange
#               by = "S")
#   ess = sim_df %>%  group_by(S) %>% filter(time == max(time)) %>% ungroup
#   perc = rep(NA,max(ess$S))
#   for (i in 1:max(ess$S)) {
#     tmp = filter(mc_df, S==i & id == 2) %>% pull(sum_Y)> ess$sum_Y[i] 
#     perc[i] = (1-round(mean(tmp,na.rm = T),2))
#   }
#   df3 = data.frame(x = dfN2 %>% with(unique(meanT)), l = perc)
#   out= 
#     ggplot(filter(dfN2, row_number() %% thin == 0), aes(x = log10(meanT), y = sum_Y))+
#     geom_vline(xintercept = log10(pull(filter(sim_df, branching==1),time)), linetype = "solid", color = "grey80")+
#     geom_quasirandom(aes(group = S), alpha = .5, color = "grey80", shape=16)+
#     geom_line(data=sim_df,aes(x = log10(time),y=sum_Y), size = 1, color = "red")+
#     geom_text(data=dfN2,aes(x=log10(meanT), label = paste("S = ",S)), y = (3.3))+
#     geom_text(data=df3,aes(x=log10(x), label = l), y = (3.2))+
#     labs(x="Evolutionary time", y = expression(paste("Total relative yield")))+
#     theme_bw() +
#     theme(panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.title=element_text(size=15),
#           axis.text =element_text(size=10))
#   return(out)
# }
# 
# 
label_s = function(x,lettr){
  out = rep(NA,length(x))
 for(i in 1:length(x)){
   out[i] = if(i==1){paste0(lettr,"=",x[i])}else{paste(x[i])}
 }
  return(out)
}
# 
# # label_s2 = as_labeller(label_s)
# 
# figure4d = function(sim_df,mc_df,thin=1){
#   
#   dfN2 <-
#     mc_df %>% 
#     filter(id == 1) %>% 
#     left_join(sim_df %>% group_by(S) %>% mutate(logt = ifelse(time>0,log(time),NA)) %>% summarise(meanT = exp(mean(range(logt,na.rm=T)))), # midrange
#               by = "S")
#   
#   
#   ess = sim_df %>%  group_by(S) %>% filter(time == max(time)) %>% ungroup
#   perc = rep(NA,max(ess$S))
#   for (i in 1:max(ess$S)) {
#     tmp = filter(mc_df, S==i & id == 1) %>% pull(sum_Y)> ess$sum_Y[i] 
#     perc[i] = (1-round(mean(tmp,na.rm = T),2))
#   }
#   df3 = data.frame(x = dfN2 %>% with(unique(meanT)), l = perc)
#   out= 
#     ggplot(filter(dfN2, row_number()%%thin == 0), aes(x = log10(meanT), y = sum_Y))+
#     geom_vline(xintercept = log10(pull(filter(sim_df, branching==1),time)), linetype = "solid", color = "grey80")+
#     geom_quasirandom(aes(group = S), alpha = .5, color = "grey50", shape=16)+
#     geom_line(data=sim_df,aes(x = log10(time),y=sum_Y), size = 1, color = "red")+
#     geom_text(data=dfN2,aes(x=log10(meanT), label = paste("S = ",S)), y = (2.1))+
#     # geom_text(data=dfN2,aes(x=log10(meanT), label = label_s(S)), y = (2.1))+
#     geom_text(data=df3,aes(x=log10(x), label = l), y = (2))+
#     labs(x="Evolutionary time", y = expression(paste("Total relative yield")))+
#     theme_bw() +
#     theme(panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.title=element_text(size=15),
#           axis.text =element_text(size=10))
#   return(out)
# }



figure4g = function(sim_df,mc_df, ID, col_name, thin=1, labs_y = NULL, lims_y = NULL, textsize = 3.88){
  
  stopifnot(all(ID %in% c(1,2)))
  stopifnot(col_name %in% colnames(mc_df) )

  dfN2 <-
    mc_df %>% 
    filter(id == ID) %>% 
    left_join(sim_df %>% group_by(S) %>% mutate(logt = ifelse(time>0,log(time),NA)) %>% summarise(meanT = exp(mean(range(logt,na.rm=T)))), # midrange
              by = "S")
  ess = sim_df %>%  mutate(branch = lead(S) - S) %>% filter(branch == 1 | is.na(branch))

  perc = matrix(NA,nrow=max(sim_df$S),ncol=length(ID))
  for(randid in 1:length(ID)){
    for (i in min(ess$S):max(ess$S)) {
      tmp = (filter(mc_df, S==i & id == randid) %>% pull(col_name) ) > filter(ess,S==i)[[col_name]] 
      out = (1-round(mean(tmp,na.rm = T),2))
      # print(out)
      perc[i,randid] = out
    }
  }
  stopifnot(nrow(ess) == max(sim_df$S)-min(sim_df$S)+1)

  y_lab = switch(col_name,
    "sum_n_stars" = "Total abundance",
    "sum_Y" = "Total relative yield",
    "log_eta" = expression(paste("log(", eta,")")),
    "eta" =    expression(paste("Robustness (",eta,")")),
    col_name
  )
  
  dat = filter(mutate(dfN2, time = log10(time), meanT = log10(meanT)), row_number()%%thin == 0)
  range_perc = function(x,add){range(x[!is.infinite(x)],na.rm=T) %>% {diff(.)*add + .[1]}}
  
  
  if(is.null(labs_y) | !(length(labs_y) == 3)){
    labs_y = c(
      range_perc(mc_df[[col_name]], 1.3),
      range_perc(mc_df[[col_name]], 1.2),
      range_perc(mc_df[[col_name]], 1.1)
    )
  }
  
  if(is.null(lims_y)){
    lims_y = c(NA, labs_y[1])
  }
  
  
  out= 
    ggplot(dat, aes_(x = "meanT", y = col_name))+
    geom_vline(xintercept = (pull(filter(sim_df %>% mutate(time=log10(time)), branching==1),time)), linetype = "solid", color = "grey80")+
    # geom_quasirandom(aes_string(x = "meanT", y = col_name,group = "S"), alpha = .5, color = if(ID == 1){"grey50"}else{"grey80"}, shape=16)+
    {if(2 %in% ID)geom_violinhalf(data=filter(dat,id == 2), aes_string(x = "meanT", y = col_name,group = "S"),  flip = FALSE, color= "grey50", size = .3, fill = "grey80",scale="width")}+
    {if(1 %in% ID)geom_violinhalf(data=filter(dat,id == 1), aes_string(x = "meanT", y = col_name,group = "S"), flip = TRUE, color="grey50", size = .3, fill = "grey50",scale = "width")}+
    geom_line(data=sim_df %>% mutate(time=log10(time)),aes_string(x = "time",y=col_name), size = 1, color = "red")+
    geom_text(data=distinct(dfN2,S,.keep_all = T) %>% mutate(y=labs_y[1]), 
              aes(x=log10(meanT), label = label_s(S,"S"),y=y), size = textsize)+
    {if(1 %in% ID)geom_text(data = data.frame(x= log10(unique(dfN2$meanT)),y=labs_y[2] ,text=perc[,1] %>% .[!is.na(.)]),
                            aes(x = x, y = y, label = label_s(format(round(text, 2), nsmall = 2),"P1")),
                            size = textsize)}+
    {if(2 %in% ID)geom_text(data = data.frame(x= log10(unique(dfN2$meanT)),y=labs_y[3] ,text=perc[,2] %>% .[!is.na(.)]),
                            aes(x = x, y = y, label = label_s(format(round(text, 2), nsmall = 2),"P2")),
                            size = textsize)}+
    labs(x="Evolutionary time", y = y_lab)+
    scale_y_continuous(limits = lims_y, breaks = c(-1000:1000))+
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title=element_text(size=15),
          axis.text =element_text(size=10))
  return(out)
}
# 


# data.frame(x = runif(100),y=runif(100,1,2)) %>% 
#   ggplot(aes(x=x,y=y))+
#   geom_point()+
#   scale_y_continuous(trans = "log", limits = function(x){c(min(x),range_perc(x,1.5))} )
  

figure4g_log = function(sim_df,mc_df, ID, col_name, thin=1, labs_y = NULL, lims_y = NULL, textsize = 3.88){
  
  stopifnot(all(ID %in% c(1,2)))
  stopifnot(col_name %in% colnames(mc_df) )
  
  dfN2 <-
    mc_df %>% 
    filter(id == ID) %>% 
    left_join(sim_df %>% group_by(S) %>% mutate(logt = ifelse(time>0,log(time),NA)) %>% summarise(meanT = exp(mean(range(logt,na.rm=T)))), # midrange
              by = "S")
  ess = sim_df %>%  mutate(branch = lead(S) - S) %>% filter(branch == 1 | is.na(branch))
  
  perc = matrix(NA,nrow=max(sim_df$S),ncol=length(ID))
  for(randid in 1:length(ID)){
    for (i in min(ess$S):max(ess$S)) {
      tmp = (filter(mc_df, S==i & id == randid) %>% pull(col_name) ) > filter(ess,S==i)[[col_name]] 
      out = (1-round(mean(tmp,na.rm = T),2))
      # print(out)
      perc[i,randid] = out
    }
  }
  stopifnot(nrow(ess) == max(sim_df$S)-min(sim_df$S)+1)
  
  y_lab = switch(col_name,
                 "sum_n_stars" = "Total abundace",
                 "sum_Y" = "Total relative yield",
                 "log_eta" = expression(paste("log(", eta,")")),
                 "eta" = expression(paste("Robustness (",eta,")")),
                 col_name
  )
  

  dat = filter(mutate(dfN2, time = log10(time), meanT = log10(meanT)), row_number()%%thin == 0)
  range_perc = function(x,add){range(x[!is.infinite(x)],na.rm=T) %>% {diff(.)*add + .[1]}}
  
  if(is.null(labs_y) | !(length(labs_y) == 3)){
    labs_y = c(
      range_perc(mc_df[[col_name]], 10),
      range_perc(mc_df[[col_name]], 5),
      range_perc(mc_df[[col_name]], 3)
    )}
  
  if(is.null(lims_y)){
    lims_y = c(NA, labs_y[1])
  }
  
  
  out= 
    ggplot(dat, aes_(x = "meanT", y = col_name))+
    geom_vline(xintercept = (pull(filter(sim_df %>% mutate(time=log10(time)), branching==1),time)), linetype = "solid", color = "grey80")+
    {if(2 %in% ID)geom_violinhalf(data=filter(dat,id == 2), aes_string(x = "meanT", y = col_name,group = "S"), size = .3, flip = F, color= "grey50", fill = "grey80",scale="width")}+
    {if(1 %in% ID)geom_violinhalf(data=filter(dat,id == 1), aes_string(x = "meanT", y = col_name,group = "S"), size = .3, flip = T, color="grey50", fill = "grey50",scale = "width")}+
    geom_line(data=sim_df %>% mutate(time=log10(time)),aes_string(x = "time",y=col_name, group ="S"), size = 1, color = "red")+
    geom_text(data=distinct(dfN2,S,.keep_all = T) %>% mutate(y= labs_y[1]), aes(x=log10(meanT), label = label_s(S,"S"),y=y), size = textsize)+
    {if(1 %in% ID)geom_text(data = data.frame(x= log10(unique(dfN2$meanT)),y=labs_y[2] ,text=perc[,1] %>% .[!is.na(.)]),
                            aes(x = x, y = y, label = label_s(format(round(text, 2), nsmall = 2),"P1")),
                            size = textsize)}+
    {if(2 %in% ID)geom_text(data = data.frame(x= log10(unique(dfN2$meanT)),y=labs_y[3] ,text=perc[,2] %>% .[!is.na(.)]),
                            aes(x = x, y = y, label = label_s(format(round(text, 2), nsmall = 2),"P2")),
                            size = textsize)}+
    labs(x="Evolutionary time", y = y_lab)+
    scale_y_continuous(trans = "log", limits = lims_y , breaks = (10^(-10:10)))+
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title=element_text(size=15),
          axis.text =element_text(size=10))
  return(out)
}
