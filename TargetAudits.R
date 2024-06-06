rm(list=ls())

library(tidyverse)
library(data.table)
library(lubridate)
library(readxl)
library(sf)
library(ggplot2)
library(nloptr)

setwd("C:/Users/lbeatty/Documents/TargetAudits")

source('funs.R')

xmin = -125 # area of interest for now encapsulate CONUS
xmax = -75
ymin = 26
ymax = 50

#general outline
# read in data and set parameters of problem (grid size, number of grid/n bins, inspection budget)
# grid up emissions data and calculate n wells and agg emissions by grid
# calculate detected emissions under uniform policy (probability of inspection is budget/n_wells)
# calculate optimal inspection probabilities by emission/n_wells bin
# calculate detected emissions at optimum
# make plots and write output

for(gridloop in c(0.2, 1)){
  for(binloop in c(5,10)){
    for(budgetloop in c(5000, 10000,50000)){
      
      #read file
      emis = read_excel('17 wells_emissions.xlsx')
      emis = st_as_sf(emis, coords=c('Longitude', 'Latitude'), crs='WGS84')
      
      ##############
      #define stuff #
      ###############
      gridsize = gridloop #for now just make grid thats lat/lons rather than dist
      
      budget = budgetloop #size of audit budget
      
      n_bins = binloop #size of grid emission and grid count bins
      #I'm going to have the regulator observe gridded emissions and gridded well counts
      #then assign an audit probability to each combination of of emission bin and count bin that
      #maximizes detected emissions
      
      #############
      #plot stuff #
      #############
      # ggplot(emis)+
      #   geom_sf(size=0.1)+
      #   theme_bw()+
      #   coord_sf(xlim=c(xmin,xmax), ylim=c(ymin,ymax))
      # ggsave(filename='source_map.jpg',
      #        height=6,
      #        width=8,
      #        device='jpeg')
      
      ################
      ## make grid ##
      ###############
      
      gridx = seq(xmin, xmax, by=gridsize)
      gridy = seq(ymin, ymax, by=gridsize)
      grid = expand.grid(gridx, gridy)
      grid = st_as_sf(grid, coords = c("Var1", "Var2"))
      grid$geometry = st_buffer(grid$geometry, dist=gridsize/2, endCapStyle="SQUARE")
      st_crs(grid)="WGS84"
      grid$ID = 1:nrow(grid)
      
      #calculate grid-level emissions
      emis_grid = st_join(emis, grid)
      emis_grid = data.table(emis_grid)
      
      emis_grid = emis_grid[,.(count=.N, grid_emis=sum(ER)), by='ID']
      emis_grid = emis_grid[!is.na(ID),]
      
      ##############
      ## Policies
      #############
      emis = st_join(emis, grid)
      emis=data.table(emis)
      emis = emis[!is.na(ID),] #drop if not within xlim, ylim
      
      #################################
      ##### uniform probability #######
      #################################
      emis_uniform = copy(emis)
      emis_uniform[,prob:=budget/nrow(emis_uniform)]
      
      #just going to calculate expectation of discovered emissions rather than simulate audit outcomes
      emis_uniform[,detected:=ER*prob]
      
      
      #####################################
      ##### target based on grid emis #####
      #####################################
      
      #make bins
      count_bins = quantile(emis_grid$count, seq(0,1,by=1/n_bins))
      emis_bins = quantile(emis_grid$grid_emis, seq(0,1,by=1/n_bins))
      
      #throw out bins if bin1=bin2
      delete_count=c()
      delete_emis = c()
      for(i in 1:(n_bins-1)){
        if(count_bins[i]==count_bins[i+1]){delete_count=c(delete_count, i)}
        if(emis_bins[i]==emis_bins[i+1]){delete_emis = c(delete_emis, i)}
      }
      if(!is.null(delete_count)){count_bins = count_bins[-delete_count]}
      if(!is.null(delete_emis)){emis_bins = emis_bins[-delete_emis]}
      
      #calculate each bin
      emis_grid[,count_bin := cut(count, count_bins, include.lowest = T, labels=1:(length(count_bins)-1))]
      emis_grid[,emis_bin := cut(grid_emis, emis_bins, include.lowest=T, labels=1:(length(emis_bins)-1))]
      
      
      #################
      # there are some grid combos that don't show up
      prob_by_bin = expand.grid(1:(length(count_bins)-1), 1:(length(emis_bins)-1))
      names(prob_by_bin)=c("count_bin", "emis_bin")
      prob_by_bin = data.table(prob_by_bin)
      prob_by_bin[, names(prob_by_bin) := lapply(.SD, as.character)]
      empties = emis_grid[prob_by_bin, on=.(count_bin, emis_bin)]
      empties = empties[is.na(ID),]
      empties[,empty:=1]
      empties = empties[,c('count_bin', 'emis_bin', 'empty')]
      prob_by_bin = empties[prob_by_bin, on=.(count_bin, emis_bin)]
      prob_by_bin = prob_by_bin[is.na(empty),]
      prob_by_bin[,empty:=NULL]
      
      ## Optimize based on grid emissions
      start=rep(0, nrow(prob_by_bin))
      sol = nloptr(start,
                   eval_f=grid_emis_objective,
                   eval_g_ineq = grid_emis_constraint,
                   lb = rep(0, nrow(prob_by_bin)),  #constrain between 0,1
                   ub = rep(1, nrow(prob_by_bin)),
                   opts = list("algorithm" = "NLOPT_LN_COBYLA",
                               "ftol_rel"=1e-4,
                               "xtol_rel"=1e-4,
                               "maxeval" = 10000))
      
      #calculate discovered emis at optimum
      emis_grid_target=copy(emis)
      prob_by_bin[,prob:=sol$solution]
      grid_probs = emis_grid[prob_by_bin, on=.(count_bin, emis_bin)]
      emis_grid_target = emis_grid_target[grid_probs, on=.(ID)]
      emis_grid_target[,detected:=ER*prob]
      detected_targeted = sum(emis_grid_target$detected)
      audits_targeted = sum(emis_grid_target$prob)
      
      
      
      #make plots of outputs
      
      #optimal audit probs by bin
      ggplot(prob_by_bin%>%mutate(count_bin=as.numeric(count_bin), emis_bin=as.numeric(emis_bin))) +
        geom_tile(aes(x=emis_bin, y=count_bin, fill=prob)) +
        scale_fill_gradient(low = "green", high = "red") +
        ggtitle("optimal audit probabilities by emissions bin and count bin")+
        coord_fixed()+
        theme_bw()
      ggsave(filename = paste('gridsize', gridsize, '/audit_probabilities_nbins', n_bins, '_budget', budget, '.jpg', sep=''),
             width=8,
             height=6,
             device='jpeg')
      
      
      ####
      ## Maps
      emis_grid = left_join(emis_grid, grid, by="ID")
      emis_grid = st_sf(emis_grid)
      
      ggplot(emis_grid)+
        geom_sf(aes(fill=grid_emis), color=NA)+
        scale_fill_viridis_b()+
        theme_bw()
      ggsave(filename=paste('gridsize', gridsize, '/grid_emissions_map_CONUS.jpg', sep=''),
             width=8,
             height=6,
             device='jpeg')
      
      ggplot(emis_grid)+
        geom_sf(aes(fill=grid_emis), color=NA)+
        scale_fill_viridis_b()+
        coord_sf(xlim=c(-104,-100), ylim=c(30,34))+
        theme_bw()
      ggsave(filename=paste('gridsize', gridsize, '/grid_emissions_map_Permian.jpg', sep=''),
             width=8,
             height=6,
             device='jpeg')
      
      
      
      ggplot(emis_grid)+
        geom_sf(aes(fill=count), color=NA)+
        scale_fill_viridis_b()+
        theme_bw()
      ggsave(filename=paste('gridsize', gridsize, '/grid_count_map_CONUS.jpg', sep=''),
             width=8,
             height=6,
             device='jpeg')
      
      ggplot(emis_grid)+
        geom_sf(aes(fill=count), color=NA)+
        scale_fill_viridis_b()+
        coord_sf(xlim=c(-104,-100), ylim=c(30,34))+
        theme_bw()
      ggsave(filename=paste('gridsize', gridsize, '/grid_count_map_Permian.jpg', sep=''),
             width=8,
             height=6,
             device='jpeg')
      
      #plot solution audit probabilities
      emis_grid_target = left_join(emis_grid, prob_by_bin, by=c("count_bin", "emis_bin"))
      ggplot(emis_grid_target)+
        geom_sf(aes(fill=prob), color=NA)+
        scale_fill_viridis_b()+
        theme_bw()
      ggsave(filename=paste('gridsize', gridsize, '/bins', n_bins, '_budget', budget, '_grid_map_prob_CONUS.jpg', sep=''),
             width=8,
             height=6,
             device='jpeg')
      
      ggplot(emis_grid_target)+
        geom_sf(aes(fill=prob), color=NA)+
        scale_fill_viridis_b()+
        coord_sf(xlim=c(-104,-100), ylim=c(30,34))+
        theme_bw()
      ggsave(filename=paste('gridsize', gridsize, '/bins', n_bins, '_budget', budget, '_grid_map_prob_PERMIAN.jpg', sep=''),
             width=8,
             height=6,
             device='jpeg')
      
      write.table(t(c(gridsize, n_bins, budget, sum(emis_uniform$detected), detected_targeted, audits_targeted, sol$message)),
                  file = 'output.csv',
                  append=T, sep=",", col.names = F)
      
      write.csv(prob_by_bin,
                file=paste('gridsize', gridsize, '/nbins', n_bins, '_budget', budget, '_prob_solution.csv', sep=''))
    }
  }
}
