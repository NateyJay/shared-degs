
get_shared_degs <- function(lfc.df, lfc_threshold=1, fdr_threshold=0.95, study_threshold=2) {
  ## This function takes a dataframe as an input where columns=experiments, rows=genes, and the values are lfc's (also allows NA or 0 for non-DEGs). 
  
  out.ls <- list()
  
  out.ls$lfc_threshold = lfc_threshold
  out.ls$fdr_threshold = fdr_threshold
  out.ls$study_threshold = study_threshold
  
  for (direction in c('up','down')) {
    
    message(direction)
    
    out.ls[[direction]] = list()
    
    ## formatting the dataframe to be 0's for not-DE and 1 for DE in a given direction
    
    de.df <- lfc.df
    de.df[is.na(de.df)] <- 0
    
    if (direction == 'up') {
      sign = 1
      de.df[de.df < lfc_threshold] <- 0
      de.df[de.df >= lfc_threshold] <- 1
      
    } else {
      sign = -1
      de.df[de.df > lfc_threshold * sign] <- 0
      de.df[de.df <= lfc_threshold * sign] <- 1
    }
    
    
    deg_count = colSums(de.df)
    
    
    get_sde.df <- function(de.df) {
      ## this is meant to aggregate the de.df data here by study. This expects columnames to be named like this: "Study.repgroup"
      s.df <- de.df
      s.df$gene <- rownames(de.df)
      s.df <- melt(s.df, id.vars = c("gene"))
      names(s.df) <- c('gene', 'experiment', 'de')
      
      s.df$study <- str_split_fixed(s.df$experiment, '\\.',2)[,1]  
      s.df$experiment <- NULL
      
      s.df <- dcast(s.df, gene ~ study, sum, value.var = 'de')
      
      row.names(s.df) <- s.df$gene
      s.df$gene <- NULL
      s.df <- s.df[rownames(de.df),] 
      
      s.df[s.df > 0] <- 1
      
      return(s.df)
      
    }
    
    
    sde.df <- get_sde.df(de.df)
    
    
    
    
    
    get_random_genes <- function(x, background, times=1) {
      # column_names = c("rep", 1:length(x))
      out = data.frame(row.names=1:times)
      
      
      for (time in 1:times) {
        
        ## picking random genes (by a number) for all of the gene sets
        y = c()
        for (i in 1:length(x)) {
          number_of_degs = x[i]
          s = sample(x = 1:length(background), size = number_of_degs,replace = F)
          y <- c(y, s)
        }
        
        
        
        tab = tabulate(table(y), nbins= length(x))
        # tab = c(background - sum(tab), tab)
        out = rbind(out, tab)
        
      }
      names(out) <- 1:length(x)
      
      out = as.data.frame(t(out))
      
      names(out) <- str_glue("rand{1:times}")
      
      
      out <- rbind(`0`=length(background) - colSums(out), out)
      
      return(out)
      
    }
    
    random.df <- get_random_genes(x=deg_count, background=rownames(lfc.df), times=10)
    
    random.df$mean <- rowMeans(random.df)
    random.df$cumulative_mean <- cumsum(random.df$mean)
    random.df$total_proportion <- random.df$cumulative_mean / nrow(lfc.df)
    
    
    shared_count = rowSums(de.df)
    
    filtered_count = shared_count
    filtered_count[rowSums(sde.df) < study_threshold] <- NA
    
    
    comb.df <- random.df[, c('mean','cumulative_mean','total_proportion')]
    
    names(comb.df) <- c('rand_mean', 'rand_cmean', 'rand_prop')
    
    
    tab = table(shared_count)
    comb.df$shared <- tab[row.names(comb.df)]
    comb.df$shared[is.na(comb.df$shared)] <- 0
    
    comb.df$include <- comb.df$rand_prop > fdr_threshold
    
    tab = table(filtered_count)
    comb.df$filtered <- tab[row.names(comb.df)]
    comb.df$filtered[is.na(comb.df$filtered)] <- 0
    comb.df$filtered[!comb.df$include] <- 0
    
    
    
    out.ls[[direction]][['threshold']] = which(comb.df$include)[1]
    out.ls[[direction]][['de.df']] <- de.df
    out.ls[[direction]][['shared_counts']] <- shared_count
    out.ls[[direction]][['sde.df']] <- sde.df
    out.ls[[direction]][['study_counts']] <- rowSums(sde.df)
    
  }
  
  
  summary.df <- data.frame(row.names = rownames(lfc.df))
  
  summary.df$up           <- out.ls$up$shared_counts[rownames(summary.df)]
  summary.df$up_studies   <- out.ls$up$study_counts[rownames(summary.df)]
  summary.df$down         <- out.ls$down$shared_counts[rownames(summary.df)]
  summary.df$down_studies <- out.ls$down$study_counts[rownames(summary.df)]
  
  
  summary.df <- summary.df[order(-summary.df$up, summary.df$down),]
  
  f_up   <- summary.df$up >= out.ls$up$threshold & summary.df$up_studies >= study_threshold
  f_down <- summary.df$down >= out.ls$down$threshold & summary.df$down_studies >= study_threshold
  
  summary.df$de_call <- '-'
  summary.df$de_call[f_up] <- 'up'
  summary.df$de_call[f_down] <- 'down'
  summary.df$de_call[f_up * f_down] <- 'up/down'
  
  out.ls$summary.df <- summary.df
  
  return(out.ls)
  
}