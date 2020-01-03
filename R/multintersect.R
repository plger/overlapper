#' multintersect
#'
#' Performs pair-wise overlaps between the elements of one or two lists
#'
#' @param ll A list of vectors to be compared
#' @param ll2 An optional second list of vectors
#' @param universe An optional vector of the universe (all terms), used to calculate overlap probabilities and enrichments. If NULL, the union of all lists will be used. Only elements of `ll` and `ll2` that are in the `universe` will be considered.
#' @param addSetSize Logical; whether to add the set sizes to set names (default TRUE)
#' @param breakNames If not NULL (default), should be an number indicating the character length threshold above which set names should be split onto two lines.
#'
#' @return Returns a list
#'
#' @export
multintersect <- function(ll, ll2=NULL, universe=NULL, addSetSize=TRUE, breakNames=NULL){
  if(!is.list(ll) || is.null(names(ll))) stop("`ll` should be a named list.")
  if(is.null(ll2)){
    symm <- TRUE
    ll2 <- ll
  }else{
    if(!is.list(ll2) || is.null(names(ll2))) 
      stop("If given, `ll2` should be a named list.")
    symm <- FALSE
  }
  if(is.null(universe)) universe <- unique(c(unlist(ll),unlist(ll2)))
  ll <- lapply(ll,y=universe,FUN=intersect)
  ll2 <- lapply(ll2,y=universe,FUN=intersect)
  ll <- ll[sapply(ll,FUN=length)>0]
  ll2 <- ll2[sapply(ll2,FUN=length)>0]
  if(addSetSize){
    names(ll) <- paste0(names(ll)," (",sapply(ll,FUN=length),")")
    names(ll2) <- paste0(names(ll2)," (",sapply(ll2,FUN=length),")")
  }
  if(!is.null(breakNames)){
    names(ll) <- breakStrings(names(ll),breakNames)
    names(ll2) <- breakStrings(names(ll2),breakNames)
  }
  n <- length(ll)
  j <- length(ll2)
  m <- matrix(0,nrow=n,ncol=j)
  colnames(m) <- names(ll2)
  rownames(m) <- names(ll)
  prob <- m
  enr <- m
  jacc <- m
  of <- m
  for(i in 1:n){
    for(j in 1:length(ll2)){
      of[i,j] <- paste(intersect(ll[[i]],ll2[[j]]),collapse=", ")
      m[i,j] <- length(intersect(ll[[i]],ll2[[j]]))
      if(names(ll)[i]==names(ll2)[j]){
        prob[i,j] <- NA
        enr[i,j] <- NA
        jacc[i,j] <- NA
      }else{
        enr[i,j] <- getEnrichment(ll[[i]],ll2[[j]],universe)
        prob[i,j] <- overlap.prob(ll[[i]],ll2[[j]],universe,lower=enr[i,j]<1)
        jacc[i,j] <- m[i,j]/length(unique(c(ll[[i]],ll2[[j]])))
      }
    }
  }
  list(ll1=ll,ll2=ll2,enr=enr,prob=prob,jacc=jacc,m=m,of=of,symm=symm)
}


#' plot.multintersect
#'
#' Plots the results of `multintersect` using an interactive heatmap
#'
#' @param res A list as produced by `multintersect`
#' @param keyCol The values to use for the heatmap's colors. Possible values are: prob, enrichment, log2Enrichment, log10Prob, overlap, jaccard
#' @param keyWrite The values to be written in the cells. Possible values are: prob, enrichment, log2Enrichment, log10Prob, overlap, jaccard
#' @param margin Either a single number indicating the size of the left and bottom margins, or a list (l,r,b,t,pad) indicating the size of each margin.
#' @param title Plot title, default ''
#' @param cluster Logical; whether to cluster rows and columns (default TRUE)
#'
#' @return Plots a plotly heatmap
#'
#' @import plotly
#' @export
plot.multintersect <- function(res, keyCol="log2Enrichment", keyWrite="overlap", margin=100, title="", cluster=TRUE){
  keyCol <- match.arg(keyCol, c("prob","enrichment","log2Enrichment","log10Prob","overlap","jaccard"))
  keyWrite <- match.arg(keyWrite, c("prob","enrichment","log2Enrichment","log10Prob","overlap","jaccard"))
  x <- switch(keyCol,
              prob=res$prob,
              enrichment=res$enr,
              log2Enrichment=log2(res$enr),
              log10Prob=-log10(res$prob),
              overlap=res$m,
              jaccard=res$jacc)
  y <- switch(keyWrite,
              prob=format(res$prob,digits=1),
              enrichment=format(res$enr,digits=2,trim=T,drop0trailing=T),
              log2Enrichment=round(log2(res$enr),1),
              log10Prob=round(-log10(res$prob),1),
              overlap=res$m,
              jaccard=round(res$jacc,3))
  y[is.infinite(y)] <- NA
  x[is.infinite(x) & x>0] <- max(x[which(!is.infinite(x))],na.rm=T)
  x[is.infinite(x) & x<0] <- min(x[which(!is.infinite(x))],na.rm=T)
  lab <- matrix(paste0( rep( gsub("\n"," ",row.names(x),fixed=T),ncol(x) ), "\n",
                        rep( gsub("\n"," ",colnames(x), fixed=T),each=nrow(x)), "\n",
                        "overlap: ", as.numeric(res$m), "\n",
                        as.character(format(res$enr,digits=2,trim=T,drop0trailing=T)),"-fold enrichment \n",
                        "p~",as.character(format(res$prob,digits=2))
  ),nrow=nrow(x),ncol=ncol(x))

  if(cluster){
    ro <- hclust(dist(x))$order
    co <- hclust(dist(t(x)))$order
    x <- x[ro,co]
    y <- y[ro,co]
    lab <- lab[ro,co]
  }

  if(res$symm){
    for(i in 1:ncol(x)) lab[i,i] <- row.names(x)[i]
  }
  if(length(margin)==1){
    margin <- list(l=margin, r=5, b=margin, t=ifelse(title=="",4,30), pad=4)
  }
  p <- plot_ly(x=colnames(x), y=row.names(x), z=x, type="heatmap", text=lab, hoverinfo = 'text') %>% colorbar(title = keyCol) %>%
    layout(title = title, xaxis = list(showgrid = FALSE), yaxis=list(showgrid = FALSE), margin = margin)
  p %>% add_annotations(x = rep(colnames(x),each=nrow(x)), y = rep(row.names(x),ncol(x)), text = y, xref="x", yref="y", showarrow=FALSE)
}


#' dotplot.multintersect
#'
#' Plots the results of `multintersect` using a dot plot
#'
#' @param m A list as produced by `multintersect`
#' @param sizeRange Size range of the dots
#' @param colors A named vector of 'low', 'mid', and 'high' colors
#' @param forceDivergent Logical; whether to force the use of a divergent color
#' palette; otherwise only low/high colors are used when there is no negative 
#' enrichment.
#'
#' @return A ggplot.
#' @import ggplot
#' @importFrom reshape2 melt
#' @export
dotplot.multintersect <- function(m, sizeRange=c(0,20), 
                                  colors=c(low="blue",mid="grey",high="yellow"),
                                  forceDivergent=FALSE){
  ml <- list( fill=log2(m$enr), val=m$m, size=-log10(m$prob) )
  ml <- lapply(ml, FUN=function(x){
    if(nrow(x)==ncol(x) && all(row.names(x)==colnames(x))){
      x[lower.tri(x,diag=TRUE)] <- NA
      x <- x[-nrow(x),-1]
    }
    x <- t(x)
    x[nrow(x):1,]
  })
  for(f in names(ml)) ml[[f]] <- melt(ml[[f]], value.name=f)
  d <- cbind(ml[[1]], do.call(cbind, lapply(ml[-1], FUN=function(x) x[,3,drop=FALSE])))
  
  p <- ggplot(d, aes(Var1, Var2, label=val)) + 
    geom_point(aes(size=size,colour=fill)) + geom_text() + 
    scale_size_continuous(range=sizeRange) + 
    theme(axis.line = element_blank(), axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(colour="log2(Enrichment)", size="-log10(p-value)")
  if(any(d$fill<0,na.rm = TRUE) || forceDivergent){
    p <- p + scale_color_gradient2(low=colors["low"], mid=colors["mid"], high=colors["high"], na.value = "white")
  }else{
    p <- p + scale_color_gradient(low=colors["low"], high=colors["high"], na.value = "white")
  }
  p
}


#' breakStrings
#'
#' breaks a string of words (or vector thereof) into two lines
#'
#' @param x a character vector
#' @param minSizeForBreak the minimum number of characters to break on two lines (default 20)
#' @param lb the line break character (default "\n")
#'
#' @return a character vector of length=length(x)
#'
#' @export
breakStrings <- function(x, minSizeForBreak=20, lb="\n"){
  sapply(x,minSizeForBreak=minSizeForBreak,lb=lb,FUN=function(x,minSizeForBreak,lb){
    if(nchar(x)<=minSizeForBreak)	return(x)
    g <- gregexpr(" ", x)[[1]]
    if(length(g)==0) return(x)
    if(length(g)==1 & all(g==-1)) return(x)
    mid <- nchar(x)/2
    mid <- g[order(abs(g-mid))[1]]
    substr(x, mid, mid) <- lb
    return(x)
  })
}

#' getEnrichment
#'
#' Returns the enrichment ratio of set2 in set1
#'
#' @param set1 the set in which to test for enrichment
#' @param set2 the set for which to test enrichment
#' @param universe either a character vector containing the universe/background (recommended), or an integer indicating the size of the universe
#'
#' @return A numeric value indicating the enrichment ratio (observed/expected)
#'
#' @export
getEnrichment <- function(set1,set2,universe){
  set1 <- unique(as.character(set1))
  set2 <- unique(as.character(set2))
  if(!is.numeric(universe) | length(universe)>1){
    universe <- unique(universe)
    set1 <- intersect(set1,universe)
    set2 <- intersect(set2,universe)
    universe <- length(universe)
  }
  ov <- sum(set1 %in% set2)
  expected <- length(set1)*length(set2)/universe
  return(ov/expected)
}


#' overlap.prob
#'
#' Calculates the probability of the observed overlap between two sets on the basis of the hypergeometric distribution.
#'
#' @param set1 a character vector
#' @param set2 a character vector
#' @param universe either a character vector containing the universe/background (recommended), or an integer indicating the size of the universe
#' @param lower logical; whether the lower tail should be considered (default FALSE). By default the function tests for an enrichment, set to TRUE to test for a depletion.
#'
#' @return the probability of the overlap
#'
#' @examples
#' overlap.prob( c("A","B","C","D"), c("B","D","G"), LETTERS )
#'
#' @export
overlap.prob <- function(set1,set2,universe,lower=F){
  set1 <- as.character(set1)
  set2 <- as.character(set2)
  if(class(universe)=="character"){
    set1 <- intersect(set1,universe)
    set2 <- intersect(set2,universe)
    universe <- length(unique(universe))
  }
  set1 <- unique(set1)
  set2 <- unique(set2)
  ov <- sum(set1 %in% set2)
  phyper(max(0,ov-1), length(set1), universe-length(set1), length(set2), lower.tail=lower)
}


#' multintersect.test
#'
#' Tests the significance of an overlap across multiple sets using 
#' permutations.
#'
#' @param ll A list of sets.
#' @param u The universe.
#' @param nperm The number of permutations to use.
#' @param two.tailed Whether the test should be two-tailed (defaults to testing
#' a greater alternative)
#'
#' @return A p-value, with precision roughly 1/nperm
#' @export
#'
#' @examples
#' # create random sets
#' ll <- lapply(1:3,FUN=function(x) sample(letters, 8))
#' multintersect.test(ll, u=letters)
multintersect.test <- function(ll, u, nperm=2000, two.tailed=FALSE){
  ll <- lapply(ll, unique)
  ll <- lapply(ll, intersect, y=u)
  ll <- ll[sapply(ll,length)>0]
  if(length(ll)<2) return(NULL)
  o <- u
  for(i in seq_along(ll)) o <- intersect(o,ll[[i]])
  if(length(o)==0 && !two.tailed) return(1)
  sizes <- sapply(ll,length)
  perms <- sapply(seq_len(nperm), FUN=function(x){
    ll2 <- lapply(sizes, FUN=function(s) sample(u,s))
    o <- u
    for(i in seq_along(ll2)) o <- intersect(o,ll2[[i]])
    length(o)
  })
  if(two.tailed) return(sum(perms >= length(o) | perms <= length(o))/nperm)
  sum(perms >= length(o))/nperm
}
