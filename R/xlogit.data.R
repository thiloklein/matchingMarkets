#' @title data.frame for exploded logit model
#'
#' @description shape a data.frame in a suitable form for the use of the exploded logit with the mlogit function in package mlogit.
#
#' @param data a data.frame
## @param x a mlogit.data or a pseries object
#' @param choice the variable indicating the choice made: it can be either a logical vector, a numerical vector with 0 where the alternative is not chosen, a factor with level 'yes' when the alternative is chosen
#' @param shape the shape of the data.frame: whether long if each row is an alternative or wide if each row is an observation
#' @param varying the indexes of the variables that are alternative specific
#' @param sep the seperator of the variable name and the alternative name (only relevant for a wide data.frame)
#' @param alt.var the name of the variable that contains the alternative index (for a long data.frame only) or the name under which the alternative index will be stored (the default name is alt)
#' @param chid.var the name of the variable that contains the choice index or the name under which the choice index will be stored
#' @param alt.levels the name of the alternatives: if null, for a wide data.frame, they are guessed from the variable names and the choice variable (both should be the same), for a long data.frame, they are guessed from the alt.var argument
#' @param id.var the name of the variable that contains the individual index if any
## @param group.var the name of the variable that contains the group index if any
#' @param opposite returns the opposite of the specified variables
#' @param drop.index should the index variables be dropped from the data.frame
#' @param ranked a logical value which is true if the response is a rank
## @param subset a logical expression which defines the subset of observations to be selected
#' @param ... further arguments passed to reshape
#'
#' @export
#'
#' @return \code{xlogit.data} returns a mlogit.data object, which is a data.frame in long format, i.e. one line for each alternative. It has a index attribute, which is a data.frame that contains the index of the choice made ('chid'), the index of the alternative ('alt') and, if any, the index of the individual ('id'). The choice variable is a boolean which indicates the choice made. This function use reshape if the data.frame is in wide format.
#' @author Thilo Klein
#' @keywords generate
#' @examples
#' \dontrun{
#' ## Game2 is a data.frame in long format for which the response is a
#' ## ranking variable (see example in package mlogit, function mlogit.data)
#' 
#' ## --- 1. baseline case (equivalent to mlogit package)
#' 
#' ## load data
#' data("Game2", package = "mlogit")
#' head(Game2, 13)
#' 
#' ## transform data
#' G <- xlogit.data(Game2, shape = "long", choice = "ch", 
#'                  alt.var = "platform", ranked = TRUE)
#' head(G, 30)
#'
#' ## --- 2. handle missing values 
#' 
#' ## introduce NAs
#' Game2$ch[Game2$ch != 1] <- NA
#'
#' ## transform data
#' G <- xlogit.data(Game2, shape = "long", choice = "ch", alt.var = "platform", ranked = TRUE)
#' head(G, 30)
#'
#' ## --- 3. handle varying choice sets 
#' 
#' ## drop choice options
#' Game2 <- Game2[-c(7:8),]
#' 
#' ## transform data
#' G <- xlogit.data(Game2, shape = "long", choice = "ch", alt.var = "platform", ranked = TRUE)
#' head(G, 30)
#' 
#' ## --- 4. run models
#' summary(mlogit::mlogit(ch ~ own | -1 + hours + age, G, reflevel = "PC"))
#' summary(mlogit::mlogit(ch ~ 0 | -1 + age | own, G))
#' summary(mlogit::mlogit(ch ~ 0 | -1 + age, G))
#' }

xlogit.data <- function(data, choice, shape = c("wide", "long"), varying = NULL, 
                           sep = ".", alt.var = NULL, chid.var = NULL, alt.levels = NULL, 
                           id.var = NULL, opposite = NULL, drop.index = FALSE, ranked = FALSE, 
                           ...) {
  if (is.null(chid.var)) {
    chid.name <- "chid"
    chid.is.variable <- FALSE
  }
  alt.name <- alt.var 
  alt.is.variable <- TRUE 
  if (!is.factor(data[[alt.name]])) 
    data[[alt.name]] <- factor(data[[alt.name]])
  alt.levels <- levels(data[[alt.name]])
  J <- length(alt.levels)
  alt <- data[[alt.name]]
  ts <- table(data$chid)
  ts <- as.numeric(ts)
  n <- length(ts)
  if (!chid.is.variable) 
    chid <- rep(1:n, ts)
  # chid <- as.factor(chid)
  # alt <- as.factor(alt)
  row.names(data) <- paste(chid, alt, sep = ".")
  chidpos <- which(names(data) == chid.name)
  altpos <- which(names(data) == alt.name)
  index <- data.frame(chid = chid, alt = alt)
  rownames(index) <- rownames(data)
  attr(data, "index") <- index
  attr(data, "class") <- c("mlogit.data", "data.frame")
  

  choicename = choice
  x<- data
  choicepos <- match(choicename, names(x))
  id <- attr(x, "index")$chid
  lev.id <- levels(id)
  theid <- as.numeric(id)
  oalt <-  attr(x, "index")$alt
  lev.alt <- levels(oalt)
  choice <- x[[choicename]]
  J <- length(unique(choice))
  d <- data.frame()
  chid <- c()
  alt <- c()
  id <- c()
  k <- 0
  for (i in unique(theid)){
    aid <- which(theid == i)
    adata <- x[aid, - choicepos]
    achoice <- choice[aid]
    JJ <- length(unique(achoice))
    aalt <- oalt[aid]
    J1 = length(achoice)
    remAlts <- rep(TRUE, J1)
    alogchoice <- achoice == 1
    d <- rbind(d, cbind(adata, alogchoice))
    Z <- sum(remAlts)
    k <- k + 1
    chid <- c(chid, rep(k, Z))
    id <- c(id, rep(i, Z))
    alt <- c(alt, aalt)
    if((JJ - 2)>0){
      for (j in 1:(JJ - 2)){
        k <- k + 1
        min.index <- achoice == j
        remAlts[min.index] <- FALSE
        Z <- sum(remAlts)
        chid <- c(chid, rep(k, Z))
        alt <- c(alt, aalt[remAlts])
        id <- c(id, rep(i, Z))
        alogchoice <- achoice[remAlts] == j + 1
        d <- rbind(d, cbind(adata[remAlts,], alogchoice))
      }
    }
  }
  colnames(d)[length(d)] <- choicename
  alt <- factor(alt, labels = lev.alt)
  index <- data.frame(chid = chid, alt = alt, id = id)
  rownames(d) <- rownames(index) <- paste(chid, as.character(alt), sep = ".")
  zzzz <- structure(d, index = index, class = c('mlogit.data', 'data.frame'))
  zzzz$ch[which(is.na(zzzz$ch))] = FALSE
  zzzz
}
