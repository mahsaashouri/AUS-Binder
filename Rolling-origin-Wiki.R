library(hts)
library(Matrix)
library(reshape2)
library(tidyverse) 
library(readr)
source('OLSmodel_se.R')
source('smatrix.R')

################
### Rolling origin
################
## reading data

unzip(zipfile = "wikipedia_data.zip", exdir = getwd())
wikipedia_data <- read_csv("wikipedia_data.csv")

## Data reshaping
# 394: length of each series and 913: number of series
wikipedia_wide1 <- wikipedia_data$views %>%
  matrix(nrow = 394, ncol = 913) %>%
  as.data.frame() %>%
  ts(frequency = 7)
colnames(wikipedia_wide1) <- unique(wikipedia_data$cat_column) %>% substr(1,14)

##################
## using hierarchies and groupings up to 2-way combinations
##################
wikigts <- gts(wikipedia_wide1, character=c(7,2,2,3),
               gnames = c("Access",
                          "Agent",
                          "Language",
                          "Purpose",
                          "Access x Agent",
                          "Access x Language",
                          "Access x Purpose",
                          "Agent x Language",
                          "Agent x Purpose",
                          "Language x Purpose"))

# Splitting data into training and test sets
wikitrain <- window(wikigts, end = c(1, 366))
wikitest <- window(wikigts, start = c(1, 367))

# Construct matrix of all time series including aggregates
ally <- aggts(wikigts)
ally.test <- aggts(wikitest)

n <- nrow(ally)
h <- 28
fc <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=5))
dimnames(fc) <- list(
  Horizon = paste0("h=",seq(h)),
  Series = colnames(ally),
  Method = c("OLS", "OLS.lwr", "OLS.upr", "OLS.se", "OLS.residual.scale")
)
## computation time
start_time <- Sys.time()
for(i in seq(NCOL(ally)))
{
  # OLS forecasts
  fit.OLS <- OLSmodel(ally[,i],7,7,h, nolag = c(1,7))
  fc[, i,"OLS"] <- fit.OLS[[1]]
  fc[, i,"OLS.lwr"] <- fit.OLS[[2]]
  fc[, i,"OLS.upr"] <- fit.OLS[[3]]
  fc[, i,"OLS.se"] <- fit.OLS[[4]]
  fc[, i,"OLS.residual.scale"] <- fit.OLS[[5]]
}
end_time <- Sys.time()
end_time - start_time

fc.OLS <- as.data.frame(fc[,,"OLS"])
fc.OLS.lwr <- as.data.frame(fc[,,"OLS.lwr"])
fc.OLS.upr <- as.data.frame(fc[,,"OLS.upr"])
fc.OLS.se <- as.data.frame(fc[,,"OLS.se"])
fc.OLS.residual.scale <- as.data.frame(fc[,,"OLS.residual.scale"])

colnames(fc.OLS) <- c(1:ncol(ally))
colnames(fc.OLS.lwr) <- c(1:ncol(ally))
colnames(fc.OLS.upr) <- c(1:ncol(ally))
colnames(fc.OLS.se) <- c(1:ncol(ally))
colnames(fc.OLS.residual.scale) <- c(1:ncol(ally))

as.matrix <- Matrix::as.matrix
t <- Matrix::t
solve <- Matrix::solve
diag <- Matrix::diag

## computing reconceliation matrix
gmat <- GmatrixG(wikigts$groups)
smatrix <- as((SmatrixM(gmat)), 'dgCMatrix')
lambda <- as(diag(rowSums(smatrix)), 'dgCMatrix')

rec.adj.lambda <- as.matrix(smatrix%*%solve(t(smatrix)%*%solve(lambda)%*%smatrix)%*%t(smatrix)%*%solve(lambda))

fc.rec <- matrix(NA, nrow = 28, ncol = ncol(ally))
for(i in 1:nrow(fc.OLS)){
  f.1 <- matrix(as.numeric(fc.OLS[i,]), ncol = 1, nrow = ncol(fc.OLS))
  fc.rec [i,] <- rec.adj.lambda %*% f.1
}
colnames(fc.rec ) <- c(1:ncol(ally))



## PI variance
fc.OLS.PI <- (fc.OLS.se)^2 + (fc.OLS.residual.scale)^2

## computing predictors (trend, dummy seasonality, lags) for each series
Xmat<-list()
freq <-7
maxlag <- 7
ally.final <- as.list(ally) 
Xmatrix<-function(X){
  X<-as.vector(X)
  intercept <- rep(1, length(X))
  trend1<-seq(NROW(X))
  trend2<-seq(NROW(X))^2
  season<-forecast::seasonaldummy(ts(X,frequency = freq))
  Xlag<-quantmod::Lag(X,k=1:maxlag)
  X_mat<-cbind.data.frame(intercept,trend1, trend2,season,Xlag)
  Xmat[[length(Xmat)+1]] <- X_mat 
}
Xmat.final <- lapply(ally.final, Xmatrix)
k <- 28
result.var <- matrix(NA, nrow = k, ncol = NCOL(ally))
## for loops for computing rolling base forecasts
for(i in 1:k){
  Xmat.final.train <- lapply(Xmat.final, function(x)x[1:((n - k) + (i - 1)),])
  Xmat.final.test <- lapply(Xmat.final, function(x)x[(n - k) + i,])
  mat <- as(Matrix::bdiag(lapply(Xmat.final.train, function(x){as.matrix(na.omit(x))})), 'dgCMatrix')
  mat.inverse <- as(Matrix::solve(t(mat)%*%mat), 'dgCMatrix')
  mat.test <- as(Matrix::bdiag(lapply(Xmat.final.test, function(x){as.matrix(na.omit(x))})), 'dgCMatrix')
  H.matrix <- as(mat.test %*% mat.inverse %*% t(mat.test), 'dgCMatrix')
  Sigma.mat <- as(Matrix::diag(fc.OLS.PI[i,]) + (Matrix::diag(fc.OLS.PI[i,]) %*% H.matrix), 'dgCMatrix')
  rec.p <- as(as.matrix(solve(t(smatrix)%*%solve(lambda)%*%smatrix)%*%t(smatrix)%*%solve(lambda)), 'dgCMatrix')
  var.for <- as.matrix((smatrix %*% rec.p) %*% Sigma.mat %*% (t(rec.p) %*% t(smatrix)))
  result.var[i,] <- as.vector(diag(var.for))
}

OLS.unrec = reshape2::melt(fc.OLS)
OLS.rec = reshape2::melt(fc.rec) 
OLS.var.rec = reshape2::melt(result.var) 
OLS.lower.unrec = reshape2::melt(fc.OLS.lwr) 
OLS.upper.unrec = reshape2::melt(fc.OLS.upr) 
OLS.upper.rec = OLS.rec$value + 1.96*sqrt(OLS.var.rec$value)
OLS.lower.rec = OLS.rec$value - 1.96*sqrt(OLS.var.rec$value)
date = rep(1:k, 1035)
Actual = reshape2::melt(ally.test)
error.rec = Actual$value- OLS.rec$value
error.unrec = Actual$value- OLS.unrec$value
fc.OLS <- cbind(OLS.unrec$value,
                OLS.rec$value, 
                Actual,  
                OLS.var.rec$value, 
                OLS.lower.unrec$value,
                OLS.upper.unrec$value, 
                OLS.lower.rec,
                OLS.upper.rec, 
                error.rec,
                error.unrec)
colnames(fc.OLS) <- c('OLS.unrec', 'OLS.rec','date','Series', 'Actual', 'OLS.var.rec', 'OLS.lower.unrec', 'OLS.upper.unrec', 
                      'OLS.lower.rec', 'OLS.upper.rec', 'error.rec', 'error.unrec')


##############
## ETS - ARIMA
##############

library(tsibble)
library(fabletools)
library(fable)
library(tidyverse)
library(lubridate)

# Convert two different date formats
make_date <- function(x) {
  numeric_dates <- !str_detect(x, "\\/")
  output <- Date(length = length(x))
  output[numeric_dates] <- as.Date(as.numeric(x[numeric_dates]), origin = "1899-12-30")
  output[!numeric_dates] <- mdy(x[!numeric_dates])
  return(output)
}

wikigts <- read_csv("wikipedia_data.csv") %>%
  select(date, views, cat_column) %>%
  mutate(date = make_date(date)) %>%
  group_by(date, cat_column) %>%
  summarise(views = mean(views)) %>%
  ungroup() %>%
  mutate(
    Access = str_sub(cat_column, 1, 7),
    Agent = str_sub(cat_column, 8, 9),
    Language = str_sub(cat_column, 10, 11),
    Purpose = str_sub(cat_column, 12, 14),
    Article = str_sub(cat_column, 15, 16),
  ) %>%
  as_tsibble(index = date, key = c(Access, Agent, Language, Purpose, Article)) %>%
  aggregate_key(
    Access + Agent + Language + Purpose +
      Access:Agent + Access:Language + Access:Purpose + Agent:Language + Agent:Purpose + Language:Purpose +
      Access:Agent:Language:Purpose:Article,
    views = sum(views)
  )


gts.rolling <- wikigts %>%
  filter(date < ymd("2017-06-29"))

gts.rolling2 <- gts.rolling %>%
  tsibble::stretch_tsibble(.init = 366 , .step = 1)

new_data <- wikigts %>%
  dplyr::filter(date >= ymd("2017-06-02")) %>%
  rename(actual = views)%>% 
  arrange(`date`) %>%
  mutate(new_index = dense_rank(date))

# ETS
start_time <- Sys.time()
fc.ets <- gts.rolling2 %>%
  model(ets = ETS(views))
end_time <- Sys.time()
end_time - start_time

m <- c(1:28)
fc.ets.rec <- data.frame(a=c(), b=c())
for(i in m){
  result <- fc.ets %>%
    filter(.id == i) %>%
    reconcile(ets_adjusted = min_trace(ets, method="wls_struct")) %>%
    forecast(h = 1) %>%
    hilo(level=95)%>%
    unpack_hilo("95%") 
  
  result <- result %>% 
    distinct(across(-views)) %>%
    mutate('h' = i)
  new_data2 <- new_data %>% 
    filter(new_index == i)
  
  fc.rec <- result %>%
    left_join(new_data2) %>%
    mutate(error = actual - .mean)
  
  fc.ets.rec <- bind_rows(fc.ets.rec, fc.rec)
}


# ARIMA
start_time <- Sys.time()
fc.arima <- gts.rolling2 %>%
  model(arima = ARIMA(views))
end_time <- Sys.time()
end_time - start_time

m <- c(1:28)
fc.arima.rec <- data.frame(a=c(), b=c())
for(i in m){
  result <- fc.arima %>%
    filter(.id == i) %>%
    reconcile(arima_adjusted = min_trace(arima, method="wls_struct")) %>%
    forecast(h = 1) %>%
    hilo(level=95)%>%
    unpack_hilo("95%") 
  
  result <- result %>% 
    distinct(across(-views)) %>%
    mutate('h' = i)
  new_data2 <- new_data %>% 
    filter(new_index == i)
  
  fc.rec <- result %>%
    left_join(new_data2) %>%
    mutate(error = actual - .mean)
  
  fc.arima.rec <- bind_rows(fc.arima.rec, fc.rec)
}

fc.ets.arima <- bind_rows (fc.arima.rec, fc.ets.rec)

##### Plotting the results

fc.OLS <- bind_rows(fc.OLS %>%
                      filter(Series == 'Total') %>%
                      mutate (Level = 'Total'),
                    fc.OLS %>% filter(Series %in% c('Access/desktop','Access/mobilea', 
                                                    'Access/mobilea')) %>%
                      mutate (Level = 'Access'), 
                    fc.OLS %>% filter(Series %in% c('Agent/sp','Agent/us')) %>%
                      mutate (Level = 'Agent'), 
                    fc.OLS %>% filter(Series %in% c('Language/de','Language/es', 
                                                    'Language/en', 'Language/zh')) %>%
                      mutate (Level = 'Language'), 
                    fc.OLS %>% filter(Series %in% c('Purpose/Blo', 'Purpose/Bus', 'Purpose/Gam',
                                                    'Purpose/Gen', 'Purpose/Lif', 'Purpose/Pho', 
                                                    'Purpose/Reu', 'Purpose/Tra', 'Purpose/Vid')) %>%
                      mutate (Level = 'Purpose'), 
                    fc.OLS %>% filter(grepl('Access x Purpose/', Series)) %>%
                      mutate (Level = 'Access x Purpose'), 
                    fc.OLS %>% filter(grepl('Agent x Purpose/', Series)) %>%
                      mutate (Level = 'Agent x Purpose'), 
                    fc.OLS %>% filter(grepl('Language x Purpose/', Series)) %>%
                      mutate (Level = 'Language x Purpose'),
                    fc.OLS %>% filter(grepl('Access x Agent/', Series)) %>%
                      mutate (Level = 'Access x Agent'),
                    fc.OLS %>% filter(grepl('Access x Language/', Series)) %>%
                      mutate (Level = 'Access x Language'),
                    fc.OLS %>% filter(grepl('Agent x Language/', Series)) %>%
                      mutate (Level = 'Agent x Language'),
                    fc.OLS %>% filter( !grepl('Access', Series) & !grepl('Agent', Series) & 
                                         !grepl('Language', Series) & !grepl('Purpose', Series) & 
                                         !grepl('Total', Series)) %>%
                      mutate (Level = 'Bottom level'))

error.wiki.OLS <-bind_rows( dplyr::select ( fc.OLS, error = error.rec, Level) %>% 
                              mutate(Rec = 'rec') %>%
                              mutate(Method = 'OLS') , dplyr::select ( fc.OLS, error = error.unrec, Level) %>% 
                              mutate(Rec = 'unrec') %>%
                              mutate(Method = 'OLS'))


fc.ets.arima <- bind_rows( fc.ets.arima %>%
                             filter(
                               is_aggregated(Access),
                               is_aggregated(Agent),
                               is_aggregated(Language),
                               is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Total'), 
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(Access),
                               is_aggregated(Agent),
                               is_aggregated(Language),
                               is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Access'), 
                           fc.ets.arima %>%
                             filter(
                               is_aggregated(Access),
                               !is_aggregated(Agent),
                               is_aggregated(Language),
                               is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Agent') ,
                           fc.ets.arima %>%
                             filter(
                               is_aggregated(Access),
                               is_aggregated(Agent),
                               !is_aggregated(Language),
                               is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Language') ,
                           fc.ets.arima %>%
                             filter(
                               is_aggregated(Access),
                               is_aggregated(Agent),
                               is_aggregated(Language),
                               !is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>%  mutate (Level = 'Purpose') ,
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(Access),
                               !is_aggregated(Agent),
                               is_aggregated(Language),
                               is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Access x Agent') ,
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(Access),
                               is_aggregated(Agent),
                               is_aggregated(Language),
                               !is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Access x Purpose') ,
                           fc.ets.arima %>%
                             filter(
                               is_aggregated(Access),
                               !is_aggregated(Agent),
                               is_aggregated(Language),
                               !is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Agent x Purpose'),
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(Access),
                               is_aggregated(Agent),
                               !is_aggregated(Language),
                               is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Access x Language'),
                           fc.ets.arima %>%
                             filter(
                               is_aggregated(Access),
                               !is_aggregated(Agent),
                               !is_aggregated(Language),
                               is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Agent x Language'),
                           fc.ets.arima %>%
                             filter(
                               is_aggregated(Access),
                               is_aggregated(Agent),
                               !is_aggregated(Language),
                               !is_aggregated(Purpose),
                               is_aggregated(Article)
                             ) %>% mutate (Level = 'Language x Purpose'), 
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(Access),
                               !is_aggregated(Agent),
                               !is_aggregated(Language),
                               !is_aggregated(Purpose),
                               !is_aggregated(Article)
                             ) %>% mutate (Level = 'Bottom level'))

fc.ets.arima <- bind_rows(fc.ets.arima %>% 
                            filter(.model %in% c('arima', 'ets')) %>%
                            mutate(Rec = 'unrec'), 
                          fc.ets.arima %>% 
                            filter(.model %in% c('arima_adjusted', 'ets_adjusted')) %>%
                            mutate(Rec = 'rec'))

fc.ets.arima <- bind_rows(fc.ets.arima %>% 
                            filter(.model %in% c('arima', 'arima_adjusted')) %>%
                            mutate(Method = 'ARIMA'), 
                          fc.ets.arima %>% 
                            filter(.model %in% c('ets', 'ets_adjusted')) %>%
                            mutate(Method = 'ETS'))



error.wiki <- bind_rows(
  dplyr::select ( fc.ets.arima, error , Level, Rec, Method), error.wiki.OLS) %>%
  filter(Level %in% c("Total", "Access", "Agent", "Language", "Purpose", "Bottom level")) %>%
  mutate( facet = factor(Level,
                         levels = c("Total", "Access", "Agent", "Language", "Purpose", "Bottom level")))


### Computing RMSE - table 20
rmse <- error.wiki %>%
  group_by(Rec, Method, facet) %>%
  summarise(
    rmse = sqrt(mean(error^2))
  ) %>%
  spread(value = rmse, key = Method) %>%
  ungroup() %>%
  select(Rec, facet, ETS, ARIMA, OLS) %>%
  mutate(
    facet = str_replace(facet, "level([0-9])", "facet \\1")
  )

### Plotting the results

m <- c(746, 450, 420, 576, 596, 631, 611, 514, 411, 452, 582, 927, 1106, 1313, 926, 422, 442, 619, 630, 605, 585, 606, 414, 405, 623, 614, 644, 653) 

forecast.wiki.OLS <- bind_rows(fc.OLS %>% filter(Series == 'Total'), fc.OLS %>% filter(Series == 'desktopusenPho') %>% filter(Actual%in% m) )
forecast.wiki.OLS <- forecast.wiki.OLS %>% select(-error.rec, -error.unrec, -Level) 

arima.unrec <- bind_rows( fc.ets.arima %>% filter(.model == 'arima' , is_aggregated(Access), 
                                                  is_aggregated(Agent), is_aggregated(Language), 
                                                  is_aggregated(Purpose), 
                                                  is_aggregated(Article) ) %>%
                            select(ARIMA.unrec = .mean, ARIMA.lower.unrec = `95%_lower`, ARIMA.upper.unrec = `95%_upper`),
                          fc.ets.arima %>% filter(.model == 'arima' , Access == 'desktop', 
                                                  Agent == 'us', Language == 'en', Purpose == 'Pho', Article == '21') %>%
                            select(ARIMA.unrec = .mean, ARIMA.lower.unrec = `95%_lower`, ARIMA.upper.unrec = `95%_upper`))

arima.rec <- bind_rows( fc.ets.arima %>% filter(.model == 'arima_adjusted' , is_aggregated(Access), 
                                                is_aggregated(Agent), is_aggregated(Language), 
                                                is_aggregated(Purpose), 
                                                is_aggregated(Article)) %>%
                          select(ARIMA.rec = .mean, ARIMA.lower.rec = `95%_lower`, ARIMA.upper.rec = `95%_upper`),
                        fc.ets.arima %>% filter(.model == 'arima_adjusted' , Access == 'desktop', 
                                                Agent == 'us', Language == 'en', Purpose == 'Pho', Article == '21') %>%
                          select(ARIMA.rec = .mean, ARIMA.lower.rec = `95%_lower`, ARIMA.upper.rec = `95%_upper`))

ets.unrec <- bind_rows( fc.ets.arima %>% filter(.model == 'ets' , is_aggregated(Access), 
                                                is_aggregated(Agent), is_aggregated(Language), 
                                                is_aggregated(Purpose), 
                                                is_aggregated(Article)) %>%
                          select(ETS.unrec = .mean, ETS.lower.unrec = `95%_lower`, ETS.upper.unrec = `95%_upper`),
                        fc.ets.arima %>% filter(.model == 'ets' , Access == 'desktop', 
                                                Agent == 'us', Language == 'en', Purpose == 'Pho', Article == '21') %>%
                          select(ETS.unrec = .mean, ETS.lower.unrec = `95%_lower`, ETS.upper.unrec = `95%_upper`))

ets.rec <- bind_rows( fc.ets.arima %>% filter(.model == 'ets_adjusted' , is_aggregated(Access), 
                                              is_aggregated(Agent), is_aggregated(Language), 
                                              is_aggregated(Purpose), 
                                              is_aggregated(Article)) %>%
                        select(ETS.rec = .mean, ETS.lower.rec = `95%_lower`, ETS.upper.rec = `95%_upper`),
                      fc.ets.arima %>% filter(.model == 'ets_adjusted' , Access == 'desktop', 
                                              Agent == 'us', Language == 'en', Purpose == 'Pho', Article == '21') %>%
                        select(ETS.rec = .mean, ETS.lower.rec = `95%_lower`, ETS.upper.rec = `95%_upper`))
forecast.wiki.data <- bind_cols (forecast.wiki.OLS, arima.unrec, arima.rec, ets.unrec, ets.rec)

forecast.wiki <- forecast.wiki.data %>%
  select(-OLS.lower.rec, -OLS.upper.rec, -OLS.lower.unrec, -OLS.upper.unrec,
         -ARIMA.lower.rec, -ARIMA.upper.rec, -ARIMA.lower.unrec, -ARIMA.upper.unrec,
         -ETS.lower.rec, -ETS.upper.rec, -ETS.lower.unrec, -ETS.upper.unrec, - OLS.var.rec) %>%
  gather(-Series, -date , key = "Method", value = "Count") %>%
  mutate(
    Rec = str_extract(Method, "[a-z]*$"),
    Rec = if_else(Rec == "ctual", "Actual", Rec),
    Model = str_extract(Method, "^[A-Z]*"),
    Model = if_else(Model == "A", "Actual", Model)
  )

## error boxplot

## remove the outliers in boxplots
boxplot.stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

## Figure 23
error.wiki %>%
  mutate(id = factor(paste(Method, Rec, sep = "."),
                     levels = c("ETS.rec", "ETS.unrec", "ARIMA.rec", "ARIMA.unrec", "OLS.rec", "OLS.unrec"),
                     labels = c("ETS.rec", "ETS.unrec", "ARIMA.rec", "ARIMA.unrec", "OLS.rec", "OLS.unrec")
  )) %>%
  ggplot(aes(x = id, y = error, fill = id)) +
  stat_summary(fun.data = boxplot.stat, geom = "boxplot", alpha = 0.5) +
  xlab("Method") +
  ylab("Error") +
  facet_wrap(~facet, ncol = 3, scales = "free_y") +
  guides(fill = guide_legend(nrow = 1, bycol = TRUE)) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "ETS.rec" = "green",
    "ETS.unrec" = "lightgreen",
    "ARIMA.rec" = "blue",
    "ARIMA.unrec" = "lightblue",
    "OLS.rec" = "pink4",
    "OLS.unrec" = "pink"
  )) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

## Sample series plot with forecasts

### Total - figure 25
ylim <- forecast.wiki %>%
  filter(Series == "Total") %>%
  pull(Count) %>%
  range()
forecast.wiki %>%
  filter(Series == "Total") %>%
  ggplot(aes(x = date, y = Count, colour = Model, linetype = Rec, size = Model)) +
  geom_line() +
  ylim(ylim) +
  xlab("Horizon") +
  ylab("Count") +
  ggtitle("Rolling origin 1-step forecasts") +
  scale_linetype_manual(name = "Reconciled", values = c(Actual = "solid", rec = "dashed", unrec = "dotted")) +
  scale_size_manual(values = c(Actual = 0.8, ETS = 0.5, ARIMA = 0.5, OLS = 0.5), guide = "none") +
  scale_color_manual(
    name = "Series",
    values = c(
      ETS = "green",
      ARIMA = "blue",
      OLS = "red",
      Actual = "black"
    )
  ) +
  theme_bw()
### desktopusenPho - figure 26
ylim <- forecast.wiki %>%
  filter(Series == "desktopusenPho") %>%
  pull(Count) %>%
  range()

forecast.wiki %>%
  filter(Series == "desktopusenPho") %>%
  ggplot(aes(x = date, y = Count, colour = Model, linetype = Rec, size = Model)) +
  geom_line() +
  ylim(ylim) +
  xlab("Horizon") +
  ylab("Count") +
  ggtitle("Rolling origin 1-step forecasts") +
  scale_linetype_manual(name = "Reconciled", values = c(Actual = "solid", rec = "dashed", unrec = "dotted")) +
  scale_size_manual(values = c(Actual = 0.8, ETS = 0.5, ARIMA = 0.5, OLS = 0.5), guide = "none") +
  scale_color_manual(
    name = "Series",
    values = c(
      ETS = "green",
      ARIMA = "blue",
      OLS = "red",
      OLSX = "yellow",
      ARIMAX = "orchid1",
      Actual = "black"
    )
  ) +
  theme_bw() 
## Sample series prediction interval plots
### Total - figure 27
forecast.wiki.data %>%
  filter(Series == "Total") %>%
  ggplot(aes(x = date, y = Actual, colour = "Actual", size = 'Actual')) +
  geom_ribbon(aes(x = date, ymax = ARIMA.upper.rec, ymin = ARIMA.lower.rec), fill = "lightskyblue", colour = "lightskyblue", alpha = .2, size = 0.5) +
  geom_ribbon(aes(x = date, ymax = ETS.upper.rec, ymin = ETS.lower.rec), fill = "lightgreen", colour = "lightgreen", alpha = .2, size = 0.5) +
  geom_ribbon(aes(x = date, ymax = OLS.upper.rec, ymin = OLS.lower.rec), fill = "pink", colour = "pink", alpha = .2, size = 0.5)  +
  geom_line(aes(y = ARIMA.rec, colour = "ARIMA.rec", size = 'ARIMA.rec')) +
  geom_line(aes(y = ETS.rec, colour = "ETS.rec", size = 'ETS.rec')) +
  geom_line(aes(y = OLS.rec, colour = "OLS.rec", size = 'OLS.rec')) +
  geom_line() +
  xlab("Horizon") +
  ylab("Count") +
  ggtitle("Rolling origin 1-step forecasts") +
  scale_colour_manual("Method",
                      breaks = c("Actual", "ARIMA.rec",  "ETS.rec", "OLS.rec"),
                      values = c("black", "blue", "green",  "red"))+
  scale_size_manual(breaks = c("Actual", "ARIMA.rec",  "ETS.rec", "OLS.rec"),
                    values = c( 0.8, 0.5,  0.5,  0.5), guide = "none") +
  theme_bw()
### desktopusenPho - figure 28
forecast.wiki.data %>%
  filter(Series == "desktopusenPho", Actual %in% m) %>%
  ggplot(aes(x = date, y = Actual, colour = "Actual", size = 'Actual')) +
  geom_ribbon(aes(x = date, ymax = ARIMA.upper.rec, ymin = ARIMA.lower.rec), fill = "lightskyblue", colour = "lightskyblue", alpha = .2, size = 0.5) +
  geom_ribbon(aes(x = date, ymax = ETS.upper.rec, ymin = ETS.lower.rec), fill = "lightgreen", colour = "lightgreen", alpha = .2, size = 0.5) +
  geom_ribbon(aes(x = date, ymax = OLS.upper.rec, ymin = OLS.lower.rec), fill = "pink", colour = "pink", alpha = .2, size = 0.5)  +
  geom_line(aes(y = ARIMA.rec, colour = "ARIMA.rec", size = 'ARIMA.rec')) +
  geom_line(aes(y = ETS.rec, colour = "ETS.rec", size = 'ETS.rec')) +
  geom_line(aes(y = OLS.rec, colour = "OLS.rec", size = 'OLS.rec')) +
  geom_line() +
  xlab("Horizon") +
  ylab("Count") +
  ggtitle("Rolling origin 1-step forecasts") +
  scale_colour_manual("Method",
                      breaks = c("Actual", "ARIMA.rec",  "ETS.rec", "OLS.rec"),
                      values = c("black", "blue", "green",  "red"))+
  scale_size_manual(breaks = c("Actual", "ARIMA.rec",  "ETS.rec", "OLS.rec"),
                    values = c( 0.8, 0.5,  0.5,  0.5), guide = "none") +
  theme_bw()




