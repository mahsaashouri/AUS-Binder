
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
aus <-ts(readr::read_csv("TourismData_v3.csv")[-c(1,2)],start=1,frequency =12)
## creating hierarchy structure
ausgts <- gts(aus, characters = list(c(1, 1, 1), 3),
              gnames = c("State", "Zone", "Region", "Purpose","State x Purpose", "Zone x Purpose"))
## forecast horizon
h <- 24
n <- nrow(aus)
## train and test sets
train_tourist <- window(ausgts,start = c(1, 1),end = c(1, (n-h)))
test_tourist <- window(ausgts,start = c(1, ((n-h)+1)),end = c(1, n))
ally <- aggts(ausgts)
ally.test <- aggts(test_tourist)
## empty array for forecast results
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
  fit.OLS <- OLSmodel(ally[,i],12,12,h, nolag = c(1,12))
  fc[, i,"OLS"] <- fit.OLS[[1]]
  fc[, i,"OLS.lwr"] <- fit.OLS[[2]]
  fc[, i,"OLS.upr"] <- fit.OLS[[3]]
  fc[, i,"OLS.se"] <- fit.OLS[[4]]
  fc[, i,"OLS.residual.scale"] <- fit.OLS[[5]]
}
end_time <- Sys.time()
end_time - start_time

fc.OLS.base <- as.data.frame(fc[,,"OLS"])
fc.OLS.lwr <- as.data.frame(fc[,,"OLS.lwr"])
fc.OLS.upr <- as.data.frame(fc[,,"OLS.upr"])
fc.OLS.se <- as.data.frame(fc[,,"OLS.se"])
fc.OLS.residual.scale <- as.data.frame(fc[,,"OLS.residual.scale"])


## PI variance
fc.OLS.PI <- (fc.OLS.se)^2 + (fc.OLS.residual.scale)^2

## making matrix computations sparse from Matrix package
as.matrix <- Matrix::as.matrix
t <- Matrix::t
solve <- Matrix::solve
diag <- Matrix::diag

## computing reconciliation matrix
gmat <- GmatrixG(ausgts$groups)
smatrix <- as((SmatrixM(gmat)), 'dgCMatrix')
lambda <- as(diag(rowSums(smatrix)), 'dgCMatrix')

rec.adj.lambda <- as.matrix(smatrix%*%solve(t(smatrix)%*%solve(lambda)%*%smatrix)%*%t(smatrix)%*%solve(lambda))

## computing reconciled forecasts
fc.rec <- matrix(NA, nrow = h, ncol = ncol(ally))
for(i in 1:nrow(fc.OLS.base)){
  f.1 <- matrix(as.numeric(fc.OLS.base[i,]), ncol = 1, nrow = ncol(fc.OLS.base))
  fc.rec [i,] <- rec.adj.lambda %*% f.1
}
colnames(fc.rec) <- colnames(ally)

##################
### Computing reconciled variance for prediction intervals
##################

## computing predictors (trend, dummy seasonality, lags) for each series
Xmat<-list()
freq <-12
maxlag <- 12
ally.final <- as.list(ally) 
Xmatrix<-function(X){
  X<-as.vector(X)
  intercept <- rep(1, length(X))
  trend<-seq(NROW(X))
  season<-forecast::seasonaldummy(ts(X,frequency = freq))
  Xlag<-quantmod::Lag(X,k=1:maxlag)
  X_mat<-cbind.data.frame(intercept,trend,season,Xlag)
  Xmat[[length(Xmat)+1]] <- X_mat 
}
Xmat.final <- lapply(ally.final, Xmatrix)

result.var <- matrix(NA, nrow = h, ncol = NCOL(ally))
## for loops for computing reconciled variance
for(i in 1:h){
  ## train and test sets matrix
  Xmat.final.train <- lapply(Xmat.final, function(x)x[1:((n - h) + (i - 1)),])
  Xmat.final.test <- lapply(Xmat.final, function(x)x[(n - h) + i,])
  mat <- as(Matrix::bdiag(lapply(Xmat.final.train, function(x){as.matrix(na.omit(x))})), 'dgCMatrix')
  ## (X'X)^-1
  mat.inverse <- as(Matrix::solve(t(mat)%*%mat), 'dgCMatrix')
  ## X*_t+h
  mat.test <- as(Matrix::bdiag(lapply(Xmat.final.test, function(x){as.matrix(na.omit(x))})), 'dgCMatrix')
  ## Formula 7 in the text
  H.matrix <- as(mat.test %*% mat.inverse %*% t(mat.test), 'dgCMatrix')
  Sigma.mat <-as( Matrix::diag(fc.OLS.PI[i,]) + (Matrix::diag(fc.OLS.PI[i,]) %*% H.matrix), 'dgCMatrix')
  rec.p <- as(as.matrix(solve(t(smatrix)%*%solve(lambda)%*%smatrix)%*%t(smatrix)%*%solve(lambda)), 'dgCMatrix')
  var.for <- as.matrix((smatrix %*% rec.p) %*% Sigma.mat %*% (t(rec.p) %*% t(smatrix)))
  result.var[i,] <- as.vector(diag(var.for))
}

## saving the forecasts, errors and prediction interval results - based and reconciled forecasts
OLS.unrec = reshape2::melt(fc.OLS.base)
OLS.rec = reshape2::melt(fc.rec) 
OLS.var.rec = reshape2::melt(result.var) 
OLS.lower.unrec = reshape2::melt(fc.OLS.lwr) 
OLS.upper.unrec = reshape2::melt(fc.OLS.upr) 
Actual =  reshape2::melt(ally.test)$value
error.rec = Actual- OLS.rec$value
error.unrec = Actual- OLS.unrec$value
OLS.upper.rec = OLS.rec$value + 1.96*sqrt(OLS.var.rec$value)
OLS.lower.rec = OLS.rec$value - 1.96*sqrt(OLS.var.rec$value)
date = rep(1:h, 555)

fc.OLS <- cbind(OLS.unrec$value,
                OLS.rec$value, 
                Actual,  
                OLS.var.rec$value, 
                OLS.lower.unrec,
                OLS.upper.unrec$value, 
                OLS.lower.rec,
                OLS.upper.rec, 
                error.rec,
                error.unrec, date)
colnames(fc.OLS) <- c('OLS.unrec', 'OLS.rec', 'Actual', 'OLS.var.rec', 'Series', 'OLS.lower.unrec', 'OLS.upper.unrec', 
                      'OLS.lower.rec', 'OLS.upper.rec', 'error.rec', 'error.unrec', 'date')

##############
## ETS - ARIMA
##############

library(fable)
library(fabletools)
library(tidyverse)
library(tsibble)
library(readr)

aus.ets.arima <- read.csv('TourismData_v3.csv', header = TRUE)[,-c(1,2)]
aus.ets.arima <-  tibble(aus.ets.arima)
aus.ets.arima$Date <-  rep(yearmonth("1998 Jan") + 0:227)

aus.ets.arima <- aus.ets.arima %>%
  pivot_longer(-Date, names_to = "group", values_to = "value") %>%
  mutate(
    State = stringr::str_sub(group, 1, 1),
    Zone = stringr::str_sub(group, 1, 2),
    Region = stringr::str_sub(group, 1, 3),
    Purpose = stringr::str_sub(group, 4, 6),
  ) %>%
  select(-group) %>%
  as_tsibble(index = Date, key=c(State, Zone, Region, Purpose))


ausgts.ets.arima <- aus.ets.arima %>%
  aggregate_key(Purpose * (State/ Zone/ Region), value = sum(value)) 

new_data <- ausgts.ets.arima %>%
  dplyr::filter(Date > yearmonth ("2014 Dec")) %>%
  rename(actual = value)

## rolling window
gts.rolling <- ausgts.ets.arima %>%
  filter(Date < yearmonth ("2016 Dec")) %>%
  stretch_tsibble(.init = 204 , .step = 1)

new_data <- ausgts.ets.arima %>%
  dplyr::filter(Date > yearmonth ("2014 Dec")) %>%
  rename(actual = value)%>% 
  arrange(`Date`) %>%
  mutate(new_index = dense_rank(Date))

# ETS
## computation time
start_time <- Sys.time()
fc.ets <- gts.rolling %>%
  model(ets = ETS(value))
end_time <- Sys.time()
end_time - start_time

m <- c(1:24)
fc.ets.rec <- data.frame(a=c(), b=c())
for(i in m){
  result <- fc.ets %>%
    filter(.id == i) %>%
    reconcile(ets_adjusted = min_trace(ets, method="wls_struct")) %>%
    forecast(h = 1) %>%
    hilo(level=95)%>%
    unpack_hilo("95%") 
  
  result <- result %>% 
    distinct(across(-value)) %>%
    mutate('h' = i)
  new_data2 <- new_data %>% 
    filter(new_index == i)
  
  fc.rec <- result %>%
    left_join(new_data2) %>%
    mutate(error = actual - .mean)
  
  fc.ets.rec <- bind_rows(fc.ets.rec, fc.rec)
}


# ARIMA
## computation time
start_time <- Sys.time()
fc.arima <- gts.rolling %>%
  model(arima = ARIMA(value))
end_time <- Sys.time()
end_time - start_time

m <- c(1:24)
fc.arima.rec <- data.frame(a=c(), b=c())
for(i in m){
  result <- fc.arima %>%
    filter(.id == i) %>%
    reconcile(arima_adjusted = min_trace(arima, method="wls_struct")) %>%
    forecast(h = 1) %>%
    hilo(level=95)%>%
    unpack_hilo("95%") 
  
  result <- result %>% 
    distinct(across(-value)) %>%
    mutate('h' = i)
  new_data2 <- new_data %>% 
    filter(new_index == i)
  
  fc.rec <- result %>%
    left_join(new_data2) %>%
    mutate(error = actual - .mean)
  
  fc.arima.rec <- bind_rows(fc.arima.rec, fc.rec)
}

## saving  ets and arima results

fc.ets.arima <- bind_rows (fc.arima.rec, fc.ets.rec)

##### Plotting the results

fc.OLS <- bind_rows(fc.OLS %>%
                      filter(Series == 'Total') %>%
                      mutate (Level = 'Total'),
                    fc.OLS %>% filter(grepl('State/', Series)) %>%
                      mutate (Level = 'State'), 
                    fc.OLS %>% filter(grepl('Zone/', Series)) %>%
                      mutate (Level = 'Zone'), 
                    fc.OLS %>% filter(grepl('Region/', Series)) %>%
                      mutate (Level = 'Region'), 
                    fc.OLS %>% filter(Series %in% c('Purpose/Bus','Purpose/Hol', 
                                                   'Purpose/Oth','Purpose/Vis')) %>%
                      mutate (Level = 'Purpose'), 
                    fc.OLS %>% filter(grepl('State x Purpose/', Series)) %>%
                      mutate (Level = 'State x Purpose'), 
                    fc.OLS %>% filter(grepl('Zone x Purpose/', Series)) %>%
                      mutate (Level = 'Zone x Purpose'), 
                    fc.OLS %>% filter( !grepl('State', Series) & !grepl('Zone', Series) & 
                                        !grepl('Region', Series) & !grepl('Purpose', Series) & 
                                        !grepl('Total', Series)) %>%
                      mutate (Level = 'Region x Purpose'))

error.fc.OLS<-bind_rows( dplyr::select ( fc.OLS, error = error.rec, Level) %>% 
                           mutate(Rec = 'rec') %>%
                           mutate(Method = 'OLS') , dplyr::select ( fc.OLS, error = error.unrec, Level) %>% 
                           mutate(Rec = 'unrec') %>%
                           mutate(Method = 'OLS'))

fc.ets.arima <- bind_rows( fc.ets.arima %>%
                             filter(
                               is_aggregated(State),
                               is_aggregated(Zone),
                               is_aggregated(Region),
                               is_aggregated(Purpose)
                             ) %>% mutate (Level = 'Total'), 
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(State),
                               is_aggregated(Zone),
                               is_aggregated(Region),
                               is_aggregated(Purpose)
                             ) %>% mutate (Level = 'State'), 
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(State),
                               !is_aggregated(Zone),
                               is_aggregated(Region),
                               is_aggregated(Purpose)
                             ) %>% mutate (Level = 'Zone') ,
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(State),
                               !is_aggregated(Zone),
                               !is_aggregated(Region),
                               is_aggregated(Purpose)
                             ) %>% mutate (Level = 'Region') ,
                           fc.ets.arima %>%
                             filter(
                               is_aggregated(State),
                               is_aggregated(Zone),
                               is_aggregated(Region),
                               !is_aggregated(Purpose)
                             ) %>%  mutate (Level = 'Purpose') ,
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(State),
                               is_aggregated(Zone),
                               is_aggregated(Region),
                               !is_aggregated(Purpose)
                             ) %>% mutate (Level = 'State x Purpose') ,
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(State),
                               !is_aggregated(Zone),
                               is_aggregated(Region),
                               !is_aggregated(Purpose)
                             ) %>% mutate (Level = 'Zone x Purpose') ,
                           fc.ets.arima %>%
                             filter(
                               !is_aggregated(State),
                               !is_aggregated(Zone),
                               !is_aggregated(Region),
                               !is_aggregated(Purpose)
                             ) %>% mutate (Level = 'Region x Purpose') )

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



error.tourism <- bind_rows(
  dplyr::select ( fc.ets.arima, error , Level, Rec, Method), error.fc.OLS) %>%
  mutate( facet = factor(Level,
                         levels = c("Total", "State", "Zone", "Region", "Purpose", "State x Purpose", "Zone x Purpose", "Region x Purpose")))

### Computing RMSE - table 3
rmse <- error.tourism %>%
  group_by(Rec, Method, facet) %>%
  summarise(
    rmse = sqrt(mean(error^2))
  ) %>%
  spread(value = rmse, key = Method) %>%
  ungroup() %>%
  select(Rec, facet, ETS, ARIMA, OLS) %>%
  mutate(facet = str_replace(facet, "level([0-9])", "facet \\1"))

### Plotting the results

forecast.tourism.OLS <- bind_rows(fc.OLS %>% filter(Series == 'Total') , fc.OLS %>% filter(Series == 'AAAVis'))
forecast.tourism.OLS <- forecast.tourism.OLS %>% select(-error.rec, -error.unrec, -Level) 

arima.unrec <- bind_rows( fc.ets.arima %>% filter(.model == 'arima' , is_aggregated(State), is_aggregated(Zone), is_aggregated(Region), is_aggregated(Purpose)) %>%
                            select(ARIMA.unrec = .mean, ARIMA.lower.unrec = `95%_lower`, ARIMA.upper.unrec = `95%_upper`),
                          fc.ets.arima %>% filter(.model == 'arima' , State == 'A', 
                                                  Zone == 'AA', Region == 'AAA', Purpose == 'Vis') %>%
                            select(ARIMA.unrec = .mean, ARIMA.lower.unrec = `95%_lower`, ARIMA.upper.unrec = `95%_upper`))

arima.rec <- bind_rows( fc.ets.arima %>% filter(.model == 'arima_adjusted' , is_aggregated(State), is_aggregated(Zone), is_aggregated(Region), is_aggregated(Purpose)) %>%
                          select(ARIMA.rec = .mean, ARIMA.lower.rec = `95%_lower`, ARIMA.upper.rec = `95%_upper`),
                        fc.ets.arima %>% filter(.model == 'arima_adjusted' , State == 'A', 
                                                Zone == 'AA', Region == 'AAA', Purpose == 'Vis') %>%
                          select(ARIMA.rec = .mean, ARIMA.lower.rec = `95%_lower`, ARIMA.upper.rec = `95%_upper`))

ets.unrec <- bind_rows( fc.ets.arima %>% filter(.model == 'ets' , is_aggregated(State), is_aggregated(Zone), is_aggregated(Region), is_aggregated(Purpose)) %>%
                          select(ETS.unrec = .mean, ETS.lower.unrec = `95%_lower`, ETS.upper.unrec = `95%_upper`), 
                        fc.ets.arima %>% filter(.model == 'ets' , State == 'A', 
                                                Zone == 'AA', Region == 'AAA', Purpose == 'Vis') %>%
                          select(ETS.unrec = .mean, ETS.lower.unrec = `95%_lower`, ETS.upper.unrec = `95%_upper`))

ets.rec <- bind_rows( fc.ets.arima %>% filter(.model == 'ets_adjusted' , is_aggregated(State), is_aggregated(Zone), is_aggregated(Region), is_aggregated(Purpose)) %>%
                        select(ETS.rec = .mean, ETS.lower.rec = `95%_lower`, ETS.upper.rec = `95%_upper`), 
                      fc.ets.arima %>% filter(.model == 'ets_adjusted' , State == 'A', 
                                              Zone == 'AA', Region == 'AAA', Purpose == 'Vis') %>%
                        select(ETS.rec = .mean, ETS.lower.rec = `95%_lower`, ETS.upper.rec = `95%_upper`))
forecast.tourism.data <- bind_cols (forecast.tourism.OLS, arima.unrec, arima.rec, ets.unrec, ets.rec)

forecast.tourism <- forecast.tourism.data %>%
  select( -OLS.lower.rec, -OLS.upper.rec, -OLS.lower.unrec, -OLS.upper.unrec,
          -ARIMA.lower.rec, -ARIMA.upper.rec, -ARIMA.lower.unrec, -ARIMA.upper.unrec,
          -ETS.lower.rec, -ETS.upper.rec, -ETS.lower.unrec, -ETS.upper.unrec, - OLS.var.rec) %>%
  gather(-Series, -date, key = "Method", value = "Count") %>%
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

###### Figure 5
error.tourism %>%
  mutate(id = factor(paste(Method, Rec, sep = "."),
                     levels = c("ETS.rec", "ETS.unrec", "ARIMA.rec", "ARIMA.unrec", "OLS.rec", "OLS.unrec"),
                     labels = c("ETS.rec", "ETS.unrec", "ARIMA.rec", "ARIMA.unrec", "OLS.rec", "OLS.unrec")
  )) %>%
  ggplot(aes(x = id, y = error, fill = id)) +
  stat_summary(fun.data = boxplot.stat, geom = "boxplot", alpha = 0.5) +
  xlab("Method") +
  ylab("Error") +
  facet_wrap(~facet, ncol = 4, scales = "free_y") +
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

### Total - figure 7
ylim <- forecast.tourism %>%
  filter(Series == "Total") %>%
  pull(Count) %>%
  range()

forecast.tourism %>%
  filter(Series == "Total") %>%
  ggplot(aes(x = date, y = Count, colour = Model, linetype = Rec, size = Model)) +
  geom_line() +
  ylim(ylim) +
  xlab("Horizon") +
  ylab("Count") +
  ggtitle("Fixed origin multi-step forecasts") +
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
### AAAVis - figure 8
ylim <- forecast.tourism %>%
  filter(Series == "AAAVis") %>%
  pull(Count) %>%
  range()

forecast.tourism %>%
  filter(Series == "AAAVis") %>%
  ggplot(aes(x = date, y = Count, colour = Model, linetype = Rec, size = Model)) +
  geom_line() +
  ylim(ylim) +
  xlab("Horizon") +
  ylab("Count") +
  ggtitle("Fixed origin multi-step forecasts") +
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

## Sample series prediction interval plots
### Total - figure 9
forecast.tourism.data %>%
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
  ggtitle("Fixed origin multi-step forecasts") +
  scale_colour_manual("Method",
                      breaks = c("Actual", "ARIMA.rec",  "ETS.rec", "OLS.rec"),
                      values = c("black", "blue", "green",  "red"))+
  scale_size_manual(breaks = c("Actual", "ARIMA.rec",  "ETS.rec", "OLS.rec"),
                    values = c( 0.8, 0.5,  0.5,  0.5), guide = "none") +
  theme_bw()
### AAAVis - figure 10
forecast.tourism.data %>%
  filter(Series == "AAAVis") %>%
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
  ggtitle("Fixed origin multi-step forecasts") +
  scale_colour_manual("Method",
                      breaks = c("Actual", "ARIMA.rec",  "ETS.rec", "OLS.rec"),
                      values = c("black", "blue", "green",  "red"))+
  scale_size_manual(breaks = c("Actual", "ARIMA.rec",  "ETS.rec", "OLS.rec"),
                    values = c( 0.8, 0.5,  0.5,  0.5), guide = "none") +
  theme_bw()



######################################
#### different reconciliation results
######################################


####################################
#### Reconciling forecasts - compute SP matrix based on the desired reconciliation type  - OLS
###################################

#mint_shrink - OLS
error.unrec.OLS <- as.matrix(ally.test) - as.matrix(fc.OLS.base)
n <- nrow(error.unrec.OLS)
covm <- crossprod(stats::na.omit(error.unrec.OLS)) / n
tar <- diag(apply(error.unrec.OLS, 2, compose(crossprod, stats::na.omit))/n)
corm <- cov2cor(covm)
xs <- scale(error.unrec.OLS, center = FALSE, scale = sqrt(diag(covm)))
xs <- xs[stats::complete.cases(xs),]
v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
diag(v) <- 0
corapn <- cov2cor(tar)
d <- (corm - corapn)^2
lambda <- sum(v)/sum(d)
lambda <- max(min(lambda, 1), 0)
W <- lambda * tar + (1 - lambda) * covm
gmat <- GmatrixG(ausgts$groups)
smatrix <- as.matrix(SmatrixM(gmat))
R <- t(smatrix)%*%solve(W)
P <- Matrix::solve(R%*%smatrix)%*%R
SP <- as.matrix(smatrix%*%P)


## multiply the based forecasts by SP matrix

fc.mint.shrink.OLS <- matrix(NA, nrow = 24, ncol = ncol(ally))

for(i in 1:nrow(fc.OLS.base)){
  f.1 <- matrix(as.numeric(fc.OLS.base[i,]), ncol = 1, nrow = ncol(fc.OLS.base))
  fc.mint.shrink.OLS [i,] <- SP %*% f.1
}
colnames(fc.mint.shrink.OLS) <- colnames(ally)

OLS.rec.mint.shrink <- reshape2::melt(fc.mint.shrink.OLS) %>%
  mutate(actual = reshape2::melt(ally.test)$value, error = reshape2::melt(ally.test)$value - value , Rec = 'mint_shrink')

## wls_var

W <- diag(diag(covm))
gmat <- GmatrixG(ausgts$groups)
smatrix <- as.matrix(SmatrixM(gmat))
R <- t(smatrix)%*%solve(W)
P <- solve(R%*%smatrix)%*%R
SP <- smatrix%*%P

fc.wls.var.OLS <- matrix(NA, nrow = 24, ncol = ncol(ally))

for(i in 1:nrow(fc.OLS.base)){
  f.1 <- matrix(as.numeric(fc.OLS.base[i,]), ncol = 1, nrow = ncol(fc.OLS.base))
  fc.wls.var.OLS [i,] <- SP %*% f.1
}
colnames(fc.wls.var.OLS) <- colnames(ally)

OLS.rec.wls.var <- reshape2::melt(fc.wls.var.OLS) %>%
  mutate(actual = reshape2::melt(ally.test)$value, error = reshape2::melt(ally.test)$value - value , Rec = 'wls_var')
## saving results
fc.OLS.mint.shrink.wls.var <- bind_rows(OLS.rec.mint.shrink, OLS.rec.wls.var)

###############################
### different reconciliation - ETS-ARIMA
##############################

### mint_shrink reconciliation
fc.ets.mint.shrink <- fc.ets %>%
  reconcile(ets_adjusted = min_trace(ets, method="mint_shrink"))%>%
  forecast(h = "2 years") 

fc.ets.mint.shrink.error <- fc.ets.mint.shrink %>%
  left_join(new_data) %>%
  mutate(error = actual - .mean)

fc.ets.mint.shrink <- fc.ets.mint.shrink.error %>%
  hilo(level=95) %>% 
  unpack_hilo("95%")

### wls_var reconciliation
fc.ets.wls.var <- fc.ets %>%
  reconcile(ets_adjusted = min_trace(ets, method="wls_var"))%>%
  forecast(h = "2 years") 

fc.ets.wls.var.error <- fc.ets.wls.var %>%
  left_join(new_data) %>%
  mutate(error = actual - .mean)

fc.ets.wls.var <- fc.ets.wls.var.error %>%
  hilo(level=95) %>% 
  unpack_hilo("95%")

#ARIMA

### mint_shrink reconciliation
fc.arima.mint.shrink <- fc.arima %>%
  reconcile(arima_adjusted = min_trace(arima, method="mint_shrink"))%>%
  forecast(h = "2 years") 

fc.arima.mint.shrink.error <- fc.arima.mint.shrink %>%
  left_join(new_data) %>%
  mutate(error = actual - .mean)

fc.arima.mint.shrink <- fc.arima.mint.shrink.error %>%
  hilo(level=95) %>% 
  unpack_hilo("95%")

### wls_var reconciliation
fc.arima.wls.var <- fc.arima %>%
  reconcile(arima_adjusted = min_trace(arima, method="wls_var"))%>%
  forecast(h = "2 years") 

fc.arima.wls.var.error <- fc.arima.wls.var %>%
  left_join(new_data) %>%
  mutate(error = actual - .mean)

fc.arima.wls.var <- fc.arima.wls.var.error %>%
  hilo(level=95) %>% 
  unpack_hilo("95%")
## mint_shrink & wls_var results
ets.arima.dif.rec <- bind_rows(bind_rows (fc.arima.mint.shrink, fc.ets.mint.shrink) %>% 
                                 distinct(across(-value)) %>%
                                 mutate(Rec = 'mint_shrink') %>%
                                 filter(.model %in% c('ets_adjusted', 'arima_adjusted')), bind_rows (fc.arima.wls.var, fc.ets.wls.var) %>% 
                                 distinct(across(-value)) %>%
                                 mutate(Rec = 'wls_var') %>%
                                 filter(.model %in% c('ets_adjusted', 'arima_adjusted')))


ets.arima.dif.rec <- bind_rows( ets.arima.dif.rec %>%
                                  filter(
                                    is_aggregated(State),
                                    is_aggregated(Zone),
                                    is_aggregated(Region),
                                    is_aggregated(Purpose)
                                  ) %>% mutate (Level = 'Total'), 
                                ets.arima.dif.rec %>%
                                  filter(
                                    !is_aggregated(State),
                                    is_aggregated(Zone),
                                    is_aggregated(Region),
                                    is_aggregated(Purpose)
                                  ) %>% mutate (Level = 'State'), 
                                ets.arima.dif.rec %>%
                                  filter(
                                    !is_aggregated(State),
                                    !is_aggregated(Zone),
                                    is_aggregated(Region),
                                    is_aggregated(Purpose)
                                  ) %>% mutate (Level = 'Zone') ,
                                ets.arima.dif.rec %>%
                                  filter(
                                    !is_aggregated(State),
                                    !is_aggregated(Zone),
                                    !is_aggregated(Region),
                                    is_aggregated(Purpose)
                                  ) %>% mutate (Level = 'Region') ,
                                ets.arima.dif.rec %>%
                                  filter(
                                    is_aggregated(State),
                                    is_aggregated(Zone),
                                    is_aggregated(Region),
                                    !is_aggregated(Purpose)
                                  ) %>%  mutate (Level = 'Purpose') ,
                                ets.arima.dif.rec %>%
                                  filter(
                                    !is_aggregated(State),
                                    is_aggregated(Zone),
                                    is_aggregated(Region),
                                    !is_aggregated(Purpose)
                                  ) %>% mutate (Level = 'State x Purpose') ,
                                ets.arima.dif.rec %>%
                                  filter(
                                    !is_aggregated(State),
                                    !is_aggregated(Zone),
                                    is_aggregated(Region),
                                    !is_aggregated(Purpose)
                                  ) %>% mutate (Level = 'Zone x Purpose') ,
                                ets.arima.dif.rec %>%
                                  filter(
                                    !is_aggregated(State),
                                    !is_aggregated(Zone),
                                    !is_aggregated(Region),
                                    !is_aggregated(Purpose)
                                  ) %>% mutate (Level = 'Region x Purpose') )


ets.arima.dif.rec <- bind_rows(ets.arima.dif.rec %>% 
                                 filter(.model %in% c('arima_adjusted')) %>%
                                 mutate(Method = 'ARIMA'), 
                               ets.arima.dif.rec %>% 
                                 filter(.model %in% c('ets_adjusted')) %>%
                                 mutate(Method = 'ETS'))


fc.OLS.mint.shrink.wls.var <- bind_rows(fc.OLS.mint.shrink.wls.var %>%
                                          filter(Var2 == 'Total') %>%
                                          mutate (Level = 'Total'),
                                        fc.OLS.mint.shrink.wls.var %>% filter(grepl('State/', Var2)) %>%
                                          mutate (Level = 'State'), 
                                        fc.OLS.mint.shrink.wls.var %>% filter(grepl('Zone/', Var2)) %>%
                                          mutate (Level = 'Zone'), 
                                        fc.OLS.mint.shrink.wls.var %>% filter(grepl('Region/', Var2)) %>%
                                          mutate (Level = 'Region'), 
                                        fc.OLS.mint.shrink.wls.var %>% filter(Var2 %in% c('Purpose/Bus','Purpose/Hol', 
                                                                                          'Purpose/Oth','Purpose/Vis')) %>%
                                          mutate (Level = 'Purpose'), 
                                        fc.OLS.mint.shrink.wls.var %>% filter(grepl('State x Purpose/', Var2)) %>%
                                          mutate (Level = 'State x Purpose'), 
                                        fc.OLS.mint.shrink.wls.var %>% filter(grepl('Zone x Purpose/', Var2)) %>%
                                          mutate (Level = 'Zone x Purpose'), 
                                        fc.OLS.mint.shrink.wls.var %>% filter( !grepl('State', Var2) & !grepl('Zone', Var2) & 
                                                                                 !grepl('Region', Var2) & !grepl('Purpose', Var2) & 
                                                                                 !grepl('Total', Var2)) %>%
                                          mutate (Level = 'Region x Purpose'))

error.OLS.dif.rec <- dplyr::select ( fc.OLS.mint.shrink.wls.var, error , Level, Rec) %>% 
  mutate(Method = 'OLS')

error.tourism.struct <- error.tourism %>%
  select(-Rec) %>%
  mutate(Rec = 'wls_struct')

error.dif.rec <- bind_rows(
  dplyr::select ( ets.arima.dif.rec, error , Level, Rec, Method), error.OLS.dif.rec, error.tourism.struct) %>%
  mutate( facet = factor(Level,
                         levels = c("Total", "State", "Zone", "Region", "Purpose", "State x Purpose", "Zone x Purpose", "Region x Purpose")))
############################## Table 16
rmse.dif.rec <- error.dif.rec %>%
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
########################## Figure 21

p1 <- forecast.tourism.data %>%
  filter(Series == "AAAVis") %>%
  ggplot(aes(x = date, y = Actual, colour = "Actual", size = 'Actual')) +
  geom_ribbon(aes(x = date, ymax = ARIMA.upper.unrec, ymin = ARIMA.lower.unrec), fill = "gray83", colour = "gray83", alpha = .2, size = 0.5) +
  geom_ribbon(aes(x = date, ymax = ARIMA.upper.rec, ymin = ARIMA.lower.rec), fill = "lightblue", colour = "lightblue", alpha = .2, size = 0.5)  +
  geom_line(aes(y = ETS.rec, colour = "ARIMA.unrec", size = 'ARIMA.unrec')) +
  geom_line(aes(y = ARIMA.rec, colour = "ARIMA.rec", size = 'ARIMA.rec')) +
  geom_line() +
  xlab("Horizon") +
  ylab("Count") +
  ggtitle("Fixed origin multi-step forecasts") +
  scale_colour_manual("Method",
                      breaks = c("Actual", "ARIMA.unrec", "ARIMA.rec"),
                      values = c("black", "gray",  "blue"))+
  scale_size_manual(breaks = c("Actual", "ARIMA.unrec", 'ARIMA.rec'),
                    values = c( 0.8, 0.5,  0.5), guide = "none") +
  theme_bw()

p2  <- forecast.tourism.data %>%
  filter(Series == "AAAVis") %>%
  ggplot(aes(x = date, y = Actual, colour = "Actual", size = 'Actual')) +
  geom_ribbon(aes(x = date, ymax = ETS.upper.unrec, ymin = ETS.lower.unrec), fill = "gray83", colour = "gray83", alpha = .2, size = 0.5) +
  geom_ribbon(aes(x = date, ymax = ETS.upper.rec, ymin = ETS.lower.rec), fill = "lightgreen", colour = "lightgreen", alpha = .2, size = 0.5)  +
  geom_line(aes(y = ETS.rec, colour = "ETS.unrec", size = 'ETS.unrec')) +
  geom_line(aes(y = ETS.rec, colour = "ETS.rec", size= 'ETS.rec')) +
  geom_line() +
  xlab("Horizon") +
  ylab("Count") +
  ggtitle("") +
  scale_colour_manual("Method",
                      breaks = c("Actual", "ETS.unrec", 'ETS.rec'),
                      values = c("black", "gray",  "green"))+
  scale_size_manual(breaks = c("Actual", "ETS.unrec", 'ETS.rec'),
                    values = c( 0.8, 0.5,  0.5), guide = "none") +
  theme_bw()

p3 <- forecast.tourism.data %>%
  filter(Series == "AAAVis") %>%
  ggplot(aes(x = date, y = Actual, colour = "Actual", size = 'Actual')) +
  geom_ribbon(aes(x = date, ymax = OLS.upper.unrec, ymin = OLS.lower.unrec), fill = "gray83", colour = "gray83", alpha = .2, size = 0.5) +
  geom_ribbon(aes(x = date, ymax = OLS.upper.rec, ymin = OLS.lower.rec), fill = "pink", colour = "pink", alpha = .2, size = 0.5)  +
  geom_line(aes(y = ETS.rec, colour = "OLS.unrec", size = 'OLS.unrec')) +
  geom_line(aes(y = OLS.rec, colour = "OLS.rec", size = 'OLS.rec')) +
  geom_line() +
  xlab("Horizon") +
  ylab("Count") +
  ggtitle("") +
  scale_colour_manual("Method",
                      breaks = c("Actual", "OLS.unrec", 'OLS.rec'),
                      values = c("black", "gray",  "red"))+
  scale_size_manual(breaks = c("Actual", "OLS.unrec", 'OLS.rec'),
                    values = c( 0.8, 0.5,  0.5), guide = "none") +
  theme_bw()

cowplot::plot_grid(p1, p2, p3)


