library(fable)
lung_deaths_agg <- as_tsibble(cbind(mdeaths, fdeaths)) %>%
    aggregate_key(key, value = sum(value))
  
result <- lung_deaths_agg %>%
    model(lm = TSLM(log(value + 1) ~ trend() + season())) %>%
    reconcile(lm_adjusted = min_trace(lm)) %>%
    forecast()
  
tail(result)
