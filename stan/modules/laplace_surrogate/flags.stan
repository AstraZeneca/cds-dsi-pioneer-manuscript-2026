// Runtime switch for marginalizing backgrounded-trial patients via the
// log-concave quadratic surrogate. 0 = off (backgrounded patients contribute
// nothing to the LFO likelihood, i.e. current behavior); 1 = on.
int<lower=0, upper=1> enable_background_surrogate;
