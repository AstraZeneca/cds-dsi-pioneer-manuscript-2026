// Fixed calendar anchors (weeks) for the surrogate quadratic. The first anchor
// MUST be 0 (baseline normalization pins g(0)=0; the surrogate intercept relies
// on it). Tune the other two to the backgrounded trials' observation window.
vector[3] surrogate_anchor_times;
