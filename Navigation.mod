param bigM;

param horizon;

param nDims;

param minmazebounds_d{d in 1..nDims};

param maxmazebounds_d{d in 1..nDims};

param minactionbounds_d{d in 1..nDims};

param maxactionbounds_d{d in 1..nDims};

param minmudbounds_d{d in 1..nDims};

param maxmudbounds_d{d in 1..nDims};

param initial_d{d in 1..nDims};

param goal_d{d in 1..nDims};

param discount;

param nCond;

var move_dt{d in 1..nDims, t in 1..horizon};

var location_dt{d in 1..nDims, t in 1..horizon+1};

var locationPlus_dt{d in 1..nDims, t in 1..horizon+1} >=0;

var locationMinus_dt{d in 1..nDims, t in 1..horizon+1} >=0;

var if_cdt{c in 1..nCond, d in 1..nDims, t in 1..horizon} binary;

var then_t{t in 1..horizon} binary;

minimize z: sum{d in 1..nDims, t in 2..horizon+1} (locationPlus_dt[d,t] + locationMinus_dt[d,t]);

subject to initialCon{d in 1..nDims}:
	location_dt[d,1] = initial_d[d];

subject to AuxGoalCon{d in 1..nDims, t in 2..horizon+1}:
	location_dt[d,t] - locationPlus_dt[d,t] + locationMinus_dt[d,t] = goal_d[d];

subject to ifCon1{d in 1..nDims, t in 1..horizon}:
	location_dt[d,t] - bigM*if_cdt[1,d,t] <= minmudbounds_d[d];

subject to ifCon2{d in 1..nDims, t in 1..horizon}:
	location_dt[d,t] - bigM*if_cdt[1,d,t] >= minmudbounds_d[d] - bigM;

subject to ifCon3{d in 1..nDims, t in 1..horizon}:
	location_dt[d,t] + bigM*if_cdt[2,d,t] >= maxmudbounds_d[d];

subject to ifCon4{d in 1..nDims, t in 1..horizon}:
	location_dt[d,t] + bigM*if_cdt[2,d,t] <= maxmudbounds_d[d] + bigM;

subject to ifThenCon1{t in 1..horizon}:
	sum{c in 1..nCond, d in 1..nDims} (if_cdt[c,d,t]) - nCond*nDims - then_t[t] <= -1;

subject to ifThenCon2{c in 1..nCond, d in 1..nDims, t in 1..horizon}:
	then_t[t] - if_cdt[c,d,t]<= 0;

subject to minMoveCon{d in 1..nDims, t in 1..horizon}:
	move_dt[d,t] >= minactionbounds_d[d];

subject to maxMoveCon{d in 1..nDims, t in 1..horizon}:
	move_dt[d,t] <= maxactionbounds_d[d];

subject to minLocationCon{d in 1..nDims, t in 1..horizon+1}:
	location_dt[d,t] >= minmazebounds_d[d];

subject to maxLocationCon{d in 1..nDims, t in 1..horizon+1}:
	location_dt[d,t] <= maxmazebounds_d[d];

subject to nextLocationCon1{d in 1..nDims, t in 1..horizon}:
	location_dt[d,t+1] - location_dt[d,t] - discount*move_dt[d,t] + bigM*then_t[t] <= bigM;

subject to nextLocationCon2{d in 1..nDims, t in 1..horizon}:
	location_dt[d,t+1] - location_dt[d,t] - discount*move_dt[d,t] - bigM*then_t[t] >= -1*bigM;

subject to nextLocationCon3{d in 1..nDims, t in 1..horizon}:
	location_dt[d,t+1] - location_dt[d,t] - move_dt[d,t] - bigM*then_t[t] <= 0;

subject to nextLocationCon4{d in 1..nDims, t in 1..horizon}:
	location_dt[d,t+1] - location_dt[d,t] - move_dt[d,t] + bigM*then_t[t] >= 0;