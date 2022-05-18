## Math6120 - Nonlinear Optimisation
## Coursework 1 - AMPL Model
## Emma Tarmey, 2940 4045


## This file specifies the mathematical model used to solve our problem


# parameters
param N; # number of subintervals to be created
param theta_max;
param height_0;
param velocity_0;
param D_0 := (velocity_0 / 2);


# sets
set TIMEINDEX  := {0..N};     # moments in time
set KTIMEINDEX := {1..(N-1)}; # moments in time except for time 0 and time N


# variables

# decision variables
var t_final >= 0;
var height   {t in TIMEINDEX} >= 0;
var velocity {t in TIMEINDEX} >= 0;
var mass     {t in TIMEINDEX} >= 0;
var thrust   {t in TIMEINDEX} >= 0;

# define time passing in the interval [t_0, t_final], points in time are indexed by the set TIMEINDEX
var REALTIME {t in TIMEINDEX}          = ((t / N) * t_final);

# defines the starting points of the velocity and mass functions
# these starting points are used via the 'let' command in the corresponding 'cw1.run' file
var velocity_function {t in TIMEINDEX} = ((REALTIME[t] / t_final) * (1 - (REALTIME[t] / t_final)));
var mass_function {t in TIMEINDEX}     = (1 + ((REALTIME[t] / t_final) * (0.6 - 1)));

# D(height, velocity) function from spec
var D {t in TIMEINDEX}   = D_0 * (velocity[t] * velocity[t]) * exp(-1 * height_0 * ( (height[t] - height[0]) / height[0] ));

# Right-hand-side of below constraints, seperated out for readability
var RHS {t in TIMEINDEX} = (((thrust[t] - D[t]) / mass[t]) - ((height[0] / height[t]) * (height[0] / height[t])));


# objective function
maximize final_height: height[ N ] ;


# equality constraints
subject to velocity_initial:
	velocity[0] = 0;

subject to height_initial:
	height[0]   = 1;

subject to mass_initial:
	mass[0]     = 1;

subject to velocity_final:
	mass[N] = 0.6;


# lower and upper bound inequality constraints
subject to velocity_lower_bound {t in TIMEINDEX}:
	velocity[t] >= 0;

subject to height_lower_bound {t in TIMEINDEX}:
	height[t]   >= height[0];

subject to mass_lower_bound {t in TIMEINDEX}:
	mass[t]     >= mass[N]; # final mass as lower bound

subject to mass_upper_bound {t in TIMEINDEX}:
	mass[t]     <= mass[0]; # initial mass as upper bound

subject to thrust_upper_bound {t in TIMEINDEX}:
	thrust[t]   <= theta_max;

subject to thrust_lower_bound {t in TIMEINDEX}:
	thrust[t]   >= 0;


# rocket motion height constraints
subject to height_change_initial:
	N * (height[1] - height[0])     =     (t_final * velocity[0]);

subject to height_change {k in KTIMEINDEX}:
	N * (height[k+1] - height[k-1]) = (2 * t_final * velocity[k]);

subject to height_change_final:
	N * (height[N] - height[N-1])   =     (t_final * velocity[N]);


# rocket motion velocity constraints
subject to velocity_change_initial:
	N * (velocity[1] - velocity[0])     =     (t_final * RHS[0]);

subject to velocity_change {k in KTIMEINDEX}:
	N * (velocity[k+1] - velocity[k-1]) = (2 * t_final * RHS[k]);

subject to velocity_change_final:
	N * (velocity[N] - velocity[N-1])   =     (t_final * RHS[N]);


# rocket motion mass constraints
subject to mass_change_initial:
	N * (mass[1] - mass[0])     = ((-2) * t_final * thrust[0]);

subject to mass_change {k in KTIMEINDEX}:
	N * (mass[k+1] - mass[k-1]) = ((-4) * t_final * thrust[k]);

subject to mass_change_final:
	N * (mass[N] - mass[N-1])   = ((-2) * t_final * thrust[N]);

