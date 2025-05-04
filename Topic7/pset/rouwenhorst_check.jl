
rho = 0.9;
sig = 10;
mu = -(1-rho)*0.5*sig^2
mc = QuantEcon.rouwenhorst(3,rho,sqrt(1-rho^2)*sig,mu)

eg = exp.(mc.state_values)
ss_e = stationary_distributions(mc)[1]

sum(yg.*ss_e)





