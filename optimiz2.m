clear 
clc
close all

%theta_d1,p_d1,theta_d2,p_d2,rho
% defining the bounds
lb = [0, 2,0,2,0]; %limitation of theta and p
ub = [100, 2,100,2,100];

%calling the function
fun = @(variables) co_Krigging_method_1(variables(1), variables(2), variables(3),variables(4), variables(5));

%optimization
[x,fval] = ga(fun,5,[],[],[],[],lb,ub)