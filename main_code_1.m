function [x,fval] = main_code_1()

data_c = load('x_c_data_values.txt');

num_dim = size(data_c,2);
num_variables = num_dim*2;

%%variables array contains n elements. The first n/2 elements are theta
%%values and the rest are P values

%calling the function
fun = @(variables) co_Krigging_method_smart_1(variables);

%boundaries
lb_theta = zeros(1,num_dim);
up_theta = 100.*ones(1,num_dim);


lb_P = 2.*ones(1,num_dim);
up_P = 2.*ones(1,num_dim);

lb_rho = 0;
up_rho = 100;

lb = [lb_theta lb_P lb_rho];
ub = [up_theta up_P up_rho];



%optimization
[x,fval] = ga(fun,num_variables+1,[],[],[],[],lb,ub);
end
