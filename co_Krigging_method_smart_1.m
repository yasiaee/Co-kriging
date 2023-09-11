function [objective_fun_1] = co_Krigging_method_smart_1(variables)
% variable1 is an array that contains theta values
% variable2 is an array that contains P values

y_vals_c = load('y_data_values_c.txt');
y_vals_e = load('y_data_values_e.txt');
y_vals_c_e = load('y_data_values_c_e.txt');

y_c = y_vals_c(:,1);
y_e = y_vals_e(:,1);
y_c_e = y_vals_c_e(:,1);

data_c = load('x_c_data_values.txt');
data_e = load('x_e_data_values.txt');
n_dim = size(data_e,2);
n_rows = size(data_e,1);

F = ones(n_rows,1);

theta_d_vals = variables(1:n_dim);
P_d_vals = variables(n_dim+1:end-1);
rho = variables(end);


d_e_e_x = zeros(n_rows,n_rows,n_dim);

for i=1: n_dim
    d_e_e_x(:,:,i) = pdist2(data_e(:,i), data_e(:,i));
end

d_e_e_sum = zeros(n_rows,n_rows,n_dim);

for i=1: n_dim
    d_e_e_sum(:,:,i) = theta_d_vals(i).*d_e_e_x(:,:,i).^P_d_vals(i);
end

d_e_e_sum_sum = sum(d_e_e_sum,3);
psi_e_e = exp(-d_e_e_sum_sum);

d = y_e - rho*y_c_e;

mu_d = (transpose(F)*inv(psi_e_e)*d)/(transpose(F)*inv(psi_e_e)*F);
sigma_d_2 = transpose(d - F*mu_d)*inv(psi_e_e)*(d -F*mu_d)/4;

objective_fun_1 = -((-n_rows/2)*log(sigma_d_2) - 0.5*log((det(psi_e_e))));


end