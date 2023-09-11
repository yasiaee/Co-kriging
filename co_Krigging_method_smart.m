function [objective_fun] = co_Krigging_method_smart(variables)
% variable1 is an array that contains theta values
% variable2 is an array that contains P values

y_vals_c = load('y_data_values_c.txt');
y_vals_e = load('y_data_values_e.txt');
y_vals_c_e = load('y_data_values_c_e.txt');

y_c = y_vals_c(:,1);
y_e = y_vals_e(:,1);
y_c_e = y_vals_c_e(:,1);


data_c = load('x_c_data_values.txt');
n_dim = size(data_c,2);
n_rows = size(data_c,1);

F = ones(n_rows,1);

theta_vals = variables(1:n_dim);
P_vals = variables(n_dim+1:end);


d_c_c_x = zeros(n_rows,n_rows,n_dim);

for i=1: n_dim
    d_c_c_x(:,:,i) = pdist2(data_c(:,i), data_c(:,i));
end

d_c_c_sum = zeros(n_rows,n_rows,n_dim);

for i=1: n_dim
    d_c_c_sum(:,:,i) = theta_vals(i).*d_c_c_x(:,:,i).^P_vals(i);
end
d_c_c_sum_sum = sum(d_c_c_sum,3);
psi_c_c = exp(-d_c_c_sum_sum);

% d_c_e = pdist2(x_c(:), x_e(:));
% psi_e_c = exp(-theta.*d_c_e.^p);
% 
% d_e_e = pdist2(x_e(:), x_e(:));
% psi_e_e = exp(-theta.*d_e_e.^p);

mu_c = (transpose(F)*inv(psi_c_c)*y_c)/(transpose(F)*inv(psi_c_c)*F);
sigma_c_2 = transpose(y_c - F*mu_c)*inv(psi_c_c)*(y_c -F*mu_c)/4;

objective_fun = -((-n_rows/2)*log(sigma_c_2) - 0.5*log(det(psi_c_c)));

end

% 
% 
% mu_first_term = (transpose(F)*inv(psi)*F);
% mu_second_term = (transpose(F)*inv(psi)*y);
% mu = inv(mu_first_term)*mu_second_term;
% 
% sigma = 0.25*( transpose(y - F*mu)*inv(psi)*(y - F*mu));
% 
% 
% objective_func = -(-0.5*(4*log(sigma) + log(det(psi))));
% 
% 
% end