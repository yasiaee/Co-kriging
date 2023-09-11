clear
close all
clc

N = 10;
fileout1 = fopen('theta.txt','w');
fileout2 = fopen('p_vals.txt','w');
fileout3 = fopen('theta_d.txt','w');
fileout4 = fopen('p_d_vals.txt','w');
fileout5 = fopen('rho.txt','w');
fileout6 = fopen('rmse.txt','w');
for loopnum=1:N
    [x,fval] = main_code();
    n_dim = length(x)/2;

    theta_vals = x(1:n_dim);
    P_vals = x(n_dim+1:end);

    [x1,fval1] = main_code_1();
    n_dim1 = (length(x1) -1)/2;

    theta_d_vals = x1(1:n_dim1);
    %theta_d_vals =1.994785633563367e-04;
    P_d_vals = x1(n_dim1+1:end-1);
    rho = x1(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Loading data

   y_vals_c = load('y_data_values_c.txt');
   y_vals_e = load('y_data_values_e.txt');
   y_vals_c_e = load('y_data_values_c_e.txt');
   
   y_c = y_vals_c(:,1);
   y_e = y_vals_e(:,1);
   y_c_e = y_vals_c_e(:,1);

   y = [y_c;y_e];

   data_e = load('x_e_data_values.txt');
   n_rows_e = size(data_e,1);
   n_dim = size(data_e,2);

   data_c = load('x_c_data_values.txt');
   n_rows_c = size(data_c,1);
%n_dim = size(data_c,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   F_1 = ones(n_rows_c,1);
   F_2 = ones(n_rows_e,1);

%%
   d_c_c_x = zeros(n_rows_c,n_rows_c,n_dim);
   for i=1: n_dim
       d_c_c_x(:,:,i) = pdist2(data_c(:,i), data_c(:,i));
   end


   d_c_c_sum = zeros(n_rows_c,n_rows_c,n_dim);

   for i=1: n_dim
       d_c_c_sum(:,:,i) = theta_vals(i).*d_c_c_x(:,:,i).^P_vals(i);
      
   end
   d_c_c_sum_sum = sum(d_c_c_sum,3);
   psi_c_c = exp(-d_c_c_sum_sum);
%%
   d_e_e_x = zeros(n_rows_e,n_rows_e,n_dim);

for i=1: n_dim
    d_e_e_x(:,:,i) = pdist2(data_e(:,i), data_e(:,i));
end

d_e_e_sum = zeros(n_rows_e,n_rows_e,n_dim);

for i=1: n_dim
    d_e_e_sum(:,:,i) = theta_vals(i).*d_e_e_x(:,:,i).^P_vals(i);
end

d_e_e_sum_sum = sum(d_e_e_sum,3);
psi_e_e = exp(-d_e_e_sum_sum);
%%
d_c_e_x = zeros(n_rows_c,n_rows_e,n_dim);

for i=1: n_dim
    d_c_e_x(:,:,i) = pdist2(data_c(:,i), data_e(:,i));
end

d_c_e_sum = zeros(n_rows_c,n_rows_e,n_dim);

for i=1: n_dim
    d_c_e_sum(:,:,i) = theta_vals(i).*d_c_e_x(:,:,i).^P_vals(i);
end

d_c_e_sum_sum = sum(d_c_e_sum,3);
psi_c_e = exp(-d_c_e_sum_sum);

%%
d_e_c_x = zeros(n_rows_e,n_rows_c,n_dim);

for i=1: n_dim
    d_e_c_x(:,:,i) = pdist2(data_e(:,i), data_c(:,i));
end

d_e_c_sum = zeros(n_rows_e,n_rows_c,n_dim);

for i=1: n_dim
    d_e_c_sum(:,:,i) = theta_vals(i).*d_e_c_x(:,:,i).^P_vals(i);
end

d_e_c_sum_sum = sum(d_e_c_sum,3);
psi_e_c = exp(-d_e_c_sum_sum);

%%
d_e_e_x_d = zeros(n_rows_e,n_rows_e,n_dim);

for i=1: n_dim
    d_e_e_x_d(:,:,i) = pdist2(data_e(:,i), data_e(:,i));
end

d_e_e_sum_d = zeros(n_rows_e,n_rows_e,n_dim);

for i=1: n_dim
    d_e_e_sum_d(:,:,i) = theta_d_vals(i).*d_e_e_x_d(:,:,i).^P_d_vals(i);
end

d_e_e_sum_sum_d = sum(d_e_e_sum_d,3);
psi_e_e_d = exp(-d_e_e_sum_sum_d);
%%
mu_c = (transpose(F_1)*inv(psi_c_c)*y_c)/(transpose(F_1)*inv(psi_c_c)*F_1);
sigma_c_2 = transpose(y_c - F_1*mu_c)*inv(psi_c_c)*(y_c -F_1*mu_c)/8;

d = y_e - rho*y_c_e;

mu_d = (transpose(F_2)*inv(psi_e_e_d)*d)/(transpose(F_2)*inv(psi_e_e_d)*F_2);
sigma_d_2 = transpose(d - F_2*mu_d)*inv(psi_e_e_d)*(d -F_2*mu_d)/8;

A = sigma_c_2*psi_c_c;
B = rho*sigma_c_2*psi_c_e;
D = rho*sigma_c_2*psi_e_c;
E = (rho^2*sigma_c_2*psi_e_e + sigma_d_2*psi_e_e_d);

C = [A B;D E];
%%


 

data_new = load('x_test.txt');
n_rows_new = size(data_new,1);
n_dim = size(data_new,2);






for i=1:n_rows_new
  

    % d_c_new_x1 = pdist2(data_c_x1(:), x_new_x1(:));
    % d_c_new_x2 = pdist2(data_c_x2(:), x_new_x2(:));
    % d_c_new_sum = -((-theta1.*d_c_new_x1.^p1) + (-theta2.*d_c_new_x2.^p2));
    % psi_c_new = exp(d_c_new_sum);

    d_c_new_x = zeros(n_rows_c,n_dim);
    for j=1: n_dim
        d_c_new_x(:,j) = pdist2(data_c(:,j), data_new(i,j));
    end
    d_c_new_sum = zeros(n_rows_c,n_dim);
    for j=1: n_dim
        d_c_new_sum(:,j) = theta_vals(j).*d_c_new_x(:,j).^P_vals(j);
    end

    d_c_new_sum_sum = sum(d_c_new_sum,2);
    psi_c_new = exp(-d_c_new_sum_sum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check after this

    % d_e_new_x1 = pdist2(data_e_x1(:), x_new_x1(:));
    % d_e_new_x2 = pdist2(data_e_x2(:), x_new_x2(:));
    % d_e_new_sum_c = -((-theta1.*d_e_new_x1.^p1) + (-theta2.*d_e_new_x2.^p2));
    % psi_e_new_c = exp(d_e_new_sum_c);

    d_e_new_x = zeros(n_rows_e,n_dim);
    for j=1: n_dim
        d_e_new_x(:,j) = pdist2(data_e(:,j), data_new(i,j));
    end
    d_e_new_sum = zeros(n_rows_e,n_dim);
    for j=1: n_dim
        d_e_new_sum(:,j) = theta_vals(j).*d_e_new_x(:,j).^P_vals(j);
    end

    d_e_new_sum_sum = sum(d_e_new_sum,2);
    psi_e_new = exp(-d_e_new_sum_sum);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check!



    % d_e_new_x1 = pdist2(data_e_x1(:), x_new_x1(:));
    % d_e_new_x2 = pdist2(data_e_x2(:), x_new_x2(:));
    % d_e_new_sum_d = -((-theta_d1.*d_e_new_x1.^p_d1) + (-theta_d2.*d_e_new_x2.^p_d2));
    % psi_e_new_d = exp(d_e_new_sum_d);


    d_e_new_x_d = zeros(n_rows_e,n_dim);
    for j=1: n_dim
        d_e_new_x_d(:,j) = pdist2(data_e(:,j), data_new(i,j));
    end
    d_e_new_sum_d = zeros(n_rows_e,n_dim);
    for j=1: n_dim
        d_e_new_sum_d(:,j) = theta_d_vals(j).*d_e_new_x(:,j).^P_d_vals(j);
    end

    d_e_new_sum_sum_d = sum(d_e_new_sum_d,2);
    psi_e_new_d = exp(-d_e_new_sum_sum_d);




    a = rho*sigma_c_2*psi_c_new;
    b = rho^2*sigma_c_2*psi_e_new + sigma_d_2*psi_e_new_d;

    c = [a;b];

    mu = (transpose(ones(n_rows_e+n_rows_c,1))*inv(C)*y)/(transpose(ones(n_rows_e+n_rows_c,1))*inv(C)*ones(n_rows_e+n_rows_c,1));

    y_predict(i) = mu + transpose(c)*inv(C)*(y - ones(n_rows_e+n_rows_c,1)*mu);
end

y_true_data = load('y_true.txt');
y_true = y_true_data';



for lpp = 1:length(theta_vals)
    fprintf(fileout1,'%d',theta_vals);
    fprintf(fileout2,'%d',P_vals);
    fprintf(fileout3,'%d',theta_d_vals);
    fprintf(fileout4,'%d',P_d_vals);
end
fprintf(fileout1,'\n');
fprintf(fileout2,'\n');
fprintf(fileout3,'\n');
fprintf(fileout4,'\n');

fprintf(fileout5,'%d\n',rho);

rmse = sqrt(sum((y_true-y_predict).^2)/length(y_predict));
fprintf(fileout6,'%d\n',rmse);
end

figure(1)
plot(data_new,y_true)
hold on
plot(data_new,y_predict)
legend('Ground Truth','Prediction')
title("3 high")