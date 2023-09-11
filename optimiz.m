clear 
clc
close all

% defining the bounds
lb = [0, 2,0,2]; %limitation of theta and p
ub = [100, 2,100,2];

%calling the function
fun = @(variables) co_Krigging_method(variables(1), variables(2),variables(3), variables(4));

%optimization
[x,fval] = ga(fun,4,[],[],[],[],lb,ub)


% x_test = linspace(0,1,50);
% 
% for i =1:50
%     y_pred(i) = Krigging_predict(x(1),x(2),x_test(i));
% end
% 
% m = linspace(0,1,50);
% y_true = (6.*m - 2).^2 .* sin(2.*(6.*m-2));
% 
% 
% figure(1)
% plot(x_test,y_pred,'r--')
% hold on
% plot(m,y_true,'k-')
% legend('Krigging','True')
% xlabel('x')
% ylabel('y')
