% clear; clc; close all;
% 
% theta = -30:0.5:40;
% 
% yk = [6.13, 4.88, 5.89, 6.29, 5.78, 5.04];
% 
% sigma = 1;
% ml_a = -(1/2)*log(2*pi) - (1/2)*(yk(1) - theta).^2;
% 
% 
% sum_int = zeros(2, size(theta, 2));
% for i = 1:2
%     sum_int(i, :) = (yk(i) - theta).^2;
% end
% 
% ml_b = -(2/2)*log(2*pi) - (1/2)*sum(sum_int, 1);
% 
% 
% sum_int = zeros(6, size(theta, 2));
% for i = 1:6
%     sum_int(i, :) = (yk(i) - theta).^2;
% end
% 
% ml_c = -(6/2)*log(2*pi) - (1/2)*sum(sum_int, 1);
% 
% 
% [max_liklihood, I] = max(ml_a)
% [max_liklihood, I] = max(ml_b)
% [max_liklihood, I] = max(ml_c)
% 
% subplot(3,1,1);
% plot(theta, ml_a, ".");
% title("N=1")
% ylabel("Liklihood")
% xlabel("theta")
% 
% 
% subplot(3,1,2);
% plot(theta, ml_b, ".");
% title("N=2")
% ylabel("Liklihood")
% xlabel("theta")
% 
% subplot(3,1,3);
% plot(theta, ml_c, ".");
% title("N=6")
% ylabel("Liklihood")
% xlabel("theta")


clear; clc; close all;
theta = 0:0.01:1;

N = 1000000;
n = 10;
k = N.*(theta);
X = zeros(1, size(k,2));
Y = hygepdf(X,N,k,n);

AOQ = theta.*Y;

subplot(2,1,1)
plot(theta, Y, '.')
title("OC")
xlabel("Theta")
ylabel("OC(theta; N,n,c)")

subplot(2,1,2)
plot(theta, AOQ, '.')
title("AOQ")
xlabel("Theta")
ylabel("AOQ(theta; N,n,c)")

