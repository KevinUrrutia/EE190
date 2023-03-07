clear; clc; close all;

% x = 1: 100;
% by = binocdf(x, 100, .4);
% y = normcdf(x,40,sqrt(24));
% 
% subplot(1,2,1)
% plot(x, by, '.');
% hold on;
% plot(x, y, '.');
% title("Binomial and Normal PDF Curves")
% xlabel("Trials")
% ylabel("Frequency")
% legend("Binomial", "Normal");
% 
% subplot(1,2,2)
% Y = abs(by - y);
% plot(x, Y, '.');
% title("Binomial and Normal Subtraction")
% xlabel("Trials")
% ylabel("Frequency")

N = 100;
p = .40;
mu = 40; 
sigma = sqrt(24);
k1 = 9;
k2 = k1 : N; 

norm_approx = normcdf(k2, mu, sigma) - normcdf(k1, mu, sigma);
de_moive_approx = normcdf(k2, -0.5 + mu, sigma) - normcdf(k1, 0.5 + mu, sigma);
binom = binocdf(k2,N, p) - binocdf(k1, N, p);
subplot(2,1,1);
plot(k2, norm_approx, ".")
hold on
plot(k2, de_moive_approx, ".")
hold on 
plot(k2, binom, ".")
xlabel("k_2");
ylabel("P(k_2)");
legend("normal approx", "de moive approx", "binomial")

subplot(2,1,2);
norm_diff = abs(binom - norm_approx);
de_moive_diff = abs(binom - de_moive_approx);
plot(k2, norm_diff, ".");
hold on
plot(k2, de_moive_diff, ".");
xlabel("k_2");
ylabel("P(k_2)");
legend("Normal difference", "De moicve difference");