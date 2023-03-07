clear; clc; close all;

mu = 0;
sigma = sqrt(2);
b = 4;

N = 100000;

D1 = clt(1, N, mu, sigma);  
D2 = clt(2, N, mu, sigma); 
D3 = clt(3, N, mu, sigma); 
D4 = clt(4, N, mu, sigma); 
D8 = clt(8, N, mu, sigma); 
D16 = clt(16, N, mu, sigma); 

x = [-5:.1:5];
y = normpdf(x,0,1); 

subplot(3,2,1)
histogram(D1, 'Normalization', 'pdf');
hold on;
plot(x, y)
title("n = 1")
xlabel("Standard Devivations")
ylabel("Frequency")

subplot(3,2,2)
histogram(D2, 'Normalization', 'pdf');
hold on;
plot(x, y)
title("n = 2")
xlabel("Standard Devivations")
ylabel("Frequency")

subplot(3,2,3)
histogram(D3, 'Normalization', 'pdf');
hold on;
plot(x, y)
title("n = 3")
xlabel("Standard Devivations")
ylabel("Frequency")

subplot(3,2,4)
histogram(D4, 'Normalization', 'pdf');
hold on;
plot(x, y)
title("n = 4")
xlabel("Standard Devivations")
ylabel("Frequency")

subplot(3,2,5)
histogram(D8, 'Normalization', 'pdf');
hold on;
plot(x, y)
title("n = 8")
xlabel("Standard Devivations")
ylabel("Frequency")

subplot(3,2,6)
histogram(D16, 'Normalization', 'pdf');
hold on;
plot(x, y)
title("n = 16")
xlabel("Standard Devivations")
ylabel("Frequency")

% function z = clt(n, N, mu, sigma)
%     z = zeros(N, 1);
%     for i = 1:N
%         x = repmat(mu, n, 1) + randl(n, 1)*sigma;
%         z(i, n) = (sum(x) - n*mu);
% %         z(i, n) =  (sum(x) - n*mu) / (sigma * sqrt(n));
%     end
% end

function z = clt(n, N, mu, sigma)
    x = repmat(mu, n, N) + randl(n, N)*sigma;
    Sn = sum(x, 1);

    z = (Sn - n*mu) / (sigma * sqrt(n));
end