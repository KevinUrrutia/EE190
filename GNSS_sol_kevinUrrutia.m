clear; clc; close all;

y = 1e7 .* [1.6019005943;
            1.6169407378;
            1.6461726771;
            1.7403040927;
            1.8006539877;
            1.8665085004]; % 6x1

P = 1e7 .* [ 0.1743114855,  0           , 1.9923893962;
             0.2588190451,  0.4482877361, 1.9318516526;
            -0.4226182617,  0.7319963015, 1.8126155741;
            -1.4142135624,  0           , 1.4142135624;
            -0.8191520443, -1.4188129598, 1.1471528727;
             0.9063077870, -1.5697711344, 0.8452365235]; %6x3

N = 10;

x = zeros(3, N); %3xN
f = zeros(size(y, 1), N); %6xN
F = zeros(6, 3); %6x3

i = 1;
while(true)
    for j = 1 : size(y, 1)
        v = (x(:, i) - P(j, :)');

        %calculate the range
        f(j, i) = (v' * v)^(1/2);

        %calculate the Jacobian
        F(j, :) = v / f(j, i)';
    end

    %calculate the new estimate of x
    x(:, i + 1) = x(:, i) + inv(F' * F) * (F' * (y - f(:, i))) ;

    %if the difference in position estimate is less than 1cm break loop
    delta_x = x(:, 1 +i) - x(:, i);
    if(((delta_x' * delta_x)^(1/2)) < 0.1)
        break;
    end

    %increment the column counter
    i = i + 1;
end