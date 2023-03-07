x = [1.77,4.85,2.54,2.63,3.71,4.36,3.91,4.34,5.42];
y = [6.90,-0.18,3.86,3.80,3.09,2.65,2.17,2.51,-0.50];

mu_x = mean(x);
mu_y = mean(y);

Pxy = (1/(length(x)))*sum((y-mu_y).*(x -mu_x));
Pxx = (1/(length(x)))*sum((x -mu_x));
Pyy = (1/(length(x)))*sum((y -mu_y));

mu_xy = mu_x + Pxy * Pyy * (y - mu_y);

P_xy = Pxx - Pxy * (1/Pyy) * Pxy;