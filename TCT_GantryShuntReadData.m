function [] = TCT_GantryShuntReadData()

%% Read in the data file
if isfile('GantryShunt.mat')
    % eg is the measurement of the gantry position.
    % ds is the measurement of the shunt machine position.
    % tm is the time vector corresponding to the above measurements
    load('GantryShunt.mat');
    nr = length(tm);
else
    error('Error: File GantryShunt.mat does not exist.');
end

%% Initialize variables
% student code here

%% Preallocate space
xsh         = zeros(3,nr);          % preallocate shunt state vector
xgh(1,:)    = ones(1,nr)*eg(1,1);   % preallocate gantry state vector with initialization
xgh(2:3,:)  = zeros(2,nr);          % preallocate gantry state vector 
dt          = diff(tm);             % time step between adjacent measurements (not constant)
phi         = zeros(3, 3, nr);
Qsh         = zeros(3,3,nr);
Qgh         = zeros(3,3,nr);
gamma       = zeros(3, nr);

%back_difference the velocity estimate 
xgh(2,2:end) = diff(eg) ./ dt; % gantry [m/s]
xsh(2,2:end) = diff(eg) ./ dt; %shunt [m/s]

%back difference the acceleration estimate
xgh(3,2:end) = diff(xgh(2,:)) ./ dt'; %gantry [m/s^2]
xsh(3,2:end) = diff(xsh(3,:)) ./ dt'; %shunt [m/s^2]

F = [0 1 0;
     0 0 1;
     0 0 0];

H = [1, 0, 0]; 

Rgh = var(eg);
Rsh = var(ds);
%phi is variable at every dt since dt is not constant
for i = 1:nr-1
    phi(:,:,i) = eye(3) + (F * dt(i)) + (F * dt(i))^2/factorial(2);
    gamma(:, i) = [0.5^2; dt(i); 1];
end

% p = polyfit(tm', xgh(1,:), 2);
% xgh(1,:) = polyval(p,tm');
% xgh(2,:) = polyval(polyder(p),tm');
% xgh(3,:) = polyval(polyder(polyder(p)), tm');
% 
% p = polyfit(tm', xsh(1,:), 2);
% xsh(1,:) = polyval(p,tm');
% xsh(2,:) = polyval(polyder(p),tm');
% xsh(3,:) = polyval(polyder(polyder(p)), tm');

%smooth position measurement to fit a second order polynomial 
for k = 4:nr - 3
    %y = [y(k* - 1), y(k*), y(k* + 1)]
    y_gantry = [eg(k - 3),eg(k - 2),eg(k - 1), eg(k), eg(k + 1), eg(k + 2),eg(k + 3)]';
    y_shunt = [ds(k - 3),ds(k - 2),ds(k - 1), ds(k), ds(k + 1),ds(k + 2),ds(k + 3)]';

    A = [1, (k - 3), 0.5*(k - 3)^2;
         1, (k - 2), 0.5*(k - 2)^2;
         1, (k - 1), 0.5*(k - 1)^2;
         1,  0,       0;
         1, (k + 1), 0.5*(k + 1)^2;
         1, (k + 2), 0.5*(k + 2)^2;
         1, (k + 3), 0.5*(k + 3)^2;];

    C_gantry = (A'*A)\A'*y_gantry;
    C_shunt = (A'*A)\A'*y_shunt;

 
    %velocity
    xgh(2,k) = C_gantry(2) + 2*C_gantry(3)*k; %gantry
    xsh(2,k) = C_shunt(2) + 2*C_shunt(3)*k; %shunt


    %acceleration
    xgh(3,k) = 2*C_gantry(3); %gantry 
    xsh(3,k) = 2*C_shunt(3); %shunt
end
omegagh = diff(xgh(3,:)) ./ dt';
omegash = diff(xsh(3,:)) ./ dt';

for i = 1:nr-1
    Qgh(:,:,i) = gamma(:,i) * gamma(:, i)' * var(omegagh);
    Qsh(:,:,i) = gamma(:,i) * gamma(:, i)' * var(omegash);
end

%% Process measurements
for i = 1:nr - 1
    if i == 1            % initialize the state and covariance for estimators
        % student code here
        
        Pgh = blkdiag(var(eg),var(xgh(2,:)),var(omegagh));
        Psh = blkdiag(var(ds),var(xsh(2,:)),var(omegash));
 
    else
      
        % student code here    % time propagation over each dt(i)
        xsh(:,i) = phi(:,:,i) * xgh(:, i) + gamma(i)*omegagh(i);
        Pgh_next = phi(:,:,i) * Pgh * phi(:,:,i)' + Qgh(:,:,i);

        xsh(:, i) = phi(:,:,i) * xsh(:, i) + gamma(i)*omegash(i);
        Psh_next = phi(:,:,i) * Psh * phi(:,:,i)' + Qsh(:,:,i);

        % student code here    % measurement update
        %calculate y hat
        y_hat_gh = H * xgh(:, i);
        y_hat_sh = H * xsh(:, i);

        %calcualte the residual
        rgh = eg(i) - y_hat_gh;
        rsh = ds(i) - y_hat_sh;

        %calculate the Kalman gain
        Kgh = Pgh * H' * inv(H * Pgh * H' + Rgh);
        xgh(:, i) = xgh(:, i) + Kgh * rgh;
        Pgh = Pgh - Kgh*(H*Pgh);

        Ksh = Psh * H' * inv(H * Psh * H' + Rsh);
        xsh(:, i) = xsh(:, i) + Ksh * rsh;
        Psh = Psh - Ksh*(H*Psh);

    end
end
pos_std_gh = ones(1, size(tm, 1)) * sqrt(Rgh);
pos_std_sh = ones(1, size(tm, 1)) * sqrt(Rsh);

neg_std_gh = ones(1, size(tm, 1)) * (-sqrt(Rgh));
neg_std_sh = ones(1, size(tm, 1)) * (-sqrt(Rsh));

na = 0;   % number of axes  
figure(1)
clf
na = na+1;
ax(na) = subplot(411);
plot(tm,xgh(1,:)','.', tm, eg,'.');grid on
ylabel('Gantry pos, m')
na = na+1;
ax(na) = subplot(412);
plot(tm,xgh(2,:),'.');grid on
ylabel('Gantry vel, m/s')
xlabel('Time, t, sec')
na = na+1;
ax(na) = subplot(413);
plot(tm,xgh(3,:)','.');grid on
ylabel('Gantry accel, m/s/s')
na = na+1;
ax(na) = subplot(414);
resg = eg - xgh(1,:)';
plot(tm,resg,'.');
hold on
plot(tm, pos_std_gh,'r.');
hold on
plot(tm, neg_std_gh,'r.');
grid on
ylabel('Gantry res, m')
xlabel('Time, t, sec')

figure(2)
clf
na = na+1;
ax(na) = subplot(411);
plot(tm,xsh(1,:)','.',tm, ds,'.');grid on
ylabel('Shunt pos, m')
na = na+1;
ax(na) = subplot(412);
plot(tm,xsh(2,:)','.');grid on
ylabel('Shunt vel, m/s')
xlabel('Time, t, sec')
na = na+1;
ax(na) = subplot(413);
plot(tm,xsh(3,:)','.');grid on
ylabel('Shunt accel, m/s/s')
na = na+1;
ax(na) = subplot(414);
ress  = ds - xsh(1,:)';
plot(tm,ress,'.');
hold on
plot(tm, pos_std_sh,'r.');
hold on
plot(tm, neg_std_sh,'r.');
grid on 
ylabel('Shunt res, m')
xlabel('Time, t, sec')

linkaxes(ax,'x')
figure(3)
maxdt = max(dt)
mindt = min(dt)
hist(dt); grid on
xlabel('Time step, dt, sec.')
ylabel('Histogram, counts')