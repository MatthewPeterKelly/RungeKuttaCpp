% MAIN_readDataLog
%
%


DATA = csvread('logFile.csv');
t = DATA(:,1);
q = DATA(:,2);
w = DATA(:,3);


figure(1); clf;

subplot(2,1,1);
plot(t,q);

subplot(2,1,2);
plot(t,w);


%%%% Now check using ode45:
tSpan = t([1,end]);
z0 = [q(1); w(1)];

% Dynamical system - driven damped pendulum
dynSys = @(t,z)( [z(2,:);  cos(t) - 0.1*z(2,:) - sin(z(1,:))] );

% simulate:
soln = ode45(dynSys,tSpan,z0);
zSoln = deval(soln,t);
qSoln = zSoln(1,:);
wSoln = zSoln(2,:);

subplot(2,1,1); hold on;
plot(t,qSoln);

subplot(2,1,2); hold on;
plot(t,wSoln);


