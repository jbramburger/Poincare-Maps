% ------------------------------------------------------------------
% SINDy method for discovering mappings in Poincaré sections 
% ------------------------------------------------------------------
% Application to the Rossler system
%
%           x' = -y - z
%           y' = x + a*y
%           z' = b + z*(x-c)
%
% Here a,b,c are real-valued parameters. Fixing a = b = 0.1 and increasing 
% c from 0 causes a sequence of period-doubling bifurcations leading to
% chaos. Attractor at certain parameter values: 
%               c | Attractor
%               -------------------
%               6 | period 2 orbit
%               9 | 2-banded chaos
%            12.6 | period 6 orbit
%              13 | 3-banded chaos
%
% This code is associated with the paper 
% "Poincaré maps for multiscale physics discovery and nonlinear Floquet
% theory" by Jason J. Bramburger and J. Nathan Kutz (Physica D, 2020). 
% This script is used to obtain the results in Section 3.5.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

%Model parameters 
a = 0.1;
b = 0.1;
c = 6;

%ODE generation parameters
m = 3; %Dimension of ODE
n = m-2; %Dimension of Poincaré section
dt = 0.005;
tspan = (0:1000000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

%Generate Trajectories
x0(1,:) = [0; -8; 0]; 
[~,xdat(1,:,:)]=ode45(@(t,x) Rossler(x,a,b,c),tspan,x0(1,:),options);
kfinal = 1; %Number of trajectories
if kfinal >= 2
    for k = 2:kfinal
        x0(k,:) = [0; -4*rand-8; 0];  
        [~,xdat(k,:,:)]=ode45(@(t,x) Rossler(x,a,b,c),tspan,x0(k,:),options);
    end
end

%% Poincare section data

%Counting parameter
count = 1;

%Initialize
Psec = [];
PsecNext = [];

%Create Poincare section data
for i = 1:kfinal
    temp = [];
    for j = 1:length(xdat(i,:,1))-1 
        if  (xdat(i,j,1) < 0 && xdat(i,j+1,1) >= 0) 
            temp(count,:) = xdat(i,j+1,2); %nth iterate
            count = count + 1;
        end
    end
    Psec = [Psec; temp(1:length(temp)-1,:)];
    PsecNext = [PsecNext; temp(2:length(temp),:)];
   	count = 1;
end

%% SINDy for Poincaré Sections

% Access SINDy directory
addpath Util

% Create the recurrence data
xt = Psec(:,1);
xtnext = PsecNext(:,1);

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5; %polynomial order 
usesine = 0; %use sine on (1) or off (0)

Theta = poolData(xt,n,polyorder,usesine);

% compute Sparse regression: sequential least squares
lambda = 0.01;      % lambda is our sparsification knob.

% apply iterative least squares/sparse regression
Xi = sparsifyDynamicsAlt(Theta,xtnext,lambda,n);
if n == 4
[yout, newout] = poolDataLIST({'x','y','z','w'},Xi,n,polyorder,usesine);
elseif n == 3
[yout, newout] = poolDataLIST({'x','y','z'},Xi,n,polyorder,usesine);
elseif n == 2
 [yout, newout] = poolDataLIST({'x','y'},Xi,n,polyorder,usesine);
elseif n == 1 
  [yout, newout] = poolDataLIST({'x'},Xi,n,polyorder,usesine);
end 

fprintf('SINDy model: \n ')
for k = 2:size(newout,2) 
    SINDy_eq = newout{1,k}; 
    SINDy_eq = [SINDy_eq  ' = '];
    new = 1;
   for j = 2:size(newout, 1)
       if newout{j,k} ~= 0 
           if new == 1 
             SINDy_eq = [SINDy_eq  num2str(newout{j,k}) newout{j,1} ];  
             new = 0;
           else 
             SINDy_eq = [SINDy_eq  ' + ' num2str(newout{j,k}) newout{j,1} ' '];
           end 
       end
   end
  fprintf(SINDy_eq)
  fprintf('\n ')
end 

%% Simulate SINDy Map

a = zeros(500,1); %SINDy map solution
a(1) = Psec(1,1);

for k = 1:499
   for j = 1:length(Xi)
      a(k+1) = a(k+1) + Xi(j)*a(k)^(j-1);
   end
   
   if usesine == 1
        a(k+1) = a(k+1) + Xi(j+1)*sin(a(k)) + Xi(j+2)*cos(a(k));
   end
end

%% Plot Results

% Figure 1: Simulations of the discovered Poincaré map
figure(1)
plot(1:100,a(1:100),'b.','MarkerSize',10)
hold on
plot(1:20,Psec(1:20,1),'r.','MarkerSize',10)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Iterates of the Discovered Poincaré Mapping','Interpreter','latex','FontSize',20,'FontWeight','Bold')
legend({'SINDy Mapping','Training Data'}, 'Interpreter','latex','FontSize',16,'Location','best')

%% Rossler right-hand-side

function dx = Rossler(x,a,b,c)

    dx = [-x(2) - x(3); x(1) + a*x(2); b + x(3)*(x(1) - c)];

end










