% ------------------------------------------------------------------
% SINDy method for discovering mappings in Poincaré sections 
% ------------------------------------------------------------------
% Application to the Hopf normal form
%
%           x' = x - omega*y - x*(x^2 + y^2)
%           y' = omega*x + y - y*(x^2 + y^2)     
%
% Here omega > 0 is a real-valued parameter.
%
% This code is associated with the paper 
% "Poincaré maps for multiscale physics discovery and nonlinear Floquet
% theory" by Jason J. Bramburger and J. Nathan Kutz (Physica D, 2020). 
% This script is used to obtain the results in Section 3.2.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

%Model parameters
T = 10*pi; % Period 

%Generate Trajectories Hopf normal form
m = 2; %Dimension of ODE
n = m-1; %Dimension of Poincaré section
dt = 0.005;
tspan = (0:10000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

%Generate Trajectories
x0(1,:)=[0.0001; 0];  % Initial condition (Guarantee one starts near origin)
x0(2,:) = [0; 0]; %Trajectory at the ubstable equilibrium
[~,xdat(1,:,:)]=ode45(@(t,x) Hopf(x,T),tspan,x0(1,:),options);
[~,xdat(2,:,:)]=ode45(@(t,x) Hopf(x,T),tspan,x0(2,:),options);
kfinal = 5; % Number of trajectories
if kfinal >= 3
    for k = 3:kfinal
        x0(k,:) = [20*rand; 0]; %Initial conditions start in section 
        [~,xdat(k,:,:)]=ode45(@(t,x) Hopf(x,T),tspan,x0(k,:),options);
    end
end

%% Poincaré section data

%Counting parameter
count = 1;

%Initialize
Psec(1) = xdat(1,1,1);
count = count + 1;

%Create Poincare section data
for i = 1:kfinal
    for j = 1:length(xdat(1,:,1))-1 
        if  (xdat(i,j,2) < 0) && (xdat(i,j+1,2) >= 0) 
            Psec(count) = xdat(i,j+1,1); %nth iterate
            PsecNext(count - 1) = xdat(i,j+1,1); %(n+1)st iterate
            count = count + 1;
        end
    end
    Psec = Psec(1:length(Psec)-1);
    count = count - 1;
end

%% SINDy for Poincaré Sections

% Access SINDy directory
addpath Util

% Create the recurrence data
xt = Psec';
xtnext = PsecNext';

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

a = zeros(length(Psec),1); %SINDy map solution
a(1) = 2*rand;

for k = 1:length(Psec)-1
   for j = 1:length(Xi)
      a(k+1) = a(k+1) + Xi(j)*a(k).^(j-1);
   end
end

%% Plot Results

% Figure 1: Simulations of the discovered Poincaré map
figure(1)
plot(1:100,a(1:100),'b.','MarkerSize',10)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Iterates of the Discovered Poincaré Mapping','Interpreter','latex','FontSize',20,'FontWeight','Bold')


%% Hopf Right-hand-side
function dx = Hopf(x,T)

    dx = [x(1) - T*x(2) - x(1)*(x(1)^2 + x(2)^2); T*x(1) + x(2) - x(2)*(x(1)^2 + x(2)^2)];

end






