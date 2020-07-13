% ------------------------------------------------------------------
% SINDy method for discovering mappings in Poincaré sections 
% ------------------------------------------------------------------
% Application to the non-autonomous RC circuit equation
%
%           x' = A*sin(omega*t)-x,      A,omega > 0 
% 
% The Poincaré section can be computed explicitly and is linear. 
%
% This code is associated with the paper 
% "Poincaré maps for multiscale physics discovery and nonlinear Floquet
% theory" by Jason J. Bramburger and J. Nathan Kutz (Physica D, 2020). 
% This script is used to obtain the results in Section 3.1.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

%Model parameters 
A = 1;
omega = 2*pi;

%Generate Trajectories Hopf normal form
m = 2; %Dimension of ODE
n = m-1; %Dimension of Poincaré section
dt = 0.01;
tspan = (0:20000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

%Generate More Trajectories
x0(1,:) = [0; 0]; %At least one solution
[~,xdat(1,:,:)]=ode45(@(t,x) RC(x,A,omega),tspan,x0(1,:),options);
kfinal = 5; %Number of trajectories
if kfinal >= 2
    for k = 2:kfinal
        x0(k,:) = [10*rand-5; 0]; %Initial conditions start in section 
        [~,xdat(k,:,:)]=ode45(@(t,x) RC(x,A,omega),tspan,x0(k,:),options);
    end
end

%% Poincaré section data

%Counting parameter
count = 1;

%Initialize
Psec(1) = xdat(1,1,1);
count = count + 1;

%Create Poincaré section data
for i = 1:kfinal
    for j = 1:length(xdat(i,:,1))-1 
        if (j == 1) && (i > 1) %Trajectories start in the section
            Psec(count) = xdat(i,j,1);
            count = count + 1; 
        elseif  (mod(xdat(i,j,2),2*pi/omega) >= 2*pi/omega-dt && mod(xdat(i,j+1,2),2*pi/omega) <= dt) 
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
xt = Psec(1:length(Psec)-1)';
xtnext = Psec(2:length(Psec))';

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5; %polynomial order 
usesine = 0; %use sine on (1) or off (0)

Theta = poolData(xt,n,polyorder,usesine);

% compute Sparse regression: sequential least squares
lambda = 0.005;      % lambda is our sparsification knob.

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
b = zeros(length(Psec),1);
a(1) = x0(1,1,1);
b(1) = 3*rand;

% Simulate section
for k = 1:length(Psec)-1
   for j = 1:length(Xi)
      a(k+1) = a(k+1) + Xi(j)*a(k)^(j-1);
      b(k+1) = b(k+1) + Xi(j)*b(k)^(j-1);
   end
end

%% Plot Solutions

% Figure 1: Continuous-time solution
figure(1)
plot(tspan,xdat(1,:,1),'b')
set(gca,'FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Solution of the ODE with $x(0) = 0$','Interpreter','latex','FontSize',20,'FontWeight','Bold')

% Figure 2: Simulations of the discovered Poincaré map
figure(2)
plot(1:100,a(1:100),'b.','MarkerSize',10)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x(n)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Iterates of the Discovered Poincaré Mapping','Interpreter','latex','FontSize',20,'FontWeight','Bold')

%% Rc Circuit right-hand-side

function dx = RC(x,A,omega)

dx = [A*sin(omega*x(2)) - x(1); 1];

end


