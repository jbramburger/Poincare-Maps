% ------------------------------------------------------------------
% SINDy method for discovering mappings in Poincaré sections 
% ------------------------------------------------------------------
% Application to the driven Brusselator
%
%           x' = a + alpha*sin(2*pi*t/T) - (b+1)*x + x^2*y
%           y' = b*x - x^2*y     
%
% Here a,b,alpha are real-valued parameters and T > 0 controls the 
% period of forcing.
%
% This code is associated with the paper 
% "Poincaré maps for multiscale physics discovery and nonlinear Floquet
% theory" by Jason J. Bramburger and J. Nathan Kutz (Physica D, 2020). 
% This script is used to obtain the results in Section 3.4.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

%Model parameters
a = 0.4;
b = 1.2;
T = 1; 
alpha = 0.1;

%ODE generation parameters
m = 3; %Dimension of ODE
n = m-1; %Dimension of Poincaré section
dt = 0.005;
tspan = (0:100000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

%Generate Trajectories
x0(1,:) = [0; 0; 0]; 
[~,xdat(1,:,:)]=ode45(@(t,x) Brusselator(x,a,b,T,alpha),tspan,x0(1,:),options);
kfinal = 2;
if kfinal >= 2
    for k = 2:kfinal
        x0(k,:) = [4*rand; 4*rand; 0]; %Initial conditions start in section 
        [~,xdat(k,:,:)]=ode45(@(t,x) Brusselator(x,a,b,T,alpha),tspan,x0(k,:),options);
    end
end

%% Poincaré section data

%Counting parameter
count = 1;

%Initialize
Psec = [];
PsecNext = [];

%Create Poincaré section data
for i = 1:kfinal
    for j = 1:length(xdat(i,:,1))-1 
        if  (mod(xdat(i,j,3),T) >= T-dt && mod(xdat(i,j+1,3),T) <= dt) %&& j >= length(xdat(i,:,1))/50) 
            temp(count,:) = xdat(i,j+1,1:2); %nth iterate
            count = count + 1;
        end
    end
    Psec = [Psec; temp(1:length(temp)-1,:)];
    PsecNext = [PsecNext; temp(2:length(temp),:)];
   	count = 1;
    temp = [];
end

%% SINDy for Poincaré Sections

% Access SINDy directory
addpath Util

% Create the recurrence data
xt = Psec;
xtnext = PsecNext;

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

a = zeros(1000,1); %SINDy map solution
b = zeros(1000,1);
a(1) = Psec(1,1);
b(1) = Psec(1,2);

for k = 1:999
    
    % Constant terms
    a(k+1) = Xi(1,1);
    b(k+1) = Xi(1,2);
    
    %Polynomial terms
   for p = 1:polyorder
       for j = 0:p
           a(k+1) = a(k+1) + Xi(1 + j + p*(p+1)/2,1)*(a(k)^(p-j))*(b(k)^j);
           b(k+1) = b(k+1) + Xi(1 + j + p*(p+1)/2,2)*(a(k)^(p-j))*(b(k)^j);
       end
   end
   
   if usesine == 1
        a(k+1) = a(k+1) + Xi((p+1)*p/2+p+2,1)*sin(a(k)) + Xi((p+1)*p/2+p+3,1)*sin(b(k))+ Xi((p+1)*p/2+p+4,1)*cos(a(k)) + Xi((p+1)*p/2+p+5,1)*cos(b(k));
        b(k+1) = b(k+1) + Xi((p+1)*p/2+p+2,2)*sin(a(k)) + Xi((p+1)*p/2+p+3,2)*sin(b(k))+ Xi((p+1)*p/2+p+4,2)*cos(a(k)) + Xi((p+1)*p/2+p+5,2)*cos(b(k));
   end
end

%% Plot Results

% Figure 1: Continuous-time solution x(t)
figure(1)
plot(tspan,xdat(1,:,1),'b','LineWidth',2)
set(gca,'FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Solution of the ODE','Interpreter','latex','FontSize',20,'FontWeight','Bold')

% Figure 2: Continuous-time solution y(t)
figure(2)
plot(tspan,xdat(1,:,2),'r','LineWidth',2)
set(gca,'FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Solution of the ODE','Interpreter','latex','FontSize',20,'FontWeight','Bold')

% Figure 3: Simulations of the discovered Poincaré map
figure(3)
plot(1:100,a(1:100),'b.','MarkerSize',10)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Iterates of the Discovered Poincaré Mapping','Interpreter','latex','FontSize',20,'FontWeight','Bold')

% Figure 4: Simulations of the discovered Poincaré map
figure(4)
plot(1:100,b(1:100),'r.','MarkerSize',10)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Iterates of the Discovered Poincaré Mapping','Interpreter','latex','FontSize',20,'FontWeight','Bold')

%% Driven Brusselator right-hand-side
function dx = Brusselator(x,a,b,T,alpha)

    dx = [a + alpha*sin(2*pi*x(3)/T) - b*x(1) + x(1)*x(1)*x(2) - x(1); b*x(1) - x(1)*x(1)*x(2); 1];

end










