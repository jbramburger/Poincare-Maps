% ------------------------------------------------------------------
% SINDy method for discovering mappings in Poincaré sections 
% ------------------------------------------------------------------
% Application to coarse-grained dynamical evolution of spiral 
% wave solutions to the lambda-omega PDE
%
%       u_t = D?u + u*(1 - u^2 + v^2) - beta*v*(u^2 + v^2)
%       v_t = D?v + v*(1 - u^2 + v^2) + beta*u*(u^2 + v^2)
%               
% Here D,beta > 0 are real-valued parameters. Data is gathered using 
% D = 0.1 and beta = 1 on a spatial domain with x,y in [-10,10] with 
% periodic boundary conditions. Simulations are performed using spectral
% methods.
%
% Solutions are simulated and projected onto their two dominant PCA modes. 
% Temporal data for the PCA modes is contained in spiral_data.mat.
%
% This code is associated with the paper 
% "Poincaré maps for multiscale physics discovery and nonlinear Floquet
% theory" by Jason J. Bramburger and J. Nathan Kutz (Physica D, 2020). 
% This script is used to obtain the results in Section 3.6.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

%Initializations
n = 2; %Using only the two principal components from u

%% Aggregate Data

load spiral_pca_series.mat

%Create section data
for j = 1:80
   xt(j,:) = [pcaSeries_u(5*(j-1)+1,:)]; %pcaSeries_v(10*(j-1)+1,:)];
   xtnext(j,:) = [pcaSeries_u(5*j+1,:)]; %pcaSeries_v(10*j+1,:)];
end

%% SINDy for Coarse-Grained Forecasting

% Access SINDy directory
addpath Util

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5; %polynomial order 
usesine = 0; %use sine on (1) or off (0)

Theta = poolData(xt,n,polyorder,usesine);

% compute Sparse regression: sequential least squares
lambda = 10;      % lambda is our sparsification knob.

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



%% Simulate Poincare Map

a = zeros(1000,1); %SINDy map solution
b = zeros(1000,1);
a(1) = xt(1,1);
b(1) = xt(1,2);

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

%% Uncertainty Quantification

sample = 5000;

A(:,1) = normrnd(a(1),0.1,sample,1);
B(:,1) = normrnd(b(1),0.1,sample,1);

for k = 1:999
    
    % Constant terms
    A(:,k+1) = Xi(1,1)*ones(sample,1);
    B(:,k+1) = Xi(1,2)*ones(sample,1);
    
    %Polynomial terms
   for p = 1:polyorder
       for j = 0:p
           A(:,k+1) = A(:,k+1) + Xi(1 + j + p*(p+1)/2,1)*(A(:,k).^(p-j)).*(B(:,k).^j);
           B(:,k+1) = B(:,k+1) + Xi(1 + j + p*(p+1)/2,2)*(A(:,k).^(p-j)).*(B(:,k).^j);
       end
   end
end

%Eliminate blow-up terms
A(any(isnan(A), 2), :) = [];
B(any(isnan(B), 2), :) = [];

% Figure 1: Simulations of the discovered mapping
figure(1)
plot(-A(:,1:50)','k.')
hold on
plot(-a(1:50),'b.','MarkerSize',40)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$\tilde{x}(n)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Iterates of the Discovered Mapping: Component 1','Interpreter','latex','FontSize',20,'FontWeight','Bold')

% Figure 1: Simulations of the discovered mapping
figure(2)
plot(B(:,1:50)','k.')
hold on
plot(b(1:50),'r.','MarkerSize',40)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$\tilde{y}(n)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Iterates of the Discovered Mapping: Component 2','Interpreter','latex','FontSize',20,'FontWeight','Bold')












