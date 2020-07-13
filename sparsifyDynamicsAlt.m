function Xi = sparsifyDynamicsAlt(Theta,dXdt,lambda,n)
% Augmentation of the code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

%SVD
[U,S,V]=svd(Theta,0); 
Xi = V*((U'*dXdt)./diag(S));

%Xi = pinv(Theta)*dXdt;
k = 1;
Xi_new = Xi;

% lambda is our sparsification knob.
while sum(sum(abs(Xi - Xi_new))) > 0  || k == 1 
    Xi = Xi_new;
    smallinds = (abs(Xi)<lambda);   % find small coefficients
    Xi_new(smallinds)=0;                % and threshold
    for ind = 1:n                   % n is state dimension
        biginds = ~smallinds(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        [U,S,V]=svd(Theta(:,biginds),0);
        Xi_new(biginds,ind) = V*((U'*dXdt(:,ind))./diag(S));
    end
    k = k + 1;
end
 Xi = Xi_new;
 
 
 
 
 
 
 
