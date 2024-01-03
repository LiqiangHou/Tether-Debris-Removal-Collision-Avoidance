function [ybar,Py]=UT(xbar,Px,tau_0,tau_f)
%Unscented transformation takes input of mean x bar and covariance Px,
%along with scaling parameters alpha, kappa and beta
alpha=1e-3;                                 %default, tunable
kappa=0;                                    %default, tunable
beta=2;                                     %default, tunable




% x bar is a row vector of dimension L
L=length(xbar)
lambda=alpha^2*(L+kappa)-L

%Finding Chi_i vectors
Chi=zeros(2*L+1,L) %create matrix of dimension 2L+1 x L
Chi(1,:)=xbar

%define matrix sqrt

matrix_sqrt=sqrtm((L+lambda)*Px)
% Obtain Chi MATRIX
for i=1:(2*L);
    if i<=L
        Chi(i+1,:)=xbar+matrix_sqrt(i,:); % CHi_i= xbar + ith row of matrix sqrt
    elseif i>L
        Chi(i+1,:)=xbar+matrix_sqrt((i-L),:);
    end
end

%Obtain Weight matrix

W0_mean=lambda/(lambda+L)

W0_co=lambda/(lambda+L)+(1-alpha^2+beta)

Wi_mean=1/(2*(L+lambda))


%propagate CHI through the nonlinear function g, currently use sine
%function as dummy

ybar=zeros(1,L);


for i=1:(2*L+1);
    ybar=ybar+Wi_mean*fitness(tau_0,tau_f,Chi(i,:));
end

%Obtaining Py

Py=zeros(L,L);

for i=1:(2*L+1);
    Py=Py+Wi_mean*((fitness(tau_0,tau_f,Chi(i,:))-ybar)*transpose(fitness(tau_0,tau_f,Chi(i,:))-ybar));

end
end
