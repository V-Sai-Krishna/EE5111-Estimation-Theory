%% Generating Y's
mu = zeros(2,1);
sigma = [[1,0];[0,2]];
n=1000;
R = mvnrnd(mu,sigma,n);
R = transpose(R);
%% Question 1 
% MLE
fprintf("Covariance Estimated using MLE:\n");
mean = sum(R,2)/n;

X = zeros(2,2);
for i=1:n

   X = X + (R(:,i)-mean)*transpose(R(:,i)-mean);

end
covariance = X/n;
fprintf("For n=%d:\n",n);
disp(covariance);
%% Question 2 and 3
mu = zeros(2,1);
sigma = [[1,0];[0,2]];
n=10;
R = mvnrnd(mu,sigma,n);
R = transpose(R);
%%% prior distribution parameters
df = 5;
Tau = [[4,0];[0,5]];
fprintf("For n=%d:\n",n);
fprintf("Covariance Estimated using Conjugate Prior:\n");
%%% posterior distribution parameters
df1 = df + n;
sum1 = zeros(2,2);
for i=1:n
    sum1 = sum1 + (R(:,i)*transpose(R(:,i)));
end
Tau1 = Tau + sum1;R = mvnrnd(mu,sigma,n);
R = transpose(R);

%%% estimate of the covariance
covariance_est_2 = iwishrnd(Tau1,df1);
disp(covariance_est_2);
%% Question 3 %%%%%%%%%%%%%%%

%%% posterior under non informative jeffrey's prior
%%% estimated sigma has inverse wishart distribution where delta_0=0 and
%%% v_0=1

fprintf("Covariance Estimated using non informative jeffrey's prior:\n")
covariance_est_3a = iwishrnd(sum1,n+1);
disp(covariance_est_3a);

%%% posterior under independent Jeffrey's prior
%%% estimated sigma has inverse wishart distribution where deta_0=0 and
%%% v_0=0

fprintf("Covariance Estimated using independent Jeffrey's prior:\n")
covariance_est_3b = iwishrnd(sum1,n);
disp(covariance_est_3b);

%% Question 4
mu = zeros(2,1);
sigma = [[1,0];[0,2]];
n=10;
R = mvnrnd(mu,sigma,n);
R = transpose(R);

% Inv Wishart Distribution %

df = 5;
Tau = [[4,0];[0,5]];
m=10000;
prior = zeros(2,2,m);
for i=1:m
    
    prior(:,:,i) = iwishrnd(Tau,df);

end

numerator = zeros(2,2);
denominator = 0;
P = zeros(m,1);

for i=1:m
    
    P(i) = sum(sum(transpose(R)/(prior(:,:,i)).*transpose(R)));
    denominator = denominator + det(prior(:,:,i))^(-n/2)*exp(-0.5*P(i))/m;
    numerator = numerator + prior(:,:,i)*det(prior(:,:,i))^(-n/2)*exp(-0.5*P(i))/m;
    
end

A = numerator/denominator;
fprintf("For n=%d:\n",n);
fprintf("Covariance Estimated using Monte Carlo:\n");
disp(A);
%% Q5
n=1000;
meany=[0,0];
sigmay=[1 0;0 2];
y=mvnrnd(meany,sigmay,n);

iter=1000;
A1=0.05;
A2=0.05;
a1=(gamrnd(0.5,1/(A1^2)))^-1;
a2 = (gamrnd(0.5,1/(A2^2)))^-1;
v=2;
d=2;
for i = 1:iter
    Tau=2*v*diag([a1,a2]) + y'*y;
    sig=iwishrnd(Tau,v+n+d-1);
    a1=1/gamrnd((v+n)/2,v*sig(1,1)+1/(A1^2));
    a2=1/gamrnd((v+n)/2,v*sig(2,2)+1/(A2^2));
end
fprintf("For n=%d:\n",n);
fprintf("Covariance Estimated using Gibbs Estimation:\n");
disp(sig);
%% Q6
n=1000;
meany=[0,0];
sigmay=[1 0;0 2];
y=mvnrnd(meany,sigmay,n);

del0=[4 0;0 5];
v=5;
diff=10^5;
d=2;
vvalues=[];
objvalues=[];
while diff>1e-10
    nterm=0;
    dterm=0;
    objterm=0;
    for i=1:d
        for j = 1:(n/2)
            objterm = objterm + log((v+n-i+1-2*j)/2);
            nterm=nterm+1/(v+n-i+1-2*j);
            dterm=dterm+1/((v+n-i+1-2*j)^2);
        end
    end
    obj=v*log((v+n)/v) + n*log((v+n)/v) - objterm;
    objvalues=[objvalues,obj];
    vvalues=[vvalues,v];
    vtemp=v-(d*log((v+n)/v)/2 - nterm)/(dterm - d*n/(2*v*(v+n)));
    diff=abs(vtemp-v);
    v=vtemp;
end

del=v*(y'*y)/n;
sig=iwishrnd(del+y'*y,v+n);
fprintf("For n=%d:\n",n);
fprintf("Covariance Estimated using Emperical Bayes:\n");
disp(sig);
%%
n=1000;
meany=[0,0,0];
sigmay=[1 0 0;0 2 0;0 0 3];
y=mvnrnd(meany,sigmay,n);

iter=1000;
A1=0.05;
A2=0.05;
A3=0.05;
a1=(gamrnd(0.5,1/(A1^2)))^-1;
a2 = (gamrnd(0.5,1/(A2^2)))^-1;
a3=(gamrnd(0.5,1/(A3^2)))^-1;
v=2;
d=3;
for i = 1:iter
    Tau=2*v*diag([a1,a2,a3]) + y'*y;
    sig=iwishrnd(Tau,v+n+d-1);
    a1=1/gamrnd((v+n)/2,v*sig(1,1)+1/(A1^2));
    a2=1/gamrnd((v+n)/2,v*sig(2,2)+1/(A2^2));
    a3=1/gamrnd((v+n)/2,v*sig(2,2)+1/(A3^2));
end
fprintf("For n=%d:\n",n);
fprintf("Covariance Estimated using Gibbs Estimation:\n");
disp(sig);