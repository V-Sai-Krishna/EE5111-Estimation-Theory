%% Common Initializations
F=zeros(512,32);
h=zeros(32,1);
p=zeros(32,1);
lambda=0.2;
sigma=[0.1,0.01];
X=zeros(512,512);

% F matrix initialized
for rindex = 1:512
    for cindex = 1:32
        F(rindex,cindex)=exp(2*pi*1j*(rindex-1)*(cindex-1)/512);
    end
end

% p[k] (will be used for h[k]) are initialized
p(1:32)=exp(-lambda*(0:31));
P=sqrt(sum(p.*p)/length(p));
[hr,hc]=size(h);
[xr,xc]=size(X);
%% Q1
H=zeros(10000,32);
H_hat=zeros(10000,32);
tmse=zeros(10000,32);
smse=zeros(10000,32);
for iterations = 1:10000
    % sets of 2 random numbers that take either -1 or +1 are created. This
    % is used to create X.
    r=2*randi([0 1],[xr,2])-1;
    X(1:(xr+1):xr*xc)=r(1:xr,1)+r(1:xr,2)*1i;
    % complex noise
    n=sigma(1)*randn(xr,1) + 1i*sigma(1)*randn(xr,1);
    % a[k] and b[k]
    ab=0.5*randn([hr,2]);
    h(1:hr)=(ab(1:hr,1)+ab(1:hr,2)*1i).*p(1:hr)/(P^2);
    % y vector and h using Least Squares are calculated.
    y = (X*F*h) + n;
    h_LSE=((X*F)'*(X*F))\((X*F)'*y);
    % each value of h and h_LSE is stored in a array.
    H(iterations,1:hr)=h(1:hr);
    H_hat(iterations,1:hr)=h_LSE(1:hr);
    tmse(iterations,1:hr)=diag(inv((X*F)'*(X*F)).*(sigma(1)^2));
    smse(iterations,1:hr)=(abs(h-h_LSE)).^2;
end    
%%
figure(1);
subplot(2,1,1);
plot(real(H(3,1:hr)),'o');
hold on;
plot(real(H_hat(3,1:hr)),'.');
hold off;
legend('Real(Actual h)','Real(Estimated h)');
title('Real part of h and estimated h for one random trial');
subplot(2,1,2)
plot(imag(H(3,1:hr)),'o');
hold on;
plot(imag(H_hat(3,1:hr)),'.');
hold off;
legend('Imag(Actual h)','Imag(Estimated h)');
title('Imaginary part of h and estimated h for one random trial');
s1=sum(H)/10000;
s2=sum(H_hat)/10000;

figure(2);
subplot(2,1,1);
plot(real(s1),'o');
hold on;
plot(real(s2),'.');
hold off;
legend('Real(E[h])','Real(h)');
title('Real part of E[h] and h averaged over 10000 trials');
subplot(2,1,2);
plot(imag(s1),'o');
hold on;
plot(imag(s2),'.');
hold off;
legend('Imaginary(E[h])','Imaginary(h)');
title('Imaginary part of E[h] and h averaged over 10000 trials');

figure(3);
m1=sum(tmse)/10000;
m2=sum(smse)/10000;
plot(real(m1),'o');
hold on;
plot(real(m2),'.');
hold off;
legend('Theoretical MSE','Simulated MSE');
title('Theoretical and Simulated MSE, avg over 10000 trials');
%% Q2);
H2=zeros(10000,32);
H2_hat=zeros(10000,32);
tmse2=zeros(10000,32);
smse2=zeros(10000,32);

i2=randperm(32,26);
AC2=zeros(26,32);
for in = 1:26
    AC2(in,i2(in))=1;
end

for iterations = 1:10000
    
    h2=zeros(32,1);
    r=2*randi([0 1],[xr,2])-1;
    X(1:(xr+1):xr*xc)=r(1:xr,1)+r(1:xr,2)*1i;
    
    n=sigma(1)*randn(xr,1) + 1i*sigma(1)*randn(xr,1);
    
    ab=0.5*randn([hr,2]);
    h2(1:hr)=(ab(1:hr,1)+ab(1:hr,2)*1i).*p(1:hr)/(P^2);
    h2(i2(1:26))=0;
    

    y = (X*F*h2) + n;
    h2_LSE=((X*F)'*(X*F))\((X*F)'*y);
    
    lagra2=(AC2*(((X*F)'*(X*F))\AC2'))\(AC2*h2_LSE);
    
    h2_c=h2_LSE-((X*F)'*(X*F))\(AC2'*lagra2);
    
    H2(iterations,1:hr)=h2(1:hr);
    H2_hat(iterations,1:hr)=h2_c(1:hr);
    tmse2(iterations,1:hr)=diag(inv((X*F)'*(X*F)).*(sigma(1)^2));
    smse2(iterations,1:hr)=(abs(h2-h2_c)).^2;
    
end   
%%
figure(2);
subplot(3,1,1);
stem(real(H2(3,1:hr)),real(H2_hat(3,1:hr)));
hold on;
stem(imag(H2(3,1:hr)),imag(H2_hat(3,1:hr)));
legend('Real','Imaginary');
title('Actual h and estimated h for one random trial');

s1=sum(H2)/10000;
s2=sum(H2_hat)/10000;
subplot(3,1,2);
stem(real(s1),real(s2));
hold on;
stem(imag(s1),imag(s2));
legend('Real','Imaginary');
title('E[h] and h averaged over 10000 trials');

m1=sum(tmse2)/10000;
m2=sum(smse2)/10000;
subplot(3,1,3);
plot(real(m1),'.');
hold on;
plot(real(m2),'o');
legend('Theoretical','Simulated');
title('Theoretical and Simulated MSE, avg over 10000 trials');
%% Q3 
X3=zeros(512,512);
H3=zeros(10000,32);
H3_hat=zeros(10000,32);
tmse3=zeros(10000,32);
smse3=zeros(10000,32);
i3=randperm(32,26);
AC3=zeros(26,32);
H32=zeros(10000,32);
H32_hat=zeros(10000,32);
tmse32=zeros(10000,32);
smse32=zeros(10000,32);
H3b_hat=zeros(10000,32);
H3b2_hat=zeros(10000,32);
H3b3_hat=zeros(10000,32);
smse3b=zeros(10000,32);
smse3b2=zeros(10000,32);
smse3b3=zeros(10000,32);
for in = 1:26
    AC3(in,i3(in))=1;
end
warning('off','all');
for iterations = 1:10000

    % Q1 with guard band
    h3=zeros(32,1);
    r=2*randi([0 1],[152,2])-1;
    % X3's diagonal elements only between index 181 and end-180 are
    % non-zero.
    X3((181-1)*xr+181:(xr+1):(512-180-1)*xr+512-180)=r(1:152,1)+r(1:152,2)*1i;
    ab=0.5*randn([hr,2]);
    h3(1:hr)=(ab(1:hr,1)+ab(1:hr,2)*1i).*p(1:hr)/(P^2);
    n=sigma(1)*randn(xr,1) + 1i*sigma(1)*randn(xr,1);
    

    y = (X3*F*h3) + n;
    h3_LSE=((X3*F)'*(X3*F))\((X3*F)'*y);
    H3(iterations,1:hr)=h3(1:hr);
    H3_hat(iterations,1:hr)=h3_LSE(1:hr);
    tmse3(iterations,1:hr)=diag(inv((X3*F)'*(X3*F)).*(sigma(1)^2));
    smse3(iterations,1:hr)=(abs(h3-h3_LSE)).^2;
    % Q2 with guard band
    h32=zeros(32,1);
    ab=0.5*randn([hr,2]);
    h32(1:hr)=(ab(1:hr,1)+ab(1:hr,2)*1i).*p(1:hr)/(P^2);
    h32(i3(1:26))=0;
    y = (X3*F*h32) + n;
    h32_LSE=((X3*F)'*(X3*F))\((X3*F)'*y);
    lagra32=(AC3*(((X3*F)'*(X3*F))\AC3'))\(AC3*h32_LSE);
    h32_c=h32_LSE-((X3*F)'*(X3*F))\(AC3'*lagra32);
    H32(iterations,1:hr)=h32(1:hr);
    H32_hat(iterations,1:hr)=h32_c(1:hr);
    tmse32(iterations,1:hr)=diag(inv((X3*F)'*(X3*F)).*(sigma(1)^2));
    smse32(iterations,1:hr)=(abs(h32-h32_c)).^2;

    % Q1 with guard band and alpha = 1, 1e-3 and 1e-8 
    alpha = [1,1e-3,1e-8];
    I=eye(32);
    y = (X3*F*h3) + n;
    h3b_LSE=((X3*F)'*(X3*F) + alpha(1)*I)\((X3*F)'*y);
    h3b2_LSE=((X3*F)'*(X3*F) + alpha(2)*I)\((X3*F)'*y);
    h3b3_LSE=((X3*F)'*(X3*F) + alpha(3)*I)\((X3*F)'*y);
    H3b_hat(iterations,1:hr)=h3b_LSE(1:hr);    
    H3b2_hat(iterations,1:hr)=h3b2_LSE(1:hr);   
    H3b3_hat(iterations,1:hr)=h3b3_LSE(1:hr);  
    smse3b(iterations,1:hr)=(abs(h3-h3b_LSE)).^2;
    smse3b2(iterations,1:hr)=(abs(h3-h3b2_LSE)).^2;
    smse3b3(iterations,1:hr)=(abs(h3-h3b3_LSE)).^2;
end
%%
figure(3);
subplot(3,1,1);
stem(real(H3(3,1:hr)),real(H3_hat(3,1:hr)));
hold on;
stem(imag(H3(3,1:hr)),imag(H3_hat(3,1:hr)));
legend('Real','Imaginary');
title('Actual h and estimated h for one random trial (Q1) with guard band');

s1=sum(H3)/10000;
s2=sum(H3_hat)/10000;
subplot(3,1,2);
stem(real(s1),real(s2));
hold on;
stem(imag(s1),imag(s2));
legend('Real','Imaginary');
title('E[h] and h averaged over 10000 trials');

m1=sum(tmse3)/10000;
m2=sum(smse3)/10000;
subplot(3,1,3);
plot(real(m1),'.');
hold on;
plot(real(m2),'o');
legend('Theoretical','Simulated');
title('Theoretical and Simulated MSE of estimation, averaged over trials (Q1)');
%%
figure(4);
subplot(3,1,1);
stem(real(H32(3,1:hr)),real(H32_hat(3,1:hr)));
hold on;
stem(imag(H32(3,1:hr)),imag(H32_hat(3,1:hr)));
legend('Real','Imaginary');
title('Actual h and estimated h for one random trial (Q2) with guard band');

s1=sum(H32)/10000;
s2=sum(H32_hat)/10000;
subplot(3,1,2);
stem(real(s1),real(s2));
hold on;
stem(imag(s1),imag(s2));
legend('Real','Imaginary');
title('E[h] and h averaged over 10000 trials');

m1=sum(tmse32)/10000;
m2=sum(smse32)/10000;
subplot(3,1,3);
plot(real(m1),'.');
hold on;
plot(real(m2),'o');
legend('Theoretical','Simulated');
title('Theoretical and Simulated MSE of estimation, averaged over trials (Q2)');
%%
figure(5);
subplot(3,1,1);
stem(real(H3(3,1:hr)),real(H3b_hat(3,1:hr)));
hold on;
stem(imag(H3(3,1:hr)),imag(H3b_hat(3,1:hr)));
legend('Real','Imaginary');
title('Actual h and estimated h for one random trial (alpha=1)');

s1=sum(H3)/10000;
s2=sum(H3b_hat)/10000;
subplot(3,1,2);
stem(real(s1),real(s2));
hold on;
stem(imag(s1),imag(s2));
legend('Real','Imaginary');
title('E[h] and h averaged over 10000 trials');

m1=sum(tmse3)/10000;
m2=sum(smse3b)/10000;
subplot(3,1,3);
plot(real(m1),'.');
hold on;
plot(real(m2),'o');
legend('Theoretical','Simulated');
title('Theoretical and Simulated MSE of estimation, averaged over trials (alpha=1)');
%%
figure(6);
subplot(3,1,1);
stem(real(H3(3,1:hr)),real(H3b2_hat(3,1:hr)));
hold on;
stem(imag(H3(3,1:hr)),imag(H3b2_hat(3,1:hr)));
legend('Real','Imaginary');
title('Actual h and estimated h for one random trial (alpha=1e-3)');

s1=sum(H3)/10000;
s2=sum(H3b2_hat)/10000;
subplot(3,1,2);
stem(real(s1),real(s2));
hold on;
stem(imag(s1),imag(s2));
legend('Real','Imaginary');
title('E[h] and h averaged over 10000 trials');

m1=sum(tmse3)/10000;
m2=sum(smse3b2)/10000;
subplot(3,1,3);
plot(real(m1),'.');
hold on;
plot(real(m2),'o');
legend('Theoretical','Simulated');
title('Theoretical and Simulated MSE of estimation, averaged over trials (alpha=1e-3)');
%%
figure(7);
subplot(3,1,1);
stem(real(H3(3,1:hr)),real(H3b3_hat(3,1:hr)));
hold on;
stem(imag(H3(3,1:hr)),imag(H3b3_hat(3,1:hr)));
legend('Real','Imaginary');
title('Actual h and estimated h for one random trial (alpha=1e-8)');

s1=sum(H3)/10000;
s2=sum(H3b3_hat)/10000;
subplot(3,1,2);
stem(real(s1),real(s2));
hold on;
stem(imag(s1),imag(s2));title('Theoretical and Simulated MSE, avg over 10000 trials');
legend('Real','Imaginary');
title('E[h] and h averaged over 10000 trials');

m1=sum(tmse3)/10000;
m2=sum(smse3b3)/10000;
subplot(3,1,3);
plot(real(m1),'.');
hold on;
plot(real(m2),'o');
legend('Theoretical','Simulated');
title('Theoretical and Simulated MSE of estimation, averaged over trials (alpha=1e-8)');
%% Q4

H4=zeros(10000,32);
H4_hat=zeros(10000,32);
tmse4=zeros(10000,32);
smse4=zeros(10000,32);

AC4=zeros(3,32);
AC4(1,1)=1;
AC4(1,2)=-1;
AC4(2,3)=1;
AC4(2,4)=-1;
AC4(3,5)=1;
AC4(3,6)=-1;
    
    
for iterations = 1:10000
    
    h4=zeros(32,1);
    r=2*randi([0 1],[xr,2])-1;
    X(1:(xr+1):xr*xc)=r(1:xr,1)+r(1:xr,2)*1i;
    
    n=sigma(1)*randn(xr,1) + 1i*sigma(1)*randn(xr,1);
    
    ab=0.5*randn([hr,2]);
    h4(1:hr)=(ab(1:hr,1)+ab(1:hr,2)*1i).*p(1:hr)/(P^2);
    
    y = (X*F*h4) + n;
    h4_LSE=((X*F)'*(X*F))\((X*F)'*y);
    
    lagra4=(AC4*(((X*F)'*(X*F))\AC4'))\(AC4*h4_LSE);
    
    h4_c=h4_LSE-((X*F)'*(X*F))\(AC4'*lagra4);
    
    H4(iterations,1:hr)=h4(1:hr);
    H4_hat(iterations,1:hr)=h4_c(1:hr);
    
    tmse4(iterations,1:hr)=diag(inv((X*F)'*(X*F)).*(sigma(1)^2));
    smse4(iterations,1:hr)=(abs(h4-h4_c)).^2;
end 
%%
figure(8);
subplot(3,1,1);
stem(real(H4(3,1:hr)),real(H4_hat(3,1:hr)));
hold on;
stem(imag(H4(3,1:hr)),imag(H4_hat(3,1:hr)));
legend('Real','Imaginary');
title('Actual h and estimated h for one random trial');

s1=sum(H4)/10000;
s2=sum(H4_hat)/10000;
subplot(3,1,2);
stem(real(s1),real(s2));
hold on;
stem(imag(s1),imag(s2));
legend('Real','Imaginary');
title('E[h] and h averaged over 10000 trials');

m1=sum(tmse4)/10000;
m2=sum(smse4)/10000;
subplot(3,1,3);
plot(real(m1),'.');
hold on;
plot(real(m2),'o');
legend('Theoretical','Simulated');
title('Theoretical and Simulated MSE of estimation, averaged over 10000 trials');
fprintf('h[1],h[2] values that were predicted from a random iteration are %f, %f\n',H4_hat(3,1),H4_hat(3,2));
fprintf('h[3],h[4] values that were predicted from a random iteration are %f, %f\n',H4_hat(3,3),H4_hat(3,4));
fprintf('h[5],h[6] values that were predicted from a random iteration are %f, %f\n',H4_hat(3,5),H4_hat(3,6));
%% Q5

h5=zeros(32,1);
i5=randperm(32,6);
r=2*randi([0 1],[xr,2])-1;
X(1:(xr+1):xr*xc)=r(1:xr,1)+r(1:xr,2)*1i;
n=sigma(1)*randn(xr,1) + 1i*sigma(1)*randn(xr,1);
ab=0.5*randn([6,2]);
h5(i5(1:6))=(ab(1:6,1)+ab(1:6,2)*1i).*p(i5(1:6))/(P^2);
y = (X*F*h5) + n;
S=[];
r5=y;
t=0;
A5=X*F;
for i = 1:6
    [val,t]=max(abs(A5'*r5));
    S=[S,t];
    As=A5(:,S);
    Ats=(As'*As)\As';
    P5=As*Ats;
    r5=(eye(512)-P5)*y;    
end  
%%
AC5=zeros(26,32);
i5=(1:1:32);
setdiff(i5,S);
for in = 1:26
    AC2(in,i5(in))=1;
end
h5_LSE=((X*F)'*(X*F))\((X*F)'*y);
lagra5=(AC5*(((X*F)'*(X*F))\AC5'))\(AC5*h5_LSE);
h5_c=h5_LSE-((X*F)'*(X*F))\(AC5'*lagra5);
%%
figure(9);
stem(real(h5(1:hr)),real(h5_LSE(1:hr)));
hold on;
stem(imag(h5(1:hr)),imag(h5_LSE(1:hr)));
legend('Real','Imaginary');
title('Actual h and estimated h for one random trial');
%%
