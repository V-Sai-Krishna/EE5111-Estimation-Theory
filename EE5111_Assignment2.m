%% Question 1

A=1;
%A=10;
t=50000;

A_est_1 = zeros(1,t);
A_est_10 = zeros(1,t);
A_est_100 = zeros(1,t);
A_est_1000 = zeros(1,t);
A_est_10000 = zeros(1,t);

for i=1:1:t
    
    x = zeros(1,1);
    ni = randn(1,1);
    xi = A + ni;
    A_est_1(i) = mean(xi);

    x=zeros(1,10);
    ni = randn(1,10);
    xi = A + ni;
    A_est_10(i) = mean(xi);


    x=zeros(1,100);
    ni = randn(1,100);
    xi = A + ni;
    A_est_100(i) = mean(xi);


    x=zeros(1,1000);
    ni = randn(1,1000);
    xi = A + ni;
    A_est_1000(i) = mean(xi);


    x=zeros(1,10000);
    ni = randn(1,10000);
    xi = A + ni;
    A_est_10000(i) = mean(xi);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_A_est_1 = mean(A_est_1)
mean_A_est_10 = mean(A_est_10)
mean_A_est_100 = mean(A_est_100)
mean_A_est_1000 = mean(A_est_1000)
mean_A_est_10000 = mean(A_est_10000)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variance_A_est_1 = var(A_est_1)
variance_A_est_10 = var(A_est_10)
variance_A_est_100 = var(A_est_100)
variance_A_est_1000 = var(A_est_1000)
variance_A_est_10000 = var(A_est_10000)

%% Question 1: E[A] and Variance(A) for different N values
fprintf("Question 1:\nN\t E[A]\t\tVariance(A)\n");
fprintf("%d\t%f\t%f\n",1,mean_A_est_1,variance_A_est_1);
fprintf("%d\t%f\t%f\n",10,mean_A_est_10,variance_A_est_10);
fprintf("%d\t%f\t%f\n",100,mean_A_est_100,variance_A_est_100);
fprintf("%d\t%f\t%f\n",1000,mean_A_est_1000,variance_A_est_1000);
fprintf("%d\t%f\t%f\n",10000,mean_A_est_10000,variance_A_est_10000);

%% Question 1: Plots (Please give appropriate inputs for different plots)

[p,x] = hist(A_est_10000); 
figure(1);
plot(x,p/sum(p),'g'); %PDF
hold on;
[p,x] = hist(A_est_100); 
plot(x,p/sum(p),'r'); %PDF
hold on;
[p,x] = hist(A_est_1000); 
plot(x,p/sum(p)); %PDF
%%

[f,x] = ecdf(sqrt(1)*(A_est_1-A)); 
figure(2);
plot(x,f,'r'); %CDF
hold on;

pd = makedist('Normal','mu',0,'sigma',1);
y = -5:.1:5;
p = cdf(pd,y);
plot(y,p,'g');
%% Question 2

N = [1,10,100,1000,10000];
A = 10;
E_A = zeros([5,1]);
iterations = 5000;
for j= 1:5
    U = zeros([N(j),iterations]);
    n = zeros([N(j),iterations]);
    x = zeros([N(j),iterations]);
    for i= 1:iterations
        %generating Laplace distribution
        U(:,i) = rand([N(j),1]) -0.5;
        n(:,i) = -1.*sign(U(:,i)).*log(1-2.*abs(U(:,i)))./(sqrt(2));
        x(:,i) = A + n(:,i);
        A_cap(i,j) = median(x(:,i));
        E_A(j) = E_A(j) + A_cap(i,j);
    end
    E_A(j) = E_A(j)/iterations;
    var(j) = (sum((A_cap(:,j)-E_A(j).*ones([iterations,1])).^2))/iterations;
end

%% Question 2: E[A] and Variance(A) for different N values
fprintf("\nQuestion 2:\nN\t E[A]\t\tVariance(A)\n");
for i = 1:5
    fprintf("%d\t%f\t%f\n",N(i),E_A(i),var(i));
end
%% Question 2: PLOTS %%

x_axis = [1:1:5];
stem(x_axis,E_A,'r','LineWidth',2);
hold on
stem(x_axis,A*ones([5,1]),'b--','LineWidth',1);
ylim(A*[0.9,1.1]);
xlabel('Number of samples N = [1,10,100,1000,10000]');
ylabel('Expectation of MLE of A and A');
title('Comparison of expectation of MLE of A and A');
legend('Expectation of MLE of A', 'A');
%%
stem(x_axis,var,'r','LineWidth',2);
xlabel('Number of samples N = [1,10,100,1000,10000]');
ylabel('Variance of MLE of A');
title('Variance vs number of samples');

%%
[c,x_axis_c] = ecdf(A_cap(:,5));
[c2,x_axis_c] = ecdf(sqrt(10000)*(A_cap(:,5)-E_A(5)*ones([iterations,1])));
gauss_fit = normpdf(x_axis_c,0,1/2);
[c3,x_axis_c] = ecdf(gauss_fit);
plot(x_axis_c,c,'g');
xlabel('A');
ylabel('CDF(A-cap)');
title('CDF of the estimate of A');
%%
plot(x_axis_c,c2,'b');
hold on
plot(x_axis_c,c3, 'r--');
xlabel('x');
ylabel('CDF(x)');
title('Normal convergence of MLE in distribution');
legend('CDF(root(N)*(A-cap-E[A]))', 'CDF(N(0,1/I(A)))');
%%
[pdf_A,x_axis_c]=hist(A_cap(:,5));
plot(x_axis_c,pdf_A/N(5),'r');
xlabel('A cap');%% Question 2: E[A] and Variance(A) for different N values
fprintf("\nQuestion 2:\nN\t E[A]\t\tVariance(A)\n");
for i = 1:5
    fprintf("%d\t%f\t%f\n",N(i),E_A(i),var(i));
end
ylabel('PDF(A cap)');
title('PDF of the MLE of A');
%%
[pdf,x_axis_c]=hist(sqrt(N(5))*(A_cap(:,5)-E_A(5)*ones([iterations,1])));
plot(x_axis_c,pdf/iterations,'r');
hold on
gauss_fit = normpdf(x_axis_c,0,1/2);
plot(x_axis_c,gauss_fit,'b');
xlabel('x');
ylabel('PDF(x)');
title('Normal convergence of MLE in distribution');
legend('PDF(root(N)*(A-cap-E[A]))', 'PDF(N(0,1/I(A)))');

%% Question 3
A=10;
N=[1,10,100,1000,10000];
Cg=1.78;
gamma=(2*Cg)^0.5;
E_A = zeros([5,1]);
var_A = zeros([5,1]);
A_est=zeros(5,10000);
iterations=10000;

for i=1:5
    for j =1:iterations
        noise=gamma*tan(pi*(rand(N(i),1)-1/2));
        X=noise+A;
        if(N(i)>1)
            c=sort(X);
            A_est0=c(floor(N(i)/2));
        else
            A_est0=X;
        end
        for k = 1:10
            fder=sum(2*(X-A_est0)./(gamma^2 + (X-A_est0).^2));
            sder= sum(2*(2*((X-A_est0).^2)./(gamma^2+(X-A_est0).^2).^2 - (gamma^2+(X-A_est0).^2).^(-1)));
            A_est0=A_est0-(fder/sder);
        end
        A_est(i,j)=A_est0;
        E_A(i)=E_A(i)+A_est0;
    end
    E_A(i)=E_A(i)/iterations;
    var_A(i) = (sum((A_est(i,:)-E_A(i)).^2))/iterations;
end
%% Question 3: E[A] and Variance(A) for different N values
fprintf("\nQuestion 3:\nN\t E[A]\t\tVariance(A)\n");
for i = 1:5
    fprintf("%d\t%f\t%f\n",N(i),E_A(i),var_A(i));
end
%% Question 3: Plots (Please give appropriate inputs for different plots)
figure(3);
[p,x] = hist(A_est(1,:)); 
plot(x,p/sum(p));%PDF
grid();
xlabel('A cap');
ylabel('PDF(A cap)');
title('PDF of the MLE of A');
figure(4);
[f,x] = ecdf(sqrt(N(5))*(A_est(5,:)-A)); 
plot(x,f,'r'); %CDF
grid();
hold on;
pd = makedist('Normal','mu',0,'sigma',gamma*(2^0.5));
y = -5:.1:5;
p = cdf(pd,y);
plot(y,p,'g');
hold off;
legend('sqrt(N)*(A−A0)', 'Gaussian(0,I(A)^-1)');
title('CDF of estimate and Gaussian(0,I(A)^-1)');