clear all; close all; clc
Nx = 499;
x = linspace(0, 1, Nx+1);
dx = x(2) - x(1);
Nt = 100000;
t = linspace(0, 0.2, Nt+1);
dt = t(2) - t(1);
F = dt/(dx)^2;
%f1 = dt/(0.001)^2;
upre = zeros(length(x), 1); %initial condition
upre(end) = 1;  %right end BC enfored on initial condition
u = zeros(length(x), 1);
data = zeros(length(u), length(t));
data(:,1) = upre(:);
for m=1:length(t)
    for i=2:(Nx)
        u(i) = upre(i) + F*(upre(i+1) + upre(i-1) -2*upre(i));
    end
    u(1) = 0; %BC1
    u(Nx+1) = 1; %BC2
    upre(:) = u(:);
    if (m==1)
        continue
    else
        data(:,m) = u(:); 
    end
end

%Here I will feed the data from finite solution to the DMD algorithm
data1 = data(:, 1:200:end);
X1 = data1;
%X1prime = data1(:, 2:end);
t1 = t(1:200:end);

%Looking at different X's and Xprime's to improve the condition #s
% X_ = X1 + eye(500,500);
% X__  = X_(:, 1:199);
% Xprime__ = X_(:,2:200);

X = X1(:,1:199);
Xprime = X1(:,2:200);
[U, Sig, V] = svd(X, 'econ');
sig = diag(Sig);
[s,d] = size(sig);
sig2 = zeros(s, d);
for h=1:s
   if sig(h) > (sig(1)*10^(-8))
       sig2(h) = sig(h);
   else
       break;
   end
end
r1 = nnz(sig2);
sig3 = sig2(1:r1);
%r = 20;
[Phi, Lambda, b] = DMDsolution(X, Xprime, r1, t);
lamb = diag(Lambda);
omega = log(lamb)/dt;
dmd_sol = Phi*diag(exp(omega*1))*b;
sol_ = Phi*(Lambda.^(100))*b;
myerr = norm(dmd_sol - data1(:,100));


dt1 = t1(2)- t1(1);
[Phi2, lambda2, b2] = DMDsolution(X, Xprime, r1, dt1);

% figure(1)
% plot(Xdmd2(:,100))

%plot for myerr
snapshot = [201 210 250 300 350 400 425 500]; %snapshots; interpolation mode has 1-200 snapshots
er = [];
% figure(1)
% plot(snapshot, er, 'r-')
% xlabel('n')
% ylabel('e')
% title('L2 error taken between finite differnece and DMD solution')
% hold on
% scatter(snapshot, er, 'o', 'MarkerFaceColor', 'b')

%Putting Error bound from literature here
format long;
%K = Xprime*(V*pinv(diag(sig3))*U');
[U_, S_, V_] = svd(Xprime);
Xprime11 = U_(:,1:r1)*S_(1:r1,1:r1)*(V_(1:r1,1:r1)');
K = Xprime*pinv(X); %this might be causing the problem
X11 = U(:,1:r1)*diag(sig3)*(V(1:r1,1:r1))';
K1 = Xprime11*pinv(X11); 
epsilon_m = norm(K - (Phi*Lambda*pinv(Phi)), 'fro')*norm(data(:,1));
cm = norm(K - (Phi*Lambda*pinv(Phi)), 'fro'); 
x_dmd2 = Phi*(b.*exp(omega*0.0796)); %time here corresponds to 200th snapshot;last snapshot in interpolation mode
em = data(:,200) - x_dmd2; 
err_bound = norm(pinv(Phi), 'fro')*(norm(em)+(2*epsilon_m));

snapshots = [200 ];
trun_er = [2.008529534303560 ];

%implementing longer version of error bound

longerbound = norm(Phi*(Lambda^1)*pinv(Phi), 'fro')*norm(em) + 1*epsilon_m*norm(Phi*pinv(Phi), 'fro');%, norm(Phi*(Lambda^1)*pinv(Phi)));

[A, B, b] = DMDsolution(X, Xprime, 5, t); 

%Here I will do the analytical solution
uderiv = zeros(length(u), 1);
uvec = zeros(length(u),1);
sum = 0;
sumderiv = 0;
for z=1:length(u)
   for n=1:5000
        sum = sum + (((-1)^n)/n)*sin(n*pi*(dx*(z-1)))*exp(-1*(n^2)*(pi^2)*0.2);%input time value explicitly here
        sumderiv = sumderiv + ((-1)^n)*(n*pi)*sin(n*pi*(dx*(z-1)))*exp(-1*(n^2)*(pi^2)*0.1998); %deriv to find error for CDM
   end
   uanlyt = (dx*(z-1)) + (2/pi)*sum;
   uanlytderiv = -2*sumderiv;
   uvec(z) = uanlyt;
   uderiv(z) = uanlytderiv;
   sumderiv = 0;
   uanlytderiv = 0;
   sum = 0;
   uanlyt = 0;
end


