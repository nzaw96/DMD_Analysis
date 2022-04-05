close all; clear all; clc
dx = 0.02;
L = 20;
x = -L/2:dx:L/2-dx;
%f = 0*x;
f = sech(x);
%plot(x, f)

data = zeros(length(f),110);
dt = 0.025;
for k=1:110
    t = k*dt; 
    [t,u] = ode45(@(t,u)rhsBurgers1b(t,u,L),[0,dt],f);
    f = u(end,:);
    figure(1)
    plot(x,real(f))
    ylabel('u')
    xlabel('x')
    axis([-10 10 -1.5 1.5])
    title('Spectral Solution')
    txt = 't = 2.75s';
    text(6, 0.9,txt,'FontSize',15)
    %pause(0.1)
    data(:,k) = real(f);
end


%model the solution by DMD and compare
 X = data(:, 1:99);
 Xprime = data(:,2:100);
 [U, Sig, V] = svd(X,'econ');
 sig = diag(Sig);
 %setting a threshold for the number of singular values
[s,d] = size(sig);
sig2 = zeros(s, d);
for h=1:s
   if sig(h) > (sig(1)*10^(-12))
       sig2(h) = sig(h);
   else
       break;
   end
end
 r = nnz(sig2);
 [Phi, Lambda, b] = DMDsolution(X, Xprime, r, t);
 lamb = diag(Lambda);
omega = log(lamb)/dt;

%Trying something with r
sig3 = zeros(99, 1);
for c=1:20
   sig3(c) = sig(c);
end
%dt = 0.025;
Ur = U(:,1:r);
Sigr = Sig(1:r,1:r);
Vr = V(:,1:r);

time_dynamics = zeros(r,100);
time = 0;
% x_soln = zeros(r,length(t));
dt = 0.025;
for n=1:100 
     time = n*dt;
     time_dynamics(:,n) = b.*exp(omega*time);
     %x_dmd = Phi*time_dynamics;
     figure(2) 
     plot(x,real(Phi*inv(Lambda)*time_dynamics(:,n))) %added inv(Lambda) in this line
     axis([-10 10 -1.5 1.5])
     xlabel('x')
     ylabel('u')
     title('DMD solution')
     txt = 'r = 10';
     text(6,1.25,txt, 'FontSize',15)
     txt1 = 't = 2.5s';
     text(6, 0.9,txt1,'FontSize',15)
%      ta = annotation('textarrow',[0.4 0.5], [0.679 0.679])
%      ta.String = 'Distortion in wave '
     %pause(0.1)
end
 figure(3)
 x_dmd = Phi*(b.*exp(omega*2.75));
 x_dmd2 = Phi*(b.*exp(omega*2.5));
 plot(x, real(x_dmd))
 title('DMD prediction at t=2.75s')
 xlabel('x')
 ylabel('u')
 axis([-10 10 -1.5 1.5])

%Error estimation from literature
udmd = Phi*(Lambda.^(110))*b; %dmd solution at 2.75secs
I = eye(1000);
B = (Phi*Lambda*pinv(Phi))-I;
K = Xprime*(V*pinv(diag(sig2))*U');
K1 = Xprime*pinv(X);
%K = Xprime*pinv(X); %this might be causing the problem
epsilon_m = norm(K - (Phi*Lambda*pinv(Phi)), 'fro')*norm(data(:,1));
cm = norm(K - (Phi*Lambda*pinv(Phi)), 'fro');
em = data(:,100) - x_dmd2;
err_bound = norm(pinv(Phi), 'fro')*(norm(em)+(10*epsilon_m));
rel_err_bound = err_bound/norm(data(:,1));

myerr = norm(data(:,end) - x_dmd, 'fro');
% syms x y
% x = linspace(1, 100);
% y = cm*(norm(em) + x*epsilon_m);
% figure(4)
% loglog(x, y, 'o')
% xlabel('n')
% ylabel('bounded error')
% %set(gca, 'YScale', 'log')

%IMPLEMENTING OPTIMIZED DMD
% Xopt = data(:, 1:100);
% tvec = linspace(0, 2.5, 100);
% imode = 1;
% [w, e, b2] = optdmd(Xopt, tvec, r, imode);
% Xout = w*diag(b2)*exp(e*tvec);
% e_diag = diag(e);
% dt = 0.025;
% omega2 = log(e)/dt;
% time_dynamics2 = zeros(r, 100);
% time2 = 0;
% figure(4)
% plot(x, real(w*(b2.*exp(omega2*4))))
% axis([-10 10 -1.5 1.5])
% for c=1:100
%     time2 = c*dt;
%     figure(3)
%     plot(x, real(Xout(:,c)))
% %     time_dynamics2 = b2.*exp(omega2*time2);
% %     figure(3)
% %     plot(x, real(w*time_dynamics2(:,c)))
%     axis([-10 10 -1.5 1.5])
%     xlabel('x')
%     ylabel('u')
%     title('Optimized DMD')
% end
%rel_er_dmd  = norm(real(x_dmd) - data(:,end))/norm(data(:,end));
%rel_er_optdmd = norm(real(w*diag(b2)*exp(e*4)) - data(:,end))/norm(data(:,end));

%exact_dmd_er = norm(real(Phi*time_dynamics)-data,'fro')/norm(data,'fro');
%opt_dmd_er = norm(Xout-data, 'fro')/norm(data,'fro');
