% BackwardEuler TE
clear all

syms nu dt u ux uxx uxxx uxxxx uxxxxx uxxxxxx dt1 dt2 dt3 dt4 u_m1 u_m2 x1 x2 x3 x4 u1 u2 u3 u4 reals %ut uxt uxxt uxxxt uxxxxt uxxxxxt uxxxxxxt c1 c2 c3 f fx fxx fxxx fxxxx fxxxxx nu ubar reals
% clc


%% analytic derivation
if 0
% Taylor series
u_m1 = u  - ux*( dt )   +  1/2*uxx*( dt )^2  -  1/6*uxxx*( dt )^3 ;% +    1/24*uxxxx*(dt1)^4;%  +  1/factorial(5)*uxxxxx*(dt1)^5;
u_m2 = u  - ux*( 2*dt )   +  1/2*uxx*( 2*dt )^2  -  1/6*uxxx*( 2*dt )^3 ;%+     1/24*uxxxx*(dt1+dt2)^4;%  +  1/factorial(5)*uxxxxx*(dt1+dt2)^5;


u_m1t = u  - ux*( (dt1+dt2)/2 )   +  1/2*uxx*( (dt1+dt2)/2 )^2  -  1/6*uxxx*( (dt1+dt2)/2 )^3 ;% +    1/24*uxxxx*(dt1)^4;%  +  1/factorial(5)*uxxxxx*(dt1)^5;
u_m2t = u  - ux*( dt1/2+dt2+dt3/2 )   +  1/2*uxx*(dt1/2+dt2+dt3/2)^2  -  1/6*uxxx*(dt1/2+dt2+dt3/2)^3 ;%+     1/24*uxxxx*(dt1+dt2)^4;%  +  1/factorial(5)*uxxxxx*(dt1+dt2)^5;


fprintf('Truncation error for time integration scheme:\n tau_h(u) = ');
fnph = u + dt1/(dt1+dt2)*(u - u_m1t);
fnmh = u_m1t + dt2/(dt2+dt3)*(u_m1t - u_m2t);
dfdt2 = (fnph-fnmh)/dt1;
pretty(simplify(expand(dfdt2-ux)));

fprintf('Equal Spacing (dt = dt1 = dt2 = dt3): \n tau_h(u) = ');
dfdt2eq = subs(dfdt2,{dt1,dt2,dt3},{dt,dt,dt});
pretty(simplify(expand(dfdt2eq-ux)));

fprintf('Direct calculation assuming equal spacing: \n tau_h(u) = ');
dfdt = (3*u - 4*u_m1 + u_m2)/(2*dt);
pretty(simplify(expand(dfdt-ux)));
end

%% Numerical Tests
% Test setup
% Burgers equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
re=32; l=8; u=2; nu = l*u/re;
c1=-u; c2 = re/(2*l);

u = @(x) c1*tanh(c2*x);
ux = @(x) -c2/c1*(u(x).^2 - c1.^2);
uxx = @(x) -c2/c1*(2*u(x).*ux(x));
uxxx = @(x) -2*c2/c1*(ux(x).^2 + u(x).*uxx(x));


n = 129;
x0 = 0;
dL = 2;
time = linspace(x0-dL,x0+dL,n);

teq = (time(2:end)+time(1:end-1))/2; %this is acutally the cell center time
dteq = time(2:end) - time(1:end-1); %this is the time spacing centered at the cell center

dt0 = (2*dL)/(n-1);


% All nodes perturbed
pert = 0.2; %fraction of mesh spacing
rng(1);
time(2:end-1) = time(2:end-1)+dt0*pert*rand(1,n-2);

% One node perturbed
nn = floor(n/2);
% time(nn) = time(nn)+dt0*pert;

t = (time(2:end)+time(1:end-1))/2; %this is acutally the cell center time
dt = time(2:end) - time(1:end-1); %this is the time spacing centered at the cell center

% Continuous and Discrete Operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = @(x) ux(x);
Lh = @(x1, x2, x3, dt1, dt2, dt3) (u(x1)-u(x2))./dt1 + (u(x1)-u(x2))./(dt1+dt2) - (dt2./dt1).*(u(x2)-u(x3))./(dt2+dt3);
Lheq = @(u1, u2, u3, dt) (3/2*u1 - 2*u2 + 1/2*u3)./dt;

% Lh_new = @(x1, x2, x3, u1, u2, u3) (u1.*(x2).^2 - u1.*(x3).^2 + u2.*(x3).^2 - u3.*(x2).^2)./( (x2).^2.*(x3) - (x2).*(x3).^2);%function of x
Lh_new = @(dx1, dx2, dx3, u1, u2, u3) -(u1-u3).*(2*(dx1+dx2))./( (dx2+dx3).*(dx1+2*dx2+dx3) ) + (u1-u2)*2.*(dx1+2*dx2+dx3)./( (dx1+dx2).*(dx2+dx3) );
% TE estimation functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact TE

TEex = @(x1, x2, x3, dt1, dt2, dt3) (u(x1)-u(x2))./dt1 + (u(x1)-u(x2))./(dt1+dt2) - (dt2./dt1).*(u(x2)-u(x3))./(dt2+dt3) - L(x1);

TEex_new = @(x1,x2,x3, dt1, dt2, dt3) Lh_new(dt1, dt2, dt3, u(x1), u(x2), u(x3)) - L(x1);
TE_new = @(x, dt1, dt2, dt3, uxxx) -uxxx.*(   ((dt1+dt2)/2) .* (dt1/2 + dt2 + dt3/2) )*1/6;

% TEex = @(x1,x2,x3, dt1, dt2, dt3) Lh(x1, x2, x3, dt1, dt2, dt3) - L(x1);

% Direct TE calculation
TE1 = @(x, ux, dt1, dt2, dt3) -ux.*( 0 );
TE2 = @(x, uxx, dt1, dt2, dt3) -uxx.*( dt1.^2*12 + dt1.*dt2*6 - dt2.^2*12 - dt2.*dt3*6 )./(48*dt1);
TE3 = @(x, uxxx, dt1, dt2, dt3) -uxxx.*( -dt1.^3*2 - dt1.^2.*dt2*2 + dt1.*dt2.^2*5 + dt1.*dt2.*dt3*3 + dt2.^3*6 + dt2.^2.*dt3*5 + dt2.*dt3.^2 )./(48*dt1);
TE = @(x, dt1, dt2, dt3, ux, uxx, uxxx) TE1(x, ux, dt1, dt2, dt3) + TE2(x, uxx, dt1, dt2, dt3) + TE3(x, uxxx, dt1, dt2, dt3);

% Equal spacing Direct TE calculation
TEeq = @(x,dt) -uxxx(x).*dt.^2*1/3;

% TE estimation methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Higher order Reconstruction
p=@(x,coef) coef(1) + coef(2)*x + coef(3)*x.^2 + coef(3)*x.^3;
A3fit = @(x1,x2,x3,x4) [   1, 0, 0,0 
                           1, -x2, x2.^2, -x2.^3 
                           1, -x3, x3.^2, -x3.^3 
                           1, -x4, x4.^2, -x4.^3 ];
A3inv = @(x1,x2,x3,x4)  [                                  1,                                               0,                                               0,                                               0
                         (x2*x3 + x2*x4 + x3*x4)/(x2*x3*x4),   (x3*x4)/(x2^2*x3 + x2^2*x4 - x2^3 - x2*x3*x4),   (x2*x4)/(x2*x3^2 + x3^2*x4 - x3^3 - x2*x3*x4),   (x2*x3)/(x2*x4^2 + x3*x4^2 - x4^3 - x2*x3*x4)
                                  (x2 + x3 + x4)/(x2*x3*x4), (x3 + x4)/(x2^2*x3 + x2^2*x4 - x2^3 - x2*x3*x4), (x2 + x4)/(x2*x3^2 + x3^2*x4 - x3^3 - x2*x3*x4), (x2 + x3)/(x2*x4^2 + x3*x4^2 - x4^3 - x2*x3*x4)
                                               1/(x2*x3*x4),         1/(x2^2*x3 + x2^2*x4 - x2^3 - x2*x3*x4),         1/(x2*x3^2 + x3^2*x4 - x3^3 - x2*x3*x4),         1/(x2*x4^2 + x3*x4^2 - x4^3 - x2*x3*x4)];
                          

i = 5:n-2;

for j=i
    coef = A3inv(0, t(j)-t(j-1), t(j)-t(j-2), t(j)-t(j-3))*[u(t(j)); u(t(j-1)); u(t(j-2)); u(t(j-3))];
    
    TE_est(j-4) = TE(t(j), dt(j), dt(j-1), dt(j-2), coef(2), 2*coef(3), 6*coef(4) );
    TEdir_est(j-4) =  Lh( t(j), t(j-1), t(j-2), dt(j), dt(j-1), dt(j-2) ) - coef(2);
    
   
    %smooth grid

    dtmean = dt(j);
    
    TEeq_estA(j-4) = -2*coef(4)*dtmean^2;
    
    phi_eq=A3fit(0, dtmean, 2*dtmean,3*dtmean)*coef;
    TEsmooth_estA(j-4) = Lheq(phi_eq(1), phi_eq(2), phi_eq(3), dtmean)-coef(2);

    dtmean = mean(dt(j-3:j));
    TEeq_estB(j-4) = -2*coef(4)*dtmean^2;
    
    phi_eq=A3fit(0, dtmean, 2*dtmean,3*dtmean)*coef;
    TEsmooth_estB(j-4) = Lheq(phi_eq(1), phi_eq(2), phi_eq(3), dtmean)-coef(2);
    
    
    TE_new_est(j-4) = TE_new(t(i), dt(j), dt(j-1), dt(j-2), 6*coef(4));
    
    %Effect of aveage dt for M2b
    for jj=1:5
        N = jj-1;
        dtmean = mean(dt(j-N:j));
        TEeq_est_meant(j-4,jj) = -2*coef(4)*dtmean^2;
    end
    

end

%% Plots


figure(1)
% subplot(1,2,1)
plot(teq(i), TEex(teq(i), teq(i-1), teq(i-2), dteq(i), dteq(i-1), dteq(i-2)) ,'k-*',...
     t(i),   TEex(t(i), t(i-1), t(i-2), dt(i), dt(i-1), dt(i-2)) ,'k-v',...
     t(i),   TE(t(i), dt(i), dt(i-1), dt(i-2), ux(t(i)), uxx(t(i)), uxxx(t(i))), 'r-o',...
     t(i),   TEeq(t(i),dt(i)), 'b-*')
legend('TE exact (dt=constant)', 'TE exact', 'TE(u) analytic','TE(u) analytic (dt=constant)')
j1=i(1);j2=i(end); terr = TEex(teq(i), teq(i-1), teq(i-2), dteq(i), dteq(i-1), dteq(i-2));
axis([x0-dL,x0+dL, 0.5*min(terr), 1.5*max(terr)]);
title(['Truncation Error (Perturbed max ',num2str(pert),'*dt)'])

figure(2)
plot(teq(i), TEex(teq(i), teq(i-1), teq(i-2), dteq(i), dteq(i-1), dteq(i-2)) ,'k-*',...
     t(i),   TEex(t(i), t(i-1), t(i-2), dt(i), dt(i-1), dt(i-2)) ,'k-v',...
     t(i),   TE_est, 'r-o',...
     t(i),   TEdir_est, 'r-*',...
     t(i),   TEeq_estA,'b-o',...
     t(i),   TEsmooth_estA,'b-*',...
     t(i),   TEeq_estB,'g-o',...
     t(i),   TEsmooth_estB,'g-*')
legend('TE exact (dt=constant)', 'TE exact', 'M1a','M1b', 'M2a', 'M3a','M2b (N=2)','M3b (N=2)')
j1=i(1);j2=i(end); terr = TEex(teq(i), teq(i-1), teq(i-2), dteq(i), dteq(i-1), dteq(i-2));
axis([x0-dL,x0+dL, 0.5*min(terr), 1.5*max(terr)]);
title(['Truncation Error Estimates (Perturbed max ',num2str(pert),'*dt)'])

figure(3)
plot(teq(i), TEex(teq(i), teq(i-1), teq(i-2), dteq(i), dteq(i-1), dteq(i-2)) ,'k-*',...
     t(i),   TEeq_est_meant(:,1), 'r-o',...
     t(i),   TEeq_est_meant(:,2), 'g-o',...
     t(i),   TEeq_est_meant(:,3), 'b-o',...
     t(i),   TEeq_est_meant(:,4), 'r-*',...
     t(i),   TEeq_est_meant(:,5), 'g-*')
legend('TE exact (dt=constant)', 'N=0','N=1','N=2','N=3','N=4')
j1=i(1);j2=i(end); terr = TEex(teq(i), teq(i-1), teq(i-2), dteq(i), dteq(i-1), dteq(i-2));
% axis([x0-dL,x0+dL, 0.5*min(terr), 1.5*max(terr)]);
title(['Truncation Error Estimates (M2b) with Variable Smoothing (Perturbed max ',num2str(pert),'*dt)'])


figure(4)
plot(teq(i), TEex(teq(i), teq(i-1), teq(i-2), dteq(i), dteq(i-1), dteq(i-2)) ,'k-*',...
     t(i),   TEex(t(i), t(i-1), t(i-2), dt(i), dt(i-1), dt(i-2)) ,'k-v',...
     t(i),   TEex_new(t(i),t(i-1),t(i-2), dt(i), dt(i-1), dt(i-2)),'r-*',...
     t(i),   TE_new_est,'g-o',...
     t(i),   TEeq_estB,'b-*')
 
legend('TE exact (dt=constant)', 'TE exact', 'TE exact (new)', 'TE estimate M1a (new)', 'TE estimate M2b')
j1=i(1);j2=i(end); terr = TEex(teq(i), teq(i-1), teq(i-2), dteq(i), dteq(i-1), dteq(i-2));
axis([x0-dL,x0+dL, 0.5*min(terr), 1.5*max(terr)]);
title(['Truncation Error (M2b) for a new derivative calculation (Perturbed max ',num2str(pert),'*dt)'])

