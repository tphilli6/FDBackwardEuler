%% 2nd order Backward Euler TE analysis
clear all

syms nu dx u ux uxx uxxx uxxxx uxxxxx uxxxxxx dx1 dx2 dx3 dx4 x1 x2 x3 x4 u1 u2 u3 u4 reals %ut uxt uxxt uxxxt uxxxxt uxxxxxt uxxxxxxt c1 c2 c3 f fx fxx fxxx fxxxx fxxxxx nu ubar reals
clc


u_m1 = u  - ux*dx   +  1/2*uxx*dx^2  -  1/6*uxxx*dx^3  +    1/24*uxxxx*dx^4  +  1/factorial(5)*uxxxxx*dx^5;
u_m2 = u  - ux*dx*2 +  4/2*uxx*dx^2  -  8/6*uxxx*dx^3  +    16/24*uxxxx*dx^4  +  32/factorial(5)*uxxxxx*dx^5;
u_m3 = u  - ux*dx*3 +  9/2*uxx*dx^2  -  27/6*uxxx*dx^3  +    81/24*uxxxx*dx^4  +  243/factorial(5)*uxxxxx*dx^5;

u_m1t = u  - ux*( (dx1+dx2)/2 )   +  1/2*uxx*( (dx1+dx2)/2 )^2  -  1/6*uxxx*( (dx1+dx2)/2 )^3 ;% +    1/24*uxxxx*(dx1)^4;%  +  1/factorial(5)*uxxxxx*(dx1)^5;
u_m2t = u  - ux*( dx1/2+dx2+dx3/2 )   +  1/2*uxx*(dx1/2+dx2+dx3/2)^2  -  1/6*uxxx*(dx1/2+dx2+dx3/2)^3 ;%+     1/24*uxxxx*(dx1+dx2)^4;%  +  1/factorial(5)*uxxxxx*(dx1+dx2)^5;
% u_m3t = u  - ux*(dx1/2+dx2+dx3+dx4/2)   +  1/2*uxx*(dx1/2+dx2++dx3+dx4/2)^2  -  1/6*uxxx*(dx1/2+dx2++dx3+dx4/2)^3;%  +    1/24*uxxxx*(dx1+dx2+dx3)^4;%  +  1/factorial(5)*uxxxxx*(dx1+dx2+dx3)^5;

% u_m1t = u  - ux*( x2 )   +  1/2*uxx*( x2 )^2  -  1/6*uxxx*( x2 )^3 ;% +    1/24*uxxxx*(dx1)^4;%  +  1/factorial(5)*uxxxxx*(dx1)^5;
% u_m2t = u  - ux*( x3 )   +  1/2*uxx*( x3 )^2  -  1/6*uxxx*( x3 )^3 ;%+     1/24*uxxxx*(dx1+dx2)^4;%  +  1/factorial(5)*uxxxxx*(dx1+dx2)^5;
% u_m3t = u  - ux*(dx1/2+dx2+dx3+dx4/2)   +  1/2*uxx*(dx1/2+dx2++dx3+dx4/2)^2  -  1/6*uxxx*(dx1/2+dx2++dx3+dx4/2)^3;%  +    1/24*uxxxx*(dx1+dx2+dx3)^4;%  +  1/factorial(5)*uxxxxx*(dx1+dx2+dx3)^5;


A = [1, -1, 1/2, -1/6, 1/24
     1, -2, 4/2, -8/6, 16/24
     1, -3, 9/2, -27/6, 81/24];
At=A(1:3,2:4)';
b = [0, 0, 1]';
cc = At\b;


dfdx = (3*u - 4*u_m1t + u_m2t)/(2*dx);
pretty(simplify(expand(dfdx)));
fnph = u + dx1/(dx1+dx2)*(u - u_m1t);
fnmh = u_m1t + dx2/(dx2+dx3)*(u_m1t - u_m2t);


% pretty(expand(fnph))
% pretty(expand(fnmh))
dfdx2 = (fnph-fnmh)/dx1;
dfdx2eq = subs(dfdx2,{dx1,dx2,dx3},{dx,dx,dx});
pretty(expand(dfdx2-ux));

pretty(expand( (u - 3*u_m1 + 3*u_m2 - u_m3)/dx^3 ))




dfdx2_nonequal = (u-u_m1t)/dx1 + (u-u_m1t)/(dx1+dx2) - dx2/dx1 * (u_m1t - u_m2t)/(dx2+dx3);
pretty(expand(dfdx2_nonequal))
pretty(expand(subs(dfdx2_nonequal,{dx1, dx2, dx3}, {dx, dx, dx}) ));


% fx1 = @(x1, x2, u1, u2, u3) -(u1.*x1.^2 - u1.*x2.^2 + u2.*x2.^2 - u3.*x1.^2)./(- x1.^2.*x2 + x1.*x2.^2)
dfdx2_nonequal = (u-u_m1t)/dx1 + (u-u_m1t)/(dx1+dx2) - dx2/dx1 * (u_m1t - u_m2t)/(dx2+dx3);
dfdx3_curvefit = -(u*x2^2 - u*x3^2 + u_m1t*x3^2 - u_m2t*x2^2)/(x2*x3*(x2-x3));
pretty(simplify(expand(dfdx3_curvefit)));
% fnmh



re=32; l=8; u=2; nu = l*u/re;
c1=-u; c2 = re/(2*l);

u = @(x) c1*tanh(c2*x);
ux = @(x) -c2/c1*(u(x).^2 - c1.^2);
uxx = @(x) -c2/c1*(2*u(x).*ux(x));
uxxx = @(x) -2*c2/c1*(ux(x).^2 + u(x).*uxx(x));

% u = @(x) c1*x.^1;
% ux = @(x) 1.*c1.*ones(size(x));
% uxx = @(x) 2.*c1.*zeros(size(x));
% uxxx = @(x) 0;



A1 = [ 1, 0, 0
         1, x1, x1.^2 
         1, x2, x2.^2 ];
A1 = [ 1, 0, 0
       1, -(dx1+dx2)/2, ((dx1+dx2)/2).^2
       1, -(dx1/2+dx2+dx3/2), ((dx1/2+dx2+dx3/2)).^2 ];
A1 = [ 1, 0, 0
     1, x2, (x2).^2
     1, x3, (x3).^2];
A1u1 = A1; A1u1(:,2) = [u1; u2; u3];

fx1a = det(A1u1)/det(A1);

fx1 = @(dx1, dx2, dx3, u1, u2, u3) (u1.*(dx1/2+dx2/2).^2 - u1.*(dx1/2+dx2+dx3/2).^2 + u2.*(dx1/2+dx2+dx3/2).^2 - u3.*(dx1/2+dx2/2).^2)./( (dx1/2+dx2/2).^2.*(dx1/2+dx2+dx3/2) - (dx1/2+dx2/2).*(dx1/2+dx2+dx3/2).^2);%function of x


A3 = [ 1, 0, 0,0 
         1, x2-x1, (x2-x1).^2, (x2-x1).^3 
         1, x3-x1, (x3-x1).^2, (x3-x1).^3 
         1, x4-x1, (x4-x1).^2, (x4-x1).^3 ];
A3 = [ 1, 0, 0,0 
       1, -(dx1+dx2)/2, ((dx1+dx2)/2).^2, -((dx1+dx2)/2).^3 
       1, -(dx1/2+dx2+dx3/2), ((dx1/2+dx2+dx3/2)).^2, -((dx1/2+dx2+dx3/2)).^3 
       1, -(dx1/2+dx2+dx3+dx4/2), ((dx1/2+dx2+dx3+dx4/2)).^2, -((dx1/2+dx2+dx3+dx4/2)).^3 ];
A3 = [ 1, 0, 0,0 
       1, -x2, x2.^2, -x2.^3 
       1, -x3, x3.^2, -x3.^3 
       1, -x4, x4.^2, -x4.^3 ];

A3lsq = @(x2,x3,x4,x5) [ 1, 0, 0, 0
                         1, -x2, x2.^2, -x2.^3 
                         1, -x3, x3.^2, -x3.^3 
                         1, -x4, x4.^2, -x4.^3
                         1, -x5, x5.^2, -x5.^3 ];


A3fit = @(x1,x2,x3,x4) [   1, 0, 0,0 
                           1, -x2, x2.^2, -x2.^3 
                           1, -x3, x3.^2, -x3.^3 
                           1, -x4, x4.^2, -x4.^3 ];
A3inv = @(x1,x2,x3,x4)  [                                  1,                                               0,                                               0,                                               0
                         (x2*x3 + x2*x4 + x3*x4)/(x2*x3*x4),   (x3*x4)/(x2^2*x3 + x2^2*x4 - x2^3 - x2*x3*x4),   (x2*x4)/(x2*x3^2 + x3^2*x4 - x3^3 - x2*x3*x4),   (x2*x3)/(x2*x4^2 + x3*x4^2 - x4^3 - x2*x3*x4)
                                  (x2 + x3 + x4)/(x2*x3*x4), (x3 + x4)/(x2^2*x3 + x2^2*x4 - x2^3 - x2*x3*x4), (x2 + x4)/(x2*x3^2 + x3^2*x4 - x3^3 - x2*x3*x4), (x2 + x3)/(x2*x4^2 + x3*x4^2 - x4^3 - x2*x3*x4)
                                               1/(x2*x3*x4),         1/(x2^2*x3 + x2^2*x4 - x2^3 - x2*x3*x4),         1/(x2*x3^2 + x3^2*x4 - x3^3 - x2*x3*x4),         1/(x2*x4^2 + x3*x4^2 - x4^3 - x2*x3*x4)];

A3u3 = A3; A3u3(:,4) = [u1; u2; u3; u4];


A1u3 = A3; A1u3(:,2) = [u1; u2; u3; u4];
fx3a = det(A1u3)/det(A3);
fx3 = @(x1,x2,x3,x4,u1,u2,u3,u4) (u1*x2^2*x3^3)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) - (u1*x2^3*x3^2)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) - (u1*x2^2*x4^3)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) + (u1*x2^3*x4^2)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) + (u1*x3^2*x4^3)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) - (u1*x3^3*x4^2)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) - (u2*x3^2*x4^3)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) + (u2*x3^3*x4^2)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) + (u3*x2^2*x4^3)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) - (u3*x2^3*x4^2)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) - (u4*x2^2*x3^3)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3) + (u4*x2^3*x3^2)/(- x2^3*x3^2*x4 + x2^3*x3*x4^2 + x2^2*x3^3*x4 - x2^2*x3*x4^3 - x2*x3^3*x4^2 + x2*x3^2*x4^3);


fxxx3a = 6*det(A3u3)/det(A3);
fxxx3 = @(x1,x2,x3,x4,u1,u2,u3,u4) -(6*u1*x2*x3^2 - 6*u1*x2^2*x3 - 6*u1*x2*x4^2 + 6*u1*x2^2*x4 + 6*u1*x3*x4^2 - 6*u1*x3^2*x4 - 6*u2*x3*x4^2 + 6*u2*x3^2*x4 + 6*u3*x2*x4^2 - 6*u3*x2^2*x4 - 6*u4*x2*x3^2 + 6*u4*x2^2*x3)/(x2*x3*x4*(x2 - x3)*(x2 - x4)*(x3 - x4));




fxxx1 = @(u1, u2, u3, u4, dx1, dx2, dx3, dx4) ((3*dx3^3*u1)/2 + (3*dx2^3*u3)/2 - (3*dx3^3*u2)/2 - (3*dx2^3*u4)/2 + (9*dx1*dx2^2*u3)/4 - (9*dx1*dx3^2*u2)/4 + (9*dx2*dx3^2*u1)/4 + (3*dx1^2*dx2*u3)/4 - (3*dx1^2*dx3*u2)/4 + (3*dx2^2*dx3*u1)/4 - (9*dx1*dx2^2*u4)/4 + 3*dx1*dx3^2*u3 - (3*dx1*dx4^2*u2)/4 - (9*dx2*dx3^2*u2)/2 + (3*dx2*dx4^2*u1)/4 - (3*dx1^2*dx2*u4)/4 + (3*dx1^2*dx3*u3)/2 - (3*dx1^2*dx4*u2)/4 - 3*dx2^2*dx3*u2 + (3*dx2^2*dx4*u1)/4 - (3*dx1*dx3^2*u4)/4 + (3*dx1*dx4^2*u3)/4 + 3*dx2*dx3^2*u3 - (3*dx2*dx4^2*u2)/2 + (3*dx3*dx4^2*u1)/4 - (3*dx1^2*dx3*u4)/4 + (3*dx1^2*dx4*u3)/4 + (9*dx2^2*dx3*u3)/2 - 3*dx2^2*dx4*u2 + (9*dx3^2*dx4*u1)/4 - (3*dx2*dx3^2*u4)/4 + (3*dx2*dx4^2*u3)/4 - (3*dx3*dx4^2*u2)/4 - (9*dx2^2*dx3*u4)/4 + (9*dx2^2*dx4*u3)/4 - (9*dx3^2*dx4*u2)/4 - 3*dx1*dx2*dx3*u2 + 6*dx1*dx2*dx3*u3 - 3*dx1*dx2*dx4*u2 - 3*dx1*dx2*dx3*u4 + 3*dx1*dx2*dx4*u3 - 3*dx1*dx3*dx4*u2 + 3*dx2*dx3*dx4*u1 + 3*dx1*dx3*dx4*u3 - 6*dx2*dx3*dx4*u2 + 3*dx2*dx3*dx4*u3)/((dx1^3*dx2^2*dx3)/64 + (dx1^3*dx2^2*dx4)/64 + (3*dx1^3*dx2*dx3^2)/64 + (dx1^3*dx2*dx3*dx4)/16 + (dx1^3*dx2*dx4^2)/64 + (dx1^3*dx3^3)/32 + (3*dx1^3*dx3^2*dx4)/64 + (dx1^3*dx3*dx4^2)/64 + (5*dx1^2*dx2^3*dx3)/64 + (5*dx1^2*dx2^3*dx4)/64 + (9*dx1^2*dx2^2*dx3^2)/32 + (3*dx1^2*dx2^2*dx3*dx4)/8 + (3*dx1^2*dx2^2*dx4^2)/32 + (19*dx1^2*dx2*dx3^3)/64 + (15*dx1^2*dx2*dx3^2*dx4)/32 + (3*dx1^2*dx2*dx3*dx4^2)/16 + (dx1^2*dx2*dx4^3)/64 + (3*dx1^2*dx3^4)/32 + (11*dx1^2*dx3^3*dx4)/64 + (3*dx1^2*dx3^2*dx4^2)/32 + (dx1^2*dx3*dx4^3)/64 + (dx1*dx2^4*dx3)/8 + (dx1*dx2^4*dx4)/8 + (33*dx1*dx2^3*dx3^2)/64 + (11*dx1*dx2^3*dx3*dx4)/16 + (11*dx1*dx2^3*dx4^2)/64 + (45*dx1*dx2^2*dx3^3)/64 + (9*dx1*dx2^2*dx3^2*dx4)/8 + (15*dx1*dx2^2*dx3*dx4^2)/32 + (3*dx1*dx2^2*dx4^3)/64 + (3*dx1*dx2*dx3^4)/8 + (11*dx1*dx2*dx3^3*dx4)/16 + (3*dx1*dx2*dx3^2*dx4^2)/8 + (dx1*dx2*dx3*dx4^3)/16 + (dx1*dx3^5)/16 + (dx1*dx3^4*dx4)/8 + (5*dx1*dx3^3*dx4^2)/64 + (dx1*dx3^2*dx4^3)/64 + (dx2^5*dx3)/16 + (dx2^5*dx4)/16 + (9*dx2^4*dx3^2)/32 + (3*dx2^4*dx3*dx4)/8 + (3*dx2^4*dx4^2)/32 + (7*dx2^3*dx3^3)/16 + (45*dx2^3*dx3^2*dx4)/64 + (19*dx2^3*dx3*dx4^2)/64 + (dx2^3*dx4^3)/32 + (9*dx2^2*dx3^4)/32 + (33*dx2^2*dx3^3*dx4)/64 + (9*dx2^2*dx3^2*dx4^2)/32 + (3*dx2^2*dx3*dx4^3)/64 + (dx2*dx3^5)/16 + (dx2*dx3^4*dx4)/8 + (5*dx2*dx3^3*dx4^2)/64 + (dx2*dx3^2*dx4^3)/64);
fxxx = @(x1, x2, x3, x4, dx1, dx2, dx3, dx4) fxxx1(u(x1), u(x2), u(x3), u(x4), dx1, dx2, dx3, dx4);
 

% fxxx = @(x1, x2, x3, x4, dx1, dx2, dx3, dx4) -(96.*dx3.^3.*u(x1) + 96.*dx2.^3.*u(x3) - 96.*dx3.^3.*u(x2) - 96.*dx2.^3.*u(x4) + 144.*dx1.*dx2.^2.*u(x3)...
%                                       - 144.*dx1.*dx3.^2.*u(x2) + 144.*dx2.*dx3.^2.*u(x1) + 48.*dx1.^2.*dx2.*u(x3) - 48.*dx1.^2.*dx3.*u(x2)...
%                                       + 48.*dx2.^2.*dx3.*u(x1) - 144.*dx1.*dx2.^2.*u(x4) + 192.*dx1.*dx3.^2.*u(x3) - 48.*dx1.*dx4.^2.*u(x2)...
%                                       - 288.*dx2.*dx3.^2.*u(x2) + 48.*dx2.*dx4.^2.*u(x1) - 48.*dx1.^2.*dx2.*u(x4) + 96.*dx1.^2.*dx3.*u(x3)...
%                                       - 48.*dx1.^2.*dx4.*u(x2) - 192.*dx2.^2.*dx3.*u(x2) + 48.*dx2.^2.*dx4.*u(x1) - 48.*dx1.*dx3.^2.*u(x4)...
%                                       + 48.*dx1.*dx4.^2.*u(x3) + 192.*dx2.*dx3.^2.*u(x3) - 96.*dx2.*dx4.^2.*u(x2) + 48.*dx3.*dx4.^2.*u(x1)...
%                                       - 48.*dx1.^2.*dx3.*u(x4) + 48.*dx1.^2.*dx4.*u(x3) + 288.*dx2.^2.*dx3.*u(x3) - 192.*dx2.^2.*dx4.*u(x2)...
%                                       + 144.*dx3.^2.*dx4.*u(x1) - 48.*dx2.*dx3.^2.*u(x4) + 48.*dx2.*dx4.^2.*u(x3) - 48.*dx3.*dx4.^2.*u(x2)...
%                                       - 144.*dx2.^2.*dx3.*u(x4) + 144.*dx2.^2.*dx4.*u(x3) - 144.*dx3.^2.*dx4.*u(x2) - 192.*dx1.*dx2.*dx3.*u(x2)...
%                                       + 384.*dx1.*dx2.*dx3.*u(x3) - 192.*dx1.*dx2.*dx4.*u(x2) - 192.*dx1.*dx2.*dx3.*u(x4) + 192.*dx1.*dx2.*dx4.*u(x3)...
%                                       - 192.*dx1.*dx3.*dx4.*u(x2) + 192.*dx2.*dx3.*dx4.*u(x1) + 192.*dx1.*dx3.*dx4.*u(x3) - 384.*dx2.*dx3.*dx4.*u(x2)...
%                                       + 192.*dx2.*dx3.*dx4.*u(x3))./((dx1 + dx2).*(dx2 + dx3).*(dx3 + dx4).*(dx1 + 2.*dx2 + dx3).*(dx2 + 2.*dx3...
%                                       + dx4).*(dx1 + 2.*dx2 + 2.*dx3 + dx4));
% %                                   
% fxxx3 = @(x1,x2,x3,x4)[0, 0, 0, 6]*( [   1, 0, 0,0 
%                                          1, x2-x1, (x2-x1).^2, (x2-x1).^3 
%                                          1, x3-x1, (x3-x1).^2, (x3-x1).^3 
%                                          1, x4-x1, (x4-x1).^2, (x4-x1).^3 ]\[u(x1); u(x2); u(x3); u(x4)]);
fxxx4 = @(x1,x2,x3,x4,dx1,dx2,dx3,dx4) [0, 0, 0, 6]*(   [  1, 0, 0,0 
                                                           1, -(dx1+dx2)/2, ((dx1+dx2)/2).^2, -((dx1+dx2)/2).^3 
                                                           1, -(dx1/2+dx2+dx3/2), ((dx1/2+dx2+dx3/2)).^2, -((dx1/2+dx2+dx3/2)).^3 
                                                           1, -(dx1/2+dx2+dx3+dx4/2), ((dx1/2+dx2+dx3+dx4/2)).^2, -((dx1/2+dx2+dx3+dx4/2)).^3 ]\[u(x1); u(x2); u(x3); u(x4)]);

                                                       
n = 21;
x0 = 1;
dL = 0.1;

xnodes = linspace(x0-dL,x0+dL,n);
xnodes0 = xnodes;
dxnodes = xnodes(2:end)-xnodes(1:end-1);
dx0 = (2*dL)/(n-1);
nn = floor(n/2);

% xnodes(2:end-1) = xnodes(2:end-1)+dx0*0.05*rand(1,n-2);

xnodes(nn) = xnodes(nn)+dx0*0.05;


x = (xnodes(2:end)+xnodes(1:end-1))/2;
dx = xnodes(2:end) - xnodes(1:end-1);



% xc = (x(2:end)+x(1:end-1))/2;
% dx = [dx0, xc(2:end)-xc(1:end-1), dx0];

L = @(x) ux(x);

% TE for non-uniform spacing
TE1 = @(x, dx1, dx2, dx3) -ux(x).*( 0 );
TE2 = @(x, dx1, dx2, dx3) -uxx(x).*( dx1.^2*12 + dx1.*dx2*6 - dx2.^2*12 - dx2.*dx3*6 )./(48*dx1);
TE3 = @(x, dx1, dx2, dx3) -uxxx(x).*( -dx1.^3*2 - dx1.^2.*dx2*2 + dx1.*dx2.^2*5 + dx1.*dx2.*dx3*3 + dx2.^3*6 + dx2.^2.*dx3*5 + dx2.*dx3.^2 )./(48*dx1);

% TE for non-uniform curvefit
% TE1 = @(x, dx1, dx2, dx3) ux(x).*0;
% TE2 = @(x, dx1, dx2, dx3) uxx(x).*0;
% TE3 = @(x, dx1, dx2, dx3) -uxxx(x).*(   ((dx1+dx2)/2) .* (dx1/2 + dx2 + dx3/2) )*1/6;

TEeq = @(x,dx) -uxxx(x).*dx^2*1/3;


TE = @(x, dx1, dx2, dx3) TE1(x, dx1, dx2, dx3) + TE2(x, dx1, dx2, dx3) + TE3(x, dx1, dx2, dx3);

TEest = @(x1, x2, x3, x4, dx1)- dx1.^2/3.*(u(x1) - 3*u(x2) + 3*u(x3) - u(x4))./dx1.^3;

Lh = @(x1, x2, x3, dx1, dx2, dx3) (u(x1)-u(x2))./dx1 + (u(x1)-u(x2))./(dx1+dx2) - (dx2./dx1).*(u(x2)-u(x3))./(dx2+dx3) - L(x1);
fx = @(dx1, dx2, dx3, u1, u2, u3) (u1-u2)./dx1 + (u1-u2)./(dx1+dx2) - (dx2./dx1).*(u2-u3)./(dx2+dx3);
Lh1 = @(x1, x2, x3, dx1, dx2, dx3) fx1( dx1, dx2, dx3, u(x1), u(x2), u(x3) ) - L(x1);

Lheq = @(x1, x2, x3, dx1, dx2, dx3) 1./dx1.*(3/2*u(x1) - 2*u(x2) + 1/2*u(x3) ) - L(x1);
% TE = @(x1, dx1, dx2, dx3) -1/3*dx1.^2.*uxxx(x1);


i = 5:n-1;
TEex = Lh(x(i), x(i-1), x(i-2), dx(i), dx(i-1), dx(i-2));
TEcalc = TE(x(i), dx(i), dx(i-1), dx(i-2) );
% TEest2 = TEest(x(i), x(i-1), x(i-2), x(i-3), dx0);
% TEest3 = -fxxx(x(i), x(i-1), x(i-2), x(i-3), dx(i), dx(i-1), dx(i-2), dx(i-3)).*dx(i).^2/3;
for j=i
    TEest3(j-4) = fxxx3(0, x(j-1)-x(j), x(j-2)-x(j), x(j-3)-x(j), u(x(j)), u(x(j-1)), u(x(j-2)), u(x(j-3)) ).*dx(j).^2/3;
    TEest1(j-4) = -fxxx(x(j), x(j-1), x(j-2), x(j-3), dx(j), dx(j-1), dx(j-2), dx(j-3) ).*dx(j).^2/3;
    
    % TEest 
    TEest4(j-4) = fx( dx(j), dx(j-1), dx(j-2), u(x(j)), u(x(j-1)), u(x(j-2)) )...
                + fx3(0, x(j-1)-x(j), x(j-2)-x(j), x(j-3)-x(j), u(x(j)), u(x(j-1)), u(x(j-2)), u(x(j-3)) );
    
    
    A3lhs = A3lsq( x(j-1)-x(j), x(j-2)-x(j), x(j-3)-x(j), x(j-4)-x(j) );
    b = flipud(u(x(j-4:j))' );
%     coef = (A3lhs'*A3lhs)\(A3lhs'*b);
    
    coef = A3inv(0, x(j-1)-x(j), x(j-2)-x(j), x(j-3)-x(j))*[u(x(j)); u(x(j-1)); u(x(j-2)); u(x(j-3))];
    dxbar = mean(dx(j-2:j));
    ufit = A3fit(0, dxbar, 2*dxbar, 3*dxbar)*coef;
    dfdx = (3/2*ufit(1) - 2*ufit(2) + 1/2*ufit(3) )./dxbar;
    TEeqErr(j-4) = -dfdx + coef(2);
    
%     fx = ux(x(j));
%     fxx = uxx(x(j));
%     fxxx = uxxx(x(j));
%     x1=x(j);
%     x2=x(j-1);
%     x3=x(j-2);
%     dx1 = dx(j);
%     dx2 = dx(j-1);
%     dx3 = dx(j-2);
%     TEest1(j-3) = (3*fxxx*x1^3 - 12*fxxx*x1^2*x2 + 3*fxxx*x1^2*x3 + 9*fxx*x1^2 + 12*fxxx*x1*x2^2 - 24*fxx*x1*x2 - 3*fxxx*x1*x3^2 + 6*fxx*x1*x3 + 18*fx*x1 - 4*fxxx*x2^3 + 12*fxx*x2^2 - 24*fx*x2 + fxxx*x3^3 - 3*fxx*x3^2 + 6*fx*x3)/(12*dx1) - fx;


end

% figure(1)
plot(xnodes0(i), Lh(xnodes0(i), xnodes0(i-1), xnodes0(i-2), dxnodes(i), dxnodes(i-1), dxnodes(i-2)) ,'k-o',...
     x(i), Lh(x(i), x(i-1), x(i-2), dx(i), dx(i-1), dx(i-2)) ,'k-v',...
     x(i), Lh1(x(i), x(i-1), x(i-2), dx(i), dx(i-1), dx(i-2)) ,'k-*',...
     x(i), TEcalc, 'r-*',...
     x(i), TEeqErr, 'b-*',...
     x(i), TEest3, 'r-o',...
     x(i), TEest4, 'g-^')

% figure(2)
% plot(xnodes0(i), Lh(xnodes0(i), xnodes0(i-1), xnodes0(i-2), dxnodes(i), dxnodes(i-1), dxnodes(i-2)) ,'k-o',...
%      x(i), Lh(x(i), x(i-1), x(i-2), dx(i), dx(i-1), dx(i-2)) ,'k-v',...
%      x(i), Lh1(x(i), x(i-1), x(i-2), dx(i), dx(i-1), dx(i-2)) ,'k-*',...
%      x(i), TEcalc, 'r-*',...
%      x(i), TEeqErr, 'b-*',...
%      x(i), TEest3, 'r-o',...
%      x(i), TEest4, 'g-^')

 legend('TE-equal','TE','TE-mod','TE(x)-mod','TEsmooth','TE(x)','Direct Error Calc')
% plot(x(i), Lh(x(i), x(i-1), x(i-2), dx(i)) ,'k-o',...
%      x(i), TE(x(i), dx(i), dx(i-1), dx(i-2) ), 'r-o' )
x1=x(i); x2=x(i-1); x3=x(i-2); x4=x(i-3);
