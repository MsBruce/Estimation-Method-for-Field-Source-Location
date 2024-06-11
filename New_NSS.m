function [ NSS ] = New_NSS(mag_grd,n,m,Xmin,Xmax,Ymin,Ymax,h)
%N-MNSS
% mag_grd: magnetic data
% n,m the row and col of the data
% Xmin Xmax Ymin Ymax, the range of the data
% h: upward height
[n0,n1,n2,n3,dx] = Calculate_m1m2_dx(n,Xmin,Xmax);% x
[m0,m1,m2,m3,dy] = Calculate_m1m2_dx(m,Ymin,Ymax); % y
TMI = zeros(m3,n3);
TMI(m1:m2,n1:n2) = mag_grd; %
TMI = expound_2D(TMI,m0,m1,m2,m3,n0,n1,n2,n3);%
[S,U,V] = get_s(m3,n3,dx,dy);%S,U,V


m = m2-m1+1; %row
n = n2-n1+1; %col

NSS = zeros(m,n); % NSS
FTMI = fft2(TMI);% The fourier transform of TMI
Sxx = zeros(m3,n3); 
Sxy = zeros(m3,n3);
Sxz = zeros(m3,n3);
Syy = zeros(m3,n3);
Syz = zeros(m3,n3);
Szz2 = zeros(m3,n3);
jie = 0;
for j = 1:m3
    for k = 1:n3
        if j==1 && k==1
            Sxx(j,k) = 0;
            Sxy(j,k) = 0;
            Sxz(j,k) = 0;
            Syy(j,k) = 0;
            Syz(j,k) = 0; Szz2(j,k) = 0;
        else
%         qt = complex(double(N0*S(j,k)),double(L0*U(j)+M0*V(k))); % 
        qt = 1;
        qt2 = 1/exp(S(j,k)*-h);
        qt2 = qt2/(S(j,k))^jie;
        
%         size(qt)
        Sxx(j,k) = -U(j)/(qt/U(j))*FTMI(j,k); % The fourier transform of Sxx
        Sxx(j,k) = Sxx(j,k)/qt2 ;%* exp(S(j,k)*h);
        
        Sxy(j,k) = -U(j)/(qt/V(k))*FTMI(j,k); % The fourier transform of Sxy
        Sxy(j,k) = Sxy(j,k)/qt2 ;% * exp(S(j,k)*h);
        
        Sxz(j,k) = complex( 0,U(j) )/(qt/S(j,k))*FTMI(j,k); % The fourier transform of Sxz
        Sxz(j,k) = Sxz(j,k)/qt2 ;% * exp(S(j,k)*h);
        
        Syy(j,k) = -V(k)/(qt/V(k))*FTMI(j,k); % The fourier transform of Syy
        Syy(j,k) = Syy(j,k)/qt2 ;% * exp(S(j,k)*h);
        
        Syz(j,k) = complex( 0,V(k) )/(qt/S(j,k))*FTMI(j,k); % The fourier transform of Syz
        Syz(j,k) = Syz(j,k)/qt2 ;% * exp(S(j,k)*h);
        end
    end
end
Sxx = real( ifft2(Sxx) );
Sxy = real( ifft2(Sxy) );
Sxz = real( ifft2(Sxz) );
Syy = real( ifft2(Syy) );
Syz = real( ifft2(Syz) );
for i = m1:m2
    for j = n1:n2
        MGT = [Sxx(i,j) Sxy(i,j) Sxz(i,j);
               Sxy(i,j) Syy(i,j) Syz(i,j);
               Sxz(i,j) Syz(i,j) -Sxx(i,j)-Syy(i,j)]; %  
        D = eig(MGT); % 
        D = sort(D,'descend'); % 
        NSS(i-m1+1,j-n1+1) = sqrt(-(D(2)*D(2))-D(1)*D(3)); % NSS
    end
end
end

function [m0,m1,m2,m3,dx] = Calculate_m1m2_dx(mpoint,Xmin,Xmax)
%
m0 = 1;
dx = (Xmax-Xmin)/(mpoint-1);
mtemp = int16(mpoint);
class(mtemp);
while mod(mtemp,2)==0 && mtemp~=0 %
    mtemp = mtemp/2;
end

if mtemp==1
   m=mpoint*2;
else
    mu=int16(fix(log(double(mpoint))/0.693147+0.05)+2);
         m=2.^mu;
end

m1 = 1+(m-mpoint)/2;
m2 = m1+mpoint-1;
m3 = m;
end


function g = expound_2D(g,m0,m1,m2,m3,n0,n1,n2,n3)

for i = n1:n2
%     g(m0,i) = (g(m1,i)+g(m2,i))/2;
    g(m0,i) = 0;
    g(m3,i) = g(m0,i);
end

for j = n1:n2
    for i = m0+1:m1-1
        g(i,j) = g(m0,j)+cos(double(pi/2.0*(m1-i)/(m1-m0)))*(g(m1,j)-g(m0,j));
    end
    for i = m2+1:m3-1
        g(i,j) = g(m2,j)+cos(double(pi/2.0*(m3-i)/(m3-m2)))*(g(m3,j)-g(m2,j));
    end
end

for i = m0:m3
%     g(i,n0) = (g(i,n1)+g(i,n2))/2;
    g(i,n0) = 0;
    g(i,n3) = g(i,n0);
end

for i = m0:m3
    for j = n0+1:n1-1
        g(i,j) = g(i,n0)+cos(double(pi/2.0*(n1-j)/(n1-n0)))*(g(i,n1)-g(i,n0));
    end
    for j = n2+1:n3-1
        g(i,j) = g(i,n2)+cos(double(pi/2.0*(n3-j)/(n3-n2)))*(g(i,n3)-g(i,n2));
    end
end

end

function [S,U,V] = get_s(m3,n3,dx,dy)
S = zeros(m3,n3);
U = zeros(m3,1);
V = zeros(n3,1);
mp2 = (m3/2);
np2 = (n3/2);
xlength = double(n3-1)*dx;
ylength = double(m3-1)*dy;
qx = 2.0*pi/ylength;
qy = 2.0*pi/xlength;

for i = 0:m3-1
    uk = double(i);
%       uk = mp2-i;
    if i >= mp2
        uk = double(i-m3);
%         uk = double(mp2-i);
    end
    U(i+1) = uk*qx;
end

for j = 0:n3-1
    vl = double(j);
%     vl = double(np2-j);
    if j >= np2
       vl = double(j-n3);
%          vl = double(np2-j);
    end
    V(j+1) = vl*qy;
end

for j = 1:n3
    for i = 1:m3
        ss = U(i)*U(i) + V(j)*V(j);
        S(i,j) = sqrt(ss);
    end
end

end

