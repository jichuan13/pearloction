clc
clear
n=30;
g=BoundaryValue(n);
rand('seed',0)
g=(1+0.01*(2*rand(n+1,1)-1)).*g;
loc=[-0.7,0.4,0.02+0.02*(2*rand(1,1)-1)];%1%
p=length(loc)/3;
new=loc;
h=2*pi/n;
n1=3;n2=30;
kk=1;
tol=norm(g);og1=g;
a=0:2*pi/30:2*pi;
r=loc(3);
jj=1;
while (jj<=p)
x2=loc((jj-1)*2+1)+r*cos(a);y2=loc((jj-1)*2+2)+r*sin(a);
hold on
plot(x2,y2,'k-')
jj=jj+1;
end


while(kk<=20000&&tol>1e-6)
now=new;
og2=og1;
g11=pois(now,n,p);
og1=g11;
GG=g-g11;
dth=2*pi/n; th=[0:dth:2*pi];
jj=1;
while (jj<=p)
x0=now((jj-1)*2+1); y0=now((jj-1)*2+2); r0=now(2*p+jj);
y1=x0+r0.*cos(th);y2=y0+r0.*sin(th);
y1d=-r0.*sin(th);y2d=r0.*cos(th);
z=cos(th)'+1i*sin(th)';x=y1+1i*y2;  
for k=1:length(th); 
    G=(log(abs(z-x(k)))-log(abs(z./abs(z)-abs(z)*x(k))))/(2*pi);
    GGG1(:,k)=G;
end
z=(1+1.e-7)*z;
for k=1:length(th); 
    G=(log(abs(z-x(k)))-log(abs(z./abs(z)-abs(z)*x(k))))/(2*pi);
    GGG2(:,k)=G; 
end
F=(GGG2-GGG1)/1.e-7; 
g1(:,(jj-1)*2+1)=F*(y2d.*1-y1d.*0)'*dth;
g1(:,(jj-1)*2+2)=F*(y2d.*0-y1d.*1)'*dth;
g1(:,2*p+jj)=F*(y2d.*cos(th)-y1d.*sin(th))'*dth;
jj=jj+1;
end

J=g1;R=-GG;
%阻尼最小二乘法
t0=svd(inv(J'*J));
lambda=(1/max(t0)+min(t0))/2;
new=now-(inv(J'*J+lambda*eye(length(now)))*(J'*R))';
%高斯牛顿法
% new=now-(inv(J'*J)*(J'*R))';
%最速下降法
 % alpha=0.6;
 % new=now-alpha*((J'*R))';
 % 
tol=norm(og1-og2)
kk=kk+1
end
new

x=new;
domain
a=0:2*pi/100:2*pi;
hold on
jj=1;
while (jj<=p)
x2=x((jj-1)*2+1)+x(2*p+jj)*cos(a);y2=x((jj-1)*2+2)+x(2*p+jj)*sin(a);
plot(x2,y2,'b-')
jj=jj+1;
end
% new =
% 
%    -0.0073    0.0015    0.3057
