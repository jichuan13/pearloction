function g11=pois(now,n,p)
h=2*pi/30;n1=20;n2=100;
jj=1;sum1=0;
while (jj<=p)
x0=now((jj-1)*2+1); y0=now((jj-1)*2+2); r0=now(2*p+jj);
for j=1:n+1 
th=(j-1)*h; c=[cos(th) sin(th)];
h1=r0/n1;h2=2*pi/n2;
[r,theta]=meshgrid(h1/2:h1:r0-h1/2,h2/2:h2:2*pi-h2/2);
x1=x0+r.*cos(theta); x2=y0+r.*sin(theta);
g1(j,1)=sum(sum(1./(2.*pi).*(1-x1.^2-x2.^2)./((c(1,1)-x1).^2+(c(1,2)-x2).^2).*r,2).*h1).*h2;
end;
sum1=sum1+g1;
jj=jj+1;
end

g11=sum1;

