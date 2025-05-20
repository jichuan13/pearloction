function og=OriginalBoundaryValue(beta0,n,loc)
h=2*pi/n;
beta1=beta0(:,[1:5]);x0=loc(1,1);y0=loc(1,2);
rmax=@(theta) beta1(1,1)+beta1(1,2)*cos(theta)+beta1(1,3)*sin(theta)+...
   beta1(1,4)*cos(2*theta)+beta1(1,5)*sin(2*theta); %polar coordination:
for j=1:n+1
    theta1=(j-1)*h;
    z1=cos(theta1);z2=sin(theta1);
    c=[z1,z2];
    g(j,1)=quad2d(@(theta,r) func3(theta,r,c,x0,y0),0,2*pi,0,rmax);
end
g1=g;

beta2=beta0(:,[6:10]);x0=loc(1,3);y0=loc(1,4);
rmax=@(theta) beta2(1,1)+beta2(1,2)*cos(theta)+beta2(1,3)*sin(theta)+...
   beta2(1,4)*cos(2*theta)+beta2(1,5)*sin(2*theta); %polar coordination:
for j=1:n+1
    theta1=(j-1)*h;
    z1=cos(theta1);z2=sin(theta1);
    c=[z1,z2];
    g(j,1)=quad2d(@(theta,r) func3(theta,r,c,x0,y0),0,2*pi,0,rmax);
end
g2=g;

beta3=beta0(:,[11:15]);x0=loc(1,5);y0=loc(1,6);
rmax=@(theta) beta3(1,1)+beta3(1,2)*cos(theta)+beta3(1,3)*sin(theta)+...
   beta3(1,4)*cos(2*theta)+beta3(1,5)*sin(2*theta); %polar coordination:
for j=1:n+1
    theta1=(j-1)*h;
    z1=cos(theta1);z2=sin(theta1);
    c=[z1,z2];
    g(j,1)=quad2d(@(theta,r) func3(theta,r,c,x0,y0),0,2*pi,0,rmax);
end
g3=g;
og=g1+g2+g3;