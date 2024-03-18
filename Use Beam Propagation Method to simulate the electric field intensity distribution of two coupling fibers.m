clear all
close all
tic
wavelength=1.55e-6;
k=2*pi/wavelength;
nf=1.45;
nc=1;
nr=1.45;
col = 1001;
x=-5e-6:(10e-6)/(col-1):5e-6;
z=0:(1e-3)/(col-1):1e-3;
% x=linspace(-5e-6,5e-6,col);
% z=linspace(0,1e-3,col);

[zz,xx]=meshgrid(z,x);
clad=nc*ones(col);
fiber=find(0.75e-6<=xx & xx<=1.25e-6 & 0.01e-6<=zz & zz<=1e-3); 
clad(fiber)=nf; 
fiber=find(-0.75e-6>=xx & xx>=-1.25e-6 & 0<=zz & zz<=1e-3); 
clad(fiber)=nf; 

subplot(1,2,1);
mesh(zz,xx,clad);
xlabel('z, meter');
ylabel('x,meter');
title('Close Proximity Optical Fibers');
ylim([-5e-6 5e-6]);
view(0,90);

number1=find(zz==0);
delx=k*(10/(col-1))*10^-6;
delz=k*(1000/(col-1))*10^-6;
b=(1/(delx^2));
Aa=(clad(number1)).^2-2*b;
temp=repmat([b],1,col-1);
Ab1=[0 temp];
Ab2=[temp 0];
Array_a=sparse(1:col,1:col,Aa,col,col);
Array_b1=sparse(2:col,1:col-1,Ab1(2:size(Ab1,2)),col,col);
Array_b2=sparse(1:col-1,2:col,Ab2(1:size(Ab2,2)-1),col,col);
A_result=Array_a+Array_b1+Array_b2;
Unit=sparse(1:col,1:col,1,col,col);
[evector,evalue,iresult]=sptarn(A_result,Unit,nc^2,nf^2); %evector is bpm input
size = col-1;
for n=1:size
    phi(:,1) = evector;
    a_p=(((clad(:,n)+clad(:,n+1))/2).^2-nr.^2-2.*b)+(i.*4.*nr.*(1/delz));
    a_m=(((clad(:,n)+clad(:,n+1))/2).^2-nr.^2-2.*b)-(i.*4.*nr.*(1/delz)); 
    A=sparse(1:col,1:col,a_p,col,col);
    B=sparse(1:col,1:col,a_m,col,col);
    AA=A+Array_b1+Array_b2; 
    BB=B+Array_b1+Array_b2;
    phi(:,n+1)=-inv(BB)*AA*phi(:,n);
end
E=abs(phi);
subplot(1,2,2);
mesh(zz,xx,E);
view(0,90);
xlabel('z, meter');
ylabel('x, meter');
ylim([-5e-6 5e-6]);
title('Mode coupling between 2 fibers');
toc