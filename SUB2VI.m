close all;
clc;
clear all;
warning off
%% System1 description
M=21;
m=0.42;
RK=0.106;
D=0.44;
L=0.3;
Jw=0.0024;
Jd=0.3388;
Jp=0.63;
g=9.8;
a23=(-M^2*L^2*g)/(M*Jp+2*(Jp+M*L^2)*(m+(Jw/RK^2)));
a43=(M^2*g*L+2*M*g*L*(m+(Jw/RK^2)))/(M*Jp+2*(Jp+M*L^2)*(m+(Jw/RK^2)));
b21=(((Jp+M*L^2)/RK)+M*L)/(M*Jp+2*(Jp+M*L^2)*(m+(Jw/RK^2)));
b41=-((((RK+L)*M)/RK))-2*(m+(Jw/RK^2))/(M*Jp+2*(Jp+M*L^2)*(m+(Jw/RK^2)));
b61=(D/(2*RK))/(Jd+(D^2/(2*RK))*(m*RK+(Jw/RK)));

A1=[1  0.1;
      0    1]
B1=[0.02596;
     0.5192]



R=1;
Q=5*eye(2);
x=[0.1 0.1]';
K0=[0 0];
G=[Q [0;0];[0 0] R];

K=K0;
zbar(1:6,1:6)=0;d_target(1:6,1)=0;
 H=[0 0 0 ;0 0 0 ;0 0 0 ];
H1=[];H2=[];H3=[];

 for i=1:400
    a1=1;

    if i>=150
        a1=0;
    end

 noise1= a1*0.001*(sin(1*i)^2*cos(9*i)+sin(2*i)^2*cos(0.1*i)+ ...
     sin(-1.2*i)^2*cos(0.5*i)+sin(i)^5+sin(1.12*i)^2+cos(2.4*i)* ...
     sin(2.4*i)^3);
ns(i)=noise1;
    u=K*x(:,i)+noise1;  
    x(:,i+1)=A1*x(:,i)+B1*u;
    d_target(1,1)=d_target(2,1);
    d_target(2,1)=d_target(3,1);
    d_target(3,1)=d_target(4,1);
    d_target(4,1)=d_target(5,1);
    d_target(5,1)=d_target(6,1);
    d_target(6,1)=[x(:,i); u]'*G*[x(:,i); u]+...
    [x(:,i+1);K*x(:,i+1)]'*H*[x(:,i+1);K*x(:,i+1)];

    zbar(:,1)=zbar(:,2);
    zbar(:,2)=zbar(:,3);
    zbar(:,3)=zbar(:,4);
    zbar(:,4)=zbar(:,5);
    zbar(:,5)=zbar(:,6);
    zbar(:,6)=[x(1,i)^2;x(1,i)*x(2,i);x(1,i)*u;x(2,i)^2;x(2,i)*u;u^2];

    if mod(i,6)==0
       if i<=150
        m1=zbar*zbar';
        size(m1);
        q=zbar*d_target;
        rank(m1);
        vH=inv(m1)*q;
        H=[vH(1,1) vH(2,1)/2 vH(3,1)/2 ; vH(2,1)/2 vH(4,1) vH(5,1)/2 ;...
            vH(3,1)/2 vH(5,1)/2  vH(6,1)];   
        Huu=H(3,3);
        Hux=H(3,1:2);
        H3=[H3,Huu];
        H1=[H1,H(3,1)]
        H2=[H2,H(3,2)]
        K=-inv(Huu)*Hux; 
       
       end

    end
 end
 figure(1);
sgtitle('State trajectories')
subplot(2,1,1);
plot(x(1,:),'linewidth',1.2);
xlabel('Time step');
ylabel('yaw angle δ (rad)');
hold on; 
subplot(2,1,2);
plot(x(2,:),'LineWidth',1)
xlabel('Time step');
ylabel('yaw angular velocity ̇𝛿 (rad/s)');
hold on; 
figure

plot(H1,'-o','Color','b','LineWidth',1)
hold on;
plot(H2,'-X','Color','r','LineWidth',1)
hold on;

plot(H3,'-diamond','Color','y','LineWidth',1)
xlabel('Iterations');
ylabel('Parameter Estimates');
title("Subsystem II")
hold on;
legend('Hux(1)','Hux(2)','Huu')
H