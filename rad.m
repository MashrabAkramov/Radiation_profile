clear all; close all; clc; clearvars;

N=400; L=4;
x=linspace(-L,L,2*N+2); dx=x(2)-x(1);
dt=0.2*dx.^2;

J=20000;
J=1;
t0=0; for j=1:J+1 t(j)=(j-1)*dt; end
aa=-1i*dt/(2*dx^2);

%% ABC
alpha = 1; g0 = 20;
sink = exp(-g0*dt*(1./cosh(alpha*(x-L)).^2 + 1./cosh(alpha*(x+L)).^2));
q=zeros(2*N+2,J+1);
% a=0.5; b=3; 
% %% The initial condition
% for n=1:2*N+2
%     if x(n)>=-a && x(n)<=a
%         q(n,1)=b;
%     end
% end

z0=0;
k0=5;
sigma=.3; A=2;
q(:,1)=A*exp(1i*k0*x-(x-z0).^2/(2*sigma)^2)/sqrt(sqrt(2*pi)*sigma);%+A*exp(-1i*k0*x-(x-z0).^2/(2*sigma)^2)/sqrt(sqrt(2*pi)*sigma);

%%
center=ones(1,2*N+2); left=ones(1,2*N+1);
AA=aa*diag(left,-1)+aa*diag(left,1)-(2*aa+1)*diag(center,0);
BB=-aa*diag(left,-1)-aa*diag(left,1)+(2*aa-1)*diag(center,0);

AA(N+1,N)=0; AA(N+1,N+1)=1; AA(N+1,N+2)=-1; 
BB(N+1,N)=0; BB(N+1,N+1)=0; BB(N+1,N+2)=0; 

AA(N+2,N+1)=1; AA(N+2,N+2)=1; AA(N+2,N+3)=0; 
BB(N+2,N+1)=1; BB(N+2,N+1)=0; BB(N+2,N+2)=0; BB(N+2,N+3)=1; 

AA=sparse(AA); BB=sparse(BB);

% for j=1:J
%     for n=1:2*N+2
%         q(n,j)=q(n,j)*exp(-2*1i*dt*q(n,j)*conj(q(n,j))/2);
%     end
%     
%     q(:,j+1)=AA\BB*q(:,j);
%     
%     for n=1:2*N+2
%         q(n,j+1)=q(n,j+1)*exp(-2*1i*dt*q(n,j+1)*conj(q(n,j+1))/2);
%     end
% %     q(:,j+1)=sink'.*q(:,j+1); 
% %     plot(x,abs(q(:,j+1)).^2); %ylim([0 12])
% %     drawnow
%     fprintf('%s%i\n','j=',j)
% end

% plot(x,abs(q(:,2000)).^2); %ylim([0 2])


time=1;
for kk=1:50:j+1
    Q(:,time)=q(:,kk);
    time=time+1;
    plot(x,abs(q(:,kk)).^2); %ylim([0 160])
    drawnow
end


figure(4)
plot(x,abs(q(:,1)).^2,'r','Linewidth',2); %ylim([0 160]); 
ylabel('|q|^2','Fontsize',16); xlabel('x','Fontsize',16)
hold on
plot(x,abs(Q(:,40)).^2,'b','Linewidth',2); %ylim([0 160])
hold on
plot(x,abs(Q(:,80)).^2,'k','Linewidth',2); %ylim([0 160])
legend('t=0','t=0.2','t=0.4')




%% Star graph
beta_m1=46/70; beta_p1=46/58; 
beta_p2=46/22; beta_m2=46/17;
beta_p3=46/16; beta_m3=46/15;

alpha_m1=sqrt(beta_m1);
alpha_p1=sqrt(beta_p1);
alpha_p2=sqrt(beta_p2);

alpha_m2=sqrt(beta_m2);
alpha_p3=sqrt(beta_p3);
alpha_m3=sqrt(beta_m3);

q_m1=sqrt(2/beta_m1)*Q(1:N+1,:)/2; 
q_m2=sqrt(2/beta_m2)*Q(1:N+1,:)/2; 
q_m3=sqrt(2/beta_m3)*Q(1:N+1,:)/2;
q_p1=sqrt(2/beta_p1)*Q(N+2:2*N+2,:)/2; 
q_p2=sqrt(2/beta_p2)*Q(N+2:2*N+2,:)/2; 
q_p3=sqrt(2/beta_p3)*Q(N+2:2*N+2,:)/2; 

q_m1(:,1)=sqrt(2/beta_m1)*q(1:N+1,1)/2; 
q_m2(:,1)=sqrt(2/beta_m2)*q(1:N+1,1)/2; 
q_m3(:,1)=sqrt(2/beta_m3)*q(1:N+1,1)/2;
q_p1(:,1)=sqrt(2/beta_p1)*q(N+2:2*N+2,1)/2; 
q_p2(:,1)=sqrt(2/beta_p2)*q(N+2:2*N+2,1)/2; 
q_p3(:,1)=sqrt(2/beta_p3)*q(N+2:2*N+2,1)/2;

x1=x(1:N+1); x2=x(N+2:2*N+2);
%% Plotting
middle=60; final=120;
font=16;
figure('DefaultAxesFontSize',14); set(gcf,'Position',[200 0 600 1600]);
subplot(3,2,1)
plot(x1,abs(q_m1(:,1)).^2,'r','Linewidth',2); title('b_{-1}','Fontsize',font); xlabel('x','Fontsize',font);
ylabel('|q_{-1}|^2','Fontsize',font);
ylim([0 4])
hold on
plot(x1,abs(q_m1(:,middle)).^2,'g','Linewidth',2);
hold on
plot(x1,abs(q_m1(:,final)).^2,'b','Linewidth',2);

subplot(3,2,2)
plot(x2,abs(q_p1(:,1)).^2,'r','Linewidth',2); title('b_{1}','Fontsize',font); xlabel('x','Fontsize',font);
ylabel('|q_{1}|^2','Fontsize',font);
ylim([0 4])
hold on
plot(x2,abs(q_p1(:,middle)).^2,'g','Linewidth',2);
hold on
plot(x2,abs(q_p1(:,final)).^2,'b','Linewidth',2);

subplot(3,2,3)
plot(x1,abs(q_m2(:,1)).^2,'r','Linewidth',2); title('b_{-2}','Fontsize',font); xlabel('x','Fontsize',font);
ylabel('|q_{-2}|^2','Fontsize',font);
ylim([0 1.5])
hold on
plot(x1,abs(q_m2(:,middle)).^2,'g','Linewidth',2);
hold on
plot(x1,abs(q_m2(:,final)).^2,'b','Linewidth',2);

subplot(3,2,4)
plot(x2,abs(q_p2(:,1)).^2,'r','Linewidth',2);  title('b_{2}','Fontsize',font); xlabel('x','Fontsize',font);
ylabel('|q_{2}|^2','Fontsize',font);
ylim([0 1.5])
hold on
plot(x2,abs(q_p2(:,middle)).^2,'g','Linewidth',2);
hold on
plot(x2,abs(q_p2(:,final)).^2,'b','Linewidth',2);

subplot(3,2,5)
plot(x1,abs(q_m3(:,1)).^2,'r','Linewidth',2);  title('b_{-3}','Fontsize',font); xlabel('x','Fontsize',font);
ylabel('|q_{-3}|^2','Fontsize',font);
ylim([0 1])
hold on
plot(x1,abs(q_m3(:,middle)).^2,'g','Linewidth',2);
hold on
plot(x1,abs(q_m3(:,final)).^2,'b','Linewidth',2); 

subplot(3,2,6)
plot(x2,abs(q_p3(:,1)).^2,'r','Linewidth',2);  title('b_{3}','Fontsize',font); xlabel('x','Fontsize',font); 
ylabel('|q_{3}|^2','Fontsize',font);
ylim([0 1])
hold on
plot(x2,abs(q_p3(:,middle)).^2,'g','Linewidth',2); 
hold on
plot(x2,abs(q_p3(:,final)).^2,'b','Linewidth',2);
legend('t=0','t=0.06','t=0.12')





