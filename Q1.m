clc; 
clear all;
close all;

N=1200;
m1=[1 1]';
m2=[7 7]';
m3=[15 1]';
s1=[12 0; 0 1];
s2=[8 3;3 2];
s3=[2 0; 0 2];



mu = m1;
Sigma =s1;
rng('default')  % For reproducibility
R1 = mvnrnd(mu,Sigma,N/3);

mu = m2;
Sigma =s2;
rng('default')  % For reproducibility
R2 = mvnrnd(mu,Sigma,N/3);


mu = m3;
Sigma =s3;
rng('default')  % For reproducibility
R3 = mvnrnd(mu,Sigma,N/3);

figure;
plot(R1(:,1),R1(:,2),'+')
hold on
plot(R2(:,1),R2(:,2),'o')
plot(R3(:,1),R3(:,2),'*')

title('equi propable classes')
legend('class 1','class 2','class 3')
hold off


% % b part
p=[0.6 0.3 00.1'];


mu = m1;
Sigma =s1;
rng('default')  % For reproducibility
R1 = mvnrnd(mu,Sigma,N*p(1));

mu = m2;
Sigma =s2;
rng('default')  % For reproducibility
R2 = mvnrnd(mu,Sigma,N*p(2));


mu = m3;
Sigma =s3;
rng('default')  % For reproducibility
R3 = mvnrnd(mu,Sigma,N*p(3));

figure;
plot(R1(:,1),R1(:,2),'+')
hold on
plot(R2(:,1),R2(:,2),'o')
plot(R3(:,1),R3(:,2),'*')
title('for a given propability classes')
legend('class 1','class 2','class 3')
hold off

