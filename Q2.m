clc; 
clear all;
close all;

N=1000;
m1=[1 1]';
m2=[14 7]';
m3=[16 1]';
s1=[5 3; 3 4];
s2=[5 3; 3 4];
s3=[5 3; 3 4];

x1=[5 2]';
x2=[17 5]';
x3=[9 2]';

pw=[1/3 1/3 1/3];




mu = m1;
Sigma =s1;
rng('default')  % For reproducibility
R1 = mvnrnd(mu,Sigma,333);

mu = m2;
Sigma =s2;
rng('default')  % For reproducibility
R2 = mvnrnd(mu,Sigma,333);


mu = m3;
Sigma =s3;
rng('default')  % For reproducibility
R3 = mvnrnd(mu,Sigma,334);

figure;
plot(R1(:,1),R1(:,2),'+')
hold on
plot(R2(:,1),R2(:,2),'o')
plot(R3(:,1),R3(:,2),'*')

title('equi propable classes')
legend('class 1','class 2','class 3')
hold off








% basian classification;

text="Basian classification";
disp(text)
disp(' ')


a1=Bayesian(x1);
text="given x1 belongs to class "+ num2str(a1);
disp(text)

a2=Bayesian(x2);
text2="given x2 belongs to class "+ num2str(a2);
disp(text2)

a3=Bayesian(x3);
text3="given x3 belongs to class "+ num2str(a3);
disp(text3)

%end

% Mahalanobis distance classification
disp(' ')
text="Mahalanobis distance classification";
disp(text)
disp(' ')

a1=md(x1);
text="given x1 belongs to class "+ num2str(a1);
disp(text)

a2=md(x2);
text2="given x2 belongs to class "+ num2str(a2);
disp(text2)

a3=md(x3);
text3="given x3 belongs to class "+ num2str(a3);
disp(text3)




% Euclidean distance classification
disp(' ')
text="Euclidean distance classification";
disp(text)
disp(' ')

a1=ed(x1);
text="given x1 belongs to class "+ num2str(a1);
disp(text)

a2=ed(x2);
text2="given x2 belongs to class "+ num2str(a2);
disp(text2)

a3=ed(x3);
text3="given x3 belongs to class "+ num2str(a3);
disp(text3)























function out=Bayesian(x)
m1=[1 1]';
m2=[14 7]';
m3=[16 1]';
s=[5 3; 3 4];

pw=[1/3 1/3 1/3];

px1_w1=1/(2*pi*det(s))*exp((-1/2)*(x-m1)'*inv(s)*(x-m1));
px1_w2=1/(2*pi*det(s))*exp((-1/2)*(x-m2)'*inv(s)*(x-m2));
px1_w3=1/(2*pi*det(s))*exp((-1/2)*(x-m3)'*inv(s)*(x-m3));

p=[px1_w1 px1_w2 px1_w3];

for i=1:3
    if p(i)==max(p)
        %text="given x belongs to class "+ num2str(i);
    %disp(text)
    out=i;
    end
end
end



function out=md(x)
m1=[1 1]';
m2=[14 7]';
m3=[16 1]';
s=[5 3; 3 4];


px1_w1=(x-m1)'*inv(s)*(x-m1);
px1_w2=(x-m2)'*inv(s)*(x-m2);
px1_w3=(x-m3)'*inv(s)*(x-m3);

p=[px1_w1 px1_w2 px1_w3];

for i=1:3
    if p(i)==min(p)
        %text="given x belongs to class "+ num2str(i);
    %disp(text)
    out=i;
    end
end
end


function out=ed(x)
m1=[1 1]';
m2=[14 7]';
m3=[16 1]';
s=[5 3; 3 4];


px1_w1=norm(x-m1);
px1_w2=norm(x-m2);
px1_w3=norm(x-m3);

p=[px1_w1 px1_w2 px1_w3]

for i=1:3
    if p(i)==min(p)
        %text="given x belongs to class "+ num2str(i);
    %disp(text)
    out=i;
    end
end
end




