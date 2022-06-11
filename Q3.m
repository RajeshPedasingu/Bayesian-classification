
clc; 
clear all;
close all;


%A = xlsread('iristrain.xlsx');
data = readtable('iristrain.xlsx');

p1=data.Var5;
p2=zeros(135,1);
t1='Iris-setosa';
t2='Iris-versicolor';
t3='Iris-virginica';

for i=1:135
    p3=cell2mat(p1(i));
    if strcmp(p3,t1)
        p2(i)=1;
    elseif strcmp(p3,t2)
        p2(i)=2;
    elseif strcmp(p3,t3)
        p2(i)=3;
    end
end

p2;


%for test data convertion 
data1 = readtable('iristest.xlsx');
pp1=data1.Var5;
pp2=zeros(15,1);
t1='Iris-setosa';
t2='Iris-versicolor';
t3='Iris-virginica';

for i=1:15
    pp3=cell2mat(pp1(i));
    if strcmp(pp3,t1)
        pp2(i)=1;
    elseif strcmp(pp3,t2)
        pp2(i)=2;
    elseif strcmp(pp3,t3)
        pp2(i)=3;
    end
end

pp2;






%prior probabilities 
pw1=sum(p2==1)/135;
pw2=sum(p2==2)/135;
pw3=sum(p2==3)/135;

%mean vector 

c1=(p2==1);
m1=sum(c1.*data.Var1)/sum(p2==1);
m2=sum(c1.*data.Var2)/sum(p2==1);
m3=sum(c1.*data.Var3)/sum(p2==1);
m4=sum(c1.*data.Var4)/sum(p2==1);

mx1=[m1 m2 m3 m4];

c1=(p2==2);
m1=sum(c1.*data.Var1)/sum(p2==2);
m2=sum(c1.*data.Var2)/sum(p2==2);
m3=sum(c1.*data.Var3)/sum(p2==2);
m4=sum(c1.*data.Var4)/sum(p2==2);

mx2=[m1 m2 m3 m4];

c1=(p2==3);
m1=sum(c1.*data.Var1)/sum(p2==3);
m2=sum(c1.*data.Var2)/sum(p2==3);
m3=sum(c1.*data.Var3)/sum(p2==3);
m4=sum(c1.*data.Var4)/sum(p2==3);

mx3=[m1 m2 m3 m4];


% co-variance by considering all 4 feture vectors at once.

cv1=[data.Var1(p2==1) data.Var2(p2==1) data.Var3(p2==1) data.Var4(p2==1)];
cv_c1=cov(cv1);

cv2=[data.Var1(p2==2) data.Var2(p2==2) data.Var3(p2==2) data.Var4(p2==2)];
cv_c2=cov(cv2);

cv3=[data.Var1(p2==3) data.Var2(p2==3) data.Var3(p2==3) data.Var4(p2==3)];
cv_c3=cov(cv3);


%%%%

d1=zeros(1,2,6);
kovar1=zeros(2,2,6);
kovar2=zeros(2,2,6);
kovar3=zeros(2,2,6);
l=1;
k1 = [data.Var1(p2==1) data.Var2(p2==1) data.Var3(p2==1) data.Var4(p2==1) ];
k2 = [data.Var1(p2==2) data.Var2(p2==2) data.Var3(p2==2) data.Var4(p2==2) ];
k3 = [data.Var1(p2==3) data.Var2(p2==3) data.Var3(p2==3) data.Var4(p2==3) ];

for i=1:4
    for j=1:4
        if i~=j
            if i<j
             mean1=[mx1(i) mx1(j)];
             mean2=[mx2(i) mx2(j)];
             mean3=[mx3(i) mx3(j)];
             covar1=cov(k1(:,i),k1(:,j));
             kovar1(:,:,l)=covar1;
             covar2=cov(k2(:,i),k2(:,j));
             kovar2(:,:,l)=covar2;
             covar3=cov(k3(:,i),k3(:,j));
             kovar3(:,:,l)=covar3;
             
             d1(1,:,l)=mean1;
             d1(2,:,l)=mean2;
             d1(3,:,l)=mean3;
             l=l+1;
            end
        end
    end
end



% all parammeters;

% mean 
% d1
% co variance
% kovar1,2,3




% basian classificatiion 
x1=[5.1 0.2]';
n1=3;

a1=Bayesian(x1,n1,d1,kovar1,kovar2,kovar3);
if a1==1
text="given x1 belongs to class "+ t1 ;
elseif a1==2
    text="given x1 belongs to class "+ t2 ;
elseif a1==3
    text="given x1 belongs to class "+ t3 ;
end 
disp(text)




x1=[5.7 4.1]';
n1=2;

a1=Bayesian(x1,n1,d1,kovar1,kovar2,kovar3);
if a1==1
text="given x2 belongs to class "+ t1 ;
elseif a1==2
    text="given x2 belongs to class "+ t2 ;
elseif a1==3
    text="given x2 belongs to class "+ t3 ;
end 
disp(text)



x1=[6.5 2]';
n1=3;

a1=Bayesian(x1,n1,d1,kovar1,kovar2,kovar3);
if a1==1
text="given x3 belongs to class "+ t1 ;
elseif a1==2
    text="given x3 belongs to class "+ t2 ;
elseif a1==3
    text="given x3 belongs to class "+ t3 ;
end 
disp(text)




% testing data

mat1=[data1.Var1 data1.Var2 data1.Var3 data1.Var4];

bb1=zeros(15,1);

n1=0;

nn1=zeros(6,1);
for i=1:4
    for j=1:4
        if i~=j 
                if i<j
                    n1=n1+1;
                    count=0;
                    count1=0;
                    count2=0;
                        for k=1:15
                            x1=[mat1(k,i) mat1(k,j)]';
                            a1=Bayesian(x1,n1,d1,kovar1,kovar2,kovar3);
                            

                             if a1==pp2(k)
                                count=count+1;
                             else
                                % [n1 x1']
                             end
                             a2=mohalanobis(x1,n1,d1,kovar1,kovar2,kovar3);
                             if a2==pp2(k)
                                count1=count1+1;
                             else
                                % [n1 x1']
                             end
                             a3=Euclidian(x1,n1,d1,kovar1,kovar2,kovar3);
                             if a3==pp2(k)
                                count2=count2+1;
                             else
                                % [n1 x1']
                             end



                        end
                       nn1(n1)=count;
                       nn2(n1)=count1;
                       nn3(n1)=count2;
                end
        end
    end

   end

nn1;
nn2;
nn3;
error_bayesian = (15-nn1)*100/15
error_mohalanobis = (15-nn2')*100/15
error_Euclidian = (15-nn3')*100/15
% 
% 
% 


% Q3_c part ----------------------------


X=[data.Var1 data.Var2 data.Var3 data.Var4];

X=[data.Var1 data.Var2];
y=categorical(data.Var5);

classifier_name = {'Bayes classification by taking 2 fetures [A and B] '};
classifier{1} = fitcnb(X,y);


x1range = min(X(:,1)):.01:max(X(:,1));
x2range = min(X(:,2)):.01:max(X(:,2));
[xx1, xx2] = meshgrid(x1range,x2range);
XGrid = [xx1(:) xx2(:)];


predictedspecies = predict(classifier{1},XGrid);

figure;
gscatter(xx1(:), xx2(:), predictedspecies,'rgb');

title(classifier_name{1})
legend off, axis tight

legend('Iris-setosa','Iris-versicolor','Iris-virginica')

 % %

X=[data.Var3 data.Var4];
y=categorical(data.Var5);

classifier_name = {'Bayes classification by taking 2 fetures [C and D] '};
classifier{1} = fitcnb(X,y);
                

x1range = min(X(:,1)):.01:max(X(:,1));
x2range = min(X(:,2)):.01:max(X(:,2));
[xx1, xx2] = meshgrid(x1range,x2range);
XGrid = [xx1(:) xx2(:)];


predictedspecies = predict(classifier{1},XGrid);

figure;
gscatter(xx1(:), xx2(:), predictedspecies,'rgb');

title(classifier_name{1})
legend off, axis tight

legend('Iris-setosa','Iris-versicolor','Iris-virginica')

% %

X=[data.Var1 data.Var4];
y=categorical(data.Var5);

classifier_name = {'Bayes classification by taking 2 fetures [A and D] '};
classifier{1} = fitcnb(X,y);
                

x1range = min(X(:,1)):.01:max(X(:,1));
x2range = min(X(:,2)):.01:max(X(:,2));
[xx1, xx2] = meshgrid(x1range,x2range);
XGrid = [xx1(:) xx2(:)];


predictedspecies = predict(classifier{1},XGrid);

figure;
gscatter(xx1(:), xx2(:), predictedspecies,'rgb');

title(classifier_name{1})
legend off, axis tight

legend('Iris-setosa','Iris-versicolor','Iris-virginica')



 



function out1 = Bayesian(x,n1,d1,kovar1,kovar2,kovar3)
mad1=d1(1,:,:);
mad2=d1(2,:,:);
mad3=d1(3,:,:);
s1=kovar1;
s2=kovar2;
s3=kovar3;

pw=[1/3 1/3 1/3];

if n1==1
    m1=mad1(:,:,1)';
    m2=mad2(:,:,1)';
    m3=mad3(:,:,1)';
    r1=s1(:,:,1);
    r2=s2(:,:,1);
    r3=s3(:,:,1);
px1_w1=1/(2*pi*det(r1))*exp((-1/2)*(x-m1)'*inv(r1)*(x-m1));
px1_w2=1/(2*pi*det(r2))*exp((-1/2)*(x-m2)'*inv(r2)*(x-m2));
px1_w3=1/(2*pi*det(r3))*exp((-1/2)*(x-m3)'*inv(r3)*(x-m3));


p=[px1_w1 px1_w2 px1_w3];


elseif n1==2
    m1=mad1(:,:,2)';
    m2=mad2(:,:,2)';
    m3=mad3(:,:,2)';
    r1=s1(:,:,2);
    r2=s2(:,:,2);
    r3=s3(:,:,2);
px1_w1=1/(2*pi*det(r1))*exp((-1/2)*(x-m1)'*inv(r1)*(x-m1));
px1_w2=1/(2*pi*det(r2))*exp((-1/2)*(x-m2)'*inv(r2)*(x-m2));
px1_w3=1/(2*pi*det(r3))*exp((-1/2)*(x-m3)'*inv(r3)*(x-m3));
p=[px1_w1 px1_w2 px1_w3];

elseif n1==3
    m1=mad1(:,:,3)';
    m2=mad2(:,:,3)';
    m3=mad3(:,:,3)';
    r1=s1(:,:,3);
    r2=s2(:,:,3);
    r3=s3(:,:,3);
px1_w1=1/(2*pi*det(r1))*exp((-1/2)*(x-m1)'*inv(r1)*(x-m1));
px1_w2=1/(2*pi*det(r2))*exp((-1/2)*(x-m2)'*inv(r2)*(x-m2));
px1_w3=1/(2*pi*det(r3))*exp((-1/2)*(x-m3)'*inv(r3)*(x-m3));
p=[px1_w1 px1_w2 px1_w3];


elseif n1==4
    m1=mad1(:,:,4)';
    m2=mad2(:,:,4)';
    m3=mad3(:,:,4)';
    r1=s1(:,:,4);
    r2=s2(:,:,4);
    r3=s3(:,:,4);
px1_w1=1/(2*pi*det(r1))*exp((-1/2)*(x-m1)'*inv(r1)*(x-m1));
px1_w2=1/(2*pi*det(r2))*exp((-1/2)*(x-m2)'*inv(r2)*(x-m2));
px1_w3=1/(2*pi*det(r3))*exp((-1/2)*(x-m3)'*inv(r3)*(x-m3));
p=[px1_w1 px1_w2 px1_w3];


elseif n1==5
    m1=mad1(:,:,5)';
    m2=mad2(:,:,5)';
    m3=mad3(:,:,5)';
    r1=s1(:,:,5);
    r2=s2(:,:,5);
    r3=s3(:,:,5);
px1_w1=1/(2*pi*det(r1))*exp((-1/2)*(x-m1)'*inv(r1)*(x-m1));
px1_w2=1/(2*pi*det(r2))*exp((-1/2)*(x-m2)'*inv(r2)*(x-m2));
px1_w3=1/(2*pi*det(r3))*exp((-1/2)*(x-m3)'*inv(r3)*(x-m3));
p=[px1_w1 px1_w2 px1_w3];


elseif n1==6
    m1=mad1(:,:,6)';
    m2=mad2(:,:,6)';
    m3=mad3(:,:,6)';
    r1=s1(:,:,6);
    r2=s2(:,:,6);
    r3=s3(:,:,6);
px1_w1=1/(2*pi*det(r1))*exp((-1/2)*(x-m1)'*inv(r1)*(x-m1));
px1_w2=1/(2*pi*det(r2))*exp((-1/2)*(x-m2)'*inv(r2)*(x-m2));
px1_w3=1/(2*pi*det(r3))*exp((-1/2)*(x-m3)'*inv(r3)*(x-m3));
p=[px1_w1 px1_w2 px1_w3];

end



for i=1:3
    %display(p)
if p(i)==max(p)
   
    %text="given x belongs to class "+ num2str(i);
    %disp(text)
out1=i;
end
end

end










% Mahalanobis distance classification
function out1 = mohalanobis(x,n1,d1,kovar1,kovar2,kovar3)
mad1=d1(1,:,:);
mad2=d1(2,:,:);
mad3=d1(3,:,:);
s1=kovar1;
s2=kovar2;
s3=kovar3;

pw=[1/3 1/3 1/3];

if n1==1
    m1=mad1(:,:,1)';
    m2=mad2(:,:,1)';
    m3=mad3(:,:,1)';
    r1=s1(:,:,1);
    r2=s2(:,:,1);
    r3=s3(:,:,1);

px1_w1=(x-m1)'*inv(r1)*(x-m1);
px1_w2=(x-m2)'*inv(r2)*(x-m2);
px1_w3=(x-m3)'*inv(r3)*(x-m3);



p=[px1_w1 px1_w2 px1_w3];


elseif n1==2
    m1=mad1(:,:,2)';
    m2=mad2(:,:,2)';
    m3=mad3(:,:,2)';
    r1=s1(:,:,2);
    r2=s2(:,:,2);
    r3=s3(:,:,2);
px1_w1=(x-m1)'*inv(r1)*(x-m1);
px1_w2=(x-m2)'*inv(r2)*(x-m2);
px1_w3=(x-m3)'*inv(r3)*(x-m3);


p=[px1_w1 px1_w2 px1_w3];

elseif n1==3
    m1=mad1(:,:,3)';
    m2=mad2(:,:,3)';
    m3=mad3(:,:,3)';
    r1=s1(:,:,3);
    r2=s2(:,:,3);
    r3=s3(:,:,3);

px1_w1=(x-m1)'*inv(r1)*(x-m1);
px1_w2=(x-m2)'*inv(r2)*(x-m2);
px1_w3=(x-m3)'*inv(r3)*(x-m3);
    
    p=[px1_w1 px1_w2 px1_w3];


elseif n1==4
    m1=mad1(:,:,4)';
    m2=mad2(:,:,4)';
    m3=mad3(:,:,4)';
    r1=s1(:,:,4);
    r2=s2(:,:,4);
    r3=s3(:,:,4);
px1_w1=(x-m1)'*inv(r1)*(x-m1);
px1_w2=(x-m2)'*inv(r2)*(x-m2);
px1_w3=(x-m3)'*inv(r3)*(x-m3);

p=[px1_w1 px1_w2 px1_w3];


elseif n1==5
    m1=mad1(:,:,5)';
    m2=mad2(:,:,5)';
    m3=mad3(:,:,5)';
    r1=s1(:,:,5);
    r2=s2(:,:,5);
    r3=s3(:,:,5);
px1_w1=(x-m1)'*inv(r1)*(x-m1);
px1_w2=(x-m2)'*inv(r2)*(x-m2);
px1_w3=(x-m3)'*inv(r3)*(x-m3);


p=[px1_w1 px1_w2 px1_w3];


elseif n1==6
    m1=mad1(:,:,6)';
    m2=mad2(:,:,6)';
    m3=mad3(:,:,6)';
    r1=s1(:,:,6);
    r2=s2(:,:,6);
    r3=s3(:,:,6);
px1_w1=(x-m1)'*inv(r1)*(x-m1);
px1_w2=(x-m2)'*inv(r2)*(x-m2);
px1_w3=(x-m3)'*inv(r3)*(x-m3);


p=[px1_w1 px1_w2 px1_w3];

end



for i=1:3
    %display(p)
if p(i)==min(p)
   
    %text="given x belongs to class "+ num2str(i);
    %disp(text)
out1=i;
end
end

end






% Euclidian distance classification
function out1 = Euclidian(x,n1,d1,kovar1,kovar2,kovar3)
mad1=d1(1,:,:);
mad2=d1(2,:,:);
mad3=d1(3,:,:);
s1=kovar1;
s2=kovar2;
s3=kovar3;

pw=[1/3 1/3 1/3];

if n1==1
    m1=mad1(:,:,1)';
    m2=mad2(:,:,1)';
    m3=mad3(:,:,1)';
    r1=s1(:,:,1);
    r2=s2(:,:,1);
    r3=s3(:,:,1);

px1_w1=norm(x-m1);
px1_w2=norm(x-m2);
px1_w3=norm(x-m3);


p=[px1_w1 px1_w2 px1_w3];


elseif n1==2
    m1=mad1(:,:,2)';
    m2=mad2(:,:,2)';
    m3=mad3(:,:,2)';
    r1=s1(:,:,2);
    r2=s2(:,:,2);
    r3=s3(:,:,2);
px1_w1=norm(x-m1);
px1_w2=norm(x-m2);
px1_w3=norm(x-m3);


p=[px1_w1 px1_w2 px1_w3];

elseif n1==3
    m1=mad1(:,:,3)';
    m2=mad2(:,:,3)';
    m3=mad3(:,:,3)';
    r1=s1(:,:,3);
    r2=s2(:,:,3);
    r3=s3(:,:,3);

px1_w1=norm(x-m1);
px1_w2=norm(x-m2);
px1_w3=norm(x-m3);
    
    p=[px1_w1 px1_w2 px1_w3];


elseif n1==4
    m1=mad1(:,:,4)';
    m2=mad2(:,:,4)';
    m3=mad3(:,:,4)';
    r1=s1(:,:,4);
    r2=s2(:,:,4);
    r3=s3(:,:,4);
px1_w1=norm(x-m1);
px1_w2=norm(x-m2);
px1_w3=norm(x-m3);

p=[px1_w1 px1_w2 px1_w3];


elseif n1==5
    m1=mad1(:,:,5)';
    m2=mad2(:,:,5)';
    m3=mad3(:,:,5)';
    r1=s1(:,:,5);
    r2=s2(:,:,5);
    r3=s3(:,:,5);
px1_w1=norm(x-m1);
px1_w2=norm(x-m2);
px1_w3=norm(x-m3);

p=[px1_w1 px1_w2 px1_w3];


elseif n1==6
    m1=mad1(:,:,6)';
    m2=mad2(:,:,6)';
    m3=mad3(:,:,6)';
    r1=s1(:,:,6);
    r2=s2(:,:,6);
    r3=s3(:,:,6);
px1_w1=norm(x-m1);
px1_w2=norm(x-m2);
px1_w3=norm(x-m3);


p=[px1_w1 px1_w2 px1_w3];

end



for i=1:3
    %display(p)
if p(i)==min(p)
   
    %text="given x belongs to class "+ num2str(i);
    %disp(text)
out1=i;
end
end

end















