function D = decisionRegion(features, targets, region, display)

% Classification using the Adaptive Kernel Density estimation method
% with automatic estimation of the smooth parameters based on the Least-Squares Cross-Validation method.
% Inputs:
% 	features: - computed ITD/IID
%	targets	- ideal binary classification
%	region - domain of interest  [min_itd max_itd min_iid max_iid]

% Outputs 
%    D - decision rule

N		= region(5); %Number of points on the grid
x		= ones(N,1) * linspace (region(1),region(2),N);
y		= linspace (region(3),region(4),N)' * ones(1,N);
V0		= zeros(N);
V1		= zeros(N);

FAST=1;

class1  = find(targets == 1);
class0  = find(targets == 0);


%class0

P0 = length(class0)/length(features);
n= length(class0);

%automatic estimation of h
if(FAST==0)

DIFX1=zeros(1,n*(n-1)/2);
DIFX2=zeros(1,n*(n-1)/2);

ind=0;

for i=1:n-1
    for j=i+1:n
        ind=ind+1;
        DIFX1(ind)= (features(1,class0(i))-features(1,class0(j)))^2;
        DIFX2(ind)= (features(2,class0(i))-features(2,class0(j)))^2;
    end
end
end

sigma=cov(features(1,class0),features(2,class0));

hstart1=n^(-1/6)*sqrt(sigma(1,1));
h_min1=1/4*hstart1;
h_max1=3/2*hstart1;

hstart2=n^(-1/6)*sqrt(sigma(2,2));
h_min2=1/4*hstart2;
h_max2=3/2*hstart2;

if (FAST==0)
options=optimset('LargeScale','off','Display','iter','TolFun',0.001);
h0= fmincon(@myfun,[hstart1 hstart2],[],[],[],[],[h_min1 ,h_min2],[h_max1,h_max2],[],options,DIFX1,DIFX2,n)
h01=max(h0(1),0.1);h02=max(h0(2),0.1);

clear DIFX1 DIFX2
else
h01=hstart1;
h02=hstart2;
end
%landa estimation

f_X=zeros(1,n);

for i = 1:n
    
   temp = (features(1,class0)-features(1,class0(i))).^2/(2*h01^2) + (features(2,class0)-features(2,class0(i))).^2/(2*h02^2);
   f_X   = f_X + exp(-temp);
   
end

f_X = 1/(2*pi*h01*h02)*f_X/n;


g=exp(sum(log(f_X))/n);

alpha=-0.5;

landa=(f_X/g).^alpha;

% density estimation

for i = 1:n,
   temp = ((x - features(1,class0(i))).^2/(2*h01^2) + (y - features(2,class0(i))).^2/(2*h02^2))/(landa(i)^2);   
   V0   = V0 + 1/(h01*h02*landa(i)^2)*exp(-temp);
end

V0 = 1/(2*pi)*V0/n;

clear landa f_X g

if display==1
    
figure
mesh(V0)

end  %display

% class1

P1 = length(class1)/length(features);
n  = length(class1);

%automatic estimation of h
if FAST==0
DIFX1=zeros(1,n*(n-1)/2);
DIFX2=zeros(1,n*(n-1)/2);

ind=0;

for i=1:n-1
    for j=i+1:n
        ind=ind+1;
        DIFX1(ind)= (features(1,class1(i))-features(1,class1(j)))^2;
        DIFX2(ind)= (features(2,class1(i))-features(2,class1(j)))^2;
    end
end
end

sigma=cov(features(1,class1),features(2,class1));

hstart1=n^(-1/6)*sqrt(sigma(1,1));
h_min1=1/4*hstart1;
h_max1=3/2*hstart1;

hstart2=n^(-1/6)*sqrt(sigma(2,2));
h_min2=1/4*hstart2;
h_max2=3/2*hstart2;

if FAST==0
options=optimset('LargeScale','off','Display','iter','TolFun',0.001);
h1= fmincon(@myfun,[hstart1 hstart2],[],[],[],[],[h_min1 ,h_min2],[h_max1,h_max2],[],options,DIFX1,DIFX2,n)
h11=max(h1(1),0.1);h12=max(h1(2),0.1);
clear DIFX1 DIFX2
else
    h11=hstart1;
    h12=hstart2;
end

%landa estimation

f_X=zeros(1,n);

for i = 1:n
    
   temp = (features(1,class1)-features(1,class1(i))).^2/(2*h11^2) + (features(2,class1)-features(2,class1(i))).^2/(2*h12^2);
   f_X   = f_X + exp(-temp);
   
end

f_X = 1/(2*pi*h11*h12)*f_X/n;

g=exp(sum(log(f_X))/n);

alpha=-0.5;

landa=(f_X/g).^alpha;

% density estimation

for i = 1:n,
   temp = ((x - features(1,class1(i))).^2/(2*h11^2) + (y - features(2,class1(i))).^2/(2*h12^2))/(landa(i)^2);   
   V1   = V1 + 1/(h11*h12*landa(i)^2)*exp(-temp);
end

V1 = 1/(2*pi)*V1/n;

if display==1
figure
mesh(V1)
end %display

clear landa f_X g

D = (V0*P0 < V1*(1-P0));

if display==1
figure
imshow(D,[])
end %display

function M = myfun(h,DIFX1,DIFX2,n)

h1=h(1);
h2=h(2);

D1=DIFX1./(h1^2);D2=DIFX2./(h2^2);D=D1+D2;D1=[];D2=[];
D=exp(-D/4)-4*exp(-D/2);
M=2*sum(D);
M=M/(4*pi*n^2*h1*h2)+1/(4*pi*n*h1*h2);

