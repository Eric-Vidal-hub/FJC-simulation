clearvars; close all;
tic
T = 1000;
N = 100;
b = 3.0;
L = T*(N+2);

filename = 'D:\Documents\QuanTEEM\study notes\Sem 1\soft matter\practicals\simulation_FJC_b=3.0_N=100_T=1000.xyz';
fid = fopen(filename);

xdata = zeros([T , N]);
ydata = zeros([T , N]);
zdata = zeros([T , N]);


for i = 1:T
    fgetl(fid); 
    fgetl(fid);
    for j=1:N+1
        line=fgetl(fid);
        data=strsplit(line);
        xdata(i,j)=str2double(data(2));
        ydata(i,j)=str2double(data(3));
        zdata(i,j)=str2double(data(4));
    end
end
fclose(fid);
%% Calculating Q and Rg

Q = zeros(T,1);
Rg = zeros(T,1);

for t=1:T
  Q(t)=sqrt((xdata(t,end)-xdata(t,1))^2+(ydata(t,end)-ydata(t,1))^2+(zdata(t,end)-zdata(t,1))^2);
  Rcm=[mean(xdata(t,:)) mean(ydata(t,:)) mean(zdata(t,:))];
  d2=(xdata(t,:)-Rcm(1)).^2+(ydata(t,:)-Rcm(2)).^2+(zdata(t,:)-Rcm(3)).^2;
  Rg(t)=sqrt(mean(d2));
end
toc
plots = figure;
figure(1);
plot(Q,'b.');
figure(2);
plot(Rg,'r.');

%% verification of theoretical results
Qtemp = zeros(T,1);
Qnum = zeros(N-10,1);
Qthr = zeros(N-10,1);

Rgtemp = zeros(T,1);
Rgnum = zeros(N-10,1);
Rgthr = zeros(N-10,1);

tic
for n = 10:N
    for t = 1:T
          Qtemp(t) = sqrt((xdata(t,n)-xdata(t,1))^2+(ydata(t,n)-ydata(t,1))^2+(zdata(t,n)-zdata(t,1))^2);
          Rcm=[mean(xdata(t,1:n)) mean(ydata(t,1:n)) mean(zdata(t,1:n))];
          d2=(xdata(t,1:n)-Rcm(1)).^2+(ydata(t,1:n)-Rcm(2)).^2+(zdata(t,1:n)-Rcm(3)).^2;
          Rgtemp(t)=sqrt(mean(d2));
    end
    Qnum(n-9) = mean(Qtemp.^2);
    Qthr(n-9) = sqrt(n*(b^2));

    Rgnum(n-9) = mean(Rgtemp.^2);
    Rgthr(n-9) = sqrt(n*(b^2)/6);
end
clear Qtemp;

titlestr = sprintf('Q model vs theory, N=%d T=%d',N,T);
figure;
subplot(2,1,1)
plot(10:N,Qnum(:,1),'o');
hold on;grid on;
plot(10:N,Qthr(:,1).^2);
legend('simulation','theory','Location','northwest')
title(titlestr,'FontSize',15)
xlabel('n','FontSize',15)
ylabel('Q^2','FontSize',15)


titlestr = sprintf('R_g model vs theory, N=%d T=%d',N,T);
subplot(2,1,2)
plot(10:N,Rgnum(:,1),'o');
hold on;grid on;
plot(10:N,Rgthr(:,1).^2);
legend('simulation','theory','Location','northwest')
title(titlestr,'FontSize',15)
xlabel('n','FontSize',15)
ylabel('R_g','FontSize',15)
toc

%% PDF of Q
n = 100;
PDF = zeros(T,1);
for t=1:T
      Q(t)=sqrt((xdata(t,n)-xdata(t,1)).^2+(ydata(t,n)-ydata(t,1)).^2+(zdata(t,n)-zdata(t,1)).^2);
end

Q = sort(Q);
for i = 1:T
    PDF(i) = 4*pi*Q(i)^2 * (3/(2*pi*n*b^2))^(3/2) * exp(-3*Q(i)^2/(2*n*b^2));
end

figure
Qhist = histogram(Q,'Normalization','pdf');
w = Qhist.BinWidth;
Qbins = Qhist.BinEdges;
Qhistvals = Qhist.Values;
bincenters = Qbins(1:end-1) + w/2;

plot(Q,PDF,'LineWidth',1.5,'Color','black')
hold on;grid on;
plot(bincenters,Qhistvals,'o','MarkerSize',5)
legend('Simulation','Model')

%% structure factor
n=100;
k = 0:0.00001:0.01;
I = zeros(length(k),1);
t = 1;
Guinier = zeros(length(k),1);

for l=1:length(k)
    Rcm=[mean(xdata(1,1:n)) mean(ydata(1,1:n)) mean(zdata(1,1:n))];
    d2=(xdata(1,1:n)-Rcm(1)).^2+(ydata(1,1:n)-Rcm(2)).^2+(zdata(1,1:n)-Rcm(3)).^2;
    Rg=mean(d2);
    Guinier(l) = ((n+1)^2) .* (1-(((k(l)*Rg).^2)/3));
    for i =1:n
        for j=1:n
            if j~=i                         % sin(k*dist)/k*dist is singular at i=j
                dist = abs((xdata(1:T,i).^2 + ydata(1:T,i).^2 + zdata(1:T,i).^2) - (xdata(1:T,j).^2 + ydata(1:T,j).^2 + zdata(1:T,j).^2));
                I(l) = I(l) + mean(sin(k(l).*dist)./(k(l).*dist));
            else
                I(l) = I(l) + 1;            % sinc(0) = 1        
            end
        end
    end
end

scatfig = figure;
plot(k,Guinier);
hold on; grid on;
plot(k,I,'x');
legend('Guinier Approximation','Numerical Results')

%% yeah idk
t=1;            % whatevrer random time frame
diff = zeros(1,N-1);
for i = 1:N-1
    xdiff = xdata(t,i) - xdata(t,i+1);
    ydiff = ydata(t,i) - ydata(t,i+1);
    zdiff = zdata(t,i) - zdata(t,i+1);
    diff(i) = (xdiff^2 + ydiff^2 + zdiff^2)^0.5;
end
mean(diff)
%% mont - e - carlo
% set up molecule
t=  1;            % pick some time step for a starting config
n= 100;            % number of bonds used in simulation
F = [0.01 0 0];   % force in Newtons [Fx Fy Fz]
polymer = [xdata(t,1:n); ydata(t,1:n); zdata(t,1:n)];
naccept = 100;
T=250;
kb=1.3806503e-23;
Q = zeros(3,naccept);
backup = 0;

% monte the carlo
i=1;
while i<=naccept
    polymer2 = polymer;
    bm = floor(rand()*(n-1)) + 1;              % random number to select bond to modify
    newx = b*rand();                           % first one is random in [0,b]
    newy = sqrt((b^2 - newx^2))*rand();        % second and third bonds limited to ensure random values don't exceed b^2
    newz = sqrt(b^2 - (newx^2 + newy^2));
    incrs = zeros(3,n-bm);                     % bond vectors relative to previous monomer instead of origin
    % fprintf('\nnewx =%d   newy =%d   newz =%d   ',newx,newy,newz)
    for j = bm:n-1      
        for e=1:3
            incrs(e,j+1) = polymer(e,j+1) - polymer(e,j);
        end
    end

    for j=bm:n
        polymer2(1,j+1) = polymer(1,j) + newx;
        polymer2(2,j+1) = polymer(2,j) + newy;
        polymer2(3,j+1) = polymer(3,j) + newz;
    end
    Q1 = ([(polymer(1,1)-polymer(1,end)); (polymer(2,1)-polymer(2,end)); (polymer(3,1)-polymer(3,end))]);
    Q2 = ([(polymer2(1,1)-polymer2(1,end)); (polymer2(2,1)-polymer2(2,end)); (polymer2(3,1)-polymer2(3,end))]);
    V1 = dot(-F,Q1);
    V2 = dot(-F,Q2);
    
    if((V1>=V2))
        fprintf('\n  new bond x=%d   old bond x =%d',polymer2(1,bm+1),polymer(1,bm+1))
        Q(:,i) = (Q2);
        i = i+1;
        polymer = polymer2;     
    else
        if(((V2-V1)<exp(-(V2-V1)/(kb*T))))
            fprintf('\n  new bond x=%d   old bond x =%d',polymer2(1,bm+1),polymer(1,bm+1))
            Q(:,i) = (Q2);
            i = i+1;
            polymer = polymer2;    
        end
    end

    backup = backup+1;
    if(backup>10000)
        i=naccept;
    end
end

figure
plot(norm(Q))
grid off