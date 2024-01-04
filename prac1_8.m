%% Import data
clearvars; close all;
tic
T = 1000;
N = 100;
b = 3.0;
L = T*(N+2);

% need to substitute file address here. for some reason the sprintf method
% is not working for me, but it is a work in progress
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
toc

%% mont - e - carlo
% set up molecule
tic
Tmax= 10;         % number of starting configurations to average over - kept lower for initial testing
n = 100;         % number of bonds used in simulation

naccept = 1000;

Temp=1;            % temp in planck units
kb=1;              % boltzmann constant in planck units

Fx = (0.0:0.1:5);                  % kept lower for initial testing
% Fx = 10;
ext = zeros(1,Tmax);               % mean extension for given force over Tmax configurations
Qxnum = zeros(1,length(Fx));       % mean of ext for a given force

% monte the carlo
for f = 1:length(Fx)    % vary force applied
    F = [Fx(f) 0 0];
    for t = 1:Tmax      % vary starting configs
        polymer = [xdata(t,1:n); ydata(t,1:n); zdata(t,1:n)];   % read starting config from data
        i  = 1;         % reset MC loop iteration
        count=1;        % sanity check variable
        Q = zeros(1,naccept);
        while i<=naccept   % MC loop
            polymer2 = polymer;             % candidate
            bm = floor(rand()*(n-2)) + 2;   % random number to select bond to modify
            newx = sqrt(b^2)*rand();
            newy = sqrt((b^2 - newx^2))*rand();
            newz = sqrt(b^2 - (newx^2 + newy^2));
            newbond = [newx; newy; newz];
            % newbond = randn(3,1);
            % newbond = b* newbond ./ norm(newbond);
            % 
            % approach 1
            % calculate small r
            % incrs = zeros(3,n);           % bond vectors relative to previous monomer instead of origin
            % for j = bm:n-1
            %     incrs(:,j+1) = polymer(:,j+1) - polymer(:,j);
            % end
            % 
            % % update bond at bm+1
            % polymer2(:,bm+1) = polymer(:,bm) + newbond(:);  
            % 
            % % update rest of molecule 
            % for j = bm+1:n-1
            %     polymer2(:,j+1) = polymer2(:,j) + incrs(:,j+1);
            % end

            % approach 2
            % update molecule from bm to end by removing old bond and
            % adding new bond at index bm
            oldbond = polymer2(:,bm+1) - polymer2(:,bm);
            polymer2(:,bm+1:end) = polymer(:,bm+1:end) + newbond - oldbond;

            % pass/fail parameters
            Q1 = [(polymer(1,end)-polymer(1,1)); (polymer(2,end)-polymer(2,1)); (polymer(3,end)-polymer(3,1))];
            Q2 = [(polymer2(1,end)-polymer2(1,1)); (polymer2(2,end)-polymer2(2,1)); (polymer2(3,end)-polymer2(3,1))];
            V1 = dot(-F,Q1);
            V2 = dot(-F,Q2);
            dv = V2-V1;

            % pass/fail check
            if(dv<=0)
                Q(i) = Q2(1);
                i = i+1;
                polymer = polymer2;     
            elseif(randn() <= exp(-dv/(kb*Temp)))
                Q(i) = Q2(1);
                i = i+1;
                polymer = polymer2; 
            end
        
            % sanity check
            count = count+1;
            if(count>1e6)
                warning('         Terminating due to sanity check')
                i=naccept+1;
            end
        end
            ext(t) = mean(Q(end-200:end));        
    end
    Qxnum(f) = mean(ext);
end
toc

% optional outputs
% fprintf('\n simulation run %d times',count)
% fprintf('\n rms of last 200 values = %d\n',rms(Q(end-200:end)))
% fprintf('\nAverage extension at temperature  %d over %d configurations = %d\n',Temp,Tmax,Qxnum)
 
figure
plot(Q)
grid on
title(sprintf('Q for Fx=%d',Fx(f)))
% fprintf('\nforce = %d avg of last 100 Q values = %d\n',Fx,mean(Q(end-100:end)))

% % theoretical results
alpha = Fx*1/(kb*Temp);
Qxmdl = 1 * (coth(alpha)-(1./alpha));

titlestr = sprintf('Force-Extension Monte Carlo n=%d, over %d conformations',n,Tmax);

figure
plot(Fx,Qxmdl);
grid on
hold on
plot(Fx,Qxnum/(n*b),'o','MarkerEdgeColor','red');
title(titlestr,'FontSize',20)
xlabel('force applied (Planck units)','FontSize',15)
ylabel('extension','FontSize',15)