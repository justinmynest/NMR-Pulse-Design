clear;
clc;

% so(3)
Ox=[0 0 0; 0 0 -1;0 1 0];
Oy=[0 0 1; 0 0 0;-1 0 0];
Oz=[0 -1 0; 1 0 0;0 0 0];

Ai = input('Enter Ai in rad, constant control Ai*Ix, = ');
Ay = 0;
dw0 = input('Enter chemical shift (+ve value )= ');

% initial state
X0 = [0; 0; 1];

% final state required
% Y = [0; 0; -1];                                 % for inversion
Y = [1; 0; 0];                                  % for excitation

n = input('Enter the number of flip angles, n = ');
% for i = 1:1:n
%     phase_raw(i) = input(sprintf('Flip phase of each flip %d = ',i));
% end
% 
% phase = abs(phase_raw);

% For AT = 3*pi the total flip angle divided into 10
% fprintf('delta t sampled at %d points for %s duration = \n', n, '3*pi');
ts= 10*pi;                                       % Total time for the evolution
delta = linspace(0,ts,n+1);                    % 11 points to get 10 divisions 

delta_t = ts/n;
fprintf('delta_t for the total time %s, divivded to %d divisions = %.4f \n', num2str(ts) ,n,delta_t)

% initial phase of each flip angle
fprintf('Initial phase anlges, ');
phase_init = zeros(1,n)

% The data is 

% sampling omega to omega_j
w_j = -dw0:dw0/200:dw0;

% number of samples of omega (bandwidth)
N = length(w_j); 

% generating each initial Unitarty matrices (U)
for j=1:1:N                                 % for chemcial shift
    for i = 1:1:n                           % for phase and U_i
        U(j,i,:,:) = expm(delta_t* (w_j(j)*Oz + (cos(phase_init(i)*pi/180))*Ox + (sin(phase_init(i)*pi/180))*Oy ));
    end
end

% dim of U is jxix3x3 to get 3x3 matrix for evaluation do reshaping
% example j=2,i=5 there is a 3x3 matrix, obtained as
B = reshape(U(12,5,:,:),[3,3]);

cost = 0;  cost_max =0;                % initial sample cost value, 
phase_new = phase_init;
phase_max = phase_init;
itr=1;                                      % for number of iterations, and nn for reverse order for Y_k
while itr<500 && cost <.99 
    % for finding new flip angles    
    for i=1:1:n                                     % to check at which angle
        V = zeros(3,1);
        for j=1:1:N                                 %  for each w_j
%             nn=n; Y_k = transpose(Y); X_k =X0;      % I for initial multipication, Y_(k+1) and X_(k-1) vectors
            nn=n; Y_k = Y; X_k =X0;
            for k=1:1:n                             % k for finding effective transformations for each w_j
                if k <i && k~=i
                    X_k = reshape(U(j,k,:,:),[3,3]) *X_k;               %3x1
                elseif k>i && k~=i 
%                     Y_k = Y_k *reshape(U(j,nn,:,:),[3,3]);              %1x3
                    Y_k = transpose(reshape(U(j,nn,:,:),[3,3]))*Y_k;              %1x3                    
                    nn = nn-1;                                          % nn for reverese order
                end
            end
            Vj = cross(X_k,Y_k);         % finding X_(k-1) x Y_(k+1), cross product
            V = V + Vj;                             % summing all the vectors from 
        end  
        % updating phase
        phase_new(i) = atan2(V(2),V(1))*180/pi;        % atan2(Y,X) returns values in the interval [-pi,pi]
                                                      % 180/pi converts rad to degree
        %update the U corresponding to new phase to all w_j
        for j=1:1:N
            U(j,i,:,:) = expm(delta_t*( w_j(j)*Oz + (cos(phase_new(i)*pi/180))*Ox + (sin(phase_new(i)*pi/180))*Oy ));
        end
        
        % cost function
        cost_N =0;
        for p=1:1:N
            U_t = eye(3);
            for q=1:1:n
                U_t = reshape(U(p,q,:,:),[3,3]) * U_t;
            end
            X_f = U_t*X0;
            cost_N = cost_N + dot(Y,X_f);                      % dot product
        end
        cost = cost_N/N;
        CC(itr) = cost;
        
        % storing the maximum cost condition of phase
        if cost>cost_max
            phase_max = phase_new;
            cost_max = cost
        end
        
        % breaking the loop is the condition is satisfied
        if cost>.99
%             continue          % skips any remaining statements in the body of the loop for the current iteration
%         else
            break             % break exits only from the loop in which it occurs
        end
    en
    d
    itr = itr +1
end

% error from the desired state
% error = norm(X_f-Y);                                % sqrt of sum of the components^2

% new phase angles
fprintf('The new phase angles are');
phase_new;
%phase angles for maximum cost
phase_max;

phase_new = phase_max

fprintf('Maximum cost obtained');
cost_max;

% cost over iterartions

F1 = figure;
plot(CC,'LineWidth',2)
ax=gca;                                    
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.LineWidth = 2;
grid on
title('Cost fuction over iterations','fontsize',18,'fontweight','bold')
xlabel('Iterations','fontsize',17,'fontweight','bold') 
ylabel('CostIterations','fontsize',17,'fontweight','bold')
% legend({'Profile_{}'},'Location','southeast','fontsize',16,'fontweight','bold')

% Evolution over the bandwidth

tsn = repmat(delta_t,1,n);                   % making sequence of sub sequence R, 4 times.
opts = odeset('RelTol', 1e-10);             % tolarance error for ode
l=1;                                        % for the lenght measure of total time evolution
for dw=-dw0-2:dw0/200:dw0+2
    M0=[0;0;1];                             % initial state
    m=1;                                    % for identifying flip angle 
    for T = tsn                   
        A1 = 1*Ai;    
        Ax = A1*(cos(phase_new(m)*pi/180)); 
        Ay = A1*(sin(phase_new(m)*pi/180));
        [t1, M] = ode45(@(t1,M)evolveconstant(t1, M, dw, Ax, Ay), [0, 1*T], M0, opts);
        M0 = M(end,:);                 % initilaizing the state after each flip angle
        Ys(l,m,1:3)=M(end,:);           % With different q we can get rotation after every angle  
        m=m+1;
    end
    Xs(l) = dw;                          % to get the bandwidth in same length
    l=l+1;
end  

F2 = figure;                                    % 2 for x axis , 3 for z axis
plot(Xs(:),Ys(:,n,1), 'LineWidth',2)                       
ax=gca;                                    
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.LineWidth = 2;
grid on
title('Excitation Profile with TOPS','fontsize',18,'fontweight','bold')
xlabel('Bandwidth in rad','fontsize',17,'fontweight','bold') 
ylabel('z component after different phases ','fontsize',17,'fontweight','bold')
zlabel('Excitation','fontsize',17,'fontweight','bold')
legend({'Profile_{}'},'Location','southeast','fontsize',16,'fontweight','bold')



% Evolution for R R_bar R_bar R 

phase_new_Rseq = repmat(phase_new,1,4);     % 4 times the phase angles sequence
tsnr = repmat(tsn,1,4);                     % making sequence of sub sequence R, 4 times.

l=1;                                % for the lenght measure of total time evolution
for dw=-dw0-2:dw0/200:dw0+2
    M0=[0;0;1];                         % initial state
    r=1;                            % for identifying flip angle
    for T = tsnr
        % checking which sub sequence
        if fix((r-1)/n) == 0            % 1-100 gives quotient 0 , if n =100 phase angles
            A1 = 1*Ai;
        elseif fix((r-1)/n) == 1        % 101-200 gives quotient 1 , 
            A1 = -1*Ai;
        elseif fix((r-1)/n) == 2        % 201-300 gives quotient 2 ,
            A1 = -1*Ai;
        elseif fix((r-1)/n) == 3
            A1 = 1*Ai;
        end
        % in each sub sequence the phase    
        Ax = A1*(cos(phase_new_Rseq(r)*pi/180)); 
        Ay = A1*(sin(phase_new_Rseq(r)*pi/180));
        [t1, M] = ode45(@(t1,M)evolveconstant(t1, M, dw, Ax, Ay), [0, 1*T], M0, opts);
%         ti(l) = t1 + p;              % adding evolved time
        M0 = M(end,:);                 % initilaizing the state after each flip angle
        Ys1(l,r,1:3)=M(end,:);           % With different q we can get rotation after every angle  
        r=r+1;
    end
    Xs1(l) = dw;                          % to get the bandwidth in same length
    l=l+1;
end  

F3 = figure;                                % for R R_bar R_bar R its 9x4 = 36 angles
plot(Xs1(:),Ys1(:,n,3), 'LineWidth',2)                         % when q=2, effective 90degree pulse, 3 for z component
hold on                                     % mine1, q=5 for 90 degree, mine2, q=3   
% plot(Xs1(:),Ys1(:,2*n,3))                      % inversions at 9,18,27,36 for sequence. 3 for z component
% plot(Xs1(:),Ys1(:,3*n,3))                      % inversions at 11,22,33,44 for mine1. 3 for z component
plot(Xs1(:),Ys1(:,4*n,3), 'LineWidth',2)                         % when q=36, the pulse sequence is over, 3 for z component
ax=gca;                                     % mine1 q=44/43 for 180 degree, mine2, q=7
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.LineWidth = 2;
grid on
title('Excitation Profile with TOPS Sequence 2','fontsize',18,'fontweight','bold')
xlabel('Bandwidth in rad','fontsize',17,'fontweight','bold') 
ylabel('z component after different angles ','fontsize',17,'fontweight','bold')
zlabel('Excitation','fontsize',17,'fontweight','bold')
legend({'Inversion_1','Inversion_2','Inversion_3','Inversion_4'},'Location','southeast','fontsize',16,'fontweight','bold')

phase = mod(phase_new,360);

% Writing %amplitude and phase to text file for experiment requirement. 
fid = fopen('TOPS1_exi_10pi_300n_3w.txt', 'wt');
for j = 1:1:length(phase)
    fprintf(fid, ' %3.3f,  %4.4f\n',[100 phase(j)]' );     % spacing between %f gives the same spacing in txt file
end
fclose(fid);

% Phase RF pulse
F4 = figure;
stairs(linspace(0,ts,n), phase, 'r', 'LineWidth',2);   
ax=gca;                                    
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.LineWidth = 2;
grid on
title('Piece-Wise Constant Phase Profile','fontsize',18,'fontweight','bold')
xlabel('Time "s" ','fontsize',17,'fontweight','bold') 
ylabel('Phase "°" ','fontsize',17,'fontweight','bold')
legend({'TOPS-1'},'Location','southeast','fontsize',16,'fontweight','bold')
