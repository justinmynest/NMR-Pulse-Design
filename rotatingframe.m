clear;
clc;
%Rotating Frame NMR Spectroscopy

n = 1; A1 = []; C1 = []; V1 = [];

A1(n) = input('Enter required amplitude of RF field (Hz)in the final frame (1/4*j. j=1,2,...), An = ');
% with Vn = 1
C1(n) = input('Enter max value chemical shift required in final frame (<< 1, eg = .01), Cn = ');
V1(n) = 1; % Setting the midpoint of effective chemical shift in second last frame (1 for Hz, 1000 for kHz)
c = input('Enter the max value actual chemical shift (without normalization), [1.2, 4.4] for MODE 4 (3 frames), [12, 27.6] for MODE 6, [59.6, 123.4] for MODE 8 \n [251.5, 507.4] for MODE 10, [1020, 2043] for MODE 12 , 5000 for MODE 14, 20000 for MODE 16\n C0 = ');                       

% [1.2, 4.4] for MODE 4 (3 frames), [12, 27.6] for MODE 6, [59.6, 123.4] for MODE 8, 
% If the chemical shift in the final frame is reduced we can reduce the no of frames

while c > C1(n)
    
    V1(n+1) = V1(n)*2;
    C1(n+1) = 2*sqrt(C1(n)*V1(n));
    A1(n+1) = V1(n) - C1(n);
    
    n = n+1;                % for next iteration and counting number of frames
end

fprintf('Number of frames (n+1) need for a bandwith of [-%f %f] = %d \n', C1(n), C1(n), n);
%n+1 = MODE#
% Reversing indexing order 

m = n; A = []; C = []; V = [];
while m>0
    A(n-m+1) = A1(m)/A1(n);                 % Normalizing everything by A0 = A1(n)
    C(n-m+1) = C1(m)/A1(n);
    V(n-m+1) = V1(m)/A1(n);
    m = m-1;
end

% RF Field Generation

T = pi/(2*A(n));                            % For Excitation
% T = pi/(1*A(n));                            % For Inversion
% T = 1*pi/(4*A(n));                            % For other flip angles pi/4, As we make random flip,  
% T = 1*pi/(6*A(n));                          % AT should be multiple of 2pi, hence here An=1/12
fprintf('Simulation Time = %f \n', T);

sympref('FloatingPointOutput', true);     % for displaying output in decimal values
syms t

u = 1; w = 0;

if n>2
    i = n; j = n;
    while j>1
        
        u = 2*cos(V(j)*t)*u;
        j = j-1;
    end
    u = A(i)*u;
    i = i-1;    
    while i>1        
        v = 1; k = i-1;        
        while k>1            
            v = v*2*cos(V(k)*t);
            k = k-1;
        end
        v = v*2*A(i)*sin(V(i)*t);
        w = w + v;
        i = i-1;
    end
    w = w + u;
else
    w = 2*A(2)*cos(V(2)*t);
end

fprintf('The RF field required is')
RF =[A(1), w]                           % Combined x and y rf fields\

% % Actual RF pulse
% f1 = figure;
% fplot(RF(2), [0 T]);                  % plots a symbolic expression
% grid on
% hold on
% fplot(RF(1), [0 T]);
% title('Actual RF Field')
% xlabel('Time "sec" ') 
% ylabel('Amplitude "rad" ')
% legend({'Oy','Ox'},'Location','southeast')

% % Absolute RF pulse
% f2 =figure;
% fplot(sqrt(RF(1)^2 + RF(2)^2), [0 T]);
% grid on
% title('Actual RMS RF Field')
% xlabel('Time "sec" ') 
% ylabel('Amplitude "rad" ')
% legend({'Amplitude'},'Location','southeast')

% % Phase RF pulse
% f3 =figure;
% fplot(angle(RF(2)/RF(1)), [0 T]);
% grid on
% title('Actual Phase Profile')
% xlabel('Time "sec" ') 
% ylabel('Phase "rad" ')
% legend({'Phase'},'Location','southeast')


% RMS value over the simulation time
fun = RF(1)^2 + RF(2)^2;                    % Integration in symbolic form
% q = int(fun,0,T);
% rms = sqrt(q/T);
% fprintf('RMS value of the field = %f \n', rms);
ht = matlabFunction(fun);                   % Convert symbolic expression to function handle
Ar = sqrt(integral(ht,0,T)/T);
Ar2 = Ar/(2*pi);                            % RMS/2pi = radians to hertz   (all calculations are done in radians )
fprintf('Operating RMS frequency, Ar2 = %f (Hz) \n', Ar2);

% Coefficients for higer power and shorter time
Ar1 = input('Enter required operating RMS frequency in Hz(no 2pi required), Arms = ');
P = Ar1/Ar2;
fprintf('Multiplying factor = %f \n', P);
Ap = A*P;
Vp = V*P;
Cp = C*P;
Tp = T/P;

% Control of spin

l=1;
for dw=-Cp(1):Cp(1)/100:Cp(1)

    M0=[0;0;1];
    opts = odeset('RelTol', 1e-10);
    [t1, M] = ode45(@(t1,M)rfevolve(t1, M, dw, Ap, Vp, n), [0, 1*Tp], M0, opts);
    X(l,1:3)=M(end,:);
    bw(l,1) = dw;
    l=l+1;

end

% 
% f4 = figure;
% plot(t1, M(:, 3))
% grid on

% Final positions of different omega along x direction
f5 = figure;
plot(bw(:,1)/(2*pi*1e3), X(:,1), 'LineWidth',2)               % converted rad to Hz by (1/2pi)
ax=gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
ax.LineWidth = 2;
grid on
title('Excitation Profile with MODE n','fontsize',16,'fontweight','bold')
xlabel('Offset (kHz) ','fontsize',16,'fontweight','bold')
ylabel('X component magnetization','fontsize',16,'fontweight','bold')
%legend({'M_x'},'Location','southeast','fontsize',16,'fontweight','bold')


u = 1; w = 0;

if n>2
    i = n; j = n;
    while j>1
        
        u = 2*cos(Vp(j)*t)*u;
        j = j-1;
    end
    u = Ap(i)*u;
    i = i-1;    
    while i>1        
        v = 1; k = i-1;        
        while k>1            
            v = v*2*cos(Vp(k)*t);
            k = k-1;
        end
        v = v*2*Ap(i)*sin(Vp(i)*t);
        w = w + v;
        i = i-1;
    end
    w = w + u;
else
    w = 2*Ap(2)*cos(Vp(2)*t);
end
RF =[Ap(1), w] 

Ra = (sqrt(RF(1)^2 + RF(2)^2))/(2*pi);                  % converted rad to Hz by (1/2pi)
Rp = mod(atan2(RF(2),RF(1))*(180/pi),360);              % rad to degree and then mod(theta,360)

% Absolute RF pulse
f6 =figure;
fplot(Ra, [0 Tp]);            
grid on
title('Fixed Amplitude RF Field')
xlabel('Time "s" ') 
ylabel('Amplitude "Hz" ')
legend({'Amplitude Profile'},'Location','southeast')

% Phase RF pulse
f7 = figure;
fplot(Rp, [0 Tp]);           
grid on
title('Fixed Phase Profile')
xlabel('Time "s" ') 
ylabel('Phase "°" ')
legend({'Phase Profile'},'Location','southeast')

f8 = figure;
subplot(2,1,1);
fplot(Ra/(1e3), [0 Tp],'LineWidth',2); 
ax=gca;
% ax.XAxis.Exponent = -3;
% xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'})      % for MODE4
% xticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5'})      % for MODE8
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
ax.LineWidth = 2;
grid on
title('RF Amplitude Profile','fontsize',16,'fontweight','bold')
xlabel('Time (s)','fontsize',16,'fontweight','bold') 
ylabel('Amplitude (kHz)','fontsize',16,'fontweight','bold')
%legend({'Amplitude Profile'},'Location','southeast','fontsize',16,'fontweight','bold')

subplot(2,1,2);
fplot(Rp, [0 Tp],'LineWidth',2);      
ax=gca;
% xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'})
% xticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5'})      % for MODE8
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
ax.LineWidth = 2;
grid on
title('RF Phase Profile','fontsize',16,'fontweight','bold')
xlabel('Time (s)','fontsize',16,'fontweight','bold') 
ylabel('Phase (°)','fontsize',16,'fontweight','bold')
%legend({'Phase Profile'},'Location','southeast','fontsize',16,'fontweight','bold')

Tsamp = (20*Vp(2)/(2*pi))^-1;          % rad to Hz to find sampling frequency
i=1;
for t = 0:Tsamp:Tp
    tr1(i) = t;
    Ra1(i) = double(subs(Ra));          % substitutes the symbolic variable value in the expression
    Rp1(i) = double(subs(Rp));          % converting symbolic expression to double precision
    i = i+1;
end
    
Ramax =  max(Ra1);
fprintf('Maximum RMS Amplitude of RF = %f \n', Ramax);
fprintf('Pulse duration = %f \n', Tp);

% Writing %amplitude and phase to text file for experiment requirement. 
% fid = fopen('MODE_14_10x.txt', 'wt');
% for j = 1:1:length(Ra1)
%     fprintf(fid, ' %3.3f,  %4.4f\n',[(Ra1(j)/Ramax)*100 Rp1(j)]' );     % spacing between %f gives the same spacing in txt file
% end
% fclose(fid);

% Alternative way to find rms amplitude of the RF pulse
Apr = 0; j=1;
for i = 1:n
    Apr = Apr + j*Ap(i)^2;
    j = j*2;
end
fprintf('RMS Amplitude of RF in Hz = %f \n', sqrt(Apr)/(2*pi));

% assume(t>0 & t<Tp)
% g = diff(Ra, t);
% solve(g == 0, t, 'MaxDegree', 4);    %option specifies the maximum degree of polynomials for which the solver tries to return explicit solutions
% extrema = vpa(ans, 6)    %approximate the exact solution numerically by using the vpa function