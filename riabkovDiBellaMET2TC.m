function riabkovDiBellaMET2TC()

clear all;

%-------------------------------------------------------------------------%
%                               Definitions                               %
%-------------------------------------------------------------------------%

% Globals

global t;     % Time vector
global dt;    % Time intervals lenghts vector
global vtac;  % Voxel time-activity curve = beta*TTAC + (1-beta)*PTAC
global ptac;  % Plasma time activity curve
global p;     % Blood input function model parameters
global alpha; % Extracellular fraction = V_E/V_Tissue
global beta;  % Tissue fraction = V_Tissue/V_voxel
global flag;  % Use blood input model ?


% Image

image = read4DDICOMFiles();
%image = ReadData3D('D:\Documents de Corentin\Images\DICOM\OsiriX Samples\CEREBRIX\PET PETCT_CTplusFET_LM_Brain (Adult)\dynamic recon 3x10min Volume (Corrected) - 7\IM-0001-0001.dcm', true);


% Rescale slope

rescaleSlope = 0.46163022518158;


% ROIs

% /!\ inversion of the y-axis compared with healthstation (y' = width-y)

targetROIBounds = [115 191 128 193 13 40]; % Tumor region where parameters are voxel-wise computed
%bloodROIBounds = [158 160 244 246 38 40]; % Blood sampling ROI, small region within the carotid if available
bloodROIBounds = [159 159 245 245 39 39];


% Signal weights

alpha = 0.5; 
beta = 0.95; 
gamma = 0.55; % Free plasmatic fraction


% Time

t = 0:10:20;
n = length(t);
dt = 10*ones(n);


%-------------------------------------------------------------------------%
%                           Blood input function                          %
%-------------------------------------------------------------------------%

% PTAC construction

ptac = zeros(n, 1);

for f=1:n
    
    mean = 0;
    count = 0;
    
    for x=bloodROIBounds(1):bloodROIBounds(2)
        
        for y=bloodROIBounds(3):bloodROIBounds(4)
            
            for z=bloodROIBounds(5):bloodROIBounds(6)
                
                mean = mean + rescaleSlope*double(image(y, x, z, f));    % /!\ x & y are inverted in Matlab !
                count = count + 1;
                
            end
        end
    end

    ptac(f) = gamma*mean/count;
    
end


% Blood input function model fit

flag = 1;

if flag == 1

    p0 = [851.1225, 21.8798, 20.8113, -4.1339, -0.1191, -0.0104];    
    p = fminsearch(@bloodCostFunction, p0);
    
end


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
 
time = 0:0.1:20;   

figure;
scatter(t, ptac, 'r');
hold on;
plot(time, interp1(t, ptac, time, 'pchip'), 'b');
hold on;
plot(time, cp(p0, time), 'g');
hold on; 
plot(time, cp(p, time), 'y');


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------%
%                             Voxel response                              %
%-------------------------------------------------------------------------%

% k constraints
 
eps = 1e-5;

%         k1  k2  k3  k4
Ainequ = [-1   0   0   0;    % k1 > 0 <=> -k1 <= -eps
           0  -1   0   0;    % k2 > 0 <=> -k2 <= -eps
           0   0  -1   0;    % k3 > 0 <=> -k3 <= -eps
           0   0   0  -1;];  % k4 > 0 <=> -k4 <= -eps
binequ = [-eps; -eps; -eps; -eps];

Aeq = [];
beq = [];


% k bounds

%kBounds = [0.0500, 0.3000;            % For 18Fluoride (vertebra) (Hawkins, 1992, p637)
%           0.0500, 0.8000;            % 1st column: lower bounds for {k1, ..., k4}
%           0.0500, 0.2000;            % 2nd column: upper bounds for {k1, ..., k4}
%           0.0001, 0.0100];
     
%kBounds = [0.1150, 0.2010;            % For FDG (Riabkov & Di Bella, 2004, p651)
%           0.0480, 0.2850;            % 1st column: lower bounds for {k1, ..., k4}
%           0.0100, 0.0960;            % 2nd column: upper bounds for {k1, ..., k4}
%           0.0028, 0.0098];

kBounds = [0.0000, 1.0000;            % Why does k parameters belong tp [0, 1] ? 
           0.0000, 1.0000;            
           0.0000, 1.0000;            
           0.0000, 1.0000];

lowerBounds = kBounds(:, 1);
upperBounds = kBounds(:, 2);


% fmincon setup

constraintTolerance = eps/2;
maxIterations = 150;
optimalityTolerance = 1e-10;                % default 1e-6
stepTolerance = 1e-6;                       % default 1e-6
k0 = 0.5*(lowerBounds + upperBounds);


options = optimoptions('fmincon', ...                                     % Parameter update between two iterations must be adapted since k parameters precision must be up to 10^-4
                       'Algorithm', 'sqp', ...
                       'CheckGradients', false, ...
                       'ConstraintTolerance', constraintTolerance, ...
                       'Diagnostics', 'off', ...
                       'Display', 'iter', ...
                       'MaxFunctionEvaluations', 2000, ...
                       'MaxIterations', maxIterations, ...
                       'OptimalityTolerance', optimalityTolerance, ...
                       'SpecifyConstraintGradient', false, ...           
                       'SpecifyObjectiveGradient', false, ...             % May be accelerated using analytic derivatives
                       'StepTolerance', stepTolerance, ...
                       'TypicalX', k0);
                   

%-------------------------------------------------------------------------%               
%-------------------------------------------------------------------------%               


x = 151;
y = 148;
z = 22;

vtac = zeros(n, 1);

for f=1:n
    
    vtac(f) = rescaleSlope*double(image(y, x, z, f));    % /!\ x & y are inverted in Matlab !
    
end

figure;
plot(time, interp1(t, vtac, time, 'pchip'));

return;

[k, finalCost, exitFlag, info] = fmincon(@voxelCostFunction, k0, Ainequ, binequ, Aeq, beq, lowerBounds, upperBounds);     

H = convolutionMatrix(kInit);

signal = beta*(H*ptac) + (1-beta)*ptac;
figure, plot(time2, signal);


H = convolutionMatrix(k);

signal = beta*(H*ptac) + (1-beta)*ptac;
figure, plot(time2, signal);
                        
return;   


%-------------------------------------------------------------------------%               
%-------------------------------------------------------------------------%               


% Coefficients computation

vtac = zeros(n, 1);
map = zeros(targetROIBounds(2)-targetROIBounds(1)+1, targetROIBounds(4)-targetROIBounds(3)+1, targetROIBounds(6)-targetROIBounds(5)+1, 4);

for x=targetROIBounds(1):targetROIBounds(2)
    
    for y=targetROIBounds(3):targetROIBounds(4)
        
        for z=targetROIBounds(5):targetROIBounds(6)
            
            for t=1:n
                
                vtac(t) = rescaleSlope*double(image(y, x, z, t));    % /!\ x & y are inverted in Matlab !
                
            end
            
            [k, finalCost, exitFlag, info] = fmincon(@voxelCostFunction, k0, Ainequ, binequ, Aeq, beq, lowerBounds, upperBounds, [], options);
            map(x-targetROIBounds(1)+1, y-targetROIBounds(3)+1, z-targetROIBounds(5)+1, :) = k;    % No x-y inversion for the k-map
                
        end
    end
end

return;

end

function cost = bloodCostFunction(p)

global t;
global ptac;

cost = norm(cp(p, t)-ptac)^2;

return;

end

function cost = voxelCostFunction(k)

global t;
global dt;
global vtac;
global ptac;
global p;
global alpha;
global beta;
global flag;

if flag == 1
    
    H = convolutionMatrix(k, alpha, t, dt);
    cost = norm(vtac-(beta*(H*ptac) + (1-beta)*ptac))^2; % Add weight ponderation => Feng, 1997
                                                         % /!\ vtac value is actually the mean tracer concentration observed during the frame, so should be the right term beta*(H*b) + (1-beta)*b
                                                         % => Take this into account during frame reconstruction and vtac extraction. Perhaps reconstruct smaller frames so that mean can be
                                                         % considered to be equal to the computed value ?
                                                         
else
    
    cost = norm(vtac-(beta*ct(k, p, alpha, t) + (1-beta)*cp(p, t)))^2;
    
end
    
return;

end

function H = convolutionMatrix(k, alpha, t, dt)

% Convolution matrix for parameters 'k', extracellular fraction 'alpha' and time vector 't' and time interval lenght vector 'dt'
% Riabkov & Di Bella, 2002, equation (9) extended to 2TC model

row = [dt*h(k, alpha, 0), zeros(1, lenght(t)-1)]; % dt is used to ensure that discrete convolutions equals continuous convolution (Riabkov & Di Bella, 2002 & 2004) => to be verified
                                                  % Could be simply adapted to variable dt
col = dt*h(k, alpha, t);

H = toeplitz(col, row);

return;

end

function cp = cp(p, t)
  
% Blood input function for parameters 'p' at time 't'
% Feng, 1997, equation (1)

%      1  2  3     4       5       6
% p    A1 A2 A3 lambda1 lambda2 lambda3 

cp = (p(1)*t-p(2)-p(3)).*exp(p(4)*t) + p(2)*exp(p(5)*t) + p(3)*exp(p(6)*t);

return;

end

function h = h(k, alpha, t)

% 2TC tissue impulse response function for parameters 'k' and extracellular fraction 'alpha' at time 't'
% Feng, 1997, equation (5) adapted to the intra/extracellular model

q = sqrt((k(2)+k(3)+k(4))^2 - 4*k(2)*k(4));
l1 = 0.5*(k(2)+k(3)+k(4) - q);
l2 = 0.5*(k(2)+k(3)+k(4) + q);
b1 = k(1)/(l2-l1) * ((1-alpha)*k(3) + alpha*(k(4)-l1));
b2 = k(1)/(l2-l1) * (alpha*(l2-k(4)) - (1-alpha)*k(3));

h = b1*exp(-l1*t) + b2*exp(-l2*t); 

return;

end

function ct = ct(k, p, alpha, t)

% 2TC analytic tissue response to blood input function for tissue parameters 'k', blood input parameters 'p' and extracellular fraction 'alpha' at time 't'
% Loeb, 2015, equation (15b) adapted to the intra/extracellular model

%      1  2  3     4       5       6
% p    A1 A2 A3 lambda1 lambda2 lambda3 

%      1  2  3  4
% k    k1 k2 k3 k4 

a = p(1:3);
lambda = p(4:6);

q = sqrt((k(2)+k(3)+k(4))^2 - 4*k(2)*k(4));
l1 = 0.5*(k(2)+k(3)+k(4) - q);
l2 = 0.5*(k(2)+k(3)+k(4) + q);
b1 = k(1)/(l2-l1) * ((1-alpha)*k(3) + alpha*(k(4)-l1));
b2 = k(1)/(l2-l1) * (alpha*(l2-k(4)) - (1-alpha)*k(3));

b = [b1, b2];
l = [l1, l2];

c = [-a(2)-a(3)-(a(1)/(b(1)-lambda(1))), a(2), a(3);
     -a(2)-a(3)-(a(1)/(b(2)-lambda(1))), a(2), a(3)];

ct = zeros(1, size(t));

for i=1:2
    
    ct = ct + ((a(1)*t*b(i))/(l(i)-lambda(1)) .* exp(-lambda(1)*t));
    
    for j=1:3
        
        ct = ct + (c(i, j)*b(i)/(l(i)-lambda(j)) * (exp(-lambda(j)*t)-exp(-l(i)*t)));
        
    end;
    
end;


return;

end