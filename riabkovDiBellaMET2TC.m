function riabkovDiBellaMET2TC()

clear all;


% Globals definition

global n;
global dt;
global tac;
global b;
global alpha;
global beta;
global gamma;


% ROIs defition

% /!\ inversion of the y-axis compared with healthstation (y' = width-y)

targetROIBounds = [115 191 128 193 13 40]; % Tumor region where parameters are voxel-wise computed
%bloodROIBounds = [158 160 244 246 38 40]; % Blood sampling ROI, small region within the carotid if available
bloodROIBounds = [159 159 245 245 39 39];


% 4D image reading

[image, n] = read4DDICOMFiles();
%image = ReadData3D('D:\Documents de Corentin\Images\DICOM\OsiriX Samples\CEREBRIX\PET PETCT_CTplusFET_LM_Brain (Adult)\dynamic recon 3x10min Volume (Corrected) - 7\IM-0001-0001.dcm', true);


% Rescale slope definition

rescaleSlope = 0.46163022518158;


% Signal weights

alpha = 0.5; % Extracellular fraction = V_E/V_Tissue
beta = 0.95; % Tissue fraction = V_Tissue/V_voxel
gamma = 0.55; % Free plasmatic fraction


% Linear time interval

dt = 0.1; % Min


% Image blood sampling

b = zeros(n, 1);

for t=1:n
    
    mean = 0;
    number = 0;
    
    for x=bloodROIBounds(1):bloodROIBounds(2)
        
        for y=bloodROIBounds(3):bloodROIBounds(4)
            
            for z=bloodROIBounds(5):bloodROIBounds(6)
                
                mean = mean + rescaleSlope*double(image(y, x, z, t));    % /!\ x & y are inverted in Matlab !
                number = number + 1;
                
            end
        end
    end

    b(t) = gamma*mean/number;
    
end


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

kBounds = [0.0, 1.0;            % For FDG (Riabkov & Di Bella, 2004, p651)
           0.0, 1.0;            % 1st column: lower bounds for {k1, ..., k4}
           0.0, 1.0;            % 2nd column: upper bounds for {k1, ..., k4}
           0.0, 1.0];

lowerBounds = kBounds(:, 1);
upperBounds = kBounds(:, 2);


% fmincon setup

constraintTolerance = eps/2;
maxIterations = 150;
optimalityTolerance = 1e-10;                % default 1e-6
stepTolerance = 1e-6;                       % default 1e-6
kInit = 0.5*(lowerBounds + upperBounds);


options = optimoptions('fmincon', ...
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
                       'TypicalX', kInit);
                   

%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%               
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%               
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%               

                   
time = 0:10:20; % Input time
time2 = 0:dt:35; % Resampled time


b = transpose(interp1(time, b, time2, 'pchip'));
figure, plot(time2, b);


x = 145;
y = 150;
z = 24;

tac = zeros(n, 1);

for t=1:n
    tac(t) = rescaleSlope*double(image(y, x, z, t));    % /!\ x & y are inverted in Matlab !
end

tac = transpose(interp1(time, tac, time2, 'pchip'));
figure, plot(time2, tac);    


n = length(time2);


[k, finalCost, exitFlag, info] = fmincon(@costFunction, kInit, Ainequ, binequ, Aeq, beq, lowerBounds, upperBounds);     

kInit
k

H = convolutionMatrix(kInit);

signal = beta*(H*b) + (1-beta)*b;
figure, plot(time2, signal);


H = convolutionMatrix(k);

signal = beta*(H*b) + (1-beta)*b;
figure, plot(time2, signal);
              
                   
return;   


%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%               
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%               
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%               


% Coefficients computation

tac = zeros(n, 1);
map = zeros(targetROIBounds(2)-targetROIBounds(1)+1, targetROIBounds(4)-targetROIBounds(3)+1, targetROIBounds(6)-targetROIBounds(5)+1, 4);

for x=targetROIBounds(1):targetROIBounds(2)
    
    for y=targetROIBounds(3):targetROIBounds(4)
        
        for z=targetROIBounds(5):targetROIBounds(6)
            
            for t=1:n
                
                tac(t) = rescaleSlope*double(image(y, x, z, t));    % /!\ x & y are inverted in Matlab !
                
            end
            
            [k, finalCost, exitFlag, info] = fmincon(@costFunction, kInit, Ainequ, binequ, Aeq, beq, lowerBounds, upperBounds, [], options);
            map(x-targetROIBounds(1)+1, y-targetROIBounds(3)+1, z-targetROIBounds(5)+1, :) = kuvs2ki(kuvs, dt);    % No x-y inversion for the k-map
                
        end
    end
end

return;

function cost = costFunction(k)

global tac;
global b;
global beta;

H = convolutionMatrix(k);

cost = norm(tac-(beta*(H*b) + (1-beta)*b))^2; % Add weight ponderation => Feng, 1997

return;

function H = convolutionMatrix(k)

% Convolution matrix
% Riabkov & Di Bella 2002, equation (9) extended to 2TC model

global n;
global dt;

row = [dt*tissueResponse(k, 0) zeros(1, n-1)];
col = zeros(1, n);

for i=1:n
    col(i) = dt*tissueResponse(k, double(i-1)*dt); % i belongs to [0, ..., n-1]
end

H = toeplitz(col, row);

return

function h = tissueResponse(k, t)

% 2TC tissue response for parameters k at time t

global alpha;

q = sqrt((k(2)+k(3)+k(4))^2 - 4*k(2)*k(4));
a1 = 0.5*(k(2)+k(3)+k(4) - q);
a2 = 0.5*(k(2)+k(3)+k(4) + q);

h = k(1)/(a2-a1) * (((1-alpha)*k(3) + alpha*(k(4)-a1))*exp(-a1*t) + (alpha*(a2-k(4)) - (1-alpha)*k(3))*exp(-a2*t)); 

return;