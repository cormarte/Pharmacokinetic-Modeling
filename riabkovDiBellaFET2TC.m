function riabkovDiBellaFET2TC()

clear all;


% Globals definition

global n;
global dt;
global tac;
global b;


% Rescale slope definition

rescaleSlope = 0.46163022518158;


% ROIs defition

% /!\ inversion of the y-axis compared with healthstation (y' = width-y)

targetROIBounds = [115 191 128 193 13 40]; % Tumor region where parameters are computed
bloodROIBounds = [158 160 244 246 38 40]; % Blood sampling ROI, small region within the carotid if available
%bloodROIBounds = [159 159 245 245 39 39];


% Linear time interval

dt = 10; 


% 4D image reading

[image, n] = read4DDICOMFiles();
%image = ReadData3D('D:\Documents de Corentin\Images\DICOM\OsiriX Samples\CEREBRIX\PET PETCT_CTplusFET_LM_Brain (Adult)\dynamic recon 3x10min Volume (Corrected) - 7\IM-0001-0001.dcm', true);

% Image blood sampling

b = zeros(n, 1);
alpha = 1.0; % Free plasmatic fraction

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

    b(t) = alpha*mean/number;
    
end

time = 0:10:20;
time2 = 0:0.1:20;
b2 = interp1(time, b, time2, 'pchip');
figure, plot(time2, b2);


% k u v s constraints
 
eps = 1e-5;

%          k  u  v  s
Ainequ = [-1  0  0  0;    % k > 0 <=>   -k <=  -eps
           0  1  0  0;    % u < 1 <=>    u <= 1-eps
           0 -1  0  1;    % u > s <=>  s-u <=  -eps
           0  0  1 -1;    % s > v <=>  v-s <=  -eps
           0  0 -1  0];   % v > 0 <=>   -v <=  -eps
binequ = [-eps; 1-eps; -eps; -eps; -eps];

Aeq = [];
beq = [];


% k u v s bounds

%kiBounds = [0.0500, 0.3000;            % For 18Fluoride (vertebra) (Hawkins, 1992, p637)
%            0.0500, 0.8000;            % 1st column: lower bounds for {k1, ..., k4}
%            0.0500, 0.2000;            % 2nd column: upper bounds for {k1, ..., k4}
%            0.0001, 0.0100];
     
kiBounds = [0.1150, 0.2010;             % For FDG (Riabkov & Di Bella, 2004, p651)
            0.0480, 0.2850;             % 1st column: lower bounds for {k1, ..., k4}
            0.0100, 0.0960;             % 2nd column: upper bounds for {k1, ..., k4}
            0.0028, 0.0098];

lowerBounds = ki2kuvs(kiBounds(:,1)*(1+eps), dt);
upperBounds = ki2kuvs(kiBounds(:,2)*(1-eps), dt);

sortedBounds = sort([lowerBounds upperBounds], 2);

lowerBounds = sortedBounds(:, 1);
upperBounds = sortedBounds(:, 2);


% fmincon setup

constraintTolerance = eps/2;
maxIterations = 150;
optimalityTolerance = 1e-10;                % default 1e-6
stepTolerance = 1e-6;                       % default 1e-6
kuvsInit = 0.5*(lowerBounds + upperBounds);

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
                       'TypicalX', kuvsInit);

               
% Coefficients computation

tac = zeros(n, 1);
map = zeros(targetROIBounds(2)-targetROIBounds(1)+1, targetROIBounds(4)-targetROIBounds(3)+1, targetROIBounds(6)-targetROIBounds(5)+1, 4);

for x=targetROIBounds(1):targetROIBounds(2)
    
    for y=targetROIBounds(3):targetROIBounds(4)
        
        for z=targetROIBounds(5):targetROIBounds(6)
            
            for t=1:n
                
                tac(t) = rescaleSlope*double(image(y, x, z, t));    % /!\ x & y are inverted in Matlab !
                
            end
            
            [kuvs, finalCost, exitFlag, info] = fmincon(@costFunction, kuvsInit, Ainequ, binequ, Aeq, beq, lowerBounds, upperBounds, [], options);
            map(x-targetROIBounds(1)+1, y-targetROIBounds(3)+1, z-targetROIBounds(5)+1, :) = kuvs2ki(kuvs, dt);    % No x-y inversion for the k-map
                
        end
    end
end

figure, imshow(map(:,:,3,1));
figure, imshow(map(:,:,3,2));
figure, imshow(map(:,:,3,3));
figure, imshow(map(:,:,3,4));

return;

function ki = kuvs2ki(kuvs, dt)

a1 = -log(kuvs(2))/dt;
a2 = -log(kuvs(3))/dt;
k1 = kuvs(1);
k2 = (-kuvs(4)*(a2-a1) + kuvs(3)*a2 - kuvs(2)*a1)/(kuvs(3)-kuvs(2));
k4 = a1*a2/k2;
k3 = -k2 + (a1+a2) - k4;

ki = [k1; k2; k3; k4];

return;

function kuvs = ki2kuvs(ki, dt)

q = sqrt((ki(2)+ki(3)+ki(4))^2-4*ki(2)*ki(4));
a1 = 0.5*(ki(2)+ki(3)+ki(4)-q);
a2 = 0.5*(ki(2)+ki(3)+ki(4)+q);
B1 = ki(1)/q * (ki(3)+ki(4)-a1);
B2 = ki(1)/q * (ki(3)+ki(4)-a2);
k = ki(1);
u = exp(-a1*dt);
v = exp(-a2*dt);
s = (B1*v-B2*u)/k;

kuvs = [k; u; v; s];

return;

function cost = costFunction(kuvs)

global n;
global dt;
global tac;
global b;

H = linearConvolutionMatrix(kuvs, dt, n);

cost = norm(tac-H*b)^2;

return;

function H = linearConvolutionMatrix(kuvs, dt, numberOfFrames)

% Convolution matrix
% Riabkov & Di Bella 2002, equation (9) extended to 2TC model

row = [tissueResponse(kuvs, dt, 0) zeros(1, numberOfFrames-1)];
col = zeros(1, numberOfFrames);

tissueResponse(kuvs, dt, 0);

for i=1:numberOfFrames
    col(i) = tissueResponse(kuvs, dt, double(i-1)); % i belongs to [0, ..., numberOfFrames-1]
end

H = toeplitz(col, row);

return

function h = tissueResponse(kuvs, dt, n)

% 2TC tissue response for parameters kuvs at time n*dt
% Riabkov & Di Bella 2004, equation (38)

h = (dt*kuvs(1)*((kuvs(2)-kuvs(4))*kuvs(2)^n + (kuvs(4)-kuvs(3))*kuvs(3)^n))/(kuvs(2)-kuvs(3));

return;

function dhdk = dTissueResponsedk(kuvs, dt, n)

% Partial derivative of tissueResponse w.r.t k

dhdk = h(kuvs, dt, n)/kuvs(1);

return;

function dhdu = dTissueResponsedu(kuvs, dt, n)

% Partial derivative of tissueResponse w.r.t u

dhdu = (dt*kuvs(1)*kuvs(2)^(n-1)*((n+1)*kuvs(2)-n*kuvs(4)) - h(kuvs, dt, n))/(kuvs(2)-kuvs(3));

return;

function dhdv = dTissureResponsedv(kuvs, dt, n)

% Partial derivative of tissueResponse w.r.t v

dhdv = (dt*kuvs(1)*kuvs(3)^(n-1)*(n*kuvs(4)-(n+1)*kuvs(3)) + h(kuvs, dt, n))/(kuvs(2)-kuvs(3));

return;

function dhds = dTissureResponseds(kuvs, dt, n)

% Partial derivative of tissueResponse w.r.t s

dhds = (dt*kuvs(1)*(kuvs(3)^n-kuvs(2)^n))/(kuvs(2)-kuvs(3));

return;


