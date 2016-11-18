%MBF myocardial blood flow estimation with kinetic modeling
%   [ Y, BP_TAC, MYO_TAC ] = MBF(D,E,T,U,VOI1,VOI2) estimates the myocardial
%      blood flow from a dynamic data set. The input variables are:
%       D: 4D data set
%       E: 4D error data set (maybe empty)
%       T: timing information for data and error
%       U: data and error unit (see Voxulus list)
%       VOI1: volume of interest for input function estimation (maybe empty)
%       VOI2: volume of interest for target region
%      The output variable are:
%       Y: MBF and associated error
%       BP_TAC: blood pool time activity curve
%       MYO_TAC: target region time activity curve
%      
%      The indexing order of the data is Z,Y,X,T and the for the VOIs Z,Y,X.
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%

function [ y, bp_tac, myo_tac ] = mbf (d, e, t, u, voi1, voi2)
  bp_tac = [];
  
  tracer = voxulus ('get', 'current', 'tracer');
  if strcmp (tracer, 'N-13')
    model = 'Cardiac_HutchinsII';
    voxulus ('select', 'model', model);
    voxulus ('set', 'model', 'k1', 'value', 1.0);
    voxulus ('set', 'model', 'k1', 'initvalue', 1.0);
    voxulus ('set', 'model', 'k1', 'min', 0.001);
    voxulus ('set', 'model', 'k1', 'max', 10.0);
    voxulus ('set', 'model', 'k2', 'value', 1.0);
    voxulus ('set', 'model', 'k2', 'initvalue', 1.0);
    voxulus ('set', 'model', 'k2', 'min', 0.001);
    voxulus ('set', 'model', 'k2', 'max', 10.0);
    voxulus ('set', 'model', 'k3', 'value', 0.0);
    voxulus ('set', 'model', 'k3', 'initvalue', 0.05);
    voxulus ('set', 'model', 'k3', 'min', 0.001);
    voxulus ('set', 'model', 'k3', 'max', 1.0);
    voxulus ('set', 'model', 'k3', 'state', 'const');
    voxulus ('set', 'model', 'tbv', 'value', 0.25);
    voxulus ('set', 'model', 'tbv', 'initvalue', 0.25);
    voxulus ('set', 'model', 'tbv', 'min', 0.001);
    voxulus ('set', 'model', 'tbv', 'max', 1.0);
    %voxulus ('set', 'model', 'tbv', 'state', 'const');
    voxulus ('set', 'model', 'unit', u);
  elseif strcmp (tracer, 'Rb-82')
    %model = 'Cardiac_2TC';
    model = 'Cardiac_HutchinsII';
    voxulus ('select', 'model', model);
    voxulus ('set', 'model', 'k1', 'value', 1.0);
    voxulus ('set', 'model', 'k1', 'initvalue', 1.0);
    voxulus ('set', 'model', 'k1', 'min', 0.001);
    voxulus ('set', 'model', 'k1', 'max', 10.0);
    voxulus ('set', 'model', 'k2', 'value', 1.0);
    voxulus ('set', 'model', 'k2', 'initvalue', 1.0);
    voxulus ('set', 'model', 'k2', 'min', 0.001);
    voxulus ('set', 'model', 'k2', 'max', 10.0);
    %voxulus ('set', 'model', 'V', 'value', 0.75);
    %voxulus ('set', 'model', 'V', 'state', 'const');
    voxulus ('set', 'model', 'k3', 'value', 0.0);
    voxulus ('set', 'model', 'k3', 'initvalue', 0.0);
    voxulus ('set', 'model', 'k3', 'min', 0.001);
    voxulus ('set', 'model', 'k3', 'max', 1.0);
    voxulus ('set', 'model', 'k3', 'state', 'const');
    %voxulus ('set', 'model', 'k4', 'value', 0.0);
    %voxulus ('set', 'model', 'k4', 'initvalue', 0.0);
    %voxulus ('set', 'model', 'k4', 'min', 0.0);
    %voxulus ('set', 'model', 'k4', 'max', 1.0);
    %voxulus ('set', 'model', 'k4', 'state', 'const');
    voxulus ('set', 'model', 'tbv', 'value', 0.5);
    voxulus ('set', 'model', 'tbv', 'initvalue', 0.5);
    voxulus ('set', 'model', 'tbv', 'min', 0.001);
    voxulus ('set', 'model', 'tbv', 'max', 1.0);
    %voxulus ('set', 'model', 'alpha', 'value', 1.0);
    %voxulus ('set', 'model', 'alpha', 'initvalue', 1.0);
    %voxulus ('set', 'model', 'alpha', 'min', 0.001);
    %voxulus ('set', 'model', 'alpha', 'max', 1.0);
    %voxulus ('set', 'model', 'alpha', 'state', 'const');
    %voxulus ('set', 'model', 'beta', 'value', 0.5);
    %voxulus ('set', 'model', 'beta', 'initvalue', 0.5);
    %voxulus ('set', 'model', 'beta', 'min', 0.001);
    %voxulus ('set', 'model', 'beta', 'max', 1.0);
    %voxulus ('set', 'model', 'beta', 'state', 'const');
    voxulus ('set', 'model', 'unit', u);
  else
    error ('unknown tracer ''%s''', tracer);
  end

  if ~isempty (voi1)
    inputList = create_inputs ('GenericInputFunction', u);
    if isscalar (e)
      [ bic, tac ] = input_estimate (d, e, t, u, voi1, inputList, 1);
    else
      [ bic, tac ] = input_estimate (d, [], t, u, voi1, inputList, 1);
    end
    save bp tac
    [ bicMin, i ] = min (bic(:,3));
    inputFct = inputList{i};
    bp_tac = tac(:,i+3);
    voxulus ('select', 'input', inputFct);
    voxulus ('select', 'model', model);
    voxulus ('set', 'model', 'input', inputFct);
    for j=1:size(inputList,1)
      if j ~= i
        voxulus ('delete', 'input', inputList{j});
      end
    end
   end

  voxulus ('set', 'data', 'DATA', d, t);
  voxulus ('set', 'data', 'DATA', 'midframe_approx', 1);
  voxulus ('set', 'data', 'DATA', 'unit', u);
  if ~isempty (e) && ~isscalar (e)
    voxulus ('set', 'data', 'ERROR', e, t);
    voxulus ('set', 'data', 'ERROR', 'midframe_approx', 1);
    voxulus ('set', 'data', 'ERROR', 'unit', u);
  end
  voxulus ('set', 'voi', 'MYO', voi2);
  %voxulus ('show', 'voi');


voxulus ('select', 'optimizer', 'Levenberg-Marquardt');
voxulus ('create', 'protocol', 'MyProtocol', model);
voxulus ('select', 'protocol', 'MyProtocol');
voxulus ('set', 'protocol', 'optimizer', 'Levenberg-Marquardt');

%
% average the data on the VOIs
%
[ myo_data, t2, u ] = voxulus ('filter', 'DATA', 'voi', 'MYO', 'mean');
voxulus ('set', 'data', 'myo_data', myo_data, t2);
voxulus ('set', 'data', 'myo_data', 'unit', u{1});
if ~isempty (e) && ~isscalar (e)
  myo_error = voxulus ('filter', 'ERROR', 'voi', 'MYO', 'mean');
  voxulus ('set', 'data', 'myo_error', myo_error, t2);
  voxulus ('set', 'data', 'myo_error', 'unit', u{1});
end

voxulus ('set', 'voi', 'dummy', logical (ones (1, 1, 1)));
if isempty (e) || isscalar (e)
  voxulus ('process', 'myo_data', 'none', 'none', 'dummy');
  myo_tac = voxulus ('get', 'tac', 'myo_data', 'none', 'dummy');
else
  voxulus ('process', 'myo_data', 'myo_error', 'none', 'dummy');
  myo_tac = voxulus ('get', 'tac', 'myo_data', 'myo_error', 'dummy');
end
[ k1, k1err, k1init, k1state, k1unit ] = voxulus ('get', 'model', 'k1');
%
% bloodflow
%
y = [ k1/1.042 k1err/1.042 ];
%
% cleanup
%
voxulus ('delete', 'protocol', 'MyProtocol', model);
voxulus ('delete', 'voi', 'dummy');
voxulus ('delete', 'voi', 'MYO');
voxulus ('delete', 'data', 'myo_data');
voxulus ('delete', 'data', 'DATA');
if ~isempty (e) && ~isscalar (e)
  voxulus ('delete', 'data', 'myo_error');
  voxulus ('delete', 'data', 'ERROR');
end
