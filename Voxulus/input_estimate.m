%INPUT_ESTIMATE estimate a collection of input functions
%   [ Y, TAC ] = INPUT_ESTIMATE (D, E, T, U, V, LST[, VERBOSITY])
%      estimates all Voxulus input functions in the list 'LST' (cell array of
%      strings). The input variables are:
%       D: 4D data set
%       E: 4D error data set (maybe empty)
%       T: timing information for data and error
%       U: data and error unit (see Voxulus list)
%       V: volume of interest on D, E
%       LST: list of functions to estimate
%      The output variables are:
%       Y: 2D array of model complexity criteria for all estimated functions
%       TAC: 2D array of corresponding time activity curves
%
%      The indexing order of D and E is Z,Y,X,T and of V is Z,Y,X.
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%

function [ y, tac ] = input_estimate (d, e, t, u, v, lst, verbosity)
  if nargin < 7
    verbosity = 0;
  end

  voxulus ('set', 'voi', 'BP', v);
  voxulus ('set', 'data', 'DATA', d, t);
  voxulus ('set', 'data', 'DATA', 'midframe_approx', 1);
  voxulus ('set', 'data', 'DATA', 'unit', u);
  [ bp_data, t2, u ] = voxulus ('filter', 'DATA', 'voi', 'BP', 'mean');
  voxulus ('set', 'data', 'bp_data', bp_data, t2);
  voxulus ('set', 'data', 'bp_data', 'unit', u{1});
  bp_error = ones (size (bp_data));
  if isscalar (e)
    bp_error = ones (size (bp_data));
    [ bpmax, tmax ] = max (bp_data);
    bp_error(:,:,:,tmax) = 1.0 / e;
    %bp_error(:,:,:,tmax-1) = (1.0 + 1.0 / e) / 2;
    %bp_error(:,:,:,tmax+1) = bp_error(:,:,:,tmax-1);
  elseif ~isempty (e)
    voxulus ('set', 'data', 'ERROR', e, t);
    voxulus ('set', 'data', 'ERROR', 'midframe_approx', 1);
    voxulus ('set', 'data', 'ERROR', 'unit', u{1});
    [ bp_error, t2, u ] = voxulus ('filter', 'ERROR', 'voi', 'BP', 'mean');
  end
  voxulus ('set', 'data', 'bp_error', bp_error, t2);
  voxulus ('set', 'data', 'bp_error', 'unit', u{1});    

  voxulus ('select', 'model', 'Generic_1TC');
  voxulus ('set', 'model', 'input', lst{1});
  voxulus ('set', 'model', 'unit', u{1});
  voxulus ('select', 'optimizer', 'Levenberg-Marquardt');
  voxulus ('create', 'protocol', 'InputEstimate', 'Generic_1TC');
  voxulus ('select', 'protocol', 'InputEstimate');
  voxulus ('set', 'protocol', 'optimizer', 'Levenberg-Marquardt');
  voxulus ('set', 'voi', 'dummy', logical (ones (1, 1, 1)));

  y = zeros (size (lst, 1), 5);
  tac = zeros (size (t, 1), size (lst, 1) + 3);
  tac(:,1) = t(:,1);
  tac(:,2) = bp_data;
  if ~isempty (e)
    tac(:,3) = bp_error;
  end
  for i = 1:size (lst,1)
    voxulus ('select', 'input', lst{i});
    voxulus ('set', 'input', 'init', 'DATA', 'BP');
    voxulus ('select', 'model', 'Generic_1TC');
    voxulus ('set', 'model', 'input', lst{i});
    voxulus ('set', 'protocol', 'reset');
    voxulus ('process', 'bp_data', 'bp_error', 'dummy', 'none');
    tmp = voxulus ('get', 'tac', 'bp_data', 'bp_error', 'dummy');
    tac(:,i+3) = tmp(:,5);
    y(i,1) = log (voxulus ('get', 'chisquare'));
    y(i,2) = voxulus ('get', 'aicc');
    y(i,3) = voxulus ('get', 'bic');
    y(i,4) = voxulus ('get', 'fpe');
    y(i,5) = voxulus ('get', 'hannan-quinn');
  end
  voxulus ('delete', 'protocol', 'InputEstimate');
  voxulus ('delete', 'voi', 'BP');
  voxulus ('delete', 'voi', 'dummy');
  voxulus ('delete', 'data', 'bp_data');
  voxulus ('delete', 'data', 'DATA');
  voxulus ('delete', 'data', 'bp_error');
  if ~isempty (e) && ~isscalar (e)
    voxulus ('delete', 'data', 'ERROR');
  end

  if verbosity
    fprintf (' N INPUTFCT   log CHI^2       AICC        BIC        FPE        H-Q\n');
    fprintf ('===================================================================\n');
    for i=1:size(lst,1)
      fprintf ('%2d %-8s ', i, lst{i});
      fprintf (' %10g', y(i,:));
      fprintf ('\n');
    end
    fprintf ('===================================================================\n');
    [ bicMin, j ] = min (y(:,3));
    fprintf ('%2d %-8s ', j, lst{j});
    fprintf (' %10g', y(j,:));
    fprintf ('\n');
  end
  if verbosity > 1
    figure (4);
    errorbar (tac(:,1), tac(:,2), tac(:,3), 'b-x', 'LineWidth', 2);
    hold on;
    plot (tac(:,1), tac(:,4:size(tac,2)), 'LineWidth', 2);
    set (gca, 'FontSize', 14);
    set (gca, 'FontWeight', 'bold');
    title ('Input Estimate');
    legend ('data');
    xlabel ('t [s]');
    ylabel (sprintf ('c [%s]', u{1}));
    pause;
  end
