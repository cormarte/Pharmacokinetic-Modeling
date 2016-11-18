units = 'Bq/ml';
slope = 0.46163022518158;

% Initialisation
verbosity = 3;
voxulus ('reset');
voxulus ('set', 'verbosity', verbosity);
voxulus ('set', 'seed', 1234);

% Data
[data, time] = voxulus('read','data', filename);
data = data*slope;
voxulus ('set', 'data', 'data', data, time, units);

% Input function
voxulus ('select', 'input', 'GenericInputFunction');
voxulus ('set', 'input', 'unit', units);
voxulus ('set', 'input', 'tau', 'value', 1);
voxulus ('set', 'input', 'tau', 'unit', 's');
voxulus ('set', 'input', 'tau', 'state', 'const');
voxulus ('set', 'input', 't0', 'value', 0);
voxulus ('set', 'input', 't0', 'unit', 's');
voxulus ('set', 'input', 't0', 'state', 'const');

voxulus ('set', 'input', 'A1', 'unit', units);
%voxulus ('set', 'input', 'A1', 'value', 665.0592);
%voxulus ('set', 'input', 'A1', 'initvalue', 665.0592);
%voxulus ('set', 'input', 'B1', 'value', 2);
%voxulus ('set', 'input', 'B1', 'state', 'const');
%voxulus ('set', 'input', 'C1', 'value', 0.2781);
%voxulus ('set', 'input', 'C1', 'initvalue', 0.2781);
voxulus ('set', 'input', 'A2', 'unit', units);
%voxulus ('set', 'input', 'A2', 'value', 59.9671);
%voxulus ('set', 'input', 'A2', 'initvalue', 59.9671);
%voxulus ('set', 'input', 'B2', 'value', 1);
%voxulus ('set', 'input', 'B2', 'state', 'const');
%voxulus ('set', 'input', 'C2', 'value', 0.0276);
%voxulus ('set', 'input', 'C2', 'initvalue', 0.0276);

% Model
voxulus ('select', 'model', 'Generic_2TC');
voxulus ('set', 'model', 'unit', units);
voxulus ('set', 'model', 'input', 'GenericInputFunction');

voxulus ('set', 'model', 'k1', 'value', 1.0);
voxulus ('set', 'model', 'k1', 'initvalue', 1.0);
voxulus ('set', 'model', 'k2', 'value', 1.0);
voxulus ('set', 'model', 'k2', 'initvalue', 1.0);
voxulus ('set', 'model', 'k3', 'value', 0.9);
voxulus ('set', 'model', 'k3', 'initvalue', 0.9);
voxulus ('set', 'model', 'k4', 'value', 0.08);
voxulus ('set', 'model', 'k4', 'initvalue', 0.08);
voxulus ('set', 'model', 'alpha', 'value', 1.0);
voxulus ('set', 'model', 'alpha', 'initvalue', 1.0);
voxulus ('set', 'model', 'alpha', 'state', 'const');
voxulus ('set', 'model', 'beta', 'value', 0.0);
voxulus ('set', 'model', 'beta', 'initvalue', 0.0);
voxulus ('set', 'model', 'beta', 'state', 'const');

voxulus ('select', 'optimizer', 'Levenberg-Marquardt');
voxulus ('create', 'protocol', 'Generic_2TC', 'Generic_2TC');

%Parameters
[ dm, t ] = voxulus ('filter', 'test_data', 'voi', 'bif_voi', 'mean');
voxulus ('set', 'data', 'input_data', dm, t);
voxulus ('set', 'data', 'input_data', 'unit', 'nCi/ml');
[ em, t ] = voxulus ('filter', 'test_error', 'voi', 'bif_voi', 'mean');
voxulus ('set', 'data', 'input_error', em, t);
voxulus ('set', 'data', 'input_error', 'unit', 'nCi/ml');

voxulus ('set', 'voi', 'input', logical(1));
voxulus ('show', 'voi');

voxulus ('process', 'input_data', 'input_error', 'input', 'none');

[ A0, A0err, A0init, A0state, A0unit ] = voxulus ('get', 'input', 'A0');
[ A1, A1err, A1init, A1state, A1unit ] = voxulus ('get', 'input', 'A1');
[ A2, A2err, A2init, A2state, A2unit ] = voxulus ('get', 'input', 'A2');
[ c0, c0err, c0init, c0state, c0unit ] = voxulus ('get', 'input', 'c0');
[ c11, c11err, c11init, c11state, c11unit ] = voxulus ('get', 'input', 'c11');
[ c21, c21err, c21init, c21state, c21unit ] = voxulus ('get', 'input', 'c21');
[ c12, c12err, c12init, c12state, c12unit ] = voxulus ('get', 'input', 'c12');
[ c22, c22err, c22init, c22state, c22unit ] = voxulus ('get', 'input', 'c22');

voxulus ('process', 'test_data', 'test_error', 'none', 'test_output_mask');
voxulus ('show', 'data');

[ k1, k1err, k1init, k1state, k1unit ] = voxulus ('get', 'model', 'k1');
[ k2, k2err, k2init, k2state, k2unit ] = voxulus ('get', 'model', 'k2');
[ k3, k3err, k3init, k3state, k3unit ] = voxulus ('get', 'model', 'k3');
[ k4, k4err, k4init, k4state, k4unit ] = voxulus ('get', 'model', 'k4');
[ beta, betaerr, betainit, betastate, betaunit ] = voxulus ('get', 'model', 'beta');

