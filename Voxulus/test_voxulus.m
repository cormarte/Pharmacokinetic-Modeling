global VOXULUS_VERSION;

VOXULUS_VERSION = 'opt';

voxulus ('reset');
voxulus ('set', 'verbosity', 3);

[ d, t ] = voxulus ('read', 'data', 'testcases/TestVoxulus/input/test_data');
voxulus ('set', 'data', 'test_data', d, t, 'kBq/ml');
[ e, t2 ] = voxulus ('read', 'data', 'testcases/TestVoxulus/input/test_error');
voxulus ('set', 'data', 'test_error', e, t2, 'kBq/ml');
%voxulus ('set', 'data', 'test_data', 'midframe_approx', 1);
%voxulus ('set', 'data', 'test_error', 'midframe_approx', 1);
voxulus ('show', 'data');

input_mask = voxulus ('read', 'voi', 'testcases/TestVoxulus/input/test_input_mask');
voxulus ('set', 'voi', 'test_input_mask', input_mask);
output_mask = voxulus ('read', 'voi', 'testcases/TestVoxulus/input/test_output_mask');
voxulus ('set', 'voi', 'test_output_mask', output_mask);
voxulus ('show', 'voi');

voxulus ('select', 'input', 'FMISO_InputFunction');
voxulus ('show', 'input');
voxulus ('set', 'input', 'tau', 'value', 1.0);
voxulus ('set', 'input', 'tau', 'unit', 's');
voxulus ('set', 'input', 'tau', 'state', 'const');
voxulus ('set', 'input', 'A0', 'value', 10.0);
voxulus ('set', 'input', 'A0', 'initvalue', 10.0);
voxulus ('set', 'input', 'A0', 'unit', 'kBq/ml');
voxulus ('set', 'input', 'c0', 'value', 0.33);
voxulus ('set', 'input', 'c0', 'initvalue', 0.33);
voxulus ('set', 'input', 'A1', 'value', 20.0);
voxulus ('set', 'input', 'A1', 'initvalue', 20.0);
voxulus ('set', 'input', 'A1', 'unit', 'kBq/ml');
voxulus ('set', 'input', 'c11', 'value', 1.0);
voxulus ('set', 'input', 'c11', 'initvalue', 1.0);
voxulus ('set', 'input', 'c21', 'value', 2.0);
voxulus ('set', 'input', 'c21', 'initvalue', 2.0);
voxulus ('set', 'input', 'A2', 'value', 500.0);
voxulus ('set', 'input', 'A2', 'initvalue', 500.0);
voxulus ('set', 'input', 'A2', 'unit', 'kBq/ml');
voxulus ('set', 'input', 'c12', 'value', 12.0);
voxulus ('set', 'input', 'c12', 'initvalue', 12.0);
voxulus ('set', 'input', 'c22', 'value', 120.0);
voxulus ('set', 'input', 'c22', 'initvalue', 120.0);
voxulus ('set', 'input', 'unit', 'kBq/ml');
voxulus ('set', 'input', 'init', 'test_data', 'test_input_mask');

voxulus ('select', 'model', 'FMISO_3TC');
voxulus ('set', 'model', 'k1', 'value', 0.2);
voxulus ('set', 'model', 'k1', 'initvalue', 0.2);
voxulus ('set', 'model', 'k1', 'unit', '1/min');
voxulus ('set', 'model', 'k2', 'value', 0.2);
voxulus ('set', 'model', 'k2', 'initvalue', 0.2);
voxulus ('set', 'model', 'k2', 'unit', '1/min');
voxulus ('set', 'model', 'k3', 'value', 0.02);
voxulus ('set', 'model', 'k3', 'initvalue', 0.02);
voxulus ('set', 'model', 'k3', 'unit', '1/min');
voxulus ('set', 'model', 'k4', 'value', 0.04);
voxulus ('set', 'model', 'k4', 'initvalue', 0.04);
voxulus ('set', 'model', 'k4', 'unit', '1/min');
voxulus ('set', 'model', 'alpha', 'value', 0.36);
voxulus ('set', 'model', 'alpha', 'state', 'const');
voxulus ('set', 'model', 'beta', 'value', 0.05);
voxulus ('set', 'model', 'beta', 'initvalue', 0.05);
voxulus ('set', 'model', 'eta', 'value', 0.5);
voxulus ('set', 'model', 'eta', 'state', 'const');
voxulus ('set', 'model', 'V', 'value', 1.0);
voxulus ('set', 'model', 'V', 'state', 'const');
voxulus ('set', 'model', 'input', 'FMISO_InputFunction');

voxulus ('select', 'optimizer', 'Levenberg-Marquardt');
voxulus ('create', 'protocol', 'FMISO_3TC', 'FMISO_3TC');

[ dm, t ] = voxulus ('filter', 'test_data', 'voi', 'test_input_mask', 'mean');
voxulus ('set', 'data', 'input_data', dm, t);
voxulus ('set', 'data', 'input_data', 'unit', 'kBq/ml');
[ em, t ] = voxulus ('filter', 'test_error', 'voi', 'test_input_mask', 'mean');
voxulus ('set', 'data', 'input_error', em, t);
voxulus ('set', 'data', 'input_error', 'unit', 'kBq/ml');

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

%
% figure 1: TAC on input mask
%
data1 = get_tac (d, input_mask);
error1 = get_tac (e, input_mask);
tac1 = voxulus ('get', 'tac', 'test_data', 'test_error', 'test_input_mask', 'tac1');
figure (1);
errorbar (t(:,1), tac1(:,3), tac1(:,4), 'bx', 'LineWidth', 2);
xlim = get (gca, 'XLim');
xlim(1) = 0;
set (gca, 'XLim', xlim);
ylim = get (gca, 'YLim');
ylim(1) = 0;
set (gca, 'YLim', ylim);
hold on;
plot (tac1(:,1), tac1(:,5), 'r-', 'LineWidth', 2);
hold off;
set (gca, 'FontSize', 14);
set (gca, 'FontWeight', 'bold');
title ('TAC on input mask');
legend ('data', 'fit');
xlabel ('t [s]');
ylabel ('C [hBq/ml]');
set_text (0.6, 0.6, sprintf ('A0 = %.4g +/- %.4g %s', A0, A0err, A0unit{1}));
set_text (0.6, 0.56, sprintf ('c0 = %.4g +/- %.4g %s', c0, c0err, c0unit{1}));
set_text (0.6, 0.52, sprintf ('A1 = %.4g +/- %.4g %s', A1, A1err, A1unit{1}));
set_text (0.6, 0.48, sprintf ('c11 = %.4g +/- %.4g %s', c11, c11err, c11unit{1}));
set_text (0.6, 0.44, sprintf ('c21 = %.4g +/- %.4g %s', c21, c21err, c21unit{1}));
set_text (0.6, 0.40, sprintf ('A2 = %.4g +/- %.4g %s', A2, A2err, A2unit{1}));
set_text (0.6, 0.36, sprintf ('c12 = %.4g +/- %.4g %s', c12, c12err, c12unit{1}));
set_text (0.6, 0.32, sprintf ('c22 = %.4g +/- %.4g %s', c22, c22err, c22unit{1}));

%
% figure 2: TAC on output mask
%
tac2 = voxulus ('get', 'tac', 'test_data', 'test_error', 'test_output_mask', 'tac2');
figure (2);
errorbar (t(:,1), tac2(:,3), tac2(:,4), 'bx', 'LineWidth', 2);
xlim = get (gca, 'XLim');
xlim(1) = 0;
set (gca, 'XLim', xlim);
set (gca, 'FontSize', 14);
set (gca, 'FontWeight', 'bold');
hold on;
plot (tac2(:,1), tac2(:,2), 'r-');
hold off;
title ('TAC on output mask');
legend ('data', 'fit');
xlabel ('t [s]');
ylabel ('C [hBq/ml]');
k = get_tac (k1, output_mask);
kerr = get_tac (k1err, output_mask);
set_text (0.6, 0.6, sprintf ('k1 = %.4g +/- %.4g %s', k, kerr, k1unit{1}));
k = get_tac (k2, output_mask);
kerr = get_tac (k2err, output_mask);
set_text (0.6, 0.56, sprintf ('k2 = %.4g +/- %.4g %s', k, kerr, k2unit{1}));
k = get_tac (k3, output_mask);
kerr = get_tac (k3err, output_mask);
set_text (0.6, 0.52, sprintf ('k3 = %.4g +/- %.4g %s', k, kerr, k3unit{1}));
k = get_tac (k4, output_mask);
kerr = get_tac (k4err, output_mask);
set_text (0.6, 0.48, sprintf ('k4 = %.4g +/- %.4g %s', k, kerr, k4unit{1}));
k = get_tac (beta, output_mask);
kerr = get_tac (betaerr, output_mask);
set_text (0.6, 0.44, sprintf ('beta = %.4g +/- %.4g %s', k, kerr, betaunit{1}));

%
% figure 3: TAC at specific voxel
%
tac3 = voxulus ('get', 'tac', 'test_data', 'test_error', [ 7, 2, 0 ], 'tac3');
figure (3);
errorbar (t(:,1), tac3(:,3), tac3(:,4), 'bx', 'LineWidth', 2);
xlim = get (gca, 'XLim');
xlim(1) = 0;
set (gca, 'XLim', xlim);
ylim = get (gca, 'YLim');
ylim(1) = 0;
set (gca, 'YLim', ylim);
set (gca, 'FontSize', 14);
set (gca, 'FontWeight', 'bold');
hold on;
plot (tac3(:,1), tac3(:,2), 'r-', 'LineWidth', 2);
hold off;
title ('TAC at voxel 7,2,0');
legend ('data', 'fit');
xlabel ('t [s]');
ylabel ('C [hBq/ml]');
set_text (0.6, 0.6, sprintf ('k1 = %.4g +/- %.4g %s', k1(1,3,8), k1err(1,3,8), k1unit{1}));
set_text (0.6, 0.56, sprintf ('k2 = %.4g +/- %.4g %s', k2(1,3,8), k2err(1,3,8), k2unit{1}));
set_text (0.6, 0.52, sprintf ('k3 = %.4g +/- %.4g %s', k3(1,3,8), k3err(1,3,8), k3unit{1}));
set_text (0.6, 0.48, sprintf ('k4 = %.4g +/- %.4g %s', k4(1,3,8), k4err(1,3,8), k4unit{1}));
set_text (0.6, 0.44, sprintf ('beta = %.4g +/- %.4g %s', beta(1,3,8), betaerr(1,3,8), betaunit{1}));

%
% figure 4: parametric maps
%
figure (4);
subplot (2,2,1);
imagesc (squeeze(k1(1,:,:)), [ 0 1 ]);
title ('k1');
subplot (2,2,2);
imagesc (squeeze(k2(1,:,:)), [ 0 1 ]);
colorbar ('EastOutside');
title ('k2');
subplot (2,2,3);
imagesc (squeeze(k3(1,:,:)), [ 0 0.1 ]);
title ('k3');
subplot (2,2,4);
imagesc (squeeze(k4(1,:,:)), [ 0 0.1 ]);
colorbar ('EastOutside');
title ('k4');

%
% figure 5: parametric error maps
%
figure (5);
subplot (2,2,1);
imagesc (squeeze(k1err(1,:,:)), [ 0 1 ]);
title ('k1 error');
subplot (2,2,2);
imagesc (squeeze(k2err(1,:,:)), [ 0 1 ]);
colorbar ('EastOutside');
title ('k2 error');
subplot (2,2,3);
imagesc (squeeze(k3err(1,:,:)), [ 0 0.1 ]);
title ('k3 error');
subplot (2,2,4);
imagesc (squeeze(k4err(1,:,:)), [ 0 0.1 ]);
colorbar ('EastOutside');
title ('k4 error');
