vt_ne_exp = readmatrix('VT_Ne.csv');
vv_ne_exp = readmatrix('VV_Ne.csv');
vv_ne_exp(:, 2) = 1 ./ vv_ne_exp(:, 2);

%%
fit_f = fittype('exp(a*x^2+b*x+c)')
f_vt = fit(vt_ne_exp(:, 1), vt_ne_exp(:, 2), fit_f, 'TolFun', 1e-30, ...
    'TolX', 1e-30, 'MaxIter', 1000, 'MaxFunEvals', 1000)
f_vv = fit(vv_ne_exp(:, 1), vv_ne_exp(:, 2), fit_f, 'TolFun', 1e-30, ...
    'TolX', 1e-30, 'MaxIter', 1000, 'MaxFunEvals', 1000)

%%
t = (0.05:0.001:0.15)';

figure
semilogy(vt_ne_exp(:, 1), vt_ne_exp(:, 2), 'b*', ...
    vv_ne_exp(:, 1), vv_ne_exp(:, 2), 'r*', t, f_vt(t), '-b', t, f_vv(t), '-r')
grid minor
axis padded
