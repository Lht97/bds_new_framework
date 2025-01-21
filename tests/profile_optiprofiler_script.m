clear all
% options.dim = 'big';
options.solver_names = {'cbds', 'newuoa'};
% feature = {'perturbed_x0_0.001', 'perturbed_x0_1', 'perturbed_x0_10', 'perturbed_x0_100', 'truncated_1', 'truncated_2', 'truncated_3', 'truncated_4', 'permuted', 'random_nan_5', 'random_nan_10', 'random_nan_20'};
% for i = 1:length(feature)
%     options.feature_name = feature{i};
%     profile_optiprofiler(options);
% end
options.dim = 'big';
%options.solver_names = {'rbds-zero-delay', 'rbds-one-delay', 'rbds-eighth-delay', 'rbds-quarter-delay', 'rbds-half-delay', 'rbds-n-minus-one-delay'};
options.feature_name = 'quantized_1';
profile_optiprofiler(options);

% options.feature_name = 'linearly_transformed';
% profile_optiprofiler(options);
% options.feature_name = 'perturbed_x0_1';
% profile_optiprofiler(options);

% options.feature_name = 'noisy_1e-1';
% profile_optiprofiler(options);

% options.feature_name = 'noisy_1e-2';
% profile_optiprofiler(options);

% options.feature_name = 'noisy_1e-3';
% profile_optiprofiler(options);

% options.feature_name = 'noisy_1e-4';
% profile_optiprofiler(options);

% options.feature_name = 'rotation_noisy_1e-1';
% profile_optiprofiler(options);

% options.feature_name = 'rotation_noisy_1e-2';
% profile_optiprofiler(options);

% options.feature_name = 'rotation_noisy_1e-3';
% profile_optiprofiler(options);

% options.feature_name = 'rotation_noisy_1e-4';
% profile_optiprofiler(options);



