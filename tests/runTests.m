% Copyright 2021 The MathWorks, Inc.

import matlab.unittest.TestRunner
import matlab.unittest.Verbosity
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat
 
% Add paths
current_path = mfilename("fullpath");
path_root = fileparts(fileparts(current_path));
path_src = fullfile(path_root, "src");
addpath(path_src);

% Create a test suite 
suite = testsuite(path_src, 'IncludeSubfolders', true);

% Create a test runner that displays test run progress at the matlab.unittest.Verbosity.Detailed level
runner = TestRunner.withTextOutput('OutputDetail',Verbosity.Detailed); 

% Create a CodeCoveragePlugin instance and add it to the test runner
sourceFolder = fullfile(path_root, "src");
reportFile = 'coverage.xml';
reportFormat = CoberturaFormat(reportFile);
% Use the CodeCoveragePlugin to collect coverage information for all the .m files in the src folder
p = CodeCoveragePlugin.forFolder(sourceFolder, ...
    'IncludingSubfolders', true, ...
    'Producing', reportFormat, ...
    'Files', fullfile(sourceFolder, '**', '*.m'), ... % Only contain .m files in the src folder
    'ExcludeFiles', {...
        fullfile(sourceFolder, '**', 'test_*.m'), ... % Exclude test files
        fullfile(sourceFolder, '**', '*_test.m'), ...
        fullfile(sourceFolder, '**', 'tests', '*.m') ...
    });
runner.addPlugin(p)
 
% Run the tests and fail the build if any of the tests fails
results = runner.run(suite);  
nfailed = nnz([results.Failed]);
assert(nfailed == 0,[num2str(nfailed) ' test(s) failed.'])

% Remove paths
rmpath(path_src);