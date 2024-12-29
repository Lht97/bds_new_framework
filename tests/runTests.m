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

% Create a test runner
runner = TestRunner.withTextOutput('OutputDetail',Verbosity.Detailed); 

% 创建覆盖率插件
sourceFolder = fullfile(path_root, "src");
reportFile = 'coverage.xml';
reportFormat = CoberturaFormat(reportFile);

% 使用 IncludeSubfolders 来包含所有子文件夹
p = CodeCoveragePlugin.forFolder(sourceFolder, ...
    'IncludingSubfolders', true, ...
    'Producing', reportFormat);

runner.addPlugin(p)
 
% Run tests
results = runner.run(suite);  
nfailed = nnz([results.Failed]);
assert(nfailed == 0,[num2str(nfailed) ' test(s) failed.'])

% Remove paths
rmpath(path_src);