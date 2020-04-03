function start_up(input_workdir)

% Add auxilliary programs to path.
%
% Input argument(input_workdir): Specify the working directory, the
% subdirectory 'Programs' of which, contains the necessary auxilliary 
% programs.
% Subfunctions: -

% ----- Add programs to path ----- %
addpath([input_workdir '\Programs\spm12']);
addpath([input_workdir '\Programs\r2agui_v27']);
addpath([input_workdir '\Programs\hiro3']);
addpath(genpath([input_workdir '\Programs\DPABI_V4.2\DPABI_V4.2_190919']));
