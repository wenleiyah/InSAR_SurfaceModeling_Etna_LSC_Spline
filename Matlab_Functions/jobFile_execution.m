function jobFile_execution(filename)

% Function to automatically execute the job file of geoSplinter on both 
% Windows and macOS.
% 
% INPUT VARIABLES:
% Job file name: filename
%
%
% Roberto Monti
% Politecnico di Milano
% Last update: November 2024


% Job file execution
if ismac
    job_execution = strcat('./geoSplinter_analysis_macOS < ./job/', filename, '.job');
elseif ispc
    job_execution = strcat('.\geoSplinter_analysis.exe < .\job\', filename, '.job');
end
    
system(job_execution)
system('exit')