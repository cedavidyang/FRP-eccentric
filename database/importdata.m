clear;
% import excel databaes to matlab data for model error analysis
[cfrpdatabase_raw, ~, ~] = xlsread('FRP confined RC columns under eccentric loading.xlsx', 'CFRP');
[gfrpdatabase_raw, ~, ~] = xlsread('FRP confined RC columns under eccentric loading.xlsx', 'GFRP');
indx_cfrp = (cfrpdatabase_raw(:, 26)~=0);
% indx_cfrp = (cfrpdatabase_raw(:, 26)==1);
cfrpdatabase = cfrpdatabase_raw(indx_cfrp, :);    % all specimens
% cfrpdatabase = cfrpdatabase_raw([2,4:6], :);   % Bisby-Ranger-2010
% cfrpdatabase = cfrpdatabase_raw([19, 20, 22], :);     %Fitzwilliam-Bisby-2010
% cfrpdatabase = cfrpdatabase_raw(29:36, :);    % Jiang-2014: 1layer
% cfrpdatabase = cfrpdatabase_raw(37:end, :);    % Jiang-2014: 2layers
indx_gfrp = (gfrpdatabase_raw(:, 23)==1);
gfrpdatabase = gfrpdatabase_raw(indx_gfrp, :);
save('frpdatabase.mat', 'cfrpdatabase', 'gfrpdatabase');