function stats=stats(var);

% 
%  function stats=stats(var);
%
% Returns structure contining statistics of vector
%
% stats.n = length(var);
% stats.max = max(var);
% stats.min = min(var);
% stats.mean = mean(var);
% stats.std = std(var);
% stats.median = median(var);
%

var=full(var(:));
  

stats.n = length(var);
stats.max = max(var);
stats.min = min(var);
stats.mean = mean(var);
stats.std = std(var);
stats.median = median(var);
  
disp(['Mean +/- SD = ' num2str(stats.mean) ' +/- ' num2str(stats.std) '.']);
disp(['Max = ' num2str(stats.max) ', min = ' num2str(stats.min) ', median = ' num2str(stats.median) ', n = ' num2str(stats.n)]);
disp(' ');
