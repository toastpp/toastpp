function [err,varargout] = privObjective (proj, data, sd, hReg, x)

err_data = full(sum(((data-proj)./sd).^2));
err_prior = 0;

if nargin >= 5
     % Add regularisation component
     err_prior = full(toastRegulValue (hReg, x));
end

err = err_data + err_prior;
if nargout >= 1
    varargout{1} = err_data;
end
if nargout >= 2
    varargout{2} = err_prior;
end
