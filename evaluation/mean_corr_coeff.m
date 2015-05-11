function result = mean_corr_coeff(corr_coefs)

% first transform the corr_coefs using z transform
    coeff_z = fisherz(corr_coefs);
% compute the mean value
    mean_z = mean(coeff_z);
% transform with back
    result = inverse_fisherz(mean_z);
end



function z=fisherz(r)
    z=0.5*log((1+r)./(1-r));
end


function r=inverse_fisherz(z)
    r = (exp(2*z) -1) ./ (exp(2*z) +1);
end