function nrmse = NRMSE_calc(x,x_est,N)
    rmse = sqrt(sum((x_est-x').^2)/N);
    nrmse = rmse/(max(x)-min(x));
end