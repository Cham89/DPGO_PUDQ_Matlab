function Sigma = generate_pudq_cov(sigma_pudq, df)
% Edge measurement noise covariance
    Sigma_w = ones(3,3) + diag(unifrnd(realmin,1,[3,1]));
    Sigma = wishrnd(sigma_pudq*Sigma_w,df)/df;
end
