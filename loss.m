function L = loss(x,y,s,param)
    K = size(x,2);
    beta1 = param(1:K);
    beta2 = param(K+1:2*K);
    lam = param(2*K+1);
    c = param(2*K+2);
    
    g = 1-exp(-lam*(s-c).^2);
    e = y - (x*beta1 + g.*x*beta2);
    L = e'*e;
end