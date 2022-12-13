function F = DiscreteCDF(x, distribution_type, parameters)
f = @(x)DiscretePDF(x,distribution_type, parameters);
x = sort(x);
n_length = length(x);
x_min=0;
switch (distribution_type) 
    case 'Bernoulli'
        x_min = 0;
    case 'bino'
        x_min = 0;
    case 'Poisson'
        x_min = 0;
    case 'geometric'
        x_min = 1;
    case 'nbin'
        x_min = 1;
    case 'hyge'
        N = parameters(1);
        M = parameters(2);
        K = parameters(3);
        x_min = max(0,K-N+M);
end


F = zeros(1,n_length);
F(1) = sum(f(x_min:x(1)));

for i=2:n_length
    v = (x(i-1)+1):x(i);
   F(i) = F(i-1) + sum(f(v));
end   
    
end

