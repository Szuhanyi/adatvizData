function f=ContinuousPdf(x, distribution_type, parameters)

x = sort(x);
n_length = length(x);

switch(distribution_type)

  case 'normal'
    mu = parameters(1);
    sigma = parameters(2);
    if (sigma <= 0)
      error('The parameter 2 must be positive');
    end
    f = (1.0 / sqrt(2.0* pi) / sigma) * exp(-(x - mu).^2 / 2.0 / sigma^2);

  case 'beta'
    % B(a,b) - beta distribution
    % uses the Euler's beta function

    a = parameters(1);
    b = parameters(2);

    if( a <= 0 || b <= 0)
      error('parameters must be positive!')
    end
    f = zeros(1,n_length);

    for i=1:n_length
      if (x(i) > 0)
        f(i) = 1 / beta(a,b) * x(i)^(a-1) * (1-x(i))^(b-1);
      else
        f(i) = 0;
      end

    end




  case 'chi2'

  case 'exp'
    % Exp(lambda)
    % in the other books
    lambda = parameters(1);
    mu = 1/lambda; % in matlab my is used, instead of lambda, for some reason
    if(mu <= 0)
      error('parameter 1 must be bigger than 0');
    end
    f = size(n_length);
    for i =1:n_length
      if (x(i) < 0)
        f(i) = 0;
      else
        f(i) = mu * exp(- mu * x(i));
      end
    end

  case 'gam'
    % Gamma distribution : G(a,b)
    % uses Euler's gamma function as well
    a = parameters(1);
    b = parameters(2);
    if ( a <= 0 || b <= 0)
      error('parameters must be positive');
    end

    f = zeros(1,n_length);
    for i=1:n_length
      if (x(i) > 0)
        f(i) = 1 / (b^a) / gamma(a) * x(i)^(a-1)*exp(-x(i)/b);
      else
        f(i) = 0;
      end

  end


  case  't'
    % student distribution
    % S(n)

    n = parameters(1);
    if(n < 1)
      error('parmeter must be bigger than 1!')
    end

    f = zeros(1,n_length);
    for i=1:n_length
      f = gamma((n+1)/2) / sqrt(n*pi) / gamma(n/2) * (1.+(x.^2)/n).^(-(n+1)/2);
    end


  case 'unif'
    %U(a,b)
    a = parameters(1);
    b = parameters(2);
    % f = zeros(n_length);

    for i=1:n_length
      if(x(i) >= a & x(i) <= b)
        f(i) = 1 / (b - a);
      else
        f(i) = 0;
      end

    end
  case 'f'

    % Snedecor-Fisher distribution
    % F(m,n)
    m = parameters(1);
    n = parameters(2);
    if( m <= 0 || n <= 0)
       error('parameters must be bigger than zero!');
    end
    f = zeros(1,n_length);
    for i = 1:n_length
      if(x(i) >= 0)
          f(i) =  1 / beta(m/2,n/2) * (m/n)^(m/2) * x(i)^(m/2-1)/(1 + m/n*x(i))^((m+n)/2);
      else
        f(i) = 0;
      end
    end


end



end
