function pair_diff = circ_dist2(X)
    % Given a 1D array of angles, find the 2D array of pairwise differences
    % Based on https://github.com/circstat/circstat-matlab/blob/master/circ_dist2.m
    % https://www.mathworks.com/matlabcentral/answers/67282-function-which-returns-the-outer-product-of-two-vectors

    x_prime = exp(sqrt(-1)*X);
    x_ones = ones(1,length(X));
    outer_product = [];
    for i=1:length(X)
        temp = [];
        for j=1:length(X)
            temp = [temp, x_prime(i)*x_ones(j)];
        end
        outer_product = [outer_product; temp];
    end
    % x = np.outer(np.exp(1j*X), np.ones(X.size)) # similar to meshgrid, but simpler to implement
    x = outer_product;
    % y = np.transpose(x)
    y = x';
    % return np.angle(x/y)
    pair_diff = angle(x./y);
end