function [ uTeo ] = fTeo( x, c, t )

    uTeo = exp(-10*(4*(x-c*t) - 1).^2)';

end

