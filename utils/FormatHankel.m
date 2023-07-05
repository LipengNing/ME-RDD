function [H1, H2] = FormatHankel(x)

n = length(x);
if(mod(n,2))%n odd but max moment order is even
    H1 = hankel(x(1:(n+1)/2),x((n+1)/2:end));
    H2 = hankel(x(2:(n+1)/2)-x(3:(n+3)/2),x((n+1)/2:end-1)-x((n+3)/2:end));
else
    H1 = hankel(x(2:(n/2+1)),x((n/2+1):end));
    H2 = hankel(x(1:n/2)-x(2:(n/2+1)),x(n/2:end-1)-x((n/2+1):end));
end
    