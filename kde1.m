% KDE function with Gaussian kernel
function [a] = kde1(m, s, xmin, xmax, xint)
x = xmin:xint:xmax;
n = length(m);
f = zeros(n,length(x));
for i = 1:n;
	f(i,:) = (1./ (s(i)*sqrt(2*pi)) .* exp (  (-((x-m(i)).^2)) ./ (2*((s(i)).^2))  ).*xint); % Gaussian
end
a = (sum(f))/n; %sum and normalize

end