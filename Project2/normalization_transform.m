function H = normalization_transform(P)
mu = mean(P, 1);
sigma = std(P, 0, 1);
H = [1./sigma(1), 0, -mu(1)./sigma(1); 0, 1./sigma(2), -mu(2)./sigma(2); 0, 0, 1];
end