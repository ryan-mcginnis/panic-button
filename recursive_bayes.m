function Po = recursive_bayes(freqs, prior, x)
% Inputs: 
% freqs - vector of frequencies (Hz)
% prior - probabilty density function of HR based on frequency analysis 
% x     - hr measurements from time-domain analysis

% Outputs:
% Po    - probability density function from recursive bayesian estimation

% Written by Ryan S. McGinnis - ryan.mcginis14@gmail.com - December 11, 2016

% Copyright (C) 2016  Ryan S. McGinnis


%Initialize prior with input
Porig = prior/sum(prior); % Scale so that each value is now 1/(number of states), so it all sums to one.
Pmax = max(Porig);
Po=Porig; %duplicate for iterative update


%Iterative bayes
% Ko=0.0278; % covariance matrix for making a gaussian (assumes SD of 10 bpm) - could scale covariance by original probability
Ko=0.0069; % covariance matrix for making a gaussian (assumes SD of 5 bpm) 
for n=1:length(x)
    [~,ind]=min(abs(freqs-x(n)));
    K = Ko * (1 + ((Pmax-Porig(ind))/Pmax)); %update the covariance based on the original pdf (from frequency-based approach)
    Pr=Po; %store the posterior to the prior.   
    %likelihood
    m=0*Pr; %preinitialize the likelihood 
    for i=1:length(Pr)
       me = freqs(i);
       m(i) = 1/sqrt((2*pi)^2*det(K)) * exp(-(x(n)-me)'*(1/K)*(x(n)-me)/2); %Compute likelihood         
       m(i) = m(i) * Pr(i); % Combine this likelihood with the prior    
    end
    Po=m/sum(m); %normalize distribution to make it a proper probability distribution.
end
end