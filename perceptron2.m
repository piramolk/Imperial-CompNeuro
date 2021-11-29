clear all

% Params
n = 100;                %number of inputs
p = 50;                 %length of a spike train
b = 1;                  %threshold
T = 50;                 %maximum number of iterations
alpha=1;                %learning rate
toll = 10^(-5);         %tollerance for stop criterion

% Init
X = round(rand(n,p));   %matrix of inputs
y_t = round(rand(1,p));   %desired output
w = ones(n,1)*0.5;      %initialise weights
y = 0;              %initialise output
error = zeros(1,T);     %initialise vector error

% Main algo
for t =1:T
    
    for mu = 1:p
        y = (w'*X(:,mu)-b)>0;       %calculate output
        w = w + alpha*(y_t(mu)-y)*X(:,mu);   %update weights
        error(t) = error(t) + ((y_t(mu) - y).^2); % save error
    end
    
    if error(t)< toll
        break                                   %stop if error is close to 0
    end
    
end

% Plot
figure;
plot(error)
xlabel('iterations')
ylabel('error')