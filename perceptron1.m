clear all;
% Params
b = 1;                                      % threshold
T=10;                                       % maximum number of iterations
alpha=1;                                    % learning rate
toll = 10^(-5);                             % tollerance for stop criterion

% Init
x = [[0;0], [0;1], [1;0], [1;1]];           % matrix of inputs
[n,p] = size(x);                            % n: dimension of inputs, p: number of patternss
figure; hold on;

% Two taks
for task = 1:2
    if task ==1; y_t = [0, 1, 1, 1];end           %desired output OR
    if task ==2; y_t = [0, 1, 1, 0];end           %desired output: XOR
    
    % Init
    w = zeros(n,1);                             %initialise weights
    y = 0;                                  %initialise output
    error = zeros(1,T);                         %initialise vector error
    
    % Main algo
    for t =1:T
        for mu = 1:p
            y = (w'*x(:,mu)-b)>0;       %calculate output
            w = w + alpha*(y_t(mu)-y)*x(:,mu);   %update weights
            error(t) = error(t) + ((y_t(mu) - y).^2); % save the weights
        end
        if error(t)< toll
            break                                   %stop if error is close to 0
        end
    end
    
    % Plot
    plot(error)
end
xlabel('iterations')
ylabel('error')
legend('OR', 'XOR')