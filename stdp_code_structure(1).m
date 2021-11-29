clear all;
dt = 1;          % Euleur integration time step [ms]
p = 60;          % number of pairings
rep = 1000;      % time between the pairing [ms], i.e. 1Hz

% parameters STDP
tau_plus = 10;  % LTP time constant [ms]
tau_minus = 20; % LTD time constant [ms]
A_plus = 1;    % LTP learning rate or amplitude
A_minus = 1;   % LTD learning rate or amplitude

lag_range = -50:5:50; % simulate for different lags between the pre and post spikes
dw_range = zeros(size(lag_range)); % save weight changes
count = 0; % counter

for lag = lag_range
    % Init
    T = rep*(p-1)+2*abs(lag)+1; % spiketrain total time length
    pre_spikes = zeros(1,T);% presynaptic spiketrain
    post_spikes = zeros(1,T);% postsynaptic spiketrain
    pre_spikes(lag+abs(lag)+1:rep:lag+abs(lag)+T) = 1;  % presynaptic spiketrain (0  if not spike, 1 if spike)
    post_spikes(abs(lag)+1:rep:T+abs(lag)) = 1;  % postsynaptic spiketrain (0  if not spike, 1 if spike)
    count = count +1; % update a counter
    x = zeros(1,T);   % presynaptic trace
    y = zeros(1,T);   % postsynaptic trace
    dw = 0;           % synaptic weight change
    
    % time iterations
    for t = 1:T
%%%%%%%%%%%%%%%  TO ADD HERE THE STDP CODE
%%%%%%%%%%%%%%% Write a code which calcules the weight change dw
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
    end
    dw_range(count) = dw;
end

figure;
plot(lag_range, dw_range)
xlabel('post spike time - pre spike time')
ylabel('weight change')