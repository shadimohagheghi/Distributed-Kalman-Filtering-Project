
function DFK()

N = 15; % number of agents

Iterations  = 100;
Samples = 50; % samples for averaging.

T = 0.1; % sampling period

%state space model

m=4; % size 

F =[ 1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
 
H = [1 0 0 0; 0 1 0 0];
G = zeros(4,1);
Q = 0;
d = -[0 0 0 1]';

%initial state vector
vx(1) = 10*cos(pi/3); % initial x-speed
vy(1) = 10*sin(pi/3); % initial y-speed
x(1)  = 5; % initial x-location
y(1)  = 20; % initial y-location

% noise variances across agents
sigma = 0.75*rand(N,1);

% Generating the Graph

algebraic_connectivity = 0;
while algebraic_connectivity < 0.01 % keep trying until a connected topology is chosen
    [A,SpectralRad,algebraic_connectivity,degree,coordinates] = graphProperties(N);
end

num_nb = zeros(N,1);
for k=1:N
    num_nb(k) = sum(A(k,:));
end

W = A;
for k = 1:N
    W(k,:) = W(k,:)/sum(W(k,:));
end

%metropolis
% for k=1:N
%     for l=1:N
%         W(k,l) = A(k,l)/max([num_nb(k), num_nb(l)]);
%     end
%     W(k,k) = 1 + W(k,k) - sum(W(k,:));
% end

% Generating data
z = zeros(2,Iterations,N); 
for i=1:N 
    z(1,1,i) = x(1) + sqrt(sigma(i))*randn; % first measurement of node i
    z(2,1,i) = y(1) + sqrt(sigma(i))*randn; % second measurement of node i
    R(:,:,i) = sigma(i)*eye(2,2); % measurement noise covariance matrix of node i
end

% calculating measurements.
for n = 1:Iterations
    vy(n) = vy(1) - 10*n*T;
    vx(n) = vx(1);
    x(n) = x(1) + vx(1)*n*T;
    y(n) = y(1) + vy(1)*n*T -0.5*10*(n*T)^2;
    for k=1:N
        z(1,n,k) = x(n) + sqrt(sigma(k))*randn;
        z(2,n,k) = y(n) + sqrt(sigma(k))*randn;
    end
end

% DKF Algorithm
MSE = zeros(1,Iterations);
x_ave  = zeros(1,Iterations);
y_ave  = zeros(1,Iterations);

for L=1:Samples
    x_local = zeros(m,N);
    x_plus = zeros(m,N);
    x_minus = zeros(m,N);
    P_plus = zeros(m,m,N);
    P_minus = zeros(m,m,N);
    
    for k=1:N
        x_minus(1:m,k) = zeros(m,1);
        P_minus(1:m,1:m,k) = eye(4,4);
    end
    
    for i=1:Iterations
        for k=1:N
            x_local(:,k) = x_minus(:,k);
            P_local = P_minus(1:m,1:m,k);
            for l=1:N
                if W(l,k) >0
                    error = z(:,i,l) - H*x_local(:,k);
                    K = P_local*(H')*inv(R(:,:,l) + H*P_local*(H'));
                    x_local(:,k) = x_local(:,k) + K*error;
                    P_local = P_local - K*H*P_local;
                end
            end
            P_plus(:,:,k) = P_local;
            P_minus(:,:,k) = F*P_plus(:,:,k)*(F') + G*Q*(G');
        end
        delta = 0.01;
        for k=1:N
            x_plus(:,k) = zeros(m,1);
            for l=1:N
                if l==k
                    x_plus(:,k)  = x_plus(:,k)  + (1+delta-degree(k)*delta)*x_local(:,k);
                else
                    x_plus(:,k)  = x_plus(:,k)  + delta*x_local(:,l);
                end
            end
            x_minus(:,k)  = F*x_plus(:,k) + d;
            
            x_ave(i) = x_ave(i) + x_plus(1,k);
            y_ave(i) = y_ave(i) + x_plus(2,k);
            
            MSE(i) = MSE(i) + (norm([x_plus(1,k)-x(i); x_plus(2,k)-y(i)],2)^2);
            
        end
    end
end

x_ave = x_ave/(Samples*N);
y_ave = y_ave/(Samples*N);

MSE = 10*log10(MSE/(Samples*N));

figure
subplot(211)
plot(x,y,'r', x_ave,y_ave,'*');
xlabel('x-coordinate');
ylabel('y-coordinate');
legend('real value','consensus estimate');
%title('algebraic_connectivity=',algebraic_connectivity, 'SpeactralRad=',SpectralRad)
'algebraic_connectivity=',algebraic_connectivity, 'SpeactralRad=',SpectralRad
axis([5 15 10 30]);
grid

subplot(212)
plot(1:Iterations,MSE,'b');
xlabel('iteration');
ylabel('MSE (dB)');
%algebraic_connectivity
axis([0 30 -15 45]);
grid

end
