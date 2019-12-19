
function [A,SpectralRad,Algebraic_Connectivity,Degree_Vector,Coordinates] = graphProperties(Num_nodes)

rand('seed',1);

N = Num_nodes;  % Number of nodes.
A = zeros(N,N); % Adjacency matrix.
L = zeros(N,N); % Laplacian matrix.

radius = 0.5;
r = radius;  %neighbors are defined by nodes being within this radius

% Generating N random (x,y) coordinates

x_coordinates = rand(1,N) + 0.1;

y_coordinates = rand(1,N) + 0.1;
Coordinates = [x_coordinates' y_coordinates'];

% Finding the adjacency matrix.

  for k=1:N
    for l=1:N
        d = sqrt((x_coordinates(1,k)-x_coordinates(1,l))^2 + (y_coordinates(1,k)-y_coordinates(1,l))^2);
        if d <= r
            A(k,l) = 1; % set entry in adjacency matrix to one if nodes k and l should be neighbors.
        end
    end
  end

% Listing the degrees of nodes
num_nb = zeros(N,1);
for k=1:N
    num_nb(k) = sum(A(k,:));
end
Degree_Vector = num_nb;  % vector of degrees for the various nodes

% Laplacian matrix L and verifying the connectivity of the graph

for k=1:N
    L(k,k) = max(0, sum(A(k,:))-1); % set diagonal entry to zero if degree-1 for node k is negative.
    for l=k+1:N
        L(k,l) = -1*A(k,l);
        L(l,k) = -1*A(l,k);
    end
end
sigma = svd(L); % vector of singular values of L.

Algebraic_Connectivity = sigma(N-1); % algebraic connectivity
SpectralRad=max(eig(A));
