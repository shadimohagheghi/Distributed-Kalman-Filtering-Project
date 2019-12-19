
% This function plots a network topology

% April 2, 2014 (A. H. Sayed, UCLA).

function graphPlot(Adjacency,Coordinates,Color)

%INPUT
% Adjacency: size NxN; A[a,b] =  1 if a and b are connected; otherwise zero.
%
% Coordinates: Nx2 matrix containing the (x,y) location coordinates of the agents in the square region [0,1.2]x[0,1.2]; 
%              each row corresponds to one agent
%
% Color: a vector of size Nx1. Location k in this vector is set to:
%        0 if the corresponding agent should have one color (yellow)
%        1 if the corresponding agent should have a second color (red)
%        2 if the corresponding agent should have a third color (green)

A = Adjacency;      % adjacency matrix
N = max(size(A));   % number of agents

x_coordinates = Coordinates(:,1); % x-coordinates of agents
y_coordinates = Coordinates(:,2); % y-coordinates of agents

figure
hold on
for k=1:N
  for l=1:N
    if A(k,l)>0
      plot([x_coordinates(k),x_coordinates(l)],[y_coordinates(k),y_coordinates(l)],'b-','LineWidth',1.5);
    end
  end
end

for k=1:N
 if Color(k) == 0 % yellow
     plot(x_coordinates(k),y_coordinates(k),'o','MarkerEdgeColor','b','MarkerFaceColor','y','MarkerSize',10);
 else
   if Color(k) == 1 % red
      plot(x_coordinates(k),y_coordinates(k),'o','MarkerEdgeColor','b','MarkerFaceColor','r','MarkerSize',10);
   else % green
      plot(x_coordinates(k),y_coordinates(k),'o','MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',10);
   end
 end
end

axis([0,1.2,0,1.2]);
axis square
grid
xlabel('x-coordinate')
ylabel('y-coordinate')

for k=1:N
  text(x_coordinates(k)+0.03,y_coordinates(k)+0.03,num2str(k),'Fontsize',7);
end  
