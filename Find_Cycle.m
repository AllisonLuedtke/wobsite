%This function takes a potential eq network adjacency matrix and finds the
%cycle

function [cycle] = Find_Cycle(A)

Pot = A; %Adjacency matrix

%How many firms?
num_firms = size(Pot,1);

g = digraph(Pot);

%Find the last edge
e = dfsearch(g, 1, 'edgetodiscovered', 'Restart', true);

%The supplier of this edge is the last firm in the cycle (the customer is
%the first)

last_node = e(2);
temp_cycle = zeros(num_firms,1);

temp_cycle(1) = last_node;

index = 1;

%next_node = find(Pot(:,last_node));
next_node = find(Pot(last_node,:));
while next_node ~= last_node
    %increment index
    temp_i = index;
    index = temp_i + 1;
    
    %store the next node in the next spot
    temp_cycle(index) = next_node;
    
    %get the supplier of the next node and made THAT the next node
    temp_node = next_node;
    next_node = find(A(temp_node,:));
   
end

%find index of last non-zero element of cycle
last = find(temp_cycle, 1, 'last');

reversed_cycle = temp_cycle(1:last);

cycle = flipud(reversed_cycle);


end