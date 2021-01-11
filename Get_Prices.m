%This function determines the input prices (price of using the output of a
% firm in your own production) that each firm pays. It takes as inputs:
%       - an adjacency matrix describing the equilibrium network, Adj
%       - a productivity vector, Prod, which lists the productivities of the
%       edges pointing to each firm
%       - the wage
%
% It returns a vector p_e which describes the price of the edge pointing to
% firm i.



function [p_e] = Get_Prices(Adj, Prod, wage)

%Get the function values
A = Adj;
z = Prod;
w = wage;

%How many firms are there?
num_firms = size(A,1);

%Create the price vector
p_e_temp = zeros(num_firms,1);

%assign alpha and wage here it'll probably be somewhere else in the big program
a = 0.33;
%w = 10;

%Get the cycles
cycles = Find_Cycle_For_Test(A);
num_cyc = size(cycles,1);

%for each cycle
for cyc = 1:num_cyc

%find the length of this cycle
cyc_row = cycles(cyc,:);
cycle = nonzeros(cyc_row);
c = length(cycle);

%construct the root node indicator vector
root_node = zeros(c,1);
%for each firm in the cycle list
for f = 1:c
    %get the firm index
    firm = cycle(f);
    %how many customers does that firm have
    cyc_customers = find(A(:,firm));
    num_customers = length(cyc_customers);
    %If that firm has more than one customer it is a root node
    if num_customers > 1
        root_node(f) = 1;
    end   
end
%root_node


%First, calculate the price for the first firm in the cycle and
%simultaneously construct the sum of alpha^k
price_base = 1;
alpha_sum = 0;
for i = 1:c
    temp_price = price_base;
    temp_alpha = alpha_sum;
    %get the index of the firm in the cycle list
    firm_index = cycle(i);
    z_i = z(firm_index);
    price_base = temp_price/((z_i)^(a^(c-i)));
    alpha_sum = temp_alpha + (a^(i - 1));
end   
%alpha_sum;

%multiply by wage function
w_expo = (1-a)*alpha_sum;
price_base_w = price_base * (w^(w_expo));
%price_base_w;
first_price = price_base_w^(1/(1 - a^c));

%put that price in the price vector for the first firm in the cycle
first_index = cycle(1);
p_e_temp(first_index) = first_price;

%find the rest of the prices of the cycle firms
%---> for each firm past the first one, calculate the MC of the prev firm
%and save that value of the price for the current firm

for j = 2:c
    %Get the index of the current firm
    firm_index = cycle(j);
    
    %Get the index of the previous firm in the cycle and its price
    prev_index = cycle(j-1);
    prev_price = p_e_temp(prev_index);
    
    %Calculate the MC of the previous firm (=price of current firm)
    MC_prev = (1/z(prev_index))*(prev_price^a)*(w^(1-a));
    
    %Store the MC of the prev firm as the price of the current firm
    p_e_temp(firm_index) = MC_prev;   
    
    
end

%Find the prices on the branches using a recursive tree search algorithm
%Search through the indicator vector, and for every firm that is a root
%node 

for r = 1:c
    %if it's a root node
    if root_node(r) == 1
        %get its index from the cycle list
        root_index = cycle(r);
        
        %get the customers of the root node, but remove the next guy on the
        %cycle
        customers = find(A(:,root_index));
        
        %cycle customer index
        %if it's the last guy on the cycle, the cycle customer you want to
        %remove is the FIRST guy on the cycle
        if r == c
            cyc_cus_index = 1;
        else
            %otherwise it's the next guy in the cycle list.
            cyc_cus_index = r+1;
        end
        
        cyc_cus = cycle(cyc_cus_index);
        %remove the cycle customer
        branch_customers = customers(customers~=cyc_cus);
        
        
        %call the recursive function for the root node
        p_e_temp = rec_branch_price(root_index,branch_customers,A,z,p_e_temp,a,w);
        
    end
    
end

    p_e = p_e_temp;
    
end

end
