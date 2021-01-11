%This function builds the list of other networks for a given coalition and
%current network by calling the recursive function, rec_other_list

%It takes as input:
%   -The current network
%   -The coalition
%   -The list of suppliers for each firm

%It returns a list of networks (the suppliers for each firm) (so one row is
%a network)

function [Other_List] = build_other_list(current_network, co, Sup)

%get function variables
current = current_network;
coalition = co;
Suppliers = Sup;

%how many firms?
num_firms = length(current);

%how many other networks are there ( = number of rows of the Other_List)
co_size = length(coalition);
num_rows = 1;
%for each firm in the coalition, take the product of the number of sups - 1
for f = 1:co_size
    %get the firm
    co_firm = coalition(f);
    %get the number of suppliers
    num_sups = nnz(Suppliers(co_firm,:));
    %if it's 1, mult by 1
    %if it's greater than 1, mult by it - 1
    if num_sups > 1
       %mult by  num_sups - 1
       temp_num_rows = num_rows;
       num_rows = temp_num_rows*(num_sups - 1);
    end   
end

%create storage matrix
temp_store = zeros(num_rows, num_firms);

%start with the current network
init_vec = current;

%call the recursive function for 0
other_mat = rec_build_other(0, coalition, current, init_vec, Suppliers, temp_store);

Other_List = other_mat;


end
