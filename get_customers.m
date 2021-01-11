function [ cus ] = get_customers( Adjacency, index )
%This function determines the number and identity of the customers
% of firm index

A = Adjacency;
i = index;

A_size = size(A);
num_firms = A_size(1);

cus = [];

count = 0;

for row = 1:num_firms
    if A(row,i) == 1
        count = count + 1;
        cus(count) = row;
    end
end
end

