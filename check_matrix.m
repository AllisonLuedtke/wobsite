%This function checks a matrix to see if any of the firms have no suppliers

%It takes as input an adjacency matrix, Adj

%It returns a binary flag that is 1 if any of the firms have no suppliers,
%flag, and the index of the first firm with no suppliers.

function [flag,index] = check_matrix(Adj)

%get function value
A = Adj;

%initialize flag and index
flag = 0;
index = 0;

%get the number of suppliers
num_sups = sum(A~=0,2);

%are any == 0?
%if any of them are 0, get the first index
if ismember(0,num_sups)
    %set the flag = 1
    flag = 1;
    %get the first index
    sup_list = find(num_sups == 0);
    index = sup_list(1);   
end

end