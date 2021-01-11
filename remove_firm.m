%This function deletes a specified firm from a Adj and a Z matrix

%It takes as inputs an Adjacency Matrix, old_A, a productivity matrix,
%old_Z, and a firm number, firm.

%It returns a new Adj, new_A and a new productivity matrix, new_Z.

%If old_A and old_Z were both n x n, then the new ones are (n-1) x (n-1)

function [new_A, new_Z] = remove_firm(old_A, old_Z, firm)

%get function values
A = old_A;
Z = old_Z;
f = firm;

%copy A
new_A = A;

%copy Z
new_Z = Z;

%remove firm f's rows in A and Z
new_A(:,f) = [];
new_Z(:,f) = [];

%remove columns
new_A(f,:) = [];
new_Z(f,:) = [];

end