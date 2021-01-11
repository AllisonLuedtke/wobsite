%This function creates an adjacency matrix, Adj, and a productivity
%matrix, Prod, and a flag, according to the following guidelines

%Adjacency Matrix:
%No self loops
%No mutliple eddges (allows for the possibility that a firm can draw
%another firm as a supplier twice, but since the z's will never be equal
%only the edge with the bigger z will be used.
%The number of suppliers is drawn from a discrete uniform distribution
%(1...num_firms-1)
%The identity of the supplier is drawn uniformly from the other firms.

%Productivity Matrix:
%If there is an edge in Adj, there is a non-zero number in Prod.
%This number is drawn from a truncated uniform distribution.

%flag
%The flag indicates that there is at least one firm with only one supplier

function [Adj, Prod, flag] = Create_A_and_Z(number_firms)

%Get function values
num_firms = number_firms;
num_firms_sq = num_firms*num_firms;

%create a sobol set from which to draw Z's
p = sobolset(1,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
%draw num_firms^2 z's to us
s_set = net(p,4*num_firms_sq);
%make them between 1 and two digits
z_set = s_set(s_set > 0.2);

%create blank A and Z
Adj = zeros(num_firms);
Prod = zeros(num_firms);

%create a blank flag;
flag = 0;

%create a matrix describing who can be a supplier for each firm
%use each row of this as the data sample set
poss_sup = zeros(num_firms, num_firms-1);
for r = 1:num_firms
   counter = 1;
   for c = 1:num_firms
      if c~=r
         poss_sup(r,counter) = c;
         temp_count = counter;
         counter = temp_count + 1;
      end
   end  
end
%poss_sup

%create a vector of discrete uniform random variables describing the number of
%suppliers for each firm. The number of suppliers is in the range
%(2...num_firms)
%temp_num_suppliers = unidrnd(num_firms - 3,num_firms,1);
%num_suppliers = temp_num_suppliers + 2*ones(num_firms,1);

%create a vector with the number of suppliers for each firm
num_suppliers = datasample(2:num_firms-1,num_firms);

rand_index = 1;

%populate the two for each firm (row)
for f = 1:num_firms
    %get the number of suppliers for that firm
    num_f_sup = num_suppliers(f);
    
    %draw that many from the available set
    %get row from sup_id_set for firm j
    sup_row = poss_sup(f,:);
    %draw num_f_sup from the possible suppliers without replacement
    sups = datasample(sup_row,num_f_sup,'Replace',false);
    
    %draw that many numbers from (1...num_firms - 1)
    %temp_sup_indices = unidrnd(num_firms-2,num_f_sup,1);
    %sup_indices = temp_sup_indices + ones(num_f_sup,1);
    
    %for each supplier
    for c = 1:num_f_sup
       %draw the next id from the list of ids
       %sup_index = sup_id_set(rand_index);
     
       %find the name of the supplier in the matrix of possible suppliers
       %supplier = poss_sup(f,sup_index);
       
       supplier = sups(c);
       
       %draw a z
       z = z_set(rand_index);
       %z = 1/denom;
       
       %increment z index
       temp_index = rand_index;
       rand_index = temp_index + 1;
       
       %store them
       %if there hasn't already been a supplier put there
       %if Adj(f,supplier) == 0
           %record the supplier and z
        %   Adj(f,supplier) = 1;
         %  Prod(f,supplier) = z;
       %else %if there has already been a supplier and thus z stored there,
            %check to see which z is biggest
        %   if z > Prod(f,supplier)
         %     Prod(f,supplier) = z; 
          % end
      % end 
      Adj(f,supplier) = 1;
      Prod(f,supplier) = z;
        
    end
    
    
    
end    
%How many suppliers did each firm show up with?
final_num_sup = sum(Adj~=0,2);
%Are any of them 1?
flag = ismember(1,final_num_sup);



end