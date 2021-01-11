%This function creates an adjacency matrix, Adj, and a productivity
%matrix, Prod, and a flag, according to the following guidelines

%Adjacency Matrix:
%No self loops
%No mutliple edges (allows for the possibility that a firm can draw
%another firm as a supplier twice, but since the z's will never be equal
%only the edge with the bigger z will be used.
%The number of suppliers is drawn from a Poisson(M) distribution and then
%shifted by 1.
%The identity of the supplier is drawn uniformly from the other firms.

%Productivity Matrix:
%If there is an edge in Adj, there is a non-zero number in Prod.
%This number is drawn from a Pareto(0.2, -1.8) distribution

%flag
%The flag indicates that there is at least one firm with only one supplier

function [Adj, Prod, flag] = PCreate_A_and_Z(number_firms)

%Get function values
num_firms = number_firms;
num_firms_sq = num_firms*num_firms;

%set M
M = 3;

%create a sobol set from which to draw Z's and uniforms for the poisson
p = sobolset(1,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
p2 = sobolset(1,'Skip',1e3,'Leap',1e2);
p2 = scramble(p2,'MatousekAffineOwen');
%draw num_firms^2 z's to us
s_set = net(p,4*num_firms_sq);
%make them Pareto
s_set_2 = s_set./0.2;
z_set = s_set_2.^(-1.8);

%sup_set
sup_set = net(p2,10*num_firms_sq);

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

%create a vector containing Poisson(M) rv's for each firm describing how
%many suppliers they have
it = 1;
num_suppliers = zeros(num_firms,1);
for g = 1:num_firms
   %set num options to 1
   num_sup_op = 1;
   prod = 1;
   poiss_flag = 1;
   while poiss_flag == 1
      %get next value from the sup_set
      u = sup_set(it);
      %incrememnt it
      temp_it = it;
      it = temp_it + 1;
      %set prod = u*prod
      temp_prod = prod;
      prod = u*temp_prod;
      %prod
      if prod < exp(-M)
          poiss_flag = 0;
      else
          %increment num_sup_op
          temp_sup_op = num_sup_op;
          num_sup_op = temp_sup_op + 1;
      end    
   end
   
   %num_sup_op is the number of supplier options for firm g
   %store it
   num_suppliers(g) = num_sup_op;
    
end
%num_suppliers



%create a vector with the number of suppliers for each firm
%num_suppliers = datasample(2:num_firms-1,num_firms);

rand_index = 1;

%populate the two for each firm (row)
for f = 1:num_firms
    %get the number of suppliers for that firm
    num_f_sup = num_suppliers(f);
    
    %get row from sup_id_set for firm f
    sup_row = poss_sup(f,:);

    

    
    %for each supplier
    for c = 1:num_f_sup
       %draw an id from the sup row
       sup_id = datasample(sup_row,1);
       
       
       %draw a z
       z = z_set(rand_index);
       %z = 1/denom;
       
       %increment z index
       temp_index = rand_index;
       rand_index = temp_index + 1;
       
       %store them
       %if there hasn't already been a supplier put there
       if Adj(f,sup_id) == 0
           %record the supplier and z
           Adj(f,sup_id) = 1;
           Prod(f,sup_id) = z;
       else %if there has already been a supplier and thus z stored there,
            %check to see which z is biggest
           if z > Prod(f,sup_id)
              Prod(f,sup_id) = z; 
           end
       end 
        
    end
    
    
    
end    
%How many suppliers did each firm show up with?
final_num_sup = sum(Adj~=0,2);
%Are any of them 1?
flag = ismember(1,final_num_sup);



end