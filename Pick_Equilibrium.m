%This function determines which of the potential equilbrium networks is the
%pairwise stable equilbrium.

%It takes as input:
%   - The list of Potential Equilbrium networks, P
%   - The Payoffs for each firm for each such network, Pay
%   - The possible suppliers for each firm, Sup

%It returns a matrix describing the equilbrium, E


function [E] = Pick_Equilibrium(P,Pay,Sup)

%Get the function values
Pot_Matrix = P; %lists the potential eq networks
Payoffs = Pay;  %lists the payoffs for each firm for each network
Suppliers = Sup; %lists the suppliers for each firm

dim = size(Pot_Matrix);

%Create a boolean indicating whether you found a profitable deviation
found_dev = 0;


%How many options for eq network are there?
num_options = dim(1);

%How many firms are there?
num_firms = dim(2);

%initialize E as negative ones
E = -1*ones(1,num_firms);

%create vectors to hold the current and other eq network options
current = zeros(num_firms,1);
other = zeros(num_firms,1);

%create vectors to hold the payoffs for the two networks
current_payoffs = zeros(num_firms,1);
other_payoffs = zeros(num_firms,1);

%Starting from the first option, check if each matrix is pairwise stable
o = 1;
while o <= num_options
    %Get the current network
    current = Pot_Matrix(o,:);
    
    %Get the payoffs for the current network
    current_payoffs = Payoffs(o,:);
    
    %set the indicator to 0
    found_dev = 0;
    
    %For each firm, check to see if there is a deviation they would take
    f = 1;
    while f <= num_firms
        %What other suppliers does this firm have
        sup = nonzeros(Suppliers(f,:));
        num_sup = length(sup);

        
        %For each supplier that isn't the one being used right now,
        %construct a row describing the alternative matrix
        s = 1;
        while s <= num_sup
            %only if the supplier is different from the current one
            %current(f)
            %other = zeros(1,num_firms);

            if sup(s) ~= current(f)
                %copy current
                other = current;
                
                %put the alternative supplier option in for firm f in the
                %alternative network
                other(f) = sup(s);
                
                %find the index of the other vector in the potential
                %network matrix
                other_index = find(ismember(Pot_Matrix,other,'rows'));
                
                %get the payoffs associated with that network
                other_payoffs = Payoffs(other_index,:);
                Payoffs(other_index,:);
                %Payoffs
                
                %Are both the current firm and the potential new
                %supplier made better off?
                if other_payoffs(f) > current_payoffs(f) && other_payoffs(sup(s)) > current_payoffs(sup(s))
                    %set the indicator to 1, break the supplier loop
                    found_dev = 1;
                    %disp('found dev')
                    break;
                else
                    %otherwise move to the next supplier
                    %disp('did not find dev')
                    temp_s_2 = s;
                    s = temp_s_2 + 1;
                end 
            else
                %If it IS the current supplier, just go to the next
                %supplier
               temp_s = s;
               s = temp_s + 1;
            end        
        end
        
        %If you found a deviation, break the firm loop, you don't need to keep
        %looking
        if found_dev == 1
            break
        else
            %otherwise go on to the next firm
            temp_f = f;
            f = temp_f + 1;    
        end
    end   
    
    %If you found a deviation, increase s becuase the one you tried wasn't
    %an eq.
if found_dev == 1
    temp_o = o;
    o = temp_o + 1;
 %Otherwise you didn't find a deviation and this is a pairwise equilbirum
   %     set E to be current and put s at the end of it's max
else
    E = current;
    break
end
    
end

end
