%This function returns the coalition-proof equilibrium network (if one
%exists and a list of -1's if one does not). It takes as input:
%   - The list of possible equilibrium matrices
%   - The list of associated payoffs
%   - The sup matrix that lists the possible suppliers for each firm

%It returns a network, e, that lists the supplier of each firm (or -1's)

function [E] = Find_CP_Eq(Pot_Mat, Pay, Sup)

%get the function variables
Pot = Pot_Mat; %The list of possible eq matrices
Payoffs = Pay; %The list of payoffs associated with each such matrix
Suppliers = Sup; %Lists the possible suppliers for each firm

%get the dimensions of the pot matrix
pot_dim = size(Pot);
num_opt = pot_dim(1);
num_firms = pot_dim(2);

%Initialize E as negative ones in case you don't find an eq
E = -1*ones(num_firms,1);

%create a counter to count the matrices visited
m = 1;
%while there are still matrices to check
%EQ MATRIX LOOP
while m <= num_opt
    
    %for each matrix, check if it is coalition proof
    %   --> check each coaliton size
    %get the current possible eq network we're trying
    poss_eq = Pot(m,:);
    
    %COALITION SIZE LOOP
    %check each coalition size
    %if any of them have a profitable deviation, set coaliton_flag = 1 and found_dev = 1 and end
    coalition_size = num_firms;
    coalition_flag = 0;
    found_dev = 0;
    while coalition_flag == 0
        %check each size
        %check the coalitions of the current coalition_size
        %coalition_size
        found_dev = check_coalition(poss_eq, m, Pot, Payoffs, Suppliers, coalition_size);
        %coalition_size
        
        %if there was a profitable deviation, give up on this network
        if found_dev == 1
            coalition_flag = 1;
        else
            %decrease the coalition size
            %if it's 1, set the flag to 1 to to end the loop
            if coalition_size == 1
                coalition_flag = 1;
            else
                temp_co_size = coalition_size;
                coalition_size = temp_co_size - 1;
                %coalition_size
            end
        end
    %END COALITION CHECKING LOOP
    end
    
    %if found_dev = 1, you found a deviation and need to increase the m
    %counter to go to the next network
    if found_dev == 1
        %go to next poss eq network
        temp_m = m;
        m = temp_m + 1;
    else
        %if you didn't find a deviation that means you got through all the
        %coalition sizes and no profitable deviations
        %   --> this is a coalition proof network
        %set E to this network
        E = poss_eq;
        %set the m counter to num_opt + 1 to end the big loop
        m = num_opt+1;  
    end 
    %If you never find a deviation, you never set E to anything different
    %so you return the -1's
    
%END OF MATRIX LOOP    
end


end