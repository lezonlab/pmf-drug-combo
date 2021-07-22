function h = HubBuilder(n,hubnum,edges,randflag,inputhubs)
    %Note, once assigning hubs, the remaing edges are assigned in a
    %Erdos-Reyni fashion. This is time consuming, so in the case that no
    %edges are hidden, we simply return a matrix of all ones
    flag = true;
    if edges == n^2
        h = ones(n,n);
        flag = false;
    end
    edges = edges - 26;
    if flag
        links = 0; %Counter for number of links in graph
        if randflag
            temp = randperm(n);
            hubs = temp(1:hubnum);
        else
            hubs = inputhubs;
        end
        %Create Hidden Ind
        Hidden_Ind = zeros(n);
        %First, ensure every node has at least one connection
        for col = 1:n %Iterate through every column
            %If the collumn has no connection yet, add one
            if nnz(Hidden_Ind(col,:))==0
                %Find all nodes that still have no connections, and select
                %one to connect it to
                check = all(Hidden_Ind==0);
                perm = find(all(Hidden_Ind==0));
                pos = randi(length(perm));
                add = perm(pos);
                while add == col %Ensure link to be added isn't on diagonal
                    pos = randi(length(perm));
                    add = perm(pos);
                end
                Hidden_Ind(col, add) = 1;
                Hidden_Ind(add, col) = 1;
                links = links + 1;
            end
        end
        %Now that every node has one connection, fully connect each hub
        for h = 1:length(hubs)
            overflow = nnz(Hidden_Ind(hubs(h),:)); %Account for links that are already connected
            Hidden_Ind(hubs(h),:) = 1;
            Hidden_Ind(:,hubs(h)) = 1;
            links = links + n - overflow - 1; %Subtracting 1 for double counting the index the row and collumn intersect on
        end
        %Remaining links are assigned randomly
        comp = Hidden_Ind;
        [row,col] = find(triu(Hidden_Ind==0));
        ix = randperm(length(row));
        total = (edges-hubnum*104);
        for i = 1:total
            %Ensure we are not adding along the diagonal
            if row(ix(i)) ~= col(ix(i))
                Hidden_Ind(row(ix(i)), col(ix(i))) = 1;
                Hidden_Ind(col(ix(i)), row(ix(i))) = 1;
                links = links + 1;
            end
        end
        h = Hidden_Ind;
        %fprintf(1,'links = %d\n', links);
    end
    
    %nnz(Hidden_Ind)
end