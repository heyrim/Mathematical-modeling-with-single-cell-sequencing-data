function Load_graph( ntest ) 
global  cDiffM  nVrtx nEdge edge vertex wUn indSig 

%% Load graph structure. Includes vertex, edge information 
    if( ntest == 0 )  
        %%%%%% Nestorowa Graph 8 
        load( 'data/Graph_Nestorowa_node8.mat' )
        vertex = centerdm(1:8,:); 
        indSig = [3:6,8];
        
        for nn = min(cluster):max(cluster)
            num(nn+1) = length(find(cluster==nn)); 
        end
        
    elseif( ntest == 21 )
        %%%%%% Nestorowa Graph 19      
        load( 'data/Graph_Nestorowa_node21.mat' ); 
        %%%% merge insignificant clusters. 
        cluster = cluster + 1; 
        cluster(cluster==20) = 8; cluster(cluster==21) = 17; cluster(cluster==22) = 2; 
        ind = find( cluster==14 ); cluster( ind(1:2:end) ) = 16; cluster( ind(2:2:end) ) = 7; 
        
        for n = 1:19; num(n) = length(find(cluster==n)); end 
        vertex = centerdm(1:19,:); 
        indSig = [19, 9, 6, 10, 18, 8, 17, 13, 5, 7, 16 ]; 

    elseif( ntest == 1 ) % Paul 
        load( 'data/Graph_Paul_G9.mat' )
        cluster(cluster==10)=6; 
        for n = 1:9; num(n) = length(find(cluster==n)); end     
        
        indSig = [2,3,4,8,9];   
    end 
    
    nVrtx = size( vertex, 1 ); 
    nEdge = size( edge, 1 );     
    
    BC.Nindout = zeros( 1, nVrtx ); BC.Nindin = zeros( 1, nVrtx ); 
    Cout = zeros(nVrtx); Cin = zeros(nVrtx);
    if( ntest == 0 ) 
        %%%%%% Graph 8 
        BC.inNout = union( [1], [] ); %3,5,6,7,8] );     
        BC.noin = [2]; 
        BC.noout = [3:8];     
    elseif( ntest == 21 )
        %%%%%% Graph 19      
        BC.noin = 15; 
        BC.noout = [19, 6, 10, 16, 17 ];     
        BC.inNout = setdiff( 1:19 , [BC.noin, BC.noout, 14 ] ); 
    elseif( ntest == 1 ) % Paul 
        BC.inNout = [1,7,5];      
        BC.noin = [6]; 
        BC.noout = [2,3,9,4,8];                 
    end 
    
%     cDiffM = diag( cDiff ); cDiffM = sparse( cDiffM ); 
    for n = 1:nEdge 
        len(n) = norm( vertex(edge(n,1),:) - vertex(edge(n,2),:) ); 
    end     
    wUn = ones(2,1)*len; 
    
%     cDiffM = sparse( diag( 1./(len.^2) / 5 ) ); 
    cDiffM = sparse( diag( 1./(len.^2) ) ); 
    
    for nn = BC.inNout 
        ind = find( edge(:,1) == nn ); 
        BC.Nindout(nn) = length(ind); 
        BC.indout(nn,1:BC.Nindout(nn)) =  ind; 
        
        ind = setdiff( ind, find(edge(:,3)==0)); 
        nind = edge(ind,2); 
        Cout(nn,nind) = num(nind)/sum(num(nind)); 
        
        ind = find( edge(:,2) == nn ); 
        BC.Nindin(nn) = length(ind); 
        BC.indin(nn,1:BC.Nindin(nn)) =  ind; 
        
        ind = setdiff( ind, find(edge(:,3)==0) ); 
        nind = edge(ind,1); 
        Cin(nind,nn) = num(nind)/sum(num(nind)); 
        
    end 
    for nn = BC.noin 
        ind = find( edge(:,1) == nn ); 
        BC.Nindout(nn) = length(ind); 
        BC.indout(nn,1:BC.Nindout(nn)) =  ind; 
        
        ind = setdiff( ind, find(edge(:,3)==0) ); 
        nind = edge(ind,2); 
        Cout(nn,nind) = num(nind)/sum(num(nind)); 
        
    end 
    %%%% no out at 3:8 
    for nn = BC.noout 
        ind = find( edge(:,2) == nn ); 
        BC.Nindin(nn) = length(ind); 
        BC.indin(nn,1:BC.Nindin(nn)) =  ind; 
        
        ind = setdiff( ind, find(edge(:,3)==0) ); 
        nind = edge(ind,1);
        Cin(nind,nn) = num(nind)/sum(num(nind)); 
        
        ind = find( edge(:,1) == nn ); 
        BC.Nindout(nn) = length(ind); 
        BC.indout(nn,1:BC.Nindout(nn)) =  ind; 
    end 


    for n = 1:nEdge 
        BC.rateout(n,1) = Cout( edge(n,1),edge(n,2) ); 
        BC.ratein( n,1) = Cin( edge(n,1),edge(n,2) ); 
    end 
    
    





end


