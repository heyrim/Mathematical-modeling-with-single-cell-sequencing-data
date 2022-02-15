
function [BC, steady, indSig] = edge_graph( nG ) 
global cAdv cDiff cDiffM cAdvM cRct cRctM addcRct N1 nVrtx nEdge edge vertex 
global wUn cDthM cAdvM2 
%% Nestorowa G8 
    if( nG == 0 )  
        %%%%%% Graph 8 
        load( 'data/Graph_Nestorowa_node8.mat' )
        vertex = centerdm(1:8,:); 
        indSig = [3:6,8];
 
    elseif( nG == 1 ) % Paul 
        load( 'data/Graph_Paul_G9.mat' )
        M(M==10)=6; 
        
        indSig = [2,3,4,8,9];   
    end 
    
    for n = 1:max(M); num(n) = length(find(M==n)); end     
    
    nVrtx = size( vertex, 1 ); 
    nEdge = size( edge, 1 );     
    
    BC.Nindout = zeros( 1, nVrtx ); BC.Nindin = zeros( 1, nVrtx ); 
    Cout = zeros(nVrtx); Cin = zeros(nVrtx);
    if( nG == 0 ) 
        %%%%%% Graph 8 
        BC.inNout = union( [1], [] ); %3,5,6,7,8] );     
        BC.noin = [2]; 
        BC.noout = [3:8];     
        
    elseif( nG == 1 ) % Paul 
        BC.inNout = [1,7,5];      
        BC.noin = [6]; 
        BC.noout = [2,3,9,4,8];                 
    end; 
    
%     cDiffM = diag( cDiff ); cDiffM = sparse( cDiffM ); 
    for n = 1:nEdge 
        len(n) = norm( vertex(edge(n,1),:) - vertex(edge(n,2),:) ); 
    end 
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
    
%     if( nargin==0 ) 
%     %%%%%% Graph 8 
% %     p2 = num( 3:8 )/ sum(num(3:8)); 
% %     wUn = zeros(2,nEdge); wUn(:,1) = 1; 
% %     wUn(:,[2,3,5,7,9,11]) = [1;1]*p2; 
% %     for n = [4,6,8,10,12] 
% %         wUn(1,n) = p2( edge(n,1)-2 ); 
% %         wUn(2,n) = p2( edge(n,2)-2 ); 
% %     end 
%     end 
    wUn = ones(2,1)*len; 

    for n = 1:nEdge 
        BC.rateout(n,1) = Cout( edge(n,1),edge(n,2) ); 
        BC.ratein( n,1) = Cin( edge(n,1),edge(n,2) ); 
    end 
    
    
%% based on steady state 
global Dx  
    cRctM  = zeros( N1, nEdge ); 
    cAdvM  = zeros( N1, nEdge ); 
    cAdvM2 = zeros( N1, nEdge ); 
    addcRct = sparse(zeros( N1, nEdge )); 
    
% load( 'steady.mat' ) 
% steady = steady(:,[1,2,3,5,7,9,11]); 

% global sig 
sig = 0.3; 
xx = [0:1/(N1(1)-1):1]; 
steady = zeros( N1, nEdge ); 
Mtmp = zeros( nVrtx ); 
for n = 1:nVrtx
  tmp = zeros( N1, nEdge ); 
  ind = find( edge(:,1) == n ); 
  tmptmp = Normal( 0, sig, xx )'; 
  tmp(:,ind) = tmptmp *ones(1,length(ind));   
    
  ind = find( edge(:,2) == n ); 
  tmptmp = Normal( 1, sig, xx )'; 
  tmp(:,ind) = tmptmp *ones(1,length(ind));   
  
  Mtmp(:,n) = plot_graph1d_half( tmp );   
end 

if( nG == 0 ) 
        ntmp = Mtmp\num(1:nVrtx)';
elseif( nG == 21 )  
        ind = [1:13,15:19]; 
        ntmp(ind) = Mtmp(ind,ind)\num(ind)';
elseif( nG == 1 ) 
        ntmp = Mtmp\num';        
end 
    
tmptmp(:,1) = Normal( 0, sig, xx )'; 
tmptmp(:,2) = Normal( 1, sig, xx )'; 
for n = 1:nEdge 
    steady(:,n) = steady(:,n) + tmptmp(:,1)*ntmp(edge(n,1))+ tmptmp(:,2)*ntmp(edge(n,2));
end 

% steady(:,1) = Normal( 0, sig, xx )*num(edge(1,1)) + Normal( 1, sig, xx )*num(edge(1,2))/2;  
% for n = 2:7 
% %   steady(:,n) = Normal( 1, sig, xx )*num(edge(n,2))/p2(n-1) + Normal( 0, sig, xx )*num(edge(n,1))/2;  
%   steady(:,n) = Normal( 1, sig, xx )*sum(num(3:8)) + Normal( 0, sig, xx )*num(edge(n,1))/2; 
% end 
% 
% steady(end,1) = -Dx(end,1:end-1)*steady(1:end-1,1) / Dx(end,end); 
% steady(1,2:7) = steady(end,1); 
% % for n = 2:7; steady(2,n) = -Dx(1,[1,3:end])*steady([1,3:end],n)/Dx(1,2); end 
% ntmp = plot_graph1d_half( steady ); 
% steady = steady / ntmp(1) * num(1); 


ntmp = Compute_nTotC( steady );  
steady = steady / ntmp; 

vv = -log(steady); 
for n = 1:nEdge 
%     cAdvM(:,n) = -cDiff(1) * Dx * vv(:,n); 
    
    cAdvM(:,n) = -cDiff(1) * Dx * vv(:,n) * cDiffM(n,n); 
    cAdvM(:,n) = cAdvM(:,n); % / norm( vertex(edge(n,1),1:2) - vertex(edge(n,2),1:2) ); 
end 

cProlif = [0.002152  0.01125  0.05658  0.1612  0.3199  0.6931]; 
cProlif = cProlif( [2 3 4 6] ); 
cAdvec = [0.77 0.7689 0.7359 0.66 0.154]; 
cAdvec = cAdvec(1:4); 

if( nG == 0 ) 
    cellind2state = [2 1 3 4 4 3 2 3]; 
    cProlif = cProlif( cellind2state ); 
    cAdvec = cAdvec( cellind2state ); 
    cRct = cProlif; 
    cAdv = cAdvec; 
    
elseif( nG == 21 )
    cellind2state([15]) = 1;     cellind2state([1]) = 2; 
    cellind2state([6,10,16]) = 4; 
    cellind2state([19,13,12,17,18,8]) = 3; 
    ind = [15,1,6,10,16,19,13,12,17,18,8]; 
    cRct(ind) = cProlif( cellind2state(ind) ); 
    cAdv(ind) = cAdvec( cellind2state(ind) ); 
    cRct([11,4]) = interp1( [0,1], cRct([15,1]), len([1,2])/sum(len(1:3)) ); 
    cRct([3]) = interp1( [0,1], cRct([1,6]), len([6])/sum(len([6,19])) );
    cRct([7]) = interp1( [0,1], cRct([13,16]), len([33])/sum(len([33,35])) );
    cRct([2]) = interp1( [0,1], cRct([3,12]), len([12])/sum(len([12,9])) );
    cRct([5]) = interp1( [0,1], cRct([17,7]), len([6])/sum(len([6,19])) );
    cRct([9]) = interp1( [0,1], cRct([3,19]), len([21])/sum(len([21,36])) );
%     for n = setdiff( [1:13,15:19], ind ) 
%         n1 = find(edge(:,2)==n); [l1,m1] = min( len(n1) ); 
%         n2 = find(edge(:,1)==n); [l2,m2] = min( len(n2) );  
%         cRct(n) = interp( [-l1,l2], cRct(edge(n1,1),edge(n2,2)), 0 ); 
%     end 
    cAdv([11,4]) = interp1( [0,1], cAdv([15,1]), len([1,2])/sum(len(1:3)) ); 
    cAdv([3]) = interp1( [0,1], cAdv([1,6]), len([6])/sum(len([6,19])) );
    cAdv([7]) = interp1( [0,1], cAdv([13,16]), len([33])/sum(len([33,35])) );
    cAdv([2]) = interp1( [0,1], cAdv([3,12]), len([12])/sum(len([12,9])) );
    cAdv([5]) = interp1( [0,1], cAdv([17,7]), len([6])/sum(len([6,19])) );
    cAdv([9]) = interp1( [0,1], cAdv([3,19]), len([21])/sum(len([21,36])) );
    
elseif( nG == 1 )
        cellind2state = [ 2 4 4 4 3 1 2 3 2 ]; 
        cProlif = cProlif( cellind2state ); 
        cAdvec = cAdvec( cellind2state ); 
        cRct = cProlif; 
        cAdv = cAdvec; 
end 

for n = 1:nEdge 
    cRctM(:,n) = linspace( cRct(edge(n,1)),  cRct(edge(n,2)), N1 );     
    cAdvM2(:,n) = linspace( cAdv(edge(n,1)),  cAdv(edge(n,2)), N1 );     
end 

cAdvM2 = 2*(1-cAdvM2).* cRctM; 
for n = 1:nEdge 
% %     cAdvM2(:,n) = cAdvM2(:,n) / norm( vertex(edge(n,1),1:2) - vertex(edge(n,2),1:2) ); 
%     cAdvM2(:,n) = cAdvM2(:,n) /len(n) / 10; 
    cAdvM2(:,n) = cAdvM2(:,n) /len(n); 
end 

for n = BC.noout
    ind = find(edge(:,2)==n); 
    cAdvM2(:,ind) = cAdvM2(:,ind) .* ( (1-xx.^2)'*ones(1,length(ind)) ); 
end 
% cAdvM2(:,[2,3,5,7,9,11]) = cAdvM2(:,[2,3,5,7,9,11]) .* ( (1-xx.^2)'*ones(1,6) ); 

cAdvM2(:,edge(:,3)==0) = 0; 


cDthM = 1./steady; 

%%%% signal - depending on the final number 
% global cK cS barv 
% % cK = 50; %12.8 * 0.1^2; 
% % disp( cK )
% cK = 20; 
% cS = 1; %1-1/(1+exp( -(sum(Un(3:end)) - barv)*cK )); 
% barv = 0.35; 

BC.Uinit(1) = steady(1,1); 

%% what is this??? 
if( nG == 0 ) 
    BC.Nnmn = [1]; 
elseif( nG == 21 ) 
    BC.Nnmn = [1]; 
elseif( nG == 1 ) 
    BC.Nnmn = [3,7,12,13];     
end
    
end 

function ntmp = Compute_nTotC( U ) 
global wUn N1 nEdge int_dx 
ntmp = 0; 
for n = 1:nEdge 
    wtmp = linspace( wUn(1,n), wUn(2,n), N1 )'; 
    ntmp = ntmp + int_dx'*(U(:,n).*wtmp);
end    
end 

function UU = plot_graph1d_half( U, nfig ) 
global N1 edge nEdge int_dx nVrtx wUn 
    Nf = (N1+1)/2; 
    
    UU = zeros( nVrtx, 1 ); 
    for n = 1:nEdge 
        wtmp = linspace( wUn(1,n), wUn(2,n), N1 )'; 
        UU(edge(n,1)) = UU(edge(n,1)) + sum(U(1:Nf, n).*int_dx(1:Nf).*wtmp(1:Nf)) ; 
        UU(edge(n,2)) = UU(edge(n,2)) + sum(U(Nf:end, n).*int_dx(Nf:end).*wtmp(Nf:end)) ; 
    end    
    
    if( nargin > 1 )
      figure(nfig); hold on;  plot( UU, '-x', 'linewidth', 1 )
    end     
end 

