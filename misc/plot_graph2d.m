function plot_graph2d( U, nfig, nfigsub ) 
global N1 edge nEdge 
    if( nargin == 1 );          figure; 
    elseif( nargin == 2 );      figure(nfig); 
    else; figure(nfig);  subplot( nfigsub(1), nfigsub(2), nfigsub(3) ); 
    end 

    imagesc( U' )
    
    set(gca,'YTick', 1:nEdge, 'YTicklabel', edge(end:-1:1,1) ); 
    yyaxis right
    set(gca,'YTick', [0:(1/12):1], 'YTicklabel', edge(end:-1:1,2) );    
    yyaxis left 
    set(gca,'XTick', [1:((N1-1)/5):N1], 'XTicklabel', [0:.2:1] )
    
end