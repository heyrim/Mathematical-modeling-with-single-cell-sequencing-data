function plot_graph3d( U, nfig, nfigsub ) 
global x N1 vertex edge nEdge  
    if( nargin == 1 );          figure; hold on; 
    elseif( nargin == 2 );      figure(nfig); hold on; 
    else; figure(nfig);  subplot( nfigsub(1), nfigsub(2), nfigsub(3) ); hold on; 
    end 
    
    set(gca, 'color', 'black' ); 
    for n = 1:nEdge 
      edgex{n} = linspace( vertex(edge(n,1),1), vertex(edge(n,2),1), N1 );
      edgey{n} = linspace( vertex(edge(n,1),2), vertex(edge(n,2),2), N1 );
    end 
    
    for m = 1:nEdge  
      surface([edgex{m};edgex{m}],[edgey{m};edgey{m}],[U(:,m)';U(:,m)'],[U(:,m)';U(:,m)'],'facecol','no','edgecol','interp','linew',2);  
    end 
    
    plot3( vertex(:,1), vertex(:,2), min(U(:))*ones(length(vertex),1), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5] ) 
    set(gca,'XTick', [], 'XTicklabel', [] )
    set(gca,'YTick', [], 'YTicklabel', [] )

    view( [0 90] );  box on;  %grid on; 
 
    set(gcf,'Renderer','opengl');     
    
end 