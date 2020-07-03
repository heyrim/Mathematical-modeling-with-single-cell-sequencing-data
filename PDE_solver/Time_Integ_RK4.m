function uSave = Time_Integ_RK4( Compute_du, BC, uinit, T, dt, Tstep ) 
global N cAdv dx cDiff V indSig cS cK barv 

Vmax = max( max(max(V{1})), max(max(V{2})) );
if( max(dx)/max(max(cAdv))/2 < dt ||  max(dx)/max(Vmax)/2 < dt || prod(dx)/max(cDiff)/2 < dt )
    disp( 'reduce dt' );
end

t = 0; 
nCell = uinit; 

nTotCell(1) = sum(sum(nCell))*prod(dx); %dxn*dxm;
nCellSave(:,1) = Compute_cluster( nCell );
uSave(:,:,1) = nCell; 

cSsave(1) = cS;


Tstart = tic; 
for nT = 1:(T/Tstep)
    for nt = 1:(Tstep/dt)
        
        %%%% RK4
        [dnCell] = Compute_du( nCell, t );
        dnCellRK = dnCell;
        nCellRK = nCell + ( dnCell )*dt/2;
        nCellRK = BC( nCellRK ); nCellRK(nCellRK<0) = 0;
        
        [dnCell] = Compute_du( nCellRK, t+dt*0.5 );
        dnCellRK = dnCellRK + dnCell*2;
        nCellRK = nCell + ( dnCell )*dt;
        nCellRK = BC( nCellRK ); nCellRK(nCellRK<0) = 0;
        
        [dnCell] = Compute_du( nCellRK, t+dt );
        dnCellRK = dnCellRK + dnCell*2;
        nCellRK = nCell + ( dnCell )*dt;
        nCellRK = BC( nCellRK ); nCellRK(nCellRK<0) = 0;
        
        [dnCell] = Compute_du( nCellRK, t+dt );
        dnCellRK = dnCellRK + dnCell;
        nCell = nCell + ( dnCellRK )*dt/6;
        nCell = BC( nCell ); nCell(nCell<0) = 0;
        
        t = t + dt;
        
        nCtmp = Compute_cluster( nCell ); 
        cS = 1/(1+exp( (sum(nCtmp(indSig)) - barv)*cK )); 
        
    end % end of time inner loop
    
    nTotCell(nT+1) = sum(sum(nCell))*prod(dx); %dxn*dxm;
    nCellSave(:,nT+1) = Compute_cluster( nCell );
    uSave(:,:,nT+1) = nCell; 
   
    cSsave(nT+1) = cS;
    
    %   if( Tstep*nT == 10 )
    %     Update_AML( nCell, ntest ) ;
    %
    %   end
    
end

toc( Tstart );

end 