%% PDF on space 
function nCell = PDFonContinuum_2019( uinit, steady, T, nfig, ntest )
global N1 N Dxn Dxxn Dxm Dxxm cDiff cDifff cRct cAdv cDth cA uS xx yy dxn dxm center dx
global cS cK V barv

%% ntest == 0 Nestorowa / 1 Paul

% prestep
if (isempty(steady))
    steady = get_usteady(ntest);
end
if( nargin < 3 ); nfig = round( rand(1)*1000 ); end
if( nargin < 4 ); T = 100;  end



%% Prolif etc
[cRct, cDth, cAdv, cA, indSig] = Compute_system( steady, ntest );
uS = steady; uS( uS<=eps ) = eps;

set_system;
% cDiff=0.1^4; cDifff =0.1^4;
% cK = 50;
barv = 0.05;


%% Potential
[V] = Compute_V( steady, ntest );
Vmax = max( max(max(V{1})), max(max(V{2})) );

%% IC
dx = [dxn dxm];
nCell = IC( uinit, 0.1, steady*0.05, ntest );
% nCell = IC( uinit, 0, steady, ntest ); %% Start from steady state
% nCell = steady/(sum(sum(steady))*prod(dx));
nCell = BC( nCell );

nTotCell(1) = sum(sum(nCell))*prod(dx); %dxn*dxm;
nCellSave(:,1) = Compute_cluster( nCell, ntest );

cS = 1-1/(1+exp( -(sum(nCellSave(indSig,1)) - barv)*cK ));
cSsave(1) = cS;

%% AML at t=0
global nAML 
if( nAML ) 
    [indSig] = Update_AML( nCell, ntest ) ;

    cS = 1-1/(1+exp( -(sum(nCellSave(indSig,1)) - barv)*cK ));
    cSsave(1) = cS;
end 


disp( strcat( 'params : ', num2str([cDiff cDifff cS cK barv]) ) );

figure(nfig); hold on; if(  N(2) == 1 ); plot( xx, nCell );
else; subplot( 3, 5, 1); imagesco( nCell' ); title( strcat( 't=', num2str(0) ) );end

%% Time integration
Tstep = 1; dt = 0.001; t = 0;
Tstart = tic;

if( max(dx)/max(max(cAdv))/2 < dt ||  max(dx)/max(Vmax)/2 < dt || prod(dx)/max(cDiff)/2 < dt )
    disp( 'reduce dt' );
end

for nT = 1:(T/Tstep)
    for nt = 1:(Tstep/dt)
        %       t = 1;
        
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
        
        %      [V, cDifff] = Update_System( V, cDifff, nCell );
        
        nCtmp = Compute_cluster( nCell, ntest );
        %       cS = 1-1/(1+exp( -(sum(nCtmp([3,4,5,6,8])) - barv)*cK ));
        cS = 1/(1+exp( (sum(nCtmp(indSig)) - barv)*cK ));
        
    end % end of time inner loop
    
    nstep = 5;
    if( mod(Tstep*nT, nstep ) == 0 )
        figure(nfig); hold on;
        if(  N(2) == 1 ); plot( xx, nCell );
        else; subplot( 3, 5, min(Tstep*nT/nstep + 1, 15)); imagesco( nCell' );
            title( strcat( 't=', num2str(nT*Tstep) ) ); end
    end
%       if( Tstep*nT == 20 )
%         disp('stop') 
%       end
    
    nTotCell(nT+1) = sum(sum(nCell))*prod(dx); %dxn*dxm;
    nCellSave(:,nT+1) = Compute_cluster( nCell, ntest );
    %   cS = 1-1/(1+exp( -(sum(nCellSave([3,4,5,6,8],nT+1)) - barv)*cK ));
    
    cSsave(nT+1) = cS;
    
    %   if( Tstep*nT == 10 )
    %     Update_AML( nCell, ntest ) ;
    %
    %   end
    
end

toc( Tstart );


%% plot results
figure(nfig+4); hold on; plot( [0:Tstep:T], nTotCell, 'black' ) 
load( 'colororder' ) 
figure(nfig+3); hold on; plot( [0:Tstep:T], nCellSave'  ) 
% %   for n = 1:8; plot( [0:Tstep:T], nCellSave(n,:), '-', 'color', co(n,:)  ); end
% figure(nfig+5); hold on; plot(  [0:Tstep:T], cSsave  )

end

function dnC = Compute_du( nC, time )
global Dxn Dxxn cDiff V Dxm Dxxm cRct cDth cAdv cA Vaml uS cS cDifff dx addcRct

React = cRct .* nC ;
Advec = cS * 2*( 1-cAdv).* React;
%     Advec = 2*( 1-cAdv*cS ).* React;
%     Advec = 2*( 1-cAdv ).* React;

%     React = (1./( 1 + exp( (nC - uS*1)) )).* React;

React = ( 1 - cDth.*nC ); React(React<-0.6925)=-0.6925;
React = cRct .*React.*nC;

dnC = cDiff *(Dxxn *nC + nC* Dxxm') ...
        - Dxn * ( V{1}.* nC + cA{1}.*Advec) - ( V{2}.* nC + cA{2}.*Advec)*Dxm' + React;
%     - Dxn * ( V{1}.* nC ) - ( V{2}.* nC )*Dxm' + React;


% %     dnC = dnC + (1./( 1 + exp( (nC - uS*1)) )).* React;
%     ntmp =  (nC - uS); ntmp(ntmp>0)=0; ntmp(ntmp<0)=1;
%     dnC = dnC + ntmp.* React;
% % %     dnC = dnC + (cRct - (1-cS)*cDth .*nC) .* nC;
% % % %     dnC = dnC + (1 - sum(sum(nC))*prod(dx))* React;

React = addcRct.*nC;
dnC = dnC - Dxn * ( Vaml{1}.*React) - ( Vaml{2}.*React) * Dxm' ; % + React;
% dnC = dnC + React;

%     dnC = dnC - Dxn * ( Vaml{1}.*nC) - ( Vaml{2}.*nC) * Dxm' ;
%     dnC = dnC + addcRct.*nC;

end

function nC = BC( nC )

% [n,m] = size( nC );
nC( [1,end], :) = 0;
% nC(1,:) = 0;
% nC(end,:) = (4*nC(end-1,:)-nC(end-2,:))/3;
if( size( nC,2 ) > 1 )
    nC( :, [1,end]) = 0;
    %     nC(:,end) = (4*nC(:,end-1)-nC(:,end-2))/3;
    % %     nC(:,1) = 0;
end
% nC(1,:) = 0;

end

function [V] = Compute_V( usol, ntest )
global cDiff N N1 xx yy cRct dx

% method 1 - with min value
%  usol(usol==0) =  min(min( usol(usol~=0) ));

ulog = -log( usol );

Dxn = Compute_Dx2nd( N(1)+1, dx(1) );

if( N(2) == 1 )
    %% 1D
    V = -cDiff * Dxn*ulog;
    
    %%%%%% including reaction does not work. SUPER UNSTABLE
    %     int_dx = ones( N1(1), 1 )*dxn; int_dx(1) = dxn/2;
    %     for n = N1(1):-1:2
    %         tmp = ulog(n) - ulog(1:n);
    %         int_dx(n) = dxn/2;
    %         tmptmp(n) = (exp(tmp).*cRct(1:n))'*int_dx(1:n);
    %     end
    %     tmptmp(1) = 0;
    %     V = V + exp(ulog).*tmptmp;
    
    % %check
    % tmp = -V .* usol + cDiff * Dxm * usol;
    
else
    
    % method 2
    [XX,YY] = meshgrid( 1:N1(2), 1:N1(1) );
    ind = find( ulog ~= inf ) ;
    
    % bdind = []; for n = 2:(N1(2)-1); bdind = [bdind; sub2ind( N1, 1,n ); sub2ind( N1, N1(1),n )]; end
    %  for n = 1:N1(1); bdind = [bdind; sub2ind( N1, n,1 ); sub2ind( N1, n,N1(2) )]; end
    % bdind = setdiff( bdind, ind );
    
    ind = find( usol > max( usol(:) )*0.0001 );
    F = scatteredInterpolant(XX(ind), YY(ind), ulog(ind), 'natural', 'none');
    ulog = F(XX,YY); 
    if( ntest == 0 ); ulogMax = 3.5; elseif( ntest == 1 ); ulogMax = 5.5; end 
    ulog(isnan(ulog)) = ulogMax; 
    ulog( ulog>ulogMax ) = ulogMax; 
    
    Dxm = Compute_Dx2nd( N(2)+1, dx(2) );
    
    % dx = 0.04;
    % Dxn = Dxn *dxn/dx; %(xx(2)-xx(1));
    % Dxm = Dxm *dxm/dx; %(yy(2)-yy(1));
    
    V{1} = -cDiff * Dxn*ulog;
    V{2} = -cDiff * ulog*Dxm';
    
end

end


function [cRct, cDth, cAdv, cAvec, indSig] = Compute_system( usol, ntest )
global Vaml addcRct N N1

if( ntest == 0 )
    [cRct, cAdv, cA] = PrePocess_Nestorowa();
    indSig = [3:6,8];
elseif( ntest == 1 )
    [cRct, cAdv, cA] = PrePocess_Paul();
    indSig = [2,3,4,8,9];
end

usol(usol<=eps)=eps;
cDth = 1./usol;

%     nn = max(max(usol));
%     cDth = cRct/nn;
%     cRct = cRct.*(usol/nn);

%     cAdv = 2*( 1-cAdv ).* cRct;

cA(isnan(cA)) = 0;
cAvec{1} = squeeze( cA(:,:,1) );
cAvec{2} = squeeze( cA(:,:,2) );

%%%%% normalization from 6 (Anna) to 0.7
%     cAvec{1} = cAvec{1}*6/0.7;
%     cAvec{2} = cAvec{2}*6/0.7;

Vaml = cAvec;
for n = 1:2; Vaml{n}(:) = 0; Vaml{n} = sparse( Vaml{n} ); end

addcRct = zeros( N1(1), N1(2) ); addcRct(:) = 0; addcRct = sparse( addcRct );

end

function nC = Compute_cluster( nCell, ntest )
global Vw dx

if( ntest == 0 )
    % load( '190226_Nestorowa_weights.mat' ) % for 101x101
    load( 'data_Matlab/190227_Nestorowa_weights.mat' )
elseif( ntest == 1 )
    load( 'data_Matlab/190517_Paul_weight.mat' ); % for 101x101
    % change to 121
    W(:, 21:121, : ) = Vw(:,:,:);
    % mean( cc(1, 6:7) )  = 0.2692
    % mean( cc(1,[7,5]) ) = 0.4340
    % mean( cc(1,[5,4]) ) = 0.5513
    W( 1:26, 1:20, 6) = 1;   W(27:43, 1:20, 7 ) = 1;  W(44:55, 1:20, 5 ) = 1;
    W(56:end,1:20, 4) = 1/2; W(56:end,1:20, 8 ) = 1/2;
    
    Vw = W;
end

nVrtx = size( Vw, 3 );
nC = zeros( nVrtx, 1 );
for n = 1:nVrtx
    nC(n) = sum(sum(nCell .* squeeze( Vw(:,:,n) )))*prod(dx); %dxn*dxm;
end


end

function indSig = Update_AML( U, ntest )
global xx yy int_dx cRct cAdv cDth Vaml V addcRct edge N1 cK barv caml

if( ntest == 0 )
    load( 'data_Matlab/181019_Data_Nestorowa.mat' )
    
    [X,Y] = meshgrid([0:.01:1],[0:.01:1]);
    [XX,YY] = meshgrid(xx,yy);
    MEP = interp2( X,Y, VV{4}, XX, YY )';
    STHSC = interp2( X,Y, VV{1}, XX, YY )';
    
    %%%%%%% signal param
    % barv = barv*10;
    % cK = 10;
    
    %%%%%%% psuedo inverse
    load( 'data_Matlab/190406_Nstrw_AML_vec.mat' );
    
    for n = 1:2; Vaml{n} = newvec(:,:,n); end
    
    %%%%%%% potential
%     % % global cDiff
%     AMLcenter = [0.6101   0.2153]; % shiifted from [1 -2]    
%     usol = Normal2D(AMLcenter, 0.04*[1,0.3;0.3,1], xx', yy' );
%     
%     Vaml = Compute_V( usol, 0 ) ;
%     
%     usol(usol>10) = 10;
%     addcRct = usol/max(usol(:)) * 0.6931;
    
    
%      STHSC(STHSC>20)=20; STHSC = STHSC/20 * (1-1/2.5);
%      cDth = cDth.*(1 - STHSC);
    STHSC(MEP>20)=20; STHSC = STHSC/20 * 4 + 1;
    cDth = cDth./STHSC;
    
    
elseif( ntest == 1 )
    
    load( 'data_Matlab/190604_Data_Paul.mat' );
    MEP(:, 21:121) = WW;
    %  STHSC(:, 21:121) = VV{6};
    
    %  load( 'data_Matlab/190604_Paul_AML_vec.mat' );
    %  load( 'data_Matlab/190612_Paul_AML_vec_straight.mat' );
    load( 'data_Matlab/190614_Paul_AML_vec.mat' );
    for n = 1:2; Vaml{n}(:,21:121) = newvec(:,:,n); end
    
    %%%%%%% potential
%     AMLcenter = [0.6  0.2];
%     usol = Normal2D(AMLcenter, 0.04*[1,0;0,0.2], xx', yy' );
%     %%%%%% vector from potential
%     Vaml = Compute_V( usol, 0 ) ;
%     
%     
%     usol(usol>10) = 10;
%     addcRct = usol/max(usol(:)) * 0.6931;
    
    
end

%  MEP(MEP>20)=20; MEP = MEP/20 * (1-0.1);
%  cDth = cDth.*(1 - MEP);
MEP(MEP>20)=20; MEP = MEP/20 * 19 + 1;
cDth = cDth./MEP;



for n = 1:2; Vaml{n} = caml *Vaml{n}; end

for n = 1:2; Vaml{n} = Vaml{n} - V{n}; end
% % ind = find( usol >= 10 );
% % for n = 1:2; Vaml{n}(ind) = Vaml{n}(ind) - V{n}(ind); end


if( ntest == 0 )
    indSig = [3,5,6,8];
elseif( ntest == 1 )
    indSig = [2,3,9];
end

end

function nCell = IC( uinit, ninit, usol, ntest )
global xx yy dxn dxm N1 N

if( ntest == 0 )
    if( isempty(uinit))
        load( 'data_Matlab/181019_Data_Nestorowa.mat' )
        
        [X,Y] = meshgrid([0:.01:1],[0:.01:1]);
        [XX,YY] = meshgrid(xx,yy);
        uinit = interp2( X,Y, VV{2}, XX, YY )';
        
    end
    
elseif( ntest == 1 )
    load( 'data_Matlab/190419_PaulData.mat' ); %dc(:,2) = (dc(:,2)+0.2)*(5/6);
    y = [-0.2:1.2/N(2):1]';
    [YY,XX] = meshgrid(y,xx);
    cluster(cluster==10) = 6; ii = find( cluster == 6 );
    [uinit,~,~] = ksdensity( [dc(ii,:)], [XX(:),YY(:)], 'bandwidth', [0.03 0.03] );
    
    %     ii = setdiff( ii, find( dc(:,1) > 0.2 ) ); ii = setdiff( ii, find( dc(:,2) > 0.2 ) );
    %    [uinit,~,~] = ksdensity( [dc(ii,:)], [XX(:),YY(:)], 'bandwidth', [0.02 0.02] );
    
    uinit = reshape( uinit, N1(1), N1(2) );
    
    %     load( 'data_Matlab/190517_PaulInit.mat' );
end

%     nCell = usol*ninit;
if( nargin == 1 )
    ninit = 1;
end

nCell = uinit/ (sum(sum(uinit))*prod([dxn dxm])) *ninit ;

if( ~isempty(usol) )
    nCell = nCell + usol;
end


end

function usol = get_usteady( ntest )
global N N1 dx xx yy
if( ntest == 0 )
    load( 'data_Matlab/181019_Data_Nestorowa.mat' )
    N = [100 125];
    cellind( cellind == 9 ) = 2;
    dc = dc_normalize(dc, 2);
    y = [0:1/N(2):1]';
elseif( ntest == 1 )
    load( 'data_Matlab/190419_PaulData.mat' );
    N = [100 120];
    y = [-0.2:1.2/N(2):1]';
end

xx = [0:1/N(1):1]';  yy = [0:1/N(2):1]'; dx = [1/N(1), 1/N(2)];
[YY,XX] = meshgrid(y,xx);


[usol,~,~] = ksdensity( [dc(:,1:2)], [XX(:),YY(:)], 'bandwidth', [0.03 0.03] );

N1 = N+1;
usol = reshape( usol, N1(1), N1(2) );


end

function Dx = Compute_Dx2nd( N1, dx )

Dx = (zeros(N1));   %%% 1D Dx
Dx(1,1:3) = [-3,4,-1]./(2*dx);
for n = 2:(N1-1)
    Dx(n, (n-1):(n+1)) = [-1,0,1]./(2*dx);
end
Dx(N1,(N1-2):N1) = [1,-4,3]./(2*dx);
Dx = sparse( Dx );
end

