function dnC = Compute_du_2Dsp( nC, time )
global Dxn Dxxn cDiff V Dxm Dxxm cRct cDth cAdv V2 Vaml cS addcRct

React = cRct .* nC ;
Advec = cS * 2*( 1-cAdv).* React;

React = ( 1 - cDth.*nC ); React(React<-0.6925)=-0.6925;
React = cRct .*React.*nC;

dnC = cDiff *(Dxxn *nC + nC* Dxxm') ...
        - Dxn * ( V{1}.* nC + V2{1}.*Advec) - ( V{2}.* nC + V2{2}.*Advec)*Dxm' + React;

% % React = addcRct.*nC;
% dnC = dnC - Dxn * ( Vaml{1}.*React) - ( Vaml{2}.*React) * Dxm' ; % + React;
% % dnC = dnC + React;


end
