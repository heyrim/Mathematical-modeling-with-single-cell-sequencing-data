function [vec, mappedX, newmapped] = Compute_Gene_Perturbation( ntest, DiffMapConst  )  

if( nargin < 1 ) 
    ntest = 0; 
    DiffMapConst = [10, 0.5, Inf]; 
elseif( nargin < 2 ) 
    global dmParams 
    DiffMapConst = [dmParams(1), 0.5, dmParams(2)]; 
end 


if( ntest == 0 ) 
    load( 'data/Nestorowa2016_scRNAseqData' )
elseif( ntest == 1 ) 
    load( 'data/Paul2015_scRNAseqData.mat' ); 
end 



%% read AML related gene alternation 
[numGene, scale] = AMLinfo( genename ); 
% [numGene, scale] = AMLinfo_YaHuei(genename) ; 

%% modify gene expression 
addscale = 0; % to amplify gene alternations, 0 for none, 4 for 2^4 times, etc. 
maxdata = min( max(scdata(:)), 16 ); % if addscale == 100; 

newdata = scdata; 
for nn = 1:length(numGene) 
    if( addscale < 100 ) 
        newdata( numGene(nn), : ) = newdata(numGene(nn),:) + scale(nn) + sign(scale(nn))* addscale; 
    else
      if( scale(nn)>0 )
        newdata( numGene(nn), : ) = maxdata; 
      else
        newdata( numGene(nn), : ) = 0;           
      end
    end
end 

clear ntmp GeneAlter 


%% Find projected direction in Diffusion Mapping space computed by raw data 
[mappedX,~,newmapped] = DiffusionMap_wNewData(scdata', newdata', DiffMapConst(1), DiffMapConst(2), DiffMapConst(3) ); 
mappedX = mappedX(:,1:2); 
newmapped = newmapped(:,1:2); 

% vector of alternation 
vec1 = newmapped(:,1)-mappedX(:,1); vec2 = newmapped(:,2)-mappedX(:,2); 
vec = [vec1,vec2]; 


% plot results 
% figure; subplot( 1, 2, 1 ); hold on; 
% plot( mappedX(:,1), mappedX(:,2), 'o', 'markersize', 3, 'markerfacecolor', [0    0.4470    0.7410] ); 
% plot( newmapped(:,1), newmapped(:,2), 'x', 'markersize', 3 ); 
% subplot( 1, 2, 2 ); quiver( mappedX(:,1), mappedX(:,2), vec1, vec2, 'black' )



end

%%%% mean of vector 0.0267 (nscale=0)
%%%%  0.0730 (nscale=4) / 0.1536 (16) / 0.1811 (max 16/min 0) 
 
function [numGene, scale] = AMLinfo( genename )

%% gene related with Leukemic Stemness (log_2 increment) 
% data from Nature 2016;540(7633):433--437. doi:10.1038/nature20598., 
% 'A 17-gene stemness score for rapid determination of risk in acute leukaemia'. 
GeneAlter = {'Cd34' 2.15
'Laptm4b' 1.8 
'Arhgap22' 1.48  
'Smim24' 1.45  
'Kiaa0125' 1.4
'Mmrn1' 1.36 
'Dnmt3b' 1.31 
'Socs2' 1.24 
'Cdk6' 1.23 
'Ngfrap1' 1.2 
'Cpxm1' 1.2
'Zbtb46' 1.19
'Dpysl3' 1.16
'Nynrin' 1.15 
'Gpr56' 1.06 % 'Adgrg1'
'Akr1c3' 1.06 
'Emp1' 1.01 };

%% Find genes that are in our single-cell dataset
for nn = 1:size( GeneAlter, 1 ) 
    n = 1; 
    while( ~strcmpi( genename{n}, GeneAlter{nn,1} ) && n<length(genename) ); n=n+1; end
    ntmp(nn,1) = n; 
    ntmp(nn,2) = GeneAlter{nn,2}; 
end 
ind = find( ntmp(:,1) < length(genename) ); 

numGene = ntmp(ind,1); 
scale = ntmp(ind,2); 

clear ntmp GeneAlter 


%% Gene related to AML (read from excel file) 
% data from Blood 2016;127(16):2018--2028. doi:10.1182/blood-2015-11-683649.
% 'GPR56 identifies primary human acute myeloid leukemia cells with high repopulating potential in vivo.'
[~,~,GeneAlter] = xlsread('AML_phenotype') ; 
GeneAlter{1,1} = 'Adgrg1'; %% Change Gpr56 to Adgrg1 -- identical gene 


%% Find genes that are in our single-cell dataset
for nn = 1:size( GeneAlter, 1 ) 
    n = 1; 
    while( ~strcmpi( genename{n}, GeneAlter{nn,1} ) && n<length(genename) ); n=n+1; end
    ntmp(nn,1) = n; 
    ntmp(nn,2) = log2(GeneAlter{nn,2}); 
end 
ind = find( ntmp(:,1) < length(genename) ); 

numGene = [numGene; ntmp(ind,1)]; 
scale   = [scale;  ntmp(ind,2)]; 


end 

function [numGene, scale] = AMLinfo_YaHuei(genename) 

% Upregulated CD16 
n=1; while( ~strcmpi( genename{n}, 'FCGR3' ) && n<length(genename) ); n=n+1; end
numGene = n; 
% Upregulated CD32 
n=1; while( ~strcmpi( genename{n}, 'FCGR2b' ) && n<length(genename) ); n=n+1; end
numGene = [numGene; n]; 

scale = [2;2]; 

end



