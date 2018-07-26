function [PathwayQuant, UniProt2Reactome] = DataProcessPathways(ProteinOutput, species, limit, uniprotfile, isPTMfile, ctrl, grouping)

%wait = waitbar(0,'Getting Reactome data...', 'CreateCancelBtn',...
 %   'setappdata(gcbf,''canceling'',1)');
%setappdata(wait,'canceling',0)

%Download pathway data from Reactome.
fprintf('Getting Reactome data...')
url = 'http://reactome.org/download/current/UniProt2Reactome.txt';

filename = 'UniProt2Reactome.tab';
try outfilename = websave(filename,url);
    UniProt2Reactome = readtable(outfilename,'FileType','text','Delimiter','tab');
    UniProt2Reactome = table2cell(UniProt2Reactome);
catch
    outfilename = websave(filename,url);
    UniProt2Reactome = readtable(outfilename,'FileType','text','Delimiter','tab');
    UniProt2Reactome = table2cell(UniProt2Reactome);
end

%Convert species name to full name
switch species
    case 'human'
        species = 'Homo sapiens';
        if grouping; genecol = 12;
        else genecol = 2;
        end
    case 'mouse'
        species = 'Mus musculus';
        if grouping; genecol = 13;
        else genecol = 2;
        end
end

if isPTMfile
    uniprotidcol = strcmp(ProteinOutput(1,:),'UNIPROT ID');
    FCbegin = 6;
    ptm = 3;
    ptm2 = 4;
else
    uniprotidcol = strcmp(ProteinOutput(2,:),'UNIPROT ID');
    FCbegin = 3;
    ptm = 0;
    ptm2 = 0;
end
uniprotidcol = find(uniprotidcol);
UniProt2Reactome = UniProt2Reactome(strcmp(UniProt2Reactome(:,6),species),:);
fprintf('Got Reactome pathway data!')

%Annotate pathway file with Identifiers used in protein quant
for i = 1:size(uniprotfile,1)
    temp = strcmp(UniProt2Reactome(:,1),uniprotfile(i,1));
    if sum(temp) > 0
        temp2 = uniprotfile(i, genecol);
        UniProt2Reactome(temp,7) = temp2(1,1);
        fprintf(['Pathways found for ', char(uniprotfile(i,1)),'\n'])
    end
    %waitbar(i/size(uniprotfile,1), wait,...
     %       'Comparing Reactome database with ProteinQuant...');
end

%Get pathways represented in ProteinOutput and assign proteins in 
%ProteinOutput to pathways
x1 = size(UniProt2Reactome,1);
x2 = size(ProteinOutput,1);
PathwayCount = zeros(x1,1);
ProteinsInPathway = zeros(x2,x1);
for i = 3:x2
    temp = strcmp(UniProt2Reactome(:,7),ProteinOutput(i,1)); %Which pathways protein i appears in?
    PathwayCount = PathwayCount + temp;                      %How many proteins in each pathway?
    ProteinsInPathway(i,:) = ProteinsInPathway(i,:) + temp'; %Record which protein is in each pathway.
    fprintf(['Protein #',num2str(i),' ', char(ProteinOutput(i,1)),' found in ', num2str(sum(temp)), ' pathway(s).\n'])
    %waitbar(i/x2, wait, 'Finding proteins in each pathway...');
end
UniProt2Reactome = UniProt2Reactome(PathwayCount >= 1,:);
ProteinsInPathway = ProteinsInPathway(:,PathwayCount >= 1);

%Collapse multiple entries for different pathways
[uniquePathways,ic,ia] = unique(UniProt2Reactome(:,2),'stable');
UniProt2Reactome2 = UniProt2Reactome(ic,:);
ProteinsInPathway2 = zeros(x2,size(uniquePathways,1));
totalPinP = zeros(size(uniquePathways,1),1);
for i = 1:size(uniquePathways,1)
    %temp = strcmp(UniProt2Reactome(:,2),uniquePathways(i));
    totalPinP(i) = numel(UniProt2Reactome(ia==i,1));
    %temp2 = UniProt2Reactome(temp,2:4);
    ProteinsInPathway2(:,i) = sum(ProteinsInPathway(:,ia==i),2);
    %UniProt2Reactome2(i,:) = temp2(1,:);
end
UniProt2Reactome2 = UniProt2Reactome2(totalPinP >= limit,:);
x1 = size(UniProt2Reactome2,1);
ProteinsInPathway2 = ProteinsInPathway2(:,totalPinP' >= limit);

GroupNum = numel(unique(ProteinOutput(2,3:uniprotidcol-1)))-ptm;
SEbegin = FCbegin + GroupNum;
if isPTMfile, SEbegin = SEbegin + 4; end
FCend = SEbegin - 1;
if isPTMfile, FCend = FCend -4; end
SEend = SEbegin + GroupNum - 1;
groups = unique(ProteinOutput(2,FCbegin:FCend),'stable');

if isPTMfile
    PathwayQuant = cell(x1+2,2+size(ProteinOutput,2));
    PathwayQuant(1:2,:) = [ProteinOutput(1:2,:),cell(2,2)]; %Headers
else
    PathwayQuant = cell(x1+2,size(ProteinOutput,2));
    PathwayQuant(1:2,:) = ProteinOutput(1:2,:); %Headers
end
%delete(wait);
PathwayQuant(2,end) = {'NUMBER OF FEATURES'};
PathwayQuant(2,end-1) = {'DEGREES OF FREEDOM'};
ii =3;
for i = 3:x1+2
    %Get list of proteins in pathway i
    ProtList = ProteinOutput(~~ProteinsInPathway2(:,i-2),1:SEend);
    
    %Create long vector table
    switch isPTMfile
        case true;
            proteins = repmat(ProtList(:,1),GroupNum,1);
            sites = repmat(ProtList(:,3),GroupNum,1);
            for iii = 1:size(ProtList,1)
                proteins(iii) = {[char(proteins{iii}),num2str(sites{iii})]};
            end
        otherwise;
            proteins = repmat(ProtList(:,1),GroupNum,1);
    end
    nprot = size(ProtList,1);
    if nprot < limit, continue; end
    
    abundances = cell2mat(ProtList(:,FCbegin:FCend));
    %abundances = bsxfun(@rdivide,abundances,max(abs(abundances),[],2));
    abundances = abundances(:);
    Groups = repmat(groups, nprot,1);
    Groups = Groups(:);
    SEs = cell2mat(ProtList(:,SEbegin:SEend));
    SEs = min(1,abs(abundances./SEs(:)));
    lmetable = table(proteins, abundances, Groups, SEs,...
        'VariableNames', {'Proteins', 'Abundances', 'Groups', 'SE'});
    lmetable.Proteins = categorical(lmetable.Proteins);
    lmetable.Groups = categorical(lmetable.Groups);
    
    %Create design matrix and response variable
    X = [dummyvar(lmetable.Groups), dummyvar(lmetable.Proteins), ones(size(abundances,1),1)];
    Xids = [repmat({'Group'},1,GroupNum),repmat({'Feature'},1,nprot),{'Intercept'}];
    
    %weight abundances by SEs
    Y = lmetable.Abundances;
    
    %Fit model
    mdl = weighted_bayeslm(X,Y,Xids,true,lmetable.SE,[]);
    
    G = strcmp(mdl.FeatureType,'Group');
    %I = strcmp(mdl.FeatureType,'Intercept');
    %results = zeros(1,sum(G,2)+1);
    results = mdl.beta_estimate(:,G);% - mdl.beta_estimate(:,I);
    %SE = zeros(1,GroupNum);
    SE = mdl.SEPred(:,G);
    %SE(1,I) = mdl.SEPred(:,I);
    P = ones(1,GroupNum);
    P(1,G) = mdl.Pvalues(:,G);
    PathwayQuant(ii,FCbegin:FCend) = num2cell(results);
    PathwayQuant(ii,SEbegin:SEend) = num2cell(SE);
    PathwayQuant(ii,SEend + 1:SEend + GroupNum) = num2cell(P);
    PathwayQuant(ii,end) = num2cell(nprot);
    PathwayQuant(ii,end-1) = num2cell(mdl.DoF);
    PathwayQuant(ii,end-2) = {unique(lmetable.Proteins)};
    PathwayQuant(ii,1) = UniProt2Reactome2(i-2,2);
    PathwayQuant(ii,2) = UniProt2Reactome2(i-2,4);
    fprintf(['#',num2str(ii),' ',char(UniProt2Reactome2(i-2,2)),' ',num2str(nprot),' ',num2str(results(:)',2),'\n'])
    ii = ii + 1;
    %waitbar(i/(x1+2), wait,'Fitting Pathway models...');
end
PathwayQuant = PathwayQuant(1:ii-1,:);
%% Do Contrasts and Empirical Bayes correction of SEMs and p vals
%Performed as detailed in Kammers et al. (2015) Detecting significant 
%changes in protein abundance. EuPA Open Proteomics 7, 11:19.
ctrl_betas = cell2mat(PathwayQuant(3:end,FCbegin + ctrl-1));
DoFs = cell2mat(PathwayQuant(3:end,end-1));
ds = zeros(GroupNum,2);
%ctrl_intbetas = cell2mat(InteractionTable(4:end,RA-1 + find(I == 1)));
for i = 0:1:GroupNum-1

    d0s0 = double.empty(0,2);
    iter = 1;
    SEs = cell2mat(PathwayQuant(3:end, SEbegin + i));
    %SEuseable = SEs < 100;
    %SEs = SEs(SEuseable,:);
    
    %Do constrasts for all groups
    PathwayQuant(3:end,FCbegin + i) = ...
        num2cell(cell2mat(PathwayQuant(3:end,FCbegin + i)) - ctrl_betas);
    
    while isempty(d0s0)
       try d0s0 = mle(SEs,'pdf',@(x,v,s)chi2pdf(x/s,v)/s,'start', [rand,rand]);
       catch
           %try d0s0 = mle(cell2mat(ProteinAbundance(3:end,3+GroupNum + i)),...
           %    'pdf',@(x,v,s)chi2pdf(x/s,v)/s,'start', [1,0.01]);
           %catch
           if iter > 10000
               warning(cat(2,'SEM estimates whacky - cannot get hyperparameters! No Empirical Bayes on column ', num2str(i),'!'));
               d0s0 = [0, 1];
               break;
           end
           iter = iter + 1;
           continue; 
       end
    end
   
    %Do correction for each protein.
    for ii = 3:size(PathwayQuant,1)
        DoF = DoFs(ii-2);
        lambda = DoF/(DoF + d0s0(1));
        oldSEM = cell2mat(PathwayQuant(ii,SEbegin + i));
        oldVar = (oldSEM/sqrt(2/((DoF + d0s0(1))/2 + 2))).^2;
        newSEM = sqrt(lambda * oldVar + (1 - lambda) * d0s0(2))...
            * sqrt(2/((DoF + d0s0(1))/2 + 2));
        PathwayQuant(ii,SEbegin + i) = {newSEM};
        FC = abs(cell2mat(PathwayQuant(ii,FCbegin+i)));
        t = FC/newSEM;
        p = min(1,(1 - tcdf(t, DoF + d0s0(1) - 1)) * 2 + 1e-15); %2-sided t-test
        PathwayQuant(ii,SEend + i + 1) = {p};
        %PathwayQuant(ii,end-1) = {DoF + d0s0(1) - 1}; %New degrees of freedom
    end
    ds(i+1,:) = d0s0;
end

%% Final step: Adjust p-values according to Benjamini & Hochberg  1995
%(doc mafdr for more info)
for i = FCbegin + (GroupNum * 2)+ptm2:FCbegin - 1 + (GroupNum * 3)+ptm2
    %if i ~= 2 + (GroupNum * 2) + find(I == 1)
    try temp_pvals = bhfdr(cell2mat(PathwayQuant(3:end,i)));%,...
           % 'BHFDR', true); %Call to modified function to avoid licence check
           if isPTMfile
               PathwayQuant(3:end,i + GroupNum) = num2cell(temp_pvals);
           else
               PathwayQuant(3:end,i + GroupNum + FCbegin) = num2cell(temp_pvals);
           end
    catch
        warning('mafdr() failed, probably a licensing thing.');
    end
end
if isPTMfile
    PathwayQuant(1,11+GroupNum) = {'SE'};
    PathwayQuant(1,11+GroupNum*2) = {'p-values'};
    PathwayQuant(1,11+GroupNum*3) = {'BH-FDR'};
    PathwayQuant = PathwayQuant(:,[1:2,6:6+GroupNum-1,6+GroupNum+4:end]);
end
end