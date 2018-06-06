%If individual donors are not used donor must be set to false. Remeber to
%order your groups such that your controls are first.
%
%If ProteinGrouping is set to true, will treat all proteins with the same
%name as a single entity using all available peptides, otherwise each one
%will be calculated separately.
%
%mods is a cell row vector of strings that specifies which mods have been
%searched for (e.g. {'Phospho', 'Oxidation'}).
%
%organism accepts either 'human' or 'mouse' and is used to pull proteome
%data from UniProt.

%nDB is number of databases used during Mascot search.
%%
function [ProteinAbundance, normedpepstonorm, InteractionTable, ModList,...
    PCA, uniprotall] = DataProcessProtein(wholepeptidelist, pepstonorm,...
    normmethod, pepmin, donorvector, ProteinGrouping, mods, organism,...
    score, regmethod, model, scorethreshold, nDB, incSubject, EB)

%% Step 0: setup variables
%setup progress bar
%wait = waitbar(0,'Getting UniProt data...', 'CreateCancelBtn',...
 %   'setappdata(gcbf,''canceling'',1)');
%setappdata(wait,'canceling',0)
fprintf('Getting UniProt data...')

%Download organism UniProt data.
switch organism
    case 'human'
        url = 'http://www.uniprot.org/uniprot/?sort=&desc=&compress=no&query=proteome:UP000005640&fil=&force=no&preview=true&format=tab&columns=id,entry%20name,protein%20names,genes,go,go(biological%20process),go(molecular%20function),go(cellular%20component),go-id,interactor,sequence,database(GeneID),reviewed';
        upcol = 12;
    case 'mouse'
        url = 'http://www.uniprot.org/uniprot/?sort=&desc=&compress=no&query=proteome:UP000000589&fil=&force=no&preview=true&format=tab&columns=id,entry%20name,protein%20names,genes,go,go(biological%20process),go(molecular%20function),go(cellular%20component),go-id,interactor,sequence,database(GeneID),database(MGI),reviewed';
        upcol = 13;
    otherwise
        error('organism name can only be human or mouse')
end
filename = 'uniprotall.tab';
try outfilename = websave(filename,url);
    uniprotall = readtable(outfilename,'FileType','text','Delimiter','tab');
    uniprotall = table2cell(uniprotall);
catch
    outfilename = websave(filename,url);
    uniprotall = readtable(outfilename,'FileType','text','Delimiter','tab');
    uniprotall = table2cell(uniprotall);
end

switch organism
    case 'human'
        if size(uniprotall,1) < 70000
            warning('UniProt file download failed, try again.')
            return;
        end
    case 'mouse'
        if size(uniprotall,1) < 50000
            warning('UniProt file download failed, try again.')
            return;
        end
end
uniprotall_reviewed = uniprotall(strcmp(uniprotall(:,end),'reviewed'),:);
fprintf('Got UniProt data!')

%Get column where Raw/Normalized abundance values start
if strcmp(normmethod, 'PG') == 1
    RAcell = strcmp(wholepeptidelist(1,:), 'Normalized abundance');
    PGendcell = strcmp(wholepeptidelist(1,:), 'Raw abundance');
    
    for i = 1:1:size(PGendcell,2)
        if PGendcell(1,i) == 1
            PGend = i;
            break;
        else PGend = size(PGendcell,2) + 1;
        end
    end
else RAcell = strcmp(wholepeptidelist(1,:), 'Raw abundance');
end

for i = 1:1:size(RAcell,2)
    if RAcell(1,i) == 1
        RA = i;
        break;
    end
end

if strcmp(normmethod, 'PG') == 1
    width = (size(wholepeptidelist(1,RA:PGend-1), 2));
else width = (size(wholepeptidelist(1,RA:end), 2));
end

%Setup parameters
if donorvector == false
    donors = RA:RA-1+width;
else
    donors = donorvector;
end

numDonors = numel(unique(donors));
peptidelistabds = wholepeptidelist;
normedpepstonorm = pepstonorm;
gene_col = 7;
pl_length = size(peptidelistabds,1);
norm_length = size(normedpepstonorm,1);

%Calculate Mascot score threshold. Take only peptides with "Score" < score.
%If score = 'BFDR', will use bonferonni calculation to ensure no false
%positive identifications. This will likely exclude A LOT of peptides.
%Using 'BHFDR' option will will use Benjamini-Hochberg adjustment of
%p-values and threshold for p<scorethreshold.

scorebf = (log10(1/(20*(norm_length - 3))) * -10) - 13;
if strcmp(score, 'BFDR')
    score = scorebf;
end

if strcmp(score, 'BHFDR')
    scores = 10.^(str2double(normedpepstonorm(4:end,8))./(-10));
    scores = [0;0;0; bhfdr(scores)];
    normedpepstonorm = normedpepstonorm(scores < scorethreshold,:);
    norm_length = size(normedpepstonorm,1);
elseif score > 0
    scored = str2double(normedpepstonorm(:,8)) > score;
    normedpepstonorm = normedpepstonorm(scored,:);
    normedpepstonorm = [pepstonorm(1:3,:); normedpepstonorm];
    norm_length = size(normedpepstonorm,1);
end
name = '';
for i = RA:RA-1+width
    if isempty(cell2mat(normedpepstonorm(2,i)))
        normedpepstonorm{2,i} = name;
    else
        name = normedpepstonorm{2,i};
    end
end
normedpepstonorm{3,gene_col} = 'GENE ID';

%% Step 1: Initial pass finds identifiers and review status of each protein.
%We'll check whether peptides are non-conflicting later in the code.
header = normedpepstonorm(1:3,:);
normedpepstonorm = [header;sortrows(normedpepstonorm(4:end,:),11)];
ii = 4;
while ii ~= norm_length + 1
    ProtNameptn = normedpepstonorm{ii,11};
    ProtFindptn = strcmp(normedpepstonorm(:,11), char(ProtNameptn));
    switch nDB
        case 1
            ProtName = char(ProtNameptn);
        otherwise
            ProtName = char(ProtNameptn(1,4:end));
    end
    UPFind = strcmp(uniprotall(:,2), ProtName);
    review_status = uniprotall(UPFind,end); %Pick reviewed protein ids
    protid = uniprotall(UPFind,[upcol,2]);
    protid = protid(strcmp(review_status,'reviewed'),:);
    try normedpepstonorm(ProtFindptn,[gene_col,11]) = repmat(protid(1,:),sum(ProtFindptn,1),1);
        review_status = 'reviewed';
    catch
        normedpepstonorm(ProtFindptn,[gene_col,11]) = {ProtName};
        review_status = 'unreviewed';
    end
    
    %NumPepsptn = numel(unique(normedpepstonorm(ProtFindptn,9)));
    normedpepstonorm(ii:end,1) = {review_status};
    fprintf(['#',num2str(ii),'-',num2str(ii+nansum(ProtFindptn, 1)),' ',ProtName,' ',char(review_status),'\n'])
    ii = ii + nansum(ProtFindptn, 1); 
    %waitbar(ii/norm_length, wait,...
     %   'Getting protein info...');
end   
normedpepstonorm = [header;sortrows(normedpepstonorm(4:end,:),gene_col)];

%% Step 2: Delete (convert to NaN) zeroes (log(0) = -Infinity), infinites.
%Get Log2(intensities) 
if strcmp(normmethod,'MSCMEF') || strcmp(normmethod,'none') || strcmp(normmethod,'nonewithimp')
    if iscellstr(peptidelistabds(4:pl_length, RA:RA-1+width))
        temp_peps = str2double(peptidelistabds(4:pl_length, RA:RA-1+width));
        %temp_impute = sqrt(nanmean(temp_peps,1))./1e9;
        %for a = 1:size(temp_peps,2)
        %    temp_peps(temp_peps(:,a) == 0,a) = temp_impute(a);
        %end
        temp_peps(isinf(temp_peps)) = NaN;
        %temp_peps(cell2mat(peptidelistabds(4:end,1)) < pepmin,:) = NaN;
        temp_peps = log2(temp_peps);
    else
        temp_peps = cell2mat(peptidelistabds(4:pl_length, RA:RA-1+width));
    end
    
    if iscellstr(normedpepstonorm(4:norm_length, RA:RA-1+width))
        temp_normed = str2double(normedpepstonorm(4:norm_length, RA:RA-1+width));
        temp_normed(temp_normed == 0) = NaN;
        temp_normed(isinf(temp_normed)) = NaN;
        %temp_normed(cell2mat(normedpepstonorm(4:end,1)) < pepmin,:) = NaN;
        temp_normed = log2(temp_normed);
    else
        temp_normed = cell2mat(normedpepstonorm(4:norm_length, RA:RA-1+width));
    end
    
elseif ~strcmp(normmethod,'MSCMEF')
    temp_peps = str2double(peptidelistabds(4:pl_length, RA:RA-1+width));
    temp_peps(isinf(temp_peps)) = NaN;
    temp_normed = str2double(normedpepstonorm(4:norm_length, RA:RA-1+width));
    temp_normed(temp_normed == 0) = NaN;
    temp_normed(isinf(temp_normed)) = NaN;
    %temp_normed(cell2mat(normedpepstonorm(4:end,1)) < pepmin,:) = NaN;
    temp_peps = log2(temp_peps);
    temp_normed = log2(temp_normed);
    
end

%% Step 3: Normalise pepstonorm by wholepeptidelist mean/median/sum/none.
%In log-space this is done by subtraction. Subtraction in log-space is
%the same as division in abundance-space
switch normmethod
    case 'mean'
        means = nanmean(temp_peps, 1);
        %overall_mean = nanmean(nanmean(temp_peps, 1),2);
        temp_normed = bsxfun(@minus, temp_normed, means);
        
    case 'median'
        medians = nanmedian(temp_peps, 1);
        %overall_median = nanmedian(nanmedian(temp_peps, 1),2);
        temp_normed = bsxfun(@minus, temp_normed, medians);
        
    case 'sum'
        sums = nansum(temp_peps, 1);
        temp_normed = bsxfun(@minus, temp_normed, sums);
        
    case 'MSCMEF'
        medians1 = nanmedian(temp_peps, 1);
        %medians2 = nanmedian(temp_normed, 1);
        temp_normed = temp_normed - repmat(medians1, norm_length-3, 1);%...
            %+ repmat(medians2, norm_length-3, 1);
    case 'none'
        warning('No normalisation used')
    case 'nonewithimp'
        warning('No normalisation used')
    otherwise
        error('normmethod not set correctly')
end

GroupList = unique(normedpepstonorm(2,RA:1:RA-1+width), 'stable');
GroupNum = numel(GroupList);
if numDonors > 1
    h = numDonors;
else
    h = numel(donors)/GroupNum;
end
%% Step 4: Find unique peptides, including PTMs.
%Re-assign peptides from unreviewed isoforms to reviewed if sequence
%matches.
norm_length = size(normedpepstonorm,1);
AllPeptideIDs = cell(norm_length,20);

for i = 4:norm_length
   AllPeptideIDs(i,1) = {cat(2, cell2mat(normedpepstonorm(i,9)),...
       cell2mat(normedpepstonorm(i,10)))};
   AllPeptideIDs(i,2:RA) = normedpepstonorm(i,1:RA-1);
   AllPeptideIDs(i,20) = normedpepstonorm(i,gene_col+(~ProteinGrouping)*4);
end
AllUniquePeptideIDs = unique(AllPeptideIDs(4:end,10), 'stable');
UniquePepNum = length(AllUniquePeptideIDs);
temp_final = zeros(UniquePepNum, width);
temp_info = cell(UniquePepNum, 18);
q = 1;
for i = 1:UniquePepNum
    p_finder = strcmp(AllPeptideIDs(:,10),char(AllUniquePeptideIDs(i,1))); %Look for all matches to sequence
    p_info = AllPeptideIDs(p_finder, [2:11,20,13:19]);
    p_values = temp_normed(p_finder(4:end),:);
    
    %Reassign peptides from unreviewed proteins if possible
    if any(strcmp(p_info(:,1),'unreviewed'),1)
        protpepfind = strfind(uniprotall_reviewed(:,11),char(p_info(1,9)));
        protindex = find(~cellfun('isempty',protpepfind));
        if ~isempty(protindex)
            new_ids = uniprotall_reviewed(protindex,[2+ProteinGrouping*10,upcol]);
            is_present = zeros(norm_length,size(new_ids,1));
            for ii = 1:size(new_ids,1)
                %Reassign to proteins already in dataset if possible
                is_present(:,ii) = strcmp(AllPeptideIDs(:,1+gene_col+(~ProteinGrouping)*4),char(new_ids(ii,1+ProteinGrouping)));
            end
            is_present = sum(is_present,1);
            [ic,ia] = max(is_present,[],2);
            if ic ~= min(is_present,[],2) || size(is_present,2) == 1; %If peptide could be shared between 2 or more equally likely reviewed proteins, skip reassignment.
                p_info(:,[11,7]) = repmat(new_ids(ia,:),size(p_info,1),1);
                p_info(:,1) = repmat({'reviewed'},size(p_info,1),1);
                AllPeptideIDs(p_finder,[12,gene_col+1]) = repmat(new_ids(ia,:),size(p_info,1),1);
            end
        end
    end
    
    p_names = p_info(:,7+(~ProteinGrouping)*4);
    if numel(unique(p_names)) > 1
        continue; %Skip if peptide used in more than 1 entity
    end
    
    %if ~strcmp(regmethod,'bayes')
     %   p_values = nanmedian(p_values,1);
      %  p_info = p_info(1,:);
    %end
     
    temp_final(q:q-1+size(p_values,1),:) = p_values;
    temp_info(q:q-1+size(p_values,1),:) = p_info;
    fprintf(['#',num2str(q),' ',char(p_info(1,9)),' ',char(p_info(1,10)),'\n'])
    q = q + size(p_values,1);
    %waitbar((i/UniquePepNum), wait,'Finding useable peptides...');
end
temp_final = temp_final(1:q-1,:);
temp_info = temp_info(1:q-1,:);
PCA.NP.missing = isnan(temp_final); %Record indices of missing values

%Impute missing values as random values drawn from a normal distribution
%calculated from the mean and std of each peptide, median shifted
%by -1.6 * std. For Bayesian Linear regression, missing values are imputed
%as part of the Gibbs sampler.
if (~strcmp(normmethod,'none') && ~strcmp(regmethod,'bayes')) && ~strcmp(regmethod,'nnnn')
    temp_final = Impute(temp_final, 'MI');
    temp_final = temp_final(~isnan(temp_final(:,1)),:);
    temp_info = temp_info(~isnan(temp_final(:,1)),:);
end

%Reassign normalised peptide intensities back to normedpepstonorm cell
%array for output.
q = size(temp_info,1);
normedpepstonorm(4:q + 3,1:18) = temp_info;
normedpepstonorm(4:q + 3,RA:RA-1+width) = num2cell(temp_final);
normedpepstonorm = normedpepstonorm(1:q + 3,:);

%% Step 5: Create table for fitlme
%Create vector of donor/sample IDs
norm_length = size(normedpepstonorm,1);
donors = repmat(donors, norm_length-3, 1);
donors = donors(:);

%Create vector of group names
Group = normedpepstonorm(2,RA:1:RA-1+width);
Group = repmat(Group, norm_length-3, 1);
Group = Group(:);

sID = repmat(RA:1:RA-1+width, norm_length-3, 1);
sID = sID(:);
Subject = cellstr(cat(2,char(Group),num2str(sID)));

%Create vector of feature IDs for all samples and groups
Features = cell(norm_length-3,1);
for i = 4:norm_length
    Features(i-3,1) = {cat(2, cell2mat(normedpepstonorm(i,9)),...
        cell2mat(normedpepstonorm(i,10)))};
end
Features = repmat(Features, width, 1);

%Get vecor of protein sequences for determining unique peptides
Seqs = repmat(normedpepstonorm(4:end,9),width,1);

%Vector of Mascot scores
Scores = repmat(str2double(normedpepstonorm(4:end,8)),width,1);
Scores = Scores./scorebf;
Scores(Scores >= 1) = 1;
%Scores = Scores.^2;
Scores(isnan(Scores)) = 0;

%Create vector of peptide intensities
meanIntensities = nanmean(temp_final,2); %Used to model missingness
meanIntensities = repmat(meanIntensities,width,1);
Intensities = temp_final(:);

%Get protein names (not used in equation, only for selecting
%peptides to fit model for each individual protein).
switch ProteinGrouping
    case true
        Proteins = normedpepstonorm(4:end,7);
        %upCol = 12;
        IDcol = 7;
    case false
        Proteins = normedpepstonorm(4:end,11);
        upcol = 2;
        IDcol = 11;
end
UniqueProteins = unique(Proteins, 'stable');
ProtNum = numel(UniqueProteins);
Proteins = repmat(Proteins, width,1);

%Create table for fitlme()
lmetable = table(Proteins, donors, Group, Features, Intensities,...
    Subject, Scores, meanIntensities, Seqs,...
    'VariableNames', {'Protein', 'Donor', 'Group', 'Feature',...
    'Intensity', 'Subject', 'Scores', 'Means', 'Sequences'});
lmetable.Donor = categorical(lmetable.Donor);

%Make arrays for preloading
ProteinAbundance = cell(length(UniqueProteins)+2,8 + GroupNum * 4);
%LM = cell(length(UniqueProteins),1);
[SortedGroupList,I] = sort(GroupList);

%Assign labels to tables
ProteinAbundance{2,1} = 'IDENTIFIER USED DURING ANALYSIS';
ProteinAbundance{2,2} = 'NUMBER OF PEPTIDES';
if strcmp(regmethod,'bayes') || strcmp(regmethod,'ols')
    ProteinAbundance(2,3:2 + (GroupNum * 3)) = repmat(SortedGroupList,1,3);
else
    ProteinAbundance(2,3:2 + (GroupNum * 3)) = repmat(GroupList,1,3);
end
ProteinAbundance{1,3} = 'EFFECT SIZE';
ProteinAbundance{1,3 + GroupNum} = 'S.E.';
ProteinAbundance{1,3 + (GroupNum * 2)} = 'P VALUE (VS. CTRL)';
ProteinAbundance(2,3 + (GroupNum * 3)) = {'UNIPROT ID'};
ProteinAbundance(2,4 + (GroupNum * 3)) = {'Protein Name'};
ProteinAbundance(2,5 + (GroupNum * 3)) = {'GENE NAME'};
if strcmp(regmethod,'bayes')
    ProteinAbundance(2,6 + (GroupNum * 3):5 + (GroupNum * 4)) = SortedGroupList;
else
    ProteinAbundance(2,6 + (GroupNum * 3):5 + (GroupNum * 4)) = GroupList;
end
ProteinAbundance{1,6 + (GroupNum * 3)} = 'FDR-ADJUSTED P VALUE (VS. CTRL)';
ProteinAbundance{1,end-2} = 'NUMBER OF FEATURE:GROUP COEFFICIENTS > SEM';
ProteinAbundance{1,end-1} = 'DEGREES OF FREEDOM';
ProteinAbundance{1,end} = 'MSE,INTERCEPT';

InteractionTable = [normedpepstonorm(:,1:RA-1), cell(q+3, GroupNum*4+2)];
InteractionTable(3,RA:RA-1 + (GroupNum * 4)) = repmat(SortedGroupList,1,4);
InteractionTable{2,RA} = 'EFFECT SIZE';
InteractionTable{2,RA + GroupNum} = 'S.E.';
InteractionTable{2,RA + (GroupNum * 2)} = 'P VALUE (VS. CTRL)';
InteractionTable{2,RA + (GroupNum * 3)} = 'FDR-ADJUSTED P VALUE (VS. CTRL)';
InteractionTable{3,end-1} = 'DEGREES OF FREEDOM';
InteractionTable{3,end} = 'MSE';

%% Step 6: Fit model for each protein
q = 1;
warning('off', 'stats:LinearModel:RankDefDesignMat');
%delete(wait); %Delete waitbar
%pepCounter = 4;
for i = 1:ProtNum
    Prottbl = lmetable(strcmp(lmetable.Protein, UniqueProteins(i)), :);
    Prottbl.Feature = categorical(Prottbl.Feature);
    ProtID = char(UniqueProteins(i));
    Y = Prottbl.Intensity;% - nanmean(Prottbl.Intensity);
    intind = find(I == 1);
    [pepC,pepia] = unique(Prottbl.Feature);
    NumPepsU = numel(unique(Prottbl.Sequences));
    NumPeps = numel(pepC);
    if NumPepsU < pepmin
       %warning('skipping protein with fewer than pepmin useable peptides');
       lmetable = lmetable(~strcmp(lmetable.Protein, UniqueProteins(i)),:);
       InteractionTable = InteractionTable(~strcmp(InteractionTable(:,IDcol), UniqueProteins(i)),:);
       normedpepstonorm = normedpepstonorm(~strcmp(normedpepstonorm(:,IDcol), UniqueProteins(i)),:);
       PCA.NP.missing = PCA.NP.missing(~strcmp(normedpepstonorm(4:end,IDcol), UniqueProteins(i)),:);
       continue;
    end 
    
    if isequal(ProtID,'') || isempty(ProtID)
        %warning('skipping unannotated protein');
        lmetable = lmetable(~strcmp(lmetable.Protein, UniqueProteins(i)),:);
        InteractionTable = InteractionTable(~strcmp(InteractionTable(:,IDcol), UniqueProteins(i)),:);
        normedpepstonorm = normedpepstonorm(~strcmp(normedpepstonorm(:,IDcol), UniqueProteins(i)),:);
        PCA.NP.missing = PCA.NP.missing(~strcmp(normedpepstonorm(4:end,IDcol), UniqueProteins(i)),:);
        continue;
    end
    %if ProteinGrouping == false && nDB > 1
     %   ProtID = ProtID(1,4:end);
   % end
    %waitmsg = cat(2,'Fitting model for ', ProtID);
    %waitbar(i/ProtNum, wait,waitmsg);
    
    switch model %Determine model parameters
        case 'full'
            %design matrix = Group + Feature + Donor + Feature:Group + Feature:Donor
            dummy = [dummyvar(categorical(Prottbl.Group)),...
                dummyvar(categorical(Prottbl.Feature))];
            dummy(:,GroupNum+1) = [];
            % Make variable IDs for easy searching
            dummyID = [repmat({'Group'},1,GroupNum), repmat({'Feature'},1,NumPeps-1)];
            
            interactions = zeros(size(Prottbl,1),1);
            if numDonors > 1
                dummy = [dummy, dummyvar(categorical(Prottbl.Donor))];
                dummy(:,GroupNum+NumPeps) = [];
                dummyID = [dummyID, repmat({'Donor'},1,numDonors-1)];
                %dummy(:,GroupNum+NumPeps+1) = 1;
                for iii = 0:1:numDonors-2
                    interactions = [interactions,...
                        bsxfun(@times, dummy(:,GroupNum+NumPeps+iii),...
                        dummyvar(categorical(Prottbl.Feature)))];
                        %[ones(size(Y,1),1),dummy(:,GroupNum+1:GroupNum+NumPeps-1)])];
                        %dummy(:,GroupNum+1:GroupNum+NumPeps-1))];
                end
                %interactions = [interactions, [ones(size(Y,1),1),dummy(:,GroupNum+1:GroupNum+NumPeps-1)]]; %???
                dummyID = [dummyID, repmat({'Feature:Donor'},1,(NumPeps)*(numDonors-1))];
            end
             
            for iii = 1:GroupNum
                %if iii == intind; 
                    %interactions = [interactions, [ones(size(Y,1),1),dummy(:,GroupNum+1:GroupNum+NumPeps-1)]];     
                %else
                if iii ~= intind
                    interactions = [interactions, bsxfun(@times, dummy(:,iii),...
                        dummyvar(categorical(Prottbl.Feature)))];
                    %[ones(size(Y,1),1),dummy(:,GroupNum+1:GroupNum+NumPeps-1)])];
                end
            end
            interactions(:,1) = [];
            dummyID = [dummyID, repmat({'Feature:Group'},1,(NumPeps)*(GroupNum-1))];
            dummy = [dummy, interactions]; 
            
            if incSubject %Add subject term if necessary
                dummy = [dummy, dummyvar(categorical(Prottbl.Subject))];
                dummy(:,end-numel(unique(Subject))) = [];
                dummyID = [dummyID, repmat({'Subject'},1,numel(unique(Subject))-1)];
            end
            
            %if numel(unique(Prottbl.Scores)) > 1 && ProteinGrouping == true %Add protein term if necessary
            %    dummy = [dummy, dummyvar(categorical(Prottbl.Scores))];
            %    dummyID = [dummyID, repmat({'Protein'},1,numel(unique(Prottbl.Scores)))];
            %end
        case 'simple' %design matrix = Group + Feature
            dummy = [dummyvar(categorical(Prottbl.Group)),...
                dummyvar(categorical(Prottbl.Feature))];
            dummy(:,GroupNum+1) = [];
            dummyID = [repmat({'Group'},1,GroupNum),...
                repmat({'Feature'},1,NumPeps-1)];
            
            if numDonors > 1
                dummy = [dummy, dummyvar(categorical(Prottbl.Donor))];
                dummy(:,GroupNum+NumPeps) = [];
                dummyID = [dummyID, repmat({'Donor'},1,numDonors-1)];
            end
            %dummy(:,find(I == 1)) = 1; %design matrix = Intensity ~ Group + Feature   
            %dummy(:,GroupNum+1) = 1;
            %dummy(:,GroupNum+NumPeps+1) = 1;
        otherwise
            error('model arguement not set properly - must be either ''simple'' or ''full''')
    end
    %dummy(:,GroupNum+1) = 1;
    %dummy(:,find(I == 1)) = [];
    %dummyID(:,find(I == 1)) = [];
    %dummyID(1,GroupNum+1) = {'Intercept'};
    %dummy(:,GroupNum+1) = 1;
    
    dummy = [dummy,ones(size(Y,1),1)]; %Add intercept term
    dummyID = [dummyID,{'Intercept'}];
    
    switch regmethod %Fit model
        case 'ME' % Mixed-effects model/similar to ridge regression
            if numel(unique(Prottbl.Feature)) == 1
                formula = 'Intensity ~ Group';
            else
                if incSubject %Add subject term if necessary
                    switch numDonors
                        case 1
                            switch model
                                case 'full'
                                    formula = 'Intensity ~ Group + (1|Subject) + (1|Feature) + (1|Feature:Group)';
                                otherwise
                                    formula = 'Intensity ~ Group + (1|Subject) + (1|Feature)';
                            end
                        otherwise
                            switch model
                                case 'full'
                                    formula = 'Intensity ~ Group + (1|Donor) + (1|Subject) + (1|Feature) + (1|Feature:Group) + (1|Feature:Donor)';
                                otherwise
                                    formula = 'Intensity ~ Group + (1|Donor) + (1|Subject) + (1|Feature)';
                            end
                    end
                else
                    switch numDonors
                        case 1
                            switch model
                                case 'full'
                                    formula = 'Intensity ~ Group + (1|Feature) + (1|Feature:Group)';
                                otherwise
                                    formula = 'Intensity ~ Group + (1|Feature)';
                            end
                        otherwise
                            switch model
                                case 'full'
                                    formula = 'Intensity ~ Group + (1|Donor) + (1|Feature) + (1|Feature:Group) + (1|Feature:Donor)';
                                otherwise
                                    formula = 'Intensity ~ Group + (1|Donor) + (1|Feature)';
                            end
                    end
                end
            end
            
            index = 1;
            MixedLME = {[]};
            while isempty(MixedLME{1,1})
                %Do fit
                try MixedLME{1,1} = fitglme(Prottbl, formula,...
                    'FitMethod', 'REMPL',...
                    'StartMethod', 'random',...
                    'OptimizerOptions', statset(statset('fitglme'),statset('RobustWgtFun','huber','Tune',1.345)));
                            
                    [~,~,S] = randomEffects(MixedLME{1,1});
                    PepGroupInteractions = S(strcmp(S.Group,'Feature:Group'),:);
                    interresults = PepGroupInteractions.Estimate;
                    interBetas = reshape(interresults, [NumPeps,GroupNum]);
                    interSEMs = reshape(PepGroupInteractions.SEPred, [NumPeps,GroupNum]);
                    interPs = reshape(PepGroupInteractions.pValue, [NumPeps,GroupNum]);
                    interresults = [interBetas, interSEMs, interPs];
                    temp = InteractionTable(strcmp(InteractionTable(:,IDcol),UniqueProteins(i)),RA:RA-1+GroupNum*3);
                    temp(pepia,:) = num2cell(interresults);
                    
                    num_interresults = sum(abs(PepGroupInteractions.Estimate) > PepGroupInteractions.SEPred,1);
                    ProteinAbundance(q+2,end-2) = {num_interresults};
                    data = dataset2cell(MixedLME{1,1}.Coefficients);
                    intercept = cell2mat(data(1,1));
                    results = [0;cell2mat(data(3:end,2))]';
                    SEMs = cell2mat(data(2:end,3))';
                    Ps = [1;cell2mat(data(3:end,6))]';
                    DF = cell2mat(data(2,5));
                    MSE = MixedLME{1,1}.SSE;
                    
                    InteractionTable(strcmp(InteractionTable(:,IDcol),UniqueProteins(i)),RA:RA-1+GroupNum*3) = temp;
                    InteractionTable(strcmp(InteractionTable(:,IDcol),UniqueProteins(i)),end-1) = {DF};
                    InteractionTable(strcmp(InteractionTable(:,IDcol),UniqueProteins(i)),end) = {MSE};
                catch
                    if index > 1000
                        msg = cat(2,'fitglme failed on ', ProtID);
                        warning(msg);
                        break;
                    end
                    index = index + 1;
                    continue;
                end
            end
            if index > 1000
                continue;
            end
        case 'bayes' %Bayesian linear regression
            [n,~] = size(dummy);
            
            %Model missing value types
            missing = isnan(Y);
            missing2 = 10.*sign(missing-0.5);
            F = [dummy(:,1:GroupNum + NumPeps-1),ones(n,1)]; %Intercept denotes intrinsic missingness probability
            missingmdl = weighted_bayeslm(F,missing2,dummyID(1,[1:GroupNum + NumPeps-1,end]),false,ones(n,1),[]);
            MNR = (missingmdl.beta_estimate > 0 & missingmdl.beta_estimate > missingmdl.beta_estimate(end));  %Missing intensities from MNR are deemed missing non-randomly
            %MNR(end) = 0;%All proteins have a non-zero chance of being missing, doesn't mean that all peptides are MNR
            Y_MNR = F(missing,MNR); 
            Y_MNR = any(Y_MNR,2);
            
            %Fit protein model
            mdl = weighted_bayeslm(dummy,Y,dummyID,true,Prottbl.Scores,Y_MNR);
            results = mdl.beta_estimate(1,1:GroupNum)';
            DF = mdl.DoF;
            MSE = mdl.residErr;
            SEMs = mdl.SEPred;
            Ps = mdl.Pvalues;
            intercept = mdl.beta_estimate(1,strcmp(dummyID,'Intercept'))';
            %if numel(unique(Subject)) / GroupNum > 1
            
            %else
             %   interresults = mdl.beta_estimate(1,strcmp(dummyID,'Feature:Group'))';
              %  num_interresults = numel(interresults(abs(interresults) > SEMs(strcmp(dummyID,'Feature:Group'))'));
            %end
            
            %Setup things for PTM quantification using interaction
            %coefficients
            interresults = mdl.beta_estimate(1,strcmp(dummyID,'Feature:Group'))';
            num_interresults = numel(interresults(abs(interresults) > SEMs(1,strcmp(dummyID,'Feature:Group'))'));
            interBetas = zeros(NumPeps,GroupNum);
            interBetas(:,[1:intind-1,intind+1:end]) = reshape(interresults, [NumPeps,GroupNum-1]);
            %interBetas = reshape(interresults, [NumPeps,GroupNum]);
            interBetas_ctrl = interBetas(:,intind);
            interBetas = bsxfun(@minus, interBetas, interBetas_ctrl); %Easier to do contrasts here
            interSEMs = zeros(NumPeps,GroupNum);
            interSEMs(:,[1:intind-1,intind+1:end]) = reshape(SEMs(1,strcmp(dummyID,'Feature:Group')), [NumPeps,GroupNum-1]);
            %interSEMs = reshape(SEMs(1,strcmp(dummyID,'Feature:Group')), [NumPeps,GroupNum]);
            interPs = (1 - tcdf(abs(interBetas)./interSEMs, DF - 1)) .* 2;%reshape(Ps(1,strcmp(dummyID,'Feature:Group')), [NumPeps,GroupNum]);
            interresults = [interBetas, interSEMs, interPs];
            temp = InteractionTable(strcmp(InteractionTable(:,IDcol),UniqueProteins(i)),RA:RA-1+GroupNum*3);
            temp(pepia,:) = num2cell(interresults);
            InteractionTable(strcmp(InteractionTable(:,IDcol),UniqueProteins(i)),RA:RA-1+GroupNum*3) = temp;
            InteractionTable(strcmp(InteractionTable(:,IDcol),UniqueProteins(i)),end-1) = {DF};
            InteractionTable(strcmp(InteractionTable(:,IDcol),UniqueProteins(i)),end) = {MSE};
            
            %Update normalised peptide intensities with imputed values
            normedpepstonorm(strcmp(normedpepstonorm(:,IDcol), UniqueProteins(i)),RA:end) = num2cell(reshape(mdl.Yimputed,[numel(Prottbl.Feature)/width,width]));
            %pepCounter = pepCounter+NumPeps;
            ProteinAbundance(q+2,end-2) = {num_interresults};
                        
        case 'ols' %Ordinary least squares
            switch numDonors
                case 1
                    formula = 'Intensity ~ Group + Feature';
                otherwise
                    formula = 'Intensity ~ Group + Donor + Feature';
            end
            InitLM = fitlm(Prottbl,formula);
            final_beta = InitLM.Coefficients.Estimate(numDonors+1:numDonors+GroupNum-1);
            intercept = InitLM.Coefficients.Estimate(1,1);
            SEMs = [InitLM.Coefficients.SE(1),InitLM.Coefficients.SE(numDonors+1:numDonors+GroupNum-1)'];
            SEMs = SEMs(1:GroupNum);
            Ps = [1,InitLM.Coefficients.pValue(numDonors+1:numDonors+GroupNum-1)'];
            Ps = Ps(1:GroupNum);
            results = [0;final_beta];
            DF = InitLM.DFE;
            MSE = InitLM.MSE;
        otherwise
            warning('regmethod not set properly, must be either ''ME'', ''bayes'' or ''ols''')
            delete(wait); %Delete waitbar
            %normedpepstonorm = [pepstonorm(1:2,:); normedpepstonorm];
            return;
    end
    
    ProteinAbundance(q+2,1) = UniqueProteins(i);
    ProteinAbundance{q+2,2} = NumPepsU;
    ProteinAbundance(q+2,3:2+GroupNum) = num2cell(results)';
    ProteinAbundance(q+2,3 + GroupNum:2 + (GroupNum * 2)) = num2cell(SEMs(1:GroupNum));
    pvals = num2cell(Ps(1:GroupNum));
    pvals(find(I==1)) = {1};
    ProteinAbundance(q+2,3 + GroupNum*2:2 + (GroupNum * 3)) = pvals;
    
    %Find protein in uniprotall array
    uniprotallfinder = strcmp(uniprotall(:,upcol), ProtID);
    info2 = unique(uniprotall(uniprotallfinder, 3));
    info3 = unique(uniprotall(uniprotallfinder, 4));
    try ProteinAbundance(q+2,end-(5+GroupNum):end-(GroupNum+3)) = ...
            [{unique(uniprotall(uniprotallfinder, 1))}, info2(1,1), ...
            info3(1,1)];
    catch
        msg = cat(2,'No data in uniprotall file for ', ProtID,...
       '. Protein is unreviewed.');
        warning(msg);
    end
    
    %Get degress of freedom, interaction betas and MSE for final model
    ProteinAbundance(q+2,end-1) = {DF};
    ProteinAbundance(q+2,end) = {[MSE,intercept]};
    %P = ProteinAbundance(q+2,:)';
    %P(2:end+1,:) = P;
    %P(1,:) = {max(cellfun('length',ProteinAbundance(:,1)))};
    %fprintf(['%*s %3d ' ,repmat('%3.2f ',1,3*3),'%3s %3s %3s ', repmat('%3.2f ',1,3),'%3d %3d %3d\n'], P{:})
    fprintf(['#',num2str(q),' ',ProtID,' ',num2str(NumPepsU),' ',num2str(results(:)',2),'\n'])
    q = q + 1;
end

ProteinAbundance = ProteinAbundance(1:q + 1, :);
t = zeros(size(InteractionTable,1),1);
for i = 1:size(InteractionTable,1)
    t(i) = ~isempty(cell2mat(InteractionTable(i,RA)));
end
InteractionTable = InteractionTable(~~t,:);
if strcmp(regmethod,'bayes') || strcmp(regmethod,'ME')
    norm_length = size(InteractionTable,1);
else
    norm_length = size(normedpepstonorm,1);
end
%ProteinAbundance(3:end,2 + GroupNum + find(I == 1)) = {0};

%% Step 7: Do Contrasts and Empirical Bayes correction of SEMs and p vals
%Performed as detailed in Kammers et al. (2015) Detecting significant 
%changes in protein abundance. EuPA Open Proteomics 7, 11:19.
if strcmp(regmethod,'bayes')
    %ProteinAbundance(3:end,2 + find(I == 1)) = {0};
    ctrl_betas = cell2mat(ProteinAbundance(3:end,2 + find(I == 1)));
end
if EB %Do Empirical Bayes
    %ctrl_intbetas = cell2mat(InteractionTable(4:end,RA-1 + find(I == 1)));
    DoFs = cell2mat(ProteinAbundance(3:end,end-1));
    ds = zeros(GroupNum,2);
    for i = 0:1:GroupNum-1
        %Get hyperparameters from chi-squared distribution fitted to current
        %SEM estimates.
        %if i == find(I == 1)-1
        %   continue; %Skip control group
        %end
        d0s0 = double.empty(0,2);
        iter = 1;
        SEs = cell2mat(ProteinAbundance(3:end,3 + GroupNum + i));
        %SEuseable = cell2mat(ProteinAbundance(3:end,3 + GroupNum + i)) < 100;
        %SEs = SEs(SEuseable,:);
        if strcmp(regmethod,'bayes')
           %Do constrasts for all groups
           ProteinAbundance(3:end,3 + i) = ...
               num2cell(cell2mat(ProteinAbundance(3:end,3 + i)) - ctrl_betas);
        %if i == 0, continue;
        end
    
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
        for ii = 3:q+1
            DoF = DoFs(ii-2);
            lambda = DoF/(DoF + d0s0(1));
            oldSEM = cell2mat(ProteinAbundance(ii,3 + GroupNum + i));
            oldVar = (oldSEM/sqrt(2/((DoF/2) + 2))).^2;
            newSEM = sqrt(lambda * oldVar + (1 - lambda) * d0s0(2))...
                * sqrt(2/((DoF + d0s0(1))/2 + 2));
            ProteinAbundance(ii,3 + GroupNum + i) = {newSEM};
            FC = abs(cell2mat(ProteinAbundance(ii,3+i)));
            t = FC/newSEM;
            p = min(1,(1 - tcdf(t, DoF + d0s0(1) - 1)) * 2 + 1e-15); %2-sided t-test, p-values cannot be 0
            ProteinAbundance(ii,3 + (GroupNum * 2) + i) = {p};
            %ProteinAbundance(ii,end-1) = {DoF + d0s0(1) - 1}; %New degrees of freedom
        end
        ds(i+1,:) = d0s0;
    end
end
ProteinAbundance{1,1} = ds;

%% Step 8: Make summary figures.
%Make protein correlation matrix with mean peptide intensities.
%PCA.MP.PreCorr = unstack(lmetable,'Intensity', 'Subject',...
%    'AggregationFunction', @LSum, 'GroupingVariables', 'Protein');
%PCA.MP.PreCorrOrder = PCA.MP.PreCorr.Properties.VariableNames(1,2:end);
%PCA.MP.PreCorr = table2cell(PCA.MP.PreCorr);
%PCA.MP.PreCorr = cell2mat(PCA.MP.PreCorr(:,2:end));
%[NormedCor(:,:,1), NormedCor(:,:,2)] = corrcoef(PCA.MP.PreCorr,'rows',...
%    'pairwise');

%Make fold change histograms
figure('Name', 'Fold change distributions')
for i = 1:GroupNum
    if i ~= find(I == 1)
        subplot(1,GroupNum,i)
        histogram(cell2mat(ProteinAbundance(3:end,2 + i)),19)
        hold on;
        xmin = nanmin(nanmin(cell2mat(ProteinAbundance(3:end,...
            3:2 + GroupNum)),[],1),[],2);
        xmax = nanmax(nanmax(cell2mat(ProteinAbundance(3:end,...
            3:2 + GroupNum)),[],1),[],2)+1;
        axis([xmin xmax 0 size(ProteinAbundance,1)/1.5]);
    end
end
hold off;

%Make 2D PCA using mean peptide intensities
%figure('Name', '2D PCA using mean peptide intensities')
%[PCA.MP.Coef,PCA.MP.Score,PCA.MP.Latent,PCA.MP.T,PCA.MP.Explained] = ...
%    pca(PCA.MP.PreCorr');
%for i = 1:size(PCA.MP.PreCorr,2)
%    scatter(PCA.MP.Score(i,1),PCA.MP.Score(i,2),'filled')
%    hold on;
%end
%hold off;

%Output correlation matrix using usable normalised peptide intensities.
try PCA.NP.PreCorr = cell2mat(normedpepstonorm(2:end,RA:end));
    PCA.NP.PreCorrOrder = normedpepstonorm(1,RA:end);
catch 
    PCA.NP.PreCorr = cell2mat(normedpepstonorm(4:end,RA:end));
    PCA.NP.PreCorrOrder = normedpepstonorm(3,RA:end);
end
%PCA.NP.PreCorr = table2cell(PCA.NP.PreCorr);
%PCA.NP.PreCorr = cell2mat(PCA.NP.PreCorr(:,2:end));
[NormedCor(:,:,1), NormedCor(:,:,2)] = corrcoef(PCA.NP.PreCorr,...
    'rows','pairwise');
figure('Name', 'Correlation Plot: normalised peptide intensities')
imagesc(NormedCor(:,:,1))
colorbar
caxis([-1,1])
n = 33;%fix(0.5*33);
r = [(0:1:n-1)/n,ones(1,n)];
g = [(0:n-1)/n, (n-1:-1:0)/n];
b = [ones(1,n),(n-1:-1:0)/n];
cmap = [r(:), g(:), b(:)];
colormap(cmap)
labs = PCA.NP.PreCorrOrder;
set(gca,'XTick', 0.5:width, 'YTick', 0.5:width, 'XTickLabel', labs,...
    'YTickLabel', labs, 'XTickLabelRotation', 45);
%{
%2D PCA using normalised peptide intensities
figure('Name', '2D PCA using normalised peptide intensities')
[PCA.NP.Coef,PCA.NP.Score,PCA.NP.Latent,PCA.NP.T,PCA.NP.Explained] =...
    pca(PCA.NP.PreCorr');
for i = 1:size(PCA.NP.Score,2)
    scatter(PCA.NP.Score(i,1),PCA.NP.Score(i,2),'filled')
    hold on;
end
hold off;
%}
figure('Name', 'Boxplots')
boxplot(cell2mat(ProteinAbundance(3:end,3:2+GroupNum)), 'OutlierSize',1)

%% Final step: Adjust p-values according to Benjamini & Hochberg  1995
%(doc mafdr for more info)
for i = 3 + (GroupNum * 2):2 + (GroupNum * 3)
    %if i ~= 2 + (GroupNum * 2) + find(I == 1)
    try temp_pvals = bhfdr(cell2mat(ProteinAbundance(3:end,i)))%,...
            %'BHFDR', true); 
        ProteinAbundance(3:end,i + GroupNum + 3) = num2cell(temp_pvals);
    catch
        warning('mafdr() failed, probably a licensing thing.');
    end
end

%% Pass over to PTM calculation code to output fold changes of PTMs.
PTMnum = 0;
if ~strcmp(mods(1),''); PTMnum = size(mods,2); end
ModList = cell(PTMnum,1);
switch regmethod
    case 'bayes'
        PrePTMList = InteractionTable;
        Groups = SortedGroupList;
        width = GroupNum*4;
    case 'ME'
        PrePTMList = InteractionTable;
        Groups = GroupList;
        width = GroupNum*4;
    otherwise
        PrePTMList = normedpepstonorm;
        
        if ~isequal(donors, RA:RA-1+width)
            %Preload donor matrix
            %DonorMatrix = zeros((norm_length-3), GroupNum, numDonors);
            %DMmeans = zeros((norm_length-3), numDonors);
            overall_mean = nanmean(PCA.NP.PreCorr, 2);
    
            %Populate donor matrix with values from temp_normed
            for i = 1:numDonors    
                DonorMatrix = PCA.NP.PreCorr(:,donorvector == donorvector(1,i));
                DMmeans = nanmean(DonorMatrix,2);
        
                %Do feature normalisation
                DonorMatrix = DonorMatrix...
                    - repmat(DMmeans,1,size(DonorMatrix,2))...
                    + repmat(overall_mean,1,size(DonorMatrix,2));
        
                %Reassign back to temp_normed
                PCA.NP.PreCorr(:,donorvector == donorvector(1,i)) = DonorMatrix;
            end
            try normedpepstonorm(2:end,RA:end) = num2cell(PCA.NP.PreCorr);
            catch 
                normedpepstonorm(4:end,RA:end) = num2cell(PCA.NP.PreCorr);
            end
        end
        Groups = GroupList;
end
for i = 1:PTMnum
    %waitmsg = cat(2,'Finding all ', mods{1,i},' PTMs.');
    %waitbar(i/PTMnum, wait,waitmsg);
    ModList{i,1} = CreatePTMList(PrePTMList, ProteinAbundance, ...
        norm_length, mods{1,i}, width, RA, IDcol, uniprotall, upcol, regmethod, Groups);
    
    interDoFs = cell2mat(ModList{i,1}(3:end,4));
    ds = zeros(GroupNum,2,PTMnum);
    for q = 0:1:GroupNum-1
        if 1+q == intind; continue; end
        %Repeat process for ModLists if using 'bayes' or 'ME'
        if (strcmp(regmethod,'bayes') || strcmp(regmethod,'ME')) && EB
            d0s0 = double.empty(0,2);
            iter = 1;
            SEs = cell2mat(ModList{i,1}(3:end,10 + GroupNum + q));
            %SEuseable = cell2mat(InteractionTable(4:end,RA + GroupNum + i)) < 100;
            %SEs = SEs(SEuseable,:);
        
            %Do constrasts for all groups
            %InteractionTable(4:end,RA + i) = ...
            %    num2cell(cell2mat(InteractionTable(4:end,RA + i)) - ctrl_intbetas);
      
            while isempty(d0s0)
                try d0s0 = mle(SEs,'pdf',@(x,v,s)chi2pdf(x/s,v)/s,'start', [rand,rand]);
                catch
                    if iter > 10000
                        warning(cat(2,'PTM SEM estimates whacky - cannot get hyperparameters! No Empirical Bayes on column ', num2str(1+q),'!'));
                        d0s0 = [0, 1];
                        break;
                    end
                    iter = iter + 1;
                    continue; 
                end
            end
    
            %Do correction for each PTM.
            for ii = 3:size(ModList{i,1},1)
                DoF = interDoFs(ii-2);
                if isempty(DoF)
                    continue;
                end
                lambda = DoF/(DoF + d0s0(1));
                oldSEM = cell2mat(ModList{i,1}(ii,10 + GroupNum + q));
                oldVar = (oldSEM/sqrt(2/((DoF/2) + 2))).^2;
                newSEM = sqrt(lambda * oldVar + (1 - lambda) * d0s0(2))...
                    * sqrt(2/((DoF + d0s0(1))/2 + 2));
                ModList{i,1}(ii,10 + GroupNum + q) = {newSEM};
                FC = abs(cell2mat(ModList{i,1}(ii,6+q)));
                t = FC/newSEM;
                p = min(1,(1 - tcdf(t, DoF + d0s0(1) - 1)) * 2 + 1e-15); %2-sided t-test, p-values cannot be 0
                ModList{i,1}(ii,10 + (GroupNum * 2) + q) = {p};
                %InteractionTable(ii,end-1) = {DoF + d0s0(1) - 1}; %New degrees of freedom
            end
            ds(q+1,:,i) = d0s0;
            try temp_pvals = bhfdr(cell2mat(ModList{i,1}(3:end,10 + (GroupNum * 2) + q)))%,...
                %'BHFDR', true);
                ModList{i,1}(3:end,q + 10 + (GroupNum * 3)) = num2cell(temp_pvals);
            catch
                warning('mafdr() failed, probably a licensing thing.');
            end
        end
    end
end
ModList{end+1,1} = ds;
%normedpepstonorm = [wholepeptidelist(1:2,:); normedpepstonorm];
end
