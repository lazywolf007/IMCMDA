clc;clear;
str = load(['.\1.miRNA-disease association\knowndiseasemirnainteraction.txt']);
[~,disease]=xlsread(['.\1.miRNA-disease association\disease name.xlsx']);
[~,miRNA]=xlsread(['.\1.miRNA-disease association\miRNA name.xlsx']);
% nd:the number of diseases
% nm:the number of miRNAs
% pp:the number of known diseae-miRNA associations
nd = max(str(:,2));
nm = max(str(:,1));
[pp,~] = size(str);
% FS:the functional similarity between m(i) and m(j)
% FSP:Functional similarity weighting matrix
% SS:the semantic similarity between d(i) and d(j).
% SSP:semantic similarity weighting matrix
%
FS = load(['.\4.miNA functional simialrity\functional similarity matrix.txt']);
FSP = load(['.\4.miNA functional simialrity\Functional similarity weighting matrix.txt']);             
SS1 = load(['.\2.disease semantic similarity 1\disease semantic similarity matrix 1.txt']);
SS2 = load(['.\3.disease semantic similarity 2\disease semantic similarity matrix 2.txt']);
SS = (SS1+SS2)/2;
SSP = load(['.\2.disease semantic similarity 1\disease semantic similarity weighting matrix1.txt']);
%interaction: adajency matrix for the disease-miRNA association network
%interaction(i,j)=1 means miRNA j is related to disease i
interaction = zeros(nd,nm);
for i = 1:pp
    interaction(str(i,2),str(i,1)) = 1;
end
[kd,km] = gaussiansimilarity(interaction,nd,nm);                   %calculate gaussiansimilarity
[sd,sm] = integratedsimilarity(FS,FSP,SS,SSP,kd,km);               % Integrated similarity for diseases and miRNAs
[score]=IMC(interaction,sd,sm,100);                                 %obtain the prediction score matrix
allresult(disease,miRNA,interaction,score);                          %write the prediction score to excel


