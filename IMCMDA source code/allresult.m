function allresult(disease,miRNA,interaction,score)
a=(interaction-1)*(-1);
score=score.*a;
[m,n]=size(score);
c=find(score~=0);
[~,k]=size(c);
p=0;
str=cell(k,3);
for i=1:m
    for j=1:n
        if score(i,j)~=0
            p=p+1;
            str(p,1)=disease(i,1);
            str(p,2)=miRNA(j,1);           
            str{p,3}=score(i,j);
        end
    end
end
str=sortrows(str,-3);                     
xlswrite('.\所有疾病的预测结果.xlsx',str);
end
