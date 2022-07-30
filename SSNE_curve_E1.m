

%The example for GSE113036 data.

clear;
clc;
close all;
[profile,pipi]= xlsread('GSE113036_data.xlsx');
[empty,network]= xlsread('network.xlsx');


profile=log(profile+1);
normal=profile(:,2:77);
case_mprofile=profile(:,2:624);

edge_source=network(:,1);
edge_target=network(:,2);
node1=unique(network(:,1));


reference_sample_num=76;
case_mprofile=fillmissing(case_mprofile,'constant',0);
es=0.00000001;
cell_num=[76,86,48,283,61,69];

tempcase(:,1,1:cell_num(1))=case_mprofile(:,1:76);     % Stage IA 
tempcase(:,2,1:cell_num(2))=case_mprofile(:,77:162);   % Stage IB
tempcase(:,3,1:cell_num(3))=case_mprofile(:,163:210);   % Stage IIA
tempcase(:,4,1:cell_num(4))=case_mprofile(:,211:493);   % Stage IIB
tempcase(:,5,1:cell_num(5))=case_mprofile(:,494:554);  % Stage IIIA
tempcase(:,6,1:cell_num(6))=case_mprofile(:,555:623); % Stage IIIB



psize=size(tempcase);
Land_entropy=zeros(length(node1),283,6);
Land_entropy_mix=zeros(length(node1),283,6);
Land_entropy_1=zeros(length(node1),283,6);
Land_entropy_1_mix=zeros(length(node1),283,6);
stage=6;


for l=1:psize(2)
    for  s=1:cell_num(l)
         for na=1:length(node1)
             edges_list=[];
             edges_list_1=[];
             [cencter_liang,cencter_idd]= find(ismember(pipi,node1(na)));
             if isempty(cencter_liang)
                 Land_entropy(na,s,l)=0;
                 Land_entropy_1(na,s,l)=0;
                 Land_entropy_mix(na,s,l)=0;
                 Land_entropy_1_mix(na,s,l)=0;
                 continue;
             else
                 [liang,nei_idd]= find(ismember(edge_source,node1(na)));
                 nei_gene=edge_target(liang);
                 e=0;
                 for n=1:length(nei_gene)
                     [nei_liang,nie_idd]= find(ismember(pipi,nei_gene(n)));
                     if ~isempty(nei_liang)
                         e=e+1;
                         edges_list(e,:)=[cencter_liang nei_liang];
                     end
                 end
                 
                 [liang_1,nei_idd_1]= find(ismember(edge_target,node1(na)));
                 nei_gene_1=edge_source(liang_1);
                 e1=0;
                 for n=1:length(nei_gene_1)
                     [nei_liang_1,nie_idd_1]= find(ismember(pipi,nei_gene_1(n)));
                     if ~isempty(nei_liang_1)
                         e1=e1+1;
                         edges_list_1(e1,:)=[cencter_liang nei_liang_1];
                     end
                 end
             end
            delt_pcc=[];
            if e<2
                Land_entropy(na,s,l)=0;
                Land_entropy_1(na,s,l)=0;
                Land_entropy_mix(na,s,l)=0;
                Land_entropy_1_mix(na,s,l)=0;
                continue;
            end
            
            nn=0;
            for i=1:e
                curr_pcc=abs(corr(normal(edges_list(i,1),:)',normal(edges_list(i,2),:)'));
                temp_add_onecase1=[normal(edges_list(i,1),:),reshape(tempcase(edges_list(i,1),l,s),1,1)];
                temp_add_onecase2=[normal(edges_list(i,2),:),reshape(tempcase(edges_list(i,2),l,s),1,1)];
                curr_pcc_add_onecase=abs(corr(temp_add_onecase1',temp_add_onecase2'));
                delt=abs(curr_pcc_add_onecase-curr_pcc);
                p_val=SNN(delt,curr_pcc,reference_sample_num);
                if p_val<0.05
                     nn=nn+1;
                     delt_pcc(nn)=delt;
                end
            end
            
            delt_pcc_1=[];
            nn_1=0;
            for i=1:e1
                curr_pcc_1=abs(corr(normal(edges_list_1(i,1),:)',normal(edges_list_1(i,2),:)'));
                temp_add_onecase_1=[normal(edges_list_1(i,1),:),reshape(tempcase(edges_list_1(i,1),l,s),1,1)];
                temp_add_onecase_2=[normal(edges_list_1(i,2),:),reshape(tempcase(edges_list_1(i,2),l,s),1,1)];
                curr_pcc_add_onecase_2=abs(corr(temp_add_onecase_1',temp_add_onecase_2'));
                delt_1=abs(curr_pcc_add_onecase_2-curr_pcc_1);
                p_val_1=SNN(delt_1,curr_pcc_1,reference_sample_num);
                if p_val_1<0.05
                     nn_1=nn_1+1;
                     delt_pcc_1(nn_1)=delt_1;
                end
            end
            
            if isempty(delt_pcc) || (length(delt_pcc)==1)
                Land_entropy(na,s,l)=0;
                Land_entropy_mix(na,s,l)=0;
            else
                delt_pcc=delt_pcc/(sum(delt_pcc)+es);
                Land_entropy(na,s,l)=(-(1/log(length(delt_pcc)))*sum(delt_pcc.*log(delt_pcc+es)))*(abs(std([normal(cencter_liang,:),reshape(tempcase(cencter_liang,l,s),1,1)])-std(normal(cencter_liang,:))));
                Land_entropy_mix(na,s,l)=(length(delt_pcc)/(length(delt_pcc)+length(delt_pcc_1)))*Land_entropy(na,s,l);
            end
            

            if isempty(delt_pcc_1) || (length(delt_pcc_1)==1)
                Land_entropy_1(na,s,l)=0;
                Land_entropy_1_mix(na,s,l)=0;
            else
                delt_pcc_1=delt_pcc_1/(sum(delt_pcc_1)+es);
                Land_entropy_1(na,s,l)=(-(1/log(length(delt_pcc_1)))*sum(delt_pcc_1.*log(delt_pcc_1+es)))*(abs(std([normal(cencter_liang,:),reshape(tempcase(cencter_liang,l,s),1,1)])-std(normal(cencter_liang,:))));
                Land_entropy_1_mix(na,s,l)=(length(delt_pcc_1)/(length(delt_pcc)+length(delt_pcc_1)))*Land_entropy_1(na,s,l);
            end
            
         end    
    end
    l
end

   
save GSE113036_entropy.mat;



mix_Land_entropy=Land_entropy_mix+Land_entropy_1_mix;

Land_entropy_size=size(mix_Land_entropy);
case_result=zeros(283,6);
for t=1:Land_entropy_size(3)
    for case_num=1:cell_num(t)
        [sort_Land_entrop,idx]=sort(mix_Land_entropy(:,case_num,t),'descend');
        case_result(case_num,t)=sum(sort_Land_entrop(1:0.05*Land_entropy_size(1)));                    
    end
    result(t)=mean(case_result(1:cell_num(t),t));
end
figure;
t=[1 2 3 4 5 6];
plot(t,result,'r','LineWidth',3);
set(gca,'XTick',1:7);
B={'d0' 'd2' 'd7' 'd14' 'd21' 'd22'};
set(gca,'XTickLabel',B);
ylim([1,9])
set(gca,'Ylim',[1 9]);
% set(gca,'YTick',[0.08:0.01:0.1]);
xlabel('Stages');
ylabel('Entropy');
%plot(t,aver_comidx,'r','LineWidth',3);
title('Average Entropy for GSE113036');


