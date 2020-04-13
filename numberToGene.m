close all;
clear all;
clc;

protein_file='3047-C-2-protein_result.xlsx';
sequence_file='3047-C-2-sequence.xlsx';
main_file='3047-C-2.csv';
output_file='3047-C-2-output-rib-fixed.xlsx';

protein_result=xlsread(protein_file);
[~,~,sequence]=xlsread(sequence_file);
[main,~,raw]=xlsread(main_file);

%find genes
sequence_list=int32([sequence{:,1}]');

for i = 1 : size(main,1)
    protein_ID=main(i,3);
    result_idx = find(protein_result(:,2)==protein_ID);
    if ~isempty(result_idx)
        NCBI_ID=int32(protein_result(result_idx(1),1));
        gene_idx = find(sequence_list==NCBI_ID);
        if ~isempty(gene_idx) && string(sequence{gene_idx,2}) == "gene"
            gene(i) = {sequence{gene_idx,3}};
        else
            gene(i) = {'NA'};
        end
    else
        gene(i) = {'NA'};
    end
    
end
gene=gene';
gene_header={'gene'};
gene_combined=cat(1,gene_header,gene);
raw_combined_genes=cat(2,raw,gene_combined);

%count modifications
accession_list=main(:,3);
unique_accession=unique(accession_list);
raw_no_header=raw(2:end,:);
j=0;
k=0;
z=1;
substring='R(+10.01)';
accession_num(1)=main(1,3);
counts(1)=0;
totals(1)=0;
for i = 1 : size(main,1)
    if i ~= 1 && main(i-1,3) ~= main(i,3)
        z=z+1;
        totals(z)=1;
        counts(z)=0;
        k=0;
        accession_num(z)=main(i,3);
    else
        totals(z)=totals(z)+1;
    end
    counts(z)=counts(z)+contains(raw_no_header(i,4),substring);
    
end
accession_num=accession_num';
counts=counts';
totals=totals';

for i = 1 : size(accession_num,1)
    protein_ID=accession_num(i);
    result_idx = find(protein_result(:,2)==protein_ID);
    if ~isempty(result_idx)
        NCBI_ID=int32(protein_result(result_idx(1),1));
        gene_idx = find(sequence_list==NCBI_ID);
        if ~isempty(gene_idx) && string(sequence{gene_idx,2}) == "gene"
            gene_an(i) = {sequence{gene_idx,3}};
        else
            gene_an(i) = {'NA'};
        end
    else
        gene_an(i) = {'NA'};
    end
    
end
gene_an=gene_an';
gene_an=strtrim(gene_an);

%assign regulon
[~,~,sigA]=xlsread('...\SigA_Regulon.xlsx');
[~,~,sigE]=xlsread('...\SigE_Regulon.xlsx');
[~,~,sigF]=xlsread('...\SigF_Regulon.xlsx');
[~,~,sigG]=xlsread('...\SigG_Regulon.xlsx');
[~,~,sigH]=xlsread('...\SigH_Regulon.xlsx');
[~,~,sigK]=xlsread('...\SigK_Regulon.xlsx');
[~,~,spo0A]=xlsread('...\Spo0A_Regulon.xlsx');

[~,~,rib]=xlsread('...\ribosomes.xlsx');

sA=zeros(size(gene_an,1),1);
sE=zeros(size(gene_an,1),1);
sF=zeros(size(gene_an,1),1);
sG=zeros(size(gene_an,1),1);
sH=zeros(size(gene_an,1),1);
sK=zeros(size(gene_an,1),1);
s0A=zeros(size(gene_an,1),1);
R=zeros(size(gene_an,1),1);
sA_mod=0;
sA_tot=0;
sE_mod=0;
sE_tot=0;
sF_mod=0;
sF_tot=0;
sG_mod=0;
sG_tot=0;
sH_mod=0;
sH_tot=0;
sK_mod=0;
sK_tot=0;
s0A_mod=0;
s0A_tot=0;
R_mod=0;
R_tot=0;
for i = 1 : size(gene_an,1)
    if string(gene_an{i})=="NA"
        sA(i)=0;
        sE(i)=0;
        sF(i)=0;
        sG(i)=0;
        sH(i)=0;
        sK(i)=0;
        s0A(i)=0;
        R(i)=0;
    else
        for j = 1 : size(rib,1)
            if contains(rib{j},gene_an{i},'IgnoreCase',true)
                R(i)=1;
                R_mod=R_mod+counts(i);
                R_tot=R_tot+totals(i);
            end
        end
        for j = 1 : size(sigA,1)
            if contains(sigA{j},gene_an{i},'IgnoreCase',true)
                sA(i)=1;
                sA_mod=sA_mod+counts(i);
                sA_tot=sA_tot+totals(i);
            end
        end
        for j = 1 : size(sigE,1)
            if contains(sigE{j},gene_an{i},'IgnoreCase',true)
                sE(i)=1;
                sE_mod=sE_mod+counts(i);
                sE_tot=sE_tot+totals(i);
            end
        end
        for j = 1 : size(sigF,1)
            if contains(sigF{j},gene_an{i},'IgnoreCase',true)
                sF(i)=1;
                sF_mod=sF_mod+counts(i);
                sF_tot=sF_tot+totals(i);
            end
        end
        for j = 1 : size(sigG,1)
            if contains(sigG{j},gene_an{i},'IgnoreCase',true)
                sG(i)=1;
                sG_mod=sG_mod+counts(i);
                sG_tot=sG_tot+totals(i);
            end
        end
        for j = 1 : size(sigH,1)
            if contains(sigH{j},gene_an{i},'IgnoreCase',true)
                sH(i)=1;
                sH_mod=sH_mod+counts(i);
                sH_tot=sH_tot+totals(i);
            end
        end
        for j = 1 : size(sigK,1)
            if contains(sigK{j},gene_an{i},'IgnoreCase',true)
                sK(i)=1;
                sK_mod=sK_mod+counts(i);
                sK_tot=sK_tot+totals(i);
            end
        end
        for j = 1 : size(spo0A,1)
            if contains(spo0A{j},gene_an{i},'IgnoreCase',true)
                s0A(i)=1;
                s0A_mod=s0A_mod+counts(i);
                s0A_tot=s0A_tot+totals(i);
            end
        end
    end
end

%exclusive regulon combinations
veg_counts=0;
veg_tot=0;
emc_counts=0;
emc_tot=0;
lmc_counts=0;
lmc_tot=0;
efs_counts=0;
efs_tot=0;
lfs_counts=0;
lfs_tot=0;
mc_counts=0;
mc_tot=0;
fs_counts=0;
fs_tot=0;
mcex_counts=0;
mcex_tot=0;
fsex_counts=0;
fsex_tot=0;
for i = 1 : size(gene_an,1)
    if sA(i) || s0A(i) || sH(i) && ~sE(i) && ~sK(i) && ~sF(i) && ~sG(i)
        veg_counts=veg_counts+counts(i);
        veg_tot=veg_tot+totals(i);
    end
    if sE(i) && sK(i)
        if sA(i) || s0A(i) || sH(i) || sF(i) || sG(i)
            mcex_counts=mcex_counts+counts(i);
            mcex_tot=mcex_tot+totals(i);
        else
            mc_counts=mc_counts+counts(i);
            mc_tot=mc_tot+totals(i);
        end
    end
    if sF(i) && sG(i)
        if sA(i) || s0A(i) || sH(i) || sE(i) || sK(i)
            fsex_counts=fsex_counts+counts(i);
            fsex_tot=fsex_tot+totals(i);
        else
            fs_counts=fs_counts+counts(i);
            fs_tot=fs_tot+totals(i);
        end
    end
    if sE(i) && ~sK(i) && ~sF(i) && ~sG(i) && ~sA(i) && ~s0A(i) && ~sH(i)
        emc_counts=emc_counts+counts(i);
        emc_tot=emc_tot+totals(i);
    end
    if ~sE(i) && sK(i) && ~sF(i) && ~sG(i) && ~sA(i) && ~s0A(i) && ~sH(i)
        lmc_counts=lmc_counts+counts(i);
        lmc_tot=lmc_tot+totals(i);
    end
    if ~sE(i) && ~sK(i) && sF(i) && ~sG(i) && ~sA(i) && ~s0A(i) && ~sH(i)
        efs_counts=efs_counts+counts(i);
        efs_tot=efs_tot+totals(i);
    end
    if ~sE(i) && ~sK(i) && ~sF(i) && sG(i) && ~sA(i) && ~s0A(i) && ~sH(i)
        lfs_counts=lfs_counts+counts(i);
        lfs_tot=lfs_tot+totals(i);
    end
end




%make output file
combined_counts=cat(2,num2cell(accession_num),gene_an,num2cell(counts),num2cell(totals));
combined_counts=cat(2,combined_counts,num2cell(sA),num2cell(sE),num2cell(sF),num2cell(sG),num2cell(sH),num2cell(s0A));
header = {'accession number','gene','mod count','total count','sA','sE','sF','sG','sH','s0A'};
combined_counts=cat(1,header,combined_counts);

regulon_top=[sA_mod,sE_mod,sF_mod,sG_mod,sH_mod,sK_mod,s0A_mod,R_mod;sA_tot,sE_tot,sF_tot,sG_tot,sH_tot,sK_tot,s0A_tot,R_tot;sA_mod/sA_tot,sE_mod/sE_tot,sF_mod/sF_tot,sG_mod/sG_tot,sH_mod/sH_tot,sK_mod/sK_tot,s0A_mod/s0A_tot,R_mod/R_tot];
regulon_header={'sA','sE','sF','sG','sH','sK','s0A','Rib'};
regulon_combined=cat(1,regulon_header,num2cell(regulon_top));

regulon_sort=[veg_counts,emc_counts,lmc_counts,efs_counts,lfs_counts,mc_counts,mc_counts+emc_counts+lmc_counts,fs_counts,fs_counts+efs_counts+lfs_counts;veg_tot,emc_tot,lmc_tot,efs_tot,lfs_tot,mc_tot,mc_tot+emc_tot+lmc_tot,fs_tot,fs_tot+efs_tot+lfs_tot;veg_counts/veg_tot,emc_counts/emc_tot,lmc_counts/lmc_tot,efs_counts/efs_tot,lfs_counts/lfs_tot,mc_counts/mc_tot,(mc_counts+emc_counts+lmc_counts)/(mc_tot+emc_tot+lmc_tot),fs_counts/fs_tot,(fs_counts+efs_counts+lfs_counts)/(fs_tot+efs_tot+lfs_tot)];
regulon_sort_header={'vegetative','early mother cell','late mother cell','early forespore','late forespore','early&late mother cell exclusive','total motehr cell','early&late forespore exclusive','total forespore'};
regulon_sort_combined=cat(1,regulon_sort_header,num2cell(regulon_sort));


xlswrite(output_file,combined_counts,1);
xlswrite(output_file,regulon_combined,2);
xlswrite(output_file,regulon_sort_combined,3);



