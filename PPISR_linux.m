function []=PPISR_linux(varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ISR methods for GWAS @AUTHOR MENG LUO  %
    % contact: czheluo@gmail.com             %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %seting the env path 
    % alias matlab='/mnt/d/linux/MATLAB2016b/bin/matlab -nodesktop -nosplash -singleCompThread -logfile `date +%Y_%m_%d-%H_%M_%S`.log -r'
    addpath([pwd,'/src/']);
    %prediction=ft_getopt(varargin, 'prediction', 'no');
    method=ft_getopt(varargin, 'method', 'randomly');%cross validation or  10 fold 
    matfile=ft_getopt(varargin, 'matfile', []);
    phefile = ft_getopt(varargin, 'phefile', 'phe.fam');
    genofile = ft_getopt(varargin, 'genofile', 'pop.traw');
    outfile = ft_getopt(varargin, 'outfile', 'pop.traw.mat');
    sample = ft_getopt(varargin, 'sample',[]);
    nSNP = ft_getopt(varargin, 'nSNP', []);
    ntrait = ft_getopt(varargin, 'ntrait', 1);
    chr = ft_getopt(varargin, 'nchr', []);
    opt_outresult = ft_getopt(varargin, 'opt_outresult', 'ISR.opt.outresult.txt');
    all_outresult = ft_getopt(varargin, 'all_outresult', 'ISR.outresult.txt');
    vcf = ft_getopt(varargin, 'vcf', []);
    bed = ft_getopt(varargin, 'bed', []);
    ncov = ft_getopt(varargin, 'ncov', []);
    IM = ft_getopt(varargin, 'IM', 1);
    sgv = ft_getopt(varargin, 'sgv', 0.05); % default bonferroni correction
    mdl = ft_getopt(varargin, 'model',1); 
    %matfile = string, the matlab data format if you already prepare you data  
    %phefile = string, can be any of file format split with "\t"(default = 'phe.fam')
    %genofile = string, .traw file format from plink (default = 'pop.traw')
    %outfile = string, save covert genotypes file name with any name you defined and save matlab format (default = 'pop.traw.mat')
    %sample = number, the number of individuals you want to analysis
    %nSNP = number, the number of SNPs.
    %ntrait = number, the number of traits.
    %chr = number, the number of chromosome.
    %opt_outresult = string, write the result to text file (default = 'ISR.opt.outresult.txt') 
    %all_outresult = string, write the result to text file (default = 'ISR.outresult.txt')  
    %vcf = string, the VCF file name.
    %bed = string, the bed file name.
    %ncov = number, the number of PCs covariates.
    %IM = impute missing genotype with mean and median value, '1' was the default method means and others was median.
    %sgv = number, the bonferroni correction for association tests results.
    %mdl = number,1 for linear model and 2 or 3 for nolinear model ; input('Using Model II(without square term 2) or Model III(with square term 3) 2/3? ');
    % Usage:
    %       matlab "PPISR_linux('phefile','../data/pop.fam','genofile','../data/pop.traw','sample',87,'nSNP',28228,'ntrait',1,'ncov',5),exit;"
    %       matlab "PPISR_linux('matfile','demo.mat','sample',798,'nSNP',92641),exit;"
           
    %default parameter
    %if ~exist('outfile')
    %    outfile = 'pop.traw.mat';
    %end
    %if ~exist('opt_outresult')
    %    opt_outresult = 'opt_outresult.txt';
    %end
    %if ~exist('all_outresult')
    %    all_outresult = 'all_outresult.txt';
    %end
    %maxNumCompThreads=1;
    % write a shell script for ISR getopt....
    %clear;clc;
    %maxNumCompThreads(1)=1;
    if ~isempty(vcf)
        system(['plink ',' --vcf ',vcf,' --recode A-transpose --out pop']);
        traw2mat(phefile,genofile,outfile,sample,nSNP,ntrait,IM);
	load(outfile);
    elseif ~isempty(bed)
        system(['plink ',' --bfile ',bed,' --recode A-transpose --out pop ']);
        traw2mat(phefile,genofile,outfile,sample,nSNP,ntrait,IM);
	load(outfile);
    elseif ~isempty(matfile)
        load(matfile)
    else
        traw2mat(phefile,genofile,outfile,sample,nSNP,ntrait,IM);
        load(outfile);
    end
    %traw2mat('../../pop.fam','../../poptraw.traw','popnew.mat',175,407045,3)
    
    Y=y;X=x;
    [~,p1]=size(Y);
    if method=='randomly'
        %clear,clc
        warning('off','all')
        diary('ISR.log'); % notes the command window diary
        yname=2;%no change the name of phenotype
        impute=2;
        alfa=0.001;
        glm=3;%the number of running default 3 
        ept0=25;%initially random chosed the number of multi-locus (probable result) default 5
        gr=2;%export the result of significant association SNP genotype (with gr=1)
        nchr=chr;%the number of Chromosome
        sgt=sgv;%  bonferroni correction
        %y=Y(:,p1);
        N = size(x,1);
        sq=randperm(N);testing=sq(1:round(N*0.2));trainning=setdiff(sq,testing);
        %for i=1:nfolds
        %which=foldid==i;
        %if verbous, disp(['Fitting fold # ' num2str(i) ' of ' num2str(nfolds)]);end
        % x=reperma
        x=X(trainning,:);y=Y(trainning,1);xte=X(testing,:);xtr=X(trainning,:);
        %find the MISSING DATA add to the program
        % [n,p]=size(x);
        % for i=1:n
         %  for j=1:p
         %      if isnan(x(i,j))
        %          disp([i,j,x(i,j)])
        %      end
        %   end
        % end
        %stop
        % ISR method
        if mdl==2
            if ~isempty(ncov)
                ISR_Epi_COV
            else
                ISR_Epi
            end
        else
            if ~isempty(ncov)
                ISR_COV
            else
                ISR_REG
            end
        end     
       %changing possible the last association results
        %writetable(opt_result,opt_outresult,'Delimiter','\t');writetable(all_result,all_outresult,'Delimiter','\t');
        diary off
        
    elseif method=='fold'
        impute=2;yname=2;
        diary(['isr',num2str(i),'.log']); % notes the command window diary
        alfa=0.001;
        glm=3;%the number of running 
        ept0=25;%initially random chosed the number of multi-locus (probable result)
        gr=2;%export the result of significant association SNP genotype (with gr=1)
        nchr=chr;%the number of Chromosome
        sgt=sgv; 
        y=Y(:,p1);
            %find the MISSING DATA add to the program
            % [n,p]=size(x);
            % for i=1:n
            %  for j=1:p
            %      if isnan(x(i,j))
            %          disp([i,j,x(i,j)])
            %      end
            %   end
            % end
            %stop
            % ISR method
        if mdl==2
             if ~isempty(ncov)
                ISR_Epi_COV
            else
                ISR_Epi
            end
        else
            if ~isempty(ncov)
                ISR_COV
            else
                ISR_REG
           end
        end
        %changing possible the last association results
        writetable(opt_result,opt_outresult,'Delimiter','\t');writetable(all_result,all_outresult,'Delimiter','\t');
        diary off 
    end
end
   
