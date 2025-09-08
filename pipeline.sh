## 1. Fastp质控
	for i in `cat sample.txt`;do
		mkdir temp/${i}/fastp
		fastp -i seq/${i}/${i}_1.fastq.gz -o temp/${i}/fastp/${i}_clean_1.fq.gz \
		-I seq/${i}/${i}_2.fastq.gz -O temp/${i}/fastp/${i}_clean_2.fq.gz \
		--json temp/${i}/fastp/${i}_fastp.json --html temp/${i}/fastp/${i}_fastp.html
	done
		
## 2. Metaphlan4物种组成分析
	conda activate metaphlan4
	for i in `cat sample.txt`;do
		mkdir temp/${i}/mp4 
		metaphlan temp/${i}/fastp/${i}_clean_1.fq.gz,temp/${i}/fastp/${i}_clean_2.fq.gz \
		--input_type fastq -o temp/${i}/mp4/${i}.txt --nproc 10 --bowtie2out temp/${i}/mp4/${i}/bowtie2out_${i}.bz2
	done
	
## 2.1 物种组成表
	for i in `cat sample.txt`;do
		merge_metaphlan_tables.py temp/${i}/mp4/${i}.txt > temp/${i}/mp4/merged_abundance_table.txt
		grep -E '(s__)|(clade_name)' temp/${i}/mp4/merged_abundance_table.txt |grep -v 't__'|sed 's/^.*s__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' \
		> temp/${i}/mp4/merged_abundance_species.txt
		grep -E '(g__)|(clade_name)' temp/${i}/mp4/merged_abundance_table.txt |grep -v 's__'|sed 's/^.*g__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' \ 
		> temp/${i}/mp4/merged_abundance_genus.txt
		grep -E '(f__)|(clade_name)' temp/${i}/mp4/merged_abundance_table.txt |grep -v 'g__'|sed 's/^.*f__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' \
		> temp/${i}/mp4/merged_abundance_family.txt
		grep -E '(o__)|(clade_name)' temp/${i}/mp4/merged_abundance_table.txt |grep -v 'f__'|sed 's/^.*o__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' \
		> temp/${i}/mp4/merged_abundance_order.txt
		grep -E '(c__)|(clade_name)' temp/${i}/mp4/merged_abundance_table.txt |grep -v 'o__'|sed 's/^.*c__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' \
		> temp/${i}/mp4/merged_abundance_class.txt
		grep -E '(p__)|(clade_name)' temp/${i}/mp4/merged_abundance_table.txt |grep -v 'c__'|sed 's/^.*p__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' \
		> temp/${i}/mp4/merged_abundance_phylum.txt
	done
		
## 3. Megahit组装Assembly
	conda activate megahit
	for i in `cat sample.txt`;do
		megahit -1 temp/${i}/fastp/${i}_clean_1.fq.gz -2 temp/${i}/fastp/${i}_clean_2.fq.gz \
		-o temp/${i}/assembly --out-prefix ${i} --presets meta-sensitive -t 10
	done

## 4. Metagenemark基因预测Gene prediction
	for i in `cat sample.txt`;do
		mkdir temp/${i}/metagenemark
		gmhmmp -m soft/metagenemark/MetaGeneMark_v1.mod \
		-o temp/${i}/metagenemark/${i}.gmm temp/${i}/assembly/${i}.contigs.fa -A temp/${i}/metagenemark/${i}_prot.fa
	done

## 5. Mmseq2构建非冗余基因集

## 5.1 分生境构建基因集
	mkdir result/geneset
	for j in {1..5};do
		mkdir result/geneset/EB00000${j}
		for i in $(cat "biome_list/EB00000${j}.txt");do
			cat temp/${i}/metagenemark/${i}_prot.fa >> result/geneset/EB00000${j}/EB00000${j}_all_prot.fa
		done
	done
	
## 5.2 Mmseq2分生境构建非冗余基因集
	conda activate mmseq2
	for j in {1..5};do
		mmseqs easy-cluster result/geneset/EB00000${j}/EB00000${j}_all_prot.fa \
		result/geneset/EB00000${j}/EB00000${j}_cluster result/geneset/EB00000${j}/EB00000${j}_cluster_tmp \
		--min-seq-id 0.95 --cov-mode 1 --cluster-mode 3 --threads 35
	done
	
	# 修改序列ID
	perl -e 'open OU1,">result/geneset/EB000001/Final_TM_geneset_prot.fa"; \
	open OU2,">result/geneset/EB000001/name.change.list"; \
	while(<>){ chomp; if(/^>/){ $id=sprintf("%010d",++$id); print OU2 "$_\tTM_$id\n"; print OU1 ">TM_$id\n"; next; } print OU1 "$_\n"; }' \
	result/geneset/EB000001/EB000001_cluster_rep_seq.fasta
	
## 5.3 分生境构建蛋白结构预测集
	conda activate mmseq2
	for j in {1..5};do
		mmseqs easy-cluster result/geneset/EB00000${j}/Final_TM_geneset_prot.fa \
		result/geneset/EB00000${j}/EB00000${j}_cluster2 result/geneset/EB00000${j}/EB00000${j}_cluster2_tmp \
		--min-seq-id 0.20 -c 0.50 --cov-mode 1 --cluster-mode 3 --threads 35
		awk '{print $1}' result/geneset/EB00000${j}/EB00000${j}_cluster2_cluster.tsv | sort | uniq -c | awk '$1>=30' \
		> result/geneset/EB00000${j}/EB00000${j}_cluster2.txt
	done
	seqkit grep -f result/geneset/EB000001/EB000001_cluster2.txt \
	result/geneset/EB000001/Final_TM_geneset_prot.fa -o result/geneset/EB000001/Final_TM_geneset_prot2.fa

## 6. 基因注释Gene annotation
	for j in {1..5};do
		mkdir result/geneset/EB00000${j}/eggnog
		mkdir result/geneset/EB00000${j}/cazy
		mkdir result/geneset/EB00000${j}/card
	done
	
## 6.1 EggNOG注释
	conda activate eggnog
	emapper.py -i result/geneset/EB000001/Final_TM_geneset_prot.fa \
	-o result/geneset/EB000001/eggnog/EB000001 --cpu 20 \
	--data_dir soft/envs/eggnog/lib/python3.12/site-packages/data --dbmem
	
## 6.2 Cazy注释
	diamond blastp --db soft/cazydb/cazy.dmnd --query result/geneset/EB000001/Final_TM_geneset_prot.fa \
	-e 0.00001 --outfmt 6 --more-sensitive --max-target-seqs 1 --threads 22 --quiet \
	--out result/geneset/EB000001/cazy/EB000001_cazy

## 6.3 Card注释
	conda activate rgi
	rgi main --input_sequence result/geneset/EB000001/Final_TM_geneset_prot.fa \
	--output_file result/geneset/EB000001/card/EB000001_card --local --clean --include_loose -t protein

## 7. ESMfold蛋白结构预测
	mkdir result/structureset
	for j in {1..5};do
		mkdir result/structureset/EB00000${j}
	done
	
	conda activate esmfold
	esm-fold -i result/geneset/EB000001/Final_TM_geneset_prot2.fa -o result/structureset/EB000001/EB000001_structureset

## 8. MetaWRAP样本分箱Binning
	conda activate metawrap
	for i in `cat sample.txt`;do
		mkdir temp/${i}/metawrap
		gzip -dc temp/${i}/fastp/${i}_clean_1.fq.gz > temp/${i}/fastp/${i}_clean_1.fq
		gzip -dc temp/${i}/fastp/${i}_clean_2.fq.gz > temp/${i}/fastp/${i}_clean_2.fq

		metawrap binning -o temp/${i}/metawrap/${i} -t 10 --metabat2 -l 1000 -a \
		temp/${i}/assembly/${i}.contigs.fa temp/${i}/fastp/${i}_clean_1.fq temp/${i}/fastp/${i}_clean_2.fq
	
		checkm lineage_wf -t 10 -x fa temp/${i}/metawrap/metabat2_bins temp/${i}/metawrap/CheckM --pplacer_threads 10 && \
		python soft/metaWRAP/bin/metawrap-scripts/summarize_checkm.py \
		temp/${i}/metawrap/CheckM/storage/bin_stats_ext.tsv ${i}|(read -r;printf %sn ; sort) > \
		temp/${i}/metawrap/CheckM/${i}_metaBat2.stats.xls
	
		awk '$2>50&&$3<10' temp/${i}/metawrap/CheckM/${i}_metaBat2.stats.xls > temp/${i}/metawrap/CheckM/${i}_metaBat2_50_10.stats.xls
		mkdir temp/${i}/metawrap/CheckM/metabat2_bins_50_10
		for k in $(awk 'NR>1 {print $1}' "temp/${i}/metawrap/CheckM/${i}_metaBat2_50_10.stats.xls");do
			cp temp/${i}/metawrap/metabat2_bins/${k}.fa temp/${i}/metawrap/CheckM/metabat2_bins_50_10/${i}_${k}.fa
		done
	done

## 9. 构建基因组集
	## 9.1 分生境构建基因组集
	for j in {1..5};do
		mkdir result/genomeset/EB00000${j}
		for i in $(cat "biome_list/EB00000${j}.txt");do
			cp temp/${i}/metawrap/CheckM/metabat2_bins_50_10/*.fa result/genomeset/EB00000${j}/EB00000${j}_all_genome
		done
		ls result/genomeset/EB00000${j}/EB00000${j}_all_genome > result/genomeset/EB00000${j}/EB00000${j}_all_genome.txt
	done
	
	## 9.2 GTDB-tk物种注释
	conda activate gtdbtk
	for j in {1..5};do
		gtdbtk classify_wf --skip_ani_screen --genome_dir result/genomeset/EB00000${j}/EB00000${j}_all_genome \
		--out_dir result/genomeset/EB00000${j}/gtdbtk --cpus 20 -x fa
	done
	
	## 9.3 分生境构建非冗余基因集dRep
	conda activate drep
	for j in {1..5};do
		dRep dereplicate result/genomeset/EB00000${j}/EB00000${j}_drep_genome -g result/genomeset/EB00000${j}/EB00000${j}_genome_path.txt \
		-p 50 --length 0 --completeness 50 --contamination 10 --S_algorithm ANImf --P_ani 0.9 --S_ani 0.95 --cov_thresh 0.3 --strain_heterogeneity_weight 0 \
		--genomeInfo result/genomeset/EB00000${j}/EB00000${j}_checkm.csv 
	done

## 10. 生物合成基因簇（BGCs）
	## 10.1 Antismash
	conda activate antismash
	for j in {1..5};do
		for i in $(cat "result/genomeset/EB00000${j}/EB00000${j}_all_genome.txt");do
			antismash --minimal --genefinding-tool prodigal --cpus 4 \
			--output-dir result/bgcset/EB00000${j}/antismash/result/${i} \
			result/genomeset/EB00000${j}/EB00000${j}_all_genome/${i} --cb-knownclusters
		done
	done
	
	for j in {1..5};do
		for i in $(cat "result/genomeset/EB00000${j}/EB00000${j}_all_genome.txt");do
			cp result/bgcset/EB00000${j}/antismash/result/${i}/*.region*.gbk result/bgcset/EB00000${j}/antismash/EB00000{j}_allgbks
		done
	done
	
	## 10.2 Bigscape
	conda activate bigscape
	for j in {1..5};do
		python soft/bigscape/BiG-SCAPE-1.1.5/bigscape.py -i result/bgcset/EB00000${j}/antismash/EB00000{j}_allgbks \
		-o result/bgcset/EB00000${j}/bigscape_result \
		--migbig --mix --clan_cutoff 0.3 0.7 --hybrids-off -c 20 --mode auto --pfam_dir soft/pfam_dir
	done

