## ExMODE -- Explore extreme environments and protect the code of life on earth
数据分析流程
### biome_list
五个生境样本列表
### metadata.xlsx
所有样本元信息表，包括BioSample、Taxonomy name、Isolation Source、Sample Title、Sample description、Project name、Project description等信息。
### sample.txt
所有样本列表
### pipeline.sh 
#### 1. Fastp质控
#### 2. Metaphlan4物种组成分析
#### 3. Megahit宏基因组组装
#### 4. Metagenemark基因预测
#### 5. MMseq2构建非冗余基因集
#### 6. 基因功能注释
#### 7. ESMfold蛋白结构预测
#### 8. MetaWRAP宏基因组分箱
#### 9. 构建非冗余基因组集、物种注释
#### 10. 生物合成基因簇（BGCs）

