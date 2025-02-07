# 细菌完成图流程

## Introduction

细菌完成图流程，用于细菌基因组组装。

## Usage

### 1、项目信息来源
项目信息通过产品线上游同事邮件通知或者信息通知，有样本需要进行分析。 

### 2、项目样本信息获取
通过以下脚本获取样本的信息单， 脚本与质粒商检流程中一样使用IT给的JAVA程序。 

```bash
sh script/pull_infor_form.sh
```

信息单会在路径中LimsData 文件夹中。
获取到信息单后，需要筛选出项目类型(Project)是细菌完成图的样本，进行后续的分析。 
分析后需要手动清除LimsData 文件夹中的信息。

### 3、检查修改下机单数据
常见问题
1. 手动修改下机单文件的下机路径, 路径修改精确到 `fastq_pass` 如下
`/BSequenator03/25011311/no_sample/20250117_1736_P2S-02106-B_PAY44743_3eb272e8/fastq_pass`
2. 检查 `Genome_size`
要求如 1.5M  或者 1500000， 客户没写单位 或者 单位不是M 都要改写。
下面为示例，`Genome_size`,  客户填写的是1790， 询问前端，是1790K，修改后如下 1.97M。

3. 检查Species_name，两侧不能有引号 ，如 "Coxiella Burnetii" , 引号要去掉。

```bash
# 下机单示例
Project           Sample_name   Barcode     Barcode_type   Sample_type   Ref   Client           Length   PLASMID_LENGTH   Demand   Demand_150   Haplotype   Path                                                                                       Detect_no      Detect_cnt   Report_path                                      Species_name        Genome_size   Data_volume
---------------   -----------   ---------   ------------   -----------   ---   --------------   ------   --------------   ------   ----------   ---------   ----------------------------------------------------------------------------------------   ------------   ----------   ----------------------------------------------   -----------------   -----------   -----------
   细菌完成图        D3-LCV        barcode15   SQK            dna           N     SD250120145031   0        0                D        N            N           /BSequenator04/25012212/no_sample/20250122_1733_P2S-02167-A_PAY44458_864cf8e3/fastq_pass   B22501200998   1            /gwservice/b/primecx/25012212/B22501200998.zip   Coxiella Burnetii   1.97M         -
   细菌完成图        D10-SCV       barcode14   SQK            dna           N     SD250120145031   0        0                D        N            N           /BSequenator04/25012212/no_sample/20250122_1733_P2S-02167-A_PAY44458_864cf8e3/fastq_pass   B22501200999   1            /gwservice/b/primecx/25012212/B22501200999.zip   Coxiella Burnetii   1.97M         -
```

### 4、项目分析
建立项目分析文件夹，以样本ID （即在信息单中 Detect_no 列的信息）作为文件夹的名字。


```bash
# 举例
B22412250031.xls

mkdir B22412250031
cd B22412250031
cp ../LimsData/path/B22412250031.xls. ./

# 
#创建 run.sh 内容如下进行分析 --csv B22412250031.xls 这个参数为输入， 输入文件为下机单文件。
source /nas02/software/conda/Miniconda3/miniconda3/bin/activate /nas02/project/huyifan/software/nextflow/v24.04.4
nextflow run /nas02/project/zhaolei/pipeline/bacteria_genome_assembly/main.nf --csv B22412250031.xls  -resume -profile sge

# 运行run.sh 即可
sh run.sh
或者
nohup sh run.sh & 

任务在内部自动qsub投递
```


### 5、释放数据
在结果文件路径： ./analysis/sample_id/  下有release.sh 
该脚本目前是只生成，不执行。
检查结果没问题，执行该脚本进行数据释放。


### 6、参考路径
/nas02/project/zhaolei/project/bacteria_genome_assembly/B22501200913
结果文件路径： ./analysis/B22501200913/
执行脚本：    run.sh
输入文件：    input.xls

