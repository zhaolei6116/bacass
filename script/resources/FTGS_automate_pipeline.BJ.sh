#!/bin/bash
echo "#########################################################################"
echo "FTGS自动化脚本:FTGS_automate_pipeline.sh"
echo "注意该脚本10分钟扫描一次下机样本信息单"
echo "自动化脚本路径: /bioinformation/Project/FTGS_sh/FTGS_automate_pipeline.sh"
#########################################################################
current_path=$(dirname "${BASH_SOURCE[0]}")
echo $current_path
DATE=$(date +%Y%m%d)
second=600                             #监控时间间隔单位为10分钟
#=========================main===========================================
function fn_showlog()
{
	local curtime;
	curtime=`date +"%Y%m%d-%H:%M:%S"`
	echo "$curtime ------ $1";
}
while true
do
if [ ! -f sampleinfo.end ];then
	export start=`date -d "240 minute ago" "+%Y-%m-%d %H:%M:%S"`
	export end=`date "+%Y-%m-%d %H:%M:%S"`
	echo "Start time:$start  End Time:$end"
	################### Get sample xls and creat conf #########################
	#rm -rf /bioinformation/Project/FTGS_bj/LimsData/*
	#/data/software/miniconda/bin/java -jar /bioinformation/Project/FTGS_bj/script/CwbioRequestDataLims.V2.jar --startTime="$start" --endTime="$end" --config="/bioinformation/Project/FTGS_bj/script/config/config.properties.BJ"
	#/data/software/miniconda/bin/java -jar /bioinformation/Project/FTGS_bj/script/CwbioRequestDataLims.V2.jar --startTime="$start" --endTime="$end" --config="/bioinformation/Project/FTGS_bj/script/config/config.properties.GZ"
	#/data/zwx/miniconda_R/bin/Rscript /bioinformation/Project/FTGS_bj/script/xlsx_transfer_conf.R "/bioinformation/Project/FTGS_bj/LimsData/*/"
	###################### analysis  TNS Sample #########################
	/data/software/miniconda/bin/perl /nas04/Project/FTGS.automate.pipe/script/Monitor.FTGS.pl --dir /nas04/Project/FTGS/conf/ --suf conf --stat /nas04/Project/FTGS/Stat.conf --cmd "/data/software/miniconda/bin/perl /nas02/project/wangkx/pipe_FTGS_update/FTGSpipe_v3.12/00.FTGSPipeline/TGSpipelineV4.pl --xml PCR_ref=/nas02/project/wangkx/pipe_FTGS_update/FTGSpipe_v3.12/XML/FTGSpipe.pcr.REF.v5.0.xml --xml PCR_noref=/nas02/project/wangkx/pipe_FTGS_update/FTGSpipe_v3.12/XML/FTGSpipe.pcr.NOREF.v4.0.xml --xml  PLA_ref=/nas02/project/wangkx/pipe_FTGS_update/FTGSpipe_v3.12/XML/FTGSpipe.plasmid.REF.v6.0.xml --xml PLA_noref=/nas02/project/wangkx/pipe_FTGS_update/FTGSpipe_v3.12/XML/FTGSpipe.plasmid.NOREF.v6.0.xml --xml Hap_noref=/nas02/project/wangkx/pipe_FTGS_update/FTGSpipe_v3.12/XML/FTGSpipe.Haplotype.NOREF.v5.0.xml --xml Only_seq=/nas02/project/wangkx/pipe_FTGS_update/FTGSpipe_v3.12/XML/FTGSpipe.OnlySeq.v2.0.xml --cfg \$CfgFile --cts /nas04/Project/FTGS/ContactsList.xls --out /nas04/Project/FTGS/\$OutDate --sh --inexe && /data/software/miniconda/bin/perl /nas02/project/wangkx/pipe_FTGS_update/FTGSpipe_v3.12/00.FTGSPipeline/distribution.v3.1.pl /nas04/Project/FTGS/\$OutDate/Tasks/\$Basename/bin/\$Basename.sh /nas04/Project/FTGS/\$OutDate/Tasks/\$Basename" --log /nas04/Project/FTGS/Monitor.log 
else
	echo `date` "数据循环结束，请查找原因" >> ctDNA_auto.log
	exit 1
fi
sleep $second
done

