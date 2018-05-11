#########################################################################
# File Name: work.sh
# Author: jdw
# mail: jiangdawei@icarbonx.com
# Created Time: Thu 01 Feb 2018 04:27:03 PM CST
#########################################################################
#!/bin/bash
#split rawdata using barcode_Spacer.txt file
#note: There is a sequence of varied length  after the 6bp barcode.
#echo "Start : split raw data"`date`
#../SplitPip  ./raw_data/test_data.fq1.gz  ./raw_data/test_data.fq2.gz   ./sample.txt  ./split_result
#echo "End : split data"`date`

#echo "Start : usearch pipline"`date`
#../AnalysisPip  \
#	./sample.list  \
#	./split_result  \
#	./analysis_result
#echo "End :usearch "`date`

#echo "Start : vsearch pipline"`date`
#../AnalysisPips  \
#	./sample.list  \
#	./split_result  \
#	./analysis_results
#echo "End :vsearch "`date`
