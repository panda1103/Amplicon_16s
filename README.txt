#
1. SplitPip
Usage: ./SplitPip [Read1_Fastq.gz] [Read2_Fastq.gz] [Barcode to Sample list] [Outdir] [minratio:default 1.0]

脚本SplitPip用于QC和样品数据拆分:

QC: 过滤掉Q20比例<70%的pair reads,具体阈值可以在SplitPip里进行修改;

数据拆分:
1)依据barcode_Spacer.txt 和 Sample.list进行数据拆分,注意该'barcode_Spacer.txt'的barcode为6bp,
而barcode长度与primer match位置是对应的,所以使用其他长度的barcode需要对barcode.py进行一定的修改;
2)barcode后,设计了spacer region区域,用于增加序列的多样性,如果你实验没有在barcode后设计该序列,
请使用barcode_noSpacer.txt格式的文件替换barcode_Spacer.txt文件,然后分析;
3)barcode与序列进行匹配时,默认是的min match ratio为1.0, 即不允许mismatch;
建议minratio不低于0.81

2. AnalysisPip
Usage: ./AnalysisPip [Sample list] [Fastq Indir] [Outdir]

1)脚本AnalysisPip用于16s标准分析,该pipline主要参照usearch进行分析,丰度< 0.5%的otu会被过滤掉;
但免费版usearch在数据量比较大时无法运行,此时建议使用AnalysisPip2进行分析
2)AnalysisPip2 相较于AnalysisPip,其主要使用vsearch进行分析,但两个软件在质量过滤,嵌合体过滤的结果各不相同;
AnalysisPip2会过滤掉reads数目<2的otu


