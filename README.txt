######################
##  16s 扩增子分析流程
######################

1. 下机数据处理
使用脚本: SplitPip

Usage: 
    ./SplitPip [Read1_Fastq.gz] [Read2_Fastq.gz] [Barcode to Sample list] [Outdir] [minratio:default 1.0]

参数说明:
    "Read1_Fastq.gz, Read2_Fastq.gz"    测序的pair reads, gz压缩格式;
    "Barcode to Sample list"            barcode与样品的对应列表;
    "Outdir"                            结果的输出路径;
    "minratio:default 1.0"              barcode 与reads匹配时,需满足的最低match ratio.

分析流程说明:
    脚本SplitPip,是将barcode.py进行包装,用于QC和样品数据拆分,其中
   
    barcode.py 参数包括:

    --fq1                       FQ1
    --fq2                       FQ2
    --barcodefile               BARCODEFILE
    --primerfile                PRIMERFILE
    --samplefile                SAMPLEFILE
    --minratio                  MINRATIO
    --low_quality_threshold     LOW_QUALITY_THRESHOLD
    --high_quality_percent      HIGH_QUALITY_PERCENT
    --outpath                   OUTPATH

    barcode.py处理数据的具体步骤如下:

    1)QC: 
        a. 统计reads的碱基质量,在pair reads中,如果Q20的比例低于70%, 则过滤掉,该阈值由参数 "high_quality_percent" 确定; 
        b. 统计reads末端质量值低于20且连续的碱基序列,截去这些序列,该阈值由参数"low_quality_threshold" 确定;

    2)primer匹配:
        a. 根据文件"primerfile.txt"中的引物,在reads 中确定引物的起始和终止位置, 并将引物序列从reads中截去;

    3)barcode匹配:
        a. 根据引物的起始位置,将该位置之前的序列作为与barcode匹配序列;
        b. 在barcode序列之后,设计了spacer region区域,该序列与barcode序列一起作为匹配依据,spacer 与barcode的对应关系由参数"barcodefile"输入;
        c. 序列的match ratio由"minratio"确定,默认为1.0, 建议不低于0.81;

2. 数据分析流程
使用脚本: AnalysisPip

Usage: 
    ./AnalysisPip [Sample list] [Fastq Indir] [Outdir]

参数说明:
    "Sample list"        需要分析的样品ID列表
    "Fastq Indir"        输入的Fastq文件路径
    "Outdir"             结果的输出路径

分析流程说明:
    脚本AnalysisPip, 是用于16s扩增子分析的pipline,主要参照 usearch10 官方流程(http://www.drive5.com/usearch/manual/);
    该pipline 由以下4个小流程包装而成:

         1) ./MergeFq [Sample list] [Faste Indir] [Merged Fasetq Outdir]
         2) ./UsearchPip [fastq dir] [Outdir]
         3) ./plot.py [-h] [--path PATH]
         4) ./pcoa.py [-h] [--distance_file DISTANCE_FILE] [--distance_type DISTANCE_TYPE] [--group_file GROUP_FILE] \
                         [--group_lable GROUP_LABLE] [--outfile_prefix OUTFILE_PREFIX]
    
    该pipline具体步骤如下:
    
    1)数据合并:根据overlap合并pair reads;
    2)数据质量过滤,并将fastq格式文件转化为fasta格式;
    3)去除冗余序列,保留unique序列;
    4)OTU聚类
    5)去除嵌合体序列
    6)生成OTU table
    7)过滤低丰度OTU
    8)计算alpha多样性,beta多样性,稀释曲线
    9)使用数据库"RDP training set v16" 进行种属鉴定(http://www.drive5.com/usearch/manual/sintax_downloads.html)
    10)计算不同分类水平下,各个门类的丰度
    11)结果可视化

3.其他
    1)免费版usearch在数据量较大时无法运行,此时建议使用AnalysisPip2进行分析
    2)AnalysisPip2 相较于AnalysisPip,其中内存受限的步骤使用vsearch将usearch替换,两个软件的使用方法大致相同;
    3)vsearch 与 usearch在嵌合体的过滤原理不同,结果有差异,不能确定哪个更为准确;
    4)AnalysisPip2会过滤掉reads数目<2的otu, 而AnalysisPip过滤标准为丰度低于0.001, 具体请根据需要做修改;


