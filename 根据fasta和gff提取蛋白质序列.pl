#!/usr/bin/perl

#打开fna文件。把对应染色体的序列放入哈希表。key是染色体名，value是碱基序列。
open(SEQ,"<GCF_000001405.39_GRCh38.p13_genomic.fna") or die "file.txt 文件无法打开，$!";
$ori=$/; 
$/=">"; #将分隔符改成>,等会就不是每次读一行了，而是每次读的内容以>分割，即每次读一个染色体

%chromNameToSeq = (); #初始化哈希表，用于存放染色体名和染色体序列的对应关系
while(<SEQ>){
	next if (/^>$/); #就是如果每次以>分隔来读入的话，那么每次读入内容的结尾就是>。由于该文件第一个字符就是>,所以第一次读入的内容就是一个>,所以需要跳过第一次。
	/(.*?)\s.*?\n(.*)/s; #读入的内容是放到变量$_中，不写对哪个匹配就默认对它匹配。该正则表达式的结尾修饰符s是指让.可以匹配换行\n，默认是不能匹配的。
	#注意正则匹配时，一个单位匹配到哪里，取决于它的下一个单位是啥以及是否贪婪。加上问号就代表不贪婪：即只要遇到它的下一个单位就停止，否则是遇到最远的下一个单位才停止
	#那么上面是先匹配任意直到遇到空格（非贪婪），那么可见匹配的就是染色体的名字；然后再匹配任意直到换行（非贪婪），就是匹配染色体名字后面的那串解释，匹配完这一整行；
	#接着是匹配任意字符（贪婪），即匹配后面所有的碱基序列（注意其中混有换行符了，等会需要去掉）
	$nameOfChromosome = $1;
	$seqOfChromsome = $2;
	$seqOfChromsome =~ s/\n//g;  #这里就是把所有的换行符\n替换成空。结尾g代表全局替换，即把所有的换行符都替换。
	$seqOfChromsome =~ tr/>//;   #这是把结尾的>删去。
	$seqOfChromsome = "\U$seqOfChromsome"; #在字符串前加上\U 即代表全转换成大写
	$chromNameToSeq{$nameOfChromosome} = $seqOfChromsome; #然后是把对应的染色体名字作为哈希表的key，对应的序列作为value（这哈希表也可以不用初始化直接用，狗得很）
}
close (SEQ);
$/=$ori; #恢复默认分隔符的值

#将标准密码子表存入哈希表condonTable1
%condonTable1=('TCA' => 'S','TCC' => 'S','TCG' => 'S','TCT' => 'S','TTC' => 'F',
                'TTT' => 'F','TTA' => 'L','TTG' => 'L','TAC' => 'Y','TAT' => 'Y',
                'TAA' => '*','TAG' => '*','TGC' => 'C','TGT' => 'C','TGA' => '*',
                'TGG' => 'W','CTA' => 'L','CTC' => 'L','CTG' => 'L','CTT' => 'L',
                'CCA' => 'P','CCC' => 'P','CCG' => 'P','CCT' => 'P','CAC' => 'H',
                'CAT' => 'H','CAA' => 'Q','CAG' => 'Q','CGA' => 'R','CGC' => 'R',
                'CGG' => 'R','CGT' => 'R','ATA' => 'I','ATC' => 'I','ATT' => 'I',
                'ATG' => 'M','ACA' => 'T','ACC' => 'T','ACG' => 'T','ACT' => 'T',
                'AAC' => 'N','AAT' => 'N','AAA' => 'K','AAG' => 'K','AGC' => 'S',
                'AGT' => 'S','AGA' => 'R','AGG' => 'R','GTA' => 'V','GTC' => 'V',
                'GTG' => 'V','GTT' => 'V','GCA' => 'A','GCC' => 'A','GCG' => 'A',
                'GCT' => 'A','GAC' => 'D','GAT' => 'D','GAA' => 'E','GAG' => 'E',
                'GGA' => 'G','GGC' => 'G','GGG' => 'G','GGT' => 'G');

#将第二套密码子表存入哈希表condonTable2
%condonTable2=('TCA' => 'S','TCC' => 'S','TCG' => 'S','TCT' => 'S','TTC' => 'F',
                'TTT' => 'F','TTA' => 'L','TTG' => 'L','TAC' => 'Y','TAT' => 'Y',
                'TAA' => '*','TAG' => '*','TGC' => 'C','TGT' => 'C','TGA' => 'W',#由*变为W
                'TGG' => 'W','CTA' => 'L','CTC' => 'L','CTG' => 'L','CTT' => 'L',
                'CCA' => 'P','CCC' => 'P','CCG' => 'P','CCT' => 'P','CAC' => 'H',
                'CAT' => 'H','CAA' => 'Q','CAG' => 'Q','CGA' => 'R','CGC' => 'R',
                'CGG' => 'R','CGT' => 'R','ATA' => 'M',#由I变为M
                'ATC' => 'I','ATT' => 'I',
                'ATG' => 'M','ACA' => 'T','ACC' => 'T','ACG' => 'T','ACT' => 'T',
                'AAC' => 'N','AAT' => 'N','AAA' => 'K','AAG' => 'K','AGC' => 'S',
                'AGT' => 'S','AGA' => '*',#由R变为*
                'AGG' => '*',#由R变为*
                'GTA' => 'V','GTC' => 'V',
                'GTG' => 'V','GTT' => 'V','GCA' => 'A','GCC' => 'A','GCG' => 'A',
                'GCT' => 'A','GAC' => 'D','GAT' => 'D','GAA' => 'E','GAG' => 'E',
                'GGA' => 'G','GGC' => 'G','GGG' => 'G','GGT' => 'G');

#打开gff文件。匹配到CDS的行，然后根据起止位置、正负链、CDS的id（相同id的CDS拼到一起），从刚才那个染色体序列哈希表中去查找对应的序列，然后把CDS名和序列存入哈希表。
open(ANNO,"<GCF_000001405.39_GRCh38.p13_genomic.gff") or die "file.txt 文件无法打开，$!";
%useWhichTable=(); #初始化一个哈希表，用于记录某个CDS名是用哪一个密码子表
%cdsNameToSeq=();  #初始化一个哈希表，用于记录CDS名与碱基序列的对应关系

while(<ANNO>){
	#chomp; #每次读一行（默认按\n分隔，所以结尾是\n, 因此需要用chomp去掉。省略参数就是默认对$_操作）
    if(/ID=cds-(.*?);/){  #如果匹配到cds的话，说明这一行是cds。小括号里就匹配-之后的任意个字符直到分号，（非贪婪匹配），因此匹配的就是cds名。
        my @line=split("\t",$_); #将匹配到的这一行用制表符分割成数组
        $useWhichTable{$1}=index($_,"transl_table=2"); #$1就是匹配到的CDS名，看它对应的是哪个密码子表，并记录。
		#上面那个index是查找含有"transl_table=2"的下标，找不到返回-1就说明是第一个表，
        if($line[6] eq "-"){
			$aPartOfSeq=reverse(substr($chromNameToSeq{$line[0]},$line[3]-1,$line[4]-$line[3]+1)); #substr就是从某个数组中，提取从第几号位置开始的n个字符。
			$aPartOfSeq=~tr/TCGA/AGCT/; #正则替换。
		}
        else{
          $aPartOfSeq=substr($chromNameToSeq{$line[0]},$line[3]-1,$line[4]-$line[3]+1);
        }
        $cdsNameToSeq{$1}.=$aPartOfSeq; #将提取到的序列拼过来。（perl中的.=就相当于字符串的+=，这也太妖了吧。。。）
    }
}
close (ANNO);

#将CDS序列存入文件cdsNameToSeq.txt
open(IDSEQ,">cdsNameToSeq.txt"); #>代表写入的意思
foreach $key (keys %cdsNameToSeq){    #keys 加哈希表 再用括号括起来，外面foreach，这种语法到时候再记一下。太weird。
	print IDSEQ ">$key\n$cdsNameToSeq{$key}\n";  #print + 文件句柄 + 内容 就是把内容输入到文件。（这种写法狗得很，让我一个java过来的很是不爽）
}
close(IDSEQ);

#将CDS序列翻译为蛋白质序列存入哈希%cdsToProtenSeq
%cdsToProtenSeq=();#初始化哈希表，存放CDS名对应的蛋白质序列

foreach $key (keys %cdsNameToSeq){ #对于每个CDS序列，在之后再弄一层循环，一个个把碱基翻译成蛋白质，然后拼在一起
    $protein=""; #初始蛋白质是空，之后一个个拼过来
    if($useWhichTable{$key}==-1){ #如果用的是第一个密码子表
		for(my $a=0;$a<length($cdsNameToSeq{$key})-2;$a+=3){    #沿着碱基序列三个三个往后走
        $protein.=$condonTable1{substr($cdsNameToSeq{$key},$a,3)};  #substr就是对于某一段字符串，把它第n个位置的3个碱基提出来
        }
    }
    else{
		for(my $a=0;$a<length($cdsNameToSeq{$key})-2;$a+=3){
        $protein.=$condonTable2{substr($cdsNameToSeq{$key},$a,3)};
        }
    }
    $cdsToProtenSeq{$key}=$protein; 
    print  ">$key\n$cdsToProtenSeq{$key}\n";  #把CDS名和对应的蛋白质序列输出看一下。注意print后跟双引号，默认结尾没有换行。所以需要手动加上一个换行。
}

#将蛋白序列存入文件cdsToProtenSeq.txt
open(IDPRO,">cdsToProtenSeq.txt");
foreach $key (keys %cdsToProtenSeq){
    print IDPRO ">$key\n$cdsToProtenSeq{$key}\n";
}
close IDPRO;
