#处理fasta文件

import linecache
import re
import sys

a= linecache.getlines('GCF_000001405.39_GRCh38.p13_genomic.fna') #返回一个字符串数组
dict={} #承装染色体和序列的映射关系

theStartPosition=[] #用于记录所有染色体起始行的行号
indexToName=[] #完成theStartPosition下标（即染色体序号）到染色体名字的映射


for i in range(0,len(a)):
	if re.match(r'>',a[i]):
		theStartPosition.append(i)
		indexToName.append(re.split('>| ',a[i])[1]) #把染色体的名字添加到那个映射数组中
		


for i in range(0,len(theStartPosition)-1):
	myList=[] #暂时承装某一条染色体的所有碱基序列
	for j in range(theStartPosition[i]+1,theStartPosition[i+1]):
		a[j] = a[j].upper(); #把每一行的字符串都转换成大写
		myList += list(a[j])[0:-1]
	dict[indexToName[i]]=myList

#下面是补上最后一个染色体
myList = []
for j in range(theStartPosition[-1]+1, len(a)):
	a[j] = a[j].upper(); #把每一行的字符串都转换成大写
	myList += list(a[j])[0:-1]

dict[indexToName[-1]]=myList


___________________________________
# -*- coding: utf-8 -*-
#处理gff文件并根据上一个脚本生成的字典，进行查找。
#在第一个脚本生成dict（染色体的序列字典）后，将gff文件进行处理并根据dict截取出gff，保存在一个字典cdsDict中
# import sys
import re

gff = open('GCF_000001405.39_GRCh38.p13_genomic.gff',mode='r') 

cdsDict = {} #用于装每条CDS对应的碱基序列。



def insertToDict(dict, CDS_Name, CDS_Sequence):
	if CDS_Name in dict:
		dict[CDS_Name]+= CDS_Sequence
	else:
		dict[CDS_Name] = CDS_Sequence


for line in gff:
	if re.match('#',line): #跳过#开头的行。
		continue
	myList = re.split('\t|;',line) #注意这里是用制表符匹配，如果用空格匹配的话则出错。因为文件本身就是制表符
	if myList[2]=='CDS':
		#start = int(myList[3]) + int(myList[7]) #注意文本中的数字也是字符串的，所以要先转换成int
		start = int(myList[3])
		end = int(myList[4])
		
		#if myList[8] in cdsDict:
		#	start = int(myList[3])
		#else:
		#	start = int(myList[3]) + int(myList[7]) #注意文本中的数字也是字符串的，所以要先转换成int
		#end = int(myList[4])
		
		if myList[6]=='+':
			insertToDict(cdsDict, myList[8], dict[myList[0]][start-1:end])
		else:
			if(start == 1):
				cds = dict[myList[0]][-end:]; #倒着切片然后reverse翻转(当包含最后一个元素时，就不能[-4:0]这么写了，要省略这个0，所以这里特殊处理)。
				cds.reverse()
			else:
				cds = dict[myList[0]][-end:-start+1]; #倒着切片然后reverse翻转。		
				cds.reverse()
			match = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
			for i in range(0,len(cds)): #按配对方式替换
				cds[i] = match[cds[i]]
			insertToDict(cdsDict, myList[8], cds)
		

gff.close()



____________________________________________________

#直接用上一问的列表cdsList，对应成蛋白质。把蛋白质id和序列保存到字典中。

codonDict = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

proteinSequenceDict = {} #用于装各个CDS的蛋白质序列（每个蛋白质序列是一个字符串）
#对于每一个CDS序列


for key in cdsDict:
	protein = "" #用于承装该序列的所有氨基酸
	#对于该序列中的每一个三连密码子
	for j in range(0,len(cdsDict[key]),3):
		try:
			codon = cdsDict[key][j]+cdsDict[key][j+1]+cdsDict[key][j+2] #将三联密码子相加变成字符串			
		except: #如果密码子不是3的倍数，最后数量差了的话，就break，后面这个就不算了。
			break
		try:
			protein += codonDict[codon]
		except:
			protein += '*' #  *表示未知蛋白
	proteinSequenceDict[key] = protein
