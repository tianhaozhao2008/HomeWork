IFS_bak=$IFS   #备份一下IFS的值。IFS是系统默认变量，储存默认的分隔符
IFS=$'\t'      #把默认分隔符改成制表符
declare -A result     #声明关联数组用于储存最后结果，即数组下标可以是字符串等。


while read line  #对于读入的每一行
do
  arr=($line)   #把这行字符串拆分成数组
  for ((i=2;i<${#arr[*]};i++)) #对于每个通路，都把arr[0]的蛋白质加进去
  do 
    #下面这句，无论result[${arr[i]}]是否有初始化，都可以直接拿来用。如果没初始化，那么字符串拼接时就默认自己是空字符串。
    result[${arr[i]}]="${result[${arr[i]}]}  ${arr[0]}"
  done
done < pathway.txt

#IFS=$IFS_bak    #用完IFS后，恢复它的值[]。这里这样写会导致bug，具体解释参考readme.txt


#最后输出我们的结果：各个通路所包含的蛋白质
for index in ${!result[*]} #遍历所有下标
do #下面的-e是可以使转义字符换行\n生效
  echo -e "$index 通路中的蛋白质包括：${result[$index]} \n"
done
