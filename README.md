# 需要用到的自定义函数
## Classification_picture  
此函数为分类图，返回样本分类图。根据偏最小二乘分析的结果对象plsda对样本进行分类。  
  plsda：偏最小二乘分析的结果对象  
  genderFc：样本分类对象  
  sample.score：plsda中得分矩阵  
  注：为得到合适的图形，不同的数据集需要调节图形范围、置信系数、甚至成像函数。通常情况下成像选用函数stat_ellipse（），但是当样本数据非常少时，不能画出置信椭圆，需选用函数geom_encircle（）
```r
Classification_picture <-function(plsda, genderFc)
{
  sample.score = plsda@scoreMN#得分矩阵
  sample.score=as.data.frame( sample.score)##转换成数据框
  #genderFc<-as.integer(genderFc)
  sample.score$gender<-genderFc##增加gender列
  # 解释率
  x_lab <- plsda @modelDF[1, "R2X"] * 100
  y_lab <- plsda @modelDF[2, "R2X"] * 100
  # 画图
  p = ggplot(sample.score, aes(p1, p2, color =factor( gender))) +
    ###在排序图中根据个体属性（gender）给样本上色
    geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    geom_point(size=2) +
    #geom_point(aes(-10, 10), color = 'white') + xlim(-10,10)+ylim(-10,10)+
    #geom_point(aes(-20, 20), color = 'white') + xlim(-20,20)+ylim(-20,20)+
    geom_point(aes(-40, 40), color = 'white') +#xlim(-40,40)+ylim(-50,50)+
    labs(x = paste0("p1(", x_lab, ")"), y = paste0("p2(", y_lab, ")"))+
    stat_ellipse(level = 0.98, linetype = 'solid', linewidth = 1, show.legend = FALSE) +
    ##数据太少时，可能画不出置信（0.95）椭圆，改用下面的函数替代stat_ellipse
    #geom_encircle(aes(group = gender), alpha = 0.9,color="black", size=2,expand=0.05, show.legend = F) +
    ###代谢物数据专用 第一次expand=0.2，aes(-40, 40)第二次expand=0.15，aes(-20, 20)第三次expand=0.05，aes(-10, 10)
    scale_color_manual(values = c('#008000','#FFA74F', '#001000', '#004000', '#002200')) +
    theme_bw() +
    theme(legend.position = c(0.8,0.9),
          legend.text = element_text(color = 'black', size = 10, face = 'plain'),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black', size = 8, face = 'plain'),#坐标
          axis.title = element_text(color = 'black', size = 10, face = 'plain'),
          axis.ticks = element_line(color = 'black'))
  p
  ggsave(p, filename = 'pls.pdf', width = 5, height = 5, device = cairo_pdf)
  return  (p)
}
 ```
## Arg_VIP_Order
该函数根据自变量的重要程度得分进行排序，达到逆序排列的自变量。  
  plsda：偏最小二乘分析的结果对象  
  dd：阈值，舍弃VIP<dd的自变量  
  vip.score：自变量的重要程度得分  
  P：自变量VIP值图形，保存在文件pls_VIP_0.pdf中  
```r
Arg_VIP_Order<- function(plsda,dd)
{
  vip.score = as.data.frame(plsda@vipVn)  ##每个指标的重要程度得分
  colnames(vip.score) = 'vip'
  ### colnames函数来指定矩阵的列名称，vip.score只有一列
  vip.score$metabolites = rownames(vip.score)
  ## vip.score增加一列，列名metabolites，值为vip.score行名
  vip.score = vip.score[order(-vip.score$vip),]
  ##按照vip.score$vip，逆序排列
  #write.csv(as.data.frame(vip.score)," PLSDA_VIP.csv",row.names =FALSE)
  #VIP值保存到文件里面
  vip.score$metabolites=factor(vip.score$metabolites,levels= vip.score$metabolites) ## 转换成因子
  p = ggplot(vip.score[vip.score$vip >=dd,], aes(metabolites, vip)) +
    geom_segment(aes(x = metabolites, xend = metabolites,
                     y = 0, yend = vip)) +
    geom_point(shape = 21, size= 5, color = '#008000' ,fill = '#008000') +
    geom_point(aes(1,2.5), color = 'white') +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = '', y = 'VIP value') +
    theme_bw() +
    theme(legend.position = 'none',
          legend.text = element_text(color = 'black',size = 12, face = 'plain'),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black',size = 15, face = 'plain'),
          axis.text.x = element_text(angle = 90),
          axis.title = element_text(color = 'black',size = 15, face = 'plain'),
          axis.ticks = element_line(color = 'black'),
          axis.ticks.x = element_blank())
  ggsave(p, filename = 'pls_VIP_0.pdf',width = 20, height = 5, device = cairo_pdf)
  return (vip.score)
}
```
## Arg_screening:
此函数用于筛选R2Y和Q2最大的自变量组合。将每个plsda中R2Y和Q2的变化曲线保存到指定文件中，并返回筛选的自变量个数（自变量是从大到小排序的）。  
  plsda：偏最小二乘分析的结果对象  
  dataMatrix：做plsda分析的数据框，用于提取自变量对应的样本数据  
  genderFc：样本的分类变量  
  dd：阈值，仅筛选VIP>1的自变量  
  zhimu：标识符，在R2Y和Q2的变化曲线上角做标记。  
  vip.score：自变量的重要程度得分  
  333_300.tiff:保存R2Y和Q2的变化曲线，分辨率为300
```r
Arg_screening<- function(plsda, dataMatrix, genderFc,dd,zhimu) ###循环验证
{
  vip.score = as.data.frame(plsda@vipVn)  #用PLSDA的VIP值
  colnames(vip.score) = 'vip'
  vip.score$metabolites = rownames(vip.score)
  vip.score<-vip.score[which(vip.score$vip>dd),]
  vip.score = vip.score[order(-vip.score$vip),]
  N= nrow(vip.score)
  R2X <- array(1: N)
  R2X<-rep(0,N)
  R2Y <- array(1: N)
  R2Y <-rep(0, N)
  Q2 <- array(1: N)
  Q2 <-rep(0, N)
  MA<-0
  MI<-0
  for(i in seq(2,N,by=1))
  {         
    otu_select <- rownames(vip.score)[1:i]
    Top25 <- dataMatrix [ ,c(otu_select)]
    plsda <-opls(Top25,genderFc) #plsda <-opls(Top25,genderFc, predI = 1, orthoI = NA)
    aa<-plsda@summaryDF
    R2X [i]<-aa[,1]
    R2Y [i]<-aa[,2]
    Q2 [i]<-aa[,3]  
    if(MA<aa[,3])
    {
      MA<-aa[,3] 
      MI<-i
    } 
  }
  colnames(vip.score) = 'vip' ### colnames函数来指定矩阵的列名称，vip.score只有一列
  vip.score$metabolites = rownames(vip.score) ## vip.score增加一列，列名metabolites，值为vip.score行名
  vip.score = vip.score[order(-vip.score$vip),]
  ##按照vip.score$vip，逆序排列
  #write.csv(as.data.frame(vip.score)," PLSDA_VIP_new.csv",row.names =FALSE)
  Q2<-as.data.frame(Q2)
  Q2$num<-rownames(Q2)
  Q2<-Q2[-1,]  #一个指标不能做模型，从前两个开始，所以删除第一个
  R2X <-as.data.frame(R2X)
  R2X $num<-rownames(R2X)
  R2X <- R2X [-1,]
  R2Y <-as.data.frame(R2Y)
  R2Y $num<-rownames(R2Y)
  R2Y <- R2Y [-1,]
  #Q2<-Q2[Q2!=0,] #删除0值
  Q2 <- Q2[!is.na(Q2$Q2), ]#删除缺失值
  R2Y <- R2Y [!is.na(R2Y $ R2Y), ] 
  R2X <- R2X [!is.na(R2X $ R2X), ] 
  tiff(file="333_300.tiff",compression="lzw",units="in",res=300,pointsize=8,height=4,width=4)##  1in=2.54cm
  par(oma=c(1,1,1,1))
  plot(Q2$num,Q2$Q2,xlab="num",ylab="value", ylim =c(0.5,1.0), col="red",pch=2,cex=1.2,lwd=2,type='l')
  axis(side=1,at=c(MI),label=c(MI),col.axis="red",line=-1)
  lines(R2Y$num,R2Y$R2Y,pch=2,cex=1.2,lwd=2,col="blue",type='l')
  legend("right", c("Q2", "R2Y"), lwd=2,col=c("red","blue"),bg ="white")
  #abline(h = MA)
  abline(v = MI) 
  pushViewport(viewport(x=0.02, y=0.97, width=1, height=1, angle=0))
  grid.text(zhimu, gp=gpar(col="black", cex=2))
  dev.off()
  return (MI)
}
```
## BoxPlot 
画指定变量Top的箱线图，返回箱线图对象  
  Top:变量数据框  
  dd：箱线图的行数  
```r
BoxPlot<- function(Top,dd) 
{
  Top<- as.data.frame(Top)
  Top<-scale(Top)
  Top<- as.data.frame(Top)
  Top$group <- genderFc
  tail(Top)
  dat2 = gather(Top,key = "gene",value = "expression",- group)
  #dat2$gene=factor(dat2$gene,ordered = TRUE,levels = paste0("gene",1:10))
  df_p_val1 <- dat2 %>% 
    group_by(gene) %>% 
    wilcox_test(formula = expression ~ group) %>% 
    add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
    add_xy_position(x='group')
  
  PP<-ggplot(data = dat2)+
    geom_boxplot(aes(x = group,y = expression,fill = group))+# color = group
    scale_fill_manual(values = c('#E21C21','#3A7CB5'))+
    stat_pvalue_manual(df_p_val1,label = '{p}',tip.length = 0)+# p.signif
    theme_bw()+
    #theme(axis.text = element_text(color = 'black'),
    #  legend.position = c(0.7,0.1),
    #   legend.direction = 'horizontal')+
    facet_wrap(~gene,nrow =dd)
  return  (PP)
}
```
## Save_picture
保存PLSDA与分类图。保存为分辨率为300的图像。  
  plsda：PLSDA对象  
  p1:分类图对象  
  zhi1: plsda图像标识  
  zhi2: p1图像标识  
  111_300.tiff： plsda图像名称  
  222_300.tiff： p1图像名称  
```r
Save_picture<-function(plsda, p1,zhi1,zhi2)
{
  tiff(file="111_300.tiff",compression="lzw",units="in",res=300,pointsize=8,height=4,width=4)##  1in=2.54cm
  par(oma=c(1,1,1,1))
  plot(plsda)
  pushViewport(viewport(x=0.02, y=0.97, width=1, height=1, angle=0))
  grid.text(zhi1, gp=gpar(col="black", cex=2))
  dev.off()
  tiff(file="222_300.tiff",compression="lzw",units="in",res=300,pointsize=8,height=4,width=4)##  1in=2.54cm
  par(oma=c(1,1,1,1))
  plot(p1)
  pushViewport(viewport(x=0.02, y=0.97, width=1, height=1, angle=0))
  grid.text(zhi2, gp=gpar(col="black", cex=2))
  dev.off()
}
```
## Merge_picture
将plsda图、分类图和筛选图合并为一幅图像，三个图像横向排列。返回合并图像对象。合并的图像存储为分辨率为300的图像。
  name1：分辨率为300的图像名称
  dd:合并图像的数量，有时候只合并plsda图、分类图，则dd=2
```r
Merge_picture <-function(name1,dd)
{
  if(dd==3)
  {
    a1<-image_read("111_300.tiff")
    a2<-image_read("222_300.tiff")
    a3<-image_read("333_300.tiff")
    p1<-image_append(c(a1,a2,a3))
    image_write(p1,name1)
  }
  if(dd==2)
  {
    a1<-image_read("111_300.tiff")
    a2<-image_read("222_300.tiff")
    p1<-image_append(c(a1,a2))
    image_write(p1,name1)
  }
}
```
# Example
使用data文件夹中的GSE90028.xls作为例子，对RPLS相关算法进行复现。  
为了便于代码的理解，现对后续代码中出现的部分变量进行说明：  
  dataMatrix1：训练集样本数据框  
  genderFc1：  训练集样本分类变量  
  dataMatrix2：验证集样本数据框  
  genderFc2：  验证集样本分类变量  
  dataMatrix： 训练集转化为秩的数据集  
  genderFc：   PLSDA分析的分类变量  
  Top42：      迭代过程中用于衔接的中间变量

## 1.Adjust the format of raw data
### raw data of training set
```r
table_test1<- read_excel ("GSE90028.xls",2)
dataMatrix1<-table_test1 
dataMatrix1 <-dataMatrix1%>%mutate(p=NULL, zhi1=NULL, zhi2=NULL, cpd_ID =NULL, HMDB =NULL)
dataMatrix1 <-t(dataMatrix1) 
dataMatrix1 <-as.data.frame(dataMatrix1 ) 
colnames(dataMatrix1 )= dataMatrix1 [1,]
dataMatrix1 <-dataMatrix1 [-1,] 
genderFc1<-as.factor(dataMatrix1 $ MYKIND)# 
n<-ncol(dataMatrix1 ) 
dataMatrix1 <-dataMatrix1 [,-n]
rowname1 <- rownames(dataMatrix1 )
colname1 <-colnames(dataMatrix1 )
dataMatrix1 <- as.data.frame(lapply(dataMatrix1 ,as.numeric)) 
rownames(dataMatrix1 )= rowname1 
colnames(dataMatrix1 )= colname1
plsda1 = opls(dataMatrix1,genderFc1)##PLSDA，, predI = 2
if(plsda1@summaryDF$pre==1)
{
  plsda1 = opls(dataMatrix1,genderFc1 , predI = 2)
}
```
### raw data of validation set
```r
table_test2<- read_excel ("GSE90028.xls",3)
dataMatrix2<-table_test2  
dataMatrix2 <-dataMatrix2%>%mutate(p=NULL, zhi1=NULL, zhi2=NULL, cpd_ID =NULL, HMDB =NULL)
dataMatrix2 <-t(dataMatrix2)  
dataMatrix2 <-as.data.frame(dataMatrix2) 
colnames(dataMatrix2 )= dataMatrix2 [1,]
dataMatrix2 <-dataMatrix2 [-1,] 
genderFc2<-as.factor(dataMatrix2 $ MYKIND)
n<-ncol(dataMatrix2 )  
dataMatrix2 <-dataMatrix2 [,-n]
rowname2 <- rownames(dataMatrix2 )
colname2 <-colnames(dataMatrix2 )
dataMatrix2 <- as.data.frame(lapply(dataMatrix2 ,as.numeric)) 
rownames(dataMatrix2 )= rowname2
colnames(dataMatrix2 )= colname2
plsda2 = opls(dataMatrix2,genderFc2)
```
## 2.Data analysis model based on RPLS_DA
### Data transformation
calculate the rank corresponding to the raw data, that is, sort the metabolite data in ascending order, and the subscript of the sorting is the rank of the raw data.
```c
aaaaaaaaaaaaaaaaaaaa
```
### Adjust the format of the rank data of training set and establish PLS model based on the rank data
```r
table_test<- read_excel ("GSE90028.xls",4)
dataMatrix<-table_test 
dataMatrix <-dataMatrix%>%mutate(p=NULL, zhi1=NULL, zhi2=NULL, cpd_ID =NULL, HMDB =NULL)
dataMatrix <-t(dataMatrix)  
dataMatrix <-as.data.frame(dataMatrix ) 
colnames(dataMatrix )= dataMatrix [1,]
dataMatrix <-dataMatrix [-1,] 
genderFc<-as.factor(dataMatrix $ MYKIND)
n<-ncol(dataMatrix )  
dataMatrix <-dataMatrix [,-n]
rowname <- rownames(dataMatrix )
colname <-colnames(dataMatrix )
dataMatrix <- as.data.frame(lapply(dataMatrix ,as.numeric)) 
rownames(dataMatrix )= rowname 
colnames(dataMatrix )= colname
#head(dataMatrix)
## Establish PLS model
plsda <-opls(dataMatrix,genderFc)
if(plsda@summaryDF$pre==1)
{
  plsda = opls(dataMatrix,genderFc , predI = 2)
}
Top42 <- dataMatrix
```
## 3.Variable screening model based on RPLS_DA
```r
dataMatrix0<- Top42
p1<- Classification_picture (plsda, genderFc) ###分类图
vip.score <- Arg_VIP_Order (plsda,dd=1) ####获取VIP值
Save_picture (plsda,p1, "a","b") ##保存PLSDA图形和分类图形
MI<- Arg_screening (plsda, dataMatrix0, genderFc,dd=1.,"c") ###根据Q2Y和R2Y曲线图，筛选自变量。
#由于Q2Y和R2Y可能不会同时达到极值，因此需要根据两条曲线的峰值人工设置筛选值，但主要参考的是Q2Y值。通常num就是Arg_screening返回的MI，也可在下面根据图形人工设置num值
#MI<-num
#将每次迭代的plsda图、分类图、筛选过程中Q2Y和R2Y曲线图合并为一个图
Merge_picture ("0000_300.tiff ",3) 
#Merge _picture ("0000_300.tiff ",2)
##获取筛选出的自变量
otu_select <- rownames(vip.score)[1: MI]
```
## 4.预测并统计预测准确率
```r
##使用raw data建立模型
Top25 <- dataMatrix1 [ ,c(otu_select)] ## Top25:筛选自变量otu_select对应的原训练集样本数据
sacurine.oplsda<-opls(Top25, genderFc1)
##使用validation set进行预测
Top250 <- dataMatrix2 [ ,c(otu_select)] ## Top250：otu_select对应的验证集样本数据
tab=table(genderFc2,predict(sacurine.oplsda, Top250))
P <- round(sum(diag(tab))/sum(tab)*100,2)
```
## 5.Iterative 
Assume the number of candidate biomarkers `otu_select` is N, the determine whether to continue the iteration based on N and the prediction rate P: if N＜Ne (the expected final number of  remaining biomarkers) or P＜Pe (the expected prediction rate), terminate the iteration (such as Ne = 10, Pe = 60%), and the markers in variable `otu_select` are the final screened biomarers. Otherwise, let `Top42 <- dataMatrix0 [ ,c(otu_select)]` and start the iteration from Step “3.Variable screening model based on RPLS_DA”.  
It should be highlighted that when the prediction rate P is lower than the expected value, if the number of remaining metabolites N is still large, the iteration can continue until (N＜Ne). This is because the prediction rate may increase again after removing redundant metabolites. Therefore, if P shows a fluctuating trend, the size of N should be comprehensively considered to select the stopping point;  if P shows a decreasing trend, the point where P is at its maximum should be chosen as the stopping point.  



