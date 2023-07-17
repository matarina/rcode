#数据下载
rm(list = ls())
library(GEOquery)
setwd('H:/liver/')
gse_number = "GSE76427"  #这里你只需要改成你要下载文献的GEO号
eSet <- getGEO(gse_number, 
               destdir = './', 
               getGPL = F)
length(eSet)
eSet = eSet[[1]]
#(1)提取表达矩阵exp
exp <- exprs(eSet)    #表达矩阵在assay data中exprs
exp[1:4,1:4]
if(T){exp = log2(exp+1)}#避免有0值的存在 加1不影响大小（观察表达矩阵中有无取过log）
#boxplot(exp)   #观察是否差不多的箱长，如果不一样
if(T){exp = limma::normalizeBetweenArrays(exp)}  
#boxplot(exp)   #这时候箱长就会一样的
#(2)提取临床信息
pd <- pData(eSet)
pd
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号  
gpl_number <- eSet@annotation
Group = rep(c("treat","control"),times = c(115,52)) 
#为了大家分析不出错，建议大家把实验组的sample调到前面，control组的sample调到后面，
#这样就可以无脑运行后面的代码了，并且匹配设置group
#设置参考水平，指定levels，对照组在前，处理组在后
Group = factor(Group,
               levels = c("control","treat"))  #这里一定要这样设置因子的位置
Group                                       #因为后面的代码是按照这样匹配的                                           


library(tinyarray)
  #假如你自己的数据，请自行给gpl_number赋予平台number
ids <- AnnoProbe::idmap('GPL10558')  #在这里就可以看到ids有两列 一个是探针 一个是symbol





#差异分析，用limma包来做
#需要表达矩阵和Group，不需要改
library(limma)
design=model.matrix(~Group) #根据group分组信息进行贝叶斯检验得出差异分析logFC以及p value
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)


#第一种
#利用tinyarray的包


#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))  #mutate函数给数据框添加一列probe_id 内容是rownames(deg)
#2.加上探针注释      
ids = ids[!duplicated(ids$symbol),]  #探针去重 自行了解       
deg <- inner_join(deg,ids,by="probe_id")    #将deg与ids合并  自动把deg中重复的probe_id给删除
nrow(deg)

res1 = deg
library(EnhancedVolcano)
keyvals <- ifelse(
  res1$logFC < -1, '#0081B4',
  ifelse(res1$logFC > 1, '#ff6361',
         'black'))
keyvals[is.na(keyvals)] <- '#B6EADA'
names(keyvals)[keyvals == '#ff6361'] <- 'high'
names(keyvals)[keyvals == '#B6EADA'] <- 'mid'
names(keyvals)[keyvals == '#0081B4'] <- 'low'
volcano = EnhancedVolcano(res1,legendPosition = 'none',
                          FCcutoff = 1,       
                          gridlines.minor=FALSE, gridlines.major=FALSE,
                          colCustom = keyvals,#set color to custom
                          selectLab = rownames(res1)[which(names(keyvals) %in% c('high', 'low'))],
                          lab =NA,title = NULL,subtitle = NULL,#legendLabels = NULL,
                          caption = '',hline = NULL,vline = NULL,
                          x = 'logFC',
                          y = 'P.Value',parseLabels = FALSE,
                          border = 'full',
                          borderWidth = 1,pointSize = 2.0,
                          axisLabSize = 6,
                          pCutoff = 0.05)

volcano
ggsave('H:/liver/volcano.pdf',plot = volcano,width = 5.8,height =5.7,units = 'in')

#get diff_genes in immue pathway 
deg = na.omit(deg)
diff = deg %>% filter(P.Value < 0.05 & abs(logFC) >0.7)
imm = read.csv('H:/liver/immue_genes.csv',header = FALSE)
imm_gene = c(imm$V1,imm$V2,imm$V3,imm$V4,imm$V5)
imm_diff = intersect(diff$symbol,imm_gene)
write.csv(imm_diff,file = 'immue_diffgene.csv',quote = FALSE)

#get high ppi socre gene 
string_interaction = read.table(file = 'string_interactions.tsv', sep = '\t', header = TRUE)
intergene  = string_interaction %>% filter(combined_score > 0.5) %>% distinct(node1, .keep_all= TRUE)
intergene$node1

#lassox_cox 
library(tidyverse)
expr2 <- mutate(as.data.frame(exp),probe_id=rownames(exp))  
idexp <- inner_join(expr2,ids,by="probe_id")   
idexp = idexp %>% column_to_rownames(var="symbol") %>% select(-(probe_id)) 
idexp = idexp %>% filter(row.names(idexp) %in% imm_diff)


library(ggcorrplot)


cor = round(cor(t(idexp),use = "complete.obs"),3)
head(cor)

pmat = round(cor_pmat(t(idexp)),)

ggcorrplot(cor,hc.order = TRUE, colors = c("#153462", "white", "#ff6361"), tl.cex = 3,
          
           outline.color = "white",ggtheme = theme_bw())

















pd = pd %>% select(`duryears_os:ch1`,`event_os:ch1`) %>% rename(osyear = 1,os = 2)
input = na.omit(merge(t(idexp),pd, by = 'row.names' ,all = FALSE))

library("glmnet") 
library("survival")
x=input[,c(2:(ncol(input)-2))]
x2 = x %>% mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))  
y=data.matrix(Surv(as.numeric(input$osyear),as.numeric(input$os)))

cvfit <- cv.glmnet(data.matrix(x2), y, family = 'cox', type.measure = 'deviance', nfolds = 10)
plot(cvfit)

lambda.1se <- cvfit$lambda.1se
lambda.1se

cvfit$lambda.min
coefficient <- coef(cvfit, s = "lambda.min")

#系数不等于0的为纳入的变量（基因）
Active.index <- which(as.numeric(coefficient) != 0)
Active.coefficient <- as.numeric(coefficient)[Active.index]
sig_gene_mult_cox <- rownames(coefficient)[Active.index]
#查看具体哪些基因
sig_gene_mult_cox
rownames(input) = input$Row.names
input = input[,-1]
data2_cox <- input %>% dplyr::select(os,osyear,all_of(sig_gene_mult_cox))

multiCox <- coxph(Surv(as.numeric(osyear),as.numeric(os)) ~ ., data =  data2_cox)

summary(multiCox)

#predict函数计算风险评分
riskScore=predict(multiCox,type="risk",newdata=data2_cox) 
riskScore<-as.data.frame(riskScore)

riskScore$sample <- rownames(riskScore)
head(riskScore,2) 



######riskScore 二分绘制KM##########
riskScore_cli <-merge(input,riskScore,by = 'row.names')
#按照中位数分为高低风险两组
riskScore_cli$riskScore2 <- ifelse(riskScore_cli$riskScore > median(riskScore_cli$riskScore),
                                   "High","Low")
#KM分析
fit <- survfit(Surv(as.numeric(osyear), as.numeric(os)) ~ riskScore2, data=riskScore_cli)
library('survminer')
plot(fit)

ggsurvplot(fit, data = riskScore_cli,
           pval = T,size = 0.7,
           risk.table = T,
           surv.median.line = "hv", #添加中位生存曲线
           palette=c("#FEC260","#153462"),  #更改线的颜色
           #legend.labs=c("High risk","Low risk"), #标签
           legend.title="RiskScore",
           title="Overall survival", #标题
           ylab="Cumulative survival (percentage)",xlab = " Time (Days)", #更改横纵坐标
           censor.shape = 124,censor.size = 2,conf.int = FALSE, #删失点的形状和大小
           break.x.by = 720#横坐标间隔
)
ggsurvplot(fit,
           conf.int = TRUE,# 显示置信区间
           linetype = "strata", # 根据性别分组自动设置曲线类型
           surv.median.line = "hv", # 设置中位生存期显示
           ggtheme = theme_bw(), # 设置ggplot2主题
           palette = c("#E7B800", "#2E9FDF"))
km_fit <- survfit(Surv(time, status) ~ sex, data = lung)
plot(km_fit, main = "KM Survival Curve for Lung Cancer Patients")

riskScore_cli$os = as.numeric(riskScore_cli$os)

riskScore_cli$osyear = as.numeric(riskScore_cli$osyear)

#ROC

library(timeROC)
result = with(riskScore_cli,
              ROC_riskscore <<- timeROC(T = osyear,
                                        delta = os,
                                        marker = riskScore,
                                        cause = 1,
                                        weighting = "marginal",
                                        times = c(365,1095,1825),
                                        ROC = TRUE,
                                        iid = TRUE)
)
plot(ROC_riskscore, time = 365, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 1080, col = "blue", add = T)
result <-with(new_dat, timeROC(T=time,
                               delta=event,
                               marker=fp,
                               cause=1,
                               times=c(365,1095,1825),
                               iid = TRUE))
#identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(365,1095,1825)),each = nrow(result$TP)))

library(ggplot2)
lasso_roc = ggplot() + 
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 0.7) + 
  scale_color_manual(name = NULL,values = c("#FEC260","#153462", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
lasso_roc
ggsave('H:/prna_transcript_analysis/annotationnewplot/fig_backup/lasso_roc.pdf',plot = lasso_roc,width = 2.8,height =2.7,units = 'in')

