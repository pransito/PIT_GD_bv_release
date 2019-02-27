## get the distribution of KFG and SOGS
KFG_SOGS = dat_match[c('HCPG','KFG','SOGS')]
KFG_SOGS[c('KFG','SOGS')] = KFG_SOGS[c('KFG','SOGS')] + randn(dim(KFG_SOGS)[1],dim(KFG_SOGS)[2]-1)*0.1
KFG_SOGS$Group= as.factor(agk.recode(as.character(KFG_SOGS$HCPG),c('HC','PG'),c('HC','GD')))

p = ggplot(KFG_SOGS, aes(KFG, SOGS))
p = p + geom_point(size = 4, shape = 17)
p = p + geom_point(aes(colour = Group), size = 4, shape = 17)
p = p + xlab('Score Short Questionnaire Pathological Gambling (KFG)')
p = p + ylab('Score South Oaks Gambling Screen (SOGS)')
p = p + geom_hline(yintercept = 4) # linetype, color, size)
p = p + geom_vline(xintercept = 15) # linetype, color, size)
p = p + ggtitle('Sorting subjects into groups according to score on the KFG')
p = p + theme(axis.text=element_text(size=14, face = "bold"),
              axis.title=element_text(size=20,face="bold"))
p = p + theme(plot.title = element_text(size=22))
p = p + theme(legend.text = element_text(size=18))
p = p + theme(legend.title= element_text(size=18))
print(p)