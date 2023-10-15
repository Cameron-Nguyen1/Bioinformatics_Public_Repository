# Regressional analysis that aims to compare slope

model = lm(Run_1~Week,df)  ##Fit data to model
model_prediction = as.data.frame(predict(model,interval="predict"))
sum1 = summary(model)

plot1 = ggplot(df, aes(x=Week,y=Run_1))+ ##Plot
ggtitle("LinReg Run_1 2012c Antibody Decay")+
geom_text(x=20,y=15,label=paste0("y = ",signif(model$coefficients[2],4),"x + ",signif(model$coefficients[1],6)))+
geom_point()+
stat_smooth(method = lm) +
geom_line(aes(y=model_prediction$lwr),col="coral2",linetype="dashed")+
geom_line(aes(y=model_prediction$upr),col="coral2",linetype="dashed")+
ylab("Log_2 Value")

model2 = lm(Run_2~Week,df)  ##Fit data to model
model2_prediction = as.data.frame(predict(model2,interval="predict"))
sum2 = summary(model2)

plot2 = ggplot(df, aes(x=Week,y=Run_2))+ ##Plot
ggtitle("LinReg Run_2 2012c Antibody Decay")+
geom_text(x=20,y=12,label=paste0("y = ",signif(model2$coefficients[2],4),"x + ",signif(model2$coefficients[1],6)))+
geom_point()+
stat_smooth(method = lm) +
#geom_line(aes(y=model2_prediction$lwr),col="coral2",linetype="dashed")+
#geom_line(aes(y=model2_prediction$upr),col="coral2",linetype="dashed")+
ylab("Log_2 Value")

coefficient1 = sum1$coefficients[2]
coefficient2 = sum2$coefficients[2]
se1 = sum1$coefficients[4]
se2 = sum2$coefficients[4]

#Modified Z-test tests equality of regression coefficients
sig = (coefficient1 - coefficient2)/(sqrt((se1^2)+(se2^2)))