library(ggplot2)


#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
  # for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
# to be summariezed
# groupnames : vector of column names to be used as
# grouping variables

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

##############################

## Loading Data
df <- read.csv("myfile.csv")

## Calculate Statistics
df2 <- data_summary(df, varname="len", groupnames=c("supp", "dose"))

# Convert dose to a factor variable
df2$dose=as.factor(df2$dose)


## Plot bar chart
p<- ggplot(df2, aes(x=dose, y=len, fill=supp)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2, position=position_dodge(.9)) 
