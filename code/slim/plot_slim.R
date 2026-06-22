# Analyze results of SLiM simulations and plot

library(ggplot2)

setwd("your_dir")

# Function to calculate fitness
fitness_calc<-function(z,opt,omega) {
  Lambda<-solve(omega)
  D2<-t(z-opt)%*%Lambda%*%(z-opt)
  w<-exp(-D2/2)
  return(as.numeric(w))
}

n_share_all=c(0,10,20,30,40,50) # Number of pleiotropic cis-elements
cor_selection_all=c(0,0.9,-0.9) # Regimes of correlational selection (selection gradient aligned with SLLR)

# Data matrix of parameter combinations to examine
par_comb_all=matrix(0,nrow=1,ncol=2)
for(i in 1:length(n_share_all)){
    for(j in 1:length(cor_selection_all)){
      row=c(n_share_all[i],cor_selection_all[j])
      par_comb_all=rbind(par_comb_all,row)
    }
}
par_comb_all=par_comb_all[2:nrow(par_comb_all),]
colnames(par_comb_all)=c("n_share","cor_selection") # Add column names for plotting
rownames(par_comb_all)=NULL # Remove row names

prefix="sim_out_all_"
for(c in 1:nrow(par_comb_all)){
  n=par_comb_all[c,1]
  cor_selection=par_comb_all[c,2]
  fn=paste(prefix,n,"_",cor_selection,".txt",sep="")
  # Decide how many lines to remove (dependent on types of genetic elements initialized)
  if(n!=0&n!=50){
    start=23 # 3 genetic element types 
  }else{
    if(n==0){
      start=22 # 2 genetic element types (unique targets of 2 regulators only)
    }else{
      start=21 # 1 genetic element type (shared targets of 2 regulators only)
    }
  }
  d<-readLines(fn) # Read file
  d=d[start:length(d)] # Remove starting lines
  fn_new=paste(prefix,n,"_",cor_selection,"_processed.txt",sep="") # New file name
  writeLines(d,fn_new) # Write processed file
}

# Re-read processed data files and analyze
out_list=list() # A list of data matrices (each matrix corresponding to a parameter combination)
for(c in 1:nrow(par_comb_all)){
  out=matrix(0,nrow=100,ncol=6) # Data matrix corresponding to this parameter combination
  
  n=par_comb_all[c,1]
  cor_selection=par_comb_all[c,2]
  if(cor_selection==0){
    opt=c(log(100),0)
  }else{
    if(cor_selection==0.9){
      opt=c(log(100),log(100))
    }else{
      opt=c(log(100),-log(100))
    }
  }
  mat_selection=rbind(c(1,cor_selection),c(cor_selection,1))
  
  out[,1]=n;out[,2]=cor_selection
  fn_new=paste(prefix,n,"_",cor_selection,"_processed.txt",sep="") # New file name
  d<-read.table(fn_new,sep="\t") # Re-read rearranged file
  
  for(i in 1:100){
    dsub=d[which(d[,1]==100*i),] # Extract rows corresponding to the focal time point
    out[i,3]=100*i # Time
    out[i,4]=mean(dsub[,3]) # Mean of log(z1) (converted to difference relative to ancestral state)
    out[i,5]=mean(dsub[,4]) # Mean of log(z2) (converted to difference relative to ancestral state)
    out[i,6]=fitness_calc(z=out[i,4:5],opt=opt,omega=mat_selection)
  }
  out_list[[c]]=out
}

out_mat=out_list[[1]]
for(c in 1:nrow(par_comb_all)){
  out_mat=rbind(out_mat,out_list[[c]])
}
colnames(out_mat)=c("n_share","cor_selection","t","z1","z2","w")
write.table(out_mat,file="sum_slim_temp.txt",sep="\t")

d<-read.table("sum_slim_temp.txt",sep="\t")
d$n_share=factor(d$n_share,levels=sort(unique(d$n_share)),ordered=TRUE)

# Subsets based on correlational selection
d1=d[which(d$cor_selection==0),]
d2=d[which(d$cor_selection==0.9),]
d3=d[which(d$cor_selection==-0.9),]

# Plot z1 (log-scale) against time
g1<-ggplot(d1,aes(x=t,y=z1,col=n_share))
g1=g1+geom_point()+geom_line(lwd=1)+ylim(-0.1,3.5)
g1=g1+theme_classic()+xlab("Time")+ylab("")
g1=g1+labs(color=NULL) # Remove legend title (to be manually re-added)
g1=g1+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0)))

g2<-ggplot(d2,aes(x=t,y=z1,col=n_share))
g2=g2+geom_point()+geom_line(lwd=1)+ylim(-0.1,3.5)
g2=g2+theme_classic()+xlab("Time")+ylab("")
g2=g2+labs(color=NULL) # Remove legend title (to be manually re-added)
g2=g2+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0)))

g3<-ggplot(d3,aes(x=t,y=z1,col=n_share))
g3=g3+geom_point()+geom_line(lwd=1)+ylim(-0.1,3.5)
g3=g3+theme_classic()+xlab("Time")+ylab("")
g3=g3+labs(color=NULL) # Remove legend title (to be manually re-added)
g3=g3+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0)))

# Plot fitness against time
h1<-ggplot(d1,aes(x=t,y=w,col=n_share))
h1=h1+geom_point()+geom_line(lwd=1)+ylim(0,0.4)
h1=h1+theme_classic()+xlab("Time")+ylab("")
h1=h1+labs(color=NULL) # Remove legend title (to be manually re-added)
h1=h1+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0)))

h2<-ggplot(d2,aes(x=t,y=w,col=n_share))
h2=h2+geom_point()+geom_line(lwd=1)+ylim(0,0.4)
h2=h2+theme_classic()+xlab("Time")+ylab("")
h2=h2+labs(color=NULL) # Remove legend title (to be manually re-added)
h2=h2+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0)))

h3<-ggplot(d3,aes(x=t,y=w,col=n_share))
h3=h3+geom_point()+geom_line(lwd=1)+ylim(0,0.4)
h3=h3+theme_classic()+xlab("Time")+ylab("")
h3=h3+labs(color=NULL) # Remove legend title (to be manually re-added)
h3=h3+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0)))

ggsave("plot_z1_temp_0.pdf",plot=g1,width=7,height=5)
ggsave("plot_z1_temp_0.9.pdf",plot=g2,width=7,height=5)
ggsave("plot_z1_temp_-0.9.pdf",plot=g3,width=7,height=5)
ggsave("plot_w_temp_0.pdf",plot=h1,width=7,height=5)
ggsave("plot_w_temp_0.9.pdf",plot=h2,width=7,height=5)
ggsave("plot_w_temp_-0.9.pdf",plot=h3,width=7,height=5)


