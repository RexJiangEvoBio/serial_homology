# Simulation of evolution under fluctuating selection
# Optimum is resampled periodically, and the main axis of distribution is aligned with SLLR if applicable
# z1 and z2 in the simulation code are positive trait values

library(mvtnorm) # To sample from multivariate normal distribution
library(resample) # To use colVars
library(ggplot2)

setwd("your_dir")

# Function to calculate fitness
fitness_calc<-function(z,opt,omega){
  Lambda<-solve(omega)
  D2<-t(z-opt)%*%Lambda%*%(z-opt)
  w<-exp(-D2/2)
  return(as.numeric(w))
}

# Function to generate a series of optima over time, with optimum resampled periodically
# Input and output phenotypes are original, not log scale
fluc_wn<-function(n_step,t_cycle,t_fluc_start,opt_ances,mat_fluc){
  out=matrix(0,nrow=n_step,ncol=length(opt_ances))
  out[1,]=log(opt_ances)
  for(t in 1:n_step){
    if((t-t_fluc_start)%%t_cycle==0){
      out[t,]=rmvnorm(1,mean=rep(0,length(opt_ances)),sigma=mat_fluc)[1,]
    }else{
      out[t,]=out[t-1,]
    }
  }
  return(exp(out))
}

# Parameters to stay the same across simulations
u=1e-7 # Mutation rate per cis-element
n_cis=50 # Total number of cis-element deployed in each serial homolog
sig_1=.1;sig_2=.1;sig_p=c(.1,.1) # SD of mutation's effect on b*log(a)
z_a=c(1,1) # Ancestral phenotype (a=1 and b=0.01 for all 50 loci)
width=c(10,10) # Width (SD) of fitness function
n_rep=100 # Number of populations
n_step=200 # Number of time steps per population
t_step=500 # Generations per time step
Ne=1e5 # Focus on a single Ne
t_cycle=20 # Number of time steps between optimum shifts

# Parameter combinations to examine
n_p_all=n_cis*c(0,0.2,0.4,0.6,0.8,1) # Number of pleiotropic elements
r_m_all=c(0.9,0.5) # Mutational correlation for pleiotropic elements
r_selection_all=c(0,0.9,-0.9) # Correlational selection (also applied to optimum movement)
sd_fluc_all=c(10,20) # Magnitude of optimum movement
# Assemble a data matrix containing all parameter combinations
par_all=rep(0,4)
for(i in 1:length(n_p_all)){
  for(j in 1:length(r_m_all)){
    for(k in 1:length(r_selection_all)){
      for(l in 1:length(sd_fluc_all)){
        row=c(n_p_all[i],r_m_all[j],r_selection_all[k],sd_fluc_all[l])
        par_all=rbind(par_all,row)
      }
    }
  }
}
par_all=par_all[-1,]
rownames(par_all)=1:nrow(par_all)
colnames(par_all)=c("n_p","r_m","r_selection","sd_fluc")

z1_var_time=matrix(0,nrow=nrow(par_all),ncol=n_step) # Phenotypic variance of SH1 across populations over time
z2_var_time=matrix(0,nrow=nrow(par_all),ncol=n_step) # Phenotypic variance of SH2 across populations over time
w_mean_time=matrix(0,nrow=nrow(par_all),ncol=n_step) # Mean fitness over time
w_var_time=matrix(0,nrow=nrow(par_all),ncol=n_step) # Variance of fitness over time
cor_time=matrix(0,nrow=nrow(par_all),ncol=n_step) # Correlation between between SHs over time
w_geomean=matrix(0,nrow=nrow(par_all),ncol=n_rep) # Geometric mean of all populations
for(c in 1:nrow(par_all)){
  
  n_p=par_all[c,1] # Read the number pleiotropic cis-elements
  n_1=n_cis-n_p;n_2=n_cis-n_p # The number of non-pleiotropic cis-elements
  V_1=u*(n_1*sig_1^2+n_p*sig_p[1]^2);V_2=u*(n_2*sig_2^2+n_p*sig_p[2]^2) # Mutational variances
  
  r_m=par_all[c,2] # Read mutational correlation for pleiotropic elements
  cov_m=u*n_p*sig_p[1]*sig_p[2]*r_m # Mutational covariance
  M=rbind(c(V_1,cov_m),c(cov_m,V_2)) # M-matrix
  G=2*Ne*M # G-matrix
  
  r_selection=par_all[c,3] # Read correlational selection
  S=rbind(c(width[1]^2,width[1]*width[2]*r_selection),c(width[1]*width[2]*r_selection,width[2]^2)) # Selection matrix (shape of adaptive landscape)
  
  var_fluc=par_all[c,4]^2
  mat_fluc=var_fluc*rbind(c(1,r_selection),c(r_selection,1))
  
  z1_time=matrix(0,nrow=n_rep,ncol=n_step) # z_1 over time
  z2_time=matrix(0,nrow=n_rep,ncol=n_step) # z_2 over time
  w_time=matrix(0,nrow=n_rep,ncol=n_step) # Fitness over time
  for(i in 1:n_rep){
    z=z_a
    # Simulate a trajectory of optimum shift
    z_opt=fluc_wn(n_step=n_step,t_cycle=t_cycle,t_fluc_start=1,opt_ances=z_a,mat_fluc=mat_fluc)
    
    for(t in 1:n_step){
      delta_m=rmvnorm(1,sigma=2*M*t_step) # Neutral evolutionary changes over the time step
      delta_s=-t_step*G%*%solve(S)%*%(log(z)-log(z_opt[t,])) # Evolutionary change due to selection
      z=exp(log(z)+delta_m[1,]+delta_s[,1])
      z1_time[i,t]=z[1];z2_time[i,t]=z[2]
      w_time[i,t]=fitness_calc(log(z),log(z_opt[t,]),S)
    }
    w_geomean[c,i]=(prod(w_time[i,]))^(1/n_step)
  }
  
  z1_var_time[c,]=colVars(log(z1_time));z2_var_time[c,]=colVars(log(z2_time))
  w_mean_time[c,]=colMeans(w_time)
  w_var_time[c,]=colVars(w_time)
  
  for(t in 1:n_step){
    cor_time[c,t]=cor(log(z1_time[,t]),log(z2_time[,t]))
  }
}

write.table(data.frame(par_all,z1_var_time),file="var_time_1_fluc_wn.txt",sep="\t")
write.table(data.frame(par_all,z2_var_time),file="var_time_2_fluc_wn.txt",sep="\t")
write.table(data.frame(par_all,w_mean_time),file="w_mean_time_fluc_wn.txt",sep="\t")
write.table(data.frame(par_all,w_var_time),file="w_var_time_fluc_wn.txt",sep="\t")
write.table(data.frame(par_all,cor_time),file="cor_time_fluc_wn.txt",sep="\t")
write.table(data.frame(par_all,w_geomean),file="w_geomean_fluc_wn.txt",sep="\t")

# Parameter combinations to plot (fluctuation speed, correlational selection)
par_plot_all=rbind(c(10,0),
                   c(10,0.9),
                   c(10,-0.9),
                   c(20,0),
                   c(20,0.9),
                   c(20,-0.9)
)
# Fitness over time
d<-read.table("w_mean_time_fluc_wn.txt",sep="\t")
d$n_p=d$n_p/max(d$n_p) # Convert to fractions
d=d[which(d$r_m==0.9),] # Focus on scenarios where the regulators have similar binding preference
for(row in 1:nrow(par_plot_all)){
  par_plot=par_plot_all[row,]
  d_sub=d[which(d$sd_fluc==par_plot[1]&d$r_selection==par_plot[2]),]
  # Rearrange into a format for ggplot
  d_new=cbind(rep(d_sub$n_p[1],ncol(d_sub)-4),1:(ncol(d_sub)-4),as.numeric(d_sub[1,5:ncol(d_sub)]))
  for(i in 2:length(unique(d_sub$n_p))){
    d_new=rbind(d_new,cbind(rep(d_sub$n_p[i],ncol(d_sub)-4),1:(ncol(d_sub)-4),as.numeric(d_sub[i,5:ncol(d_sub)])))
  }
  colnames(d_new)=c("n_p","t","w")
  d_new=data.frame(d_new)
  
  d_new$n_p=factor(d_new$n_p,levels=sort(unique(d_new$n_p)),ordered=TRUE)
  g<-ggplot(d_new,aes(x=t,y=w,colour=n_p))
  g=g+geom_point()+geom_line(lwd=1)+ylim(0,1)#+geom_hline(yintercept=log(10),color="red")
  g=g+theme_classic()+xlab("")+ylab("")
  g=g+labs(color=NULL) # Remove legend title (to be manually re-added)
  g=g+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0))) # Replace dots in the legend with boxes
  fn=paste("plot_w_fluc_wn_",par_plot[2],"_",par_plot[1],".pdf",sep="")
  ggsave(fn,plot=g,width=7,height=5)
}

# Fitness variance over time
d<-read.table("w_var_time_fluc_wn.txt",sep="\t")
d$n_p=d$n_p/max(d$n_p) # Convert to fractions
d=d[which(d$r_m==0.9),] # Focus on scenarios where the regulators have similar binding preference
for(row in 1:nrow(par_plot_all)){
  par_plot=par_plot_all[row,]
  d_sub=d[which(d$sd_fluc==par_plot[1]&d$r_selection==par_plot[2]),]
  # Rearrange into a format for ggplot
  d_new=cbind(rep(d_sub$n_p[1],ncol(d_sub)-4),1:(ncol(d_sub)-4),as.numeric(d_sub[1,5:ncol(d_sub)]))
  for(i in 2:length(unique(d_sub$n_p))){
    d_new=rbind(d_new,cbind(rep(d_sub$n_p[i],ncol(d_sub)-4),1:(ncol(d_sub)-4),as.numeric(d_sub[i,5:ncol(d_sub)])))
  }
  colnames(d_new)=c("n_p","t","vw")
  d_new=data.frame(d_new)
  
  # Set upper limit of plot based on fluctuation speed
  ub=max(d[which(d$sd_fluc==par_plot[1]),5:ncol(d)])*1.1
  
  d_new$n_p=factor(d_new$n_p,levels=sort(unique(d_new$n_p)),ordered=TRUE)
  g<-ggplot(d_new,aes(x=t,y=vw,colour=n_p))
  g=g+geom_point()+geom_line(lwd=1)+ylim(0,ub)#+geom_hline(yintercept=log(10),color="red")
  g=g+theme_classic()+xlab("")+ylab("")
  g=g+labs(color=NULL) # Remove legend title (to be manually re-added)
  g=g+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0))) # Replace dots in the legend with boxes
  fn=paste("plot_vw_fluc_wn_",par_plot[2],"_",par_plot[1],".pdf",sep="")
  ggsave(fn,plot=g,width=7,height=5)
}

# Phenotypic variance over time and phylogenetic signal
d1<-read.table("var_time_1_fluc_wn.txt",sep="\t")
d2<-read.table("var_time_2_fluc_wn.txt",sep="\t")
d1$n_p=d1$n_p/max(d1$n_p) # Convert to fractions
d2$n_p=d2$n_p/max(d2$n_p) # Convert to fractions
# Focus on scenarios where the regulators have similar binding preference
d1=d1[which(d1$r_m==0.9),]
d2=d2[which(d2$r_m==0.9),]
physig_out=matrix(0,nrow=nrow(d1),ncol=2);colnames(physig_out)=c("physig_z1","physig_z2")
for(i in 1:nrow(d1)){
  physig_out[i,1]=cor(1:(ncol(d1)-4),as.numeric(d1[i,5:ncol(d1)]))
  physig_out[i,2]=cor(1:(ncol(d1)-4),as.numeric(d2[i,5:ncol(d1)]))
}
physig_out=data.frame(d1[,1:4],physig_out)
write.table(physig_out,file="physig_fluc_wn.txt",sep="\t")
physig_out$sd_fluc=factor(physig_out$sd_fluc,levels=sort(unique(physig_out$sd_fluc)),ordered=TRUE)
for(r in c(0,0.9,-0.9)){
  d_sub=physig_out[which(physig_out$r_selection==r),]
  
  g1=ggplot(d_sub,aes(x=n_p,y=physig_z1,colour=sd_fluc))
  g1=g1+geom_point()+geom_line(lwd=1)
  g1=g1+ylim(0,1)+theme_classic()+xlab("")+ylab("")+scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))
  g1=g1+labs(color=NULL) # Remove legend title (to be manually re-added)
  g1=g1+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0))) # Replace dots in the legend with boxes
  fn1=paste("plot_physig_fluc_z1_",r,".pdf",sep="")
  ggsave(fn1,plot=g1,width=5,height=4)
  
  g2=ggplot(d_sub,aes(x=n_p,y=physig_z2,colour=sd_fluc))
  g2=g2+geom_point()+geom_line(lwd=1)
  g2=g2+ylim(0,1)+theme_classic()+xlab("")+ylab("")+scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))
  g2=g2+labs(color=NULL) # Remove legend title (to be manually re-added)
  g2=g2+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0))) # Replace dots in the legend with boxes
  fn2=paste("plot_physig_fluc_z2_",r,".pdf",sep="")
  ggsave(fn2,plot=g2,width=5,height=4)
}

# Phenotypic correlation over time
d<-read.table("cor_time_fluc_wn.txt",sep="\t")
d$n_p=d$n_p/max(d$n_p) # Convert to fractions
d=d[which(d$r_m==0.9),] # Focus on scenarios where the regulators have similar binding preference
for(row in 1:nrow(par_plot_all)){
  par_plot=par_plot_all[row,]
  d_sub=d[which(d$sd_fluc==par_plot[1]&d$r_selection==par_plot[2]),]
  # Rearrange into a format for ggplot
  d_new=cbind(rep(d_sub$n_p[1],ncol(d_sub)-4),1:(ncol(d_sub)-4),as.numeric(d_sub[1,5:ncol(d_sub)]))
  for(i in 2:length(unique(d_sub$n_p))){
    d_new=rbind(d_new,cbind(rep(d_sub$n_p[i],ncol(d_sub)-4),1:(ncol(d_sub)-4),as.numeric(d_sub[i,5:ncol(d_sub)])))
  }
  colnames(d_new)=c("n_p","t","corr")
  # Keep time points that are multiples of pre-set window size
  d_new=data.frame(d_new)
  d_new=d_new[which(d_new$t>=10),] # Starting from the 10th step, omitting early steps when there was little variation
  
  d_new$n_p=factor(d_new$n_p,levels=sort(unique(d_new$n_p)),ordered=TRUE)
  g<-ggplot(d_new,aes(x=t,y=corr,colour=n_p))
  g=g+geom_point()+geom_line(lwd=1)+xlim(0,200)+ylim(-1,1)#+geom_hline(yintercept=log(10),color="red")
  g=g+theme_classic()+xlab("")+ylab("")
  g=g+labs(color=NULL) # Remove legend title (to be manually re-added)
  g=g+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0))) # Replace dots in the legend with boxes
  fn=paste("plot_cor_fluc_wn_",par_plot[2],"_",par_plot[1],".pdf",sep="")
  ggsave(fn,plot=g,width=7,height=5)
}

# Analyze geometric mean fitness
d<-read.table("w_geomean_fluc_wn.txt",sep="\t")
d$n_p=d$n_p/max(d$n_p) # Convert to fractions
d=d[which(d$r_m==0.9),] # Focus on scenarios where the regulators have similar binding preference
signif=matrix(0,nrow=nrow(par_plot_all),ncol=2)
for(row in 1:nrow(par_plot_all)){
  par_plot=par_plot_all[row,]
  d_sub=d[which(d$sd_fluc==par_plot[1]&d$r_selection==par_plot[2]),]
  # Rearrange into a format for ggplot
  d_new=cbind(rep(d_sub$n_p[1],ncol(d_sub)-4),as.numeric(d_sub[1,5:ncol(d_sub)]))
  for(i in 2:length(unique(d_sub$n_p))){
    d_new=rbind(d_new,cbind(rep(d_sub$n_p[i],ncol(d_sub)-4),as.numeric(d_sub[i,5:ncol(d_sub)])))
  }
  colnames(d_new)=c("n_p","geow")
  d_new=data.frame(d_new)
  signif[row,1]=cor(d_new$n_p,d_new$geow)
  signif[row,2]=cor.test(d_new$n_p,d_new$geow)$p.value
  
  g<-ggplot(d_new,aes(x=n_p,y=geow))
  g=g+geom_point()+geom_smooth(method="loess")+scale_x_continuous(breaks =c(0,.2,.4,.6,.8,1))+ylim(0,1)
  g=g+theme_classic()+xlab("")+ylab("")
  fn=paste("plot_geo_fluc_wn_",par_plot[2],"_",par_plot[1],".pdf",sep="")
  ggsave(fn,plot=g,width=7,height=5)
}
signif=cbind(par_plot_all,signif)
colnames(signif)=c("fluc_rate","selection_cor","cor","pv")
write.table(data.frame(signif),file="w_geomean_fluc_wn_cor.txt")

