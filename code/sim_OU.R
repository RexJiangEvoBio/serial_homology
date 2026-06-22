# Simulate the evolution of population mean phenotypes under selection 
# Evolution modeled as an OU process
# Mutational parameters are derived from a developmental model
# z1 and z2 in the simulation code are positive trait values

library(mvtnorm) # To sample from multivariate normal distribution
library(resample) # To use colVars
library(ggplot2)

setwd("your_dir")

# Function to calculate fitness
fitness_calc<-function(z,opt,omega) {
  Lambda<-solve(omega)
  D2<-t(z-opt)%*%Lambda%*%(z-opt)
  w<-exp(-D2/2)
  return(as.numeric(w))
}

# Parameters to stay the same across simulations
u=1e-7 # Mutation rate per cis-element
n_cis=50 # Total number of cis-element deployed in each serial homolog
sig_1=.1;sig_2=.1;sig_p=c(.1,.1) # SD of mutation's effect on b*log(a)
z_a=c(1,1) # Ancestral phenotype (a=1 and b=0.01 for all 50 loci)
width=c(10,10) # Width (SD) of fitness function
n_rep=20 # Number of populations
n_step=100 # Number of time steps per population
t_step=500 # Generations per time step

# Parameter combinations to examine
Ne_all=10^c(3,4,5) # Effective population sizes
n_p_all=n_cis*c(0,0.2,0.4,0.6,0.8,1) # Number of pleiotropic elements
r_m_all=c(0.9,0.5) # Mutational correlation for pleiotropic elements
z_opt_all=list(c(10*z_a[1],z_a[2]),10*z_a,c(10*z_a[1],0.1*z_a[2]),z_a) # Optimum (z_1 under directional selection, concordant selection, discordant selection, stabilizing selection)
r_selection_all=c(0,0.5,0.9,-0.5,-0.9) # Correlational selection
# Assemble a data matrix containing all parameter combinations
par_all=rep(0,6)
for(i in 1:length(Ne_all)){
  for(j in 1:length(n_p_all)){
    for(k in 1:length(r_m_all)){
      for(l in 1:length(z_opt_all)){
        for(m in 1:length(r_selection_all)){
          row=c(Ne_all[i],n_p_all[j],r_m_all[k],z_opt_all[[l]][1],z_opt_all[[l]][2],r_selection_all[m])
          par_all=rbind(par_all,row)
        }
      }
    }
  }
}
par_all=par_all[-1,]
rownames(par_all)=1:nrow(par_all)
colnames(par_all)=c("Ne","n_p","r_m","z1_opt","z2_opt","r_selection")

# Data matrices to store output (all using log-transformed phenotypes)
z1_mean_time=matrix(0,nrow=nrow(par_all),ncol=n_step) # Mean z1 across populations over time
z2_mean_time=matrix(0,nrow=nrow(par_all),ncol=n_step) # Mean z2 across populations over time
z1_var_time=matrix(0,nrow=nrow(par_all),ncol=n_step) # Variance of z1 across populations over time
z2_var_time=matrix(0,nrow=nrow(par_all),ncol=n_step) # Variance of z2 across populations over time
out_end_sum=matrix(0,nrow=nrow(par_all),ncol=4);colnames(out_end_sum)=c("ln_z1","ln_z2","var1","var_2") # End-point means and variances
# Simulate
for(c in 1:nrow(par_all)){
  Ne=par_all[c,1] # Read Ne
  
  n_p=par_all[c,2] # Read the number pleiotropic cis-elements
  n_1=n_cis-n_p;n_2=n_cis-n_p # The number of non-pleiotropic cis-elements
  V_1=u*(n_1*sig_1^2+n_p*sig_p[1]^2);V_2=u*(n_2*sig_2^2+n_p*sig_p[2]^2) # Mutational variances
  
  r_m=par_all[c,3] # Read mutational correlation for pleiotropic elements
  cov_m=u*n_p*sig_p[1]*sig_p[2]*r_m # Mutational covariance
  M=rbind(c(V_1,cov_m),c(cov_m,V_2)) # M-matrix
  G=2*Ne*M # G-matrix
  
  z_opt=par_all[c,4:5] # Read optimal phenotype
  r_selection=par_all[c,6] # Read correlational selection
  S=rbind(c(width[1]^2,width[1]*width[2]*r_selection),c(width[1]*width[2]*r_selection,width[2]^2)) # Selection matrix (shape of adaptive landscape)
  
  z1_time=matrix(0,nrow=n_rep,ncol=n_step) # z_1 over time
  z2_time=matrix(0,nrow=n_rep,ncol=n_step) # z_2 over time
  z_end=matrix(0,nrow=n_rep,ncol=2) # End-point phenotypes
  for(i in 1:n_rep){
    z=z_a
    for(t in 1:n_step){
      delta_m=rmvnorm(1,sigma=2*M*t_step) # Neutral evolutionary changes over the time step
      delta_s=-t_step*G%*%solve(S)%*%(log(z)-log(z_opt)) # Evolutionary change due to selection
      z=exp(log(z)+delta_m[1,]+delta_s[,1])
      z1_time[i,t]=z[1];z2_time[i,t]=z[2]
    }
    z_end[i,]=z
  }
  
  z1_mean_time[c,]=colMeans(log(z1_time));z2_mean_time[c,]=colMeans(log(z2_time))
  z1_var_time[c,]=colVars(log(z1_time));z2_var_time[c,]=colVars(log(z2_time))
  out_end_sum[c,]=c(colMeans(log(z_end)),colVars(log(z_end)))
  
}

write.table(data.frame(par_all,out_end_sum),file="out_end_sum.txt",sep="\t")
write.table(data.frame(par_all,z1_mean_time),file="mean_time_1.txt",sep="\t")
write.table(data.frame(par_all,z2_mean_time),file="mean_time_2.txt",sep="\t")
write.table(data.frame(par_all,z1_var_time),file="var_time_1.txt",sep="\t")
write.table(data.frame(par_all,z2_var_time),file="var_time_2.txt",sep="\t")

# Plot z1 (log-scale) against time for selected scenarios of directional selection
d<-read.table("mean_time_1.txt",sep="\t")
d=d[which(d$Ne==1e5),] # Pick the Ne to focus on
d$n_p=d$n_p/max(d$n_p) # Convert to fractions
# Extract the subset to plot (z1 optimum, z2 optimum, correlational selection, mutational correlation at pleiotropic loci)
par_plot_all=rbind(c(10,1,0,0.9),
                   c(10,10,0,0.9),
                   c(10,0.1,0,0.9),
                   c(10,10,0.9,0.9),
                   c(10,0.1,-0.9,0.9),
                   c(10,1,0,0.5),
                   c(10,10,0,0.5),
                   c(10,0.1,0,0.5),
                   c(10,10,0.9,0.5),
                   c(10,0.1,-0.9,0.5)
                   )

for(row in 1:nrow(par_plot_all)){
  par_plot=par_plot_all[row,]
  d_sub=d[which(d$z1_opt==par_plot[1]&d$z2_opt==par_plot[2]&d$r_selection==par_plot[3]&d$r_m==par_plot[4]),]
  # Rearrange into a format for ggplot
  d_new=cbind(rep(d_sub$n_p[1],ncol(d_sub)-6),1:(ncol(d_sub)-6),as.numeric(d_sub[1,7:ncol(d_sub)]))
  for(i in 2:length(unique(d_sub[,2]))){
    d_new=rbind(d_new,cbind(rep(d_sub$n_p[i],ncol(d_sub)-6),1:(ncol(d_sub)-6),as.numeric(d_sub[i,7:ncol(d_sub)])))
  }
  colnames(d_new)=c("n_p","t","z1")
  d_new=data.frame(d_new)

  d_new$n_p=factor(d_new$n_p,levels=sort(unique(d_new$n_p)),ordered=TRUE)
  g<-ggplot(d_new,aes(x=t,y=z1,colour=n_p))
  g=g+geom_point()+geom_line(lwd=1)+ylim(-0.05,2.5)+geom_hline(yintercept=log(10),color="red")
  g=g+theme_classic()+xlab("Time")+ylab("")
  g=g+labs(color=NULL) # Remove legend title (to be manually re-added)
  g=g+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0))) # Replace dots in the legend with boxes

  fn=paste("plot_z1_",log10(par_plot[1]),"_",log10(par_plot[2]),"_",par_plot[3],"_",par_plot[4],".pdf",sep="")
  ggsave(fn,plot=g,width=7,height=5)
}

# Plot variance across populations against time (mainly for scenarios of stabilizing selection)
d<-read.table("var_time_1.txt",sep="\t")
d=d[which(d$z1_opt==1&d$z2_opt==1),] # Stabilizing selection
d$n_p=d$n_p/max(d$n_p) # Convert to fractions
# Extract the subset to plot (correlational selection, mutational correlation at pleiotropic loci)
par_plot_all=rep(0,3)
for(i in 1:length(unique(d$Ne))){
  for(j in 1:length(unique(d$r_selection))){
    for(k in 1:length(unique(d$r_m))){
      par_plot_all=rbind(par_plot_all,c(unique(d$Ne)[i],unique(d$r_selection)[j],unique(d$r_m)[k]))
    }
  }
}
par_plot_all=par_plot_all[-1,]

for(row in 1:nrow(par_plot_all)){
  par_plot=par_plot_all[row,]
  d_sub=d[which(d$Ne==par_plot[1]&d$r_selection==par_plot[2]&d$r_m==par_plot[3]),]
  # Rearrange into a format for ggplot
  d_new=cbind(rep(d_sub$n_p[1],ncol(d_sub)-6),1:(ncol(d_sub)-6),as.numeric(d_sub[1,7:ncol(d_sub)]))
  for(i in 2:length(unique(d_sub$n_p))){
    d_new=rbind(d_new,cbind(rep(d_sub$n_p[i],ncol(d_sub)-6),1:(ncol(d_sub)-6),as.numeric(d_sub[i,7:ncol(d_sub)])))
  }
  colnames(d_new)=c("n_p","t","var1")
  d_new=data.frame(d_new)

  d_new$n_p=factor(d_new$n_p,levels=sort(unique(d_new$n_p)),ordered=TRUE)
  g<-ggplot(d_new,aes(x=t,y=log(var1),colour=n_p))
  g=g+geom_point()+geom_line(lwd=1)#+ylim(-0.05,2.5)#+geom_hline(yintercept=log(10),color="red")
  g=g+theme_classic()+xlab("Time")+ylab("")
  g=g+labs(color=NULL) # Remove legend title (to be manually re-added)
  g=g+guides(color=guide_legend(override.aes=list(shape=15,size=6,linetype=0))) # Replace dots in the legend with boxes
  fn=paste("plot_var1_0_0_",par_plot[1],"_",par_plot[2],".pdf",sep="")
  ggsave(fn,plot=g,width=7,height=5)
}

