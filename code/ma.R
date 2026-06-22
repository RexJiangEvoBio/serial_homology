# Sanity check: see if the distribution of mutant phenotypes are as expected using MA simulations
# Simulation is agnostic to ploidy (given no dominance & ancestral phenotypic value is always 0)

setwd("your_dir")

library(mvtnorm) # To sample from multivariate normal distribution
library(ggplot2)

# Simulate a single trait under
n_cis=50 # Number of cis-elements affecting the trait
a_ances=rep(1,n_cis) # cis-trans binding affinity
b_ances=0.1 # Each cis-element's effect on the phenotype
z_ances=sum(b_ances*log(a_ances)) # Ancestral phenotype
sig=1 # SD of mutation's effect on log-transformed binding affinity
n_line=500 # Number of MA lines to simulate
n_mut=100 # Number of mutations per MA line
z_out=rep(0,n_line) # Vector to store mutant phenotypes at the end
for(i in 1:n_line){
  a=a_ances # Initialize the mutant binding affinities
  for(j in 1:n_mut){
    affect=sample(1:n_cis,1) # Randomly pick a cis-element to be the target of the mutation
    a[affect]=a[affect]*exp(rnorm(1,mean=0,sd=sig)) # Binding affinity in mutant
  }
  z_out[i]=sum(b_ances*log(a)) # Mutant phenotype at the end
}
hist(z_out,breaks=50) # Check distribution
var_exp=n_mut*b_ances^2*sig^2 # Expected phenotypic variance (number of mutations)
var(z_out) # Compare observed and expected variances

# Simulate 2 traits (z1, z2)
n_cis=50 # Number of cis-elements affecting each trait
n_p=30 # Number of pleiotropic cis-elements
n_np_1=n_cis-n_p;n_np_2=n_cis-n_p # Number of non-pleiotropic cis-elements affecting each trait
a_ances_1=rep(1,n_np_1);a_ances_2=rep(1,n_np_2);a_ances_p=matrix(1,nrow=2,ncol=n_p) # Ancestral binding affinity
b_ances=.1 # Each cis-element's effect on the phenotype
sig=1 # SD of mutation's effect on log-transformed binding affinity
r_m=0.9;m=sig*(rbind(c(1,r_m),c(r_m,1))) # Correlation between pleiotropic mutation's effect on a1 and a2
n_line=500 # Number of MA lines to simulate
n_mut=100 # Number of mutations per MA line
z_out=matrix(0,nrow=n_line,ncol=2)
for(i in 1:n_line){
  # Initialize the mutant binding affinities
  a_1=a_ances_1
  a_2=a_ances_2
  a_p=a_ances_p
  for(j in 1:n_mut){
    type=sample(1:3,1,prob=c(n_np_1,n_np_2,n_p)/n_cis) # Decide which type of cis-element is affected
    if(type==3){ # Mutation is in a pleiotropic cis-element
      affect=sample(1:n_p,1) # Randomly pick a cis-element to be the target of the mutation
      effect=rmvnorm(1,mean=c(0,0),sigma=m) # Mutation's effect
      a_p[,affect]=a_p[,affect]*exp(effect[1,]) # Mutant binding affinity
    }else{ # Mutation is in a non-pleiotropic cis-element
      if(type==1){ # z1 is affected
        affect=sample(1:n_np_1,1) # Randomly pick a cis-element to be the target of the mutation
        a_1[affect]=a_1[affect]*exp(rnorm(1,mean=0,sd=sig)) # Mutant binding affinity
      }else{ # z2 is affected
        affect=sample(1:n_np_2,1) # Randomly pick a cis-element to be the target of the mutation
        a_2[affect]=a_2[affect]*exp(rnorm(1,mean=0,sd=sig)) # Mutant binding affinity
      }
    }
  }
  # Mutant phenotype at the end
  z_out[i,1]=sum(b_ances*c(log(a_1),log(a_p[1,])))
  z_out[i,2]=sum(b_ances*c(log(a_2),log(a_p[2,])))
}
plot(z_out[,1],z_out[,2]) # Check distribution
cov_exp=(r_m*n_p/n_cis)*n_mut*b_ances^2*sig^2 # Expected covariance between z1 and z2
cov_exp;cov(z_out) # Compare observed and expected covariances
cor_exp=r_m*n_p/n_cis # Expected correlation between z1 and z2
cor_exp;cor(z_out) # Compare observed and expected correlations

# Go through mutational parameters
# Assemble a data matrix containing combinations to examine
n_p_all=c(0,10,20,30,40,50) # Number of pleiotropic cis-elements
r_m_all=c(.9,.5) # Correlation between pleiotropic mutation's effect on a1 and a2
par_all=c(0,0)
for(i in 1:length(r_m_all)){
  for(j in 1:length(n_p_all)){
    par_all=rbind(par_all,c(n_p_all[j],r_m_all[i]))
  }
}
par_all=par_all[-1,]

n_cis=50 # Number of cis-elements affecting each trait
b_ances=.1 # Each cis-element's effect on the phenotype
sig=1 # Ancestral binding affinity
n_line=200 # Number of MA lines to simulate
n_mut=5 # Number of mutations per cis-element per MA line
z_out_all=list() # List to store output data frames
for(c in 1:nrow(par_all)){
  n_p=par_all[c,1];r_m=par_all[c,2]
  n_np_1=n_cis-n_p;n_np_2=n_cis-n_p
  n_mut_total=n_mut*(n_p+n_np_1+n_np_2) # Total number of mutations (proportional to the total number of cis-elements present)
  a_ances_1=rep(1,n_np_1);a_ances_2=rep(1,n_np_2);a_ances_p=matrix(1,nrow=2,ncol=n_p)
  m=sig*(rbind(c(1,r_m),c(r_m,1)))
  z_out=matrix(0,nrow=n_line,ncol=2)
  for(i in 1:n_line){
    # Initialize the mutant binding affinities
    a_1=a_ances_1
    a_2=a_ances_2
    a_p=a_ances_p
    for(j in 1:n_mut_total){
      if(n_p>0){ # There exist pleiotropic cis-elements
        if(n_p==n_cis){ # All cis-elements are pleiotropic
          type=3 # Mutation is pleiotropic
        }else{
          type=sample(1:3,1,prob=c(n_np_1,n_np_2,n_p)/(n_np_1+n_np_2+n_p)) # Decide which type of cis-element is affected
        }
      }else{
        type=sample(1:2,1,prob=c(n_np_1,n_np_2)/(n_np_1+n_np_2)) # Decide which type of cis-element is affected
      }
      if(type==3){ # Mutation is pleiotropic
        affect=sample(1:n_p,1) # Randomly pick a cis-element to be the target of the mutation
        effect=rmvnorm(1,mean=c(0,0),sigma=m) # Mutation's effect sampled from a bivariate distribution
        a_p[,affect]=a_p[,affect]*exp(effect[1,]) # Mutant binding affinity
      }else{ # Mutation isn't pleiotropic
        if(type==1){
          affect=sample(1:n_np_1,1) # Randomly pick a cis-element to be the target of the mutation
          a_1[affect]=a_1[affect]*exp(rnorm(1,mean=0,sd=sig)) # Mutant binding affinity
        }else{
          affect=sample(1:n_np_2,1) # Randomly pick a cis-element to be the target of the mutation
          a_2[affect]=a_2[affect]*exp(rnorm(1,mean=0,sd=sig)) # Mutant binding affinity
        }
      }
    }
    # Mutant phenotype at the end
    z_out[i,1]=sum(b_ances*c(log(a_1),log(a_p[1,])))
    z_out[i,2]=sum(b_ances*c(log(a_2),log(a_p[2,])))
  }
  z_out_all[[c]]=data.frame(rep(n_p,n_line),rep(r_m,n_line),z_out)
}

# Combine into a single dataset
d=z_out_all[[1]]
for(c in 2:nrow(par_all)){
  d=rbind(d,z_out_all[[c]])
}
colnames(d)=c("np","rm","z1","z2")
write.table(d,file="z_out_MA.txt",sep="\t")

# Calculate mutational (co)variances from mutant phenotypes
m_out=matrix(0,nrow=nrow(par_all),ncol=8)
for(c in 1:nrow(par_all)){
  dsub=d[which(d$np==par_all[c,1]&d$rm==par_all[c,2]),]
  m_out[c,1]=var(dsub$z1) # z1 variance
  m_out[c,2]=var(dsub$z2) # z2 variance
  m_out[c,3]=cov(dsub$z1,dsub$z2) # Covariance
  m_out[c,4]=cor(dsub$z1,dsub$z2) # Correlation
  m_out[c,5]=cor.test(dsub$z1,dsub$z2)$p.value # Significance of correlation
  
  m_out[c,6]=n_mut*n_cis*b_ances^2*sig^2 # Expected variance (the same for z1 and z2)
  m_out[c,7]=n_mut*par_all[c,1]*par_all[c,2]*b_ances^2*sig^2 # Expected covariance
  m_out[c,8]=(par_all[c,1]/n_cis)*par_all[c,2] # Expected correlation
}
m_out=data.frame(par_all,m_out)
colnames(m_out)=c("np","rm","v1","v2","cov","cor","p_cor","v_exp","cov_exp","cor_exp");rownames(m_out)=NULL
write.table(m_out,file="m_mat_out.txt",sep="\t")

# Plot (multiple scatter plots in one figure)
g<-ggplot(d,aes(x=z1,y=z2))
g=g+geom_point()+facet_grid(rows=vars(np),cols=vars(rm),axes="all")
g=g+geom_smooth(method="lm")+theme_classic()
g=g+theme(strip.background=element_blank(),strip.text.y=element_text(angle=0))
ggsave("plot_MA.pdf",plot=g,width=10,height=20)


