## R version of the community model 
library(deSolve)
library(tidyverse)
library(igraph)
library(ggnetwork)
library(GGally)
library(intergraph)
community_pop <- 1000
init_conds <- c(S = community_pop-1, I = 1, R = 0)
print(init_conds)
parms <- data.frame(bet =0.00035 , gam = 1/14) 
print(parms)
times <- seq(from = 1, to = 144, by = 1)
print(times)


SIR_community_model <- function(t, x, parms) {
    S <- x[1]
    I <- x[2]
    R <- x[3]

    with(parms, {
    dS <- -bet*S*I
    dI <- bet*S*I - gam*I
    dR <- gam*I

    dt <- c(dS,dI,dR)
    return(list(dt))
})}

Community_output <- data.frame(lsoda(y = init_conds, func = SIR_community_model,times = times, parms=parms))

Community_output %>% pivot_longer(cols = !time) %>% arrange(desc(time))%>% 
  ggplot(aes(x=time,y =value, color = name))+geom_line()+
  scale_color_manual(name=NULL,values=c("blue","red","purple"),breaks = c("S","I","R"))+
  theme_classic()+
  labs(x= "Time (days)", y = "Number of individuals")+ggtitle("Community outbreak of pathogen")
ggplot(Community_output)+geom_line(aes(x=time,y = S), color = "blue")+
 geom_line(aes(x=time, y=I),color = "red") +
 geom_line(aes(x=time, y = R), color = "purple")+scale_color_manual()+theme_classic()

ggplot(Community_output)+geom_line(aes(x=time, y = I/(S+I+R)))
##Okay next part is to use the solved numeric output from the community model for the initial conditions for the building model
N <-20 #number of rooms
m <- 4 #number of equations for each room
Maxtime <- 8 # max time we are looking at (8 for an 8 hour work day)
delt <- 0.1 # proportionality constant ( what proportion of individuals from the community are in the building of interest)
day <- 31
Building_ICs <- Community_output[day,]

Sb <- Community_output$S[day]*delt
Ib <- Community_output$I[day]*delt
Rb <- Community_output$R[day]*delt
Rb
Sr <- runif(N, min = 0, max = 1)
Ir <- runif(N, min = 0, max = 1)
Rr <- runif(N, min = 0, max = 1)

Sr <- (Sr/sum(Sr))*Sb
Ir <- (Ir/sum(Ir))*Ib
Rr <- (Rr/sum(Rr))*Rb
Pr <- c(rep(0,N))
#initial conditions for each room 
Init_conds <-c(S=Sr, I = Ir, R = Rr, P = Pr)


#create a random network with N patches and fixed connectivity 
set.seed(13156166)
Net <- sample_gnp(N, 0.5, directed = FALSE, loops = FALSE)
plot(Net)
Dis <- distances(Net)
Dis
Net<-set_vertex_attr(Net, "Infectious", value = Ir)
Net<-set_vertex_attr(Net, "Room", value = seq(1:N))
plot(Net)
ggnetwork(Net)
set.seed(20)
ggplot(Net, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = Infectious), color = "#C21807")+
  theme_blank()+
  scale_size(name = "Number of\nInfectious people\nin rooms")+ggtitle("Network representation of a building")
#create transition matrix

ind_mov <- data.frame(matrix(runif(N^2), nrow = N))
diag(ind_mov) <- 0
ind_mov

theta_mov <- data.frame(matrix(runif(N^2), nrow = N))
diag(theta_mov) <- 0

t_mov <- 0.1
#fraction of people that actually move
i_mov <- 0.2


#making it so only a proportion of people actually move
#first have to make sure that people leaving rooms is not more than 1 -- normalized by rows
#then multiplied by the proportion we want to actually mov
ind_mov_norm <- t(apply(ind_mov, 1, function(x) x / sum(x)))
ind_mov <- i_mov*ind_mov_norm
ind_mov


theta_mov_norm <- t(apply(theta_mov, 1, function(x) x / sum(x)))
theta_mov <- t_mov*theta_mov_norm
theta_mov

#Next we need to make sure that people can only move between rooms that are actually connected
#will need to look at the distance matrix of the network for that.

for(i in 1:ncol(mov)){
  for( j in 1:nrow(mov)){
    if(Dis[i,j] > 1){
      mov[i,j] <- 0
    }
  }
}

##Okay, now we have our network set up, we have a transition matrix, we have the initial conditions 
#for how many people are in each room of each class. Now comes the time to solve the ODE

#First, define necessary parameters
s <- 10 #shedding rate
a <- 0.005 #absorption
#theta <- 0 #particles coming in through other rooms
d <- 0.6 # decay of particles
K <- 100 #carrying capacity

times <- seq(from = 0, to = Maxtime, by = 0.1)

parms <- c(s=s,a=a, d=d, K = K)

flux_of_people <- function(N, ind_mov, State){
  all_room_change <- c(seq(N))
  for(j in 1:N){
    flux_in_temp <- 0
    flux_out_temp <- 0
    for(i in 1:N){
      flux_in_temp <- flux_in_temp + (ind_mov[i,j]*State[i])
      flux_out_temp <- flux_out_temp + (ind_mov[j,i]*State[j])
    }
    all_room_change[j] <- flux_in_temp - flux_out_temp
  }
  test <- c(M = as.vector(all_room_change))
  return(test)
}



Other_room_particle_sums <- function(State,theta){
  particle_room_change <- c(seq(N))
  for(j in 1:N){
    particles_in_temp <- 0
    for(i in 1:N){
      if(i == j){
        particles_in_temp <- particles_in_temp
      }else
      particles_in_temp <- particles_in_temp + (theta[i,j]*State[i])
      
    }
    all_room_change[j] <- particles_in_temp
  }
  particle_change <- c(M = as.vector(all_room_change))
  return(particle_change)
}

Particle_model <- function(t, x, parms){
  ncompartment <- 4
  n_rooms <- length(x)/ncompartment
  S <- as.matrix(x[1:n_rooms])
  I <- as.matrix(x[(n_rooms+1):(2*n_rooms)])
  R <- as.matrix(x[(2*n_rooms+1):(3*n_rooms)])
  P <- as.matrix(x[(3*n_rooms+1):(4*n_rooms)])

         dS <- as.matrix(flux_of_people(N, ind_mov, State = S))
         dI <- as.matrix(flux_of_people(N, ind_mov, State = I))
         dR <- as.matrix(flux_of_people(N, ind_mov, State = R))
         #last step will be the particle EQ
         dP <- as.matrix(s*as.vector(I)) - as.matrix(a*P)*(1-(as.matrix(P/K)))*(S+I+R)+ as.matrix(Other_room_particle_sums(State = P,theta_mov)) - as.matrix(d*P)
       
         dt <- c(dS,dI,dR,dP)
         return(list(dt))

}

#sum of outgoing thetas must be less than d

Building_output <- data.frame(lsoda(y = Init_conds, func = Particle_model,times = times, parms = parms))
ggplot(Building_output)+geom_line(aes(x=time, y = P1))+geom_line(aes(x=time, y = P2))+ geom_line(aes(x=time, y =P3))+geom_line(aes(x = time, y = P1+P2+P3+P4+P5+P6+P7+P8+P9+P10))

ggplot(Building_output)+geom_line(aes(x=time, y = P1+P2+P3))+geom_line(aes(x=time, y = P1))+ geom_line(aes(x=time, y =P2))+geom_line(aes(x = time, y = P3))

#### data viz ####
Building_data <- Building_output
Community_data <- Community_output
head(Community_data)

#Prop_across_time <- Community_data 


category <- c("Begining", "Peak", "End")
Community_updated<-Community_data %>% 
  mutate(POP = S+I+R , Prop_S = S/POP, Prop_I = I/POP, Prop_R = R/POP) 
peak_time <- Community_updated%>%  filter( Prop_I ==max(Prop_I)) %>% select(time)
   
Community_bar_plot<-Community_updated %>%  filter(time == 1 | time == as.numeric(peak_time) | time == 100) %>%
  mutate(Category = category) %>% 
  select(Prop_S,Prop_I,Prop_R, Category) %>% 
  pivot_longer(cols = !Category)
Community_bar_plot
Community_bar_plot$Category <- factor(Community_bar_plot$Category, levels = c("Begining", "Peak", "End"))
Community_bar_plot$name <- factor(Community_bar_plot$name, levels = c("Prop_S", "Prop_I", "Prop_R"))

ggplot(Community_bar_plot,aes(x =Category, y = value, fill = name))+
  geom_col(position = "dodge")+
  scale_fill_manual(name="Infection status",
                    values=c("#2a52be","#C21807","#7bb661"),
                    breaks = c("Prop_S","Prop_I","Prop_R"),
                    labels = c("Susceptible","Infectious","Recovered"))+
  theme_minimal()+
  labs(x = NULL,
       y = "Proportion", 
       title = "Proportion of individuals that are Susceptible, Infected, or Recovered\nthroughout an outbreak")


clean_data <- Building_data %>% pivot_longer(cols = !time,
                              names_to = c("State", "Room"),
                              names_pattern = "([A-Za-z]+)(\\d+)",
                              values_to = "Number")
clean_data %>% filter(State == "S") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of Susceptible individuals")+ggtitle("Susceptible individuals across rooms")

clean_data %>% filter(State == "I") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of Infectious individuals")+ggtitle("Infectious individuals across rooms")
clean_data %>% filter(State == "R") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of Recovered individuals")+ggtitle("Recovered individuals across rooms")
clean_data %>% filter(State == "P") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of infectious particles")+ggtitle("Infectious particles in rooms within a building")

