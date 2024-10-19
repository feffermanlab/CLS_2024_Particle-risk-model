

## Step by step for creating risk model code

#For this document I hope to go through my code for the risk model and double check that each piece does what it is supposed to do. The main document has functions first and then parameter declaration and then running the simulations. For this document I will go in order of the pieces needed, so I will start withe the community level model, then to adjacency matrix, transition, and then lastly all the code need to solve the system of equations


library(deSolve)
library(tidyverse)
library(igraph)
library(ggnetwork)
library(GGally)
library(intergraph)
library(ggthemes)
library(paletteer)
library(pals)
library(ggrepel)
library(ggpubr)



##Community model parameters

community_pop <- 100000 #pop size of the larger community
init_conds <- c(S = community_pop-1, I = 1, R = 0) 
print(init_conds)
parms <- data.frame(bet =0.0000035 , gam = 1/21) 
print(parms)
times <- seq(from = 1, to = 144, by = 1)
print(times)



#r communtiy model function
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


#Community model output**
  


Community_output <- data.frame(lsoda(y = init_conds, func = SIR_community_model,times = times, parms=parms))

Community_output %>% pivot_longer(cols = !time) %>% arrange(desc(time))%>% 
  ggplot(aes(x=time,y =value, color = name))+geom_line()+
  scale_color_manual(name=NULL,values=c("blue","red","purple"),breaks = c("S","I","R"))+
  theme_classic()+
  labs(x= "Time (days)", y = "Number of individuals")+ggtitle("Community outbreak of pathogen")



#Make building**
  


small_bld_3_rooms <- matrix(c(c(0,1,1),
                              c(1,0,1),
                              c(1,1,0)),nrow = 3,ncol = 3)
small_bld_5_rooms <- matrix(c(c(0,1,1,1,1),
                              c(1,0,1,1,0),
                              c(1,1,0,1,0),
                              c(1,1,1,0,0),
                              c(1,0,0,0,0)),nrow=5, ncol = 5)

church_adjacency_matrix <- matrix(c(c(0,rep(1,3),rep(0,8),rep(1,4),rep(0,30-16)), #column 1 - main area / Hallway
                                    c(1,0,0,1,rep(0,12),rep(1,4),rep(0,30-20)), # column 2 - Hallway
                                    c(1,0,0,1,rep(0,4),rep(1,4),rep(0,30-12)), # column 3 - Hallway
                                    c(1,1,1,0,1,0,1,1,rep(0,30-8)), #column 4 Main room
                                    c(rep(0,3),1,0,1,0,1,rep(0,17),1,rep(0, 30-26)), # # column 5 Hallway
                                    c(rep(0,4),1,0,1,rep(0,15),rep(1,5),0,1,1), # Column 6 Hallway
                                    c(rep(0,3),1,0,1,0,1,rep(0,30-8)), # Column 7 Hallway
                                    c(rep(0,3),1,1,0,1,rep(0,21-8),1,1,rep(0,30-22)), # Column 8 Hallway
                                    rep(c(0,0,1,rep(0,30-3)),4), # columns 9-12
                                    rep(c(1,rep(0,30-1)),4), # columns 13-16
                                    rep(c(0,1,rep(0,30-2)),4), # columns 17-20
                                    rep(c(rep(0,7),1,rep(0, 30-8)),2), # Columns 21 and 22
                                    rep(c(rep(0,5),1,rep(0,30-6)),3), # columns 23-25
                                    c(rep(0,4),1,1,rep(0,30-6)), # Column 26
                                    c(rep(0,5),1,rep(0,21),1,rep(0,30-28)), # column 27
                                    c(rep(0,26),1,rep(0,30-27)), # column 28
                                    c(rep(0,5),1,rep(0,30-6)), # column 29 
                                    c(rep(0,5),1,1,rep(0,30-7))# column 30
                                    
),nrow=30, ncol = 30)

small_bld_3_rooms_graph<- graph_from_adjacency_matrix(small_bld_3_rooms, mode = "undirected")
#take a look
plot(small_bld_3_rooms_graph)
church_graph<- graph_from_adjacency_matrix(church_adjacency_matrix, mode = "undirected")
#take a look
plot(church_graph)


#r set carrying capacities
graph_to_use <- small_bld_3_rooms
adjacency_matrix_to_use <- church_adjacency_matrix 
N_rooms <- ncol(adjacency_matrix_to_use) #number of rooms
# Church_C<- c(200, #column 1 - main area / Hallway
#        100, # column 2 - Hallway
#        100, # column 3 - Hallway
#        400, #column 4 Main room
#        200
# ) 
small_bld_3_rooms_C <- c(100,200,300)
Church_C <- c(20, #column 1 - main area / Hallway
              5, # column 2 - Hallway
              5, # column 3 - Hallway
              362, #column 4 Main room
              5, # # column 5 Hallway
              7, # Column 6 Hallway - slightly bigger hallway
              5, # Column 7 Hallway
              5, # Column 8 Hallway
              50, # Column 9 Large classroom
              30, # column 10 classroom
              30, #column 11 classroom
              2, # column 12 janitors closet
              30, # column 13 classroom
              5, # column 14 - bathroom
              30, # column 15 classroom
              5, # column 16 - bathroom
              30, # column 17 classroom
              2, # column 18 - closet
              30, # column 19 classroom
              50, # column 20 large classroom
              2, # column 21 dressing room
              2, # column 22 dressing room
              4, # column 23 office
              4, # column 24 office
              4, # column 25 office
              4, # column 26 office
              5, # column 27 office - large
              2, # column 28 small bathroom
              4, # column 29 office
              4 # column 30 office
)


#Now that we have building defined by the adjacency matrix, we need to specify the conditions in the building
#r Building paramaters}
Max_Building_Capacity <- sum(Church_C) 
Prop_full <- 0.8 # how full do we want our building capacity to be
Adj_Max_Building_Capacity <- round(Prop_full*Max_Building_Capacity,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
delt <- Adj_Max_Building_Capacity/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)



#Next, we need to put people in our building. These will be our initial conditions for the building level model.


#r initialize building}

day <- 31 #What day to extract results from the community level model
N_rooms <- length(Church_C)

Bld_setup_func_v2 <- function(Community_output,day, delt,N_rooms, C_x, N_total){
  Building_ICs <- Community_output[day,] #Retrieves the number of S, I, and R individuals in the community at a particular day.
  
  N_x  <- c(rep(0, N_rooms)) #initialize a vector
  for (i in 1:length(C_x)) {
    if (sum(N_x) >= N_total) break  # Stop if all individuals are distributed
    
    # Maximum we can place in this room, constrained by remaining individuals and room capacity
    max_for_room <- min(C_x[i], N_total - sum(N_x))
    
    # Randomly assign a number between 0 and max_for_room
    N_x[i] <- sample(0:max_for_room, 1)
  }
  
  # If after the first pass, there are still individuals left, distribute the remaining in a while loop
  while (sum(N_x) < N_total) {
    # Randomly select a room index
    room_index <- sample(1:length(N_x), 1)
    
    # Check if the selected room can accommodate more individuals
    if (N_x[room_index] < C_x[room_index] && sum(N_x) < N_total) {
      N_x[room_index] <- N_x[room_index] + 1
    }
    
    # Break if we've distributed all individuals or hit the building limit
    if (sum(N_x) >= N_total || sum(N_x) >= N_total) {
      break
    }
  }
  
  #Now we need to find how many S, I, and R individuals we should have in the building based on the proportionality constant solved for above
  Sb <- Community_output$S[day]*delt # number of susceptible in building
  Ib <- Community_output$I[day]*delt # number of Infectious in building
  Rb <- Community_output$R[day]*delt # number of Recovered in building
  
  # calculate the proportion of S,I, and R thatshould be in the building
  Sb_prop <- Sb/N_total
  Ib_prop <- Ib/N_total
  Rb_prop <- Rb/N_total
  
  # make the number of S, I, and R in each room is proportional to that in the building 
  S_x <- N_x*Sb_prop
  I_x <- N_x*Ib_prop
  R_x <- N_x*Rb_prop
  
  S_x_prop <- S_x/N_total
  I_x_prop <- I_x/N_total
  R_x_prop <- R_x/N_total
  
  # we start with no particles in the building
  P_x <- c(rep(0,N_rooms))
  
  return(data.frame(S=S_x_prop,I=I_x_prop,R=R_x_prop, P=P_x, N_x = N_x))
}

Church_setup <- Bld_setup_func_v2(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms,C_x = Church_C,N_total = Adj_Max_Building_Capacity)
Church_setup
r1_pop <-(Church_setup$S[1]+Church_setup$I[1]+Church_setup$R[1])*Adj_Max_Building_Capacity
r2_pop <-(Church_setup$S[2]+Church_setup$I[2]+Church_setup$R[2])*Adj_Max_Building_Capacity
r3_pop<- (Church_setup$S[3]+Church_setup$I[3]+Church_setup$R[3])*Adj_Max_Building_Capacity
r1_pop+r2_pop+r3_pop



#Checking initial conditions**
#   Now lets check that our building setup function did what I wanted it to. 
# Each room's proportion should sum to 1. And if we convert back to numbers using N_x the numbers should equal

#Church_setup %>% mutate(total_prop = S_prop+I_prop+R_prop)


#r props in building match community}
#prop of susceptible in building
sum(Church_setup$S)
#prop of susceptible in community
Community_output$S[day]/(Community_output$S[day]+Community_output$I[day]+Community_output$R[day])

#prop of infectious in building
sum(Church_setup$I)
#prop of infectious in community
Community_output$I[day]/(Community_output$S[day]+Community_output$I[day]+Community_output$R[day])

#prop of recovered in building
sum(Church_setup$R)
#prop of susceptible in community
Community_output$R[day]/(Community_output$S[day]+Community_output$I[day]+Community_output$R[day])


#Now people and particles need to move (Transition matrices)**


Create_T_Matrix <-function(adjacency_matrix_to_use){
  set.seed(123145) # <- easier for debugging
  T_mov <- data.frame(matrix(runif(N_rooms^2), nrow = N_rooms)) #populates a square matrix/dataframe with random numbers between 0 and 1 for the number of rooms that our building has.
  diag(T_mov) <- 0 #set diagonal to 0
  T_mov <- adjacency_matrix_to_use*T_mov #restrict the movement according to our network/adjacency matrix
  T_mov_norm <- t(apply(T_mov, 1, function(x) x / sum(x))) # normalize so that there aren't more people moving than what can (rows should sum to 1)
T_mov <-T_mov_norm
return(T_mov)
}

Church_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use)
sum(Church_T_mov[1,])
sum(Church_T_mov[2,])
sum(Church_T_mov[3,])
# sum(Church_T_mov[4,])
# sum(Church_T_mov[5,])
r1_in <- sum(Church_T_mov[,1])
r1_out <- sum(Church_T_mov[1,])
r1_change <- r1_in-r1_out
r2_in <- sum(Church_T_mov[,2])
r2_out <- sum(Church_T_mov[2,])
r2_change <- r2_in-r2_out
r3_in <- sum(Church_T_mov[,3])
r3_out <- sum(Church_T_mov[3,])
r3_change <- r3_in - r3_out
# We can use the same function for the particle matrix
Church_theta_mov <- Create_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use)

#Now set parameters for the model**
  
  
#r parameters}
parms <-data.frame(s=100,a=5, d=3,lam = 1)
# s = shedding, a = absorption, d = decay, lam = scalar for room capacities
#N_b <- Adj_Max_Building_Capacity #Building population size
Maxtime <- 24*3
times <- seq(from = 0, to = Maxtime, by = 0.2)
m <- 4 # number of equations per room


#For the model we are going to have to sum over some of the vectors/matrices so we will put those steps in functions. This is needed for the movement of individuals and particles. Particles and people will have to have a different function because there re capacities on rooms for people but not particles.

#r movement functions}

flux_in_people <- function(N_rooms, Transition_matrix,State,Room_pops,Carrying_capacity,t){
  # Transition_matrix <- matrix(c(0,0.5,1,0.7,0,0,0.3,0.5,0), nrow = 3,ncol = 3)
  # State <- Church_setup$S_prop
  # Room_pops <- Church_setup$N_x
  # Carrying_capacity <- Church_C
  all_room_change <- c(seq(N_rooms))
  for(x in 1:N_rooms){ #flow in to room x from other rooms (j)
    flux_in_temp <- 0
    for(j in 1:N_rooms){
      flux_in_temp <- flux_in_temp + State[j]*Transition_matrix[j,x]*(1-(Room_pops[x]/Carrying_capacity[x]))
    }
    all_room_change[x] <- flux_in_temp
  }
  
  test <- c(M = as.vector(all_room_change))
  return(test)
}

flux_out_people <- function(N_rooms, Transition_matrix,State,Room_pops,Carrying_capacity){
  all_room_change <- c(seq(N_rooms))
  for(x in 1:N_rooms){
    flux_out_temp <- 0
    for(j in 1:N_rooms){
      flux_out_temp <- flux_out_temp + State[x]*Transition_matrix[x,j]*(1-(Room_pops[j]/Carrying_capacity[j]))
    }
    all_room_change[x] <- flux_out_temp
  }
  test <- c(M = as.vector(all_room_change))
  return(test)
}

flux_in_particles <- function(N_rooms, Transition_matrix,State){
  all_room_change <- c(seq(N_rooms))
  for(x in 1:N_rooms){
    flux_in_temp <- 0
    for(j in 1:N_rooms){
      flux_in_temp <- flux_in_temp + State[j]*Transition_matrix[j,x]
    }
    all_room_change[x] <- flux_in_temp
  }
  
  test <- c(M = as.vector(all_room_change))
  return(test)
}

flux_out_particles <- function(N_rooms, Transition_matrix,State){
  all_room_change <- c(seq(N_rooms))
  for(x in 1:N_rooms){
    flux_out_temp <- 0
    for(i in 1:N_rooms){
      flux_out_temp <- flux_out_temp + State[x]*Transition_matrix[x,i]
    }
    all_room_change[x] <- flux_out_temp
  }
  test <- c(M = as.vector(all_room_change))
  return(test)
}


#Function for the whole model**

Particle_model_v2 <- function(t, x, parms,T_mov, theta_mov, adjacency_matrix_to_use,C_x,N_b){
  
  # x <- Church_Init_conds
  # T_mov <- Church_T_mov
  # 
  # theta_mov <- Church_theta_mov
  # adjacency_matrix_to_use <- small_bld_3_rooms
  # C_x <- small_bld_3_rooms_C
  # N_b <- N_b
  ncompartment <- 5
  n_rooms <- length(x)/ncompartment
  S <- as.matrix(x[1:n_rooms])
  I <- as.matrix(x[(n_rooms+1):(2*n_rooms)])
  R <- as.matrix(x[(2*n_rooms+1):(3*n_rooms)])
  P <- as.matrix(x[(3*n_rooms+1):(4*n_rooms)])
  N_x <- as.matrix(x[(4*n_rooms+1):(5*n_rooms)])
  
  
  with(parms,{
    
    dS <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=S,Room_pops = (S+I+R)*N_x,Carrying_capacity = C_x)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=S,Room_pops = (S+I+R)*N_x,Carrying_capacity = C_x))
    # s1_in <- S[1]*T_mov[1,1]*(1-(((S[1]+I[1]+R[1])*N_b)/C_x[1]))
    # s2_in <- S[2]*T_mov[2,1]*(1-(((S[1]+I[1]+R[1])*N_b)/C_x[1]))
    # s3_in <- S[3]*T_mov[3,1]*(1-(((S[1]+I[1]+R[1])*N_b)/C_x[1]))
    # 
    # s1_out <- S[1]*T_mov[1,1]*(1-(((S[1]+I[1]+R[1])*N_b)/C_x[1]))
    # s2_out <- S[1]*T_mov[1,2]*(1-(((S[2]+I[2]+R[2])*N_b)/C_x[2]))
    # s3_out <- S[1]*T_mov[1,3]*(1-(((S[3]+I[3]+R[3])*N_b)/C_x[3]))
    # 
    # Total_Sx <- s1_in+s2_in+s3_in - (s1_out+s2_out+s3_out)
    dI <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=I,Room_pops = (S+I+R)*N_x,Carrying_capacity = C_x)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=I,Room_pops = (S+I+R)*N_x,Carrying_capacity = C_x))
    
    dR <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=R,Room_pops = (S+I+R)*N_x,Carrying_capacity = C_x)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=R,Room_pops = (S+I+R)*N_x,Carrying_capacity = C_x))
    #last step will be the particle EQ
    
    dP <- s*as.matrix(I)*as.matrix((S+I+R)*N_x) - as.matrix(a*P/(lam*C_x*((S+I+R)*N_x)))-d*as.matrix(P) + as.matrix(as.matrix(flux_in_particles(N_rooms, theta_mov,State = P)) - as.matrix(flux_out_particles(N_rooms, theta_mov,State = P)))
    
    dN_x <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=N_x,Room_pops = N_x,Carrying_capacity = C_x)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=N_x,Room_pops = N_x,Carrying_capacity = C_x))
    # 
    dt <- c(dS,dI,dR,dP,dN_x)
    
    return(list(dt))})
  
}
```



```{r}
Church_Init_conds_v1 <-c(S=Church_setup$S, I = Church_setup$I, R = Church_setup$R, P = Church_setup$P)

Church_Init_conds_v2 <-c(S=Church_setup$S, I = Church_setup$I, R = Church_setup$R, P = Church_setup$P,N_x =Church_setup$N_x)

Church_output_v1 <- data.frame(lsoda(y = Church_Init_conds_v1, func = Particle_model_v1,times = times,
                                     parms = parms,
                                     adjacency_matrix_to_use=church_adjacency_matrix,
                                     theta_mov =Church_theta_mov,
                                     T_mov = Church_T_mov, 
                                     C_x=Church_C,N_b = N_b))
Church_output_v2 <- data.frame(lsoda(y = Church_Init_conds_v2, func = Particle_model_v2,times = times,
                                     parms = parms,
                                     adjacency_matrix_to_use=church_adjacency_matrix,
                                     theta_mov =Church_theta_mov,
                                     T_mov = Church_T_mov, 
                                     C_x=Church_C,N_b = Adj_Max_Building_Capacity))

```


```{r}
Church_data_clean <- Church_output_v1%>% pivot_longer(cols = !time,
                                                      names_to = c("State", "Room"),
                                                      names_pattern = "([A-Za-z]+)(\\d+)",
                                                      values_to = "Number")

# Church_data_ratios <- Church_data_clean %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% group_by(time,Room)  %>% 
#   mutate(N_x = S+I+R) %>% group_by(time,Room) %>% 
#   mutate(K_x = ((parms$s)*I)/((parms$a)*N_x), prop_to_K = P/K_x)
# Church_P_x_K_x_plot<-ggplot(Church_data_ratios,aes(x =time, y =prop_to_K,group= Room,color=Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")
# Church_P_x_K_x_plot
# Church_P_x_plot <- ggplot(Church_data_ratios, aes(x=time, y = P, color = Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")
# Church_P_x_plot

p1_v1 <- Church_data_clean %>% filter(State == "S") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Proportion of Susceptible individuals")+ggtitle("Proportion of Susceptible individuals across rooms")

p2_v1 <- Church_data_clean %>% filter(State == "I") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Proportion of Infectious individuals")+ggtitle("Proportion of Infectious individuals across rooms")
p3_v1 <- Church_data_clean %>% filter(State == "R") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Proportion of Recovered individuals")+ggtitle("Proportion of Recovered individuals across rooms")
p4_v1 <- Church_data_clean %>% filter(State == "P") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of infectious particles")+ggtitle("Infectious particles in rooms within a building")
p5_v1 <- Church_data_clean %>% filter(State == "x") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of individuals")+ggtitle("Number of individuals in rooms within a building")

# Church_data_clean %>% filter(State != "P") %>% 
#   group_by(time,Room)%>% summarise(total_room_pop=sum(Number)) %>% group_by(time) %>% summarise(tot_bld_pop =sum(total_room_pop)) %>% 
#   ggplot( aes(x=time, y = tot_bld_pop))+
#   geom_line()+theme_classic()+
#   labs(x = "Time (hours)", y= "Number of people")+ggtitle("number of people in building")

p6_v1 <- Church_data_clean %>% filter(State == "S"| State == "I" | State == "R") %>% ungroup() %>% group_by(time, Room) %>% summarise(Total_prop=sum(Number)) %>%
  ggplot(aes(x=time, y=Total_prop, color=Room))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in rooms within a building")

p7_v1 <- Church_data_clean %>% filter(State == "S"| State == "I" | State == "R") %>% ungroup() %>% group_by(time) %>% summarise(Total_prop=sum(Number)) %>%
  ggplot(aes(x=time, y=Total_prop))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in the building")


p8_v1 <- Church_data_clean %>% filter(State == "S"| State == "I" | State == "R") %>% ungroup() %>% group_by(time, Room) %>% summarise(Total_prop=sum(Number)*Adj_Max_Building_Capacity) %>%
  ggplot(aes(x=time, y=Total_prop, color=Room))+geom_line()+
  labs(x = "Time (hours)", y= "Number of individuals")+ggtitle("Number of individuals in the building")

```

```{r}
Church_data_clean <- Church_output_v2%>% pivot_longer(cols = !time,
                                                      names_to = c("State", "Room"),
                                                      names_pattern = "([A-Za-z]+)(\\d+)",
                                                      values_to = "Number")

# Church_data_ratios <- Church_data_clean %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% group_by(time,Room)  %>% 
#   mutate(N_x = S+I+R) %>% group_by(time,Room) %>% 
#   mutate(K_x = ((parms$s)*I)/((parms$a)*N_x), prop_to_K = P/K_x)
# Church_P_x_K_x_plot<-ggplot(Church_data_ratios,aes(x =time, y =prop_to_K,group= Room,color=Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")
# Church_P_x_K_x_plot
# Church_P_x_plot <- ggplot(Church_data_ratios, aes(x=time, y = P, color = Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")
# Church_P_x_plot

p1_v2 <- Church_data_clean %>% filter(State == "S") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Proportion of Susceptible individuals")+ggtitle("Proportion of Susceptible individuals across rooms")

p2_v2 <- Church_data_clean %>% filter(State == "I") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Proportion of Infectious individuals")+ggtitle("Proportion of Infectious individuals across rooms")
p3_v2 <- Church_data_clean %>% filter(State == "R") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Proportion of Recovered individuals")+ggtitle("Proportion of Recovered individuals across rooms")
p4_v2 <- Church_data_clean %>% filter(State == "P") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of infectious particles")+ggtitle("Infectious particles in rooms within a building")
p5_v2 <- Church_data_clean %>% filter(State == "x") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of individuals")+ggtitle("Number of individuals in rooms within a building")

# Church_data_clean %>% filter(State != "P") %>% 
#   group_by(time,Room)%>% summarise(total_room_pop=sum(Number)) %>% group_by(time) %>% summarise(tot_bld_pop =sum(total_room_pop)) %>% 
#   ggplot( aes(x=time, y = tot_bld_pop))+
#   geom_line()+theme_classic()+
#   labs(x = "Time (hours)", y= "Number of people")+ggtitle("number of people in building")

p6_v2 <- Church_data_clean %>% filter(State == "S"| State == "I" | State == "R") %>% ungroup() %>% group_by(time, Room) %>% summarise(Total_prop=sum(Number)) %>%
  ggplot(aes(x=time, y=Total_prop, color=Room))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in rooms within a building")

p7_v2 <- Church_data_clean %>% filter(State == "S"| State == "I" | State == "R") %>% ungroup() %>% group_by(time) %>% summarise(Total_prop=sum(Number)) %>%
  ggplot(aes(x=time, y=Total_prop))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in the building")


p8_v2 <- Church_data_clean %>% filter(State == "x") %>% ungroup() %>% group_by(time, Room) %>% summarise(Total_prop=sum(Number))%>% arrange(desc(Total_prop))%>%
  ggplot(aes(x=time, y=Total_prop, color=Room))+geom_line()+
  labs(x = "Time (hours)", y= "Number of individuals")+ggtitle("Number of individuals in the building")


```


```{r}
ggarrange(p1_v1,p1_v2,ncol = 2)
ggarrange(p2_v1,p2_v2,ncol = 2)
ggarrange(p3_v1,p3_v2,ncol = 2)
ggarrange(p4_v1,p4_v2,ncol = 2)
ggarrange(p5_v1,p5_v2,ncol = 2)
ggarrange(p6_v1,p6_v2,ncol = 2)
ggarrange(p7_v1,p7_v2,ncol = 2)
ggarrange(p8_v1,p8_v2,ncol = 2)






