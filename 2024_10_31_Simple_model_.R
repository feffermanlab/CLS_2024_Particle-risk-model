
## Step by step for creating risk model code


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


# delt is the proportionality constant
# C_x is the room capacity for each room
# N_total is the total number of individuals that we want in the building

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
  
  # calculate the proportion of S,I, and R that should be in the building
  Sb_prop <- Sb/N_total
  Ib_prop <- Ib/N_total
  Rb_prop <- Rb/N_total
  
  
  # we start with no particles in the building
  P_x <- c(rep(0,N_rooms))
  
  return(data.frame(Sb = Sb, Ib = Ib, Rb = Rb,Sb_prop = Sb_prop,Ib_prop = Ib_prop,Rb_prop= Ib_prop ,P=P_x, N_x = N_x))
}


Create_T_Matrix <-function(adjacency_matrix_to_use, N_rooms){
  set.seed(123145) # <- easier for debugging
  T_mov <- data.frame(matrix(runif(N_rooms^2), nrow = N_rooms)) #populates a square matrix/dataframe with random numbers between 0 and 1 for the number of rooms that our building has.
  
  T_mov <- adjacency_matrix_to_use*T_mov #restrict the movement according to our network/adjacency matrix
  T_mov_norm <- t(apply(T_mov, 1, function(x) x / sum(x))) # normalize so that there aren't more people moving than what can (rows should sum to 1)
  T_mov <-T_mov_norm
  return(T_mov)
}

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

#


#**Function for the whole model**

Particle_model_v3 <- function(t, x, parms,T_mov, theta_mov, adjacency_matrix_to_use,C_x,N_b,N_rooms,Ib_prop){
  # x <- Church_Init_conds
  # T_mov <- Church_T_mov
  # theta_mov <- Church_theta_mov
  # adjacency_matrix_to_use <- small_bld_3_rooms
  # C_x <- small_bld_3_rooms_C
  # N_b <- N_b
  N_rooms <- N_rooms
  ncompartment <- 2
  n_rooms <- length(x)/ncompartment
  N_x <- as.matrix(x[1:n_rooms])
  P <- as.matrix(x[(n_rooms+1):(2*n_rooms)])
  
  with(parms,{
    
    dN_x <- as.matrix((flux_in_people(N_rooms=N_rooms, Transition_matrix =T_mov, State=N_x,Room_pops = N_x,Carrying_capacity = C_x)) - flux_out_people(N_rooms=N_rooms, Transition_matrix = T_mov, State=N_x,Room_pops = N_x,Carrying_capacity = C_x))
    
    dP <- s*Ib_prop*as.matrix(N_x) - as.matrix(a*as.matrix(N_x)*as.matrix(P)*(1/(C_x*bet_bar))*(as.matrix(P)/bet_hat))-d*as.matrix(P) + as.matrix(as.matrix(flux_in_particles(N_rooms=N_rooms, theta_mov,State = P)) - as.matrix(flux_out_particles(N_rooms=N_rooms, theta_mov,State = P)))
    
    
    dt <- c(dN_x,dP)
    
    return(list(dt))})
}



#**Community model parameters**
#  #{r Community model parameters}
community_pop <- 100000 #pop size of the larger community
init_conds <- c(S = community_pop-1, I = 1, R = 0) 
print(init_conds)
parms <- data.frame(bet =0.0000035 , gam = 1/21) 
print(parms)
times <- seq(from = 1, to = 144, by = 1)
print(times)


##
#**Community modelfunction/equations**
 # #{r communtiy model function}
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


#**Community model output**
  
#  

Community_output <- data.frame(lsoda(y = init_conds, func = SIR_community_model,times = times, parms=parms))

Community_output %>% pivot_longer(cols = !time) %>% arrange(desc(time))%>% 
  ggplot(aes(x=time,y =value, color = name))+geom_line()+
  scale_color_manual(name=NULL,values=c("blue","red","purple"),breaks = c("S","I","R"))+
  theme_classic()+
  labs(x= "Time (days)", y = "Number of individuals")+ggtitle("Community outbreak of pathogen")




  
#### Adjacency matrices ####

small_bld_3_rooms <- matrix(c(c(1,1,1),
                              c(1,1,1),
                              c(1,1,1)),nrow = 3,ncol = 3)

small_bld_5_rooms <- matrix(c(c(0,1,1,1,1),
                              c(1,0,1,1,0),
                              c(1,1,0,1,0),
                              c(1,1,1,0,0),
                              c(1,0,0,0,0)),nrow=5, ncol = 5)
diag(small_bld_5_rooms) <-1

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


Office_adjacency_matrix <- matrix(c(c(0,1,rep(0,4),1,0,1,rep(0,34-9)), # Hallway 1
                                    c(1,0,1,rep(0,31)), # Hallway 2
                                    c(0,1,0,1,rep(0,5),1,1,1,rep(0,34-12)), # Hallway/room 3
                                    c(0,0,1,rep(0,6),rep(1,9),0,0,rep(1,3),rep(0,34-23)), # Hallway 4
                                    c(rep(0,3),1,0,1,1,rep(0,13),rep(1,5),rep(0,34-25)), # Hallway 5
                                    c(rep(0,4),1,0,1,1, rep(0,17),rep(1,4), rep(0,34-29)), # Hallway 6
                                    c(1,rep(0,3),1,1,0,rep(0,34-7)), # Hallway/room 7
                                    c(rep(0,5),1,rep(0,22),rep(1,3),rep(0,34-31)), # Hallway 8
                                    c(1,rep(0,30),1,1,1), # Hallway 9
                                    rep(c(0,0,1,1,rep(0,34-4)),3), # Rooms 10-12
                                    rep(c(0,0,0,1,rep(0,34-4)),5), # Rooms 13-17
                                    c(0,0,0,1,rep(0,14),1,1,rep(0,34-20)), # Room 18
                                    rep(c(rep(0,17),1,rep(0,34-18)),2), # Rooms 19 and 20
                                    rep(c(rep(0,3),1,1,rep(0,34-5)),3), # Rooms 21-23
                                    rep(c(rep(0,4),1,rep(0,34-5)),2), # Rooms 24 and 25
                                    rep(c(rep(0,5),1,rep(0,34-6)),3), # Rooms 26-28
                                    c(rep(0,5),1,0,1,rep(0,34-8)), # Room 29
                                    rep(c(rep(0,7),1,rep(0,34-8)),2), # Rooms 30 and 31
                                    rep(c(rep(0,8),1,rep(0,34-9)),3) # Rooms 32-34
                                    
),
nrow=34,ncol = 34)

Movie_adjacency_matrix <- matrix(c(c(0,1,1,rep(0,5),1,rep(0,5),1,0,1,rep(0,34-17)), # Room 1 hallway
                                   c(1,0,1,1,1,rep(0,6),rep(1,3),rep(0,3),rep(1,5),0,rep(1,4),rep(0,34-27)), # Room 2 - Hallway
                                   c(1,rep(0,13),1,1,1,rep(0,34-17)), # Room 3 - Hallway
                                   c(0,1,0,0,1,0,1,rep(0,18),rep(1,3),rep(0,34-28)), # Room 4 - Hallway
                                   c(0,1,0,1,0,1,rep(0,19),1,rep(0,4),1,rep(0,34-31)), # Room 5 - Hallway
                                   c(rep(0,4),1,rep(0,25),1,rep(0,34-31)), # Room 6 - Hallway
                                   c(rep(0,3),1,rep(0,3),1,rep(0,15),1,rep(0,3),1,rep(0,34-28)), # Room 7 - Hallway
                                   c(rep(0,6),1,rep(0,16),1,rep(0,5),1,rep(0,34-30)), # Room 8 - Hallway
                                   c(1, rep(0,8),1,rep(0,3),1,rep(0, 34-14)), # Room 9 main room
                                   c(rep(0,8),1,0,1,rep(0, 34-11)), # Room 10 Admin
                                   c(rep(0,9),1,rep(0, 34-10)), # Room 11 - admin
                                   rep(c(0,1, rep(0,34-2)),2), # Room 12 and 13 - facilities management
                                   c(0,1,rep(0,6),1,rep(0,34-9)), # Room 14
                                   c(1,rep(0, 34-1)), # Room 15 - Theater
                                   c(0,0,1,rep(0,34-3)), # Room 16 - Storage room
                                   c(1,0,1,rep(0,34-3)), # Room 17
                                   rep(c(0,1,rep(0,34-2)),4), # Rooms 18-21 restrooms and theater (19)
                                   c(0,1,rep(0,20),1,rep(0,34-23)), # Room 22 - cleaning
                                   c(rep(0,21),1, rep(0,34-22)), # Room 23 - cleaning
                                   c(0,1,rep(0,4),1,1,rep(0,20),1,rep(0,34-29)), # Room 24 - Theater
                                   c(0,1,rep(0,34-2)), # Room 25 - facilities
                                   c(0,1,0,1,1,rep(0,34-5)), # Room 26 - Theater
                                   c(0,1,0,1,rep(0,34-4)), # Room 27 - Theater
                                   c(rep(0,3),1,0,0,1,rep(0,34-7)), # Room 28 - Theater
                                   c(rep(0,23),1,rep(0,34-24)), # Room 29 Facilities
                                   c(rep(0,7),1,rep(0,34-8)), # Room 30 - storage
                                   c(rep(0,4),1,rep(0,34-5)), # Room 31 - Theater
                                   c(rep(0,5),1,rep(0,26),1,rep(0,34-33)), # Room 32 - maintenance
                                   c(rep(0,31),1,0,1), # Room 33 - Employee room
                                   c(rep(0,32),1,0) # Room 34 - employee bathroom
), 

nrow = 34,ncol = 34)

Movie_graph<- graph_from_adjacency_matrix(Movie_adjacency_matrix, mode = "undirected")

University_adjacency_matrix<- matrix(c(rep(c(rep(0,34),1,rep(0,231-35)),6),#Basement floor -rooms 1-6
                                       rep(c(rep(0,35),1,rep(0,231-36)),4), #rooms 7-10
                                       c(rep(0,35),1,0,0,1,1,rep(0,231-40)), #room 11
                                       c(rep(0,36),1,rep(0,231-37)), #room 12
                                       c(rep(0,13),1,rep(0,22),1,rep(0,231-37)), #room 13
                                       c(rep(0,12),1,0,1,0,1,rep(0,20),1, rep(0,231-38)), #room 14
                                       c(rep(0,13),1,0,1,0,1,1,rep(0,18),1,rep(0,231-38)), #room 15
                                       c(rep(0,14),1,rep(0,4),1,1, rep(0,16),1,rep(0,231-38)), #room 16
                                       c(rep(0,13),1,rep(0,231-14)), #room 17
                                       rep(c(rep(0,14),1,rep(0,231-15)),2), #rooms 18-19
                                       rep(c(rep(0,15),1,rep(0,231-16)),2), #rooms 20-21
                                       c(rep(0,39),1,rep(0,231-40)), #room 22
                                       c(rep(0,39),1,1,rep(0, 231-41)), #room 23
                                       rep(c(rep(0,40), 1, rep(0,231-41)),2), #rooms 24-25
                                       c(rep(0,26),1,rep(0,13),1,1,rep(0,231-42)), #room 26
                                       c(rep(0,25),1,rep(0,16),1,rep(0,231-43)), #room 27
                                       c(rep(0,29),1,rep(0,11),1,rep(0,231-42)), # room 28
                                       c(rep(0,30),1,rep(0,10),1,rep(0,231-42)), #room 29
                                       c(rep(0,27),1,rep(0,231-28)), #room 30  
                                       c(rep(0,28),1,rep(0,231-29)), #room 31
                                       c(rep(0,41),1,rep(0,231-42)), #room 32
                                       c(rep(0,41),1,1,rep(0,231-43)), #room 33
                                       c(rep(0,42),1,rep(0,231-43)), #room 34
                                       #BASEMENT LEVEL HALLWAYS
                                       c(rep(1,6),rep(0,29),1,rep(0,7),1,rep(0,231-44)), #room 35 - (hallway)
                                       c(rep(0,6),rep(1,5),rep(0,25),1,0,1,rep(0,231-39)), #room 36 - (hallway)
                                       c(rep(0,11),1,1,rep(0,22),1,0,1,rep(0,6),1,rep(0,231-45)), #room 37 - (hallway)
                                       c(rep(0,13),rep(1,3),rep(0,20),1,0,0,1,rep(0,231-40)), #room 38 - (hallway)
                                       c(rep(0,10),1,rep(0,24),1,0,0,0,1,rep(0,231-40)), #room 39 - (hallway)
                                       c(rep(0,10),1,rep(0,10),1,1,rep(0,14),1,1,0,1,0,0,0,0,1,rep(0,231-46)), #room 40 - hallway
                                       c(rep(0,22),rep(1,4),rep(0,13),1,0,1,rep(0,231-42)), #room 41 - (hallway)
                                       c(rep(0,25),1,0,1,1,0,0,1,1,rep(0,9),1, rep(0,231-43)), #room 42 - (hallway)
                                       c(rep(0,32),1,1,rep(0,7),1,rep(0,4),1,rep(0,231-47)), #room 43 - (hallway)
                                       #STAIRS AND ELEVATORS * will fill in as I go, want to check the network of each floor
                                       c(rep(0,34),1,rep(0,188),1,rep(0,231-224)), #room 44 - (stairwell 1)
                                       c(rep(0,36),1,rep(0,14),1,rep(0,90),1,0,1,rep(0,64),1,0,1,rep(0,231-212)), #room 45 - (elevator 1)
                                       c(rep(0,39),1,rep(0,46),1,rep(0,62),1,rep(0,66),1,rep(0,231-217)), #room 46 - (elevator 2)
                                       c(rep(0,42),1,rep(0,183),1,rep(0,231-227)),  #room 47 - (stairwell 2)
                                       ################# LEVEL 1 ##
                                       rep(c(rep(0,80),1,rep(0,231-81)),2), # rooms 48 and 49
                                       c(rep(0,80),1,0,1,rep(0,231-83)), # room 50
                                       c(rep(0,80),1,rep(0,231-81)), #room 51
                                       c(rep(0,44),1,rep(0,7),1,rep(0,11),1,rep(0,15),1,1,rep(0,231-82)), # room 52
                                       c(rep(0,51),1,rep(0,231-52)), # room 53
                                       rep(c(rep(0,82),1,rep(0,231-83)),7), # rooms 54-60,
                                       c(rep(0,83),1,rep(0,231-84)), # room 61
                                       c(rep(0,66),1,rep(0,17),1,0,1,rep(0,231-87)), #room 62
                                       c(rep(0,84),1,rep(0,231-85)), # room 63
                                       c(rep(0,66),1,rep(0,17),1,rep(0,231-85)), # room 64
                                       c(rep(0,51),1,rep(0,32),1,rep(0,231-85)), # room 65
                                       c(rep(0,84),1,rep(0,231-85)), # room 66
                                       c(rep(0,61),1,0,1,rep(0,231-64)), #room 67
                                       rep(c(rep(0,86),1,rep(0,4),1,rep(0,231-92)),2), # room 68 and 69
                                       rep(c(rep(0,89),1,rep(0,231-90)),3), # room 70:72
                                       rep(c(rep(0,90),1,rep(0,231-91)),3), # room 73-75
                                       rep(c(rep(0,91),1,rep(0,231-92)),4), # room 76:79
                                       c(rep(0,85),1,0,1,0,0,1,rep(0,231-91)), # room 80
                                       c(rep(0,47),rep(1,5),rep(0,29),1,rep(0,141),1,rep(0,231-224)), # room 81 - (hallway)
                                       c(rep(0,51),1,rep(0,28),1,0,1,1,1,rep(0,231-85)), # room 82 - (Hallway)
                                       c(rep(0,49),1,0,0,0,rep(1,7),rep(0,21),1,0,1,rep(0,231-84)), # room 83 - (Hallway)
                                       c(rep(0,60),1,rep(0,20),1,1,0,1,1,1,1,1,rep(0,231-89)), # room 84 - (Hallway)
                                       c(rep(0,61),rep(1,5),rep(0,15),1,0,1,rep(0,231-84)), # room 85 - (Hallway)
                                       c(rep(0,79),1,0,0,0,1,0,0,0,1,rep(0,231-88)), # room 86 - (Hallway)
                                       c(rep(0,45),1,rep(0,15),1,rep(0,5),1,1,rep(0,14),1,rep(0,4),1,1,rep(0,231-90)), # room 87 - (Hallway)
                                       c(rep(0,79),1,rep(0,3),1,0,1,rep(0,4),1,rep(0,231-91)), # room 88 - (Hallway)
                                       c(rep(0,83),1,0,0,1,rep(0,3),1,rep(0,231-91)), # room 89 - (Hallway)
                                       c(rep(0,69),rep(1,3),rep(0,14),1, rep(0,231-87)), # room 90 - Hallway
                                       c(rep(0,72),1,1,1,rep(0,4),1, rep(0,7),1,1,0,0,1,1, rep(0,231-93)), # room 91 - (Hallway)
                                       c(rep(0,67),1,1,rep(0,6),rep(1,4),rep(0,147),1,rep(0,231-227)), # room 92 - (Hallway)
                                       ########################## STAIR ON LEVEL 1 room - 93###
                                       c(rep(0,90),1,rep(0,138),1,0), # room 93 ( Staircase)
                                       ########################## BEGIN FLOOR 2 ###
                                       rep(c(rep(0,140),1,rep(0,231-141)),3), # rooms 94:96 
                                       rep(c(rep(0,141),1,rep(0,231-142)),4), # rooms 97:100
                                       rep(c(rep(0,143),1,rep(0,231-144)),4), #rooms 101: 104
                                       c(rep(0,141),1,1,rep(0,231-143)), # room 105
                                       c(rep(0,142),1,rep(0,231-143)), #room 106
                                       c(rep(0,144),1,rep(0,231-145)), # room 107 
                                       c(rep(0,112),1,1,rep(0,30),1,1,1,rep(0,231-147)), # room 108
                                       c(rep(0,144),1,1,1,rep(0,231-147)), #room 109
                                       c(rep(0,144),1,1,rep(0,231-146)), # room 110
                                       c(rep(0,145),1,rep(0,231-146)), # room 111
                                       c(rep(0,146),1, rep(0,231-147)), # room 112
                                       rep(c(rep(0,107),1,rep(0,231-108)),2), # rooms 113 and 114
                                       rep(c(rep(0,147),1,rep(0,231-148)),2), # rooms 115 and 116
                                       c(rep(0,149),1,rep(0,231-150)), # room 117
                                       rep(c(rep(0,148),1,rep(0,231-149)),3), # room 118:120
                                       c(rep(0,121),1,rep(0,25),1,rep(0,6),1, rep(0,231-155)), # room 121
                                       c(rep(0,120),1,rep(0,32),1,rep(0,231-154)), # room 122
                                       rep(c(rep(0,154),1,rep(0,231-155)),3), # rooms 123:125
                                       rep(c(rep(0,152),1,rep(0,231-153)),2), # rooms 126 and 127 
                                       rep(c(rep(0,153),1,rep(0,231-154)),5), # rooms 128:132 
                                       rep(c(rep(0,152),1,1,rep(0,231-154)),2), # rooms 133 and 134
                                       c(rep(0,135),1,0,1,1,rep(0,11),1,rep(0,231-151)), # room 135
                                       c(rep(0,134),1,0,0,0,1,rep(0,231-139)), # room 136
                                       c(rep(0,137),1,1,rep(0,231-139)), # room 137
                                       c(rep(0,134),1,0,1,rep(0,231-137)), # room 138
                                       c(rep(0,134),1,1,1,rep(0,231-137)), # room 139
                                       c(rep(0,150),1,rep(0,231-151)), # room 140
                                       ############ start of level 2 Hallways ###
                                       c(rep(0,93),1,1,1,rep(0,45),1,rep(0,82),1,rep(0,231-225)), # room 141 - (Hallway)
                                       c(rep(0,96),rep(1,4),rep(0,4),1,rep(0,35),1,0,1,1,rep(0,231-144)), # room 142 - (Hallway)
                                       c(rep(0,44),1, rep(0,59),1,1,rep(0,35),1,0,0,1,rep(0,231-145)), # room 143 - (Hallway)
                                       c(rep(0,100),rep(1,4),rep(0,45),1,1,rep(0,231-151)), # room 144 - (Hallway)
                                       c(rep(0,44),1,rep(0,61),rep(1,4),rep(0,32),1,rep(0,231-143)), # room 145 - (Hallway)
                                       c(rep(0,107),rep(1,4),rep(0,231-111)), # room 146 - (Hallway)
                                       c(rep(0,107),1,1,0,0,1,rep(0,35),1,rep(0,231-148)), # room 147 - (Hallway)
                                       c(rep(0,114),1,1,rep(0,4),1,rep(0,25),1,0,0,1,rep(0,231-150)), # room 148
                                       c(rep(0,117), 1,1,1,rep(0,29),1,rep(0,231-150)), # room 149 - (Hallway)
                                       c(rep(0,45),1,rep(0,70),1,rep(0,26),1,rep(0,3),1,1,0,1,1,rep(0,77),1,0), # room 150 - (Hallway)
                                       c(rep(0,134),1,0,0,0,1,1,rep(0,3),1,rep(0,5),1,0,0,1,rep(0,231-153)), # room 151 - (Hallway)
                                       c(rep(0,149),1,0,0,1,rep(0,231-153)), # room 152 - (Hallway)
                                       c(rep(0,125),1,1,rep(0,5),1,1,rep(0,4),1,rep(0,11),1,1,rep(0,77),1,0),# room 153 - (Hallway)
                                       c(rep(0,121),1,rep(0,4),rep(1,8),rep(0,18),1,0,1,rep(0,72),1,rep(0,231-228)), # room 154 - (Hallway)
                                       c(rep(0,120),1,0,1,1,1,rep(0,28),1,rep(0,231-154)), # room 155 - (Hallway),
                                       ##################### BEGIN LEVEL 3 ### !!!!!!!!!!!!
                                       rep(c(rep(0,207),1,rep(0,231-208)),2), # room 156 and 157
                                       c(rep(0,207),1,1,rep(0,231-209)), # room 158
                                       rep(c(rep(0,208),1,rep(0,231-209)),4), # rooms 159-162
                                       rep(c(rep(0,210),1,rep(0,231-211)),4), # rooms 163-166
                                       c(rep(0,208),1,1,rep(0,231-210)), # room 167
                                       c(rep(0,211),1,rep(0,231-212)), # room 168
                                       c(rep(0,173),1,1,rep(0,36),1,1,1,1,rep(0,231-215)), # room 169
                                       c(rep(0,211),1,1,1,rep(0,231-214)), # room 170
                                       c(rep(0,211),1,1,rep(0,231-213)), # room 171
                                       c(rep(0,212),1,rep(0,231-213)), # room 172
                                       c(rep(0,213),1,rep(0,231-214)), # room 173
                                       rep(c(rep(0,168),1,rep(0,231-169)),2), # rooms 174 and 175
                                       rep(c(rep(0,214),1,rep(0,231-215)),2), # rooms 176 and 177
                                       c(rep(0,214),1,rep(0,6),1,1,rep(0,231-223)), # room 178
                                       c(rep(0,222),1,rep(0,231-223)), # room 179
                                       rep(c(rep(0,221),1,rep(0,231-222)),3), # rooms 180-182
                                       c(rep(0,215),1,1,rep(0,231-217)), # room 183
                                       rep(c(rep(0,215),1,rep(0,231-216)),3), # rooms 184-186
                                       c(rep(0,219),1,rep(0,231-220)), # room 187
                                       c(rep(0,219),1,1,rep(0,231-221)), # room 188
                                       rep(c(rep(0,220),1,rep(0,231-221)),4), # room 189-192
                                       c(rep(0,220),1,rep(0,231-221)), # room 193
                                       rep(c(rep(0,219),1,1,rep(0,231-221)),2), # rooms 194 and 195
                                       c(rep(0,217),1,rep(0,231-218)), # room 196
                                       c(rep(0,197),1,1,1,rep(0,17),1,rep(0,231-218)), # room 197
                                       c(rep(0,196),1,rep(0,231-197)), # room 198 
                                       c(rep(0,196),1,0,0,rep(1,5),rep(0,231-204)), # room 199
                                       c(rep(0,196),1,0,1,rep(0,231-199)), # room 200
                                       rep(c(rep(0,198),1,rep(0,231-199)),4), # rooms 201-204
                                       rep(c(rep(0, 217),1,rep(0,231-218)),3), # rooms 205-207
                                       ######## LEVEL 3 HALLWAYS ###
                                       c(rep(0,43),1,rep(0,111),1,1,1,rep(0,50),1,rep(0,231-209)),# room 208 - (Hallway)
                                       c(rep(0,157),rep(1,5),rep(0,4),1,rep(0,40),1,0,1,1,rep(0,231-211)), # room 209 - (Hallway)
                                       c(rep(0,44),1,rep(0,121),1, rep(0,41),1,0,0,1,rep(0,231-212)), # room 210 - (Hallway)
                                       c(rep(0,162),rep(1,4),rep(0,42),1,rep(0,7),1,1,rep(0,231-218)), # room 211 - (Hallway)
                                       c(rep(0,44),1,rep(0,122),rep(1,4),rep(0,38),1,rep(0,231-210)), # room 212 - (Hallway)
                                       c(rep(0,168),rep(1,4),rep(0,231-172)), # room 213 - (Hallway)
                                       c(rep(0,168),1,1,0,0,1,rep(0,41),1,rep(0,231-215)), # room 214 - (Hallway)
                                       c(rep(0,168),1, rep(0,6),1,1,1,rep(0,35),1,0,0,1,rep(0,231-217)), # room 215 - (Hallway)
                                       c(rep(0,182),rep(1,4),rep(0,30),1,rep(0,231-217)), # room 216 - (Hallway)
                                       c(rep(0,45),1,rep(0,136),1,rep(0,27),1,rep(0,3),1,1,0,0,1,rep(0,231-220),1), # room 217 - (Hallway)
                                       c(rep(0,195),1,1,rep(0,7),rep(1,3),rep(0,3),1,rep(0,8),1,rep(0,231-220)), # room 218 - (Hallway)
                                       c(rep(0,216),1,0,0,1,rep(0,231-220)), # room 219 - (Hallway)
                                       c(rep(0,92),1,rep(0,93),1,1,rep(0,5),1,1,0,0,0,1,rep(0,18),1,1,0,1,rep(0,231-221)), # room 220 - (Hallway)
                                       c(rep(0,46),1,rep(0,140),rep(1,8),rep(0,24),1,0,1,1,rep(0,231-223)), # room 221 - (Hallway)
                                       c(rep(0,177),1,0,1,1,1,rep(0,38),1,rep(0,231-221)), # room 222 - (Hallway)
                                       c(rep(0,177),1,1,rep(0,41),1,rep(0,231-221)), # room 223 - (Hallway)
                                       c(rep(0,43),1,rep(0,36),1,rep(0,143),1, rep(0,231-225)),# room 224 NW stairwell floor 1 to floor 2 (connects to 44 stair below, 81 - hallway on floor, 225 stair above)
                                       c(rep(0,140),1,rep(0,82),1,0,1, rep(0, 231-226)),# room 225 NW stairwell floor 2 to floor 3 (connects to 224 stair below, 141 - hallway on floor, 226 stair above)
                                       c(rep(0,207),1,rep(0, 16),1,rep(0,231-225)),# room 226 NW stairwell floor 3 (connects to 225 stair below, 208 - hallway on floor)
                                       c(rep(0,46),1,rep(0,44),1,rep(0,135),1, rep(0, 231-228)),# room 227 SE stairwell floor 1 to floor 2 (connects to 47 stair below, 92 - hallway on floor, 228 stair above)
                                       c(rep(0,153),1, rep(0,72),1,0,1,rep(0,231-229)),# room 228 SE stairwell floor 2 to floor 3 (connects to 227 stair below, 154 - hallway on floor, 229 stair above)
                                       c(rep(0,220),1,rep(0,6),1, rep(0,231-228)),# room 229 SE stairwell floor 3 (connects to 228 stair below, 221 - hallway on floor)
                                       c(rep(0,92),1,rep(0,56),1,0,0,1, rep(0,77),1),# room 230 Central stairwell floor 2 to floor 3 (connects to 93 stair below, 150 and 153 - hallway on floor, 231 stair above)
                                       c(rep(0,216),1, rep(0,12),1,0)# room 231 Central stairwell floor3 (connects to 230 stair below, 217 - hallway on floor)
                              )
,nrow = 231,ncol = 231)
diag(church_adjacency_matrix) <- 1
diag(Office_adjacency_matrix) <- 1
diag(Movie_adjacency_matrix) <- 1
diag(University_adjacency_matrix) <- 1

small_bld_3_rooms_graph<- graph_from_adjacency_matrix(small_bld_3_rooms, mode = "undirected")
#take a look
plot(small_bld_3_rooms_graph)

church_graph <- graph_from_adjacency_matrix(church_adjacency_matrix, mode = "undirected")

Office_graph <- graph_from_adjacency_matrix(Office_adjacency_matrix, mode = "undirected")
plot(Office_graph)


#{r set carrying capacities}

small_bld_3_rooms_C <- c(100,200,300)

Church_C <- c(161, #column 1 - main area / Hallway
              161, # column 2 - Hallway
              161, # column 3 - Hallway
              362, #column 4 Main room
              5, # # column 5 Hallway
              7, # Column 6 Hallway - slightly bigger hallway
              5, # Column 7 Hallway
              5, # Column 8 Hallway
              24, # Column 9 Large classroom
              7.7, # column 10 classroom
              7.7, #column 11 classroom
              1, # column 12 janitors closet
              11.8, # column 13 classroom
              3, # column 14 - bathroom
              11.8, # column 15 classroom
              3, # column 16 - bathroom
              7.7, # column 17 classroom
              1, # column 18 - closet
              7.7, # column 19 classroom
              23.75, # column 20 large classroom
              1, # column 21 dressing room
              1, # column 22 dressing room
              1, # column 23 office
              1, # column 24 office
              1, # column 25 office
              1, # column 26 office
              2, # column 27 office - large
              1, # column 28 small bathroom
              1, # column 29 office
              1 # column 30 office
)

Office_C <- c(
  45,# room 1 - Hallway/lobby
  46,# room 2 - Hallway/lobby
  46, # room 3 - Large Common room
  41,# room 4 - Hallway 
  27,# room 5 - Hallway
  18,# room 6 - Hallway
  45,# room 7 - Intermediate room
  9,# room 8 - Hallway- small
  9,# room 9 - Hallway to restrooms and stairway
  10,# room 10 - Conference room - large
  2,# room 11 - Storage room
  8,# room 12 - Conference room - large
  3,# room 13 - Office room
  6,# room 14 - Conference room - small
  3,# room 15 - Office room 
  3,# room 16 - Office room
  2,# room 17 - Kitchenette
  5,# room 18 - Large communal Office space
  3,# room 19 - small communal Office space
  1,# room 20 - Storage room
  rep(3,11), # rooms 21-31 - Offices
  4, # room 32 - Restroom
  4, # room 33 - Restroom
  1 # room 34 - Storage room

)

Movie_C <- c(
  303, # Room 1 hallway
  200, # Room 2 - Hallway
  100, # Room 3 - Hallway 
  71, # Room 4 - Hallway
  13, # Room 5 - Hallway
  13, # Room 6 - Hallway 
  103, # Room 7 - Hallway 
  173, # Room 8 - Hallway
  68, # Room 9 main room
  3, # Room 10 Admin
  1, # Room 11 - admin
  2, # Room 12 - facilities management
  1, # Room 13 - facilities management
  10, # Room 14 - Kitchen
  100, # Room 15 - Theater
  2, # Room 16 - Storage room
  102, # Room 17 - Theater
  12, # Room 18 - restroom
  127, # Room 19 - Theater
  12, # Room 20 - restroom
  1, # Room 21 - family bathroom
  1, # Room 22 - cleaning
  1, # Room 23 - cleaning
  127, # Room 24 - Theater
  1, # Room 25 - facilities
  49, # Room 26 - Theater
  49, # Room 27 - Theater
  73, # Room 28 - Theater
  1, # Room 29 Facilities
  2, # Room 30 - storage
  63, # Room 31 - Theater
  2, # Room 32 - maintenance
  12, # Room 33 - Employee room
  1 # Room 34 - employee bathroom
  
)

University_C <- c(rep(1,5),# rooms 1-5 storage / facilities
                  12, # room 6 research lab
                  rep(3,4), # rooms 7-10 offices
                  10, # room 11 - building services / equipment
                  1, # room 12 - storage
                  rep(4,4), # rooms 13-16 research labs
                  rep(1,5), # rooms 17-21 research lab storage
                  1, # room 22 storage
                  2, # room 23 lab storage
                  1, # room 24 storage
                  1, # room 25 storage - facilities
                  2, # room 26 facilities
                  2, # room 27 facilities
                  rep(4,2), # rooms 28 and 29 - bathrooms
                  rep(1,2), # rooms 30 and 31 - bathroom storage
                  10, # room 32 large unused space - wharehouse?
                  2, # room 33 facilities
                  10, # room 34 - large facilities
                  51, # rooms 35 hallway
                  27,# room 36 hallway
                  15, # room 37 - small hallway outside elevator
                  12, # rooms 38-42 hallways
                  37,# room 39
                  17,# room 40
                  29,# room 41
                  49,# room 42
                  63, # room 43 - small hallway by stairs
                  51, # room 44 stair
                  8, # room 45 elevator
                  8, # room 46 elevator
                  63, # room 47 stair
                  1, # room 48 storage
                  4, # room 49 Teaching lab small room
                  1, # room 50 connecting room to offices/ small common room
                  2, # room 51 storage - teaching lab
                  80, # room 52 Entryway
                  2, # room 53 Building services (elevator?)
                  rep(1,7), # rooms 54-60 offices
                  30, # room 61 - small teaching lab
                  30, # room 62 large teaching lab
                  4, # room 63 small teaching lab
                  30, # room 64 large teaching lab
                  30, # room 65 medium teaching lab
                  2, # room 66 storage
                  4, # room 67 storage
                  50, # room 68 
                  18, # room 69
                  rep(1,3), # rooms 70-72 storage/ building services
                  rep(4,2), # rooms 73 and 74 - bathrooms
                  1, # room 75 storage
                  rep(1,4), # rooms 76-79 small classrooms -- update -- i think they are offices
                  225, # room 80 large lecture hall
                  58,# room 81 - hallway
                  70,# room 82 - hallway
                  7,# room 83 - hallway
                  115,# room 84 - hallway
                  81, # room 85 smaller hallway to teaching labs
                  171,# room 86 
                  131,# room 87
                  225,# room 88
                  8,# room 89
                  3, # room 90 small hallway to storage rooms
                  226, # room 91 hallway with common area 
                  115, # room 92 hallway
                  176, # room 93 - stair way
                  1, # room 94 storage
                  1, # room 95 office
                  rep(1,9), # rooms 96-104 offices
                  rep(30,2), # rooms 105 and 106 small classroom
                  1, # room 107 small research room
                  18, # room 108 research lab
                  44, # room 109 research lab / desks
                  1,# room 110
                  1,# room 111
                  1,# room 112
                  4,# room 113
                  4,# room 114
                  2,# room 115
                  rep(1,2), # rooms 116 and 117 - small storage 
                  rep(1,3), # rooms 118-120 facilities
                  16, # rooms 121 large research lab
                  8, # room 122
                  10,# room 123
                  1, # room 124
                  5, # room 125
                  rep(1,7), # rooms 126-132 offices
                  rep(4,2), # rooms 133 and 134 bathrooms
                  10, # room 135 small common area
                  rep(1,3), # rooms 136-138 small offices
                  40, # room 139 large conference room
                  1, # room 140 facilities
                  109, # rooms 141 hallway
                  107, # room 142 hallway
                  102, # room 143 hallway
                  4, # room 144 hallway
                  42,# room 145 hallway
                  10,# room 146 hallway
                  10,# room 147 hallway 
                  12,# room 148 hallway
                  3,# room 149 hallway
                  36,# room 150 hallway
                  28,# room 151 hallway
                  38,# room 152 hallway
                  40, # room 153 hallway with common area
                  43, # room 154
                  24, # room 155
                  1, # room 156 facilities
                  rep(1,10), # room 157-166 offices
                  30, # room 167 large conference room
                  1, # room 168 research lab storage
                  18, # room 169 large research lab
                  37, # room 170 research desk/lab space
                  rep(1,6 ), # room 171 - 176 research storage 
                  2, # room 177 large storage
                  16, # room 178 large research lab
                  rep(10,2), # room 179 and 180 classroom
                  rep(1,6), # rooms 181 - 186 storage rooms
                  rep(1,7), # rooms 187 - 193 offices
                  rep(4, 2), # rooms 194 and 195 bathrooms
                  1, # room 196 office
                  5, # room 197 admin open space
                  1, # room 198 office
                  5, # room 199 admin open space
                  rep(1,7), # room 200 - 206 office space
                  1, # room 207 storage
                  32, # room 208 hallway
                  60, # room 209 
                  55, # room 210 small hallway
                  64, # room 211 hallway
                  25, # room 212 hallway
                  7, # room 213 hallway 
                  6, # room 214 hallway
                  8, # room 215 hallway
                  3, # room 216 hallway
                  86,# room 217 hallway 
                  4, # room 218 hallway
                  17, # room 219 hallway
                  13, # room 220 hallway
                  7, # room 221 hallway
                  18, # room 222 hallway
                  16, # room 223 hallway
                  58,# room 224 stairwell
                  109,# room 225 stairwell
                  32,# room 226 stairwell
                  115,# room 227 stairwell 
                  43,# room 228 stairwell
                  7,# room 229 stairwell
                  106,# room 230 stairwell
                  86# room 231 stairwell
) #  vector of the carrying capacities for each room - people


length(University_C)

#Now that we have building defined by the adjacency matrix, we need to specify the conditions in the building
#{r Building paramaters}

Max_Building_Capacity_small_bld <- sum(small_bld_3_rooms_C) 
Max_Building_Capacity_Church <- sum(Church_C) 
Max_Building_Capacity_Office <- sum(Office_C) 
Max_Building_Capacity_Movie <- sum(Movie_C) 
Max_Building_Capacity_University <- sum(University_C) 

Prop_full <- 0.8 # how full do we want our building capacity to be
Adj_Max_Building_Capacity_small_bld <- round(Prop_full*Max_Building_Capacity_small_bld,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Church <- round(Prop_full*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office <- round(Prop_full*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie <- round(Prop_full*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University <- round(Prop_full*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_small_bld <- Adj_Max_Building_Capacity_small_bld/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Church <- Adj_Max_Building_Capacity_Church/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office <- Adj_Max_Building_Capacity_Office/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie <- Adj_Max_Building_Capacity_Movie/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University <- Adj_Max_Building_Capacity_University/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)


day <- 31 

#What day to extract results from the community level model
#N_rooms <- length(Church_C)


Small_bld_setup <- Bld_setup_func_v2(Community_output = Community_output,day = day,delt = delt_small_bld,N_rooms =length(small_bld_3_rooms_C),C_x = small_bld_3_rooms_C,N_total = Adj_Max_Building_Capacity_small_bld)
Small_bld_setup

Church_setup <- Bld_setup_func_v2(Community_output = Community_output,day = day,delt = delt_Church, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church)
Office_setup <- Bld_setup_func_v2(Community_output = Community_output,day = day,delt = delt_Office, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office)
Movie_setup <- Bld_setup_func_v2(Community_output = Community_output,day = day,delt = delt_Movie, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie)
University_setup <- Bld_setup_func_v2(Community_output = Community_output,day = day,delt = delt_University, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University)

#


#**Checking initial conditions**
#   Now lets check that our building setup function did what I wanted it to. 
# Each room's proportion should sum to 1. And if we convert back to numbers using N_x the numbers should equal
#{r props sum to 1}
#Small_bld_setup %>% summarise(total_prop = sum(S)+sum(I)+sum(R),build_pop = sum(N_x))
#
#props sum to 1 so that is good

#{r props in building match community}
#prop of susceptible in building
# sum(Small_bld_setup$S)
# #prop of susceptible in community
# Community_output$S[day]/(Community_output$S[day]+Community_output$I[day]+Community_output$R[day])
# 
# #prop of infectious in building
# sum(Small_bld_setup$I)
# #prop of infectious in community
# Community_output$I[day]/(Community_output$S[day]+Community_output$I[day]+Community_output$R[day])
# 
# #prop of recovered in building
# sum(Small_bld_setup$R)
# #prop of susceptible in community
# Community_output$R[day]/(Community_output$S[day]+Community_output$I[day]+Community_output$R[day])
# 

#
#Okay, great all of that checks out

#**Now people and particles need to move (Transition matrices)**



#Church_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = church_adjacency_matrix, N_rooms = length())
Small_bld_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = small_bld_3_rooms,N_rooms = length(small_bld_3_rooms_C))

Church_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = church_adjacency_matrix, N_rooms = length(Church_C))
Office_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = Office_adjacency_matrix, N_rooms = length(Office_C))
Movie_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = Movie_adjacency_matrix, N_rooms = length(Movie_C))
University_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = University_adjacency_matrix, N_rooms = length(University_C))

sum(Small_bld_T_mov[1,])
sum(Small_bld_T_mov[2,])
sum(Small_bld_T_mov[3,])
# sum(Church_T_mov[4,])
# sum(Church_T_mov[5,])
r1_in <- sum(Small_bld_T_mov[,1])
r1_out <- sum(Small_bld_T_mov[1,])
r1_change <- r1_in-r1_out
r2_in <- sum(Small_bld_T_mov[,2])
r2_out <- sum(Small_bld_T_mov[2,])
r2_change <- r2_in-r2_out
r3_in <- sum(Small_bld_T_mov[,3])
r3_out <- sum(Small_bld_T_mov[3,])
r3_change <- r3_in - r3_out
# We can use the same function for the particle matrix
#Church_theta_mov <- Create_T_Matrix(adjacency_matrix_to_use = church_adjacency_matrix)
Small_bld_theta_mov <- Create_T_Matrix(adjacency_matrix_to_use = small_bld_3_rooms,N_rooms = length(small_bld_3_rooms_C))

Church_theta_mov <- Create_T_Matrix(adjacency_matrix_to_use = church_adjacency_matrix, N_rooms = length(Church_C))
Office_theta_mov <- Create_T_Matrix(adjacency_matrix_to_use = Office_adjacency_matrix, N_rooms = length(Office_C))
Movie_theta_mov <- Create_T_Matrix(adjacency_matrix_to_use = Movie_adjacency_matrix, N_rooms = length(Movie_C))
University_theta_mov <- Create_T_Matrix(adjacency_matrix_to_use = University_adjacency_matrix, N_rooms = length(University_C))

#

#**Now set parameters for the model**
  
  
  #{r parameters}
parms <-data.frame(s=2000,a=396, d=1.81,bet_bar = 25569.9504, bet_hat =67320)
# s = shedding, a = absorption, d = decay, lam = scalar for room capacities
#N_b <- Adj_Max_Building_Capacity #Building population size
Maxtime <- 24*3
times <- seq(from = 0, to = Maxtime, by = 0.2)
m <- 5 # number of equations per room
#


Small_bld_Init_conds_v3 <-c(N_x =Small_bld_setup$N_x, P = Small_bld_setup$P)

Church_Init_conds <- c(N_x = Church_setup$N_x, P = Church_setup$P)
Office_Init_conds <- c(N_x = Office_setup$N_x, P = Office_setup$P)
Movie_Init_conds <- c(N_x = Movie_setup$N_x, P = Movie_setup$P)
University_Init_conds <- c(N_x = University_setup$N_x, P = University_setup$P)


Small_bld_output_v3 <- data.frame(lsoda(y = Small_bld_Init_conds_v3, func = Particle_model_v3,times = times,
                                        parms = parms,
                                        adjacency_matrix_to_use=small_bld_3_rooms,
                                        theta_mov =Small_bld_theta_mov,
                                        T_mov = Small_bld_T_mov, 
                                        C_x=small_bld_3_rooms_C,
                                        N_b = Adj_Max_Building_Capacity_small_bld,
                                        N_rooms = length(small_bld_3_rooms_C),
                                        Ib_prop = Small_bld_setup$Ib_prop[1]))


data_clean <- Small_bld_output_v3%>% pivot_longer(cols = !time,
                                                  names_to = c("State", "Room"),
                                                  names_pattern = "(P|N_x)(\\d+)",
                                                  values_to = "Number")


data_clean %>% filter(State == "P") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of infectious particles")

data_clean %>% filter(State == "N_x") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of people")

data_clean %>% filter(State == "N_x") %>% ungroup() %>% group_by(time) %>% summarise(Total_prop=sum(Number))%>%
  ggplot(aes(x=time, y=Total_prop))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in rooms within a building")
temp <- data_clean
temp_C_x <- data.frame(Room = as.character(seq(1:length(small_bld_3_rooms_C))),C_x = c(small_bld_3_rooms_C))
temp_data <- left_join(temp,temp_C_x, by = "Room")

ratio_data <- temp_data %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Small_bld_setup$Ib_prop[1],Numerator = s*N_x*Ib_prop*(C_x*bet_bar)*bet_hat , Denominator =(a*N_x*P*(P)), ratio = Numerator/Denominator ) 
#temp_data %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
#+   group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Small_bld_setup$Ib_prop[1],shed = s*N_x*Ib_prop*C_x*bet_bar , Denominator =(a*P), ratio = shed/Denominator ) 
#> #(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
#(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
ggplot(ratio_data, aes(x = time, y = ratio, color = Room))+geom_line()
#


Church_output <- data.frame(lsoda(y = Church_Init_conds, func = Particle_model_v3,times = times,
                                        parms = parms,
                                        adjacency_matrix_to_use=church_adjacency_matrix,
                                        theta_mov = Church_theta_mov,
                                        T_mov = Church_T_mov, 
                                        C_x= Church_C,
                                        N_b = Adj_Max_Building_Capacity_Church,
                                        N_rooms = length(Church_C),
                                        Ib_prop = Church_setup$Ib_prop[1]))
#Church_data <- write.table(Church_output,file = "Data/Church_output_Nov24_parms1.text",sep = " ")
#Church_data <- read.table(file = "Data/Church_output_Nov24_parms1.text",sep = " ")


Church_data_clean <- Church_output %>% pivot_longer(cols = !time,
                                                  names_to = c("State", "Room"),
                                                  names_pattern = "(P|N_x)(\\d+)",
                                                  values_to = "Number")


Church_data_clean %>% filter(State == "P") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of infectious particles")

Church_data_clean %>% filter(State == "N_x") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of people")

Church_data_clean %>% filter(State == "N_x") %>% ungroup() %>% group_by(time) %>% summarise(Total_prop=sum(Number))%>%
  ggplot(aes(x=time, y=Total_prop))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in rooms within a building")
temp_Church <- Church_data_clean
temp_C_x_Church <- data.frame(Room = as.character(seq(1:length(Church_C))),C_x = c(Church_C))
temp_data_Church<- left_join(temp_Church,temp_C_x_Church, by = "Room")

ratio_data_Church <- temp_data_Church %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Church_setup$Ib_prop[1],Part_dens_room=(P/(C_x*bet_bar)),Part_dens_breath = P/bet_hat,Numerator = s*N_x*Ib_prop, Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_hat)), ratio = Numerator/Denominator ) 
#temp_data %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
#+   group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Small_bld_setup$Ib_prop[1],shed = s*N_x*Ib_prop*C_x*bet_bar , Denominator =(a*P), ratio = shed/Denominator ) 
#> #(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
#(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
ggplot(ratio_data_Church, aes(x = time, y = ratio, color = Room))+geom_line()+xlim(48,72)#+ylim(0,200)
#

Office_output <- data.frame(lsoda(y = Office_Init_conds, func = Particle_model_v3,times = times,
                                  parms = parms,
                                  adjacency_matrix_to_use=Office_adjacency_matrix,
                                  theta_mov = Office_theta_mov,
                                  T_mov = Office_T_mov, 
                                  C_x= Office_C,
                                  N_b = Adj_Max_Building_Capacity_Office,
                                  N_rooms = length(Office_C),
                                  Ib_prop = Office_setup$Ib_prop[1]))

# Office_data <- write.table(Office_output,file = "Data/Office_output_Nov24_parms1.text",sep = " ")
# Office_data <- read.table(file = "Data/Office_output_Nov24_parms1.text",sep = " ")

Office_data_clean <- Office_output %>% pivot_longer(cols = !time,
                                                    names_to = c("State", "Room"),
                                                    names_pattern = "(P|N_x)(\\d+)",
                                                    values_to = "Number")


Office_data_clean %>% filter(State == "P") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of infectious particles")

Office_data_clean %>% filter(State == "N_x") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of people")

Office_data_clean %>% filter(State == "N_x") %>% ungroup() %>% group_by(time) %>% summarise(Total_prop=sum(Number))%>%
  ggplot(aes(x=time, y=Total_prop))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in rooms within a building")
temp_Office <- Office_data_clean
temp_C_x_Office <- data.frame(Room = as.character(seq(1:length(Office_C))),C_x = c(Office_C))
temp_data_Office<- left_join(temp_Office,temp_C_x_Office, by = "Room")

ratio_data_Office <- temp_data_Office %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Small_bld_setup$Ib_prop[1],Numerator = s*N_x*Ib_prop*(C_x*bet_bar)*bet_hat , Denominator =(a*N_x*P*(P)), ratio = Numerator/Denominator ) 
#temp_data %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
#+   group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Small_bld_setup$Ib_prop[1],shed = s*N_x*Ib_prop*C_x*bet_bar , Denominator =(a*P), ratio = shed/Denominator ) 
#> #(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
#(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
ggplot(ratio_data_Office, aes(x = time, y = ratio, color = Room))+geom_line()
#

Movie_output <- data.frame(lsoda(y = Movie_Init_conds, func = Particle_model_v3,times = times,
                                  parms = parms,
                                  adjacency_matrix_to_use=Movie_adjacency_matrix,
                                  theta_mov = Movie_theta_mov,
                                  T_mov = Movie_T_mov, 
                                  C_x= Movie_C,
                                  N_b = Adj_Max_Building_Capacity_Movie,
                                  N_rooms = length(Movie_C),
                                  Ib_prop = Movie_setup$Ib_prop[1]))
# Movie_data <- write.table(Movie_output,file = "Data/Movie_output_Nov24_parms1.text",sep = " ")
# Movie_data <- read.table(file = "Data/Movie_output_Nov24_parms1.text",sep = " ")


Movie_data_clean <- Movie_output %>% pivot_longer(cols = !time,
                                                    names_to = c("State", "Room"),
                                                    names_pattern = "(P|N_x)(\\d+)",
                                                    values_to = "Number")


Movie_data_clean %>% filter(State == "P") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of infectious particles")

Movie_data_clean %>% filter(State == "N_x") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of people")

Movie_data_clean %>% filter(State == "N_x") %>% ungroup() %>% group_by(time) %>% summarise(Total_prop=sum(Number))%>%
  ggplot(aes(x=time, y=Total_prop))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in rooms within a building")
temp_Movie <- Movie_data_clean
temp_C_x_Movie <- data.frame(Room = as.character(seq(1:length(Movie_C))),C_x = c(Movie_C))
temp_data_Movie<- left_join(temp_Movie,temp_C_x_Movie, by = "Room")

ratio_data_Movie <- temp_data_Movie %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Small_bld_setup$Ib_prop[1],Numerator = s*N_x*Ib_prop*(C_x*bet_bar)*bet_hat , Denominator =(a*N_x*P*(P)), ratio = Numerator/Denominator ) 
#temp_data %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
#+   group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Small_bld_setup$Ib_prop[1],shed = s*N_x*Ib_prop*C_x*bet_bar , Denominator =(a*P), ratio = shed/Denominator ) 
#> #(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
#(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
ggplot(ratio_data_Movie, aes(x = time, y = ratio, color = Room))+geom_line()
#

University_output <- data.frame(lsoda(y = University_Init_conds, func = Particle_model_v3,times = times,
                                  parms = parms,
                                  adjacency_matrix_to_use=University_adjacency_matrix,
                                  theta_mov = University_theta_mov,
                                  T_mov = University_T_mov, 
                                  C_x= University_C,
                                  N_b = Adj_Max_Building_Capacity_University,
                                  N_rooms = length(University_C),
                                  Ib_prop = University_setup$Ib_prop[1]))

University_data <- write.table(University_output,file = "Data/University_output_Nov24_parms1.text",sep = " ")
University_data <- read.table(file = "Data/University_output_Nov24_parms1.text",sep = " ")

University_data_clean <- University_output %>% pivot_longer(cols = !time,
                                                    names_to = c("State", "Room"),
                                                    names_pattern = "(P|N_x)(\\d+)",
                                                    values_to = "Number")


University_data_clean %>% filter(State == "P") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of infectious particles")

University_data_clean %>% filter(State == "N_x") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number")+ggtitle("Number of people")

University_data_clean %>% filter(State == "N_x") %>% ungroup() %>% group_by(time) %>% summarise(Total_prop=sum(Number))%>%
  ggplot(aes(x=time, y=Total_prop))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in rooms within a building")
temp_University <- University_data_clean
temp_C_x_University <- data.frame(Room = as.character(seq(1:length(University_C))),C_x = c(University_C))
temp_data_University<- left_join(temp_University,temp_C_x_University, by = "Room")

ratio_data_University <- temp_data_University %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = University_setup$Ib_prop[1],Numerator = s*N_x*Ib_prop*(C_x*bet_bar)*bet_hat , Denominator =(a*N_x*P*(P)), ratio = Numerator/Denominator ) 
#temp_data %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
#+   group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Small_bld_setup$Ib_prop[1],shed = s*N_x*Ib_prop*C_x*bet_bar , Denominator =(a*P), ratio = shed/Denominator ) 
#> #(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
#(parms$s*N_x*Small_bld_setup$Ib_prop[1]*N_x*small_bld_3_rooms_C[Room]*parms$lam)
ggplot(ratio_data_University, aes(x = time, y = ratio, color = Room))+geom_line()
#



