
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
library(kableExtra)


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

Prop_full <- data.frame(ten = 0.1,twenty = 0.2, thirty = 0.3, forty = 0.4, fifty= 0.5, sixty = 0.6, seventy = 0.7, eighty = 0.8, ninety = 0.9, hundred = 1.0) # how full do we want our building capacity to be
#Adj_Max_Building_Capacity_small_bld <- round(Prop_full$ten[1]*Max_Building_Capacity_small_bld,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

#ten percent capacity
Adj_Max_Building_Capacity_Church_ten <- round(Prop_full$ten[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_ten <- round(Prop_full$ten[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_ten <- round(Prop_full$ten[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_ten <- round(Prop_full$ten[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_ten<- Adj_Max_Building_Capacity_Church_ten/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_ten <- Adj_Max_Building_Capacity_Office_ten/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_ten <- Adj_Max_Building_Capacity_Movie_ten/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_ten <- Adj_Max_Building_Capacity_University_ten/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)


#twenty percent capacity
Adj_Max_Building_Capacity_Church_twenty <- round(Prop_full$twenty[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_twenty <- round(Prop_full$twenty[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_twenty <- round(Prop_full$twenty[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_twenty <- round(Prop_full$twenty[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_twenty<- Adj_Max_Building_Capacity_Church_twenty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_twenty <- Adj_Max_Building_Capacity_Office_twenty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_twenty <- Adj_Max_Building_Capacity_Movie_twenty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_twenty <- Adj_Max_Building_Capacity_University_twenty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

#thirty percent capacity
Adj_Max_Building_Capacity_Church_thirty <- round(Prop_full$thirty[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_thirty <- round(Prop_full$thirty[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_thirty <- round(Prop_full$thirty[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_thirty <- round(Prop_full$thirty[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_thirty<- Adj_Max_Building_Capacity_Church_thirty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_thirty <- Adj_Max_Building_Capacity_Office_thirty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_thirty <- Adj_Max_Building_Capacity_Movie_thirty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_thirty <- Adj_Max_Building_Capacity_University_thirty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

#forty percent capacity
Adj_Max_Building_Capacity_Church_forty<- round(Prop_full$forty[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_forty <- round(Prop_full$forty[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_forty <- round(Prop_full$forty[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_forty<- round(Prop_full$forty[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_forty<- Adj_Max_Building_Capacity_Church_forty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_forty <- Adj_Max_Building_Capacity_Office_forty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_forty <- Adj_Max_Building_Capacity_Movie_forty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_forty <- Adj_Max_Building_Capacity_University_forty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

#fifty percent capacity
Adj_Max_Building_Capacity_Church_fifty <- round(Prop_full$fifty[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_fifty <- round(Prop_full$fifty[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_fifty <- round(Prop_full$fifty[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_fifty <- round(Prop_full$fifty[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_fifty<- Adj_Max_Building_Capacity_Church_fifty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_fifty <- Adj_Max_Building_Capacity_Office_fifty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_fifty <- Adj_Max_Building_Capacity_Movie_fifty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_fifty <- Adj_Max_Building_Capacity_University_fifty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

#sixty percent capacity
Adj_Max_Building_Capacity_Church_sixty <- round(Prop_full$sixty[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_sixty <- round(Prop_full$sixty[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_sixty <- round(Prop_full$sixty[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_sixty <- round(Prop_full$sixty[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_sixty<- Adj_Max_Building_Capacity_Church_sixty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_sixty <- Adj_Max_Building_Capacity_Office_sixty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_sixty <- Adj_Max_Building_Capacity_Movie_sixty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_sixty <- Adj_Max_Building_Capacity_University_sixty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

#seventy percent capacity
Adj_Max_Building_Capacity_Church_seventy <- round(Prop_full$seventy[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_seventy <- round(Prop_full$seventy[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_seventy <- round(Prop_full$seventy[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_seventy <- round(Prop_full$seventy[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_seventy<- Adj_Max_Building_Capacity_Church_seventy/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_seventy <- Adj_Max_Building_Capacity_Office_seventy/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_seventy <- Adj_Max_Building_Capacity_Movie_seventy/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_seventy <- Adj_Max_Building_Capacity_University_seventy/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

#eighty percent capacity
Adj_Max_Building_Capacity_Church_eighty <- round(Prop_full$eighty[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_eighty <- round(Prop_full$eighty[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_eighty <- round(Prop_full$eighty[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_eighty <- round(Prop_full$eighty[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_eighty <- Adj_Max_Building_Capacity_Church_eighty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_eighty <- Adj_Max_Building_Capacity_Office_eighty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_eighty <- Adj_Max_Building_Capacity_Movie_eighty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_eighty <- Adj_Max_Building_Capacity_University_eighty/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

#ninety percent capacity
Adj_Max_Building_Capacity_Church_ninety <- round(Prop_full$ninety[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_ninety <- round(Prop_full$ninety[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_ninety <- round(Prop_full$ninety[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_ninety <- round(Prop_full$ninety[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_ninety<- Adj_Max_Building_Capacity_Church_ninety/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_ninety <- Adj_Max_Building_Capacity_Office_ninety/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_ninety <- Adj_Max_Building_Capacity_Movie_ninety/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_ninety <- Adj_Max_Building_Capacity_University_ninety/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

#one-hundred percent capacity
Adj_Max_Building_Capacity_Church_hundred <- round(Prop_full$hundred[1]*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Office_hundred <- round(Prop_full$hundred[1]*Max_Building_Capacity_Office,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_Movie_hundred <- round(Prop_full$hundred[1]*Max_Building_Capacity_Movie,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity
Adj_Max_Building_Capacity_University_hundred <- round(Prop_full$hundred[1]*Max_Building_Capacity_University,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

delt_Church_hundred<- Adj_Max_Building_Capacity_Church_hundred/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Office_hundred <- Adj_Max_Building_Capacity_Office_hundred/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_Movie_hundred <- Adj_Max_Building_Capacity_Movie_hundred/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)
delt_University_hundred <- Adj_Max_Building_Capacity_University_hundred/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

max_I <-Community_output %>% filter(I == max(I))

beg_day <- 25 # estimated inflection point 
mid_day <- max_I$time
end_day <- 70 # lower than peak but still much greater than the begining
## have ten different levels that buildings are full.
# Also going to consider when we are at the begining, middle and peak of an outbreak
# this should change the proportion of infecteds that are in the building

#Small_bld_setup <- Bld_setup_func_v2(Community_output = Community_output,day = day,delt = delt_small_bld,N_rooms =length(small_bld_3_rooms_C),C_x = small_bld_3_rooms_C,N_total = Adj_Max_Building_Capacity_small_bld)
#Small_bld_setup

Church_setup_ten_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_ten, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_ten)
Office_setup_ten_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_ten, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_ten)
Movie_setup_ten_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_ten, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_ten)
University_setup_ten_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_ten, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_ten)

Church_setup_twenty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_twenty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_twenty)
Office_setup_twenty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_twenty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_twenty)
Movie_setup_twenty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_twenty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_twenty)
University_setup_twenty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_twenty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_twenty)

Church_setup_thirty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_thirty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_thirty)
Office_setup_thirty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_thirty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_thirty)
Movie_setup_thirty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_thirty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_thirty)
University_setup_thirty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_thirty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_thirty)

Church_setup_forty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_forty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_forty)
Office_setup_forty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_forty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_forty)
Movie_setup_forty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_forty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_forty)
University_setup_forty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_forty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_forty)

Church_setup_fifty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_fifty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_fifty)
Office_setup_fifty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_fifty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_fifty)
Movie_setup_fifty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_fifty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_fifty)
University_setup_fifty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_fifty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_fifty)

Church_setup_sixty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_sixty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_sixty)
Office_setup_sixty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_sixty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_sixty)
Movie_setup_sixty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_sixty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_sixty)
University_setup_sixty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_sixty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_sixty)

Church_setup_seventy_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_seventy, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_seventy)
Office_setup_seventy_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_seventy, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_seventy)
Movie_setup_seventy_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_seventy, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_seventy)
University_setup_seventy_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_seventy, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_seventy)

Church_setup_eighty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_eighty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_eighty)
Office_setup_eighty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_eighty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_eighty)
Movie_setup_eighty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_eighty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_eighty)
University_setup_eighty_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_eighty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_eighty)

Church_setup_ninety_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_ninety, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_ninety)
Office_setup_ninety_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_ninety, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_ninety)
Movie_setup_ninety_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_ninety, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_ninety)
University_setup_ninety_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_ninety, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_ninety)

Church_setup_hundred_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Church_hundred, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_hundred)
Office_setup_hundred_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Office_hundred, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_hundred)
Movie_setup_hundred_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_Movie_hundred, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_hundred)
University_setup_hundred_beg <- Bld_setup_func_v2(Community_output = Community_output,day = beg_day,delt = delt_University_hundred, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_hundred)

#### now peak infecteds ####
Church_setup_ten_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_ten, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_ten)
Office_setup_ten_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_ten, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_ten)
Movie_setup_ten_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_ten, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_ten)
University_setup_ten_mid <- Bld_setup_func_v2(Community_output = Community_output, day = mid_day,delt = delt_University_ten, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_ten)

Church_setup_twenty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_twenty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_twenty)
Office_setup_twenty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_twenty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_twenty)
Movie_setup_twenty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_twenty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_twenty)
University_setup_twenty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_University_twenty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_twenty)

Church_setup_thirty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_thirty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_thirty)
Office_setup_thirty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_thirty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_thirty)
Movie_setup_thirty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_thirty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_thirty)
University_setup_thirty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_University_thirty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_thirty)

Church_setup_forty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_forty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_forty)
Office_setup_forty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_forty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_forty)
Movie_setup_forty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_forty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_forty)
University_setup_forty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_University_forty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_forty)

Church_setup_fifty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_fifty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_fifty)
Office_setup_fifty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_fifty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_fifty)
Movie_setup_fifty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_fifty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_fifty)
University_setup_fifty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_University_fifty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_fifty)

Church_setup_sixty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_sixty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_sixty)
Office_setup_sixty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_sixty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_sixty)
Movie_setup_sixty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_sixty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_sixty)
University_setup_sixty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_University_sixty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_sixty)

Church_setup_seventy_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_seventy, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_seventy)
Office_setup_seventy_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_seventy, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_seventy)
Movie_setup_seventy_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_seventy, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_seventy)
University_setup_seventy_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_University_seventy, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_seventy)

Church_setup_eighty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_eighty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_eighty)
Office_setup_eighty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_eighty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_eighty)
Movie_setup_eighty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_eighty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_eighty)
University_setup_eighty_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_University_eighty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_eighty)

Church_setup_ninety_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_ninety, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_ninety)
Office_setup_ninety_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_ninety, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_ninety)
Movie_setup_ninety_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_ninety, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_ninety)
University_setup_ninety_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_University_ninety, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_ninety)

Church_setup_hundred_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Church_hundred, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_hundred)
Office_setup_hundred_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Office_hundred, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_hundred)
Movie_setup_hundred_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_Movie_hundred, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_hundred)
University_setup_hundred_mid <- Bld_setup_func_v2(Community_output = Community_output,day = mid_day,delt = delt_University_hundred, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_hundred)


### Now at the end of an outbreak
Church_setup_ten_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_ten, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_ten)
Office_setup_ten_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_ten, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_ten)
Movie_setup_ten_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_ten, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_ten)
University_setup_ten_end <- Bld_setup_func_v2(Community_output = Community_output, day = end_day,delt = delt_University_ten, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_ten)

Church_setup_twenty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_twenty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_twenty)
Office_setup_twenty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_twenty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_twenty)
Movie_setup_twenty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_twenty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_twenty)
University_setup_twenty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_University_twenty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_twenty)

Church_setup_thirty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_thirty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_thirty)
Office_setup_thirty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_thirty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_thirty)
Movie_setup_thirty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_thirty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_thirty)
University_setup_thirty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_University_thirty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_thirty)

Church_setup_forty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_forty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_forty)
Office_setup_forty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_forty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_forty)
Movie_setup_forty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_forty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_forty)
University_setup_forty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_University_forty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_forty)

Church_setup_fifty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_fifty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_fifty)
Office_setup_fifty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_fifty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_fifty)
Movie_setup_fifty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_fifty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_fifty)
University_setup_fifty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_University_fifty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_fifty)

Church_setup_sixty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_sixty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_sixty)
Office_setup_sixty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_sixty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_sixty)
Movie_setup_sixty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_sixty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_sixty)
University_setup_sixty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_University_sixty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_sixty)

Church_setup_seventy_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_seventy, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_seventy)
Office_setup_seventy_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_seventy, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_seventy)
Movie_setup_seventy_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_seventy, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_seventy)
University_setup_seventy_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_University_seventy, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_seventy)

Church_setup_eighty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_eighty, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_eighty)
Office_setup_eighty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_eighty, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_eighty)
Movie_setup_eighty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_eighty, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_eighty)
University_setup_eighty_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_University_eighty, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_eighty)

Church_setup_ninety_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_ninety, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_ninety)
Office_setup_ninety_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_ninety, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_ninety)
Movie_setup_ninety_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_ninety, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_ninety)
University_setup_ninety_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_University_ninety, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_ninety)

Church_setup_hundred_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Church_hundred, N_rooms =length(Church_C),C_x = Church_C, N_total = Adj_Max_Building_Capacity_Church_hundred)
Office_setup_hundred_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Office_hundred, N_rooms =length(Office_C),C_x = Office_C, N_total = Adj_Max_Building_Capacity_Office_hundred)
Movie_setup_hundred_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_Movie_hundred, N_rooms =length(Movie_C),C_x = Movie_C, N_total = Adj_Max_Building_Capacity_Movie_hundred)
University_setup_hundred_end <- Bld_setup_func_v2(Community_output = Community_output,day = end_day,delt = delt_University_hundred, N_rooms =length(University_C),C_x = University_C, N_total = Adj_Max_Building_Capacity_University_hundred)

#**Now people and particles need to move (Transition matrices)**



#Church_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = church_adjacency_matrix, N_rooms = length())
Small_bld_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = small_bld_3_rooms,N_rooms = length(small_bld_3_rooms_C))

Church_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = church_adjacency_matrix, N_rooms = length(Church_C))
Office_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = Office_adjacency_matrix, N_rooms = length(Office_C))
Movie_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = Movie_adjacency_matrix, N_rooms = length(Movie_C))
University_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = University_adjacency_matrix, N_rooms = length(University_C))


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
Maxtime <- 24*5
times <- seq(from = 0, to = Maxtime, by = 0.2)
m <- 5 # number of equations per room
#


#Small_bld_Init_conds_v3 <-c(N_x =Small_bld_setup$N_x, P = Small_bld_setup$P)
#### Initial conditions all prop full options -- begining of outbreak
Church_Init_conds_ten_beg <- c(N_x = Church_setup_ten_beg$N_x, P = Church_setup_ten_beg$P)
Office_Init_conds_ten_beg <- c(N_x = Office_setup_ten_beg$N_x, P = Office_setup_ten_beg$P)
Movie_Init_conds_ten_beg <- c(N_x = Movie_setup_ten_beg$N_x, P = Movie_setup_ten_beg$P)
University_Init_conds_ten_beg <- c(N_x = University_setup_ten_beg$N_x, P = University_setup_ten_beg$P)

Church_Init_conds_twenty_beg <- c(N_x = Church_setup_twenty_beg$N_x, P = Church_setup_twenty_beg$P)
Office_Init_conds_twenty_beg <- c(N_x = Office_setup_twenty_beg$N_x, P = Office_setup_twenty_beg$P)
Movie_Init_conds_twenty_beg <- c(N_x = Movie_setup_twenty_beg$N_x, P = Movie_setup_twenty_beg$P)
University_Init_conds_twenty_beg <- c(N_x = University_setup_twenty_beg$N_x, P = University_setup_twenty_beg$P)

Church_Init_conds_thirty_beg <- c(N_x = Church_setup_thirty_beg$N_x, P = Church_setup_thirty_beg$P)
Office_Init_conds_thirty_beg <- c(N_x = Office_setup_thirty_beg$N_x, P = Office_setup_thirty_beg$P)
Movie_Init_conds_thirty_beg <- c(N_x = Movie_setup_thirty_beg$N_x, P = Movie_setup_thirty_beg$P)
University_Init_conds_thirty_beg <- c(N_x = University_setup_thirty_beg$N_x, P = University_setup_thirty_beg$P)

Church_Init_conds_forty_beg <- c(N_x = Church_setup_forty_beg$N_x, P = Church_setup_forty_beg$P)
Office_Init_conds_forty_beg <- c(N_x = Office_setup_forty_beg$N_x, P = Office_setup_forty_beg$P)
Movie_Init_conds_forty_beg <- c(N_x = Movie_setup_forty_beg$N_x, P = Movie_setup_forty_beg$P)
University_Init_conds_forty_beg <- c(N_x = University_setup_forty_beg$N_x, P = University_setup_forty_beg$P)

Church_Init_conds_fifty_beg <- c(N_x = Church_setup_ten_beg$N_x, P = Church_setup_fifty_beg$P)
Office_Init_conds_fifty_beg <- c(N_x = Office_setup_ten_beg$N_x, P = Office_setup_fifty_beg$P)
Movie_Init_conds_fifty_beg <- c(N_x = Movie_setup_ten_beg$N_x, P = Movie_setup_fifty_beg$P)
University_Init_conds_fifty_beg <- c(N_x = University_setup_ten_beg$N_x, P = University_setup_fifty_beg$P)

Church_Init_conds_sixty_beg <- c(N_x = Church_setup_sixty_beg$N_x, P = Church_setup_sixty_beg$P)
Office_Init_conds_sixty_beg <- c(N_x = Office_setup_sixty_beg$N_x, P = Office_setup_sixty_beg$P)
Movie_Init_conds_sixty_beg <- c(N_x = Movie_setup_sixty_beg$N_x, P = Movie_setup_sixty_beg$P)
University_Init_conds_sixty_beg <- c(N_x = University_setup_sixty_beg$N_x, P = University_setup_sixty_beg$P)

Church_Init_conds_seventy_beg <- c(N_x = Church_setup_seventy_beg$N_x, P = Church_setup_seventy_beg$P)
Office_Init_conds_seventy_beg <- c(N_x = Office_setup_seventy_beg$N_x, P = Office_setup_seventy_beg$P)
Movie_Init_conds_seventy_beg <- c(N_x = Movie_setup_seventy_beg$N_x, P = Movie_setup_seventy_beg$P)
University_Init_conds_seventy_beg <- c(N_x = University_setup_seventy_beg$N_x, P = University_setup_seventy_beg$P)

Church_Init_conds_eighty_beg <- c(N_x = Church_setup_eighty_beg$N_x, P = Church_setup_eighty_beg$P)
Office_Init_conds_eighty_beg <- c(N_x = Office_setup_eighty_beg$N_x, P = Office_setup_eighty_beg$P)
Movie_Init_conds_eighty_beg <- c(N_x = Movie_setup_eighty_beg$N_x, P = Movie_setup_eighty_beg$P)
University_Init_conds_eighty_beg <- c(N_x = University_setup_eighty_beg$N_x, P = University_setup_eighty_beg$P)

Church_Init_conds_ninety_beg <- c(N_x = Church_setup_ninety_beg$N_x, P = Church_setup_ninety_beg$P)
Office_Init_conds_ninety_beg <- c(N_x = Office_setup_ninety_beg$N_x, P = Office_setup_ninety_beg$P)
Movie_Init_conds_ninety_beg <- c(N_x = Movie_setup_ninety_beg$N_x, P = Movie_setup_ninety_beg$P)
University_Init_conds_ninety_beg <- c(N_x = University_setup_ninety_beg$N_x, P = University_setup_ninety_beg$P)

Church_Init_conds_hundred_beg <- c(N_x = Church_setup_hundred_beg$N_x, P = Church_setup_hundred_beg$P)
Office_Init_conds_hundred_beg <- c(N_x = Office_setup_hundred_beg$N_x, P = Office_setup_hundred_beg$P)
Movie_Init_conds_hundred_beg <- c(N_x = Movie_setup_hundred_beg$N_x, P = Movie_setup_hundred_beg$P)
University_Init_conds_hundred_beg <- c(N_x = University_setup_hundred_beg$N_x, P = University_setup_hundred_beg$P)

#### Initial conditions all prop full options -- end of outbreak
Church_Init_conds_ten_mid <- c(N_x = Church_setup_ten_mid$N_x, P = Church_setup_ten_mid$P)
Office_Init_conds_ten_mid <- c(N_x = Office_setup_ten_mid$N_x, P = Office_setup_ten_mid$P)
Movie_Init_conds_ten_mid <- c(N_x = Movie_setup_ten_mid$N_x, P = Movie_setup_ten_mid$P)
University_Init_conds_ten_mid <- c(N_x = University_setup_ten_mid$N_x, P = University_setup_ten_mid$P)

Church_Init_conds_twenty_mid <- c(N_x = Church_setup_twenty_mid$N_x, P = Church_setup_twenty_mid$P)
Office_Init_conds_twenty_mid <- c(N_x = Office_setup_twenty_mid$N_x, P = Office_setup_twenty_mid$P)
Movie_Init_conds_twenty_mid <- c(N_x = Movie_setup_twenty_mid$N_x, P = Movie_setup_twenty_mid$P)
University_Init_conds_twenty_mid <- c(N_x = University_setup_twenty_mid$N_x, P = University_setup_twenty_mid$P)

Church_Init_conds_thirty_mid <- c(N_x = Church_setup_thirty_mid$N_x, P = Church_setup_thirty_mid$P)
Office_Init_conds_thirty_mid <- c(N_x = Office_setup_thirty_mid$N_x, P = Office_setup_thirty_mid$P)
Movie_Init_conds_thirty_mid <- c(N_x = Movie_setup_thirty_mid$N_x, P = Movie_setup_thirty_mid$P)
University_Init_conds_thirty_mid <- c(N_x = University_setup_thirty_mid$N_x, P = University_setup_thirty_mid$P)

Church_Init_conds_forty_mid <- c(N_x = Church_setup_forty_mid$N_x, P = Church_setup_forty_mid$P)
Office_Init_conds_forty_mid <- c(N_x = Office_setup_forty_mid$N_x, P = Office_setup_forty_mid$P)
Movie_Init_conds_forty_mid <- c(N_x = Movie_setup_forty_mid$N_x, P = Movie_setup_forty_mid$P)
University_Init_conds_forty_mid <- c(N_x = University_setup_forty_mid$N_x, P = University_setup_forty_mid$P)

Church_Init_conds_fifty_mid <- c(N_x = Church_setup_ten_mid$N_x, P = Church_setup_fifty_mid$P)
Office_Init_conds_fifty_mid <- c(N_x = Office_setup_ten_mid$N_x, P = Office_setup_fifty_mid$P)
Movie_Init_conds_fifty_mid <- c(N_x = Movie_setup_ten_mid$N_x, P = Movie_setup_fifty_mid$P)
University_Init_conds_fifty_mid <- c(N_x = University_setup_ten_mid$N_x, P = University_setup_fifty_mid$P)

Church_Init_conds_sixty_mid <- c(N_x = Church_setup_sixty_mid$N_x, P = Church_setup_sixty_mid$P)
Office_Init_conds_sixty_mid <- c(N_x = Office_setup_sixty_mid$N_x, P = Office_setup_sixty_mid$P)
Movie_Init_conds_sixty_mid <- c(N_x = Movie_setup_sixty_mid$N_x, P = Movie_setup_sixty_mid$P)
University_Init_conds_sixty_mid <- c(N_x = University_setup_sixty_mid$N_x, P = University_setup_sixty_mid$P)

Church_Init_conds_seventy_mid <- c(N_x = Church_setup_seventy_mid$N_x, P = Church_setup_seventy_mid$P)
Office_Init_conds_seventy_mid <- c(N_x = Office_setup_seventy_mid$N_x, P = Office_setup_seventy_mid$P)
Movie_Init_conds_seventy_mid <- c(N_x = Movie_setup_seventy_mid$N_x, P = Movie_setup_seventy_mid$P)
University_Init_conds_seventy_mid <- c(N_x = University_setup_seventy_mid$N_x, P = University_setup_seventy_mid$P)

Church_Init_conds_eighty_mid <- c(N_x = Church_setup_eighty_mid$N_x, P = Church_setup_eighty_mid$P)
Office_Init_conds_eighty_mid <- c(N_x = Office_setup_eighty_mid$N_x, P = Office_setup_eighty_mid$P)
Movie_Init_conds_eighty_mid <- c(N_x = Movie_setup_eighty_mid$N_x, P = Movie_setup_eighty_mid$P)
University_Init_conds_eighty_mid <- c(N_x = University_setup_eighty_mid$N_x, P = University_setup_eighty_mid$P)

Church_Init_conds_ninety_mid <- c(N_x = Church_setup_ninety_mid$N_x, P = Church_setup_ninety_mid$P)
Office_Init_conds_ninety_mid <- c(N_x = Office_setup_ninety_mid$N_x, P = Office_setup_ninety_mid$P)
Movie_Init_conds_ninety_mid <- c(N_x = Movie_setup_ninety_mid$N_x, P = Movie_setup_ninety_mid$P)
University_Init_conds_ninety_mid <- c(N_x = University_setup_ninety_mid$N_x, P = University_setup_ninety_mid$P)

Church_Init_conds_hundred_mid <- c(N_x = Church_setup_hundred_mid$N_x, P = Church_setup_hundred_mid$P)
Office_Init_conds_hundred_mid <- c(N_x = Office_setup_hundred_mid$N_x, P = Office_setup_hundred_mid$P)
Movie_Init_conds_hundred_mid <- c(N_x = Movie_setup_hundred_mid$N_x, P = Movie_setup_hundred_mid$P)
University_Init_conds_hundred_mid <- c(N_x = University_setup_hundred_mid$N_x, P = University_setup_hundred_mid$P)

#### Initial conditions all prop full options -- end of outbreak
Church_Init_conds_ten_end <- c(N_x = Church_setup_ten_end$N_x, P = Church_setup_ten_end$P)
Office_Init_conds_ten_end <- c(N_x = Office_setup_ten_end$N_x, P = Office_setup_ten_end$P)
Movie_Init_conds_ten_end <- c(N_x = Movie_setup_ten_end$N_x, P = Movie_setup_ten_end$P)
University_Init_conds_ten_end <- c(N_x = University_setup_ten_end$N_x, P = University_setup_ten_end$P)

Church_Init_conds_twenty_end <- c(N_x = Church_setup_twenty_end$N_x, P = Church_setup_twenty_end$P)
Office_Init_conds_twenty_end <- c(N_x = Office_setup_twenty_end$N_x, P = Office_setup_twenty_end$P)
Movie_Init_conds_twenty_end <- c(N_x = Movie_setup_twenty_end$N_x, P = Movie_setup_twenty_end$P)
University_Init_conds_twenty_end <- c(N_x = University_setup_twenty_end$N_x, P = University_setup_twenty_end$P)

Church_Init_conds_thirty_end <- c(N_x = Church_setup_thirty_end$N_x, P = Church_setup_thirty_end$P)
Office_Init_conds_thirty_end <- c(N_x = Office_setup_thirty_end$N_x, P = Office_setup_thirty_end$P)
Movie_Init_conds_thirty_end <- c(N_x = Movie_setup_thirty_end$N_x, P = Movie_setup_thirty_end$P)
University_Init_conds_thirty_end <- c(N_x = University_setup_thirty_end$N_x, P = University_setup_thirty_end$P)

Church_Init_conds_forty_end <- c(N_x = Church_setup_forty_end$N_x, P = Church_setup_forty_end$P)
Office_Init_conds_forty_end <- c(N_x = Office_setup_forty_end$N_x, P = Office_setup_forty_end$P)
Movie_Init_conds_forty_end <- c(N_x = Movie_setup_forty_end$N_x, P = Movie_setup_forty_end$P)
University_Init_conds_forty_end <- c(N_x = University_setup_forty_end$N_x, P = University_setup_forty_end$P)

Church_Init_conds_fifty_end <- c(N_x = Church_setup_ten_end$N_x, P = Church_setup_fifty_end$P)
Office_Init_conds_fifty_end <- c(N_x = Office_setup_ten_end$N_x, P = Office_setup_fifty_end$P)
Movie_Init_conds_fifty_end <- c(N_x = Movie_setup_ten_end$N_x, P = Movie_setup_fifty_end$P)
University_Init_conds_fifty_end <- c(N_x = University_setup_ten_end$N_x, P = University_setup_fifty_end$P)

Church_Init_conds_sixty_end <- c(N_x = Church_setup_sixty_end$N_x, P = Church_setup_sixty_end$P)
Office_Init_conds_sixty_end <- c(N_x = Office_setup_sixty_end$N_x, P = Office_setup_sixty_end$P)
Movie_Init_conds_sixty_end <- c(N_x = Movie_setup_sixty_end$N_x, P = Movie_setup_sixty_end$P)
University_Init_conds_sixty_end <- c(N_x = University_setup_sixty_end$N_x, P = University_setup_sixty_end$P)

Church_Init_conds_seventy_end <- c(N_x = Church_setup_seventy_end$N_x, P = Church_setup_seventy_end$P)
Office_Init_conds_seventy_end <- c(N_x = Office_setup_seventy_end$N_x, P = Office_setup_seventy_end$P)
Movie_Init_conds_seventy_end <- c(N_x = Movie_setup_seventy_end$N_x, P = Movie_setup_seventy_end$P)
University_Init_conds_seventy_end <- c(N_x = University_setup_seventy_end$N_x, P = University_setup_seventy_end$P)

Church_Init_conds_eighty_end <- c(N_x = Church_setup_eighty_end$N_x, P = Church_setup_eighty_end$P)
Office_Init_conds_eighty_end <- c(N_x = Office_setup_eighty_end$N_x, P = Office_setup_eighty_end$P)
Movie_Init_conds_eighty_end <- c(N_x = Movie_setup_eighty_end$N_x, P = Movie_setup_eighty_end$P)
University_Init_conds_eighty_end <- c(N_x = University_setup_eighty_end$N_x, P = University_setup_eighty_end$P)

Church_Init_conds_ninety_end <- c(N_x = Church_setup_ninety_end$N_x, P = Church_setup_ninety_end$P)
Office_Init_conds_ninety_end <- c(N_x = Office_setup_ninety_end$N_x, P = Office_setup_ninety_end$P)
Movie_Init_conds_ninety_end <- c(N_x = Movie_setup_ninety_end$N_x, P = Movie_setup_ninety_end$P)
University_Init_conds_ninety_end <- c(N_x = University_setup_ninety_end$N_x, P = University_setup_ninety_end$P)

Church_Init_conds_hundred_end <- c(N_x = Church_setup_hundred_end$N_x, P = Church_setup_hundred_end$P)
Office_Init_conds_hundred_end <- c(N_x = Office_setup_hundred_end$N_x, P = Office_setup_hundred_end$P)
Movie_Init_conds_hundred_end <- c(N_x = Movie_setup_hundred_end$N_x, P = Movie_setup_hundred_end$P)
University_Init_conds_hundred_end <- c(N_x = University_setup_hundred_end$N_x, P = University_setup_hundred_end$P)





church_graph<- graph_from_adjacency_matrix(church_adjacency_matrix, mode = "undirected")

#take a look
plot(church_graph)
Office_graph<- graph_from_adjacency_matrix(Office_adjacency_matrix, mode = "undirected")
plot(Office_graph)
Movie_graph<- graph_from_adjacency_matrix(Movie_adjacency_matrix, mode = "undirected")
University_graph<- graph_from_adjacency_matrix(University_adjacency_matrix, mode = "undirected")


UTK_colors_fun <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (UTK_colors_fun)
  
  unname(UTK_colors[cols])
}

UTK_colors <- c(turq="#00746F",
                orange ="#E65933",
                dark_blue ="#006C93",
                cream ="#F2F0E8",
                medium_blue ="#517C96",
                burgundy = "#8D2048",
                purple ="#754A7E",
                yellow = "#FED535",
                grey = "#A7A9AC",
                sage = "#579584",
                skyblue= "#B9E1E2",
                brown = "#705550",
                pink = "#EE3E80",
                green = "#ABC178",
                lightblue = "#2197A9",
                yellowgreen = "#EBEA64")

colors_to_use <- c("#44aa99","#ddcc77","#332288","#cc6677","#117733","#88ccee","#aa4499") #but we need 10 colors

my_palette <-c("#fed98e","#fe9929","#d95f0e","#993404")
brightened_variation <- c("#ffeda0", "#feb24c", "#fd8d3c", "#f03b20")
inspired_by_nature <- c("#fdae61", "#f46d43", "#d73027", "#a50026")

V(church_graph)$capacity <- Church_C
V(church_graph)$Room <- seq(1:length(Church_C))
set.seed(154780)
church_network <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes( color = UTK_colors_fun("lightblue"),aes(size = capacity))+
  geom_nodetext_repel(aes(label=Room))+
  theme_blank()+
  scale_size_continuous(
    name = "Room capacity",
    breaks = seq(1, max(V(church_graph)$capacity), by = 50),  # Set custom breaks
    limits = c(1, max(V(church_graph)$capacity)  )            # Define limits for the legend
  )+
  labs(title = "Network representation of a church")+
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust=1)
  )
church_network
#code to save figure
png(filename = "Figures/Church_network.png",units="in", width=6, height=4, res=300)
church_network
dev.off()

V(Office_graph)$capacity <- Office_C
V(Office_graph)$Room <- seq(1:length(Office_C))
set.seed(1544780)
Office_network <- ggplot(Office_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes( color = UTK_colors_fun("lightblue"),aes(size = capacity))+
  geom_nodetext_repel(aes(label=Room))+
  theme_blank()+
  scale_size_continuous(
    name = "Room capacity",
    breaks = seq(1, max(V(Office_graph)$capacity), by = 5),  # Set custom breaks
    limits = c(1, max(V(Office_graph)$capacity)  )            # Define limits for the legend
  )+
  labs(title = "Network representation of an office")+
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust=1)
  )
Office_network

#code to save figure
png(filename = "Figures/Office_network.png",units="in", width=6, height=4, res=300)
Office_network
dev.off()


V(Movie_graph)$capacity <- Movie_C
V(Movie_graph)$Room <- seq(1:length(Movie_C))
set.seed(154780)
Movie_network <- ggplot(Movie_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes( color = UTK_colors_fun("lightblue"),aes(size = capacity))+
  geom_nodetext_repel(aes(label=Room))+
  theme_blank()+
  scale_size_continuous(
    name = "Room capacity",
    breaks = seq(1, max(V(Movie_graph)$capacity), by = 50),  # Set custom breaks
    limits = c(1, max(V(Movie_graph)$capacity)  )            # Define limits for the legend
  )+
  labs(title = "Network representation of a movie theater")+
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust=1)
  )
Movie_network
#code to save figure
png(filename = "Figures/Movie_network.png",units="in", width=6, height=4, res=300)
Movie_network
dev.off()



V(University_graph)$capacity <- University_C
V(University_graph)$Room <- seq(1:length(University_C))
set.seed(154780)
University_network <- ggplot(University_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes( color = UTK_colors_fun("lightblue"),aes(size = capacity))+
  geom_nodetext_repel(aes(label=Room))+
  theme_blank()+
  scale_size_continuous(
    name = "Room capacity",
    breaks = seq(1, max(V(University_graph)$capacity), by = 50),  # Set custom breaks
    limits = c(1, max(V(University_graph)$capacity)  )            # Define limits for the legend
  )+
  labs(title = "Network representation of a university building")+
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust=1)
  )
University_network
#code to save figure
png(filename = "Figures/University_network.png",units="in", width=6, height=4, res=300)
University_network
dev.off()

#### Church analysis ####
##### Read in ALL Church the data #####
file_directory <- "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Data"

# List all files with similar naming, using a pattern
file_list_church_beg <- list.files(path = file_directory, 
                        pattern = "Church_output_Nov24_.*_beg\\.text$", 
                        full.names = TRUE)

# Read all files into a list of data frames
data_list_Church_beg <- lapply(file_list_church_beg, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_Church_beg <- str_extract(basename(file_list_church_beg), "_Nov24_([a-z]+)_beg") %>% 
  str_remove_all("_Nov24_|_beg")
# Combine data frames with the extracted IDs
Church_beg_combined_data <- bind_rows(setNames(data_list_Church_beg, file_ids_Church_beg), .id = "file_source")

# Optionally combine all data frames into one

Church_beg_data_clean <- Church_beg_combined_data%>% pivot_longer(cols = -c(time, file_source),
                                                   names_to = c("State", "Room"),
                                                   names_pattern = "(P|N_x)(\\d+)",
                                                   values_to = "Number")
Church_beg_data_factored <-Church_beg_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines<-Church_beg_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

church_net_with_num_particles_beg<-ggarrange(church_network,Num_part_in_rooms_lines)

png(filename = "Figures/church_net_with_num_particles_beg.png", width = 10, height = 5, units = "in", res = 300)
church_net_with_num_particles_beg
dev.off()



###### Mid outbreak/peak ######

# List all files with similar naming, using a pattern
file_list_church_mid <- list.files(path = file_directory, 
                                   pattern = "Church_output_Nov24_.*_mid\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_Church_mid <- lapply(file_list_church_mid, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_Church_mid <- str_extract(basename(file_list_church_mid), "_Nov24_([a-z]+)_mid") %>% 
  str_remove_all("_Nov24_|_mid")
# Combine data frames with the extracted IDs
combined_data_mid <- bind_rows(setNames(data_list_Church_mid, file_ids_Church_mid), .id = "file_source")

# Optionally combine all data frames into one

Church_mid_data_clean <- combined_data_mid%>% pivot_longer(cols = -c(time, file_source),
                                                   names_to = c("State", "Room"),
                                                   names_pattern = "(P|N_x)(\\d+)",
                                                   values_to = "Number")
Church_mid_data_factored <-Church_mid_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines_mid<-Church_data_clean_mid %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

church_net_with_num_particles_mid<-ggarrange(church_network,Num_part_in_rooms_lines_mid)

png(filename = "Figures/church_net_with_num_particles_mid.png", width = 10, height = 5, units = "in", res = 300)
church_net_with_num_particles_mid
dev.off()


###### end of outbreak ######
# List all files with similar naming, using a pattern
file_list_church_end <- list.files(path = file_directory, 
                                   pattern = "Church_output_Nov24_.*_end\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_Church_end <- lapply(file_list_church_end, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_Church_end <- str_extract(basename(file_list_church_end), "_Nov24_([a-z]+)_end") %>% 
  str_remove_all("_Nov24_|_end")
# Combine data frames with the extracted IDs
Church_end_combined_data <- bind_rows(setNames(data_list_Church_end, file_ids_Church_end), .id = "file_source")

# Optionally combine all data frames into one

Church_end_data_clean <- Church_end_combined_data%>% pivot_longer(cols = -c(time, file_source),
                                                                  names_to = c("State", "Room"),
                                                                  names_pattern = "(P|N_x)(\\d+)",
                                                                  values_to = "Number")
Church_end_data_factored <-Church_end_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines_end<-Church_end_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

church_net_with_num_particles_end<-ggarrange(church_network,Num_part_in_rooms_lines_end)

png(filename = "Figures/church_net_with_num_particles_end.png", width = 10, height = 5, units = "in", res = 300)
church_net_with_num_particles_end
dev.off()


##### making first bar plots - Church - beginning of outbreak at different capacities #####

Church_beg_equil_data<-Church_beg_data_factored %>% filter(time == max(time))

Church_beg_bar_plot <-Church_beg_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col()+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


Church_mid_equil_data<-Church_mid_data_factored %>% filter(time == max(time))

Church_mid_bar_plot <- Church_mid_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col()+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


Church_end_equil_data<-Church_end_data_factored %>% filter(time == max(time))

Church_end_bar_plot <- Church_end_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col(position = position_dodge(width = 0.8))+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()

  #ggsave("plot.png", width = 10, height = 6)
ggarrange(Church_beg_bar_plot, Church_mid_bar_plot, Church_end_bar_plot, nrow = 3)

ggsave("Church_epi_threshold.png", width = 10, height = 10)

##### HVAC threshold bar plots #####
Church_Hvac_thresh_data_beg <- Church_beg_data_clean %>% filter(time == max(time))

temp_C_x_Church <- data.frame(Room = as.character(seq(1:length(Church_C))),C_x = c(Church_C))
Church_Hvac_thresh_data_beg <- left_join(Church_Hvac_thresh_data_beg,temp_C_x_Church, by = "Room")
Church_Hvac_thresh_data_beg<- Church_Hvac_thresh_data_beg%>% 
    mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Church_Hvac_thresh_data_beg <- Church_Hvac_thresh_data_beg %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Church_setup_eighty_beg$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#Church_beg_equil_data<-Church_beg_data_factored %>% filter(time == max(time))

Church_HVAC_bar_plot_beg <-Church_Hvac_thresh_data_beg %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


Church_Hvac_thresh_data_mid <- Church_mid_data_clean %>% filter(time == max(time))

temp_C_x_Church <- data.frame(Room = as.character(seq(1:length(Church_C))),C_x = c(Church_C))
Church_Hvac_thresh_data_mid <- left_join(Church_Hvac_thresh_data_mid,temp_C_x_Church, by = "Room")
Church_Hvac_thresh_data_mid<- Church_Hvac_thresh_data_mid%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Church_Hvac_thresh_data_mid <- Church_Hvac_thresh_data_mid %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Church_setup_eighty_mid$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#Church_mid_equil_data<-Church_mid_data_factored %>% filter(time == max(time))

Church_HVAC_bar_plot_mid <-Church_Hvac_thresh_data_mid %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()



Church_Hvac_thresh_data_end <- Church_end_data_clean %>% filter(time == max(time))

temp_C_x_Church <- data.frame(Room = as.character(seq(1:length(Church_C))),C_x = c(Church_C))
Church_Hvac_thresh_data_end <- left_join(Church_Hvac_thresh_data_end,temp_C_x_Church, by = "Room")
Church_Hvac_thresh_data_end<- Church_Hvac_thresh_data_end%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Church_Hvac_thresh_data_end <- Church_Hvac_thresh_data_end %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Church_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#Church_end_equil_data<-Church_end_data_factored %>% filter(time == max(time))

Church_HVAC_bar_plot_end <- Church_Hvac_thresh_data_end %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


#ggsave("plot.png", width = 10, height = 6)
ggarrange(Church_HVAC_bar_plot_beg, Church_HVAC_bar_plot_mid, Church_HVAC_bar_plot_end, nrow = 3)

ggsave("Church_HVAC_threshold.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 15, height = 12)


##### Extracting ratio value ranges for table #####

# Church_Hvac_thresh_data_beg
# Church_Hvac_thresh_data_mid
# Church_Hvac_thresh_data_end

Church_ratio_data_list <- list(Beginning = Church_Hvac_thresh_data_beg, Middle = Church_Hvac_thresh_data_mid, End = Church_Hvac_thresh_data_end)

# Combine and add a column identifying the source
Church_ratio_data_combined <- bind_rows(Church_ratio_data_list, .id = "Outbreak_point")

Church_ratio_table_data<-Church_ratio_data_combined %>%
  ungroup() %>%
  group_by(Outbreak_point, file_source) %>%
  reframe(ratio_range = range(ratio)) %>%
  summarize(ratio_range = paste0(round(min(ratio_range),digits = 2), " - ", round(max(ratio_range),digits=2)), .by = c(Outbreak_point, file_source)) %>%
  pivot_wider(names_from = file_source, values_from = ratio_range)

Church_ratio_table_data %>%
  kbl(format = "latex", booktabs = TRUE, caption = "Church Ratio Table") %>%
  save_kable("Church_ratio_table.tex") 

##### Infection_risk VS room capacity - Church #####
Infection_risk_VS_room_cap_beg<-Church_Hvac_thresh_data_beg %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Church_beg_Infection_risk_V_Cx_beg <- ggplot(Infection_risk_VS_room_cap_beg, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Infection_risk_VS_room_cap_mid<-Church_Hvac_thresh_data_mid %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Church_mid_Infection_risk_V_Cx <- ggplot(Infection_risk_VS_room_cap_mid, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))


Infection_risk_VS_room_cap_end<-Church_Hvac_thresh_data_end %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Church_end_Infection_risk_V_Cx <- ggplot(Infection_risk_VS_room_cap_end, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(Church_beg_Infection_risk_V_Cx_beg,Church_mid_Infection_risk_V_Cx,Church_end_Infection_risk_V_Cx, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("Church_Part_v_room_cap.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)

##### Infection_risk VS Number of individuals in rooms - Church #####
Infection_risk_VS_Nx_beg<-Church_Hvac_thresh_data_beg %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Church_beg_Infection_risk_V_Nx_beg <- ggplot(Infection_risk_VS_Nx_beg, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Infection_risk_VS_Nx_mid<-Church_Hvac_thresh_data_mid %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Church_mid_Infection_risk_V_Nx <- ggplot(Infection_risk_VS_Nx_mid, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))


Infection_risk_VS_Nx_end<-Church_Hvac_thresh_data_end %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Church_end_Infection_risk_V_Nx <- ggplot(Infection_risk_VS_Nx_end, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(Church_beg_Infection_risk_V_Nx_beg,Church_mid_Infection_risk_V_Nx,Church_end_Infection_risk_V_Nx, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("Church_Part_v_Nx.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)


##### How full rooms are - Church #####

Church_Room_prop_full_beg<-Church_Hvac_thresh_data_beg %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred"))) %>% 
  mutate(Room_prop_full = N_x/C_x)
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
 ggplot(Church_Room_prop_full_beg, aes(x = Room, y = Room_prop_full, fill = file_source)) +
  geom_col()+facet_wrap(~file_source)
  # scale_color_manual(
  #   values = extended_palette, 
  #   name = "Building occupancy\nlevel (%)"
  # ) +  # Apply the extended palette
  # labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  # theme_minimal()+theme(plot.title.position = "panel",
  #                       plot.title = element_text(hjust=0.5))
 Church_Room_prop_full_mid<-Church_Hvac_thresh_data_mid %>% 
   mutate(file_source= str_to_sentence(file_source)) %>% 
   mutate(file_source = factor(file_source, 
                               levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                          "Sixty", "Seventy", "Eighty", "Ninety", "Hundred"))) %>% 
   mutate(Room_prop_full = N_x/C_x)
 #the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots
 
 extended_palette <- colorRampPalette(inspired_by_nature)(10)
 
 # Plot
 ggplot(Church_Room_prop_full_mid, aes(x = Room, y = Room_prop_full, fill = file_source)) +
   geom_col()+facet_wrap(~file_source)
 
 Church_Room_prop_full_end<-Church_Hvac_thresh_data_end %>% 
   mutate(file_source= str_to_sentence(file_source)) %>% 
   mutate(file_source = factor(file_source, 
                               levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                          "Sixty", "Seventy", "Eighty", "Ninety", "Hundred"))) %>% 
   mutate(Room_prop_full = N_x/C_x)
 #the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots
 
 extended_palette <- colorRampPalette(inspired_by_nature)(10)
 
 # Plot
 ggplot(Church_Room_prop_full_end, aes(x = Room, y = Room_prop_full, fill = file_source)) +
   geom_col()+facet_wrap(~file_source)
 ##### Network with prop full and number of individuals - Church #####
 
 Church_Room_prop_full_beg <- Church_Room_prop_full_beg %>% filter(file_source == "Fifty")
 V(church_graph)$capacity <- Church_C
 V(church_graph)$Room <- seq(1:length(Church_C))
 V(church_graph)$Particles <- Church_particles_for_network_beg$Number
 V(church_graph)$Prop_full <-Church_Room_prop_full_beg$Room_prop_full
 V(church_graph)$N_x <-Church_Room_prop_full_beg$N_x
 #head(Church_beg_equil_data)
 set.seed(154780)
 network_w_prop_full_church_beg <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
   geom_edges(color = "grey75")+
   geom_nodes(aes(size = capacity, color = Prop_full))+
   geom_nodetext_repel(aes(label=Room))+
   scale_size_continuous(name = "Capacity") +  # Customize the size legend
   scale_color_gradientn(colors = inspired_by_nature,
                         name = "Proportion full",
                         limits = c(0, 1))+ # Define specific points on the scale
   theme_blank()+
   labs(title = "Beginning of outbreak - How full are rooms")+
   theme(plot.title.position = "plot",
         plot.title = element_text(hjust=1)
   )
 ggsave(filename = "Network_with_Room_prop_full_beg.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")
 
 
 network_w_N_x_church_beg <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
   geom_edges(color = "grey75")+
   geom_nodes(aes(size = capacity, color = N_x))+
   geom_nodetext_repel(aes(label=Room))+
   scale_size_continuous(name = "Capacity") +  # Customize the size legend
   scale_color_gradientn(colors = inspired_by_nature,
                         name = "Proportion full",
                         limits = c(0, max(Church_Room_prop_full_beg$N_x)))+ # Define specific points on the scale
   theme_blank()+
   labs(title = "Beginning of outbreak - How full are rooms")+
   theme(plot.title.position = "plot",
         plot.title = element_text(hjust=1)
   )
 ggsave(filename = "Network_with_N_x_beg.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")
 
 
 Church_Room_prop_full_mid <- Church_Room_prop_full_mid %>% filter(file_source == "Fifty")
 V(church_graph)$capacity <- Church_C
 V(church_graph)$Room <- seq(1:length(Church_C))
 V(church_graph)$Particles <- Church_particles_for_network_mid$Number
 V(church_graph)$Prop_full <-Church_Room_prop_full_mid$Room_prop_full
 V(church_graph)$N_x <-Church_Room_prop_full_mid$N_x
 #head(Church_mid_equil_data)
 set.seed(154780)
 network_w_prop_full_church_mid <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
   geom_edges(color = "grey75")+
   geom_nodes(aes(size = capacity, color = Prop_full))+
   geom_nodetext_repel(aes(label=Room))+
   scale_size_continuous(name = "Capacity") +  # Customize the size legend
   scale_color_gradientn(colors = inspired_by_nature,
                         name = "Proportion full",
                         limits = c(0, 1))+ # Define specific points on the scale
   theme_blank()+
   labs(title = "Middle of outbreak - How full are rooms")+
   theme(plot.title.position = "plot",
         plot.title = element_text(hjust=1)
   )
 ggsave(filename = "Network_with_Room_prop_full_mid.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")
 set.seed(154780)
 network_w_N_x_church_mid <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
   geom_edges(color = "grey75")+
   geom_nodes(aes(size = capacity, color = N_x))+
   geom_nodetext_repel(aes(label=Room))+
   scale_size_continuous(name = "Capacity") +  # Customize the size legend
   scale_color_gradientn(colors = inspired_by_nature,
                         name = "Number of\nindividuals",
                         limits = c(0,  max(Church_Room_prop_full_mid$N_x)))+ # Define specific points on the scale
   theme_blank()+
   labs(title = "Middle of outbreak - How full are rooms")+
   theme(plot.title.position = "plot",
         plot.title = element_text(hjust=1)
   )
 ggsave(filename = "Network_with_N_x_mid.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")
 
 
 Church_Room_prop_full_end <- Church_Room_prop_full_end %>% filter(file_source == "Fifty")
 V(church_graph)$capacity <- Church_C
 V(church_graph)$Room <- seq(1:length(Church_C))
 V(church_graph)$Particles <- Church_particles_for_network_end$Number
 V(church_graph)$Prop_full <-Church_Room_prop_full_end$Room_prop_full
 V(church_graph)$N_x <-Church_Room_prop_full_end$N_x
 #head(Church_end_equil_data)
 set.seed(154780)
 network_w_prop_full_church_end <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
   geom_edges(color = "grey75")+
   geom_nodes(aes(size = capacity, color = Prop_full))+
   geom_nodetext_repel(aes(label=Room))+
   scale_size_continuous(name = "Capacity") +  # Customize the size legend
   scale_color_gradientn(colors = inspired_by_nature,
                         name = "Proportion full",
                         limits = c(0, 1))+ # Define specific points on the scale
   theme_blank()+
   labs(title = "End outbreak - How full are rooms")+
   theme(plot.title.position = "plot",
         plot.title = element_text(hjust=1)
   )
 ggsave(filename = "Network_with_Room_prop_full_beg.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")
 
 set.seed(154780)
 network_w_N_x_church_end <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
   geom_edges(color = "grey75")+
   geom_nodes(aes(size = capacity, color = N_x))+
   geom_nodetext_repel(aes(label=Room))+
   scale_size_continuous(name = "Capacity") +  # Customize the size legend
   scale_color_gradientn(colors = inspired_by_nature,
                         name = "Number of\nindividuals",
                         limits = c(0,  max(Church_Room_prop_full_end$N_x)))+ # Define specific points on the scale
   theme_blank()+
   labs(title = "End of outbreak - How full are rooms")+
   theme(plot.title.position = "plot",
         plot.title = element_text(hjust=1)
   )
 ggsave(filename = "Network_with_N_x_end.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")
 
 ################# Note to self for tomorrow ######################
 # So based on the bar plots, it does look like with increased occupancy, how full rooms are, increases
 # The small rooms still fill up really fast though (these only have a capacity of like 1 individual)
 # Looking at the raw number of people, most of the building is still within those large rooms. It might be good to sanity check with a different building
 # The church is quite small and it is easy to get to small rooms
 # does the way I set my initial conditions matter? It shouldn't, that is the whole point of being at steady state.
 # but it does mean that rooms with one individual are probably already filled from the beginning and you won't see much change. 
 # but it should be at equilibrium so it shouldn't make a difference anyways.
 # I do think it would be worth while to look at other buildings that have more large rooms and different sized rooms 
 # looking at the different outbreak times/ proportion of infected doesn't make much sense here because that's not what were interested in
 # Were interested in differing numbers of individuals, so it may actually be better to look at different occupancy levels rather than proportions of infectious individuals in the building.
 
 ggarrange(network_w_prop_full_church_beg,network_w_prop_full_church_mid,network_w_prop_full_church_end, ncol = 3,common.legend = TRUE,legend = "right")
 ggsave(filename = "Network_with_Room_prop_full_All.png",
        path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures",
        width = 9, height = 5,dpi = 300)
 

##### Particles versus degree #####

Church_temp_adj_matrix <- church_adjacency_matrix
diag(Church_temp_adj_matrix) <- 0
Church_temp_graph <- graph_from_adjacency_matrix(Church_temp_adj_matrix,mode = "undirected")

Church_degree_data_beg <- Church_beg_data_clean %>% filter(time == max(time))

temp_C_x_Church <- data.frame(Room = as.character(seq(1:length(Church_C))),C_x = c(Church_C),Degree = degree(temp_graph))
Church_degree_data_beg <- left_join(Church_degree_data_beg,temp_C_x_Church, by = "Room")
Church_degree_data_beg<- Church_degree_data_beg%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Church_degree_data_beg <- Church_degree_data_beg %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Church_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

Church_beg_Part_v_degree <-ggplot(Church_degree_data_beg,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Church_degree_data_mid <- Church_mid_data_clean %>% filter(time == max(time))

temp_C_x_Church <- data.frame(Room = as.character(seq(1:length(Church_C))),C_x = c(Church_C),Degree = degree(temp_graph))
Church_degree_data_mid <- left_join(Church_degree_data_mid,temp_C_x_Church, by = "Room")
Church_degree_data_mid<- Church_degree_data_mid%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Church_degree_data_mid <- Church_degree_data_mid %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Church_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

Church_mid_Part_v_degree <- ggplot(Church_degree_data_mid,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Church_degree_data_end <- Church_end_data_clean %>% filter(time == max(time))

temp_C_x_Church <- data.frame(Room = as.character(seq(1:length(Church_C))),C_x = c(Church_C),Degree = degree(temp_graph))
Church_degree_data_end <- left_join(Church_degree_data_end,temp_C_x_Church, by = "Room")
Church_degree_data_end<- Church_degree_data_end%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Church_degree_data_end <- Church_degree_data_end %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Church_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

Church_end_Part_v_degree <-ggplot(Church_degree_data_end,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(Church_beg_Part_v_degree,Church_mid_Part_v_degree,Church_end_Part_v_degree, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("Church_Part_v_degree.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)



##### Network with num particles - Church #####

Church_particles_for_network_beg <- Church_beg_equil_data %>% filter(file_source == "fifty")
V(church_graph)$capacity <- Church_C
V(church_graph)$Room <- seq(1:length(Church_C))
V(church_graph)$Particles <- Church_particles_for_network_beg$Number
#head(Church_beg_equil_data)
  set.seed(154780)
network_w_particles_church_beg <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 60000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 10000, 20000, 30000, 40000,50000) )+ # Define specific points on the scale
  theme_blank()+
  labs(title = "Beginning of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_beg.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")


Church_particles_for_network_mid <- Church_mid_equil_data %>% filter(file_source == "fifty")
V(church_graph)$capacity <- Church_C
V(church_graph)$Room <- seq(1:length(Church_C))
V(church_graph)$Particles <- Church_particles_for_network_mid$Number
#head(Church_beg_equil_data)
set.seed(154780)
network_w_particles_church_mid<- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 60000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 10000, 20000, 30000, 40000,50000) )  +
  theme_blank()+
  labs(title = "Middle/peak of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_mid.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")

Church_particles_for_network_end <- Church_end_equil_data %>% filter(file_source == "fifty")
V(church_graph)$capacity <- Church_C
V(church_graph)$Room <- seq(1:length(Church_C))
V(church_graph)$Particles <- Church_particles_for_network_end$Number
#head(Church_beg_equil_data)
set.seed(154780)
network_w_particles_church_end <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 60000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 10000, 20000, 30000, 40000,50000) )  +
  theme_blank()+
  labs(title = "End of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_end.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")


ggarrange(network_w_particles_church_beg,network_w_particles_church_mid,network_w_particles_church_end, ncol = 3,common.legend = TRUE,legend = "right")
ggsave(filename = "Network_with_num_particles_All.png",
       path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures",
       width = 9, height = 5,dpi = 300)
#### Office analysis ####
##### Read in ALL the Office data #####

##### Beginning of outbreak data ######
file_directory <- "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Data"

# List all files with similar naming, using a pattern
file_list_Office_beg <- list.files(path = file_directory, 
                                   pattern = "Office_output_Nov24_.*_beg\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_Office_beg <- lapply(file_list_Office_beg, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_Office_beg <- str_extract(basename(file_list_Office_beg), "_Nov24_([a-z]+)_beg") %>% 
  str_remove_all("_Nov24_|_beg")
# Combine data frames with the extracted IDs
Office_beg_combined_data <- bind_rows(setNames(data_list_Office_beg, file_ids_Office_beg), .id = "file_source")

# Optionally combine all data frames into one

Office_beg_data_clean <- Office_beg_combined_data%>% pivot_longer(cols = -c(time, file_source),
                                                                  names_to = c("State", "Room"),
                                                                  names_pattern = "(P|N_x)(\\d+)",
                                                                  values_to = "Number")
Office_beg_data_factored <-Office_beg_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines<-Office_beg_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

Office_net_with_num_particles_beg<-ggarrange(Office_network,Num_part_in_rooms_lines)

png(filename = "Figures/Office_net_with_num_particles_beg.png", width = 10, height = 5, units = "in", res = 300)
Office_net_with_num_particles_beg
dev.off()


##### Mid outbreak/peak data #####

# List all files with similar naming, using a pattern
file_list_Office_mid <- list.files(path = file_directory, 
                                   pattern = "Office_output_Nov24_.*_mid\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_Office_mid <- lapply(file_list_Office_mid, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_Office_mid <- str_extract(basename(file_list_Office_mid), "_Nov24_([a-z]+)_mid") %>% 
  str_remove_all("_Nov24_|_mid")
# Combine data frames with the extracted IDs
Office_mid_combined <- bind_rows(setNames(data_list_Office_mid, file_ids_Office_mid), .id = "file_source")

# Optionally combine all data frames into one

Office_mid_data_clean <- Office_mid_combined%>% pivot_longer(cols = -c(time, file_source),
                                                           names_to = c("State", "Room"),
                                                           names_pattern = "(P|N_x)(\\d+)",
                                                           values_to = "Number")
Office_mid_data_factored <-Office_mid_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines_mid<-Office_mid_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

Office_net_with_num_particles_mid<-ggarrange(Office_network,Num_part_in_rooms_lines_mid)

png(filename = "Figures/Office_net_with_num_particles_mid.png", width = 10, height = 5, units = "in", res = 300)
Office_net_with_num_particles_mid
dev.off()


###### end of outbreak ######
# List all files with similar naming, using a pattern
file_list_Office_end <- list.files(path = file_directory, 
                                   pattern = "Office_output_Nov24_.*_end\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_Office_end <- lapply(file_list_Office_end, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_Office_end <- str_extract(basename(file_list_Office_end), "_Nov24_([a-z]+)_end") %>% 
  str_remove_all("_Nov24_|_end")
# Combine data frames with the extracted IDs
Office_end_combined_data <- bind_rows(setNames(data_list_Office_end, file_ids_Office_end), .id = "file_source")

# Optionally combine all data frames into one

Office_end_data_clean <- Office_end_combined_data%>% pivot_longer(cols = -c(time, file_source),
                                                                  names_to = c("State", "Room"),
                                                                  names_pattern = "(P|N_x)(\\d+)",
                                                                  values_to = "Number")
Office_end_data_factored <-Office_end_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines_end<-Office_end_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

Office_net_with_num_particles_end<-ggarrange(Office_network,Num_part_in_rooms_lines_end)

png(filename = "Figures/Office_net_with_num_particles_end.png", width = 10, height = 5, units = "in", res = 300)
Office_net_with_num_particles_end
dev.off()


##### making first bar plots - Office - beginning of outbreak at different capacities #####

Office_beg_equil_data<-Office_beg_data_factored %>% filter(time == max(time))

Office_beg_bar_plot <-Office_beg_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col()+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


Office_mid_equil_data<-Office_mid_data_factored %>% filter(time == max(time))

Office_mid_bar_plot <- Office_mid_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col()+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


Office_end_equil_data<-Office_end_data_factored %>% filter(time == max(time))

Office_end_bar_plot <- Office_end_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col(position = position_dodge(width = 0.8))+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()

#ggsave("plot.png", width = 10, height = 6)
ggarrange(Office_beg_bar_plot, Office_mid_bar_plot, Office_end_bar_plot, nrow = 3)

ggsave("Office_epi_threshold.png", width = 10, height = 10)


##### HVAC threshold bar plots - Office #####
Office_Hvac_thresh_data_beg <- Office_beg_data_clean %>% filter(time == max(time))

temp_C_x_Office <- data.frame(Room = as.character(seq(1:length(Office_C))),C_x = c(Office_C))
Office_Hvac_thresh_data_beg <- left_join(Office_Hvac_thresh_data_beg,temp_C_x_Office, by = "Room")
Office_Hvac_thresh_data_beg<- Office_Hvac_thresh_data_beg%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Office_Hvac_thresh_data_beg <- Office_Hvac_thresh_data_beg %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Office_setup_eighty_beg$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#Office_beg_equil_data<-Office_beg_data_factored %>% filter(time == max(time))

Office_HVAC_bar_plot_beg <-Office_Hvac_thresh_data_beg %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


Office_Hvac_thresh_data_mid <- Office_mid_data_clean %>% filter(time == max(time))

temp_C_x_Office <- data.frame(Room = as.character(seq(1:length(Office_C))),C_x = c(Office_C))
Office_Hvac_thresh_data_mid <- left_join(Office_Hvac_thresh_data_mid,temp_C_x_Office, by = "Room")
Office_Hvac_thresh_data_mid<- Office_Hvac_thresh_data_mid%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Office_Hvac_thresh_data_mid <- Office_Hvac_thresh_data_mid %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Office_setup_eighty_mid$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#Office_mid_equil_data<-Office_mid_data_factored %>% filter(time == max(time))

Office_HVAC_bar_plot_mid <-Office_Hvac_thresh_data_mid %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()



Office_Hvac_thresh_data_end <- Office_end_data_clean %>% filter(time == max(time))

temp_C_x_Office <- data.frame(Room = as.character(seq(1:length(Office_C))),C_x = c(Office_C))
Office_Hvac_thresh_data_end <- left_join(Office_Hvac_thresh_data_end,temp_C_x_Office, by = "Room")
Office_Hvac_thresh_data_end<- Office_Hvac_thresh_data_end%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Office_Hvac_thresh_data_end <- Office_Hvac_thresh_data_end %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Office_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#Office_end_equil_data<-Office_end_data_factored %>% filter(time == max(time))

Office_HVAC_bar_plot_end <- Office_Hvac_thresh_data_end %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


#ggsave("plot.png", width = 10, height = 6)
ggarrange(Office_HVAC_bar_plot_beg, Office_HVAC_bar_plot_mid, Office_HVAC_bar_plot_end, nrow = 3)

ggsave("Office_HVAC_threshold.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 15, height = 12)
##### Extracting ratio value ranges for table #####

# Office_Hvac_thresh_data_beg
# Office_Hvac_thresh_data_mid
# Office_Hvac_thresh_data_end

Office_ratio_data_list <- list(Beginning = Office_Hvac_thresh_data_beg, Middle = Office_Hvac_thresh_data_mid, End = Office_Hvac_thresh_data_end)

# Combine and add a column identifying the source
Office_ratio_data_combined <- bind_rows(Office_ratio_data_list, .id = "Outbreak_point")

Office_ratio_table_data<-Office_ratio_data_combined %>%
  ungroup() %>%
  group_by(Outbreak_point, file_source) %>%
  reframe(ratio_range = range(ratio)) %>%
  summarize(ratio_range = paste0(round(min(ratio_range),digits = 2), " - ", round(max(ratio_range),digits=2)), .by = c(Outbreak_point, file_source)) %>%
  pivot_wider(names_from = file_source, values_from = ratio_range)

Office_ratio_table_data %>%
  kbl(format = "latex", booktabs = TRUE, caption = "Office Ratio Table") %>%
  save_kable("Office_ratio_table.tex") 

##### Infection_risk VS room capacity - Office #####
Infection_risk_VS_room_cap_beg<-Office_Hvac_thresh_data_beg %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Office_beg_Infection_risk_V_Cx_beg <- ggplot(Infection_risk_VS_room_cap_beg, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Infection_risk_VS_room_cap_mid<-Office_Hvac_thresh_data_mid %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Office_mid_Infection_risk_V_Cx <- ggplot(Infection_risk_VS_room_cap_mid, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))


Infection_risk_VS_room_cap_end<-Office_Hvac_thresh_data_end %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Office_end_Infection_risk_V_Cx <- ggplot(Infection_risk_VS_room_cap_end, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(Office_beg_Infection_risk_V_Cx_beg,Office_mid_Infection_risk_V_Cx,Office_end_Infection_risk_V_Cx, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("Office_Part_v_room_cap.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)


##### Infection_risk VS Number of individuals in rooms - Office #####
Infection_risk_VS_Nx_beg<-Office_Hvac_thresh_data_beg %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Office_beg_Infection_risk_V_Nx_beg <- ggplot(Infection_risk_VS_Nx_beg, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Infection_risk_VS_Nx_mid<-Office_Hvac_thresh_data_mid %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Office_mid_Infection_risk_V_Nx <- ggplot(Infection_risk_VS_Nx_mid, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))


Infection_risk_VS_Nx_end<-Office_Hvac_thresh_data_end %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Office_end_Infection_risk_V_Nx <- ggplot(Infection_risk_VS_Nx_end, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(Office_beg_Infection_risk_V_Nx_beg,Office_mid_Infection_risk_V_Nx,Office_end_Infection_risk_V_Nx, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("Office_Part_v_Nx.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)


##### Particles versus degree - Office #####

Office_temp_adj_matrix <- Office_adjacency_matrix
diag(Office_temp_adj_matrix) <- 0
Office_temp_graph <- graph_from_adjacency_matrix(Office_temp_adj_matrix,mode = "undirected")

Office_degree_data_beg <- Office_beg_data_clean %>% filter(time == max(time))

temp_C_x_Office <- data.frame(Room = as.character(seq(1:length(Office_C))),C_x = c(Office_C),Degree = degree(Office_temp_graph))
Office_degree_data_beg <- left_join(Office_degree_data_beg,temp_C_x_Office, by = "Room")
Office_degree_data_beg<- Office_degree_data_beg%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Office_degree_data_beg <- Office_degree_data_beg %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Office_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

Office_beg_Part_v_degree <-ggplot(Office_degree_data_beg,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Office_degree_data_mid <- Office_mid_data_clean %>% filter(time == max(time))

temp_C_x_Office <- data.frame(Room = as.character(seq(1:length(Office_C))),C_x = c(Office_C),Degree = degree(Office_temp_graph))
Office_degree_data_mid <- left_join(Office_degree_data_mid,temp_C_x_Office, by = "Room")
Office_degree_data_mid<- Office_degree_data_mid%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Office_degree_data_mid <- Office_degree_data_mid %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Office_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

Office_mid_Part_v_degree <- ggplot(Office_degree_data_mid,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Office_degree_data_end <- Office_end_data_clean %>% filter(time == max(time))

temp_C_x_Office <- data.frame(Room = as.character(seq(1:length(Office_C))),C_x = c(Office_C),Degree = degree(Office_temp_graph))
Office_degree_data_end <- left_join(Office_degree_data_end,temp_C_x_Office, by = "Room")
Office_degree_data_end<- Office_degree_data_end%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Office_degree_data_end <- Office_degree_data_end %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Office_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

Office_end_Part_v_degree <-ggplot(Office_degree_data_end,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(Office_beg_Part_v_degree,Office_mid_Part_v_degree,Office_end_Part_v_degree, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("Office_Part_v_degree.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)



##### Network with Num particles - Office #####
Office_particles_for_network_beg <- Office_beg_equil_data %>% filter(file_source == "fifty")
V(Office_graph)$capacity <- Office_C
V(Office_graph)$Room <- seq(1:length(Office_C))
V(Office_graph)$Particles <- Office_particles_for_network_beg$Number
#head(Office_beg_equil_data)
set.seed(154780)
network_w_particles_Office_beg <- ggplot(Office_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 20000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 4000, 8000, 12000, 16000,20000) )+ # Define specific points on the scale
  theme_blank()+
  labs(title = "Beginning of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_beg_Office.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")


Office_particles_for_network_mid <-Office_mid_equil_data %>% filter(file_source == "fifty")
V(Office_graph)$capacity <- Office_C
V(Office_graph)$Room <- seq(1:length(Office_C))
V(Office_graph)$Particles <- Office_particles_for_network_mid$Number
#head(Office_beg_equil_data)
set.seed(154780)
network_w_particles_Office_mid<- ggplot(Office_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 20000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 4000, 8000, 12000, 16000,20000)  )  +
  theme_blank()+
  labs(title = "Middle/peak of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_mid_Office.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")

Office_particles_for_network_end <- Office_end_equil_data %>% filter(file_source == "fifty")
V(Office_graph)$capacity <- Office_C
V(Office_graph)$Room <- seq(1:length(Office_C))
V(Office_graph)$Particles <- Office_particles_for_network_end$Number
#head(Office_beg_equil_data)
set.seed(154780)
network_w_particles_Office_end <- ggplot(Office_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 20000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 4000, 8000, 12000, 16000,20000)  )  +
  theme_blank()+
  labs(title = "End of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_end_Office.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")


ggarrange(network_w_particles_Office_beg,network_w_particles_Office_mid,network_w_particles_Office_end, ncol = 3,common.legend = TRUE,legend = "right")
ggsave(filename = "Network_with_num_particles_All_Office.png",
       path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures",
       width = 9, height = 5,dpi = 300)
#### Movie Analysis ####
##### Read in ALL the Movie data #####
file_directory <- "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Data"

# List all files with similar naming, using a pattern
file_list_Movie_beg <- list.files(path = file_directory, 
                                   pattern = "Movie_output_Nov24_.*_beg\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_Movie_beg <- lapply(file_list_Movie_beg, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_Movie_beg <- str_extract(basename(file_list_Movie_beg), "_Nov24_([a-z]+)_beg") %>% 
  str_remove_all("_Nov24_|_beg")
# Combine data frames with the extracted IDs
Movie_beg_combined_data <- bind_rows(setNames(data_list_Movie_beg, file_ids_Movie_beg), .id = "file_source")

# Optionally combine all data frames into one

Movie_beg_data_clean <- Movie_beg_combined_data%>% pivot_longer(cols = -c(time, file_source),
                                                                  names_to = c("State", "Room"),
                                                                  names_pattern = "(P|N_x)(\\d+)",
                                                                  values_to = "Number")
Movie_beg_data_factored <-Movie_beg_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines<-Movie_beg_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

Movie_net_with_num_particles_beg<-ggarrange(Movie_network,Num_part_in_rooms_lines)

png(filename = "Figures/Movie_net_with_num_particles_beg.png", width = 10, height = 5, units = "in", res = 300)
Movie_net_with_num_particles_beg
dev.off()


##### Mid outbreak/peak #####

# List all files with similar naming, using a pattern
file_list_Movie_mid <- list.files(path = file_directory, 
                                   pattern = "Movie_output_Nov24_.*_mid\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_Movie_mid <- lapply(file_list_Movie_mid, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_Movie_mid <- str_extract(basename(file_list_Movie_mid), "_Nov24_([a-z]+)_mid") %>% 
  str_remove_all("_Nov24_|_mid")
# Combine data frames with the extracted IDs
Movie_mid_combined <- bind_rows(setNames(data_list_Movie_mid, file_ids_Movie_mid), .id = "file_source")

# Optionally combine all data frames into one

Movie_mid_data_clean <- Movie_mid_combined%>% pivot_longer(cols = -c(time, file_source),
                                                             names_to = c("State", "Room"),
                                                             names_pattern = "(P|N_x)(\\d+)",
                                                             values_to = "Number")
Movie_mid_data_factored <-Movie_mid_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines_mid<-Movie_mid_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

Movie_net_with_num_particles_mid<-ggarrange(Movie_network,Num_part_in_rooms_lines_mid)

png(filename = "Figures/Movie_net_with_num_particles_mid.png", width = 10, height = 5, units = "in", res = 300)
Movie_net_with_num_particles_mid
dev.off()


##### end of outbreak #####
# List all files with similar naming, using a pattern
file_list_Movie_end <- list.files(path = file_directory, 
                                   pattern = "Movie_output_Nov24_.*_end\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_Movie_end <- lapply(file_list_Movie_end, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_Movie_end <- str_extract(basename(file_list_Movie_end), "_Nov24_([a-z]+)_end") %>% 
  str_remove_all("_Nov24_|_end")
# Combine data frames with the extracted IDs
Movie_end_combined_data <- bind_rows(setNames(data_list_Movie_end, file_ids_Movie_end), .id = "file_source")

# Optionally combine all data frames into one

Movie_end_data_clean <- Movie_end_combined_data%>% pivot_longer(cols = -c(time, file_source),
                                                                  names_to = c("State", "Room"),
                                                                  names_pattern = "(P|N_x)(\\d+)",
                                                                  values_to = "Number")
Movie_end_data_factored <-Movie_end_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines_end<-Movie_end_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

Movie_net_with_num_particles_end<-ggarrange(Movie_network,Num_part_in_rooms_lines_end)

png(filename = "Figures/Movie_net_with_num_particles_end.png", width = 10, height = 5, units = "in", res = 300)
Movie_net_with_num_particles_end
dev.off()


##### making first bar plots - Movie - beginning of outbreak at different capacities #####

Movie_beg_equil_data<-Movie_beg_data_factored %>% filter(time == max(time))

Movie_beg_bar_plot <-Movie_beg_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col()+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


Movie_mid_equil_data<-Movie_mid_data_factored %>% filter(time == max(time))

Movie_mid_bar_plot <- Movie_mid_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col()+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


Movie_end_equil_data<-Movie_end_data_factored %>% filter(time == max(time))

Movie_end_bar_plot <- Movie_end_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col(position = position_dodge(width = 0.8))+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()

#ggsave("plot.png", width = 10, height = 6)
ggarrange(Movie_beg_bar_plot, Movie_mid_bar_plot, Movie_end_bar_plot, nrow = 3)

ggsave("Movie_epi_threshold.png", width = 10, height = 10)


##### HVAC threshold bar plots - Movie #####
Movie_Hvac_thresh_data_beg <- Movie_beg_data_clean %>% filter(time == max(time))

temp_C_x_Movie <- data.frame(Room = as.character(seq(1:length(Movie_C))),C_x = c(Movie_C))
Movie_Hvac_thresh_data_beg <- left_join(Movie_Hvac_thresh_data_beg,temp_C_x_Movie, by = "Room")
Movie_Hvac_thresh_data_beg<- Movie_Hvac_thresh_data_beg%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Movie_Hvac_thresh_data_beg <- Movie_Hvac_thresh_data_beg %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Movie_setup_eighty_beg$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#Movie_beg_equil_data<-Movie_beg_data_factored %>% filter(time == max(time))

Movie_HVAC_bar_plot_beg <-Movie_Hvac_thresh_data_beg %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


Movie_Hvac_thresh_data_mid <- Movie_mid_data_clean %>% filter(time == max(time))

temp_C_x_Movie <- data.frame(Room = as.character(seq(1:length(Movie_C))),C_x = c(Movie_C))
Movie_Hvac_thresh_data_mid <- left_join(Movie_Hvac_thresh_data_mid,temp_C_x_Movie, by = "Room")
Movie_Hvac_thresh_data_mid<- Movie_Hvac_thresh_data_mid%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Movie_Hvac_thresh_data_mid <- Movie_Hvac_thresh_data_mid %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Movie_setup_eighty_mid$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#Movie_mid_equil_data<-Movie_mid_data_factored %>% filter(time == max(time))

Movie_HVAC_bar_plot_mid <-Movie_Hvac_thresh_data_mid %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()



Movie_Hvac_thresh_data_end <- Movie_end_data_clean %>% filter(time == max(time))

temp_C_x_Movie <- data.frame(Room = as.character(seq(1:length(Movie_C))),C_x = c(Movie_C))
Movie_Hvac_thresh_data_end <- left_join(Movie_Hvac_thresh_data_end,temp_C_x_Movie, by = "Room")
Movie_Hvac_thresh_data_end<- Movie_Hvac_thresh_data_end%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Movie_Hvac_thresh_data_end <- Movie_Hvac_thresh_data_end %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Movie_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#Movie_end_equil_data<-Movie_end_data_factored %>% filter(time == max(time))

Movie_HVAC_bar_plot_end <- Movie_Hvac_thresh_data_end %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


#ggsave("plot.png", width = 10, height = 6)
ggarrange(Movie_HVAC_bar_plot_beg, Movie_HVAC_bar_plot_mid, Movie_HVAC_bar_plot_end, nrow = 3)

ggsave("Movie_HVAC_threshold.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 15, height = 12)



##### Infection_risk VS room capacity - Movie #####
Infection_risk_VS_room_cap_beg<-Movie_Hvac_thresh_data_beg %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Movie_beg_Infection_risk_V_Cx_beg <- ggplot(Infection_risk_VS_room_cap_beg, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Infection_risk_VS_room_cap_mid<-Movie_Hvac_thresh_data_mid %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Movie_mid_Infection_risk_V_Cx <- ggplot(Infection_risk_VS_room_cap_mid, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))


Infection_risk_VS_room_cap_end<-Movie_Hvac_thresh_data_end %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Movie_end_Infection_risk_V_Cx <- ggplot(Infection_risk_VS_room_cap_end, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(Movie_beg_Infection_risk_V_Cx_beg,Movie_mid_Infection_risk_V_Cx,Movie_end_Infection_risk_V_Cx, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("Movie_Part_v_room_cap.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)

##### Infection_risk VS Number of individuals in rooms - Movie #####
Infection_risk_VS_Nx_beg<-Movie_Hvac_thresh_data_beg %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Movie_beg_Infection_risk_V_Nx_beg <- ggplot(Infection_risk_VS_Nx_beg, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Infection_risk_VS_Nx_mid<-Movie_Hvac_thresh_data_mid %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Movie_mid_Infection_risk_V_Nx <- ggplot(Infection_risk_VS_Nx_mid, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))


Infection_risk_VS_Nx_end<-Movie_Hvac_thresh_data_end %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
Movie_end_Infection_risk_V_Nx <- ggplot(Infection_risk_VS_Nx_end, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(Movie_beg_Infection_risk_V_Nx_beg,Movie_mid_Infection_risk_V_Nx,Movie_end_Infection_risk_V_Nx, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("Movie_Part_v_Nx.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)


##### Particles versus degree - Movie #####

Movie_temp_adj_matrix <- Movie_adjacency_matrix
diag(Movie_temp_adj_matrix) <- 0
Movie_temp_graph <- graph_from_adjacency_matrix(Movie_temp_adj_matrix,mode = "undirected")

Movie_degree_data_beg <- Movie_beg_data_clean %>% filter(time == max(time))

temp_C_x_Movie <- data.frame(Room = as.character(seq(1:length(Movie_C))),C_x = c(Movie_C),Degree = degree(Movie_temp_graph))
Movie_degree_data_beg <- left_join(Movie_degree_data_beg,temp_C_x_Movie, by = "Room")
Movie_degree_data_beg<- Movie_degree_data_beg%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Movie_degree_data_beg <- Movie_degree_data_beg %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Movie_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

Movie_beg_Part_v_degree <-ggplot(Movie_degree_data_beg,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Movie_degree_data_mid <- Movie_mid_data_clean %>% filter(time == max(time))

temp_C_x_Movie <- data.frame(Room = as.character(seq(1:length(Movie_C))),C_x = c(Movie_C),Degree = degree(Movie_temp_graph))
Movie_degree_data_mid <- left_join(Movie_degree_data_mid,temp_C_x_Movie, by = "Room")
Movie_degree_data_mid<- Movie_degree_data_mid%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Movie_degree_data_mid <- Movie_degree_data_mid %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Movie_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

Movie_mid_Part_v_degree <- ggplot(Movie_degree_data_mid,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Movie_degree_data_end <- Movie_end_data_clean %>% filter(time == max(time))

temp_C_x_Movie <- data.frame(Room = as.character(seq(1:length(Movie_C))),C_x = c(Movie_C),Degree = degree(Movie_temp_graph))
Movie_degree_data_end <- left_join(Movie_degree_data_end,temp_C_x_Movie, by = "Room")
Movie_degree_data_end<- Movie_degree_data_end%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
Movie_degree_data_end <- Movie_degree_data_end %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = Movie_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

Movie_end_Part_v_degree <-ggplot(Movie_degree_data_end,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(Movie_beg_Part_v_degree,Movie_mid_Part_v_degree,Movie_end_Part_v_degree, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("Movie_Part_v_degree.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)


##### Network with num particles fig - Movie #####
Movie_particles_for_network_beg <- Movie_beg_equil_data %>% filter(file_source == "fifty")
V(Movie_graph)$capacity <- Movie_C
V(Movie_graph)$Room <- seq(1:length(Movie_C))
V(Movie_graph)$Particles <- Movie_particles_for_network_beg$Number
#head(Movie_beg_equil_data)
set.seed(154780)
network_w_particles_Movie_beg <- ggplot(Movie_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 100000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 40000, 80000, 120000, 160000,200000) )+ # Define specific points on the scale
  theme_blank()+
  labs(title = "Beginning of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_beg_Movie.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")


Movie_particles_for_network_mid <-Movie_mid_equil_data %>% filter(file_source == "fifty")
V(Movie_graph)$capacity <- Movie_C
V(Movie_graph)$Room <- seq(1:length(Movie_C))
V(Movie_graph)$Particles <- Movie_particles_for_network_mid$Number
#head(Movie_beg_equil_data)
set.seed(154780)
network_w_particles_Movie_mid<- ggplot(Movie_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 100000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 40000, 80000, 120000, 160000,200000) )  +
  theme_blank()+
  labs(title = "Middle/peak of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_mid_Movie.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")

Movie_particles_for_network_end <- Movie_end_equil_data %>% filter(file_source == "fifty")
V(Movie_graph)$capacity <- Movie_C
V(Movie_graph)$Room <- seq(1:length(Movie_C))
V(Movie_graph)$Particles <- Movie_particles_for_network_end$Number
#head(Movie_beg_equil_data)
set.seed(154780)
network_w_particles_Movie_end <- ggplot(Movie_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 100000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 40000, 80000, 120000, 160000,200000) )  +
  theme_blank()+
  labs(title = "End of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_end_Movie.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")


ggarrange(network_w_particles_Movie_beg,network_w_particles_Movie_mid,network_w_particles_Movie_end, ncol = 3,common.legend = TRUE,legend = "right")
ggsave(filename = "Network_with_num_particles_All_Movie.png",
       path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures",
       width = 9, height = 5,dpi = 300)
#### University Analysis ####
##### Read in University the University data #####
file_directory <- "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Data"

# List all files with similar naming, using a pattern
file_list_University_beg <- list.files(path = file_directory, 
                                   pattern = "University_output_Nov24_.*_beg\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_University_beg <- lapply(file_list_University_beg, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_University_beg <- str_extract(basename(file_list_University_beg), "_Nov24_([a-z]+)_beg") %>% 
  str_remove_all("_Nov24_|_beg")
# Combine data frames with the extracted IDs
University_beg_combined_data <- bind_rows(setNames(data_list_University_beg, file_ids_University_beg), .id = "file_source")

# Optionally combine all data frames into one

University_beg_data_clean <- University_beg_combined_data%>% pivot_longer(cols = -c(time, file_source),
                                                                  names_to = c("State", "Room"),
                                                                  names_pattern = "(P|N_x)(\\d+)",
                                                                  values_to = "Number")
University_beg_data_factored <-University_beg_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines<-University_beg_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

University_net_with_num_particles_beg<-ggarrange(University_network,Num_part_in_rooms_lines)

png(filename = "Figures/University_net_with_num_particles_beg.png", width = 10, height = 5, units = "in", res = 300)
University_net_with_num_particles_beg
dev.off()


##### Mid outbreak/peak #####

# List all files with similar naming, using a pattern
file_list_University_mid <- list.files(path = file_directory, 
                                   pattern = "University_output_Nov24_.*_mid\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_University_mid <- lapply(file_list_University_mid, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_University_mid <- str_extract(basename(file_list_University_mid), "_Nov24_([a-z]+)_mid") %>% 
  str_remove_all("_Nov24_|_mid")
# Combine data frames with the extracted IDs
University_mid_combined <- bind_rows(setNames(data_list_University_mid, file_ids_University_mid), .id = "file_source")

# Optionally combine all data frames into one

University_mid_data_clean <- University_mid_combined%>% pivot_longer(cols = -c(time, file_source),
                                                             names_to = c("State", "Room"),
                                                             names_pattern = "(P|N_x)(\\d+)",
                                                             values_to = "Number")
University_mid_data_factored <-University_mid_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines_mid<-University_mid_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

University_net_with_num_particles_mid<-ggarrange(University_network,Num_part_in_rooms_lines_mid)

png(filename = "Figures/University_net_with_num_particles_mid.png", width = 10, height = 5, units = "in", res = 300)
University_net_with_num_particles_mid
dev.off()


##### end of outbreak #####
# List all files with similar naming, using a pattern
file_list_University_end <- list.files(path = file_directory, 
                                   pattern = "University_output_Nov24_.*_end\\.text$", 
                                   full.names = TRUE)

# Read all files into a list of data frames
data_list_University_end <- lapply(file_list_University_end, function(file) {
  read.table(file, header = TRUE, sep = " ")  # Adjust `sep` and `header` as needed
})

file_ids_University_end <- str_extract(basename(file_list_University_end), "_Nov24_([a-z]+)_end") %>% 
  str_remove_all("_Nov24_|_end")
# Combine data frames with the extracted IDs
University_end_combined_data <- bind_rows(setNames(data_list_University_end, file_ids_University_end), .id = "file_source")

# Optionally combine all data frames into one

University_end_data_clean <- University_end_combined_data%>% pivot_longer(cols = -c(time, file_source),
                                                                  names_to = c("State", "Room"),
                                                                  names_pattern = "(P|N_x)(\\d+)",
                                                                  values_to = "Number")
University_end_data_factored <-University_end_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 

Num_part_in_rooms_lines_end<-University_end_data_clean %>% 
  filter(State == "P") %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) %>%
  group_by(file_source, Room) %>% 
  mutate(mid_time = time[ceiling(n() / 2)],   # Middle x-value for each Room
         mid_number = Number[ceiling(n() / 2)]) %>%  # Middle y-value for each Room
  ggplot(aes(x = time, y = Number, group = Room, color = as.factor(Room))) +
  geom_line() +
  facet_wrap(~file_source) +  # Ensure that facet_wrap respects the order of file_source
  geom_text(
    data = . %>% 
      group_by(file_source) %>% 
      filter(Number == max(Number)) %>%  # Only keep the row with the max Number for each Room
      ungroup(),  # Ungroup after filtering
    aes(x = mid_time, y = mid_number, label = paste("Room", Room)),
    color = "black",  # Label color
    size = 3          # Label size
  ) +
  theme_classic() +
  labs(x = "Time (hours)", y = "Number", color = "Room") +
  ggtitle("Number of infectious particles")

University_net_with_num_particles_end<-ggarrange(University_network,Num_part_in_rooms_lines_end)

png(filename = "Figures/University_net_with_num_particles_end.png", width = 10, height = 5, units = "in", res = 300)
University_net_with_num_particles_end
dev.off()


##### making first bar plots - University - beginning of outbreak at different capacities #####

University_beg_equil_data<-University_beg_data_factored %>% filter(time == max(time))

University_beg_bar_plot <-University_beg_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col()+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


University_mid_equil_data<-University_mid_data_factored %>% filter(time == max(time))

University_mid_bar_plot <- University_mid_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col()+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


University_end_equil_data<-University_end_data_factored %>% filter(time == max(time))

University_end_bar_plot <- University_end_equil_data %>% 
  filter(State == "P") %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = Number))+
  geom_col(position = position_dodge(width = 0.8))+geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source, ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()

#ggsave("plot.png", width = 10, height = 6)
ggarrange(University_beg_bar_plot, University_mid_bar_plot, University_end_bar_plot, nrow = 3)

ggsave("University_epi_threshold.png", width = 10, height = 10)

University_particles_for_network_beg <- University_beg_equil_data %>% filter(file_source == "fifty")
V(University_graph)$capacity <- University_C
V(University_graph)$Room <- seq(1:length(University_C))
V(University_graph)$Particles <- University_particles_for_network_beg$Number
#head(University_beg_equil_data)
set.seed(154780)
network_w_particles_University_beg <- ggplot(University_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 100000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 40000, 80000, 120000, 160000,200000) )+ # Define specific points on the scale
  theme_blank()+
  labs(title = "Beginning of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_beg_University.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")


University_particles_for_network_mid <-University_mid_equil_data %>% filter(file_source == "fifty")
V(University_graph)$capacity <- University_C
V(University_graph)$Room <- seq(1:length(University_C))
V(University_graph)$Particles <- University_particles_for_network_mid$Number
#head(University_beg_equil_data)
set.seed(154780)
network_w_particles_University_mid<- ggplot(University_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 100000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 40000, 80000, 120000, 160000,200000) )  +
  theme_blank()+
  labs(title = "Middle/peak of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_mid_University.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")

##### HVAC threshold bar plots - University #####
University_Hvac_thresh_data_beg <- University_beg_data_clean %>% filter(time == max(time))

temp_C_x_University <- data.frame(Room = as.character(seq(1:length(University_C))),C_x = c(University_C))
University_Hvac_thresh_data_beg <- left_join(University_Hvac_thresh_data_beg,temp_C_x_University, by = "Room")
University_Hvac_thresh_data_beg<- University_Hvac_thresh_data_beg%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
University_Hvac_thresh_data_beg <- University_Hvac_thresh_data_beg %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = University_setup_eighty_beg$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#University_beg_equil_data<-University_beg_data_factored %>% filter(time == max(time))

University_HVAC_bar_plot_beg <-University_Hvac_thresh_data_beg %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


University_Hvac_thresh_data_mid <- University_mid_data_clean %>% filter(time == max(time))

temp_C_x_University <- data.frame(Room = as.character(seq(1:length(University_C))),C_x = c(University_C))
University_Hvac_thresh_data_mid <- left_join(University_Hvac_thresh_data_mid,temp_C_x_University, by = "Room")
University_Hvac_thresh_data_mid<- University_Hvac_thresh_data_mid%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
University_Hvac_thresh_data_mid <- University_Hvac_thresh_data_mid %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = University_setup_eighty_mid$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#University_mid_equil_data<-University_mid_data_factored %>% filter(time == max(time))

University_HVAC_bar_plot_mid <-University_Hvac_thresh_data_mid %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()



University_Hvac_thresh_data_end <- University_end_data_clean %>% filter(time == max(time))

temp_C_x_University <- data.frame(Room = as.character(seq(1:length(University_C))),C_x = c(University_C))
University_Hvac_thresh_data_end <- left_join(University_Hvac_thresh_data_end,temp_C_x_University, by = "Room")
University_Hvac_thresh_data_end<- University_Hvac_thresh_data_end%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
University_Hvac_thresh_data_end <- University_Hvac_thresh_data_end %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = University_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 


#University_end_equil_data<-University_end_data_factored %>% filter(time == max(time))

University_HVAC_bar_plot_end <- University_Hvac_thresh_data_end %>% 
  filter(file_source == "twenty"| file_source == "forty"| file_source == "sixty" | file_source == "eighty" | file_source == "hundred") %>% 
  group_by(file_source) %>% 
  ggplot(aes(x = Room, y = ratio))+
  geom_col()+#geom_hline(yintercept = 16000,color = "red3",linetype = 2)+
  facet_wrap(~file_source,scales = "free", ncol = 5,labeller = labeller(file_source = function(x) tools::toTitleCase(x)))+
  labs(y = "Number of infectious particles")+
  theme(axis.text.x=element_text(angle = 0))+coord_flip()


#ggsave("plot.png", width = 10, height = 6)
ggarrange(University_HVAC_bar_plot_beg, University_HVAC_bar_plot_mid, University_HVAC_bar_plot_end, nrow = 3)

ggsave("University_HVAC_threshold.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 15, height = 12)


##### Infection_risk VS room capacity - University #####
Infection_risk_VS_room_cap_beg<-University_Hvac_thresh_data_beg %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
University_beg_Infection_risk_V_Cx_beg <- ggplot(Infection_risk_VS_room_cap_beg, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Infection_risk_VS_room_cap_mid<-University_Hvac_thresh_data_mid %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
University_mid_Infection_risk_V_Cx <- ggplot(Infection_risk_VS_room_cap_mid, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))


Infection_risk_VS_room_cap_end<-University_Hvac_thresh_data_end %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
University_end_Infection_risk_V_Cx <- ggplot(Infection_risk_VS_room_cap_end, aes(x = C_x, y = P, color = file_source)) +
  geom_point() +
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Room capacity") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(University_beg_Infection_risk_V_Cx_beg,University_mid_Infection_risk_V_Cx,University_end_Infection_risk_V_Cx, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("University_Part_v_room_cap.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)

##### Infection_risk VS Number of individuals in rooms - University #####
Infection_risk_VS_Nx_beg<-University_Hvac_thresh_data_beg %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
University_beg_Infection_risk_V_Nx_beg <- ggplot(Infection_risk_VS_Nx_beg, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

Infection_risk_VS_Nx_mid<-University_Hvac_thresh_data_mid %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
University_mid_Infection_risk_V_Nx <- ggplot(Infection_risk_VS_Nx_mid, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))


Infection_risk_VS_Nx_end<-University_Hvac_thresh_data_end %>% 
  mutate(file_source= str_to_sentence(file_source)) %>% 
  mutate(file_source = factor(file_source, 
                              levels = c("Ten", "Twenty", "Thirty", "Forty", "Fifty", 
                                         "Sixty", "Seventy", "Eighty", "Ninety", "Hundred")))
#the data manipulation I did for the hvac threshold data is exactly what I need for these plots so just saving it to it's own data frame to use for these plots

extended_palette <- colorRampPalette(inspired_by_nature)(10)

# Plot
University_end_Infection_risk_V_Nx <- ggplot(Infection_risk_VS_Nx_end, aes(x = N_x, y = P, color = file_source)) +
  geom_jitter(alpha = 0.8,width = 5)+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Number of individuals in rooms") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(University_beg_Infection_risk_V_Nx_beg,University_mid_Infection_risk_V_Nx,University_end_Infection_risk_V_Nx, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("University_Part_v_Nx.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)


##### Particles versus degree - University #####

University_temp_adj_matrix <- University_adjacency_matrix
diag(University_temp_adj_matrix) <- 0
University_temp_graph <- graph_from_adjacency_matrix(University_temp_adj_matrix,mode = "undirected")

University_degree_data_beg <- University_beg_data_clean %>% filter(time == max(time))

temp_C_x_University <- data.frame(Room = as.character(seq(1:length(University_C))),C_x = c(University_C),Degree = degree(University_temp_graph))
University_degree_data_beg <- left_join(University_degree_data_beg,temp_C_x_University, by = "Room")
University_degree_data_beg<- University_degree_data_beg%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
University_degree_data_beg <- University_degree_data_beg %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = University_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

University_beg_Part_v_degree <-ggplot(University_degree_data_beg,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Beginning of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

University_degree_data_mid <- University_mid_data_clean %>% filter(time == max(time))

temp_C_x_University <- data.frame(Room = as.character(seq(1:length(University_C))),C_x = c(University_C),Degree = degree(University_temp_graph))
University_degree_data_mid <- left_join(University_degree_data_mid,temp_C_x_University, by = "Room")
University_degree_data_mid<- University_degree_data_mid%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
University_degree_data_mid <- University_degree_data_mid %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = University_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

University_mid_Part_v_degree <- ggplot(University_degree_data_mid,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "Middle of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

University_degree_data_end <- University_end_data_clean %>% filter(time == max(time))

temp_C_x_University <- data.frame(Room = as.character(seq(1:length(University_C))),C_x = c(University_C),Degree = degree(University_temp_graph))
University_degree_data_end <- left_join(University_degree_data_end,temp_C_x_University, by = "Room")
University_degree_data_end<- University_degree_data_end%>% 
  mutate(file_source = factor(file_source, 
                              levels = c("ten", "twenty", "thirty", "forty", "fifty", 
                                         "sixty", "seventy", "eighty", "ninety", "hundred")),  # Custom order for file_source
         Room = as.character(Room),  # Convert Room to a character vector
         Room = factor(Room, levels = as.character(1:length(unique(Room))))) 
University_degree_data_end <- University_degree_data_end %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% 
  group_by(time, Room)%>%mutate(s =(parms$s),a = parms$a,bet_bar = parms$bet_bar,bet_hat= parms$bet_hat, Ib_prop = University_setup_eighty_end$Ib_prop[1],Numerator = s*N_x*Ib_prop , Denominator =(a*N_x*(P/(C_x*bet_bar))*(P/bet_bar)), ratio = Numerator/Denominator ) 

University_end_Part_v_degree <-ggplot(University_degree_data_end,aes(x = Degree, y = P,color=file_source))+
  geom_point()+
  scale_color_manual(
    values = extended_palette, 
    name = "Building occupancy\nlevel (%)"
  ) +  # Apply the extended palette
  labs(title = "End of outbreak",y= "Infectious particles", x = "Node degree") +
  theme_minimal()+theme(plot.title.position = "panel",
                        plot.title = element_text(hjust=0.5))

ggarrange(University_beg_Part_v_degree,University_mid_Part_v_degree,University_end_Part_v_degree, 
          ncol = 3,common.legend = TRUE,legend = "right")
ggsave("University_Part_v_degree.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 9, height = 3)


##### Network with particle numbers #####
University_particles_for_network_end <- University_end_equil_data %>% filter(file_source == "fifty")
V(University_graph)$capacity <- University_C
V(University_graph)$Room <- seq(1:length(University_C))
V(University_graph)$Particles <- University_particles_for_network_end$Number
#head(University_beg_equil_data)
set.seed(154780)
network_w_particles_University_end <- ggplot(University_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = capacity, color = Particles))+
  geom_nodetext_repel(aes(label=Room))+
  scale_size_continuous(name = "Capacity") +  # Customize the size legend
  scale_color_gradientn(colors = inspired_by_nature,
                        name = "Particles",
                        limits = c(0, 100000),  # Set the range for the scale (e.g., 0 to 100)
                        breaks = c(0, 40000, 80000, 120000, 160000,200000) )  +
  theme_blank()+
  labs(title = "End of outbreak")+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust=1)
  )
ggsave(filename = "Network_with_num_particles_end_University.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures")


ggarrange(network_w_particles_University_beg,network_w_particles_University_mid,network_w_particles_University_end, ncol = 3,common.legend = TRUE,legend = "right")
ggsave(filename = "Network_with_num_particles_All_University.png",
       path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures",
       width = 9, height = 5,dpi = 300)


#### Graphs for comparing all buildings ####

#going to use the data frames with degree data, add column for saying what building it is.

Church_degree_data_beg <- Church_degree_data_beg %>% mutate(Building = "Church")
Church_degree_data_mid <- Church_degree_data_mid %>% mutate(Building = "Church")
Church_degree_data_end <- Church_degree_data_end %>% mutate(Building = "Church")

#Office
Office_degree_data_beg <- Office_degree_data_beg %>% mutate(Building = "Office")
Office_degree_data_mid <- Office_degree_data_mid %>% mutate(Building = "Office")
Office_degree_data_end <- Office_degree_data_end %>% mutate(Building = "Office")

#Movie
Movie_degree_data_beg <- Movie_degree_data_beg %>% mutate(Building = "Movie")
Movie_degree_data_mid <- Movie_degree_data_mid %>% mutate(Building = "Movie")
Movie_degree_data_end <- Movie_degree_data_end %>% mutate(Building = "Movie")

#University
University_degree_data_beg <- University_degree_data_beg %>% mutate(Building = "University")
University_degree_data_mid <- University_degree_data_mid %>% mutate(Building = "University")
University_degree_data_end <- University_degree_data_end %>% mutate(Building = "University")


data_list_beg <- list(Church = Church_degree_data_beg, Office= Office_degree_data_beg, Movie = Movie_degree_data_beg, University = University_degree_data_beg)

# Combine and add a column identifying the source
combined_data_beg <- bind_rows(data_list_beg, .id = "Building")

data_list_mid <- list(Church = Church_degree_data_mid, Office= Office_degree_data_mid, Movie = Movie_degree_data_mid, University = University_degree_data_mid)

# Combine and add a column identifying the source
combined_data_mid <- bind_rows(data_list_mid, .id = "Building")

data_list_end <- list(Church = Church_degree_data_end, Office= Office_degree_data_end, Movie = Movie_degree_data_end, University = University_degree_data_end)

# Combine and add a column identifying the source
combined_data_end <- bind_rows(data_list_end, .id = "Building")

combined_data_beg %>% 
  group_by(Building,file_source) %>% 
  summarise(
    Total_Rooms = n_distinct(Room),                       # Total number of unique rooms
    Rooms_with_high_P = n_distinct(Room[P > 16000]),      # Count of rooms where P > 16000
    Proportion_high_P = Rooms_with_high_P / Total_Rooms   # Proportion of rooms with P > 16000
  ) %>% 
  ggplot()+geom_col(aes(x = Building, y = Proportion_high_P), position = "dodge")+
  facet_wrap(~file_source)+theme(axis.text.x =element_text(angle = 90))

ggsave("Building_comparison_prop_risky_beg.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 15, height = 12)




combined_data_mid %>% 
  group_by(Building,file_source) %>% 
  summarise(
    Total_Rooms = n_distinct(Room),                       # Total number of unique rooms
    Rooms_with_high_P = n_distinct(Room[P > 16000]),      # Count of rooms where P > 16000
    Proportion_high_P = Rooms_with_high_P / Total_Rooms   # Proportion of rooms with P > 16000
  ) %>% 
  ggplot()+geom_col(aes(x = Building, y = Proportion_high_P), position = "dodge")+
  facet_wrap(~file_source)+theme(axis.text.x =element_text(angle = 90))

ggsave("Building_comparison_prop_risky_mid.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 15, height = 12)


combined_data_end %>% 
  group_by(Building,file_source) %>% 
  summarise(
    Total_Rooms = n_distinct(Room),                       # Total number of unique rooms
    Rooms_with_high_P = n_distinct(Room[P > 16000]),      # Count of rooms where P > 16000
    Proportion_high_P = Rooms_with_high_P / Total_Rooms   # Proportion of rooms with P > 16000
  ) %>% 
  ggplot()+geom_col(aes(x = Building, y = Proportion_high_P), position = "dodge")+
  facet_wrap(~file_source)+theme(axis.text.x =element_text(angle = 90))

ggsave("Building_comparison_prop_risky_end.png", path = "/Users/courtney/GitHub/CLS_2024_Particle-risk-model/Figures", width = 15, height = 12)




