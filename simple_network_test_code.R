#{r libraries}
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

# Community model - simple SIR model
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

# This function sets the initial conditions for the building
#including randomly distributing individuals throughout the building/across rooms

#We assume that the proportion of S,I, and R in the community is proportional to 
#the proportion of S, I, and R individuals in the building
Bld_setup_func <- function(Community_output,day, delt,N_rooms){
  Building_ICs <- Community_output[day,] #Retrieves the number of S, I, and R individuals in the community at a particular day.
  
  #S_c + I_c + R_c = N_c | sum of the number of susceptible, infectious and recovered in the community is the total population of the community
  # we want the proportion of S, I, and R in the community to match the proportion of S, I, and R in the building:
  # (S_c + I_c + R_c)*delt = building_pop
  # we set building_pop based on how full we want the building to be: this is our adjusted_building capacity set above
  #i.e. Prop_Full*Building_Max_Capacity= Adjusted_building_capacity.
  # then building_pop = Adjusted_building_capacity so that we can get our desired amount of individuals in the building. 
  # and then we solve for delt above: delt = adjusted_building_pop/community_pop 
  
  Sb <- Community_output$S[day]*delt # number of susceptible in building
  Ib <- Community_output$I[day]*delt # number of Infectious in building
  Rb <- Community_output$R[day]*delt # number of Recovered in building
  
  #Sb, Ib, Rb are the number of Susceptible, Infected and Recovered individuals that will be in the building
  
  #initial conditions for each room. Randomly distribute individuals throughout rooms
  
  # first assign a random number between 0 and 1 for each room,
  #broken up by S, I, and R
  S_x <- c(runif(N_rooms, min = 0, max = 1))
  I_x <- c(runif(N_rooms, min = 0, max = 1))
  R_x <- c(runif(N_rooms, min = 0, max = 1))
  
  #normalize and then assign the correct number of S, I, and R based on Sb, Ib, and Rb
  S_x <- ((S_x/sum(S_x))*Sb)
  I_x <- ((I_x/sum(I_x))*Ib)
  R_x <- ((R_x/sum(R_x))*Rb)
  
  
  S_prop <- ((S_x/sum(S_x))*Sb)/(Sb+Ib+Rb)
  I_prop <- ((I_x/sum(I_x))*Ib)/(Sb+Ib+Rb)
  R_prop <- ((R_x/sum(R_x))*Rb)/(Sb+Ib+Rb)
  
  #this model is of the proportion of S, I and R in each room so now we need to switch from numbers to proportions.
  # but we need to keep the total number of individuals in each room since that is in our model
  Init_conds_nums <- data.frame(S_num = S_x, I_num = I_x, R_num = R_x)
  Init_conds_nums <- Init_conds_nums %>% mutate(N_x = S_num + I_num + R_num)
  Init_conds_props <- Init_conds_nums %>% 
    mutate(S_prop = S_num/N_x, I_prop = I_num/N_x, R_prop = R_num/N_x)%>%
    select(S_prop,I_prop,R_prop,N_x)
  
  # we start with no particles in the building
  P_x <- c(rep(0,N_rooms))
  
  Init_conds_props <- Init_conds_props %>% mutate(S_prop = S_prop,I_prop =I_prop,R_prop= R_prop, P_x =P_x,N_x=N_x)
  #Init_conds <-c(S=S_x/sum(Sb+Ib+Rb), I = I_x/sum(Sb+Ib+Rb), R = R_x/sum(Sb+Ib+Rb), P = P_x)
  
  return(data.frame(S=S_prop,I=I_prop,R=R_prop, P=P_x, N_x = S_x+I_x+R_x))
}


# creates a square matrix with entries being random numbers from 0 to 1 and 
#then eliminates transitions that don't exist based on the adjacency matrix 
#can be used for both people transition matrices and particle transitions matrices
# right now the people and particles use the same adjacency matrix -- this could be 
# changed by giving the function a different adjacency matrix - this would be good if
# we want to take into account ventilation routes.
Create_T_Matrix <-function(adjacency_matrix_to_use){
  #set.seed(123145) # <- easier for debugging
  T_mov <- data.frame(matrix(runif(N_rooms^2), nrow = N_rooms)) #populates a square matrix/dataframe with random numbers between 0 and 1 for the number of rooms that our building has.
  diag(T_mov) <- 0 #set diagonal to 0
  T_mov <- adjacency_matrix_to_use*T_mov #restrict the movement according to our network/adjacency matrix
  T_mov_norm <- t(apply(T_mov, 1, function(x) x / sum(x))) # normalize so that there aren't more people moving than what can (rows should sum to 1)
  T_mov <-T_mov_norm
  return(T_mov)
}



# 
# Create_Particle_T_Matrix<-function(adjacency_matrix_to_use,prop){
#   #set.seed(12312145) # <- easier for finding debugging
#   theta_mov <- data.frame(matrix(runif(N_rooms^2), nrow = N_rooms))
#   theta_mov <- adjacency_matrix_to_use*theta_mov
#   #diag(theta_mov) <- 0
#   theta_mov_prop <- prop #proportion of particles that move
#   theta_mov_norm <- t(apply(theta_mov, 1, function(x) x / sum(x)))
#   theta_mov <- theta_mov_prop*theta_mov_norm
#   
#   
#   return(theta_mov)
# }

# Functions for finding the flux in and flux out of a room - specifically for people - different for
# people/particles because only so many people can fit in a room however particles don't have a limit.
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
    for(i in 1:N_rooms){
      flux_in_temp <- flux_in_temp + State[i]*Transition_matrix[i,x]
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


Particle_model <- function(t, x, parms,T_mov, theta_mov, adjacency_matrix_to_use,C_x){
  ncompartment <- 5
  n_rooms <- length(x)/ncompartment
  S <- as.matrix(x[1:n_rooms])
  I <- as.matrix(x[(n_rooms+1):(2*n_rooms)])
  R <- as.matrix(x[(2*n_rooms+1):(3*n_rooms)])
  P <- as.matrix(x[(3*n_rooms+1):(4*n_rooms)])
  N_x <- as.matrix(x[(4*n_rooms+1):(5*n_rooms)])
  
  with(parms,{
    
    
    dS <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=S,Room_pops = N_x,Carrying_capacity = C_x)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=S,Room_pops = N_x,Carrying_capacity = C_x))
    
    dI <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=I,Room_pops = N_x,Carrying_capacity = C_x)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=I,Room_pops = N_x,Carrying_capacity = C_x))
    
    dR <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=R,Room_pops = N_x,Carrying_capacity = C_x)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=R,Room_pops = N_x,Carrying_capacity = C_x))
    #last step will be the particle EQ
    
    dP <- s*as.matrix(I)*as.matrix(N_x) - as.matrix(a*P/(lam*C_x*N_x))-d*as.matrix(P) + as.matrix(as.matrix(flux_in_particles(N_rooms, theta_mov,State = P)) - as.matrix(flux_out_particles(N_rooms, theta_mov,State = P)))
    
    dN_x <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=N_x,Room_pops = N_x,Carrying_capacity = C_x)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=N_x,Room_pops = N_x,Carrying_capacity = C_x))
    
    dt <- c(dS,dI,dR,dP,dN_x)
    return(list(dt))})
  
}
#   NEW   as.matrix(s*I) + as.matrix(flux_in_particles(N_rooms, theta_mov,State = P)) - as.matrix(flux_out_particles(N_rooms, theta_mov,State = P)) - d*as.matrix(P)

# OLD s*as.matrix(I) - a*as.matrix(P)*(S+I+R)+ as.matrix(as.matrix(flux_in_particles(N_rooms, theta_mov,State = P)) - as.matrix(flux_out_particles(N_rooms, theta_mov,State = P))) - d*as.matrix(P)


######################### Community simulation} ##############################
community_pop <- 100000
init_conds <- c(S = community_pop-1, I = 1, R = 0)
print(init_conds)
parms <- data.frame(bet =0.0000035 , gam = 1/21) 
print(parms)
times <- seq(from = 1, to = 144, by = 1)
print(times)
Community_output <- data.frame(lsoda(y = init_conds, func = SIR_community_model,times = times, parms=parms))
Community_output %>% pivot_longer(cols = !time) %>% arrange(desc(time))%>% 
  ggplot(aes(x=time,y =value, color = name))+geom_line()+
  scale_color_manual(name=NULL,values=c("blue","red","purple"),breaks = c("S","I","R"))+
  theme_classic()+
  labs(x= "Time (days)", y = "Number of individuals")+ggtitle("Community outbreak of pathogen")

############# Global variables #########
m <-5 #number of equations per room
day <- 31 # what day from the community model do we want to model and base our proportion of Susceptible, Infected, and Recovered individuals do we want to look at
Prop_full <- 0.8
#day_start <- 0

Maxtime <- 24*3
times <- seq(from = 0, to = Maxtime, by = 0.2)
prop_to_mov <- 0.9

# units of particle equation are not in raw number of particles but instead in volumetric air. 
# units are Liters/time 
#### Decay ####
# particles are estimated to fall within 20 minutes so 
# so d = 1/20

#### Absorption ####
# estimated 10L/min 
# going to use minutes
# going to use nanoparticles since that is what shedding is in

#### Shedding ####
#seems to be 3 nanoLiters/min
# if in L so going to convert to liters
# 10^9 nanoliters in 1 L
# so 3*10^-9 L per minute
# if hours that becomes 1.8*10^-7
#"correct units" = Liters = (s=1.8*10^-7,a=600, d=3)
#NanoLiters = (180,a =6*10^11, d=3)
#parms_Liters <- data.frame(s=1.8*10^-7,a=600, d=3)
#parms_nanoLiters <- data.frame(180,a =6*10^11, d=3)
parms <-data.frame(s=100,a=5, d=3,lam = 1)

####################  First example - church #############################

church_adjacency_matrix <- matrix(c(c(0,1,1,1,1),
                                    c(1,0,1,1,0),
                                    c(1,1,0,1,0),
                                    c(1,1,1,0,0),
                                    c(1,0,0,0,0)),nrow=5, ncol = 5)
church_graph<- graph_from_adjacency_matrix(church_adjacency_matrix, mode = "undirected")
plot(church_graph)
#church_adjacency_matrix
#church_init_prob <- c(0,0.2,0,0.8,rep(0,nrow(church_adjacency_matrix)-4))

graph_to_use <- church_graph
adjacency_matrix_to_use <- church_adjacency_matrix 
N_rooms <- ncol(adjacency_matrix_to_use) #number of rooms
Church_C<- c(100, #column 1 - main area / Hallway
       100, # column 2 - Hallway
       100, # column 3 - Hallway
       400, #column 4 Main room
       500
       ) 

#outside should only be able to hold the total capacity of the building 
Max_Building_Capacity <- sum(Church_C)
Adj_Max_Building_Capacity <- Prop_full*Max_Building_Capacity
delt <- Adj_Max_Building_Capacity/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# ### checking that initial conditions are what we want
# Building_prop_S <- sum(Church_setup$S)
# Building_prop_S
# 
# Community_prop_S <- Community_output$S[day]/(Community_output$S[day]+Community_output$I[day]+Community_output$R[day])
# Community_prop_S
# 
# Building_prop_I <- sum(Church_setup$I)
# Building_prop_I
# 
# Community_prop_I <- Community_output$I[day]/(Community_output$S[day]+Community_output$I[day]+Community_output$R[day])
# Community_prop_I
# 
# Building_prop_R <- sum(Church_setup$R)
# Building_prop_R
# 
# Community_prop_R <- Community_output$R[day]/(Community_output$S[day]+Community_output$I[day]+Community_output$R[day])
# Community_prop_R

# so far so good


church_graph <-church_graph %>% set_vertex_attr( "Infectious", value = Church_setup$I_x) %>% 
  set_vertex_attr("Room", value = seq(1:N_rooms))# + 
Church_Init_conds <-c(S=Church_setup$S, I = Church_setup$I, R = Church_setup$R, P = Church_setup$P, N_x =Church_setup$N_x)

# church_network <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(size = Infectious), color = "#C21807")+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_size(name = "Number of\nInfectious\npeople in\nrooms")+
#   labs(title = "Network representation of a church")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))
# church_network

# code to save figure
#png(filename = "Church_network.png",units="in", width=6, height=4, res=300)
#church_network
#dev.off()

#Initialize Transition matrices, both people and particles
Church_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use)
#Church_T_mov <- adjacency_matrix_to_use*0.5

Church_theta_mov <- Create_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use)




Church_output <- data.frame(lsoda(y = Church_Init_conds, func = Particle_model,times = times,
                                  parms = parms,
                                  adjacency_matrix_to_use=adjacency_matrix_to_use,
                                  theta_mov =Church_theta_mov,
                                  T_mov = Church_T_mov, 
                                  C_x=Church_C))

Church_data_clean <- Church_output%>% pivot_longer(cols = !time,
                                                   names_to = c("State", "Room"),
                                                   names_pattern = "([A-Za-z]+)(\\d+)",
                                                   values_to = "Number")

Church_data_clean %>% filter(State == "S") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Proportion of Susceptible individuals")+ggtitle("Proportion of Susceptible individuals across rooms")

Church_data_clean %>% filter(State == "I") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Proportion of Infectious individuals")+ggtitle("Proportion of Infectious individuals across rooms")
Church_data_clean %>% filter(State == "R") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Proportion of Recovered individuals")+ggtitle("Proportion of Recovered individuals across rooms")
Church_data_clean %>% filter(State == "P") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of infectious particles")+ggtitle("Infectious particles in rooms within a building")
Church_data_clean %>% filter(State == "x") %>% 
  group_by(State,Room)%>% 
  ggplot( aes(x=time, y = Number, group =Room, color = Room))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of individuals")+ggtitle("Number of individuals in rooms within a building")

# Church_data_clean %>% filter(State != "P") %>% 
#   group_by(time,Room)%>% summarise(total_room_pop=sum(Number)) %>% group_by(time) %>% summarise(tot_bld_pop =sum(total_room_pop)) %>% 
#   ggplot( aes(x=time, y = tot_bld_pop))+
#   geom_line()+theme_classic()+
#   labs(x = "Time (hours)", y= "Number of people")+ggtitle("number of people in building")

Church_data_clean %>% filter(State == "S"| State == "I" | State == "R") %>% ungroup() %>% group_by(time, Room) %>% summarise(Total_prop=sum(Number)) %>%
  ggplot(aes(x=time, y=Total_prop, color=Room))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in rooms within a building")

Church_data_clean %>% filter(State == "S"| State == "I" | State == "R") %>% ungroup() %>% group_by(time) %>% summarise(Total_prop=sum(Number)) %>%
  ggplot(aes(x=time, y=Total_prop))+geom_line()+
  labs(x = "Time (hours)", y= "Proportion of individuals")+ggtitle("Proportion of individuals in the building")


Church_data_clean %>% filter(State== "x")%>% ungroup() %>% group_by(time) %>% summarise(Total_pop=sum(Number)) %>%
  ggplot(aes(x=time, y=Total_pop))+geom_line()+
  labs(x = "Time (hours)", y= "Number of individuals")+ggtitle("Number of individuals in the building")

Church_data_ratios <- Church_data_clean %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% group_by(time,Room)  %>% 
  mutate(N_x = S+I+R) %>% group_by(time,Room) %>% 
  mutate(K_x = ((parms$s)*I)/((parms$a)*N_x), prop_to_K = P/K_x)
Church_P_x_K_x_plot<-ggplot(Church_data_ratios,aes(x =time, y =prop_to_K,group= Room,color=Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")
Church_P_x_K_x_plot
Church_P_x_plot <- ggplot(Church_data_ratios, aes(x=time, y = P, color = Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")
Church_P_x_plot


Particles_in_rooms_data<- Church_data_clean %>% filter(State == "P") %>% group_by(Room,time) %>% summarise(Room = as.numeric(Room),Room_Pop =sum(Number))
Particles_in_rooms_data$Room<-as.factor(Particles_in_rooms_data$Room)

label_data <-  Particles_in_rooms_data %>%
  distinct(Room, .keep_all = TRUE)

church_model_output_fig<-ggplot(Particles_in_rooms_data,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
  labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
  theme_linedraw()+
  geom_text_repel(data=label_data,aes(label = Room), show.legend = FALSE)+
  theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
        plot.title.position = "plot",
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10),
        legend.key.size = unit(1,units = "cm"))

church_model_output_fig
# png(filename = "Church_model_output.png",units="in", width=6, height=4, res=300)
# church_model_output_fig
# dev.off()
############### plotting particles
ggplot(Particles_in_rooms_data, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
  geom_line() +
  geom_text(data = Particles_in_rooms_data %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
  labs(title = "Particle accumulation over time", x = "Time", y = "Number of infectious Particles", color = "Rooms") +
  theme_minimal() +
  theme(legend.position = "right")

################## plotting people
Total_people_In_rooms<- Church_data_clean %>% filter(State != "P") %>% group_by(Room,time) %>% summarise(Room = as.numeric(Room),Room_Pop =sum(Number))

Total_people_In_rooms$Room<-as.factor(Total_people_In_rooms$Room)

label_data <-  Total_people_In_rooms %>%
  distinct(Room, .keep_all = TRUE)

ggplot(Total_people_In_rooms, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
  geom_line() +
  geom_text(data = Total_people_In_rooms %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
  labs(title = "People in rooms over time", x = "Time", y = "Number of people in rooms", color = "Rooms") +
  theme_minimal() +
  theme(legend.position = "right")+ylim(0,100)
################# Facet_wrap - Particles

ggplot(Total_people_In_rooms, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
  geom_line() +
  facet_wrap(~Room)+
  #geom_text(data = Total_people_In_rooms %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
  labs(title = "People in rooms over time", x = "Time", y = "Number of people in rooms", color = "Rooms") +
  theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
        plot.title.position = "plot",
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10),
        legend.key.size = unit(1,units = "cm"))+
  theme_linedraw()+xlim(0,10)

################# Facet_wrap - Particles

ggplot(Particles_in_rooms_data,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
  labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
  facet_wrap(~Room)+
  theme_linedraw()+
  #geom_text_repel(data=Particles_in_rooms_data,aes(label = Room), show.legend = FALSE)+
  theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
        plot.title.position = "plot",
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10),
        legend.key.size = unit(1,units = "cm"))#+xlim(0,24*60)
# calculating levels of risk throughout the building. Potential calculations are:
# Number of infectious particles at the end of the day in each room
# Average number of infectious particles throughout the day, in each room
# Cumulative number of infectious particles in each room throughout the day
# Average risk in the building - Total number of infectious particles in the whole building divided by the number of rooms OR
# Number of infectious particles in each room divided by total in building then take the average of that quantity
# Total number of infectious particles in each room divided by their saturation point - this





