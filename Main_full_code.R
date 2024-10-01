#Code being added to github and now going to practice better version control

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
library(ideanet)
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
  Building_ICs <- Community_output[day,] #Retrieves the number of S, I, and R individuals in the population at a particular day.
  
  #delt is set below in the parameter declaration by taking into account 
  #the max capacity of the building and how full we set the building to be
  Sb <- Community_output$S[day]*delt 
  Ib <- Community_output$I[day]*delt
  Rb <- Community_output$R[day]*delt
  
  #Sb, Ib, Rb are the number of Susceptible, Infected and Recovered individuals that will be in the building
  
  #initial conditions for each room. Randomly distribute individuals throughout rooms
  
  # first asign a random number between o and one for each room,
  #broken up by S, I, and R
  S_x <- c(runif(N_rooms, min = 0, max = 1))
  I_x <- c(runif(N_rooms, min = 0, max = 1))
  R_x <- c(runif(N_rooms, min = 0, max = 1))
  
  #normalize and then assign the correct amount of S, I, and R based on Sb, Ib, and Rb
  S_x <- (S_x/sum(S_x))*Sb
  I_x <- (I_x/sum(I_x))*Ib
  R_x <- (R_x/sum(R_x))*Rb
  
  # we start with no particles in the building
  P_x <- c(rep(0,N_rooms))
  
  Init_conds <-c(S=S_x, I = I_x, R = R_x, P = P_x)
  
  return(data.frame(S_x=S_x,I_x=I_x,R_x=R_x,P_x=P_x))
}


# creates a square matrix with entries being random numbers from 0 to 1 and 
#then eliminates transitions that don't exist based on the adjacency matrix 
#can be used for both people transition matrices and particle transitions matrices
# right now the people and particles use the same adjacency matrix -- this could be 
# changed by giving the function a different adjacency matrix - this would be good if
# we want to take into account ventilation routes.
Create_T_Matrix <-function(adjacency_matrix_to_use){
  #set.seed(123145) # <- easier for finding debugging
  T_mov <- data.frame(matrix(runif(N_rooms^2), nrow = N_rooms))
  T_mov <- adjacency_matrix_to_use*T_mov
  T_mov_norm <- t(apply(T_mov, 1, function(x) x / sum(x)))
  T_mov <-T_mov_norm
  return(T_mov)
}

#function that uses the adjacency matrix to calculate distances and then sets transition rates
#such that moving closer to the outside room has higher rates. 
End_of_day_T_function <- function(adjacency_matrix_to_use){
  #set.seed(12312145)
  Random_T_mov <- Create_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use)
  End_day_T_mov <- Random_T_mov
  distance_vector <- distances(graph_to_use, to = as.character(nrow(adjacency_matrix_to_use)))
  for(i in 1:N_rooms){
    for(j in 1:N_rooms){
      if (distance_vector[j] < distance_vector[i]) {
        End_day_T_mov[i,j] <- 10
      }else{
        End_day_T_mov[i,j] <- 0.5
      }
    }
  }
  End_day_T_mov <- End_day_T_mov*adjacency_matrix_to_use
  return(End_day_T_mov)
}

# 
Create_Particle_T_Matrix<-function(adjacency_matrix_to_use,prop){
  #set.seed(12312145) # <- easier for finding debugging
  theta_mov <- data.frame(matrix(runif(N_rooms^2), nrow = N_rooms))
  theta_mov <- adjacency_matrix_to_use*theta_mov
  #diag(theta_mov) <- 0
  theta_mov_prop <- prop #proportion of particles that move
  theta_mov_norm <- t(apply(theta_mov, 1, function(x) x / sum(x)))
  theta_mov <- theta_mov_prop*theta_mov_norm
  
  
  return(theta_mov)
}

# Functions for finding the flux in and flux out of a room - specifically for people - different for
# people/particles because only so many people can fit in a room however particles don't have a limit.
flux_in_people <- function(N_rooms, Transition_matrix,State,Room_pops,Carrying_capacity,t){
  all_room_change <- c(seq(N_rooms))
  for(x in 1:N_rooms){
    flux_in_temp <- 0
    for(i in 1:N_rooms){
      flux_in_temp <- flux_in_temp + State[i]*Transition_matrix[i,x]*(1-(Room_pops[x]/Carrying_capacity[x]))
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
    for(i in 1:N_rooms){
      flux_out_temp <- flux_out_temp + State[x]*Transition_matrix[x,i]*(1-(Room_pops[i]/Carrying_capacity[i]))
    }
    all_room_change[x] <- flux_out_temp
  }
  test <- c(M = as.vector(all_room_change))
  return(test)
}


# Functions for finding the flux in and the flux out for particles
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


Particle_model <- function(t, x, parms,T_mov, theta_mov, adjacency_matrix_to_use,C){
  ncompartment <- 4
  n_rooms <- length(x)/ncompartment
  S <- as.matrix(x[1:n_rooms])
  I <- as.matrix(x[(n_rooms+1):(2*n_rooms)])
  R <- as.matrix(x[(2*n_rooms+1):(3*n_rooms)])
  P <- as.matrix(x[(3*n_rooms+1):(4*n_rooms)])
  
  with(parms,{
    
    
    dS <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=S,Room_pops = (S+I+R),Carrying_capacity = C)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=S,Room_pops = (S+I+R),Carrying_capacity = C))
    dI <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=I,Room_pops = (S+I+R),Carrying_capacity = C)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=I,Room_pops = (S+I+R),Carrying_capacity = C))
    dR <- as.matrix((flux_in_people(N_rooms, Transition_matrix =T_mov, State=R,Room_pops = (S+I+R),Carrying_capacity = C)) - flux_out_people(N_rooms, Transition_matrix = T_mov, State=R,Room_pops = (S+I+R),Carrying_capacity = C))
    #last step will be the particle EQ
    dP <- s*as.matrix(I) - a*as.matrix(P)*(S+I+R)+ 
      as.matrix(as.matrix(flux_in_particles(N_rooms, theta_mov,State = P)) - as.matrix(flux_out_particles(N_rooms, theta_mov,State = P))) - d*as.matrix(P)
    dt <- c(dS,dI,dR,dP)
    return(list(dt))})
  
}


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
m <-4 #number of equations per room
day <- 31 # what day from the community model do we want to model and base our proportion of Susceptible, Infected, and Recovered individuals do we want to look at
Prop_full <- 0.8
day_start <- 0
day_duration <- 8
Maxtime <- 24*7
times <- seq(from = 0, to = Maxtime, by = 0.2)

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
parms <-data.frame(s=100,a=5, d=3)

####################  First example - church #############################

church_adjacency_matrix <- matrix(c(c(0,rep(1,3),rep(0,8),rep(1,4),rep(0,31-17),1), #column 1 - main area / Hallway
                                    c(1,0,0,1,rep(0,12),rep(1,4),rep(0,31-20)), # column 2 - Hallway
                                    c(1,0,0,1,rep(0,4),rep(1,4),rep(0,31-12)), # column 3 - Hallway
                                    c(1,1,1,0,1,0,1,1,rep(0,31-8)), #column 4 Main room
                                    c(rep(0,3),1,0,1,0,1,rep(0,17),1,rep(0, 31-26)), # # column 5 Hallway
                                    c(rep(0,4),1,0,1,rep(0,15),rep(1,5),0,1,1,0), # Column 6 Hallway
                                    c(rep(0,3),1,0,1,0,1,rep(0,31-8)), # Column 7 Hallway
                                    c(rep(0,3),1,1,0,1,rep(0,21-8),1,1,rep(0,31-22)), # Column 8 Hallway
                                    rep(c(0,0,1,rep(0,31-3)),4), # columns 9-12
                                    rep(c(1,rep(0,31-1)),4), # columns 13-16
                                    rep(c(0,1,rep(0,31-2)),4), # columns 17-20
                                    rep(c(rep(0,7),1,rep(0, 31-8)),2), # Columns 21 and 22
                                    rep(c(rep(0,5),1,rep(0,31-6)),3), # columns 23-25
                                    c(rep(0,4),1,1,rep(0,31-6)), # Column 26
                                    c(rep(0,5),1,rep(0,21),1,rep(0,31-28)), # column 27
                                    c(rep(0,26),1,rep(0,31-27)), # column 28
                                    c(rep(0,5),1,rep(0,31-6)), # column 29 
                                    c(rep(0,5),1,1,rep(0,31-7)), # column 30
                                    c(1,rep(0,30)) # column 31 -- Outside
),nrow=31, ncol = 31)
church_graph<- graph_from_adjacency_matrix(church_adjacency_matrix, mode = "undirected")
plot(church_graph)
#church_adjacency_matrix
#church_init_prob <- c(0,0.2,0,0.8,rep(0,nrow(church_adjacency_matrix)-4))

graph_to_use <- church_graph
adjacency_matrix_to_use <- church_adjacency_matrix 
N_rooms <- ncol(adjacency_matrix_to_use) #number of rooms
Church_C <- c(20, #column 1 - main area / Hallway
              5, # column 2 - Hallway
              5, # column 3 - Hallway
              200, #column 4 Main room
              5, # # column 5 Hallway
              7, # Column 6 Hallway - slightly bigger hallway
              5, # Column 7 Hallway
              5, # Column 8 Hallway
              15, # comn 9
              15, # column 10
              15, #column 11
              2, # column 12 janitors closet
              15, # column 13
              5, # column 14 - bathroom
              15, # column 15 
              5, # column 16 - bathroom
              15, # column 17
              2, # column 18 - closet
              15, # column 19
              15, # column 20
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
) # column 31 -- Outside


Church_C <- c(Church_C,sum(Church_C))
#outside should only be able to hold the total capacity of the building 
Building_max <- Church_C[N_rooms] #last room in C is the "outside" room
Max_capacity <- Prop_full*Building_max
delt <- Max_capacity/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
church_graph <-church_graph %>% set_vertex_attr( "Infectious", value = Church_setup$I_x) %>% 
  set_vertex_attr("Room", value = seq(1:N_rooms))# + 
Church_Init_conds <-c(S=Church_setup$S_x, I = Church_setup$I_x, R = Church_setup$R_x, P = Church_setup$P_x)

church_network <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = Infectious), color = "#C21807")+
  geom_nodetext_repel(aes(label=Room))+
  theme_blank()+
  scale_size(name = "Number of\nInfectious\npeople in\nrooms")+
  labs(title = "Network representation of a church")+
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust=1))
church_network
# code to save figure
#png(filename = "Church_network.png",units="in", width=6, height=4, res=300)
#church_network
#dev.off()

#### Church Network metrics ####
# temp_church_adj_matrix <- church_adjacency_matrix
# colnames(temp_church_adj_matrix) <- c(paste("V",1:ncol(temp_church_adj_matrix),sep = ""))
# 
# church_metrics <- netwrite(data_type = c("adjacency_matrix"),
#                            adjacency_matrix = temp_church_adj_matrix,
#                            directed = FALSE,
#                            net_name = "church_net")
# 
# list2env(church_metrics, .GlobalEnv)
# node_measures
# church_network_measures_plot<-church_metrics$node_measure_plot
# 
# png(filename = "church_network_measures_plot.png",units="in", width=10, height=10, res=300)
# church_network_measures_plot
# dev.off()
# 
# table_data<- church_metrics$system_level_measures
# table_data <- table_data %>% select()
# church_table <- kable(church_metrics$system_level_measures, row.names = FALSE)
# church_table <- kable(church_metrics$system_level_measures, format = "latex", booktabs = TRUE, row.names = FALSE,longtable=TRUE) %>% 
#   kable_styling(latex_options = c("scale_down"))
# 
# church_latex_table <- capture.output(print(church_table))
# writeLines(church_latex_table, "church_table.tex")


#Initialize Transition matrices, both people and particles
Church_T_mov <- Create_Particle_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use,prop = 0.9)
#Church_T_mov <- adjacency_matrix_to_use*0.5


Church_theta_mov <- Create_Particle_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use,prop = 0.5)

Church_end_day_T_mov <- End_of_day_T_function(adjacency_matrix_to_use = adjacency_matrix_to_use)
End_of_day_T_mov_to_use <- Church_end_day_T_mov

Church_output <- data.frame(lsoda(y = Church_Init_conds, func = Particle_model,times = times, parms = parms,adjacency_matrix_to_use=adjacency_matrix_to_use,theta_mov =Church_theta_mov,T_mov = Church_T_mov,C =Church_C))

clean_data <- Church_output %>% pivot_longer(cols = !time,
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

clean_data %>% filter(State != "P") %>% 
  group_by(time,Room)%>% summarise(total_room_pop=sum(Number)) %>% group_by(time) %>% summarise(tot_bld_pop =sum(total_room_pop)) %>% 
  ggplot( aes(x=time, y = tot_bld_pop))+
  geom_line()+theme_classic()+
  labs(x = "Time (hours)", y= "Number of infectious particles")+ggtitle("Infectious particles in rooms within a building")


#debugging_15 <- clean_data %>% filter(Room == 19 ) %>% filter(State == "P") 

#ggplot(debugging_15, aes(x=time, y = Number, color = State))+geom_line()

Church_data_clean <- Church_output%>% pivot_longer(cols = !time,
                                                   names_to = c("State", "Room"),
                                                   names_pattern = "([A-Za-z]+)(\\d+)",
                                                   values_to = "Number")

Church_data_ratios <- Church_data_clean %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% group_by(time,Room)  %>% 
  mutate(N_x = S+I+R) %>% group_by(time,Room) %>% 
  mutate(K_x = ((parms$s)*I)/((parms$a)*N_x), prop_to_K = P/K_x)
Church_P_x_K_x_plot<-ggplot(Church_data_ratios,aes(x =time, y =prop_to_K,group= Room,color=Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")

Church_P_x_plot <- ggplot(Church_data_ratios, aes(x=time, y = P, color = Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")
Church_P_x_K_x_plot

#### OLD figure making (church ) ####
# 
# Particles_in_rooms_data<- Build_data_mod %>% filter(State == "P") %>% group_by(Room,time) %>% summarise(Room = as.numeric(Room),Room_Pop =sum(Number))
# Particles_in_rooms_data$Room<-as.factor(Particles_in_rooms_data$Room)
# 
# label_data <-  Particles_in_rooms_data %>%
#   distinct(Room, .keep_all = TRUE)
# church_model_output_fig<-ggplot(Particles_in_rooms_data,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   theme_linedraw()+
#   geom_text_repel(data=label_data,aes(label = Room), show.legend = FALSE)+
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# 
# church_model_output_fig
# # png(filename = "Church_model_output.png",units="in", width=6, height=4, res=300)
# # church_model_output_fig
# # dev.off()
# ############### plotting particles
# ggplot(Particles_in_rooms_data, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
#   geom_line() +
#   geom_text(data = Particles_in_rooms_data %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   labs(title = "Particle accumulation over time", x = "Time", y = "Number of infectious Particles", color = "Rooms") +
#   theme_minimal() +
#   theme(legend.position = "right")
# 
# ################## plotting people
# Total_people_In_rooms<- Build_data_mod %>% filter(State != "P") %>% group_by(Room,time) %>% summarise(Room = as.numeric(Room),Room_Pop =sum(Number))
# 
# Total_people_In_rooms$Room<-as.factor(Total_people_In_rooms$Room)
# 
# label_data <-  Total_people_In_rooms %>%
#   distinct(Room, .keep_all = TRUE)
# 
# ggplot(Total_people_In_rooms, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
#   geom_line() +
#   geom_text(data = Total_people_In_rooms %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   labs(title = "People in rooms over time", x = "Time", y = "Number of people in rooms", color = "Rooms") +
#   theme_minimal() +
#   theme(legend.position = "right")+ylim(0,100)
# ################# Facet_wrap - Particles
# 
# ggplot(Total_people_In_rooms, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
#   geom_line() +
#   facet_wrap(~Room)+
#   #geom_text(data = Total_people_In_rooms %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   labs(title = "People in rooms over time", x = "Time", y = "Number of people in rooms", color = "Rooms") +
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))+
#   theme_linedraw()
# 
# ################# Facet_wrap - Particles
# 
# ggplot(Particles_in_rooms_data,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   facet_wrap(~Room)+
#   theme_linedraw()+
#   #geom_text_repel(data=Particles_in_rooms_data,aes(label = Room), show.legend = FALSE)+
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# # calculating levels of risk throughout the building. Potential calculations are:
# # Number of infectious particles at the end of the day in each room
# # Average number of infectious particles throughout the day, in each room
# # Cumulative number of infectious particles in each room throughout the day
# # Average risk in the building - Total number of infectious particles in the whole building divided by the number of rooms OR
# # Number of infectious particles in each room divided by total in building then take the average of that quantity
# # Total number of infectious particles in each room divided by their saturation point - this
# 
# 
# #Number of infectious particles in each room over time for the last day - time >=48
# N_Inf_part_throughout_day <- filter(Particles_in_rooms_data, time >= 48)
# N_Inf_part_throughout_day$Room <- as.factor(N_Inf_part_throughout_day$Room)
# 
# # All on one
# ggplot(N_Inf_part_throughout_day,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   #facet_wrap(~Room)+
#   theme_linedraw()+
#   geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# #Each room with their own graph
# Church_N_inf_part_by_room_throughout_day<- ggplot(N_Inf_part_throughout_day,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   facet_wrap(~Room)+
#   theme_linedraw()+
#   #geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "none",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# 
# Church_N_inf_part_by_room_throughout_day
# 
# png(filename = "Church_N_inf_part_by_room_throughout_day.png",units="in", width=10, height=10, res=300)
# Church_N_inf_part_by_room_throughout_day
# dev.off()
# #Now taking into account the saturation point
# sat_points <- data.frame(Room = as.factor(1:N_rooms), Sat_point = K)
# N_Inf_part_rooms_sat_point <- right_join(N_Inf_part_throughout_day,sat_points, by = join_by(Room)) %>%
#   mutate(Prop_to_Sat_point = Room_Pop/Sat_point)
# Church_Ratio_inf_part_by_room_throughout_day <-ggplot(N_Inf_part_rooms_sat_point, aes(x = time, y= Prop_to_Sat_point,group=Room,color= Room, label= Room))+geom_line()+
#   labs(x = "Time (hours)", y= "Ratio",title = "Ratio of infectious particles in room to the threshold number of particles of each room")+
#   facet_wrap(~Room)+
#   theme_linedraw()+
#   #geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "none",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# Church_Ratio_inf_part_by_room_throughout_day
# png(filename = "Church_Ratio_inf_part_by_room_throughout_day.png",units="in", width=10, height=10, res=300)
# Church_Ratio_inf_part_by_room_throughout_day
# dev.off()
# 
# # Number of rooms over sat point by 10 percent
# N_rooms_over_sat_point_10 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.10)))
# N_rooms_over_sat_point_10 <- N_rooms_over_sat_point_10 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_10, aes(x=time, y =n))+ geom_line(color = "turquoise")
# 
# # Number of rooms over sat point by 30 percent
# N_rooms_over_sat_point_30 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.30)))
# N_rooms_over_sat_point_30 <- N_rooms_over_sat_point_30 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_30, aes(x=time, y =n))+ geom_line( color = "purple")
# 
# # Number of rooms over sat point by 50 percent
# N_rooms_over_sat_point_50 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.50)))
# N_rooms_over_sat_point_50 <- N_rooms_over_sat_point_50 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_50, aes(x=time, y =n))+ geom_line(color = "pink")
# 
# # Number of rooms over sat point by 70 percent
# N_rooms_over_sat_point_70 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.70)))
# N_rooms_over_sat_point_70 <- N_rooms_over_sat_point_70 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_70, aes(x=time, y =n))+ geom_line(color  = "orange")
# 
# # Number of rooms over sat point by 90 percent
# N_rooms_over_sat_point_90 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.90)))
# N_rooms_over_sat_point_90 <- N_rooms_over_sat_point_90 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_90, aes(x=time, y =n))+ geom_line(color = "green4")
# 
# # Number of rooms over sat point by 150 percent
# N_rooms_over_sat_point_150 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*1.5)))
# N_rooms_over_sat_point_150 <- N_rooms_over_sat_point_150 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_150, aes(x=time, y =n))+ geom_line(color = "red4")
# 
# # Number of rooms over sat point by 200 percent
# N_rooms_over_sat_point_200 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*2)))
# N_rooms_over_sat_point_200 <- N_rooms_over_sat_point_200 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_200, aes(x=time, y =n))+ geom_line(color= "blue4")
# 
# data_to_plot <- as.data.frame(cbind(time =N_rooms_over_sat_point_10$time,
#                                     Over_10 = N_rooms_over_sat_point_10$n,
#                                     Over_30 = N_rooms_over_sat_point_30$n,
#                                     Over_50 = N_rooms_over_sat_point_50$n,
#                                     Over_70 = N_rooms_over_sat_point_70$n,
#                                     Over_90 = N_rooms_over_sat_point_90$n,
#                                     Over_150 = N_rooms_over_sat_point_150$n,
#                                     Over_200 = N_rooms_over_sat_point_200$n))
# data_to_plot <- pivot_longer(data_to_plot, cols= starts_with("Over"),names_to = "Threshold", values_to = "Number")
# data_to_plot_V2 <- data_to_plot %>% mutate(Percent_over = (Number/N_rooms)*100)
# neworder <- c("Over_10","Over_30","Over_50", "Over_70", "Over_90","Over_150", "Over_200")
# data_to_plot_V2 <- arrange(transform(data_to_plot_V2,
#                                      Threshold=factor(Threshold,levels=neworder)),Threshold)
# 
# Church_Percent_above_sat <- ggplot(data_to_plot_V2, aes(x= time, y = Percent_over, color = Threshold))+geom_line()+
#   facet_wrap(~Threshold)+
#   theme_linedraw()+
#   labs(x = "Time (Hours)", y = "Percent of\nrooms over\nsaturation point", title = "Percentage of rooms over saturation points at different thresholds")+
#   ylim(0,50)+
#   scale_color_manual(name = "Percent over\nthreshold",
#                      labels = c("10%","30%","50%","30%","50%","70%","90%","150%","200%"),
#                      values = c("blue4","red4","green4","orange","pink","purple","turquoise"))+
#   #geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# 
# Church_Percent_above_sat
# 
# 
# png(filename = "Church_Percent_above_sat.png",units="in", width=10, height=10, res=300)
# Church_Percent_above_sat
# dev.off()
# 
# 
# # starting with just the number of particles in each room at the beginning of the last day.
# N_Inf_part_rooms_beg <- filter(Particles_in_rooms_data, time == 48)
# N_Inf_part_rooms_beg$Room <- as.numeric(N_Inf_part_rooms_beg$Room)
# N_Inf_part_rooms_beg <- N_Inf_part_rooms_beg%>% arrange(.by_group = TRUE)
# ggplot(N_Inf_part_rooms_beg, aes(x = Room, y = Room_Pop))+ geom_col()+
#   labs(x="Room",y="Number of infectious particles", title = "Number of infectious particles in each room at the end of the day")
# 
# 
# #Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# room_particle_graph <-church_graph %>% set_vertex_attr( "Infectious", value =N_Inf_part_rooms_beg$Room_Pop ) %>% 
#   set_vertex_attr("Room", value = seq(1:N_rooms))# + 
# Church_Init_conds <-c(S=Church_setup$S_x, I = Church_setup$I_x, R = Church_setup$R_x, P = Church_setup$P_x)
# 
# Church_Network_With_N_inf <- ggplot(room_particle_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(color = Infectious),size =2)+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_color_continuous(name = "Number of\nInfectious\nParticles")+
#   labs(title = "Number of infectious particles in each room at the beginning of the day")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))
# Church_Network_With_N_inf
# 
# png(filename = "Church_Network_With_N_inf.png",units="in", width=10, height=10, res=300)
# Church_Network_With_N_inf
# dev.off()
# # Proportion of infectious particles at the end of the last day.
# N_Inf_part_rooms_end <- filter(Particles_in_rooms_data, time == 56) %>% ungroup() %>%mutate(Total = sum(Room_Pop),Prop=Room_Pop/sum(Room_Pop))
# N_Inf_part_rooms_end$Room <- as.factor(N_Inf_part_rooms_end$Room)
# ggplot(N_Inf_part_rooms_end, aes(x = Room, y = Prop))+ geom_col()+
#   labs(x="Room",y="Proportion of infectious particles", title = "Proportion of infectious particles in each room at the end of the day")
# 
# 
# #Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# room_particle_graph <-church_graph %>% set_vertex_attr( "Infectious", value =N_Inf_part_rooms_end$Prop ) %>% 
#   set_vertex_attr("Room", value = seq(1:N_rooms))# + 
# 
# ggplot(room_particle_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(color = Infectious))+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_color_continuous(name = "Number of\nInfectous\nParticles")+
#   labs(title = "Raw number of infectious particles in each room")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))
# 
# # Average number of particles throughout the day (but on the last day)
# N_Inf_part_rooms_day_avg <- filter(Particles_in_rooms_data, time >= 48 & time <= 56)%>% summarise(Total_inf_part_room=sum(Room_Pop)) %>% mutate(Day_avg_inf_part = Total_inf_part_room/8)
# #%>% ungroup() %>%mutate(Total = sum(Room_Pop),Prop=Room_Pop/sum(Room_Pop))
# N_Inf_part_rooms_day_avg$Room <- as.factor(N_Inf_part_rooms_day_avg$Room)
# ggplot(N_Inf_part_rooms_end, aes(x = Room, y = Prop))+ geom_col()+
#   labs(x="Room",y="Proportion of infectious particles", title = "Proportion of infectious particles in each room at the end of the day")
# 
# 
# #Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# room_particle_graph <-church_graph %>% set_vertex_attr( "Infectious", value =N_Inf_part_rooms_day_avg$Prop ) %>% 
#   set_vertex_attr("Room", value = seq(1:N_rooms))# + 
# 
# ggplot(room_particle_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(color = Infectious))+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_color_continuous(name = "Number of\nInfectous\nParticles")+
#   labs(title = "Raw number of infectious particles in each room")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))



#### Office example ####
Office_adjacency_matrix <- matrix(c(c(0,1,rep(0,4),1,0,1,rep(0,36-9),1), # Hallway 1
                                    c(1,0,1,rep(0,32),1,1), # Hallway 2
                                    c(0,1,0,1,rep(0,5),1,1,1,rep(0,37-12)), # Hallway/room 3
                                    c(0,0,1,rep(0,6),rep(1,9),0,0,rep(1,3),rep(0,37-23)), # Hallway 4
                                    c(rep(0,3),1,0,1,1,rep(0,13),rep(1,5),rep(0,37-25)), # Hallway 5
                                    c(rep(0,4),1,0,1,1, rep(0,17),rep(1,4), rep(0,37-29)), # Hallway 6
                                    c(1,rep(0,3),1,1,0,rep(0,37-7)), # Hallway/room 7
                                    c(rep(0,5),1,rep(0,22),rep(1,3),rep(0,37-31)), # Hallway 8
                                    c(1,rep(0,30),1,1,1,1,0,0), # Hallway 9
                                    rep(c(0,0,1,1,rep(0,37-4)),3), # Rooms 10-12
                                    rep(c(0,0,0,1,rep(0,37-4)),5), # Rooms 13-17
                                    c(0,0,0,1,rep(0,14),1,1,rep(0,37-20)), # Room 18
                                    rep(c(rep(0,17),1,rep(0,37-18)),2), # Rooms 19 and 20
                                    rep(c(rep(0,3),1,1,rep(0,37-5)),3), # Rooms 21-23
                                    rep(c(rep(0,4),1,rep(0,37-5)),2), # Rooms 24 and 25
                                    rep(c(rep(0,5),1,rep(0,37-6)),3), # Rooms 26-28
                                    c(rep(0,5),1,0,1,rep(0,37-8)), # Room 29
                                    rep(c(rep(0,7),1,rep(0,37-8)),2), # Rooms 30 and 31
                                    rep(c(rep(0,8),1,rep(0,37-9)),3), # Rooms 32-34
                                    c(rep(0,8),1,rep(0,36-9),1), # Room 35
                                    c(0,1,rep(0,34),1), # Room 36
                                    c(1,1,rep(0,32),1,1,0)), # Room 37 --- :OUTSIDE: ---
                                  nrow=37,ncol = 37)

Office_graph<- graph_from_adjacency_matrix(Office_adjacency_matrix, mode = "undirected")
plot(Office_graph)
graph_to_use <- Office_graph
adjacency_matrix_to_use <- Office_adjacency_matrix
N_rooms <- ncol(adjacency_matrix_to_use) #number of rooms
Office_C <- c(
  200,# room 1 - Hallway/lobby
  200,# room 2 - Hallway/lobby
  30, # room 3 - Large Common room
  20,# room 4 - Hallway 
  20,# room 5 - Hallway
  20,# room 6 - Hallway
  2,# room 7 - Intermediate room
  10,# room 8 - Hallway- small
  10,# room 9 - Hallway to restrooms and stairway
  30,# room 10 - Conference room - large
  2,# room 11 - Storage room
  30,# room 12 - Conference room - large
  3,# room 13 - Office room
  10,# room 14 - Conference room - small
  3,# room 15 - Office room 
  3,# room 16 - Office room
  4,# room 17 - Kitchenette
  20,# room 18 - Large communal Office space
  3,# room 19 - small communal Office space
  2,# room 20 - Storage room
  rep(3,11), # rooms 21-31 - Offices
  5, # room 32 - Restroom
  5, # room 33 - Restroom
  2, # room 34 - Storage room
  20, # room 35 - Stairs
  10 # room 36 - Elevators
  # room 37 - "Outside"
)
Office_C <- c(Office_C,sum(Office_C))
#outside should only be able to hold the total capacity of the building 
Building_max <- Office_C[N_rooms] #last room in C is the "outside" room
Max_capacity <- Prop_full*Building_max
delt <- Max_capacity/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

Office_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
Office_Init_conds <-c(S=Office_setup$S_x, I = Office_setup$I_x, R = Office_setup$R_x, P = Office_setup$P_x)
Office_graph <- graph_from_adjacency_matrix(Office_adjacency_matrix,mode= "undirected") %>% 
  set_vertex_attr( "Infectious", value = Office_setup$I_x) %>% 
  set_vertex_attr( "Room", value = seq(1:N_rooms))

graph_to_use <- Office_graph

Office_network_fig <- ggplot(graph_to_use, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = Infectious), color = "#C21807")+
  geom_nodetext_repel(aes(label=Room))+
  theme_blank()+
  scale_size(name = "Number of\nInfectious\npeople in\nrooms")+
  labs(title = "Network representation of an Office")+
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust=1))
Office_network_fig
# png(filename = "Office_network_fig.png",units="in", width=6, height=4, res=300)
# Office_network_fig
# dev.off()


#Initialize Transition matrices, both people and particles
Office_T_mov <- Create_Particle_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use, prop = 0.9)
Office_theta_mov <- Create_Particle_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use, prop =0.5)

Office_end_day_T_mov <- End_of_day_T_function(adjacency_matrix_to_use = adjacency_matrix_to_use)
End_of_day_T_mov_to_use <- Office_end_day_T_mov

#parms <- data.frame(s=10,a=0.000005, d=.06)

Office_output <-  data.frame(lsoda(y = Office_Init_conds, func = Particle_model,times = times, parms = parms,adjacency_matrix_to_use=adjacency_matrix_to_use,theta_mov = Office_theta_mov,T_mov = Office_T_mov, C = Office_C))

Office_data_clean <- Office_output%>% pivot_longer(cols = !time,
                                                   names_to = c("State", "Room"),
                                                   names_pattern = "([A-Za-z]+)(\\d+)",
                                                   values_to = "Number")

Office_data_ratios <- Office_data_clean %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% group_by(time,Room)  %>% 
  mutate(N_x = S+I+R) %>% group_by(time,Room) %>% 
  mutate(K_x = (parms$s*I)/(parms$a*N_x), prop_to_K = P/K_x)
Office_P_x_to_K_x_plot <- ggplot(Office_data_ratios,aes(x =time, y =prop_to_K,group= Room,color=Room))+geom_line()+labs(title = "Office - Risk - 37 rooms")

Office_P_x_plot <-ggplot(Office_data_ratios, aes(x=time, y = P, color = Room))+geom_line()+labs(title = "Office - Risk - 37 rooms")
Office_P_x_to_K_x_plot

#### Old Office figure making ####
# Office_Particle_data<- Build_data_mod %>% filter(State == "P") %>% group_by(Room,time) %>% summarise(Room = as.numeric(Room),Room_Pop =sum(Number))
# Office_Particle_data$Room<-as.factor(Office_Particle_data$Room)
# 
# label_data <-  Office_Particle_data %>%
#   distinct(Room, .keep_all = TRUE)
# 
# Office_model_output_fig<-ggplot(Office_Particle_data,aes(x=time,y=Room_Pop,group=Room,color=as.numeric(Room)))+geom_line()+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   theme_linedraw()+
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# 
# Office_model_output_fig
# 
# # png(filename = "Office_model_output.png",units="in", width=6, height=4, res=300)
# # Office_model_output_fig
# # dev.off()
# 
# 
# ############### plotting particles
# ggplot(Office_Particle_data, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
#   geom_line() +
#   geom_text(data = Office_Particle_data %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   labs(title = "Particle accumulation over time", x = "Time", y = "Number of infectious Particles", color = "Rooms") +
#   theme_minimal() +
#   theme(legend.position = "right")
# 
# ################## plotting people
# Office_People_data<- Build_data_mod %>% filter(State != "P") %>% group_by(Room,time) %>% summarise(Room = as.numeric(Room),Room_Pop =sum(Number))
# 
# Office_People_data$Room<-as.factor(Office_People_data$Room)
# 
# label_data <-  Office_People_data %>%
#   distinct(Room, .keep_all = TRUE)
# 
# ggplot(Office_People_data, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
#   geom_line() +
#   geom_text(data = Office_People_data %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   labs(title = "People in rooms over time", x = "Time", y = "Number of people in rooms", color = "Rooms") +
#   theme_minimal() +
#   theme(legend.position = "right")+ylim(0,100)
# ################# Facet_wrap - People
# 
# ggplot(Office_People_data, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
#   geom_line() +
#   facet_wrap(~Room)+
#   #geom_text(data = Office_People_data %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   labs(title = "People in rooms over time", x = "Time", y = "Number of people in rooms", color = "Rooms") +
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))+
#   theme_linedraw()
# 
# ################# Facet_wrap - Particles
# 
# ggplot(Office_Particle_data,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   facet_wrap(~Room)+
#   theme_linedraw()+
#   #geom_text_repel(data=Office_Particle_data,aes(label = Room), show.legend = FALSE)+
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# # calculating levels of risk throughout the building. Potential calculations are:
# # Number of infectious particles at the end of the day in each room
# # Average number of infectious particles throughout the day, in each room
# # Cumulative number of infectious particles in each room throughout the day
# # Average risk in the building - Total number of infectious particles in the whole building divided by the number of rooms OR
# # Number of infectious particles in each room divided by total in building then take the average of that quantity
# # Total number of infectious particles in each room divided by their saturation point - this
# 
# #Number of infectious particles in each room over time for the last day - time >=48
# N_Inf_part_throughout_day <- filter(Office_Particle_data, time >= 48)
# N_Inf_part_throughout_day$Room <- as.factor(N_Inf_part_throughout_day$Room)
# 
# # All on one
# ggplot(N_Inf_part_throughout_day,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   #facet_wrap(~Room)+
#   theme_linedraw()+
#   geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# #Each room with their own graph
# office_N_inf_part_by_room_throughout_day<- ggplot(N_Inf_part_throughout_day,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   facet_wrap(~Room)+
#   theme_linedraw()+
#   #geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "none",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# office_N_inf_part_by_room_throughout_day
# png(filename = "office_N_inf_part_by_room_throughout_day.png",units="in", width=10, height=10, res=300)
# office_N_inf_part_by_room_throughout_day
# dev.off()
# #Now taking into account the saturation point
# sat_points <- data.frame(Room = as.factor(1:N_rooms), Sat_point = K)
# N_Inf_part_rooms_sat_point <- right_join(N_Inf_part_throughout_day,sat_points, by = join_by(Room)) %>%
#   mutate(Prop_to_Sat_point = Room_Pop/Sat_point)
# office_Ratio_inf_part_by_room_throughout_day <-ggplot(N_Inf_part_rooms_sat_point, aes(x = time, y= Prop_to_Sat_point,group=Room,color= Room, label= Room))+geom_line()+
#   labs(x = "Time (hours)", y= "Ratio",title = "Ratio of infectious particles in room to the threshold number of particles of each room")+
#   facet_wrap(~Room)+
#   theme_linedraw()+
#   #geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "none",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# office_Ratio_inf_part_by_room_throughout_day
# png(filename = "office_Ratio_inf_part_by_room_throughout_day.png",units="in", width=10, height=10, res=300)
# office_Ratio_inf_part_by_room_throughout_day
# dev.off()
# 
# # Number of rooms over sat point by 10 percent
# N_rooms_over_sat_point_10 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.10)))
# N_rooms_over_sat_point_10 <- N_rooms_over_sat_point_10 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_10, aes(x=time, y =n))+ geom_line(color = "turquoise")
# 
# # Number of rooms over sat point by 30 percent
# N_rooms_over_sat_point_30 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.30)))
# N_rooms_over_sat_point_30 <- N_rooms_over_sat_point_30 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_30, aes(x=time, y =n))+ geom_line( color = "purple")
# 
# # Number of rooms over sat point by 50 percent
# N_rooms_over_sat_point_50 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.50)))
# N_rooms_over_sat_point_50 <- N_rooms_over_sat_point_50 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_50, aes(x=time, y =n))+ geom_line(color = "pink")
# 
# # Number of rooms over sat point by 70 percent
# N_rooms_over_sat_point_70 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.70)))
# N_rooms_over_sat_point_70 <- N_rooms_over_sat_point_70 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_70, aes(x=time, y =n))+ geom_line(color  = "orange")
# 
# # Number of rooms over sat point by 90 percent
# N_rooms_over_sat_point_90 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.90)))
# N_rooms_over_sat_point_90 <- N_rooms_over_sat_point_90 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_90, aes(x=time, y =n))+ geom_line(color = "green4")
# 
# # Number of rooms over sat point by 150 percent
# N_rooms_over_sat_point_150 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*1.5)))
# N_rooms_over_sat_point_150 <- N_rooms_over_sat_point_150 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_150, aes(x=time, y =n))+ geom_line(color = "red4")
# 
# # Number of rooms over sat point by 200 percent
# N_rooms_over_sat_point_200 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*2)))
# N_rooms_over_sat_point_200 <- N_rooms_over_sat_point_200 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_200, aes(x=time, y =n))+ geom_line(color= "blue4")
# 
# data_to_plot <- as.data.frame(cbind(time =N_rooms_over_sat_point_10$time,
#                                     Over_10 = N_rooms_over_sat_point_10$n,
#                                     Over_30 = N_rooms_over_sat_point_30$n,
#                                     Over_50 = N_rooms_over_sat_point_50$n,
#                                     Over_70 = N_rooms_over_sat_point_70$n,
#                                     Over_90 = N_rooms_over_sat_point_90$n,
#                                     Over_150 = N_rooms_over_sat_point_150$n,
#                                     Over_200 = N_rooms_over_sat_point_200$n))
# data_to_plot <- pivot_longer(data_to_plot, cols= starts_with("Over"),names_to = "Threshold", values_to = "Number")
# data_to_plot_V2 <- data_to_plot %>% mutate(Percent_over = (Number/N_rooms)*100)
# neworder <- c("Over_10","Over_30","Over_50", "Over_70", "Over_90","Over_150", "Over_200")
# data_to_plot_V2 <- arrange(transform(data_to_plot_V2,
#                                      Threshold=factor(Threshold,levels=neworder)),Threshold)
# 
# office_Percent_above_sat <- ggplot(data_to_plot_V2, aes(x= time, y = Percent_over, color = Threshold))+geom_line()+
#   facet_wrap(~Threshold)+
#   theme_linedraw()+
#   labs(x = "Time (Hours)", y = "Percent of\nrooms over\nsaturation point", title = "Percentage of rooms over saturation points at different thresholds")+
#   #ylim(0,50)+
#   scale_color_manual(name = "Percent over\nthreshold",
#                      labels = c("10%","30%","50%","30%","50%","70%","90%","150%","200%"),
#                      values = c("blue4","red4","green4","orange","pink","purple","turquoise"))+
#   #geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# 
# 
# png(filename = "office_Percent_above_sat.png",units="in", width=10, height=10, res=300)
# office_Percent_above_sat
# dev.off()
# 
# 
# # starting with just the number of particles in each room at the beginning of the last day.
# N_Inf_part_rooms_beg <- filter(Office_Particle_data, time == 48)
# N_Inf_part_rooms_beg$Room <- as.numeric(N_Inf_part_rooms_beg$Room)
# N_Inf_part_rooms_beg <- N_Inf_part_rooms_beg%>% arrange(.by_group = TRUE)
# ggplot(N_Inf_part_rooms_beg, aes(x = Room, y = Room_Pop))+ geom_col()+
#   labs(x="Room",y="Number of infectious particles", title = "Number of infectious particles in each room at the end of the day")
# 
# 
# #Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# room_particle_graph <-Office_graph %>% set_vertex_attr( "Infectious", value =N_Inf_part_rooms_beg$Room_Pop ) %>% 
#   set_vertex_attr("Room", value = seq(1:N_rooms))# + 
# Church_Init_conds <-c(S=Church_setup$S_x, I = Church_setup$I_x, R = Church_setup$R_x, P = Church_setup$P_x)
# 
# Office_Network_With_N_inf <- ggplot(room_particle_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(color = Infectious),size = 2)+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_color_continuous(name = "Number of\nInfectious\nParticles")+
#   labs(title = "Number of infectious particles in each room at the begininning of the day")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))
# Office_Network_With_N_inf
# 
# png(filename = "Office_Network_With_N_inf.png",units="in", width=10, height=10, res=300)
# Office_Network_With_N_inf
# dev.off()
# # Proportion of infectious particles at the end of the last day.
# N_Inf_part_rooms_end <- filter(Particles_in_rooms_data, time == 56) %>% ungroup() %>%mutate(Total = sum(Room_Pop),Prop=Room_Pop/sum(Room_Pop))
# N_Inf_part_rooms_end$Room <- as.factor(N_Inf_part_rooms_end$Room)
# ggplot(N_Inf_part_rooms_end, aes(x = Room, y = Prop))+ geom_col()+
#   labs(x="Room",y="Proportion of infectious particles", title = "Proportion of infectious particles in each room at the end of the day")
# 
# 
# #Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# room_particle_graph <-church_graph %>% set_vertex_attr( "Infectious", value =N_Inf_part_rooms_end$Prop ) %>% 
#   set_vertex_attr("Room", value = seq(1:N_rooms))# + 
# 
# ggplot(room_particle_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(color = Infectious))+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_color_continuous(name = "Number of\nInfectous\nParticles")+
#   labs(title = "Raw number of infectious particles in each room")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))
# 
# # Average number of particles throughout the day (but on the last day)
# N_Inf_part_rooms_day_avg <- filter(Particles_in_rooms_data, time >= 48 & time <= 56)%>% summarise(Total_inf_part_room=sum(Room_Pop)) %>% mutate(Day_avg_inf_part = Total_inf_part_room/8)
# #%>% ungroup() %>%mutate(Total = sum(Room_Pop),Prop=Room_Pop/sum(Room_Pop))
# N_Inf_part_rooms_day_avg$Room <- as.factor(N_Inf_part_rooms_day_avg$Room)
# ggplot(N_Inf_part_rooms_end, aes(x = Room, y = Prop))+ geom_col()+
#   labs(x="Room",y="Proportion of infectious particles", title = "Proportion of infectious particles in each room at the end of the day")
# 
# 
# #Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# room_particle_graph <-church_graph %>% set_vertex_attr( "Infectious", value =N_Inf_part_rooms_day_avg$Prop ) %>% 
#   set_vertex_attr("Room", value = seq(1:N_rooms))# + 
# 
# ggplot(room_particle_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(color = Infectious))+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_color_continuous(name = "Number of\nInfectous\nParticles")+
#   labs(title = "Raw number of infectious particles in each room")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))


####{ Movie theater} ####
# movie_theater_adjacency_matrix <- matrix(c(c(0,1,1,rep(0,3),1,0,0,1,rep(0,22)), # Column 1
#                                            c(1,0,rep(1,3),rep(0,5),rep(1,15),rep(0,7)), # Column 2
#                                            c(1,1,rep(0,11),1,rep(0,18)), # Column 3
#                                            c(0,1,rep(0,2),1,rep(0,16),rep(1,7),rep(0,4)), # Column 4
#                                            c(rep(0,5),1, rep(0,21),1,rep(0,4)), # Column 5
#                                            c(rep(0,4),1,rep(0,24),1,1,1), # Column 6
#                                            c(1,rep(0,6),1,0,0,1,rep(0,21)), # Column 7
#                                            c(rep(0,6),1,0,1,rep(0,23)), #Column 8
#                                            c(rep(0,7),1,rep(0,24)), # Column 9
#                                            c(1,1,1,rep(0,29)), # Column 10
#                                            c(0,1,rep(0,4),1,rep(0,25)), # Column 11
#                                            rep(c(0,1,rep(0,30)),2), #Columns 12 and 13
#                                            c(0,0,1,rep(0,29)), # Column 14
#                                            c(1,1,1,rep(0,29)), #Column 15
#                                            rep(c(0,1,rep(0,30)),4), # Columns 16-19
#                                            rep(c(0,1,0,1,1,rep(0,27)),4), # Columns 20-23 
#                                            c(0,1,rep(0,22),1, rep(0,7)), # Column 24
#                                            c(rep(0,23),1,rep(0,8)), # Column 25
#                                            rep(c(0,0,0,1,rep(0,28)),3), # Columns 26-28
#                                            c(rep(0,4),1,rep(0,27)), # Column 29 
#                                            c(rep(0,5),1,rep(0,26)), # Column 30
#                                            c(rep(0,5),1, rep(0,25),1), # Column 31
#                                            c(rep(0,30),1,0) # Column 32
# ), nrow=32, ncol = 32)

movie_theater_adjacency_matrix <- matrix(c(c(0,1,1,rep(0,5),1,rep(0,5),1,0,1,rep(0,35-17)), # Room 1 hallway
                                           c(1,0,1,1,1,rep(0,6),rep(1,3),rep(0,3),rep(1,5),0,rep(1,4),rep(0,35-27)), # Room 2 - Hallway
                                           c(1,rep(0,13),1,1,1,rep(0,35-17)), # Room 3 - Hallway
                                           c(0,1,0,0,1,0,1,rep(0,18),rep(1,3),rep(0,35-28)), # Room 4 - Hallway
                                           c(0,1,0,1,0,1,rep(0,19),1,rep(0,4),1,rep(0,35-31)), # Room 5 - Hallway
                                           c(rep(0,4),1,rep(0,25),1,rep(0,35-31)), # Room 6 - Hallway
                                           c(rep(0,3),1,rep(0,3),1,rep(0,15),1,rep(0,3),1,rep(0,35-28)), # Room 7 - Hallway
                                           c(rep(0,6),1,rep(0,16),1,rep(0,5),1,rep(0,35-30)), # Room 8 - Hallway
                                           c(1, rep(0,8),1,rep(0,3),1,rep(0, 35-14)), # Room 9 main room
                                           c(rep(0,8),1,0,1,rep(0, 35-11)), # Room 10 Admin
                                           c(rep(0,9),1,rep(0, 35-10)), # Room 11 - admin
                                           rep(c(0,1, rep(0,35-2)),2), # Room 12 and 13 - facilities management
                                           c(0,1,rep(0,6),1,rep(0,34-9),1), # Room 14
                                           c(1,rep(0, 35-1)), # Room 15 - Theater
                                           c(0,0,1,rep(0,35-3)), # Room 16 - Storage room
                                           c(1,0,1,rep(0,35-3)), # Room 17
                                           rep(c(0,1,rep(0,35-2)),4), # Rooms 18-21 restrooms and theater (19)
                                           c(0,1,rep(0,20),1,rep(0,35-23)), # Room 22 - cleaning
                                           c(rep(0,21),1, rep(0,35-22)), # Room 23 - cleaning
                                           c(0,1,rep(0,4),1,1,rep(0,20),1,rep(0,35-29)), # Room 24 - Theater
                                           c(0,1,rep(0,35-2)), # Room 25 - facilities
                                           c(0,1,0,1,1,rep(0,35-5)), # Room 26 - Theater
                                           c(0,1,0,1,rep(0,35-4)), # Room 27 - Theater
                                           c(rep(0,3),1,0,0,1,rep(0,35-7)), # Room 28 - Theater
                                           c(rep(0,23),1,rep(0,35-24)), # Room 29 Facilities
                                           c(rep(0,7),1,rep(0,35-8)), # Room 30 - storage
                                           c(rep(0,4),1,rep(0,35-5)), # Room 31 - Theater
                                           c(rep(0,5),1,rep(0,26),1,rep(0,35-33)), # Room 32 - maintenance
                                           c(rep(0,31),1,0,1,0), # Room 33 - Employee room
                                           c(rep(0,32),1,0,0), # Room 34 - employee bathroom
                                           c(rep(0,8),1,rep(0,35-9))), # Room 35 - Outside
                                         
                                         nrow = 35,ncol = 35)

movie_theater_graph<- graph_from_adjacency_matrix(movie_theater_adjacency_matrix, mode = "undirected")

graph_to_use <- movie_theater_graph
adjacency_matrix_to_use <- movie_theater_adjacency_matrix
N_rooms <- ncol(adjacency_matrix_to_use) #number of rooms
Movie_C <- c(
  400, # Room 1 hallway
  300, # Room 2 - Hallway
  100, # Room 3 - Hallway 
  100, # Room 4 - Hallway
  100, # Room 5 - Hallway
  100, # Room 6 - Hallway 
  100, # Room 7 - Hallway 
  100, # Room 8 - Hallway
  500, # Room 9 main room
  6, # Room 10 Admin
  4, # Room 11 - admin
  3, # Room 12 - facilities management
  3, # Room 13 - facilities management
  25, # Room 14 - Kitchen
  100, # Room 15 - Theater
  3, # Room 16 - Storage room
  100, # Room 17 - Theater
  30, # Room 18 - restroom
  100, # Room 19 - Theater
  30, # Room 20 - restroom
  2, # Room 21 - family bathroom
  3, # Room 22 - cleaning
  3, # Room 23 - cleaning
  100, # Room 24 - Theater
  3, # Room 25 - facilities
  100, # Room 26 - Theater
  100, # Room 27 - Theater
  100, # Room 28 - Theater
  3, # Room 29 Facilities
  3, # Room 30 - storage
  100, # Room 31 - Theater
  3, # Room 32 - maintenance
  10, # Room 33 - Employee room
  2 # Room 34 - employee bathroom
  
)
Movie_C <- c(Movie_C,sum(Movie_C))
#outside should only be able to hold the total capacity of the building 
Building_max <- Movie_C[N_rooms] #last room in C is the "outside" room
Max_capacity <- Prop_full*Building_max
delt <- Max_capacity/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

movie_theater_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)

movie_theater_Init_conds <-c(S=movie_theater_setup$S_x, I = movie_theater_setup$I_x, R = movie_theater_setup$R_x, P = movie_theater_setup$P_x)

movie_theater_graph<-graph_from_adjacency_matrix(movie_theater_adjacency_matrix,mode= "undirected")%>% 
  set_vertex_attr( "Infectious", value = movie_theater_setup$I_x) %>% 
  set_vertex_attr( "Room", value = seq(1:N_rooms))
# graph_to_use <- movie_theater_graph
# layout <- layout_with_fr(movie_theater_graph)
# plot(movie_theater_graph, layout = layout)

movie_theater_network_fig <- ggplot(movie_theater_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = Infectious), color = "#C21807")+
  geom_nodetext_repel(aes(label=Room))+
  theme_blank()+
  scale_size(name = "Number of\nInfectious\npeople in\nrooms")+
  labs(title = "Network representation of a movie theater")+
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust=1))
movie_theater_network_fig
#png(filename = "movie_theater_network_fig.png",units="in", width=6, height=4, res=300)
#movie_theater_network_fig
#dev.off()

#create transition matrices for people and particles
#T_mov is the transition matrix for individuals moving between rooms
movie_theater_T_mov <- Create_Particle_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use,prop = 0.9)
movie_theater_theta_mov <- Create_Particle_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use,prop =0.5)
movie_theater_end_day_T_mov <- End_of_day_T_function(adjacency_matrix_to_use = adjacency_matrix_to_use)
End_of_day_T_mov_to_use <- movie_theater_end_day_T_mov




movie_theater_output <- data.frame(lsoda(y = movie_theater_Init_conds, func = Particle_model,times = times, parms = parms,adjacency_matrix_to_use= adjacency_matrix_to_use,theta_mov =movie_theater_theta_mov,T_mov=movie_theater_T_mov, C=Movie_C))


Movie_data_clean <- movie_theater_output%>% pivot_longer(cols = !time,
                                                         names_to = c("State", "Room"),
                                                         names_pattern = "([A-Za-z]+)(\\d+)",
                                                         values_to = "Number")


Movie_data_ratios <- Movie_data_clean %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% group_by(time,Room)  %>% 
  mutate(N_x = S+I+R) %>% group_by(time,Room) %>% 
  mutate(K_x = (parms$s*I)/(parms$a*N_x), prop_to_K = P/K_x)
Movie_P_x_to_K_x_plot <- ggplot(Movie_data_ratios,aes(x =time, y =prop_to_K,group= Room,color=Room))+geom_line()+labs(title = "Movie - Risk - 35 rooms")

Movie_P_x_plot <- ggplot(Movie_data_ratios, aes(x=time, y = P, color = Room))+geom_line()+labs(title = "Movie - Risk - 35 rooms")
Movie_P_x_to_K_x_plot

#### Old Movie theater figure making ####
# Movie_Particle_data<- Build_data_mod %>% filter(State == "P") %>% group_by(Room,time) %>% summarise(Room = as.numeric(Room),Room_Pop =sum(Number))
# Movie_Particle_data$Room<-as.factor(Movie_Particle_data$Room)
# label_data <-  Movie_Particle_data %>%
#   distinct(Room, .keep_all = TRUE)
# 
# movie_theater_model_output_fig<-ggplot(Movie_Particle_data,aes(x=time,y=Room_Pop,group=Room,color=as.numeric(Room)))+geom_line()+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   theme_linedraw()+
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# 
# movie_theater_model_output_fig
# 
# # png(filename = "movie_theater_model_output_fig",units="in", width=6, height=4, res=300)
# # movie_theater_model_output_fig
# # dev.off()
# 
# 
# ############### plotting particles
# ggplot(Movie_Particle_data, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
#   geom_line() +
#   geom_text(data = Movie_Particle_data %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   labs(title = "Particle accumulation over time", x = "Time", y = "Number of infectious Particles", color = "Rooms") +
#   theme_minimal() +
#   theme(legend.position = "right")
# 
# ################## plotting people
# Movie_People_data<- Build_data_mod %>% filter(State != "P") %>% group_by(Room,time) %>% summarise(Room = as.numeric(Room),Room_Pop =sum(Number))
# 
# Movie_People_data$Room<-as.factor(Movie_People_data$Room)
# 
# label_data <-  Movie_People_data %>%
#   distinct(Room, .keep_all = TRUE)
# 
# ggplot(Movie_People_data, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
#   geom_line() +
#   geom_text(data = Movie_People_data %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   labs(title = "People in rooms over time", x = "Time", y = "Number of people in rooms", color = "Rooms") +
#   theme_minimal() +
#   theme(legend.position = "right")+ylim(0,100)
# ################# Facet_wrap - People
# 
# ggplot(Movie_People_data, aes(x = time, y = Room_Pop, color = Room, group = Room)) +
#   geom_line() +
#   facet_wrap(~Room)+
#   #geom_text(data = Movie_People_data %>% filter(time == last(time)), aes(label = Room, x = time +2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   labs(title = "People in rooms over time", x = "Time", y = "Number of people in rooms", color = "Rooms") +
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))+
#   theme_linedraw()
# 
# ################# Facet_wrap - Particles
# 
# ggplot(Movie_Particle_data,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   facet_wrap(~Room)+
#   theme_linedraw()+
#   #geom_text_repel(data=Movie_Particle_data,aes(label = Room), show.legend = FALSE)+
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# 
# 
# #Number of infectious particles in each room over time for the last day - time >=48
# N_Inf_part_throughout_day <- filter(Movie_Particle_data, time >= 48)
# N_Inf_part_throughout_day$Room <- as.factor(N_Inf_part_throughout_day$Room)
# 
# # All on one
# ggplot(N_Inf_part_throughout_day,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   #facet_wrap(~Room)+
#   theme_linedraw()+
#   geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# #Each room with their own graph
# movie_N_inf_part_by_room_throughout_day<- ggplot(N_Inf_part_throughout_day,aes(x=time,y=Room_Pop,group=Room,color= Room, label= Room))+geom_line()+ #+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   facet_wrap(~Room)+
#   theme_linedraw()+
#   #geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "none",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# movie_N_inf_part_by_room_throughout_day
# png(filename = "movie_N_inf_part_by_room_throughout_day.png",units="in", width=10, height=10, res=300)
# movie_N_inf_part_by_room_throughout_day
# dev.off()
# #Now taking into account the saturation point
# sat_points <- data.frame(Room = as.factor(1:N_rooms), Sat_point = K)
# N_Inf_part_rooms_sat_point <- right_join(N_Inf_part_throughout_day,sat_points, by = join_by(Room)) %>%
#   mutate(Prop_to_Sat_point = Room_Pop/Sat_point)
# movie_Ratio_inf_part_by_room_throughout_day <-ggplot(N_Inf_part_rooms_sat_point, aes(x = time, y= Prop_to_Sat_point,group=Room,color= Room, label= Room))+geom_line()+
#   labs(x = "Time (hours)", y= "Ratio",title = "Ratio of infectious particles in room to the threshold number of particles of each room")+
#   facet_wrap(~Room)+
#   theme_linedraw()+
#   #geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "none",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# movie_Ratio_inf_part_by_room_throughout_day
# png(filename = "movie_Ratio_inf_part_by_room_throughout_day.png",units="in", width=10, height=10, res=300)
# movie_Ratio_inf_part_by_room_throughout_day
# dev.off()
# 
# # Number of rooms over sat point by 10 percent
# N_rooms_over_sat_point_10 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.10)))
# N_rooms_over_sat_point_10 <- N_rooms_over_sat_point_10 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_10, aes(x=time, y =n))+ geom_line(color = "turquoise")
# 
# # Number of rooms over sat point by 30 percent
# N_rooms_over_sat_point_30 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.30)))
# N_rooms_over_sat_point_30 <- N_rooms_over_sat_point_30 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_30, aes(x=time, y =n))+ geom_line( color = "purple")
# 
# # Number of rooms over sat point by 50 percent
# N_rooms_over_sat_point_50 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.50)))
# N_rooms_over_sat_point_50 <- N_rooms_over_sat_point_50 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_50, aes(x=time, y =n))+ geom_line(color = "pink")
# 
# # Number of rooms over sat point by 70 percent
# N_rooms_over_sat_point_70 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.70)))
# N_rooms_over_sat_point_70 <- N_rooms_over_sat_point_70 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_70, aes(x=time, y =n))+ geom_line(color  = "orange")
# 
# # Number of rooms over sat point by 90 percent
# N_rooms_over_sat_point_90 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*.90)))
# N_rooms_over_sat_point_90 <- N_rooms_over_sat_point_90 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_90, aes(x=time, y =n))+ geom_line(color = "green4")
# 
# # Number of rooms over sat point by 150 percent
# N_rooms_over_sat_point_150 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*1.5)))
# N_rooms_over_sat_point_150 <- N_rooms_over_sat_point_150 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_150, aes(x=time, y =n))+ geom_line(color = "red4")
# 
# # Number of rooms over sat point by 200 percent
# N_rooms_over_sat_point_200 <- N_Inf_part_rooms_sat_point %>% group_by(time) %>%   count(N_over =Room_Pop >= (Sat_point+(Sat_point*2)))
# N_rooms_over_sat_point_200 <- N_rooms_over_sat_point_200 %>% filter(N_over == TRUE)
# ggplot(N_rooms_over_sat_point_200, aes(x=time, y =n))+ geom_line(color= "blue4")
# 
# data_to_plot <- as.data.frame(cbind(time =N_rooms_over_sat_point_10$time,
#                                     Over_10 = N_rooms_over_sat_point_10$n,
#                                     Over_30 = N_rooms_over_sat_point_30$n,
#                                     Over_50 = N_rooms_over_sat_point_50$n,
#                                     Over_70 = N_rooms_over_sat_point_70$n,
#                                     Over_90 = N_rooms_over_sat_point_90$n,
#                                     Over_150 = N_rooms_over_sat_point_150$n,
#                                     Over_200 = N_rooms_over_sat_point_200$n))
# data_to_plot <- pivot_longer(data_to_plot, cols= starts_with("Over"),names_to = "Threshold", values_to = "Number")
# data_to_plot_V2 <- data_to_plot %>% mutate(Percent_over = (Number/N_rooms)*100)
# neworder <- c("Over_10","Over_30","Over_50", "Over_70", "Over_90","Over_150", "Over_200")
# data_to_plot_V2 <- arrange(transform(data_to_plot_V2,
#                                      Threshold=factor(Threshold,levels=neworder)),Threshold)
# 
# movie_Percent_above_sat <- ggplot(data_to_plot_V2, aes(x= time, y = Percent_over, color = Threshold))+geom_line()+
#   facet_wrap(~Threshold)+
#   theme_linedraw()+
#   labs(x = "Time (Hours)", y = "Percent of\nrooms over\nsaturation point", title = "Percentage of rooms over saturation points at different thresholds")+
#   #ylim(0,50)+
#   scale_color_manual(name = "Percent over\nthreshold",
#                      labels = c("10%","30%","50%","30%","50%","70%","90%","150%","200%"),
#                      values = c("blue4","red4","green4","orange","pink","purple","turquoise"))+
#   #geom_text(data = N_Inf_part_throughout_day %>% filter(time == last(time)), aes(label = Room, x = time+2, y = Room_Pop,color = Room), show.legend = FALSE) +  # Adding labels with ggrepel
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# 
# movie_Percent_above_sat
# png(filename = "movie_Percent_above_sat.png",units="in", width=10, height=10, res=300)
# movie_Percent_above_sat
# dev.off()
# 
# 
# # starting with just the number of particles in each room at the beginning of the last day.
# N_Inf_part_rooms_beg <- filter(Movie_Particle_data, time == 48)
# N_Inf_part_rooms_beg$Room <- as.numeric(N_Inf_part_rooms_beg$Room)
# N_Inf_part_rooms_beg <- N_Inf_part_rooms_beg%>% arrange(.by_group = TRUE)
# ggplot(N_Inf_part_rooms_beg, aes(x = Room, y = Room_Pop))+ geom_col()+
#   labs(x="Room",y="Number of infectious particles", title = "Number of infectious particles in each room at the end of the day")
# 
# 
# #Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# room_particle_graph <-movie_theater_graph %>% set_vertex_attr( "Infectious", value =N_Inf_part_rooms_beg$Room_Pop ) %>% 
#   set_vertex_attr("Room", value = seq(1:N_rooms))# + 
# Church_Init_conds <-c(S=Church_setup$S_x, I = Church_setup$I_x, R = Church_setup$R_x, P = Church_setup$P_x)
# 
# movie_Network_With_N_inf <- ggplot(room_particle_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(color = Infectious),size=2)+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_color_continuous(name = "Number of\nInfectious\nParticles")+
#   labs(title = "Number of infectious particles in each room at the beginning of the day")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))
# movie_Network_With_N_inf
# 
# png(filename = "movie_Network_With_N_inf.png",units="in", width=10, height=10, res=300)
# movie_Network_With_N_inf
# dev.off()
# # Proportion of infectious particles at the end of the last day.
# N_Inf_part_rooms_end <- filter(Particles_in_rooms_data, time == 56) %>% ungroup() %>%mutate(Total = sum(Room_Pop),Prop=Room_Pop/sum(Room_Pop))
# N_Inf_part_rooms_end$Room <- as.factor(N_Inf_part_rooms_end$Room)
# ggplot(N_Inf_part_rooms_end, aes(x = Room, y = Prop))+ geom_col()+
#   labs(x="Room",y="Proportion of infectious particles", title = "Proportion of infectious particles in each room at the end of the day")
# 
# 
# #Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# room_particle_graph <-church_graph %>% set_vertex_attr( "Infectious", value =N_Inf_part_rooms_end$Prop ) %>% 
#   set_vertex_attr("Room", value = seq(1:N_rooms))# + 
# 
# ggplot(room_particle_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(color = Infectious))+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_color_continuous(name = "Number of\nInfectous\nParticles")+
#   labs(title = "Raw number of infectious particles in each room")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))
# 
# # Average number of particles throughout the day (but on the last day)
# N_Inf_part_rooms_day_avg <- filter(Particles_in_rooms_data, time >= 48 & time <= 56)%>% summarise(Total_inf_part_room=sum(Room_Pop)) %>% mutate(Day_avg_inf_part = Total_inf_part_room/8)
# #%>% ungroup() %>%mutate(Total = sum(Room_Pop),Prop=Room_Pop/sum(Room_Pop))
# N_Inf_part_rooms_day_avg$Room <- as.factor(N_Inf_part_rooms_day_avg$Room)
# ggplot(N_Inf_part_rooms_end, aes(x = Room, y = Prop))+ geom_col()+
#   labs(x="Room",y="Proportion of infectious particles", title = "Proportion of infectious particles in each room at the end of the day")
# 
# 
# #Church_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)
# room_particle_graph <-church_graph %>% set_vertex_attr( "Infectious", value =N_Inf_part_rooms_day_avg$Prop ) %>% 
#   set_vertex_attr("Room", value = seq(1:N_rooms))# + 
# 
# ggplot(room_particle_graph, aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey75")+
#   geom_nodes(aes(color = Infectious))+
#   geom_nodetext_repel(aes(label=Room))+
#   theme_blank()+
#   scale_color_continuous(name = "Number of\nInfectous\nParticles")+
#   labs(title = "Raw number of infectious particles in each room")+
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust=1))
# 



#### University Building ####
University_adjacency_matrix<- matrix(c(rep(c(rep(0,34),1,rep(0,224-35)),6),#Basement floor -rooms 1-6
                                       rep(c(rep(0,35),1,rep(0,224-36)),4), #rooms 7-10
                                       c(rep(0,35),1,0,0,1,1,rep(0,224-40)), #room 11
                                       c(rep(0,36),1,rep(0,224-37)), #room 12
                                       c(rep(0,13),1,rep(0,22),1,rep(0,224-37)), #room 13
                                       c(rep(0,12),1,0,1,0,1,rep(0,20),1, rep(0,224-38)), #room 14
                                       c(rep(0,13),1,0,1,0,1,1,rep(0,18),1,rep(0,224-38)), #room 15
                                       c(rep(0,14),1,rep(0,4),1,1, rep(0,16),1,rep(0,224-38)), #room 16
                                       c(rep(0,13),1,rep(0,224-14)), #room 17
                                       rep(c(rep(0,14),1,rep(0,224-15)),2), #rooms 18-19
                                       rep(c(rep(0,15),1,rep(0,224-16)),2), #rooms 20-21
                                       c(rep(0,39),1,rep(0,224-40)), #room 22
                                       c(rep(0,39),1,1,rep(0, 224-41)), #room 23
                                       rep(c(rep(0,40), 1, rep(0,224-41)),2), #rooms 24-25
                                       c(rep(0,26),1,rep(0,13),1,1,rep(0,224-42)), #room 26
                                       c(rep(0,25),1,rep(0,16),1,rep(0,224-43)), #room 27
                                       c(rep(0,29),1,rep(0,11),1,rep(0,224-42)), # room 28
                                       c(rep(0,30),1,rep(0,10),1,rep(0,224-42)), #room 29
                                       c(rep(0,27),1,rep(0,224-28)), #room 30  
                                       c(rep(0,28),1,rep(0,224-29)), #room 31
                                       c(rep(0,41),1,rep(0,224-42)), #room 32
                                       c(rep(0,41),1,1,rep(0,224-43)), #room 33
                                       c(rep(0,42),1,rep(0,224-43)), #room 34
                                       #BASEMENT LEVEL HALLWAYS
                                       c(rep(1,6),rep(0,29),1,rep(0,7),1,rep(0,224-44)), #room 35 - (hallway)
                                       c(rep(0,6),rep(1,5),rep(0,25),1,0,1,rep(0,224-39)), #room 36 - (hallway)
                                       c(rep(0,11),1,1,rep(0,22),1,0,1,rep(0,6),1,rep(0,224-45)), #room 37 - (hallway)
                                       c(rep(0,13),rep(1,3),rep(0,20),1,0,0,1,rep(0,224-40)), #room 38 - (hallway)
                                       c(rep(0,10),1,rep(0,24),1,0,0,0,1,rep(0,224-40)), #room 39 - (hallway)
                                       c(rep(0,10),1,rep(0,10),1,1,rep(0,14),1,1,0,1,0,0,0,0,1,rep(0,224-46)), #room 40 - hallway
                                       c(rep(0,22),rep(1,4),rep(0,13),1,0,1,rep(0,224-42)), #room 41 - (hallway)
                                       c(rep(0,25),1,0,1,1,0,0,1,1,rep(0,9),1, rep(0,224-43)), #room 42 - (hallway)
                                       c(rep(0,32),1,1,rep(0,7),1,rep(0,4),1,rep(0,224-47)), #room 43 - (hallway)
                                       #STAIRS AND ELEVATORS * will fill in as I go, want to check the network of each floor
                                       c(rep(0,34),1,rep(0,45),1,rep(0,59),1,rep(0,66),1,rep(0,224-208)), #room 44 - (stairwell 1)
                                       c(rep(0,36),1,rep(0,14),1,rep(0,90),1,0,1,rep(0,64),1,0,1,rep(0,224-212)), #room 45 - (elevator 1)
                                       c(rep(0,39),1,rep(0,46),1,rep(0,62),1,rep(0,66),1,rep(0,224-217)), #room 46 - (elevator 2)
                                       c(rep(0,42),1,rep(0,48),1,rep(0,61),1,rep(0,66),1,rep(0,224-221)),  #room 47 - (stairwell 2)
                                       ################# LEVEL 1 ##
                                       rep(c(rep(0,80),1,rep(0,224-81)),2), # rooms 48 and 49
                                       c(rep(0,80),1,0,1,rep(0,224-83)), # room 50
                                       c(rep(0,80),1,rep(0,224-81)), #room 51
                                       c(rep(0,44),1,rep(0,7),1,rep(0,11),1,rep(0,15),1,1,rep(0,223-82),1), # room 52
                                       c(rep(0,51),1,rep(0,224-52)), # room 53
                                       rep(c(rep(0,82),1,rep(0,224-83)),7), # rooms 54-60,
                                       c(rep(0,83),1,rep(0,224-84)), # room 61
                                       c(rep(0,66),1,rep(0,17),1,0,1,rep(0,224-87)), #room 62
                                       c(rep(0,84),1,rep(0,224-85)), # room 63
                                       c(rep(0,66),1,rep(0,17),1,rep(0,224-85)), # room 64
                                       c(rep(0,51),1,rep(0,32),1,rep(0,224-85)), # room 65
                                       c(rep(0,84),1,rep(0,224-85)), # room 66
                                       c(rep(0,61),1,0,1,rep(0,224-64)), #room 67
                                       rep(c(rep(0,86),1,rep(0,4),1,rep(0,224-92)),2), # room 68 and 69
                                       rep(c(rep(0,89),1,rep(0,224-90)),3), # room 70:72
                                       rep(c(rep(0,90),1,rep(0,224-91)),3), # room 73-75
                                       rep(c(rep(0,91),1,rep(0,224-92)),4), # room 76:79
                                       c(rep(0,85),1,0,1,0,0,1,rep(0,224-91)), # room 80
                                       c(rep(0,43),1,rep(0,3),rep(1,5),rep(0,29),1,rep(0,223-82),1), # room 81 - (hallway)
                                       c(rep(0,51),1,rep(0,28),1,0,1,1,1,rep(0,224-85)), # room 82 - (Hallway)
                                       c(rep(0,49),1,0,0,0,rep(1,7),rep(0,21),1,0,1,rep(0,224-84)), # room 83 - (Hallway)
                                       c(rep(0,60),1,rep(0,20),1,1,0,1,1,1,1,1,rep(0,224-89)), # room 84 - (Hallway)
                                       c(rep(0,61),rep(1,5),rep(0,15),1,0,1,rep(0,224-84)), # room 85 - (Hallway)
                                       c(rep(0,79),1,0,0,0,1,0,0,0,1,rep(0,223-88),1), # room 86 - (Hallway)
                                       c(rep(0,45),1,rep(0,15),1,rep(0,5),1,1,rep(0,14),1,rep(0,4),1,1,rep(0,223-90),1), # room 87 - (Hallway)
                                       c(rep(0,79),1,rep(0,3),1,0,1,rep(0,4),1,rep(0,224-91)), # room 88 - (Hallway)
                                       c(rep(0,83),1,0,0,1,rep(0,3),1,rep(0,224-91)), # room 89 - (Hallway)
                                       c(rep(0,69),rep(1,3),rep(0,14),1, rep(0,224-87)), # room 90 - Hallway
                                       c(rep(0,72),1,1,1,rep(0,4),1, rep(0,7),1,1,0,0,1,1, rep(0,223-93),1), # room 91 - (Hallway)
                                       c(rep(0,46),1,rep(0,20),1,1,rep(0,6),rep(1,4),rep(0,223-79),1), # room 92 - (Hallway)
                                       ########################## STAIR ON LEVEL 1 room - 93###
                                       c(rep(0,90),1,rep(0,58),1,0,0,1,rep(0,63),1,0,0,1,rep(0,224-220)), # room 93 ( Staircase)
                                       ########################## BEGIN FLOOR 2 ###
                                       rep(c(rep(0,140),1,rep(0,224-141)),3), # rooms 94:96 
                                       rep(c(rep(0,141),1,rep(0,224-142)),4), # rooms 97:100
                                       rep(c(rep(0,143),1,rep(0,224-144)),4), #rooms 101: 104
                                       c(rep(0,141),1,1,rep(0,224-143)), # room 105
                                       c(rep(0,142),1,rep(0,224-143)), #room 106
                                       c(rep(0,144),1,rep(0,224-145)), # room 107 
                                       c(rep(0,112),1,1,rep(0,30),1,1,1,rep(0,224-147)), # room 108
                                       c(rep(0,144),1,1,1,rep(0,224-147)), #room 109
                                       c(rep(0,144),1,1,rep(0,224-146)), # room 110
                                       c(rep(0,145),1,rep(0,224-146)), # room 111
                                       c(rep(0,146),1, rep(0,224-147)), # room 112
                                       rep(c(rep(0,107),1,rep(0,224-108)),2), # rooms 113 and 114
                                       rep(c(rep(0,147),1,rep(0,224-148)),2), # rooms 115 and 116
                                       c(rep(0,149),1,rep(0,224-150)), # room 117
                                       rep(c(rep(0,148),1,rep(0,224-149)),3), # room 118:120
                                       c(rep(0,121),1,rep(0,25),1,rep(0,6),1, rep(0,224-155)), # room 121
                                       c(rep(0,120),1,rep(0,32),1,rep(0,224-154)), # room 122
                                       rep(c(rep(0,154),1,rep(0,224-155)),3), # rooms 123:125
                                       rep(c(rep(0,152),1,rep(0,224-153)),2), # rooms 126 and 127 
                                       rep(c(rep(0,153),1,rep(0,224-154)),5), # rooms 128:132 
                                       rep(c(rep(0,152),1,1,rep(0,224-154)),2), # rooms 133 and 134
                                       c(rep(0,135),1,0,1,1,rep(0,11),1,rep(0,224-151)), # room 135
                                       c(rep(0,134),1,0,0,0,1,rep(0,224-139)), # room 136
                                       c(rep(0,137),1,1,rep(0,224-139)), # room 137
                                       c(rep(0,134),1,0,1,rep(0,224-137)), # room 138
                                       c(rep(0,134),1,1,1,rep(0,224-137)), # room 139
                                       c(rep(0,150),1,rep(0,224-151)), # room 140
                                       ############ start of level 2 Hallways ###
                                       c(rep(0,43),1,rep(0,49),1,1,1,rep(0,45),1,rep(0,224-142)), # room 141 - (Hallway)
                                       c(rep(0,96),rep(1,4),rep(0,4),1,rep(0,35),1,0,1,1,rep(0,224-144)), # room 142 - (Hallway)
                                       c(rep(0,44),1, rep(0,59),1,1,rep(0,35),1,0,0,1,rep(0,224-145)), # room 143 - (Hallway)
                                       c(rep(0,100),rep(1,4),rep(0,45),1,1,rep(0,224-151)), # room 144 - (Hallway)
                                       c(rep(0,44),1,rep(0,61),rep(1,4),rep(0,32),1,rep(0,224-143)), # room 145 - (Hallway)
                                       c(rep(0,107),rep(1,4),rep(0,224-111)), # room 146 - (Hallway)
                                       c(rep(0,107),1,1,0,0,1,rep(0,35),1,rep(0,224-148)), # room 147 - (Hallway)
                                       c(rep(0,114),1,1,rep(0,4),1,rep(0,25),1,0,0,1,rep(0,224-150)), # room 148
                                       c(rep(0,117), 1,1,1,rep(0,29),1,rep(0,224-150)), # room 149 - (Hallway)
                                       c(rep(0,45),1,rep(0,46),1,rep(0,23),1,rep(0,26),1,rep(0,4),1,0,1,1,rep(0,224-152)), # room 150 - (Hallway)
                                       c(rep(0,134),1,0,0,0,1,1,rep(0,3),1,rep(0,5),1,0,0,1,rep(0,224-153)), # room 151 - (Hallway)
                                       c(rep(0,149),1,0,0,1,rep(0,224-153)), # room 152 - (Hallway)
                                       c(rep(0,125),1,1,rep(0,5),1,1,rep(0,4),1,rep(0,11),1,1,rep(0,224-152)),# room 153 - (Hallway)
                                       c(rep(0,46),1,rep(0,74),1,rep(0,4),rep(1,8),rep(0,18),1,0,1,rep(0,224-155)), # room 154 - (Hallway)
                                       c(rep(0,120),1,0,1,1,1,rep(0,28),1,rep(0,224-154)), # room 155 - (Hallway),
                                       ##################### BEGIN LEVEL 3 ### !!!!!!!!!!!!
                                       rep(c(rep(0,207),1,rep(0,224-208)),2), # room 156 and 157
                                       c(rep(0,207),1,1,rep(0,224-209)), # room 158
                                       rep(c(rep(0,208),1,rep(0,224-209)),4), # rooms 159-162
                                       rep(c(rep(0,210),1,rep(0,224-211)),4), # rooms 163-166
                                       c(rep(0,208),1,1,rep(0,224-210)), # room 167
                                       c(rep(0,211),1,rep(0,224-212)), # room 168
                                       c(rep(0,173),1,1,rep(0,36),1,1,1,1,rep(0,224-215)), # room 169
                                       c(rep(0,211),1,1,1,rep(0,224-214)), # room 170
                                       c(rep(0,211),1,1,rep(0,224-213)), # room 171
                                       c(rep(0,212),1,rep(0,224-213)), # room 172
                                       c(rep(0,213),1,rep(0,224-214)), # room 173
                                       rep(c(rep(0,168),1,rep(0,224-169)),2), # rooms 174 and 175
                                       rep(c(rep(0,214),1,rep(0,224-215)),2), # rooms 176 and 177
                                       c(rep(0,214),1,rep(0,6),1,1,0), # room 178
                                       c(rep(0,222),1,0), # room 179
                                       rep(c(rep(0,221),1,0,0),3), # rooms 180-182
                                       c(rep(0,215),1,1,rep(0,224-217)), # room 183
                                       rep(c(rep(0,215),1,rep(0,224-216)),3), # rooms 184-186
                                       c(rep(0,219),1,rep(0,224-220)), # room 187
                                       c(rep(0,219),1,1,rep(0,224-221)), # room 188
                                       rep(c(rep(0,220),1,rep(0,224-221)),4), # room 189-192
                                       c(rep(0,220),1,rep(0,224-221)), # room 193
                                       rep(c(rep(0,219),1,1,rep(0,224-221)),2), # rooms 194 and 195
                                       c(rep(0,217),1,rep(0,224-218)), # room 196
                                       c(rep(0,197),1,1,1,rep(0,17),1,rep(0,224-218)), # room 197
                                       c(rep(0,196),1,rep(0,224-197)), # room 198 
                                       c(rep(0,196),1,0,0,rep(1,5),rep(0,224-204)), # room 199
                                       c(rep(0,196),1,0,1,rep(0,224-199)), # room 200
                                       rep(c(rep(0,198),1,rep(0,224-199)),4), # rooms 201-204
                                       rep(c(rep(0, 217),1,rep(0,224-218)),3), # rooms 205-207
                                       ######## LEVEL 3 HALLWAYS ###
                                       c(rep(0,43),1,rep(0,111),1,1,1,rep(0,50),1,rep(0,224-209)),# room 208 - (Hallway)
                                       c(rep(0,157),rep(1,5),rep(0,4),1,rep(0,40),1,0,1,1,rep(0,224-211)), # room 209 - (Hallway)
                                       c(rep(0,44),1,rep(0,121),1, rep(0,41),1,0,0,1,rep(0,224-212)), # room 210 - (Hallway)
                                       c(rep(0,162),rep(1,4),rep(0,42),1,rep(0,7),1,1,rep(0,224-218)), # room 211 - (Hallway)
                                       c(rep(0,44),1,rep(0,122),rep(1,4),rep(0,38),1,rep(0,224-210)), # room 212 - (Hallway)
                                       c(rep(0,168),rep(1,4),rep(0,224-172)), # room 213 - (Hallway)
                                       c(rep(0,168),1,1,0,0,1,rep(0,41),1,rep(0,224-215)), # room 214 - (Hallway)
                                       c(rep(0,168),1, rep(0,6),1,1,1,rep(0,35),1,0,0,1,rep(0,224-217)), # room 215 - (Hallway)
                                       c(rep(0,182),rep(1,4),rep(0,30),1,rep(0,224-217)), # room 216 - (Hallway)
                                       c(rep(0,45),1,rep(0,46),1,rep(0,89),1,rep(0,27),1,rep(0,3),1,1,0,0,1,rep(0,224-219)), # room 217 - (Hallway)
                                       c(rep(0,195),1,1,rep(0,7),rep(1,3),rep(0,3),1,rep(0,8),1,rep(0,224-220)), # room 218 - (Hallway)
                                       c(rep(0,216),1,0,0,1,rep(0,224-220)), # room 219 - (Hallway)
                                       c(rep(0,92),1,rep(0,93),1,1,rep(0,5),1,1,0,0,0,1,rep(0,18),1,1,0,1,rep(0,224-221)), # room 220 - (Hallway)
                                       c(rep(0,46),1,rep(0,140),rep(1,8),rep(0,24),1,0,1,1,0), # room 221 - (Hallway)
                                       c(rep(0,177),1,0,1,1,1,rep(0,38),1,rep(0,224-221)), # room 222 - (Hallway)
                                       c(rep(0,177),1,1,rep(0,41),1,rep(0,224-221)), # room 223 - (Hallway)
                                       c(rep(0,51),1,rep(0,28),1,rep(0,4),1,1,0,0,0,1,1,rep(0,224-92)) # room 224 - Outside
                                       
),nrow = 224,ncol = 224)
adjacency_matrix_to_use <- University_adjacency_matrix
N_rooms <- ncol(adjacency_matrix_to_use) #number of rooms

University_C <- c(rep(3,5),# rooms 1-5 storage / facilities
                  100, # room 6 research lab
                  rep(4,4), # rooms 7-10 offices
                  100, # room 11 - building services / equipment
                  3, # room 12 - storage
                  rep(100,4), # rooms 13-16 research labs
                  rep(20,5), # rooms 17-21 research lab storage
                  3, # room 22 storage
                  40, # room 23 lab space
                  3, # room 24 storage
                  rep(3,3), # rooms 25-27 - facilities
                  rep(25,2), # rooms 28 and 29 - bathrooms
                  rep(3,2), # rooms 30 and 31 - bathroom storage
                  200, # room 32 large unused space
                  3, # room 33 facilities
                  10, # room 34 - large facilities
                  rep(100,2), # rooms 35-36 - hallways
                  50, # room 37 - small hallway outside elevator
                  rep(100,5), # rooms 38-42 hallways
                  50, # room 43 - small hallway by stairs
                  100, # room 44 stair
                  10, # room 45 elevator
                  10, # room 46 elevator
                  100, # room 47 stair
                  3, # room 48 storage
                  40, # room 49 Teaching lab small room
                  100, # room 50 connecting room to offices/ small common room
                  10, # room 51 storage - teaching lab
                  100, # room 52 Entryway
                  3, # room 53 Building services (elevator?)
                  rep(4,7), # rooms 54-60 offices
                  25, # room 61 - small teaching lab
                  50, # room 62 large teaching lab
                  25, # room 63 small teaching lab
                  50, # room 64 large teaching lab
                  35, # room 65 medium teaching lab
                  rep(5,2), # rooms 66 and 67 teaching lab storage
                  rep(60,2), # rooms 68 and 69 large teaching rooms
                  rep(3,3), # rooms 70-72 storage/ building services
                  rep(8,2), # rooms 73 and 74 - bathrooms
                  3, # room 75 storage
                  rep(30,4), # rooms 76-79 small classrooms
                  200, # room 80 large lecture hall
                  rep(20,4), # rooms 81-84
                  10, # room 85 smaller hallway to teaching labs
                  rep(20,4), # rooms 86-89
                  10, # room 90 small hallway to storage rooms
                  200, # room 91 hallway with common area 
                  20, # room 92 hallway
                  48, # room 93 - stair way
                  3, # room 94 storage
                  4, # room 95 office
                  rep(10,9), # rooms 96-104 offices
                  rep(30,2), # rooms 105 and 106 small classroom
                  10, # room 107 small research room
                  100, # room 108 research lab
                  75, # room 109 research lab / desks
                  rep(10,6), # rooms 110-115 research storage / side rooms
                  rep(3,2), # rooms 116 and 117 - small storage 
                  rep(5,3), # rooms 118-120 facilities
                  100, # rooms 121 large research lab
                  rep(50,2), # rooms 122 and 123 small research labs
                  rep(20,2), # rooms 124 and 125 research storage rooms
                  rep(4,7), # rooms 126-132 offices
                  rep(8,2), # rooms 133 and 134 bathrooms
                  75, # room 135 small common area
                  rep(4,3), # rooms 136-138 small offices
                  50, # room 139 large conference room
                  3, # room 140 facilities
                  20, # rooms 141 hallway
                  50, # room 142 hallway
                  20, # room 143 hallway
                  50, # room 144 hallway
                  rep(20,5), # rooms 145-149 hallways 
                  rep(75,3), # rooms 150-152 hallways
                  200, # room 153 hallway with common area
                  50, # room 154
                  20, # room 155
                  3, # room 156 facilities
                  rep(4,10), # room 157-166 offices
                  100, # room 167 large conference room
                  10, # room 168 research lab storage
                  100, # room 169 large research lab
                  75, # room 170 research desk/lab space
                  rep(50,7 ), # room 171 - 177 research storage 
                  100, # room 178 large research lab
                  rep(100,2), # room 179 and 180 classroom
                  rep(10,6), # rooms 181 - 186 storage rooms
                  rep(4,7), # rooms 187 - 193 offices
                  rep(20, 2), # rooms 194 and 195 bathrooms
                  4, # room 196 office
                  50, # room 197 admin open space
                  4, # room 198 office
                  50, # room 199 admin open space
                  rep(4,7), # room 200 - 206 office space
                  3, # room 207 storage
                  50, # room 208 hallway
                  100, # room 209 
                  50, # room 210 small hallway
                  100, # room 211 hallway
                  rep(25,5), # rooms 212-216 small hallways
                  rep(100,3), # rooms 217-219
                  200, # room 220 large hallway/common area
                  100, # room 221 hallway
                  rep(50,2) # room 224 outside
) #  vector of the carrying capacities for each room - people
University_C <- c(University_C,sum(University_C))
#outside should only be able to hold the total capacity of the building 
Building_max <- University_C[N_rooms] #last room in C is the "outside" room
Max_capacity <- Prop_full*Building_max
delt <- Max_capacity/community_pop # proportionality constant ( what proportion of individuals from the community are in the building of interest)

University_setup <- Bld_setup_func(Community_output = Community_output,day = day,delt = delt,N_rooms =N_rooms)

University_Init_conds <-c(S=University_setup$S_x, I = University_setup$I_x, R = University_setup$R_x, P = University_setup$P_x)
University_graph<-graph_from_adjacency_matrix(University_adjacency_matrix,mode= "undirected")%>% 
  set_vertex_attr( "Infectious", value = University_setup$I_x) %>% 
  set_vertex_attr( "Room", value = seq(1:N_rooms))
layout <- layout_with_fr(University_graph)
plot(University_graph, layout = layout)
University_init_prob <- c(0,0,0,1,rep(0,nrow(University_adjacency_matrix)-4))
graph_to_use <- University_graph
University_Net<-set_vertex_attr(graph_to_use, "Infectious", value = University_setup$I_x)
#University__Net<-set_vertex_attr(graph_to_use, "Room", value = seq(1:N_rooms))
plot(University_Net)
University_network_fig<-ggplot(University_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes(aes(size = Infectious), color = "#C21807")+
  geom_nodetext_repel(aes(label=Room))+
  theme_blank()+
  scale_size(name = "Number of\nInfectious\npeople in\nrooms")+
  labs(title = "Network representation of a movie theater")+
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust=0))
University_network_fig
# png(filename = "University_network_fig.png",units="in", width=10, height=5, res=300)
# University_network_fig
# dev.off()

#create transition matrices for people and particles
#T_mov is the transition matrix for individuals moving between rooms

University_T_mov <- Create_Particle_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use,prop = 0.9)


University_theta_mov <-Create_Particle_T_Matrix(adjacency_matrix_to_use = adjacency_matrix_to_use, prop=0.5)

University_end_day_T_mov <- End_of_day_T_function(adjacency_matrix_to_use = adjacency_matrix_to_use)
End_of_day_T_mov_to_use <- University_end_day_T_mov



University_output <- data.frame(lsoda(y = University_Init_conds, func = Particle_model,times = times, parms = parms,adjacency_matrix_to_use=adjacency_matrix_to_use,T_mov=University_T_mov,theta_mov =University_theta_mov, C= University_C))
University_data_clean <- University_output%>% pivot_longer(cols = !time,
                                                           names_to = c("State", "Room"),
                                                           names_pattern = "([A-Za-z]+)(\\d+)",
                                                           values_to = "Number")
test<- Build_data_mod %>% filter(State == "P") %>% group_by(Room,time) %>% summarise(Room_Pop =sum(Number))
# University_model_output_fig<-ggplot(test,aes(x=time,y=Room_Pop,group=Room,color=as.numeric(Room)))+geom_line()+scale_color_stepsn(name="Room number",n.breaks=24,colors = viridis(24))+
#   labs(x = "Time (hours)", y= "Number of\ninfectious\nparticles",title = "Number of infectious particles across rooms thoughout a day")+
#   theme_linedraw()+
#   theme(legend.position = "right",axis.title.y = element_text(angle = 0,vjust = 0.5, hjust =0),
#         plot.title.position = "plot",
#         axis.title = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         axis.text = element_text(size = 11),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(1,units = "cm"))
# 
# University_model_output_fig

# png(filename = "University_model_output.png",units="in", width=6, height=4, res=300)
# University_model_output_fig
# dev.off()
University_data_ratios <- University_data_clean %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% group_by(time,Room)  %>% 
  mutate(N_x = S+I+R) %>% group_by(time,Room) %>% 
  mutate(K_x = (parms$s*I)/(parms$a*N_x), prop_to_K = P/K_x)
University_P_x_to_K_x_plot <- ggplot(University_data_ratios,aes(x =time, y =prop_to_K,group= Room,color=Room))+geom_line()+theme(legend.position = "none")+labs(title = "Univerisity - Risk - 224 rooms")

University_P_x_plot <- ggplot(University_data_ratios, aes(x=time, y = P, color = Room))+geom_line()+theme(legend.position = "none")+labs(title = "Univerisity - Risk - 224 rooms")


#### All plots ####
All_buildings_plot_risk_and_particles<- ggarrange(Church_P_x_K_x_plot+theme(legend.position = "none"),Church_P_x_plot+theme(legend.position = "none"),
                                                  Office_P_x_to_K_x_plot+theme(legend.position = "none"), Office_P_x_plot+theme(legend.position = "none"),
                                                  Movie_P_x_to_K_x_plot+theme(legend.position = "none"), Movie_P_x_plot+theme(legend.position = "none"),
                                                  University_P_x_to_K_x_plot,University_P_x_plot, ncol = 2,nrow = 4)

All_buildings_plot_risk_and_particles








