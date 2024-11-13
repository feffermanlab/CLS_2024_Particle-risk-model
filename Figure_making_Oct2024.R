#### Figure making ####
library(ggplot2)
library(tidyverse)
library(tidygraph)
library(igraph)
library(NatParksPalettes)
library(ggpubr)
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



Office_adjacency_matrix <- matrix(c(c(0,1,rep(0,4),1,0,1,rep(0,36-9)), # Hallway 1
                                    c(1,0,1,rep(0,32),1), # Hallway 2
                                    c(0,1,0,1,rep(0,5),1,1,1,rep(0,36-12)), # Hallway/room 3
                                    c(0,0,1,rep(0,6),rep(1,9),0,0,rep(1,3),rep(0,36-23)), # Hallway 4
                                    c(rep(0,3),1,0,1,1,rep(0,13),rep(1,5),rep(0,36-25)), # Hallway 5
                                    c(rep(0,4),1,0,1,1, rep(0,17),rep(1,4), rep(0,36-29)), # Hallway 6
                                    c(1,rep(0,3),1,1,0,rep(0,36-7)), # Hallway/room 7
                                    c(rep(0,5),1,rep(0,22),rep(1,3),rep(0,36-31)), # Hallway 8
                                    c(1,rep(0,30),1,1,1,1,0), # Hallway 9
                                    rep(c(0,0,1,1,rep(0,36-4)),3), # Rooms 10-12
                                    rep(c(0,0,0,1,rep(0,36-4)),5), # Rooms 13-17
                                    c(0,0,0,1,rep(0,14),1,1,rep(0,36-20)), # Room 18
                                    rep(c(rep(0,17),1,rep(0,36-18)),2), # Rooms 19 and 20
                                    rep(c(rep(0,3),1,1,rep(0,36-5)),3), # Rooms 21-23
                                    rep(c(rep(0,4),1,rep(0,36-5)),2), # Rooms 24 and 25
                                    rep(c(rep(0,5),1,rep(0,36-6)),3), # Rooms 26-28
                                    c(rep(0,5),1,0,1,rep(0,36-8)), # Room 29
                                    rep(c(rep(0,7),1,rep(0,36-8)),2), # Rooms 30 and 31
                                    rep(c(rep(0,8),1,rep(0,36-9)),3), # Rooms 32-34
                                    c(rep(0,8),1,rep(0,36-9)), # Room 35
                                    c(0,1,rep(0,34)) # Room 36
),
nrow=36,ncol = 36)

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

movie_theater_graph<- graph_from_adjacency_matrix(movie_adjacency_matrix, mode = "undirected")

University_adjacency_matrix<- matrix(c(rep(c(rep(0,34),1,rep(0,223-35)),6),#Basement floor -rooms 1-6
                                       rep(c(rep(0,35),1,rep(0,223-36)),4), #rooms 7-10
                                       c(rep(0,35),1,0,0,1,1,rep(0,223-40)), #room 11
                                       c(rep(0,36),1,rep(0,223-37)), #room 12
                                       c(rep(0,13),1,rep(0,22),1,rep(0,223-37)), #room 13
                                       c(rep(0,12),1,0,1,0,1,rep(0,20),1, rep(0,223-38)), #room 14
                                       c(rep(0,13),1,0,1,0,1,1,rep(0,18),1,rep(0,223-38)), #room 15
                                       c(rep(0,14),1,rep(0,4),1,1, rep(0,16),1,rep(0,223-38)), #room 16
                                       c(rep(0,13),1,rep(0,223-14)), #room 17
                                       rep(c(rep(0,14),1,rep(0,223-15)),2), #rooms 18-19
                                       rep(c(rep(0,15),1,rep(0,223-16)),2), #rooms 20-21
                                       c(rep(0,39),1,rep(0,223-40)), #room 22
                                       c(rep(0,39),1,1,rep(0, 223-41)), #room 23
                                       rep(c(rep(0,40), 1, rep(0,223-41)),2), #rooms 24-25
                                       c(rep(0,26),1,rep(0,13),1,1,rep(0,223-42)), #room 26
                                       c(rep(0,25),1,rep(0,16),1,rep(0,223-43)), #room 27
                                       c(rep(0,29),1,rep(0,11),1,rep(0,223-42)), # room 28
                                       c(rep(0,30),1,rep(0,10),1,rep(0,223-42)), #room 29
                                       c(rep(0,27),1,rep(0,223-28)), #room 30  
                                       c(rep(0,28),1,rep(0,223-29)), #room 31
                                       c(rep(0,41),1,rep(0,223-42)), #room 32
                                       c(rep(0,41),1,1,rep(0,223-43)), #room 33
                                       c(rep(0,42),1,rep(0,223-43)), #room 34
                                       #BASEMENT LEVEL HALLWAYS
                                       c(rep(1,6),rep(0,29),1,rep(0,7),1,rep(0,223-44)), #room 35 - (hallway)
                                       c(rep(0,6),rep(1,5),rep(0,25),1,0,1,rep(0,223-39)), #room 36 - (hallway)
                                       c(rep(0,11),1,1,rep(0,22),1,0,1,rep(0,6),1,rep(0,223-45)), #room 37 - (hallway)
                                       c(rep(0,13),rep(1,3),rep(0,20),1,0,0,1,rep(0,223-40)), #room 38 - (hallway)
                                       c(rep(0,10),1,rep(0,24),1,0,0,0,1,rep(0,223-40)), #room 39 - (hallway)
                                       c(rep(0,10),1,rep(0,10),1,1,rep(0,14),1,1,0,1,0,0,0,0,1,rep(0,223-46)), #room 40 - hallway
                                       c(rep(0,22),rep(1,4),rep(0,13),1,0,1,rep(0,223-42)), #room 41 - (hallway)
                                       c(rep(0,25),1,0,1,1,0,0,1,1,rep(0,9),1, rep(0,223-43)), #room 42 - (hallway)
                                       c(rep(0,32),1,1,rep(0,7),1,rep(0,4),1,rep(0,223-47)), #room 43 - (hallway)
                                       #STAIRS AND ELEVATORS * will fill in as I go, want to check the network of each floor
                                       c(rep(0,34),1,rep(0,45),1,rep(0,59),1,rep(0,66),1,rep(0,223-208)), #room 44 - (stairwell 1)
                                       c(rep(0,36),1,rep(0,14),1,rep(0,90),1,0,1,rep(0,64),1,0,1,rep(0,223-212)), #room 45 - (elevator 1)
                                       c(rep(0,39),1,rep(0,46),1,rep(0,62),1,rep(0,66),1,rep(0,223-217)), #room 46 - (elevator 2)
                                       c(rep(0,42),1,rep(0,48),1,rep(0,61),1,rep(0,66),1,rep(0,223-221)),  #room 47 - (stairwell 2)
                                       ################# LEVEL 1 ##
                                       rep(c(rep(0,80),1,rep(0,223-81)),2), # rooms 48 and 49
                                       c(rep(0,80),1,0,1,rep(0,223-83)), # room 50
                                       c(rep(0,80),1,rep(0,223-81)), #room 51
                                       c(rep(0,44),1,rep(0,7),1,rep(0,11),1,rep(0,15),1,1,rep(0,223-82)), # room 52
                                       c(rep(0,51),1,rep(0,223-52)), # room 53
                                       rep(c(rep(0,82),1,rep(0,223-83)),7), # rooms 54-60,
                                       c(rep(0,83),1,rep(0,223-84)), # room 61
                                       c(rep(0,66),1,rep(0,17),1,0,1,rep(0,223-87)), #room 62
                                       c(rep(0,84),1,rep(0,223-85)), # room 63
                                       c(rep(0,66),1,rep(0,17),1,rep(0,223-85)), # room 64
                                       c(rep(0,51),1,rep(0,32),1,rep(0,223-85)), # room 65
                                       c(rep(0,84),1,rep(0,223-85)), # room 66
                                       c(rep(0,61),1,0,1,rep(0,223-64)), #room 67
                                       rep(c(rep(0,86),1,rep(0,4),1,rep(0,223-92)),2), # room 68 and 69
                                       rep(c(rep(0,89),1,rep(0,223-90)),3), # room 70:72
                                       rep(c(rep(0,90),1,rep(0,223-91)),3), # room 73-75
                                       rep(c(rep(0,91),1,rep(0,223-92)),4), # room 76:79
                                       c(rep(0,85),1,0,1,0,0,1,rep(0,223-91)), # room 80
                                       c(rep(0,43),1,rep(0,3),rep(1,5),rep(0,29),1,rep(0,223-82)), # room 81 - (hallway)
                                       c(rep(0,51),1,rep(0,28),1,0,1,1,1,rep(0,223-85)), # room 82 - (Hallway)
                                       c(rep(0,49),1,0,0,0,rep(1,7),rep(0,21),1,0,1,rep(0,223-84)), # room 83 - (Hallway)
                                       c(rep(0,60),1,rep(0,20),1,1,0,1,1,1,1,1,rep(0,223-89)), # room 84 - (Hallway)
                                       c(rep(0,61),rep(1,5),rep(0,15),1,0,1,rep(0,223-84)), # room 85 - (Hallway)
                                       c(rep(0,79),1,0,0,0,1,0,0,0,1,rep(0,223-88)), # room 86 - (Hallway)
                                       c(rep(0,45),1,rep(0,15),1,rep(0,5),1,1,rep(0,14),1,rep(0,4),1,1,rep(0,223-90)), # room 87 - (Hallway)
                                       c(rep(0,79),1,rep(0,3),1,0,1,rep(0,4),1,rep(0,223-91)), # room 88 - (Hallway)
                                       c(rep(0,83),1,0,0,1,rep(0,3),1,rep(0,223-91)), # room 89 - (Hallway)
                                       c(rep(0,69),rep(1,3),rep(0,14),1, rep(0,223-87)), # room 90 - Hallway
                                       c(rep(0,72),1,1,1,rep(0,4),1, rep(0,7),1,1,0,0,1,1, rep(0,223-93)), # room 91 - (Hallway)
                                       c(rep(0,46),1,rep(0,20),1,1,rep(0,6),rep(1,4),rep(0,223-79)), # room 92 - (Hallway)
                                       ########################## STAIR ON LEVEL 1 room - 93###
                                       c(rep(0,90),1,rep(0,58),1,0,0,1,rep(0,63),1,0,0,1,rep(0,223-220)), # room 93 ( Staircase)
                                       ########################## BEGIN FLOOR 2 ###
                                       rep(c(rep(0,140),1,rep(0,223-141)),3), # rooms 94:96 
                                       rep(c(rep(0,141),1,rep(0,223-142)),4), # rooms 97:100
                                       rep(c(rep(0,143),1,rep(0,223-144)),4), #rooms 101: 104
                                       c(rep(0,141),1,1,rep(0,223-143)), # room 105
                                       c(rep(0,142),1,rep(0,223-143)), #room 106
                                       c(rep(0,144),1,rep(0,223-145)), # room 107 
                                       c(rep(0,112),1,1,rep(0,30),1,1,1,rep(0,223-147)), # room 108
                                       c(rep(0,144),1,1,1,rep(0,223-147)), #room 109
                                       c(rep(0,144),1,1,rep(0,223-146)), # room 110
                                       c(rep(0,145),1,rep(0,223-146)), # room 111
                                       c(rep(0,146),1, rep(0,223-147)), # room 112
                                       rep(c(rep(0,107),1,rep(0,223-108)),2), # rooms 113 and 114
                                       rep(c(rep(0,147),1,rep(0,223-148)),2), # rooms 115 and 116
                                       c(rep(0,149),1,rep(0,223-150)), # room 117
                                       rep(c(rep(0,148),1,rep(0,223-149)),3), # room 118:120
                                       c(rep(0,121),1,rep(0,25),1,rep(0,6),1, rep(0,223-155)), # room 121
                                       c(rep(0,120),1,rep(0,32),1,rep(0,223-154)), # room 122
                                       rep(c(rep(0,154),1,rep(0,223-155)),3), # rooms 123:125
                                       rep(c(rep(0,152),1,rep(0,223-153)),2), # rooms 126 and 127 
                                       rep(c(rep(0,153),1,rep(0,223-154)),5), # rooms 128:132 
                                       rep(c(rep(0,152),1,1,rep(0,223-154)),2), # rooms 133 and 134
                                       c(rep(0,135),1,0,1,1,rep(0,11),1,rep(0,223-151)), # room 135
                                       c(rep(0,134),1,0,0,0,1,rep(0,223-139)), # room 136
                                       c(rep(0,137),1,1,rep(0,223-139)), # room 137
                                       c(rep(0,134),1,0,1,rep(0,223-137)), # room 138
                                       c(rep(0,134),1,1,1,rep(0,223-137)), # room 139
                                       c(rep(0,150),1,rep(0,223-151)), # room 140
                                       ############ start of level 2 Hallways ###
                                       c(rep(0,43),1,rep(0,49),1,1,1,rep(0,45),1,rep(0,223-142)), # room 141 - (Hallway)
                                       c(rep(0,96),rep(1,4),rep(0,4),1,rep(0,35),1,0,1,1,rep(0,223-144)), # room 142 - (Hallway)
                                       c(rep(0,44),1, rep(0,59),1,1,rep(0,35),1,0,0,1,rep(0,223-145)), # room 143 - (Hallway)
                                       c(rep(0,100),rep(1,4),rep(0,45),1,1,rep(0,223-151)), # room 144 - (Hallway)
                                       c(rep(0,44),1,rep(0,61),rep(1,4),rep(0,32),1,rep(0,223-143)), # room 145 - (Hallway)
                                       c(rep(0,107),rep(1,4),rep(0,223-111)), # room 146 - (Hallway)
                                       c(rep(0,107),1,1,0,0,1,rep(0,35),1,rep(0,223-148)), # room 147 - (Hallway)
                                       c(rep(0,114),1,1,rep(0,4),1,rep(0,25),1,0,0,1,rep(0,223-150)), # room 148
                                       c(rep(0,117), 1,1,1,rep(0,29),1,rep(0,223-150)), # room 149 - (Hallway)
                                       c(rep(0,45),1,rep(0,46),1,rep(0,23),1,rep(0,26),1,rep(0,4),1,0,1,1,rep(0,223-152)), # room 150 - (Hallway)
                                       c(rep(0,134),1,0,0,0,1,1,rep(0,3),1,rep(0,5),1,0,0,1,rep(0,223-153)), # room 151 - (Hallway)
                                       c(rep(0,149),1,0,0,1,rep(0,223-153)), # room 152 - (Hallway)
                                       c(rep(0,125),1,1,rep(0,5),1,1,rep(0,4),1,rep(0,11),1,1,rep(0,223-152)),# room 153 - (Hallway)
                                       c(rep(0,46),1,rep(0,74),1,rep(0,4),rep(1,8),rep(0,18),1,0,1,rep(0,223-155)), # room 154 - (Hallway)
                                       c(rep(0,120),1,0,1,1,1,rep(0,28),1,rep(0,223-154)), # room 155 - (Hallway),
                                       ##################### BEGIN LEVEL 3 ### !!!!!!!!!!!!
                                       rep(c(rep(0,207),1,rep(0,223-208)),2), # room 156 and 157
                                       c(rep(0,207),1,1,rep(0,223-209)), # room 158
                                       rep(c(rep(0,208),1,rep(0,223-209)),4), # rooms 159-162
                                       rep(c(rep(0,210),1,rep(0,223-211)),4), # rooms 163-166
                                       c(rep(0,208),1,1,rep(0,223-210)), # room 167
                                       c(rep(0,211),1,rep(0,223-212)), # room 168
                                       c(rep(0,173),1,1,rep(0,36),1,1,1,1,rep(0,223-215)), # room 169
                                       c(rep(0,211),1,1,1,rep(0,223-214)), # room 170
                                       c(rep(0,211),1,1,rep(0,223-213)), # room 171
                                       c(rep(0,212),1,rep(0,223-213)), # room 172
                                       c(rep(0,213),1,rep(0,223-214)), # room 173
                                       rep(c(rep(0,168),1,rep(0,223-169)),2), # rooms 174 and 175
                                       rep(c(rep(0,214),1,rep(0,223-215)),2), # rooms 176 and 177
                                       c(rep(0,214),1,rep(0,6),1,1), # room 178
                                       c(rep(0,222),1), # room 179
                                       rep(c(rep(0,221),1,0),3), # rooms 180-182
                                       c(rep(0,215),1,1,rep(0,223-217)), # room 183
                                       rep(c(rep(0,215),1,rep(0,223-216)),3), # rooms 184-186
                                       c(rep(0,219),1,rep(0,223-220)), # room 187
                                       c(rep(0,219),1,1,rep(0,223-221)), # room 188
                                       rep(c(rep(0,220),1,rep(0,223-221)),4), # room 189-192
                                       c(rep(0,220),1,rep(0,223-221)), # room 193
                                       rep(c(rep(0,219),1,1,rep(0,223-221)),2), # rooms 194 and 195
                                       c(rep(0,217),1,rep(0,223-218)), # room 196
                                       c(rep(0,197),1,1,1,rep(0,17),1,rep(0,223-218)), # room 197
                                       c(rep(0,196),1,rep(0,223-197)), # room 198 
                                       c(rep(0,196),1,0,0,rep(1,5),rep(0,223-204)), # room 199
                                       c(rep(0,196),1,0,1,rep(0,223-199)), # room 200
                                       rep(c(rep(0,198),1,rep(0,223-199)),4), # rooms 201-204
                                       rep(c(rep(0, 217),1,rep(0,223-218)),3), # rooms 205-207
                                       ######## LEVEL 3 HALLWAYS ###
                                       c(rep(0,43),1,rep(0,111),1,1,1,rep(0,50),1,rep(0,223-209)),# room 208 - (Hallway)
                                       c(rep(0,157),rep(1,5),rep(0,4),1,rep(0,40),1,0,1,1,rep(0,223-211)), # room 209 - (Hallway)
                                       c(rep(0,44),1,rep(0,121),1, rep(0,41),1,0,0,1,rep(0,223-212)), # room 210 - (Hallway)
                                       c(rep(0,162),rep(1,4),rep(0,42),1,rep(0,7),1,1,rep(0,223-218)), # room 211 - (Hallway)
                                       c(rep(0,44),1,rep(0,122),rep(1,4),rep(0,38),1,rep(0,223-210)), # room 212 - (Hallway)
                                       c(rep(0,168),rep(1,4),rep(0,223-172)), # room 213 - (Hallway)
                                       c(rep(0,168),1,1,0,0,1,rep(0,41),1,rep(0,223-215)), # room 214 - (Hallway)
                                       c(rep(0,168),1, rep(0,6),1,1,1,rep(0,35),1,0,0,1,rep(0,223-217)), # room 215 - (Hallway)
                                       c(rep(0,182),rep(1,4),rep(0,30),1,rep(0,223-217)), # room 216 - (Hallway)
                                       c(rep(0,45),1,rep(0,46),1,rep(0,89),1,rep(0,27),1,rep(0,3),1,1,0,0,1,rep(0,223-219)), # room 217 - (Hallway)
                                       c(rep(0,195),1,1,rep(0,7),rep(1,3),rep(0,3),1,rep(0,8),1,rep(0,223-220)), # room 218 - (Hallway)
                                       c(rep(0,216),1,0,0,1,rep(0,223-220)), # room 219 - (Hallway)
                                       c(rep(0,92),1,rep(0,93),1,1,rep(0,5),1,1,0,0,0,1,rep(0,18),1,1,0,1,rep(0,223-221)), # room 220 - (Hallway)
                                       c(rep(0,46),1,rep(0,140),rep(1,8),rep(0,24),1,0,1,1), # room 221 - (Hallway)
                                       c(rep(0,177),1,0,1,1,1,rep(0,38),1,rep(0,223-221)), # room 222 - (Hallway)
                                       c(rep(0,177),1,1,rep(0,41),1,rep(0,223-221)) # room 223 - (Hallway)
                                       # room 223 - Outside
                                       
),nrow = 223,ncol = 223)


#### Building carrying capacities ####
#r set carrying capacities
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
)

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
                  100,
                  100,
                  100# room 221 hallway
) #  vector of the carrying capacities for each room - people





small_bld_3_rooms_graph<- graph_from_adjacency_matrix(small_bld_3_rooms, mode = "undirected")
#take a look
plot(small_bld_3_rooms_graph)
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

V(church_graph)$capacity <- Church_C
V(church_graph)$Room <- seq(1:length(Church_C))
set.seed(154780)
church_network <- ggplot(church_graph, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey75")+
  geom_nodes( color = UTK_colors_fun("turq"),aes(size = capacity))+
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
  geom_nodes( color = UTK_colors_fun("turq"),aes(size = capacity))+
  geom_nodetext_repel(aes(label=Room))+
  theme_blank()+
  scale_size_continuous(
    name = "Room capacity",
    breaks = seq(1, max(V(Office_graph)$capacity), by = 25),  # Set custom breaks
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
  geom_nodes( color = UTK_colors_fun("turq"),aes(size = capacity))+
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
  geom_nodes( color = UTK_colors_fun("turq"),aes(size = capacity))+
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

##### Model output ####
Max_Building_Capacity_Church <- sum(Church_C) 
# how full do we want our building capacity to be
Adj_Max_Building_Capacity_Church <- round(Prop_full*Max_Building_Capacity_Church,0) #Adjusted capacity - number of individuals that will be in the building such that building is at 'Prop_full' capacity

Church_data <- read.table(file = "Data/Church_output_Oct20_parms1.text",sep = " ")
Church_data_clean <- Church_data%>% pivot_longer(cols = !time,
                                                   names_to = c("State", "Room"),
                                                   names_pattern = "([A-Za-z]+)(\\d+)",
                                                   values_to = "Number")

Church_data_ratios <- Church_data_clean %>% pivot_wider(names_from = c(State),values_from = c(Number)) %>% group_by(time,Room)  %>% 
  mutate(N_x = (S+I+R)*Adj_Max_Building_Capacity_Church) #%>% group_by(time) %>% summarise(Total_prop=sum(N_x), x =sum(x)*Adj_Max_Building_Capacity_Church)
# %>% group_by(time,Room) %>% 
#   mutate(K_x = ((parms$s)*I)/((parms$a)*N_x), prop_to_K = P/K_x)
# Church_P_x_K_x_plot<-ggplot(Church_data_ratios,aes(x =time, y =prop_to_K,group= Room,color=Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")
ggplot(Church_data_ratios,aes(x=time,color = Room))+geom_line(aes(x = time, y = N_x/x))#+geom_line(aes(x = time, y = N_x))
Church_P_x_plot <- ggplot(Church_data_ratios, aes(x=time, y = P, color = Room))+geom_line()+labs(title = "Church - Risk - 31 rooms")
Church_P_x_K_x_plot



