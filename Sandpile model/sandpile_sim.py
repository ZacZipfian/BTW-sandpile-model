import numpy as np
import pandas as pd
import random

#Calculates Radius
def radius_shortcut(vector_j, near):
    if near == True:
        rad = ((vector_j[2] - vector_j[0])**2 + (vector_j[3] - vector_j[1])**2)**(0.5)
    else:
        rad = ((vector_j[4] - vector_j[0])**2 + (vector_j[5] - vector_j[1])**2)**(0.5)

    return rad[0]


# Updates sand matrix / event matrix / vector of radii (floats)
def cascader_updated4(N, sand_matrix, prev_updt, event_add, m, k, rolling_radius, time_counter ):
    #matrix manipulation
    inter_matrix = np.zeros((N+2,N+2), dtype=int)
    topple_map = sand_matrix // 4 #Clever
    UP = np.s_[0:N, 1:N+1]
    DOWN = np.s_[2:N+2, 1:N+1]
    LEFT = np.s_[1:N+1, 0:N]
    RIGHT = np.s_[1:N+1, 2:N+2]
    MIDDLE = np.s_[1:N+1, 1:N+1]
    inter_matrix[UP] += topple_map
    inter_matrix[DOWN] += topple_map
    inter_matrix[LEFT] += topple_map
    inter_matrix[RIGHT] += topple_map
    inter_matrix = inter_matrix[MIDDLE]
    prev_updt += topple_map

    #update sand matrix
    sand_matrix += inter_matrix
    sand_matrix -= 4*topple_map

    #update event vector
    event_add[0][1] = time_counter #life-time
    event_add[0][2] += np.sum(topple_map) #size
    event_add[0][3] = np.count_nonzero(prev_updt == 1)# unique area
    event_add[0][4] = (np.sum(sand_matrix) + 1)  #mass

    #calculate radius
    coordinate_vec = np.zeros((4,1))
    rowscols = np.argwhere(topple_map>0)
    coordinate_vec[0] = rowscols[0][0] 
    coordinate_vec[1] = rowscols[0][1] 
    coordinate_vec[2] = rowscols[len(rowscols)-1][0] 
    coordinate_vec[3] = rowscols[len(rowscols)-1][1] 

    rolling_rad1 = radius_shortcut([m, k, coordinate_vec[0], coordinate_vec[1], coordinate_vec[2], coordinate_vec[3]], True)
    rolling_rad2 = radius_shortcut([m, k, coordinate_vec[0], coordinate_vec[1], coordinate_vec[2], coordinate_vec[3]],  False)
    interum_radi = 0.0
    if rolling_rad1 > rolling_rad2:
        interum_radi = rolling_rad1
    else:
        interum_radi = rolling_rad2

    if interum_radi > rolling_radius:
        rolling_radius = interum_radi
        
    return rolling_radius 


# Sand pile simulator 
def sand_pile_simulation4(N,t):
    csv_name = "SP"+str(N)+","+str(t)
    with open(str(csv_name)+'.csv', 'w',newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter =',')
        writer.writerow(['time','life-time','size','area','mass','radius'])

    sand_matrix = np.zeros((N,N), dtype=int)
    #radius_floats = np.zeros((1,1))[0]
    #event_matrix = np.zeros((1,5), dtype=int) 
    prev_updt = np.zeros((N,N), dtype=int)

    for i in range(t):
        check = np.any(sand_matrix >= 4)

        if check == False:
           m = random.randint(0, (N-1)) 
           k = random.randint(0, (N-1)) 
           sand_matrix[m,k] += 1 

        else:
            event_add = np.zeros((1,5), dtype=int) #[0:"time", 1:"life-time", 2:"size", 3:"area", 4:"mass"]
            event_add[0][0] = i
            rolling_radius = 0.0
            time_counter = 0

            while check == True:
                time_counter += 1
                rolling_radius = cascader_updated4(N, sand_matrix, prev_updt, event_add, m, k, rolling_radius, time_counter )#
                check = np.any(sand_matrix >= 4)

            #event_matrix = np.concatenate((event_matrix, event_add))
            #radius_floats = np.concatenate((radius_floats, [rolling_radius]))
            with open(str(csv_name)+'.csv', 'a',newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter =',')
                writer.writerow([event_add[0][0],event_add[0][1],event_add[0][2],event_add[0][3],event_add[0][4],rolling_radius])

            m = random.randint(0, (N-1)) 
            k = random.randint(0, (N-1)) 
            sand_matrix[m,k] += 1 
            prev_updt = np.zeros((N,N), dtype=int)

    return sand_matrix 


#Choose variables and run simulation
N = 40
t = 20000
sand_matrix= sand_pile_simulation4(N,t)
