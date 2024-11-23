#RECTANGLE SEMI MODEL

import numpy as np
import multiprocessing
import csv
import time


def euclidean_distance(A, B):
    return abs(A[0] - B[0]), abs(A[1] - B[1])

#Defining grid system

def initialize_grid_system(L, k1, k2):
    num_rows, num_cols = L // k2, L // k1
    grid_data = {i: [] for i in range(num_rows * num_cols)}
    return num_rows, num_cols, grid_data

#Root finding and union

def find_root(site, ptr):
    root = site
    while ptr[int(root)] >= 0:
        root = ptr[int(root)]
    while site != root:
        parent = ptr[int(site)]
        ptr[int(site)] = root
        site = parent
    return root


def union_objects(site1, site2, ptr):
    root1, root2 = find_root(site1, ptr), find_root(site2, ptr)
    
    if root1 != root2:
        if ptr[int(root1)] < ptr[int(root2)]:  # Union by size
            ptr[int(root1)] += ptr[int(root2)]
            ptr[int(root2)] = root1
        else:
            ptr[int(root2)] += ptr[int(root1)]
            ptr[int(root1)] = root2
    return ptr

#Calculation of neighboring grids

def get_grid_neighbors(label, L, k1, k2):
    cols = L // k1  
    rows=L//k2
    neighbors = []
    
    is_left_boundary = (label % cols == 0)
    is_right_boundary = ((label + 1) % cols == 0)

    
    for dy in range(-1, 2):  
        for dx in range(-1, 2):  
           
            if is_left_boundary and dx == -1:
                continue  
            if is_right_boundary and dx == 1:
                continue  

            new_label = label + dy * cols + dx

            if 0 <= new_label < cols * (rows):
                neighbors.append(new_label)
    
    return neighbors

def boundary(L):
    
    up=np.zeros(L//k1)
    down=np.zeros(L//k1)
    up[0:L//k1]=0
    down[0:L//k1]=1
    
    return up,down

# Checking for percolation in vertical direction (Rectangles of length k1 and width k2 are objects to be dropped)

def simulate(L, k1, k2, num_rectangles,seedvalue):
    rng=np.random.Generator(np.random.MT19937(seedvalue))
    num_rows, num_cols, grid_data = initialize_grid_system(L, k1, k2)
    ver_perc_prob = np.zeros(num_rectangles)
    ptr = np.zeros(3*num_rectangles, dtype=int)  # Initializing pointer
    up,down=boundary(L)
    ptr[0]=-L//k1  # Setting pointer for boundary rectangles
    ptr[1]=-L//k1

    initlab=1  # Initial label set to 1 later incremented at each step

    for i in range(num_rectangles):
        x, y = rng.uniform(0, L), rng.integers(0, L)+0.5  # Generating x,y coordinates
        initlab+=1
        
        gridx = min(max(int(x // k1), 0), num_cols - 1)
        y1=min(max(int(y)//k2, 0), num_rows - 1)
        gridy = int(y1)
        grid_label = gridy * num_cols + gridx
        grid_data[grid_label].append((initlab, (x, y)))
        ptr[initlab]=-1
        neighbor_grid_labels = get_grid_neighbors(grid_label, L, k1, k2)
            
        neighbors = [(idx, coord) for n_label in neighbor_grid_labels for idx, coord in grid_data[n_label]]
        for label, neighbor_coord in neighbors:
            if label!=initlab:
                dx, dy = euclidean_distance((x,y), neighbor_coord)
                if dx<=k1 and dy<=k2:           #Overlapping condition
                    ptr = union_objects(initlab, label, ptr) #Union
           
        if y<=k2/2:
                ptr=union_objects(1,initlab,ptr)
        elif abs(L-(y))<=k2/2:
                ptr=union_objects(0,initlab,ptr)
        
       
            

        if x < k1 / 2:
            initlab+=1
            gridx = min(max(int((x+L) // k1), 0), num_cols - 1)
            y1=min(max(int(y)//k2, 0), num_rows - 1)
            gridy = int(y1)
            grid_label = gridy * num_cols + gridx
            grid_data[grid_label].append((initlab, (x+L, y)))
            neighbor_grid_labels = get_grid_neighbors(grid_label, L, k1, k2)
            ptr[initlab]=-1
            neighbors = [(idx, coord) for n_label in neighbor_grid_labels for idx, coord in grid_data[n_label]]
            for label, neighbor_coord in neighbors:
                if label!=initlab:
                    dx, dy = euclidean_distance((x+L,y), neighbor_coord)
                    if dx<=k1 and dy<=k2:
                        ptr = union_objects(initlab, label, ptr)
            if y<=k2/2:
                ptr=union_objects(1,initlab,ptr)
            elif abs(L-(y))<=k2/2:
                ptr=union_objects(0,initlab,ptr)
        
            

        if L - x < k1 / 2:
            initlab+=1
            gridx = min(max(int((x-L) // k1), 0), num_cols - 1)
            y1=min(max(int(y)//k2, 0), num_rows - 1)
            gridy = int(y1)
            grid_label = gridy * num_cols + gridx
            grid_data[grid_label].append((initlab, (x-L, y)))
            neighbor_grid_labels = get_grid_neighbors(grid_label, L, k1, k2)
            ptr[initlab]=-1
            neighbors = [(idx, coord) for n_label in neighbor_grid_labels for idx, coord in grid_data[n_label]]
            for label, neighbor_coord in neighbors:
                if label!=initlab:
                    dx, dy = euclidean_distance((x-L,y), neighbor_coord)
                    if dx<=k1 and dy<=k2:
                        ptr = union_objects(initlab, label, ptr)
            if y<=k2/2:
                ptr=union_objects(1,initlab,ptr)
            elif abs(L-(y))<=k2/2:
                ptr=union_objects(0,initlab,ptr)
        
            
        
        if y < k2 / 2:
            initlab+=1
            gridx = min(max(int(x // k1), 0), num_cols - 1)
            y1=min(max(int(y+L)//k2, 0), num_rows - 1)
            gridy = int(y1)
            grid_label = gridy * num_cols + gridx
            grid_data[grid_label].append((initlab, (x, y+L)))
            ptr[initlab]=-1
            neighbor_grid_labels = get_grid_neighbors(grid_label, L, k1, k2)
            
            neighbors = [(idx, coord) for n_label in neighbor_grid_labels for idx, coord in grid_data[n_label]]
            for label, neighbor_coord in neighbors:
                if label!=initlab:
                    dx, dy = euclidean_distance((x,y+L), neighbor_coord)
                    if dx<=k1 and dy<=k2:
                        ptr = union_objects(initlab, label, ptr)
            
            if y+L<=k2/2:
                ptr=union_objects(1,initlab,ptr)
            elif abs(L-(y+L))<=k2/2:
                ptr=union_objects(0,initlab,ptr)
        
            
            
        
        if L - y < k2 / 2:
            initlab+=1
            gridx = min(max(int((x) // k1), 0), num_cols - 1)
            y1=min(max(int(y-L)//k2, 0), num_rows - 1)
            gridy = int(y1)
            grid_label = gridy * num_cols + gridx
            grid_data[grid_label].append((initlab, (x, y-L)))
            neighbor_grid_labels = get_grid_neighbors(grid_label, L, k1, k2)
            ptr[initlab]=-1
            neighbors = [(idx, coord) for n_label in neighbor_grid_labels for idx, coord in grid_data[n_label]]
            for label, neighbor_coord in neighbors:
                if label!=initlab:
                    dx, dy = euclidean_distance((x,y-L), neighbor_coord)
                    if dx<=k1 and dy<=k2:
                        ptr = union_objects(initlab,label, ptr)
            
            if y-L<=k2/2:
                ptr=union_objects(1,initlab,ptr)
            elif abs(L-(y-L))<=k2/2:
                ptr=union_objects(0,initlab,ptr)
        
        
        
        if x < k1/ 2 and y < k2 / 2:
            initlab+=1
            gridx = min(max(int((x+L) // k1), 0), num_cols - 1)
            y1=min(max(int(y+L)//k2, 0), num_rows - 1)
            gridy = int(y1)
            grid_label = gridy * num_cols + gridx
            grid_data[grid_label].append((initlab, (x+L, y+L)))
            ptr[initlab]=-1
            neighbor_grid_labels = get_grid_neighbors(grid_label, L, k1, k2)
            
            neighbors = [(idx, coord) for n_label in neighbor_grid_labels for idx, coord in grid_data[n_label]]
            for label, neighbor_coord in neighbors:
                if label!=initlab:
                    dx, dy = euclidean_distance((x+L,y+L), neighbor_coord)
                    if dx<=k1 and dy<=k2:
                        ptr = union_objects(initlab,label, ptr)
           
            if y+L<=k2/2:
                ptr=union_objects(1,initlab,ptr)
            elif abs(L-(y+L))<=k2/2:
                ptr=union_objects(0,initlab,ptr)

            
        
        
      
        if x < k1 / 2 and L - y < k2/ 2:
            initlab+=1
            gridx = min(max(int((x+L) // k1), 0), num_cols - 1)
            y1=min(max(int(y-L)//k2, 0), num_rows - 1)
            gridy = int(y1)
            grid_label = gridy * num_cols + gridx
            grid_data[grid_label].append((initlab, (x+L, y-L)))
            neighbor_grid_labels = get_grid_neighbors(grid_label, L, k1, k2)
            ptr[initlab]=-1
            neighbors = [(idx, coord) for n_label in neighbor_grid_labels for idx, coord in grid_data[n_label]]
            for label, neighbor_coord in neighbors:
                if label!=initlab:
                    dx, dy = euclidean_distance((x+L,y-L), neighbor_coord)
                    if dx<=k1 and dy<=k2:
                        ptr = union_objects(initlab, label, ptr)
           
            if y-L<=k2/2:
                ptr=union_objects(1,initlab,ptr)
            elif abs(L-(y-L))<=k2/2:
                ptr=union_objects(0,initlab,ptr)
        
        

        if L - x < k1/ 2 and y < k2 / 2:
            initlab+=1
            gridx = min(max(int((x-L) // k1), 0), num_cols - 1)
            y1=min(max(int(y+L)//k2, 0), num_rows - 1)
            gridy = int(y1)
            grid_label = gridy * num_cols + gridx
            grid_data[grid_label].append((initlab, (x-L, y+L)))
            neighbor_grid_labels = get_grid_neighbors(grid_label, L, k1, k2)
            ptr[initlab]=-1
            neighbors = [(idx, coord) for n_label in neighbor_grid_labels for idx, coord in grid_data[n_label]]
            for label, neighbor_coord in neighbors:
                if label!=initlab:
                    dx, dy = euclidean_distance((x-L,y+L), neighbor_coord)
                    if dx<=k1 and dy<=k2:
                        ptr = union_objects(initlab,label, ptr)

            if y+L<=k2/2:
                ptr=union_objects(1,initlab,ptr)
            elif abs(L-(y+L))<=k2/2:
                ptr=union_objects(0,initlab,ptr)

            
        
        
        
        if L - x < k1 / 2 and L - y < k2 / 2:
            initlab+=1
            gridx = min(max(int((x-L) // k1), 0), num_cols - 1)
            y1=min(max(int(y-L)//k2, 0), num_rows - 1)
            gridy = int(y1)
            grid_label = gridy * num_cols + gridx
            grid_data[grid_label].append((initlab, (x-L, y-L)))
            neighbor_grid_labels = get_grid_neighbors(grid_label, L, k1, k2)
            ptr[initlab]=-1
            neighbors = [(idx, coord) for n_label in neighbor_grid_labels for idx, coord in grid_data[n_label]]
            for label, neighbor_coord in neighbors:
                if label!=initlab:
                    dx, dy = euclidean_distance((x-L,y-L), neighbor_coord)
                    if dx<=k1 and dy<=k2:
                        ptr = union_objects(initlab, label, ptr)
            
            if y-L<=k2/2:
                ptr=union_objects(1,initlab,ptr)
            elif abs(L-(y-L))<=k2/2:
                ptr=union_objects(0,initlab,ptr)
            
        
       
        
            
             
            

        down[0:L//k1]=find_root(1,ptr)
        up[0:L//k1]=find_root(0,ptr)
              
               
            

        if up[0]==down[0]:
            ver_perc_prob[i]=1
        

        if ver_perc_prob[i]==1:
           
            break
        
    return i+1     # return point at which percolation occurs for first time in vertical direction



    
L=512
k1=2
k2=1
num_rectangles=200000
seedvalue=1000

#Set seedvalues as per requirement

V=simulate(L,k1,k2,num_rectangles,seedvalue)
  
print(f"System percolates at Step: {V}")