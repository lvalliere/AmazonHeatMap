import csv
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt  

# reads the edge list data from a CSV 
filename = "Amazon_L1_20180101_2100.txt"
amazon_surface = ChemicalCase(filename)
                
for i in range(1,914):
    # For all edges coming in and out of a reaction node, set to half the rxn timescale 
    for u,v in B.in_edges("R"+str(i)):
        B[u][v]['timescale'] = 0.5 / (amazon_surface.reaction_rates[i-1] + 1e-20) # 1e-20 to avoid divide by 0
    for u,v in B.out_edges("R"+str(i)):
        B[u][v]['timescale'] = 0.5 / (amazon_surface.reaction_rates[i-1] + 1e-20)


rxn_node_list = ['R' + str(i) for i in range(1,914)]
spc_node_list = [node for node in list(B.nodes) if node not in rxn_node_list]



# Calculate shortest path lengths
num_nodes = len(spc_node_list)
distance_matrix = np.full((num_nodes, num_nodes), np.nan)

# Calculate shortest path lengths
for i, source_node in enumerate(spc_node_list):
    for j, target_node in enumerate(spc_node_list):
        if source_node != target_node and nx.has_path(B, source_node, target_node):
            path = nx.shortest_path(B, source=source_node, target=target_node, weight='timescale')
            timescale = nx.path_weight(B, path, weight='timescale')
            distance_matrix[i, j] = timescale

# Plot the distance matrix
plt.imshow(distance_matrix, cmap='viridis', interpolation='nearest')
plt.colorbar(label="Chemical Pathway Timescale [seconds/molec/cc]")
plt.xlabel("Target Species Index")
plt.ylabel("Source Species Index")
plt.title("Shortest Path Lengths Between Species")
plt.show()

