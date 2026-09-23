import networkx as nx
import matplotlib.pyplot as plt

def generate_hasse_diagram(n):
    # Create a directed graph
    G = nx.DiGraph()
    # key = 1

    # Add nodes for each subgroup
    subgroups = [i+1 for i in range(n) if n % (i + 1) == 0]
    for subgroup in subgroups:
        G.add_node(subgroup,p_value=0)

    # Add edges for inclusion relations
    for i in subgroups:
        for j in subgroups:
            if i != j and (i) % (j) == 0:
                flag = 1
                s = subgroups.copy()
                s.remove(i)
                s.remove(j)
                for k in s:
                    if (k) % (j) == 0 and (i) % (k) == 0:
                        flag = 0
                if flag == 1:
                    G.add_edge(j, i)

    # Draw the Hasse diagram
    #pos = subgroups #nx.spring_layout(G)
    pos = nx.shell_layout(G)
   # pos = nx.bfs_tree(G, source=1)
    nx.draw(G, pos, with_labels=True, node_size=700, node_color='lightblue', font_size=10, font_weight='bold', arrows=True)
    plt.title(f"Hasse Diagram of Cyclic Group of Order {n}")
    plt.show()

    ## new_var = G.get_edge_data(1,n)
    ## print(G[1])


# Example usage: Generate Hasse diagram for cyclic group of order 12

# generate_hasse_diagram(24)
