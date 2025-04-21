import networkx as nx
import numpy as np
import pandas as pd
import math

# alphaCore ranking of nodes in a complex, directed network
#
# Iteratively computes a node ranking based on a feature set derived from
# edge attributes and optionally static node features using the
# mahalanobis data depth function at the origin.
#
# @param graph A networkx directed graph.
# @param features A list of selected node features; default is ["all"]
# @param stepSize Defines the stepsize of each iteration as percentage of node count
#   Higher values (>0.1) speed up execution but lower ranking resolution.
#   Lower values (<0.1) provide finer ranking but increase runtime.
# @param startEpsi The epsilon to start with. Removes all nodes with depth>epsilon at start
#   Higher values (>1.0) remove more nodes early, emphasizing denser cores.
#   Lower values (<1.0) refine ranking progressively.
# @param expoDecay Dynamically reduces the step size, to have high cores with few nodes if true
#   Recommended when ranking precision is more important than speed.
# @return A dataframe of columns nodeID, alpha value, and batchID

def alphaCore(graph, features=["all"], stepSize = 0.1, startEpsi = 1, expoDecay = False):
    #1 Extract numerical node features from the graph; if unavailable, use computeNodeFeatures.
    data = extractFeatures(graph, features)

    #2 compute cov matrix to be used for all remainder of depth calculations
    matrix = data.drop("nodeID", axis=1)  # convert dataframe to numeric matrix by removing first column containing nodeID
    cov = np.cov(matrix.values.T)
    
    #3 calculate the Mahalanobis depth and add it to the respective row of the dataframe
    data['mahal'] = calculateMahalFromCenter(data, 0, cov)
    
    #4 
    epsi = startEpsi
    #5
    node = []
    alphaVals = []
    #6
    batch = []
    #7
    alpha = 1 - epsi
    #8
    alphaPrev = alpha
    #9
    batchID = 0
    
    #10 
    while graph.number_of_nodes() > 0:
        #11
        while True:
            depthFound = False  # to simulate do-while loop; used to check if there exists a node with depth >= epsi on current iteration
            #12 
            for row in data.itertuples():
                if row.mahal >= epsi:
                    depthFound = True
                    #13
                    node.append(row.nodeID)  # set node core
                    alphaVals.append(alphaPrev)
                    #14
                    batch.append(batchID)
                    #15
                    graph.remove_node(row.nodeID)
            #16
            batchID += 1
            #19 while condition of do-while loop of #11
            if graph.number_of_nodes() == 0 or not depthFound:
                break
            #17
            data = extractFeatures(graph, features)  # recompute node properties
            #18
            data['mahal'] = calculateMahalFromCenter(data, 0, cov)  # recompute depth
        #20 
        alphaPrev = alpha
        #21 
        if expoDecay and graph.number_of_nodes() > 0:  # exponential decay
            localStepSize = math.ceil(graph.number_of_nodes() * stepSize)
            data = data.sort_values(ascending=False, by=['mahal'])
            epsi = data.iloc[localStepSize - 1]['mahal']
        else:  # step decay
            epsi -= stepSize
        #22
        alpha = 1 - epsi
    #23
    return pd.DataFrame({'nodeID': node, 'alpha': alphaVals, 'batchID': batch})


# Extract numerical node features from the graph.
# 
# @param graph A NetworkX directed graph.
# @param features A list of feature names to extract. If None, defaults to ["all"].
# @return A DataFrame containing extracted numerical node features.
def extractFeatures(graph, features=["all"]):
    #Collect all node attributes and identify numerical features
    all_features = set()
    numeric_features = set()
    selected_features = set()
    data = {}

    for node_id, attributes in graph.nodes(data=True):
        for key, value in attributes.items():
            all_features.add(key)  # Track all available features
            if isinstance(value, (int, float)):  
                numeric_features.add(key)  # Track only numeric features

    #Determine the features to extract
    if features == ["all"]:
        if not numeric_features:
            print("No numerical node features found. Reverting to default AlphaCore node features.")
            return computeNodeFeatures(graph)
        selected_features = numeric_features  # Use all numeric features
    else:
        for feature in features:
            if feature not in all_features:
                print(f"Debug: Available features: {all_features}")
                raise ValueError(f"Feature '{feature}' not found in graph nodes. Please check your graph.")
            if feature in numeric_features:
                selected_features.add(feature)

        if not selected_features:
            print("None of the selected features contain numerical values. Reverting to default AlphaCore node features.")
            return computeNodeFeatures(graph)


    #Extract feature values for each node
    for node_id, attributes in graph.nodes(data=True):
        node_features = {}
        for feature in selected_features:
            if feature in attributes:
                node_features[feature] = attributes.get(feature)  
        data[node_id] = node_features

    #Convert to DataFrame
    df = pd.DataFrame.from_dict(data, orient="index").reset_index().rename(columns={"index": "nodeID"})
    return df



# Computes the node features of a given directed graph and returns a dataframe containing the features of each node
#
# @param graph A networkx directed graph
# @return A dataframe containing the computed node features with each row as a new entry and columns as different features
def computeNodeFeatures(graph):

    # additional or different node features can be computed and added as a column to the dataframe in the same manner as below

    nodeID = []
    inDegree = []
    outDegree = []
    inStrength = []
    outStrength = []
    for node in graph:
        nodeID.append(node)
        inDegree.append(graph.in_degree(node))
        outDegree.append(graph.out_degree(node))
        inStrength.append(graph.in_degree(node, "value"))
        outStrength.append(graph.out_degree(node, "value"))

    # currently adding inDegree, outDegree, inStrength, and outStrength to dataframe
    df = pd.DataFrame({"nodeID": nodeID, "inDegree": inDegree, "inStrength": inStrength, "outDegree": outDegree,
                       "outStrength": outStrength})
    return df



# Computes the mahalanobis depth from a given center of each row of a given set of data and returns it as an array
#
# @param data Dataframe where each row is a new entry and each column after the first (nodeID) is a type of data
# @param center A center value calculated with respect to when computing mahalanobis depth
# @param cov The covariance matrix of the data matrix
# @return An array containing the mahalanobis depth of each row entry of a given set of data
def calculateMahalFromCenter(data, center, cov):
    matrix = data.drop("nodeID", axis=1)  # convert dataframe to numeric matrix by removing first column containing nodeID
    x_minus_center = matrix.values - center
    x_minus_center_transposed = x_minus_center.T
    # Try computing the inverse of the covariance matrix; if it fails, fall back to the pseudo-inverse for stability.
    try: 
        inv_cov = np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        print("Covariance matrix is not invertible, using Moore-Penrose pseudo-inverse instead.")  
        inv_cov = np.linalg.pinv(cov)
    left = np.dot(x_minus_center, inv_cov)
    mahal = np.dot(left, x_minus_center_transposed)
    mahal_diag = np.diagonal(mahal)   # diagonal contains the depth vaulues corresponding to each row from matrix
    return np.reciprocal(1 + np.maximum(mahal_diag, 0))  # Ensure stability
