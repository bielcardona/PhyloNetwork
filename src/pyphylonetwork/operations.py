import copy

def push_and_hang(net,u,newtaxa):
    """
    Builds a new Phylogenetic Tree from net and u by hanging a leaf labeled newtaxa 
    between u and its parent.
       |          |
       |    =>    |\
       u          u \newtaxa
      / \        / \
    """
    parents=net.predecessors(u)
    newnet=copy.deepcopy(net)
    w=newnet._generate_new_id()
    for parent in parents:
        newnet.remove_edge(parent,u)
        newnet.add_edge(parent,w)
    newnet.add_edge(w,u)
    newleaf=newnet._generate_new_id()
    newnet.add_edge(w,newleaf)
    newnet._labels[newleaf]=newtaxa
    #newnet.forget()
    newnet.cache = {}
    return newnet

def hold_and_hang(net,u,newtaxa):
    """
    Builds a new Phylogenetic Tree from net and u by hanging a leaf labeled newtaxa 
    from u
       |    =>    |
       u          u
      /|         /|\
                    \newtaxa
    """
    newnet=copy.deepcopy(net)
    newleaf = newnet._generate_new_id()
    newnet.add_edge(u,newleaf)
    newnet._labels[newleaf]=newtaxa
    #newnet.forget()
    newnet.cache = {}
    return newnet

def push_and_label(net,u,newtaxa):
    """
    Builds a new Phylogenetic Tree from net and u by inserting an elementary node
    between u and its parent and labeling it with newtaxa
       |          |
       |          newtaxa
       |    =>    |
       u          u
      / \        / \
    """
    parents=net.predecessors(u)
    newnet=copy.deepcopy(net)
    newnode = newnet._generate_new_id()
    for parent in parents:
        newnet.remove_edge(parent,u)
        newnet.add_edge(parent,newnode)
    newnet.add_edge(newnode,u)
    newnet._labels[newnode]=newtaxa
    #newnet.forget()
    newnet.cache = {}
    return newnet

def hold_and_label(net,u,newtaxa):
    """
    Builds a new Phylogenetic Tree from net and u by labeling u with newtaxa
       |          |
       |    =>    |
       u          u=newtaxa
      / \        / \
    """
    if net.is_labelled(u):
        return
    newnet=copy.deepcopy(net)
    newnet._labels[u]=newtaxa
    #newnet.forget()
    newnet.cache = {}
    return newnet
