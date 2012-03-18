from utils import total_cmp
from .exceptions import TaxaException



def mu_distance(net1,net2):
    """
    Compute the mu distance between two phylogenetic networks.
    """
#        if net1.taxa() != net2.taxa():
#            return
    nodes1=net1.sorted_nodes()[:]
    nodes2=net2.sorted_nodes()[:]
    d=0
    while(len(nodes1)>0 and len(nodes2)>0):
        x1=nodes1[0]
        x2=nodes2[0]
        comp=total_cmp(net1.mu(x1),net2.mu(x2))
        if comp==-1:
            del nodes1[0]
            d+=1
        elif comp==1:
            del nodes2[0]
            d+=1
        else:
            del nodes1[0]
            del nodes2[0]
    return d+len(nodes1)+len(nodes2)

#def is_isomorphic_fast(net1,net2):
#    return net1.mu_string() == net2.mu_string()

def nodal_distance_splitted(net1,net2,p=1,take_root=False, check=False):
    """
    Computes the nodal distance splitted between two phylogenetic networks.
    If check = True, then it checks if the two networks have the same taxa. Otherwise it will only check if the number of labels is equal.
    """
    
    if check:
        if not net1.taxa() == net2.taxa():
            raise TaxaException("Networks over different set of taxa")
         
    try:
        mat=net1.nodal_matrix()-net2.nodal_matrix()
    except:
        raise TaxaException("The number of possible labels is not equal")
    if p==1:
        return sum(abs(mat.flatten()))
    else:
        if take_root:
            return sum(abs(mat.flatten())**p)**(1.0/p)
        else:
            return sum(abs(mat.flatten())**p)

def nodal_distance_unsplitted(net1,net2,p=1,take_root=False, check=False):
    """
    Computes the nodal distance unsplitted between two phylogenetic networks.
    If check = True, then it checks if the two networks have the same taxa. Otherwise it will only check if the number of labels is equal.
    """
    
    if check:
        if not net1.taxa() == net2.taxa():
	    raise TaxaException("Networks over different set of taxa")
	  
    mat1=net1.nodal_matrix()
    mat1=mat1+mat1.transpose()
    mat2=net2.nodal_matrix()
    mat2=mat2+mat2.transpose()
    try:
        mat=net1.nodal_matrix()-net2.nodal_matrix()
    except:
        raise TaxaException("The number of possible labels is not equal")
    if p==1:
        return sum(abs(mat.flatten()))/2
    else:
        if take_root:
            return (sum(abs(mat.flatten())**p)/2)**(1.0/p)
        else:
            return (sum(abs(mat.flatten())**p)/2)

def cophenetic_distance(net1,net2,p=1,take_root=False, check=False):
    """
    Computes the nodal distance unsplitted between two phylogenetic networks.
    If check = True, then it checks if the two networks have the same taxa. Otherwise it will only check if the number of labels is equal.
    """
    
    if check:
        if not net1.taxa() == net2.taxa():
	    raise TaxaException("Networks over different set of taxa")
	  
    mat1=net1.cophenetic_matrix()
    mat2=net2.cophenetic_matrix()
    try:
        mat=mat1-mat2
    except:
        raise TaxaException("The number of possible labels is not equal")
    if p==1:
        return sum(abs(mat.flatten()))/2
    else:
        if take_root:
            return (sum(abs(mat.flatten())**p))**(1.0/p)
        else:
            return (sum(abs(mat.flatten())**p))

def transposition_distance(net1,net2):
    """
    Computes the transposition distance between two phylogenetic networks.
    """
    
#        pi1=permutations.Permutation(net1.matching_permutation())
#        pi2=permutations.Permutation(net2.matching_permutation())
    pi1=net1.matching_permutation()
    pi2=net2.matching_permutation()
    pi=pi2**(-1)*pi1
    cicles=pi.cycles()
    dist=0
    for cicle in cicles:
        dist+=len(cicle)-1
    return dist/2

def RF_distance(net1,net2):
    """
    Computes the RF distance between two phylogenetic networks.
    """
    
    return len(set(net1.cluster_representation())^set(net2.cluster_representation()))

def nested_label_distance(net1,net2):
    """
    Computes the nested label distance between two phylogenetic networks.
    """
    
    return len(net1.nested_label_representation() ^ net2.nested_label_representation())