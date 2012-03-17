'''
Created on Dec 24, 2011

@author: cardona
'''

from networkx import DiGraph, is_directed_acyclic_graph, dfs_successors 
#from networkx.exception import NetworkXException, NetworkXError
from networkx import single_source_shortest_path_length,all_pairs_shortest_path_length,dijkstra_path_length

import numpy,pyparsing,copy

from .eNewick import eNewickParser
from .utils import total_cmp
import permutations
from .memoize import memoize_method

class PhyloNetwork(DiGraph):
    """
    Main class for Phylogenetic Networks (with nested taxa).
    """

    def __init__(self, data=None, name='', eNewick=None, ignore_prefix=None, id_offset=0):
        # initialization here
        DiGraph.__init__(self,data)
        self.name=name
        self._labels={}
        self._lastlabel=id_offset
        if eNewick != None:
            self.from_eNewick(eNewick,ignore_prefix=ignore_prefix)
    
    @memoize_method
    def is_phylogenetic_network(self):
        if not is_directed_acyclic_graph(self):
            return False
        return True
    
    @memoize_method
    def taxa(self):
        """
        Returns the taxa (set of labels) of self.
        """
        taxa = list(set(self._labels.values()))
        taxa.sort()
        return taxa

    def label(self,node):
        """
        Returns the label of node, or None if not labelled.
        """
        return self._labels.get(node)

    @memoize_method
    def node_by_taxa(self,taxa):
        """
        Returns the node labelled by taxa or None if no node is labelled by taxa.
        Important: If more than one node is labelleb by taxa, then a random one is
        returned.
        """
        for node in self.labelled_nodes():
            if self.label(node) == taxa:
                return node
        return None

    @memoize_method
    def all_nodes_by_taxa(self,taxa):
        """
        Returns all nodes labelled by taxa.
        """
        return [node for node in self.labelled_nodes() if self.label(node)==taxa]

    def is_tree_node(self,u):
        """
        Returns True if u is a tree node, False otherwise.
        """
        return self.in_degree(u)<=1
    
    def is_hybrid_node(self,u):
        """
        Returns True if u is a tree node, False otherwise.
        """
        return self.in_degree(u)>1

    def is_leaf(self,u):
        """
        Returns True if u is a leaf, False otherwise.
        """
        return self.out_degree(u)==0
    
    def is_root(self,u):
        """
        Returns True if u is a root, False otherwise.
        """
        return self.in_degree(u)==0

    def is_elementary_node(self,u):
        return ((self.in_degree(u)<=1) and (self.out_degree(u)==1))
    
    def is_labelled(self,u):
        """
        Returns True if u is a labelled node, False otherwise.
        """
        return u in self._labels
        
    @memoize_method
    def leaves(self):
        """
        Returns the set of leaves of self.
        """
        leaves = filter(self.is_leaf, self.nodes())
        leaves.sort()
        return leaves

    @memoize_method
    def roots(self):
        """
        Returns the set of roots of self.
        """
        roots = filter(self.is_root, self.nodes())
        roots.sort()
        return roots

    @memoize_method
    def labelled_nodes(self):
        """
        Returns de set of labelled nodes.
        """
        return self._labels.keys()
    
    @memoize_method
    def unlabelled_nodes(self):
        return list(set(self.nodes())-set(self.labelled_nodes()))
            
    @memoize_method
    def interior_nodes(self):
        return list(set(self.nodes())-set(self.leaves()))
    
    @memoize_method
    def elementary_nodes(self):
        return filter(self.is_elementary_node, self.nodes())
            
    @memoize_method
    def depth(self,u):
        return min([dijkstra_path_length(self,root,u) for root in self.roots()])
        
    @memoize_method
    def height(self,u):
        """
        Returns the height of u. 
        """
        if self.is_leaf(u):
            return 0
        else:
            return max(map(self.height, self.successors(u)))+1

    @memoize_method
    def mu(self,u):
        if self.is_leaf(u):
            mu = numpy.zeros(len(self.taxa()),int)
        else:
            mu = sum(map(self.mu,self.successors(u)))
        if self.is_labelled(u):
            pos=self.taxa().index(self.label(u))
            mu[pos] += 1
        return mu

    @memoize_method
    def mu_string(self):
        return '-'.join([str(self.mu(u)) for u in self.sorted_nodes()])

    @memoize_method
    def sorted_nodes(self):
        nodes = self.nodes()[:]
        nodes.sort(cmp=lambda u,v:total_cmp(self.mu(u),self.mu(v)))
        return nodes
        
    def getlabel(self):
        try:
            self._lastlabel += 1
        except:
            self._lastlabel = 1
        return '_%d' % (self._lastlabel)
            
            
    def walk(self,parsed,ignore_prefix=None):
        if isinstance(parsed, pyparsing.ParseResults):
            
            if 'tag' in parsed:
                internal_label='#'+str(parsed['tag'])
            else:
                internal_label=self.getlabel()
            if 'length' in parsed:
                pass
            self.add_node(internal_label)
            if 'label' in parsed:
                self._labels[internal_label]=parsed['label']
            for child in parsed:
                child_label=self.walk(child,ignore_prefix=ignore_prefix)
                if child_label:
                    self.add_edge(internal_label,child_label)
            return internal_label

    def from_eNewick(self,string,ignore_prefix=None):
        #parsed=eNewickParser(string)[0]
        try:
            parsed=eNewickParser(string)[0]
        except pyparsing.ParseException:
            raise 'Malformed eNewick'
            #raise JoQueSe
            return False
        self.walk(parsed,ignore_prefix=ignore_prefix)
        self.cache = {}

    def eNewick_node(self,u,visited):
        if self.is_leaf(u):
            #return self._labels[u]
            return self._labels.get(u,'')
        if u in visited:
            return u
        visited.append(u)
        children=map(lambda x:self.eNewick_node(x,visited),self.successors(u))
        internal=','.join(children)
        mylabel=self.label(u) or ''
        if self.is_hybrid_node(u):
            mylabel+=u
        return '('+internal+')'+mylabel

    def __str__(self):
        return self.eNewick()
    
    def eNewick(self):
        visited=[]
        string = ''
        for root in self.roots():
            string += self.eNewick_node(root,visited)+';'
        return string
    
    @memoize_method
    def descendant_nodes(self,u):
        return sum(dfs_successors(self,u).values(),[]) + [u]
        #return dfs_successors(self,u)

    @memoize_method
    def descendant_taxa(self,u):
        return [self.label(desc) for desc in self.descendant_nodes(u) if self.is_labelled(desc)]

    @memoize_method
    def strict_descendant_nodes(self,u):
        if self.is_root(u):
            return self.descendant_nodes(u)
        pruned = copy.deepcopy(self)
        pruned.cache = {}
        pruned.remove_node(u)
        desc_pruned = []
        for root in self.roots():
            desc_pruned.extend(dfs_successors(pruned,root))
        return [desc for desc in self.descendant_nodes(u) if not desc in desc_pruned]

    @memoize_method
    def strict_descendant_taxa(self,u):
        return [self.label(desc) for desc in self.strict_descendant_nodes(u) if self.is_labelled(desc)]

    @memoize_method
    def ancestors(self,taxon):
        return [u for u in self.sorted_nodes() if taxon in self.descendant_taxa(u)]

    @memoize_method
    def strict_ancestors(self,taxon):
        return [u for u in self.sorted_nodes() if taxon in self.strict_descendant_taxa(u)]

    @memoize_method
    def CSA(self,tax1,tax2):
        return [u for u in self.ancestors(tax1) if
                (u in self.ancestors(tax2)) and
                ((u in self.strict_ancestors(tax1)) or
                 (u in self.strict_ancestors(tax2)))]

    @memoize_method
    def LCSA(self,tax1,tax2):
        csa=self.CSA(tax1,tax2)
        #print self,tax1,tax2,csa
        csa.sort(lambda x,y:cmp(self.height(x),self.height(y)))
        return csa[0]        

    @memoize_method
    def nodal_matrix(self):
        n=len(self.taxa())
        matrix=numpy.zeros((n,n),int)
        dicdist=all_pairs_shortest_path_length(self)
        for i in range(n):
            ti=self.taxa()[i]
            for j in range(i,n):
                tj=self.taxa()[j]
                lcsa=self.LCSA(ti,tj)
                matrix[i,j]=dicdist[lcsa][self.node_by_taxa(ti)]
                matrix[j,i]=dicdist[lcsa][self.node_by_taxa(tj)]
        return matrix

    def nodal_area(self):
        mat=self.nodal_matrix()
        #mat=mat+mat.transpose()
        return sum(abs(mat.flatten()))

    @memoize_method
    def cophenetic_matrix(self):
        n=len(self.taxa())
        matrix=numpy.zeros((n,n),int)
        for i in range(n):
            ti=self.taxa()[i]
            for j in range(i,n):
                tj=self.taxa()[j]
                lcsa=self.LCSA(ti,tj)
                matrix[i,j]=self.depth(lcsa)
        return matrix

    def common_taxa(self,net2):
        common=[]
        taxa1=self.taxa()
        taxa2=net2.taxa()
        for taxon in taxa1:
            if taxon in taxa2:
                common.append(taxon)
        return common
            
    def common_taxa_leaves(self,net2):
        common=[]
        taxa1=filter(lambda l:self.is_leaf(self.node_by_taxa(l)),self.taxa())
        taxa2=filter(lambda l:net2.is_leaf(net2.node_by_taxa(l)),net2.taxa())
        for taxon in taxa1:
            if taxon in taxa2:
                common.append(taxon)
        return common
            
    
    def topological_restriction(self,subtaxa):
        restricted=copy.deepcopy(self)
        for taxon in restricted.taxa():
            if not taxon in subtaxa:
                for u in restricted.all_nodes_by_taxa(taxon):
                    del restricted._labels[u]
        restricted.cache = {}

        while True:
            leaves_to_delete=[u for u in restricted.nodes() if \
                                  restricted.is_leaf(u) and \
                                  not u in restricted._labels ]
            if not leaves_to_delete:
                break
            else:
                restricted.remove_nodes_from(leaves_to_delete)
        for u in restricted.nodes():
            if restricted.is_elementary_node(u) and \
                    not restricted.is_labelled(u):
                for parent in restricted.predecessors(u):
                    for child in restricted.successors(u):
                        restricted.add_edge(parent,child)
                restricted.remove_node(u)
        restricted.cache = {}
        return restricted

    @memoize_method
    def has_nested_taxa(self):
        for node in self.labelled_nodes():
            if not self.is_leaf(node):
                return True
        return False

    def matching_representation(self):
        try:
            return self._matching_representation
        except:
            self._matching_representation={}
            for u in self.nodes():
                self._matching_representation[u]=0
            for u in self.leaves():
                pos=self.taxa().index(self.label(u))
                self._matching_representation[u]=pos+1
            h=1
            i=len(self.leaves())+1
            thislevel=filter(lambda u:self.height(u)==h,self.nodes())
            while thislevel:
                minims={}
                for u in thislevel:
                    minims[u]=min([self._matching_representation[v] for \
                                       v in self.successors(u)])
                thislevel.sort(lambda x,y: minims[x]-minims[y])
                for k in range(len(thislevel)):
                    self._matching_representation[thislevel[k]]=i
                    i += 1
                h += 1
                thislevel=filter(lambda u:self.height(u)==h,self.nodes())
            return self._matching_representation
                          
    def matching_permutation(self):
        try:
            return self._matching_permutation
        except:
            permutation={}
            representation=self.matching_representation()
            for u in self.nodes():
                children=self.successors(u)
                children.sort(cmp=lambda u,v:representation[u]-representation[v])
                for i in range(len(children)):
                    permutation[representation[children[i]]]=representation[children[(i+1) % len(children)]]
            self._matching_permutation=permutations.Permutation(permutation)
            return self._matching_permutation

    def cluster(self,u):
        cl=[]
        dictio=single_source_shortest_path_length(self,u)
        for node in dictio.keys():
            if self.label(node):
                cl.append(self.label(node))
        cl.sort()
        cl=tuple(cl)
        return cl

    def cluster_representation(self):
        cls=map(self.cluster,self.nodes())
        cls.sort()
        return cls

    def nested_label(self,node):
        try:
            return self._nested_label[node]
        except:
            if not hasattr(self,'_nested_label'):
                self._nested_label = {}
            if self.is_leaf(node):
                return self.label(node)
            else:
                desc_labels = [self.nested_label(u) for u in self.successors(node)]
                desc_labels.sort()
                self._nested_label[node] = '{'+','.join(desc_labels)+'}'
                return self._nested_label[node]

    def nested_label_representation(self):
        nls=map(self.nested_label,self.nodes())
        return set(nls)

    
