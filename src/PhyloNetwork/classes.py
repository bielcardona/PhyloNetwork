'''
Created on Dec 24, 2011

@author: cardona
'''

from networkx import DiGraph, is_directed_acyclic_graph, dfs_successors 
from networkx.exception import NetworkXException, NetworkXError
from networkx import single_source_shortest_path_length,all_pairs_shortest_path_length,dijkstra_path_length

import numpy,pyparsing,copy

from .eNewick import eNewickParser
from .utils import total_cmp
import permutations
from .memoize import memoize

class PhyloNetwork(DiGraph):
    """
    Main class for Phylogenetic Networks (with nested taxa).
    """

    def __init__(self, data=None, name='', eNewick=None, ignore_prefix=None):
        # initialization here
        DiGraph.__init__(self,data)
        self.name=name
        self._labels={}
        if eNewick != None:
            self.from_eNewick(eNewick,ignore_prefix=ignore_prefix)
#        self.recompute()

    #===========================================================================
    # def forget(self):
    #    """
    #    Clears precomputed values of self.
    #    """
    #    forgetables=['_heights','_mus','_taxa','_leaves','_roots',
    #                 '_labelled_nodes','_sorted_nodes','_mu_string',
    #                 '_descendant_nodes','_descendant_taxa',
    #                 '_strict_descendant_nodes','_strict_descendant_taxa',
    #                 '_ancestors','_strict_ancestors','_nodal_matrix',
    #                 '_matching_representation','_matching_permutation',
    #                 '_nested_label'
    #                 ]
    #    for forgetable in forgetables:
    #        try:
    #            del self.__dict__[forgetable]
    #        except:
    #            pass
    #===========================================================================
    
    @memoize
    def taxa(self):
        """
        Returns the taxa (set of labels) of self.
        """
        taxa = self._labels.values()
        taxa.sort()
        return taxa

        #=======================================================================
        # try:
        #    return self._taxa
        # except:
        #    self._taxa = self._labels.values()
        #    self._taxa.sort()
        #    return self._taxa
        #=======================================================================

    def label(self,node):
        """
        Returns the label of node, or None if not labelled.
        """
        return self._labels.get(node)
        #=======================================================================
        # try:
        #    return self._labels[node]
        # except:
        #    return None
        #=======================================================================

    @memoize
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
        #=======================================================================
        # try:
        #    return self._node_by_taxa[taxa]
        # except:
        #    if not hasattr(self,'_node_by_taxa'):
        #        self._node_by_taxa={}
        #    self._node_by_taxa[taxa]=None
        #    for node in self.labelled_nodes():
        #        if self.label(node)==taxa:
        #            self._node_by_taxa[taxa]=node
        #            break
        #    return self._node_by_taxa[taxa]
        #=======================================================================

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
        
    @memoize
    def leaves(self):
        """
        Returns the set of leaves of self.
        """
        leaves = filter(self.is_leaf, self.nodes())
        leaves.sort()
        return leaves
        #=======================================================================
        # try:
        #    return self._leaves
        # except:
        #    self._leaves = filter(self.is_leaf, self.nodes())
        #    self._leaves.sort()
        #    return self._leaves
        #=======================================================================

    @memoize
    def roots(self):
        """
        Returns the set of roots of self.
        """
        roots = filter(self.is_root, self.nodes())
        roots.sort()
        return roots
        #=======================================================================
        # try:
        #    return self._roots
        # except:
        #    self._roots = filter(self.is_root, self.nodes())
        #    self._roots.sort()
        #    return self._roots
        #=======================================================================

    @memoize
    def labelled_nodes(self):
        """
        Returns de set of labelled nodes.
        """
        return self._labels.keys()
        #=======================================================================
        # try:
        #    return self._labelled_nodes
        # except:
        #    self._labelled_nodes=self._labels.keys()
        #    return self._labelled_nodes
        #=======================================================================
            
    @memoize
    def depth(self,u):
        return min([dijkstra_path_length(self,root,u) for root in self.roots()])
        
    @memoize
    def height(self,u):
        """
        Returns the height of u. 
        """
        if self.is_leaf(u):
            return 0
        else:
            return max(map(self.height, self.successors(u)))+1
        #=======================================================================
        # try:
        #    return self._heights[u]
        # except:
        #    if not hasattr(self,'_heights'):
        #        self._heights = {}
        #    if self.is_leaf(u):
        #        self._heights[u] = 0
        #    else:
        #        self._heights[u] = max(map(self.height, self.successors(u)))+1
        #    return self._heights[u]
        #=======================================================================

    @memoize
    def mu(self,u):
        if self.is_leaf(u):
            mu = numpy.zeros(len(self.taxa()),int)
        else:
            mu = sum(map(self.mu,self.successors(u)))
        if self.is_labelled(u):
            pos=self.taxa().index(self.label(u))
            mu[pos] += 1
        return mu
        #=======================================================================
        # try:
        #    return self._mus[u]
        # except:
        #    if not hasattr(self,'_mus'):
        #        self._mus = {}
        #    if self.is_leaf(u):
        #        self._mus[u]=numpy.zeros(len(self.taxa()),int)
        #    else:
        #        self._mus[u]=sum(map(self.mu,self.successors(u)))
        #    if self.is_labelled(u):
        #        pos=self.taxa().index(self.label(u))
        #        self._mus[u][pos]=1                
        #    return self._mus[u]
        #=======================================================================

    @memoize
    def mu_string(self):
        return '-'.join([str(self.mu(u)) for u in self.sorted_nodes()])
        #=======================================================================
        # try:
        #    return self._mu_string
        # except:
        #    self._mu_string=\
        #    '-'.join([str(self.mu(u)) for u in self.sorted_nodes()])
        #    return self._mu_string
        #=======================================================================

    @memoize
    def sorted_nodes(self):
        nodes = self.nodes()[:]
        nodes.sort(cmp=lambda u,v:total_cmp(self.mu(u),self.mu(v)))
        return nodes
        #=======================================================================
        # try:
        #    return self._sorted_nodes
        # except:
        #    self._sorted_nodes=self.nodes()[:]
        #    self._sorted_nodes.sort(
        #        cmp=lambda u,v:total_cmp(self.mu(u),self.mu(v)))
        #    return self._sorted_nodes
        #=======================================================================
        
    def getlabel(self):
        try:
            self._lastlabel += 1
        except:
            self._lastlabel = 1
        return '_%d' % (self._lastlabel)
            
    #===========================================================================
    # def gethyblabel(self):
    #    try:
    #        self._lasthyblabel += 1
    #    except:
    #        self._lasthyblabel = 1
    #    return '_%d' % (self._lasthyblabel)
    #===========================================================================
            
    def walk(self,parsed,ignore_prefix=None):
        if isinstance(parsed, pyparsing.ParseResults):
            
            #===================================================================
            # if 'label' in parsed:
            #    internal_label=parsed['label']
            #    if ignore_prefix and internal_label[0]==ignore_prefix:
            #        internal_label=self.getlabel()
            #    else:
            #        self._labels[internal_label]=internal_label
            # else:
            #    internal_label=self.getlabel()
            #===================================================================
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
        try:
            parsed=eNewickParser(string)[0]
        except pyparsing.ParseException:
            print 'Malformed eNewick'
            #raise JoQueSe
            return False
        self.walk(parsed,ignore_prefix=ignore_prefix)
        self.cache = {}

    def eNewick_node(self,u,visited):
        if self.is_leaf(u):
            return self._labels[u]
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
    
    @memoize
    def descendant_nodes(self,u):
        return dfs_successors(self,u)
        #=======================================================================
        # try:
        #    return self._descendant_nodes[u]
        # except:
        #    if not hasattr(self,'_descendant_nodes'):
        #        self._descendant_nodes={}
        #    self._descendant_nodes[u]=dfs_successors(self,u)
        #    return self._descendant_nodes[u]
        #=======================================================================

    @memoize
    def descendant_taxa(self,u):
        return [self.label(desc) for desc in self.descendant_nodes(u) if self.is_labelled(desc)]
        #=======================================================================
        # try:
        #    return self._descendant_taxa[u]
        # except:
        #    if not hasattr(self,'_descendant_taxa'):
        #        self._descendant_taxa={}
        #    self._descendant_taxa[u]=[self.label(desc)
        #                           for desc in self.descendant_nodes(u)
        #                           if self.is_labelled(desc)]
        #    return self._descendant_taxa[u]
        #=======================================================================

    @memoize
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
        #=======================================================================
        # try:
        #    return self._strict_descendant_nodes[u]
        # except:
        #    if not hasattr(self,'_strict_descendant_nodes'):
        #        self._strict_descendant_nodes={}
        #    if self.is_root(u):
        #        self._strict_descendant_nodes[u]=self.descendant_nodes(u)
        #        return self._strict_descendant_nodes[u]
        #    pruned=copy.deepcopy(self)
        #    pruned.remove_node(u)
        #    desc_pruned=[]
        #    for root in self.roots():
        #        desc_pruned.extend(dfs_successors(pruned,root))
        #    self._strict_descendant_nodes[u]=\
        #        [desc for desc in self.descendant_nodes(u)
        #         if not desc in desc_pruned]
        #    return self._strict_descendant_nodes[u]
        #=======================================================================

    @memoize
    def strict_descendant_taxa(self,u):
        return [self.label(desc) for desc in self.strict_descendant_nodes(u) if self.is_labelled(desc)]
        #=======================================================================
        # try:
        #    return self._strict_descendant_taxa[u]
        # except:
        #    if not hasattr(self,'_strict_descendant_taxa'):
        #        self._strict_descendant_taxa={}
        #    self._strict_descendant_taxa[u]=\
        #        [self.label(desc)
        #         for desc in self.strict_descendant_nodes(u)
        #         if self.is_labelled(desc)]
        #    return self._strict_descendant_taxa[u]
        #=======================================================================

    @memoize
    def ancestors(self,taxon):
        return [u for u in self.sorted_nodes() if taxon in self.descendant_taxa(u)]
        #=======================================================================
        # try:
        #    return self._ancestors[taxon]
        # except:
        #    if not hasattr(self,'_ancestors'):
        #        self._ancestors={}
        #    self._ancestors[taxon]=[u for u in self.sorted_nodes() if
        #                            taxon in self.descendant_taxa(u)]
        #    return self._ancestors[taxon]
        #=======================================================================

    @memoize
    def strict_ancestors(self,taxon):
        return [u for u in self.sorted_nodes() if taxon in self.strict_descendant_taxa(u)]
        #=======================================================================
        # try:
        #    return self._strict_ancestors[taxon]
        # except:
        #    if not hasattr(self,'_strict_ancestors'):
        #        self._strict_ancestors={}
        #    self._strict_ancestors[taxon]=\
        #        [u for u in self.ancestors(taxon) if
        #         taxon in self.strict_descendant_taxa(u)]
        #    return self._strict_ancestors[taxon]
        #=======================================================================

    @memoize
    def CSA(self,tax1,tax2):
        return [u for u in self.ancestors(tax1) if
                (u in self.ancestors(tax2)) and
                ((u in self.strict_ancestors(tax1)) or
                 (u in self.strict_ancestors(tax2)))]

    @memoize
    def LCSA(self,tax1,tax2):
        csa=self.CSA(tax1,tax2)
        csa.sort(lambda x,y:cmp(self.height(x),self.height(y)))
        return csa[0]        

    @memoize
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

        #=======================================================================
        # try:
        #    return self._nodal_matrix
        # except:
        #    n=len(self.taxa())
        #    matrix=numpy.zeros((n,n),int)
        #    dicdist=all_pairs_shortest_path_length(self)
        #    for i in range(n):
        #        ti=self.taxa()[i]
        #        for j in range(i,n):
        #            tj=self.taxa()[j]
        #            lcsa=self.LCSA(ti,tj)
        #            matrix[i,j]=dicdist[lcsa][self.node_by_taxa(ti)]
        #            matrix[j,i]=dicdist[lcsa][self.node_by_taxa(tj)]
        #    self._nodal_matrix=matrix
        #    return self._nodal_matrix
        #=======================================================================


    def nodal_area(self):
        mat=self.nodal_matrix()
        #mat=mat+mat.transpose()
        return sum(abs(mat.flatten()))

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
                u=restricted.node_by_taxa(taxon)
                del restricted._labels[u]
        restricted.forget()

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
#                parent=restricted.predecessors(u)[0]
#                child=restricted.successors(u)[0]
#                restricted.add_edge(parent,child)
                restricted.remove_node(u)
        restricted.forget()
        return restricted

    def has_nested_taxa(self):
        taxa=self.taxa()
        for taxon in taxa:
            if not self.is_leaf(self._node_by_taxa[taxon]):
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
#            thislevel=self.leaves()
#            h=0
#            i=len(self.leaves())+1
#            while thislevel:
#                for u in thislevel:
#                    for v in self.predecessors(u):
#                        if self.height(v)==h+1 \
#                           and self._matching_representation[v]==0:
#                           self._matching_representation[v]=i
#                           i += 1
#                h += 1
#                thislevel=filter(lambda u:self.height(u)==h,self.nodes())
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

    
