from phylonetwork import PhyloNetwork, PhyloTree

network1 = PhyloNetwork(eNewick="((1,2),(3,4)5);")

def test001():
    assert network1.taxa() == ['1','2','3','4','5']

def test002():
    assert network1.leaves() == ['_3', '_4', '_6', '_7']

def test003():
    assert network1.label('_6') == '3'

def test004():
    assert network1.is_phylogenetic_network() == True

def test005():
    assert network1.node_by_taxa(network1.label('_6')) == '_6'

network2 = PhyloNetwork(eNewick="((1,2),(3,4)1);")

def test006():
    assert network2.nodes_by_taxa('1') == {'_5','_3'}

def test007():
    assert network2.nodes_by_taxa('xxx') == set()

network3 = PhyloNetwork(eNewick="((3),(1,2));")
network4 = PhyloNetwork(eNewick="((4,5#1)2,(#1,6)3)1;")

def test008():
    assert list(map(network3.is_tree_node, network3.nodes())) == [True, True, True, True, True, True]

def test009():
    assert list(map(network4.is_tree_node, network3.nodes())) == [True, True, True, True, True, False]

def test010():
    assert network4.is_tree_node('xxx') == False

def test011():
    assert list(map(network3.is_hybrid_node, network3.nodes())) == [False, False, False, False, False, False]

def test012():
    print(list(map(network4.is_hybrid_node, network3.nodes())))
    assert list(map(network4.is_hybrid_node, network4.nodes())) == [False, False, False, True, False, False]

def test013():
    assert network4.is_hybrid_node('xxx') == False

def test014():
    assert not network3.is_leaf('_1') and network3.is_leaf('_6') and not network3.is_leaf('non-existing node')

def test015():
    assert not network3.is_root('_3') and network3.is_root('_1')

network5 = PhyloNetwork(eNewick="((1)E,2);")

def test016():
    assert network5.elementary_nodes() == ['_2']

network6 = PhyloNetwork(eNewick="((3),(1,2))4;")

def test017():
    assert list(map(network6.is_labelled, network6.nodes())) == [True, False, True, False, True, True]

def test018():
    network = PhyloNetwork(eNewick="((3),(1,(5,6)2))4;")
    assert network.nodes() == ['_1', '_2', '_3', '_4', '_5', '_6', '_7', '_8']
    assert network.leaves() == ['_3', '_5', '_7', '_8']
    assert list(map(network.is_leaf, network.leaves())) == [True, True, True, True]

def test019():
    network = PhyloNetwork(eNewick="(1,2,3)ROOT;")
    assert network.nodes() == ['_1', '_2', '_3', '_4']
    assert network.roots() == ['_1']
    assert network.label('_1') == 'ROOT'

def test020():
    network = PhyloNetwork(eNewick="((A,B,C),1);")
    assert network.nodes()==['_1', '_2', '_3', '_4', '_5', '_6']
    assert network.labelled_nodes()==['_3', '_4', '_5', '_6']
    assert list(map(network.label, network.labelled_nodes())) == ['A', 'B', 'C', '1']
    assert sorted(network.unlabelled_nodes()) == ['_1', '_2']
    assert list(map(network.label, network.unlabelled_nodes()))==[None, None]

def test021():
    network = PhyloNetwork(eNewick="((1,2),(3,4));")
    assert network.nodes()
    assert sorted(network.interior_nodes())==['_1', '_2', '_5']
    assert list(map(network.is_leaf, network.interior_nodes()))==[False, False, False]

def test022():
    network = PhyloNetwork(eNewick="((((LEAF#1))),#1);")
    assert sorted(network.nodes()) == sorted(['#1', '_4', '_3', '_2', '_1'])
    assert list(map(network.depth, network.nodes()))==[0, 1, 2, 3, 1]
    assert network.depth('non-existing node') is None

def test023():
    network = PhyloNetwork(eNewick="((((LEAF#1))),#1);")
    assert list(map(network.height, network.nodes()))==[4, 3, 2, 1, 0]
    assert network.height('non-existing node') is None

def test024():
    from numpy import array, array_equal
    network = PhyloNetwork(eNewick="((((LEAF1, LEAF2#1)), #1)INT,#1);")
    assert network.taxa()==['INT', 'LEAF1', 'LEAF2']
    assert network.roots() == ['_1']
    assert array_equal(network.mu('_1'), array([1, 1, 3]))
    assert network.successors('_1')==['_2', '#1']
    assert array_equal(network.mu('#1'), array([0, 0, 1]))
    assert array_equal(network.mu('_2'), array([1, 1, 2]))

def test025():
    network = PhyloNetwork(eNewick="((((LEAF1, LEAF2#1)), #1)INT,#1);")
    assert network.sorted_nodes()==['#1', '_5', '_3', '_4', '_2', '_1']
    assert network.mu_string()=='[0 0 1]-[0 1 0]-[0 1 1]-[0 1 1]-[1 1 2]-[1 1 3]'

def test026():
    network = PhyloNetwork(eNewick="((,(3,4)#1)2,#1)1;")
    assert sorted(network.nodes()) == sorted(['_5', '_4', '_3', '_2', '_1', '#1'])
    assert network.node_by_taxa('2') == '_2'
    assert sorted(network.descendant_nodes('_2')) == sorted(['_5', '_4', '#1', '_3', '_2'])
    assert sorted(network.strict_descendant_nodes('_2'))==['_2', '_3']
    assert sorted(network.descendant_taxa('_2'))==sorted(['4', '3', '2'])
    assert network.strict_descendant_taxa('_2')==['2']

def test027():
    network = PhyloNetwork(eNewick="((3,4)1,(5,6,7)2);")
    assert network.ancestors('3')==['_3', '_2', '_1']

def test028():
    network = PhyloNetwork(eNewick='((, (3, 4)#1)2,#1)1;')
    assert '_2' in network.ancestors('3')
    assert not '_2' in network.strict_ancestors('3')

def test029():
    network = PhyloNetwork(eNewick='((, (3, 4)#1)2,#1)1;')
    assert network.CSA('3', '4')==['#1', '_1']
    assert network.LCSA('3', '4')=='#1'

def test030():
    network = PhyloNetwork(eNewick="(((1)), 2);")
    assert network.CSA('1', '2')==['_1']
    assert network.LCSA('1', '2')=='_1'

def test031():
    from numpy import array, array_equal
    network = PhyloNetwork(eNewick="(((1,2), 3), 4);")
    assert array_equal(network.nodal_matrix(),
        array([[0, 1, 2, 3],
               [1, 0, 2, 3],
               [1, 1, 0, 2],
               [1, 1, 1, 0]]))
    assert network.nodal_area() == 19

def test032():
    from numpy import array, array_equal

    network = PhyloNetwork(eNewick="(((1,2), 3), 4);")
    assert array_equal(network.cophenetic_matrix(),
        array([[3, 2, 1, 0],
               [0, 3, 1, 0],
               [0, 0, 2, 0],
               [0, 0, 0, 1]]))

def test033():
    network = PhyloNetwork(eNewick="((1,2), 3)4;")
    network2 = PhyloNetwork(eNewick="(1,(4,5)2);")
    assert network.common_taxa(network2) == ['1', '2', '4']
    assert network.common_taxa_leaves(network2) == ['1']

def test034():
    network = PhyloNetwork(eNewick="((1,2), (3,4));")
    assert network.taxa() == ['1', '2', '3', '4']
    assert network.roots() == ['_1']
    subnetwork = network.topological_restriction(['1', '2'])
    assert subnetwork.taxa() == ['1', '2']
    assert subnetwork.nodes() == ['_2', '_3', '_4']
    assert subnetwork.roots() == ['_2']
    subnetwork = network.topological_restriction(['1', '3'])
    assert '_1' in subnetwork.roots()

def test035():
    network = PhyloNetwork(eNewick="((2,3), (4,(5,6)))1;")
    assert network.nodes()==sorted(['_9', '_8', '_7', '_6', '_5', '_4', '_3', '_2', '_1'])
    assert network.nested_label('_7') == '{5,6}'
    assert network.nested_label('_6') == '4'
    assert network.nested_label('_5') == '{4,{5,6}}'

def test036():
    network = PhyloNetwork(eNewick="((1,2), (3,4));")
    assert network.nested_label_representation() == {'{3,4}', '1', '3', '2', '4', '{1,2}', '{{1,2},{3,4}}'}

def test037():
    from numpy import array_equal, array
    network = PhyloTree(eNewick="(((1,2),3),4);")
    print(network.nodal_matrix())
    assert array_equal(network.nodal_matrix(),
        array([[0, 1, 2, 3],
               [1, 0, 2, 3],
               [1, 1, 0, 2],
               [1, 1, 1, 0]]))
