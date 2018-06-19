from unittest import TestCase

from .classes import PhyloNetwork

class PhylonetworkTest(TestCase):
  def setUp(self):
    self.network = PhyloNetwork(eNewick='((,(3,4)#1)2,#1)1;')

  def test_descendant_nodes_returns_list_of_nodes(self):
    self.assertEqual(self.network.descendant_nodes('_2'),
                     ['_5', '_4', '#1' '_3', '_2'])
