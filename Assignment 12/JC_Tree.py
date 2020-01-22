
TOKENS = {
    "OPEN_SUBTREE" : "(",
    "CLOSE_SUBTREE" : ")",
    "SIBLING" : ",",
    "DISTANCE_SEPERATOR" : ":",
    "TERMINATING SYMBOL" : ";"
}


class TreeNode:
    """
    Tree nodes representing a phylogenetic tree.
    The tree should be build by calling the constructor without arguments to generate the root node
    and then calling add_child to expand it.

    Args:
        name(str, optional): Name of the node. If name is not given it will be automatically generated.
        parent(TreeNode, optional): Node which is the parent of this node. Defaults to None.
        distance(int, optional): Distance from this node to its parent. Defaults to 0.
        sequence(object, optional): Sequence which can be assigned to this node (in our implementation this is represented by a vector with integers in [0,3]). Defaults to None.
        depth(int, optional): Depth of this node. Defaults to 0.
    """
    cur_node_number = 0
    def __init__(self, name = None, parent = None, distance = 0, sequence = None, depth = 0):
        if name is None: 
            name = "unnamed node {}".format(TreeNode.cur_node_number)
            TreeNode.cur_node_number += 1
        self.parent = parent
        self.distance = distance
        self.sequence = sequence
        self.name = name
        self.depth = depth
        self.children = []

    def add_child(self, name = None, distance = 0, sequence = None):
        """
        Function to add a child to this node.

        Args:
            name(str, optional): Name of the node to add.
            distance(int, optional): Distance between this node and the node to add.
            sequence(object, optional): Sequence which is assigned to the new child.

        """
        new_child = TreeNode(name, self, distance, sequence, self.depth + 1)
        self.children.append(new_child)
        return new_child

    def __getitem__(self, i):
        return self.children[i]

    def __iter__(self):
        self.iter_index = 0
        return self
    
    def __next__(self):
        if self.iter_index < len(self.children):
            result = self.children[self.iter_index]
            self.iter_index += 1
            return result
        else:
            raise StopIteration

    def __str__(self):
        return "{}{} : {}{}".format("\t" * self.depth, self.name, self.distance, "".join(["\n{}".format(str(child)) for child in self]))


        

class Newick_Parser:
    """
    Class to parse newick strings.
    
    """
    def __init__(self):
        self._reset_state()
        self.cur_node_number = 0
  
    def _get_node_name(self, name):
        """
        Function which automatically generates a name if none is given. Otherwise just returns the given name.

        Args:
            name(str): Name to check.

        Returns:
            Automatically generated name if given name was empty, otherwise just returns the given name.

        """
        if not name:
            self.cur_node_number += 1
            return "node {}".format(self.cur_node_number - 1)
        else:
            return name

    def _reset_state(self):
        """
        Resets the state which is currently read, name or distance. This should be called when a token signifiying the end of a node is parsed.
        
        """
        self.name = []
        self.distance = []
        self.state = "Name"
            
    def parse_newick(self, newick_string):
        """
        Parses the given newick string to a TreeNode representing the root of a phylogenetic tree. Sequences will not yet be written (as those are not part of the newick string).
        
        Args:
            newick_string(str): String to be parsed.

        Returns:
            Phylogentic tree represented by the given newick string.

        """
        root_node = TreeNode(name = "root")
        cur_node = root_node
        cur_tree = cur_node
        self._reset_state()
        for character in newick_string:
            if character == TOKENS["OPEN_SUBTREE"]:
                cur_tree = cur_node = cur_node.add_child()

            elif character == TOKENS["CLOSE_SUBTREE"]:
                cur_node.name = self._get_node_name("".join(self.name))
                cur_node.distance = float("".join(self.distance))
                self._reset_state()
                cur_tree = cur_node = cur_node.parent

            elif character == TOKENS["SIBLING"]:
                cur_node.name = self._get_node_name("".join(self.name))
                cur_node.distance = float("".join(self.distance))
                self._reset_state()
                cur_node = cur_node.parent.add_child()

            elif character == TOKENS["DISTANCE_SEPERATOR"]:
                self.state = "Distance"

            elif character == TOKENS["TERMINATING SYMBOL"]:
                break

            else:
                if self.state == "Distance":
                    self.distance.append(character) 
                elif self.state == "Name":
                    self.name.append(character)
        return root_node
                

if __name__ == "__main__":
    reader = Newick_Parser()
    print(reader.parse_newick("(((f1:0.25,f2:0.25):0.25,(f3:0.05,f4:0.05):0.45):0.5,((f5:0.4,f6:0.4):0.2,(f7:0.3,f8:0.3):0.3):0.4);"))