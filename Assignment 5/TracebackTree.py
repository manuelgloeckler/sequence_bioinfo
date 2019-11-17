

#constants
MATCH = 0
X_GAP = 1
Y_GAP = 2

class TracebackNode:
    def __init__(self, value, parent = None):
        self.parent = parent
        self.value = value
        self.children = []
        if not self.parent is None:
            self.parent.add_child(self)
    
    def add_child(self, child):
        self.children.append(child)


def recreate_alignment(seq1, seq2, leaf_node):
    align1 = ""
    align2 = ""
    node = leaf_node
    i1 = 0
    i2 = 0
    while node:
        if node.value == MATCH:
            align1 += seq1[i1]
            align2 += seq2[i2]
            i1 += 1
            i2 += 1
        if node.value == X_GAP:
            align1 += "-"
            align2 += seq2[i2]
            i2 += 1
        if node.value == Y_GAP:
            align1 += seq1[i1]
            align2 += "-"
            i1 += 1
        node = node.parent
    return align1, align2