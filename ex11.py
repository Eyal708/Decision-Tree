"""
Name: Eyal Haluts
ID: 318456290
CSE username: eyal_haluts
"""

import itertools

class Node:
    def __init__(self, data, positive_child=None, negative_child=None):
        self.data = data
        self.positive_child = positive_child
        self.negative_child = negative_child

    def is_leaf(self):
        """
        :return: True if the Node is a leaf, False otherwise
        """
        return self.negative_child is None and self.positive_child is None

    def are_equals(self, other):
        """
        :param other: Node object
        :return: True if both Nodes are equal, meaning they they are roots
        of trees with identical leaves. False otherwise
        """
        if self.is_leaf() and other.is_leaf():
            return self.data == other.data
        if not self.is_leaf() and not other.is_leaf():
            return self.positive_child.are_equals(other.positive_child) and \
                   self.negative_child.are_equals(other.negative_child)
        if self.is_leaf() and not other.is_leaf():
            return self.are_equals(other.positive_child) and \
                   self.are_equals(other.negative_child)
        if other.is_leaf() and not self.is_leaf():
            return self.positive_child.are_equals(other) and \
                   self.negative_child.are_equals(other)

    def replace_node(self, other):
        """
        sets all attributes of Node self with the attributes of Node other
        other: Node object
        :return: None
        """
        self.data = other.data
        self.positive_child = other.positive_child
        self.negative_child = other.negative_child


class Record:
    def __init__(self, illness, symptoms):
        self.illness = illness
        self.symptoms = symptoms


def parse_data(filepath):
    with open(filepath) as data_file:
        records = []
        for line in data_file:
            words = line.strip().split()
            records.append(Record(words[0], words[1:]))
        return records


class Diagnoser:
    def __init__(self, root: Node):
        self.root = root

    def diagnose(self, symptoms):
        """
        :param symptoms: list of of symptoms to diagnose : list(string)
        :return: diagnosis according to the Diagnoser object: str
        """
        return self.diagnose_helper(symptoms, self.root)

    def diagnose_helper(self, symptoms, node):
        """
        :param symptoms: list of of symptoms to diagnose : list[string]
        :param node: the root/sub root of the tree : Node
        :return: the diagnosis according to the Diagnoser: string
        """
        if node.is_leaf():
            return node.data
        if node.data in symptoms:
            return self.diagnose_helper(symptoms, node.positive_child)
        else:
            return self.diagnose_helper(symptoms, node.negative_child)

    def calculate_success_rate(self, records):
        """
        :param records: list of Record objects : list[Record]
        :return: rate in which the right diagnosis was done by the Diagnoser:
        float
        """
        success_counter = 0
        if len(records) == 0:
            raise ValueError("Records list must contain at least one record")
        else:
            for record in records:
                if self.diagnose(record.symptoms) == record.illness:
                    success_counter += 1
            return success_counter / len(records)

    def all_illnesses(self):
        """
        :return: a sorted list of illnesses in the tree of the Diagnoser object
        (list is sorted by frequency of appearance- highest to lowest): List
        """
        illness_dict = {}
        self.all_illnesses_helper(self.root, illness_dict)
        # return the dictionary as a sorted list(highest value to lowest)
        return sorted(list(illness_dict), key=lambda key: illness_dict[key],
                      reverse=True)

    def all_illnesses_helper(self, node: Node, illness_dict):
        """
        :param node: root/ sub root of the Diagnoser's tree :Node
        :param illness_dict: a dictionary[(key)illness: (value)appearances]
        :return: a dictionary with al illnesses in the tree as keys and
        appearances as values:Dict
        """
        if node.is_leaf():
            if node.data is not None:
                if node.data in illness_dict:
                    # illness already appeared, add to number of appearances
                    illness_dict[node.data] += 1
                else:
                    # add a key with illness name
                    illness_dict[node.data] = 1
        else:
            self.all_illnesses_helper(node.positive_child, illness_dict)
            self.all_illnesses_helper(node.negative_child, illness_dict)

    def paths_to_illness(self, illness):
        """
        :param illness: string with a name of an illness or None.
        :return: a list of lists of all paths from the root of the Diagnoser
        to the illness:list[list]. an empty list if the illness was not found.
        """
        results_lst = []
        self.paths_to_illness_helper(self.root, illness, [], results_lst)
        return results_lst

    def paths_to_illness_helper(self, node, illness, path, results_lst):
        """
        :param node: a root/ sub root of the tree
        :param illness: a string with a name of an illness or None
        :param path: the current path taken from the root
        :param results_lst: a list of lists of all current paths found
        that lead to illness
        :return: a list containing all the possible paths to the illness:
        list[list]
        """
        if node.is_leaf():
            if node.data == illness:
                results_lst.append(path)

        else:
            self.paths_to_illness_helper(node.negative_child, illness,
                                         path + [False], results_lst)
            self.paths_to_illness_helper(node.positive_child, illness,
                                         path + [True], results_lst)

    def minimize(self, remove_empty=False):
        """
        :param remove_empty: boolean expression, if True- replace Nodes that
        lead to None with their child who doesn't lead to None(if one exists,
        if not replace with one of the children. if False, minmize the tree
        as musch as possible to lead to same diagnosis.
        :return: None
        """
        self.minimize_helper(self.root, remove_empty)

    def minimize_helper(self, node, remove_empty):
        """
        recursively minimizes the tree.
        :param node: the root/sub root of the tree
        :param remove_empty: boolean expression
        :return: None
        """
        if node.is_leaf():
            return
        self.minimize_helper(node.negative_child, remove_empty)
        self.minimize_helper(node.positive_child, remove_empty)
        if node.positive_child.are_equals(node.negative_child):
            node.replace_node(node.positive_child)
        elif remove_empty:
            if not Diagnoser(
                    node.positive_child).all_illnesses():  # empty list
                node.replace_node(node.negative_child)
            elif not Diagnoser(node.negative_child).all_illnesses():
                node.replace_node(node.positive_child)


def build_tree(records, symptoms):
    """
    :param records: a list of Record objects
    :param symptoms: a list of symptoms: list[string]
    :return: Diagnoser object built according to the list of symptoms. the
    leaves of tree are the illnesses from the record objects that most likely
    fit the symptoms in the tree.
    """
    for symptom in symptoms:
        if not isinstance(symptom, str):
            raise TypeError("All symptoms must be strings!")
    if len(symptoms) == 0:  # tree will only have one leaf
        tree_root = Node(one_leaf_tree(records), None, None)
        return Diagnoser(tree_root)
    tree_root = build_symp_tree(symptoms)
    leaves_dict = find_records_to_leaves(tree_root, records)
    for leaf in leaves_dict:
        ill_to_set = most_common_word(leaves_dict[leaf])
        leaf.data = ill_to_set
    return Diagnoser(tree_root)


def one_leaf_tree(records):
    """
    :param records: a list of Record objects
    :return: the data of a root with only one leaf- the most common illness
    from the records list.
    """
    illness_lst = []
    for record in records:
        illness_lst.append(record.illness)
    return most_common_word(illness_lst)


def most_common_word(lst):
    """
    :param lst: a list of objects
    :return: the most common object in the list
    """
    if not lst:
        return None
    items_dict = {}
    for item in lst:
        if item not in items_dict:
            items_dict[item] = 0
        else:
            items_dict[item] += 1
    return max(items_dict, key=lambda key: items_dict[key])


def find_records_to_leaves(tree_root, records):
    """
    :param tree_root: a root of a binary tree
    :param records: a list of record objects
    :return: a dictionary. keys: leaves that are at the end of a path of
    at least one record object.values: a list of all illnesses names that got
    to the leaf(determined by the the symptoms attribute of each Record obj)
    """
    leaves_dict = {}  # a dictionary. keys: all leaves that are at the end
    # of a path that fits a record's illness and symptoms. values: list of
    # strings of illnesses that got to the leaf.
    for record in records:
        if not isinstance(record, Record):
            raise TypeError("All records must be instances of class Record!")
        node = tree_root
        while not node.is_leaf():
            if node.data in record.symptoms:
                node = node.positive_child
            else:
                node = node.negative_child
            if node.is_leaf():
                if node not in leaves_dict:
                    leaves_dict[node] = [record.illness]
                else:
                    leaves_dict[node].append(record.illness)
    return leaves_dict


def build_symp_tree(symptoms_lst):
    """
    :param symptoms_lst: a list of symptoms: list[string]
    :return: the root of a binary tree containing
    all the symptoms from the list
    """
    if len(symptoms_lst) == 1:
        return Node(symptoms_lst[0], Node(None), Node(None))
    return Node(symptoms_lst[0], build_symp_tree(symptoms_lst[1:]),
                build_symp_tree(symptoms_lst[1:]))


def optimal_tree(records, symptoms, depth):
    """
    :param records: a list of Record objects
    :param symptoms: a list of symptoms: list[string]
    :param depth: an integer representing the depth of the tree.
    :return: Diagnoser object, who's root is the root of the tree with the
    best success rate, from all the trees with depth given.
    """
    if has_duplicates(symptoms):
        raise TypeError("symptoms can not have duplicates!")
    if depth < 0 or depth > len(symptoms):
        raise TypeError("depth must be between 0 and number of symptoms!")
    optimal_rate = 0
    best_tree = None
    for combination in itertools.combinations(symptoms, depth):
        tree_to_create = build_tree(records, combination)
        success_rate = tree_to_create.calculate_success_rate(records)
        if success_rate >= optimal_rate:
            best_tree = tree_to_create
            optimal_rate = success_rate
    return best_tree


def has_duplicates(lst):
    """
    :param lst: a list
    :return: True if the list has duplicate items, false otherwise
    """
    no_duplicates = set(lst)
    if len(no_duplicates) != len(lst):
        return True
    return False


