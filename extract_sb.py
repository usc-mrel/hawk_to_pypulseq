import xml.etree.ElementTree as ET
import sys

def extract_element(root, name):
    for i in range(0, len(root)):
        if root[i].text == name:
            return root[i+1].text
    return -1

def extract_element_recursive(parent, index, name):
    # some parameters have a "position" element, which is a child of a fellow parent element.
    # this makes the parsing a bit hairy and requires a double loop which is not very nice.
    # TODO: make this more elegant.
    if len(parent[index]) == 0:
        # nothing left to recurse.
        if parent[index].text == name:
            for (i, child) in enumerate(parent):
                if len(child) > 0:
                    for (j, child2) in enumerate(child):
                        if child2.text == "position":
                            return [parent[index+1].text, child[j+1].text]
            return [parent[index+1].text, -1]
        return [-1, -1];
    return [max([float(extract_element_recursive(parent[index], i, name)[0]) for i in range(0, len(parent[index]))]),
            max([float(extract_element_recursive(parent[index], i, name)[1]) for i in range(0, len(parent[index]))])]

def extract_sb_parameters(root, elements_to_extract):
    extracted_elements = []
    positions = []
    for element in elements_to_extract:
        [extract, position] = extract_element_recursive(root, 0, element)
        extracted_elements.append(extract)
        positions.append(position)

    # also return positions of the elements if they exist.
    return extracted_elements, positions
