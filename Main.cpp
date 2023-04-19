#include "Mesh.h"

int main() {
	Mesh mesh;
	vector<Point> nodes = mesh.fill_nodes();
	vector<Element> elements = mesh.fill_elements();
	mesh.output("..\\nodes.txt", nodes);
	mesh.output("..\\elements.txt", elements);
	vector<Point> newNodes = mesh.edit_nodes(nodes);
	vector<Element> newElements = mesh.edit_elements(elements);
	mesh.output("..\\nodes1.txt", newNodes);
	mesh.output("..\\elements1.txt", newElements);


	return 0;
}