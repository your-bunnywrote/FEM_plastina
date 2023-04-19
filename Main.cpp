#include "Mesh.h"

int main() {
	Mesh mesh;
	vector<Point> nodes = mesh.fill_nodes();
	vector<Element> elements = mesh.fill_elements();
	mesh.output("..\\nodes.txt", nodes);
	mesh.output("..\\elements.txt", elements);
	mesh.edit_nodes(nodes);
	mesh.edit_elements(elements);
	mesh.output("..\\nodes1.txt", nodes);
	mesh.output("..\\elements1.txt", elements);


	return 0;
}