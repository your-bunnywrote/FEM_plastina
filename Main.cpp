#include "Mesh.h"

int main() {
	Mesh mesh;
	mesh.output("nodes.txt", mesh.nodes);
	mesh.output("elements.txt", mesh.elements);
	return 0;
}